"""Build a self-contained, standard Jupyter notebook with embedded helpers."""
from pathlib import Path
import nbformat as nb

root=Path(__file__).parent
source=(root/'pointing_workflow.py').read_text()
cells=[]
def md(s): cells.append(nb.v4.new_markdown_cell(s.strip().replace(chr(92)*2,chr(92))))
def code(s): cells.append(nb.v4.new_code_cell(s.strip()))

md('''# Approximate pointing geometry and continuum-ratio diagnostics
**MAUVE / PHANGS-MUSE | 2026-10-24 | Setonix**

This exploratory notebook infers approximate pointing regions from repeated
structures in the input cube variance, then compares integrated fitted-continuum
spectra between neighbouring regions in the **7000 and 8900 runs**.

The motivating observation is reduced pointing structure in NGC4689's inferred
stellar M/L after increasing SFH MDEG from 12 to 16. Here we test continuum
consistency; these regional comparisons alone do **not** establish sky-subtraction
causality or prove a particular polynomial degree is unbiased.

## Contents
1. Configuration, inputs and AO exclusion
2. Geometry algorithm and candidate-status rules
3. FITS discovery, extraction and constant spatial support
4. Output and plotting helpers
5. Inventory and provenance
6. NGC4689 / NGC4321 demonstrations and single-pointing controls
7. Full population run, status audit and ratio gallery
8. Interpretation and follow-up checks

All analysis functions are included below: no external helper module is required.
Use a scientific Python environment with numpy, scipy, astropy, matplotlib,
pandas and pyyaml. Run the full population on a compute node. Inputs are read-only;
only this notebook's dedicated output directory is written.''')

md('''## 1. Configuration and AO safety
The roots below label the standard 7000 and 8900 runs. Each product's `CONFIG`
provides its actual input cube, wavelength range, AO reader and continuum MDEG.
Special SFH16/SNR test runs are excluded from the baseline population.

**AO exclusion: observed-frame 5770-6050 Angstrom, inclusive.** This matches
`specMask_KIN_narrow10_AO` / `specMask_GAS_narrow10_AO` in the local nGIST setup
(entry `5910 280 sky AO_global_NaN_5770.00_6050.00`). NGC4321 uses
`MUSE_WFMAON`; its cube header reports `WFM-AO-N`.

Detection windows are 5100, 5500, 6400 and 6800 A, each +/-5 A, so none intersects
the AO interval. AO ratio curves leave a gap; no interpolation across it is used.

`SPECTRAL_STRIDE=4` samples every fourth native channel (normally 5 A), adequate
for the broad continuum-shape diagnostic. This is channel sampling, not smoothing;
set it to 1 for all native channels. Every valid interior spaxel is integrated.
No random spatial subsampling is used. `GALAXIES=None` processes all discovered
baseline galaxies, including combined cube IDs as named in CONFIG.''')
code('''from pathlib import Path
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')

PRODUCT_ROOT = Path('/scratch/pawsey1308/mauve/products')
OUTPUT = Path('/home/rhuang/20261024_pointing_continuum_diagnostics')
GALAXIES = None
SPECTRAL_STRIDE = 4
SEAM_PROMINENCE_DEX = 0.055
OUTPUT.mkdir(parents=True, exist_ok=True)

def show_png(path, width=1100):
    path=Path(path)
    if not path.exists():
        print('Figure unavailable:', path.name)
        return
    if '_record_display_png' in globals():
        _record_display_png(path)
    else:
        from IPython.display import display, Image
        display(Image(filename=str(path), width=width))

def show_table(frame):
    if '_record_table' in globals():
        _record_table(frame)
    else:
        from IPython.display import display
        display(frame)
print('Outputs:', OUTPUT, '| spectral stride:', SPECTRAL_STRIDE)''')

md('''## 2. Geometry: what is measured and what is assumed
For each window we take the median of valid `STAT` planes, using finite `DATA`,
finite positive variance and `DQ==0` when available. This is a median **per-channel
variance map**, not an uncertainty on a wavelength-integrated measurement.

1. Fit the external valid-data convex hull with a minimum-area rectangle to
   estimate a common orientation, without assuming a 2x3 grid.
2. Form median log-variance profiles along its two axes. Subtract a broad
   Gaussian spatial trend; narrow depressions **or enhancements** can indicate
   joins. Candidate cuts must recur in at least three of four windows.
3. Use an approximately 60 arcsec single-pointing scale as a loose geometric
   prior. IC3392 and NGC4606 provide empirical controls. This is not a precise
   field-stop model or proof of exposure identity.
4. Refine horizontal joins separately within columns, allowing staggered mosaics
   such as NGC4321. Remove empty grid cells. The detector still assumes broadly
   rectangular regions at a common orientation; arbitrary mosaics need review.
5. Define **interiors** by trimming 4 arcsec from candidate joins and the outer
   rim. Define separate **shared footprints** by expanding cells 4 arcsec into
   neighbouring cells. Overlap pixels belong to both neighbours; crossings can
   belong to four. This chosen half-width is approximate, not measured exposure
   support. The diagnostic map makes both definitions visible.

Statuses are `single_candidate` (skip ratios), `multiple_candidate` (exploratory
ratios), or `review_required` (do not force an answer). A large mosaic with no
detected cuts is not silently called single. Even `multiple_candidate` requires
visual inspection; it is not an automatic science-quality approval. R labels are
geometric labels and cannot be assigned to specific OB IDs from these files.''')
first=source.index('def sample_cube')
second=source.index('def map_interiors')
third=source.index('def plot_geometry')
code(source[:first])

md('''## 3. Discovery, cache provenance and fixed-area spectral integration
Discovery uses live product `CONFIG` files, accepts combined galaxy IDs, checks
the actual continuum wavelength limit, and reports absent products explicitly.
Only one input mosaic is segmented for both runs; a run using a different input
cube is flagged for review. AO flags come from CONFIG and input metadata.

Small map caches are keyed by input path, file size, modification time and window
settings. Spectrum caches also include geometry and channel stride. Delete the
output `cache/` after modifying helper algorithms beyond these settings.

Spatial masks are mapped through celestial WCS. Integer pixel offsets are
supported; non-integer reprojection is rejected for explicit review. For each
region and run, a spaxel is retained only if finite at **every sampled supported
wavelength**. First, a deterministic preflight of up to 1024 evenly spaced
interior pixels per region requires at least 50 percent finite model coverage
in every region for a channel to be used. This excludes unsupported interpolated
endpoints and leaves them blank, rather than allowing one all-NaN channel to
reject every otherwise good spaxel. The coverage fractions and supported-channel
mask are saved in the spectrum cache. Summing the same support at each wavelength avoids false ratio
features caused by changing extraction area. AO wavelengths are excluded from
that support test. At least 100 retained spaxels are required. Counts and retained
fractions are recorded; they can differ between the two runs, so inspect them.

Memory use is bounded by processing one region at a time. Large compressed
continuum cubes are flagged rather than silently expanded in memory.''')
code(source[first:second])
code(source[second:third])

md(r'''## 4. Ratio definition, plots and checkpoints
For region $i$, the extracted continuum is the sum over its fixed valid interior
support, $F_i(\lambda)=\sum_{p\in I_i} C_p(\lambda)$. Let
$m_i=\mathrm{median}_{5400\leq\lambda\leq5600}F_i(\lambda)$.
The plotted shape ratio is

$$Q_{ij}(\lambda)=\frac{F_i(\lambda)/m_i}{F_j(\lambda)/m_j}.$$

The normalization removes the overall brightness scale; it does not remove
population or dust gradients. Non-positive/non-finite denominators and invalid
normalizations produce NaNs, not fabricated ratios. `ratio_summary.csv` also
records the median raw integrated-flux ratio. NPZ files retain raw spectra and
the normalized curve. These are fitted model continua, not independent raw
observations or uncertainty-weighted sky-residual estimates.

Only neighbours sharing a substantial edge are compared; corner-only contacts
are excluded. Both runs use the same region geometry. Side-by-side panels show
7000 and 8900, with the actual observed wavelength axis preserved and AO gaps
shaded. Plot limits adapt to each curve; inspect axis ranges before comparing
amplitudes. Percentile summaries are descriptive, not formal errors.

Each galaxy writes geometry PNG/FITS/JSON, ratio PNG/NPZ and checkpointed status
tables. Errors are visible and do not masquerade as single-pointing detections.''')
code(source[third:])

md('''## 5. Inventory and sanity checks
Inspect the CONFIG-derived run limits and continuum/SFH polynomial degrees.
Changing **SFH** MDEG alone does not change a pre-existing **CONT** cube; the
NGC4689 SFH16 test had CONT disabled and is intentionally not a second continuum
measurement in this baseline comparison.''')
code('''inventory = discover_products(PRODUCT_ROOT)
show_table(pd.DataFrame(inventory))
assert not allowed_wavelength(np.array([5770.,5890.,6050.]), True).any()
assert all(not (5770 <= w <= 6050) for w in WINDOWS)
print('AO-window guards passed; detection windows:', WINDOWS)''')

md('''## 6. Demonstrations: two mosaics and two controls
Run the four demonstrations first. NGC4689 should resemble a six-region mosaic;
NGC4321 has an irregular/staggered footprint and 11 pointing groups in the input
header. Those counts are validation context, not detector inputs. IC3392 and
NGC4606 should be single-pointing candidates and skip spectral ratios.

For every target, examine all four variance maps, interior labels and shared
membership. A visually wrong boundary must be corrected or the galaxy excluded,
even when the automatic status is `multiple_candidate`.''')
code('''demo_status, demo_metrics = run_population(
    PRODUCT_ROOT, OUTPUT, galaxies=['NGC4689','NGC4321','IC3392','NGC4606'],
    stride=SPECTRAL_STRIDE, prominence=SEAM_PROMINENCE_DEX)
show_table(pd.DataFrame(demo_status))
for galaxy in ['NGC4689','NGC4321','IC3392','NGC4606']:
    print(galaxy)
    show_png(OUTPUT/f'{galaxy}_geometry.png')''')

md('''## 7. Population execution and audit
This pass discovers all baseline galaxies unless `GALAXIES` is set above. The
demonstrations reuse their cached maps/spectra. Missing products, unsupported
geometry and extraction errors are retained in the status table. Outputs are
exploratory until their geometry maps have been inspected. This may take
substantial filesystem I/O; the supplied Slurm runner uses a compute node.

Checkpoint files: `inventory.csv`, `status.csv`, `ratio_summary.csv`. When running
interactively, interrupting the kernel leaves already completed galaxy outputs
and caches intact. Re-running uses the same deterministic masks and support.''')
code('''all_status, all_metrics = run_population(
    PRODUCT_ROOT, OUTPUT, galaxies=GALAXIES, stride=SPECTRAL_STRIDE,
    prominence=SEAM_PROMINENCE_DEX)
status_frame = pd.DataFrame(all_status)
show_table(status_frame)
show_table(pd.DataFrame(all_metrics))
if not status_frame.empty:
    print(status_frame.groupby(['stage','status']).size().to_string())
    problems=status_frame[status_frame['status'].isin(['error','review_required','insufficient_support'])]
    print('Rows needing attention:', len(problems))
    if not problems.empty: show_table(problems)''')

md('''### Per-galaxy gallery: geometry followed by both-run ratios
Each candidate mosaic appears below with its inferred regions and adjacent-pair
ratios. Single-pointing candidates have no ratio plot. Missing or failed run
panels remain labelled, rather than being silently omitted. Full-resolution
figures are saved in the output directory for zooming.''')
code('''for galaxy in sorted({r['galaxy'] for r in all_status
                      if r.get('stage')=='geometry' and r.get('status')=='multiple_candidate'}):
    print('GALAXY:', galaxy)
    show_png(OUTPUT/f'{galaxy}_geometry.png')
    show_png(OUTPUT/f'{galaxy}_ratios_7000_8900.png')''')

md('''## 8. Interpretation and follow-up
- A repeated variance grid supports approximate overlap geometry. It does not
  recover individual input spectra from the merged mosaic.
- Wavy normalized regional ratios may reflect reduction differences, population
  gradients, dust or the fitting model. They are not a unique sky diagnostic.
- Compare the same curve across 7000/8900 only over their common wavelength
  range, accounting for spatial-support and continuum-fit differences.
- Check retained fractions, geometry images, alternative seam thresholds and
  2/4/6 arcsec interior margins before interpreting a feature.
- For a stronger discontinuity test, compare closely spaced strips on either
  side of a join, control for galaxy gradients, and inspect observed spectra
  and data-minus-fit residuals. These are future diagnostics, not implemented
  evidence in this notebook.
- The shared-footprint masks permit overlap visualisation; adding the same
  merged overlap spectrum into both extracted regions cannot reproduce two
  independent pointing measurements and would artificially couple their ratios.
  For that reason the default ratios use only trimmed interiors.

**Reproducibility:** retain this executed notebook, inventory, geometry JSON/FITS,
status tables, spectrum caches and input products. No original science products
are modified. Review all flagged cases before any population-level conclusions.''')

notebook=nb.v4.new_notebook(cells=cells,metadata=dict(kernelspec=dict(display_name='Python 3 (nGIST scientific environment)',language='python',name='python3'),language_info=dict(name='python',version='3.11')))
nb.validate(notebook)
for i,cell in enumerate(notebook.cells,1):
    if cell.cell_type=='code': compile(cell.source,f'notebook cell {i}','exec')
target=root/'20261024_pointing_continuum_diagnostics.ipynb'
nb.write(notebook,target)
print(target,len(cells),'cells')
