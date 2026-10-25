# Setonix all-galaxy continuum-ratio notebook

Primary notebook: `/home/rhuang/20261024_pointing_continuum_diagnostics.ipynb`

Remote outputs and execution helpers:
`/home/rhuang/20261024_pointing_continuum_diagnostics/`

Local source/mirror: this directory. `build_notebook.py` embeds all helpers from
`pointing_workflow.py` into a standard 19-cell notebook. No external helper import
is required when using the notebook interactively. Configuration, equations,
AO safety, assumptions, status rules, demonstrations and interpretation are
documented in the notebook itself.

## Run on Setonix

```sh
sbatch /home/rhuang/20261024_pointing_continuum_diagnostics/run_notebook.slurm
```

The Slurm script uses the existing nGIST container read-only, one CPU, 32 GB,
and a six-hour limit. The remote environment lacks nbformat/nbclient, so the
supplied `execute_notebook.py` runs cells directly and records streams, tables,
explicit figure displays and exceptions into the standard notebook JSON. It
saves after each cell. A failed cell propagates a nonzero exit to Slurm.

## Defaults and products

- Baseline CONFIG discovery: both `v3tk_v7.6.8_7000` and `v3tk_v7.6.8`.
- Special SFH16/SNR tests excluded. CONFIG is checked against nominal run limits.
- Detection windows: 5100, 5500, 6400, 6800 A. AO interval 5770-6050 A excluded.
- Channel stride 4 (typically 5 A); all retained interior spaxels integrated.
- Geometry allows missing cells and column-specific staggered horizontal joins.
- Shared footprints retain multiple memberships; ratios use trimmed interiors.
- Normalize each integrated spectrum over observed-frame 5400-5600 A.
- Spatial support is fixed across commonly supported sampled wavelengths. A
  deterministic preflight drops channels with <50% finite model coverage in
  any region, then requires each retained pixel finite across remaining channels.
  This handles unsupported NaN endpoints without varying area with wavelength.
- Per-galaxy geometry PNG/FITS/JSON, adjacent-pair ratios PNG/NPZ, inventory CSV,
  status CSV and descriptive ratio-summary CSV. Caches carry source fingerprints.

Geometry statuses: `single_candidate` skips ratios; `multiple_candidate` permits
exploratory ratios but still requires visual inspection; `review_required` does
not force segmentation. Run errors and insufficient spectral support are explicit.

## Verification and limits

Seven focused tests cover AO exclusion, shared two/four-way overlaps, single
pointings, staggered joins, ratio normalization/invalid denominators, constant
spatial support, and globally unsupported endpoint handling. NGC4689 and NGC4321
variance maps were visually inspected; detector finds 6 and 11 regions, respectively.
IC3392 and NGC4606 were classified as single candidates.

The first execution exposed the NaN endpoint issue and was stopped; the corrected
version changes its spectrum-cache key, retaining prior diagnostic files without
reusing invalid caches. Read `status.csv` and Slurm logs for the actual population
execution state rather than treating notebook creation or job submission as proof
of completed results.

No input science products are modified. Approximate masks are not exact exposure
support, verified OB identities, or independent pointing spectra. Regional ratios
also reflect populations, dust, and fitting-model differences. Retained pixels can
differ between runs. Validate geometry and support before scientific interpretation.
