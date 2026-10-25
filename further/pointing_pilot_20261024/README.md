# NGC4689 approximate pointing-region pilot, 2026-10-24

The cube variance supports a rotated 2-column, 3-row layout. These are inferred
geometric regions, not recovered individual observations or verified OB IDs.

## Inputs and method

- Setonix: `/scratch/pawsey1308/mauve/cubes/v3tk/`, NGC4689, IC3392 and NGC4606
  `*_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk.fits`.
- Read nine nearest wavelength planes at each of 5100, 5500, 6000, 6800 A.
  Retain finite DATA, finite positive STAT, DQ==0; take per-pixel medians.
  These are median per-channel variances, not variance of a combined spectrum.
- Fit a minimum-area external valid-data bounding rectangle to NGC4689.
- In its rotated coordinates, median log-variance profiles across the footprint
  show one u trough and two v troughs. Detect prominent narrow depressions
  relative to a 35-pixel Gaussian trend; check all four windows independently.
  This pilot explicitly assumes one common orientation and rectangular layout.
- Region labels 1-6 run bottom to top, left to right in the displayed image.
  Exclude 20 pixels (4 arcsec) either side of each trough center and from the
  outer bounds. This is a conservative chosen margin, not an uncertainty fit.

## Outputs

- `variance_windows.png`: variance maps for target and single-pointing controls.
- `NGC4689_inferred_pointings.png`: inferred centers, interior masks and profiles.
- `NGC4689_approximate_interiors.fits`: 2D labels with celestial WCS; zero means
  excluded/invalid, 1-6 identify approximate interiors.
- `geometry.json`: transform, trough locations per window and region sizes.
- `NGC4689_approximate_shared_footprints.fits`: separate binary R1-R6 extensions
  with overlapping membership; MEMBERS counts 1, 2 or 4 candidate pointings.
  Overlap bands use the union of half-depth variance trough widths across the
  four windows. These are approximate support estimates, not exact exposure edges.
- `NGC4689_shared_footprints.png`: footprint contours and membership counts.
- `sampled_maps.npz`: small extracted maps plus input headers.

Reproduce local analysis with:

```sh
MPLCONFIGDIR=/private/tmp/ngc4689-mpl XDG_CACHE_HOME=/private/tmp/ngc4689-cache /opt/miniconda3/envs/ICRAR/bin/python analyse_maps.py
```

`extract_maps.py` runs on Setonix in the existing nGIST scientific container;
it emits the NPZ to stdout without modifying remote files. Full science products
were not downloaded or rewritten. Runtime dependencies: numpy, scipy, astropy,
matplotlib. The analysis requires the extracted NPZ beside the script.

## Interpretation and next step

Repeated narrow variance depressions and the six approximately square tiles
are strong evidence for overlap geometry. This does not prove sky residuals
or uniquely recover the contributing exposures. The header's six OB groups
provide consistent supporting provenance; R1-R6 cannot be matched to OB IDs.

The interior masks make exploratory regional continuum ratios possible, using
fixed valid spatial support across wavelength and a stated normalization band.
Different regions contain different stellar populations, dust and brightness;
their ratios need not be flat or unity even for perfect reduction. Merged-cube
overlap spectra cannot be separated back into the individual input spectra.
Use spatially matched adjacent boundary strips and observed-continuum/residual
diagnostics to test for discontinuities, alongside broad interior comparisons.
The fitted CONT model alone may already absorb continuum mismatch.

Overlap pixels belong to both neighbouring pointing footprints (and four at
crossings). The separate interior labels intentionally omit those pixels for
extraction; they must not be interpreted as exclusive full-footprint ownership.

Do not generalize this detector to all galaxies without checking rotations,
irregular/missing exposures, variance patterns, and boundary-mask sensitivity.
No all-galaxy or dual-wavelength-run notebook has been created yet.
