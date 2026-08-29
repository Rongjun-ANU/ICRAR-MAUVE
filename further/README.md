# MAUVE further post-processing

This folder contains the active MAUVE second-project post-processing scripts
that run on nGIST products under `v3tk_v7.6.8/{GALID}`.

## PHANGS-native cube inputs

`Mass.py` and `proxy_EWHa.py` read full datacubes from a local filesystem path.
By default the shell wrappers point `--cube-root` at:

```text
/arc/projects/mauve/cubes/v3tk
```

The three MAUVE galaxies that overlap the PHANGS-MUSE DR1 native public release
use these cube filenames:

```text
NGC4254_PHANGS_DATACUBE_native.fits
NGC4321_PHANGS_DATACUBE_native.fits
NGC4535_PHANGS_DATACUBE_native.fits
```

These are expected to be staged locally, for example on Setonix scratch, before
full Mass or proxy-EW processing. Do not point these scripts at a `vos:` URI for
full-cube processing. Use CANFAR tools such as `vcat --head` or cutouts only for
header-only or small inspection jobs; full FITS integration should use a local
staged file.

## Running selected galaxies

Both wrappers accept positional GALID arguments. With no arguments they run the
default sample, including the three PHANGS-native galaxies. With arguments they
run only the requested subset:

```bash
./Mass.sh NGC4254
./Mass.sh NGC4254 NGC4321 NGC4535

./proxy_EWHa.sh NGC4254
./proxy_EWHa.sh NGC4254 NGC4321 NGC4535
```

The Python scripts can also be run directly with explicit roots:

```bash
python Mass.py -g NGC4254 --root . --fallback-root . --cube-root /arc/projects/mauve/cubes/v3tk
python proxy_EWHa.py -g NGC4254 --root . --fallback-root . --cube-root /arc/projects/mauve/cubes/v3tk
```

## SFR, strong-line metallicity, and photoionisation-grid methods

`SFR+Z.py` reads the raw gas-bin maps plus the matching
`*_spatial_binning_maps_further.fits` and spatial mask. It writes
`*_gas_bin_maps_further.fits` back to the selected product directory. The
default extinction law remains Calzetti (2000).

Run one galaxy with the ICRAR environment as follows:

```bash
/opt/miniconda3/envs/ICRAR/bin/python SFR+Z.py \
  -g NGC4254 --root . --fallback-root . --product-subdir v3tk_v7.6.8
```

The pyqz, NebulaBayes, and JY22 products are enabled by three configuration
toggles near the beginning of `SFR+Z.py`. The implementation uses these local
project files and the released Peng2026 grid:

```text
SFR+Z.py
model_grid_compat.py
model_grid_diagnostics.py
Peng2026/photoionization_models/photoionization_grid_interpolated.fits
```

`model_grid_compat.py` is required only for pyqz 0.8.4.0 and NebulaBayes 0.9.9
under the current Python, NumPy, SciPy, and Matplotlib stack. It verifies exact
supported package versions (and exact pyqz source hashes), creates a disposable
pyqz runtime overlay, and supplies narrow NebulaBayes API aliases. It does not
modify `site-packages`. JY22 does not use that compatibility layer; its inference
is implemented in `model_grid_diagnostics.py` and fails closed if the released
grid contract is malformed.

### Shared scientific selection

All methods use the same existing HII and SF masks and the six fitted lines
`HB4861`, `OIII5006`, `HA6562`, `NII6583`, `SII6716`, and `SII6730`. pyqz and
NebulaBayes receive the existing dust-corrected fluxes and propagated errors.
JY22 instead starts from the corresponding raw fluxes and independent raw
errors so it can construct the published raw close-pair ratios and
their shared-denominator covariance. Model names map `OIII5006 -> OIII5007`,
`NII6583 -> NII6584`, `SII6716 -> SII6717`, and `SII6730 -> SII6731` where the
published grids use the alternate wavelength convention. A spectrum is
eligible only when all six inputs for that method are finite and strictly
positive inside the existing HII or SF mask. pyqz/NebulaBayes apply this check
to corrected maps; JY22 applies it independently to raw maps. This is an
input-domain check, not an additional S/N, EW, or three-dimensional selection.

The script evaluates each ascending spatial `BINID` once on the union of the
HII and SF masks and broadcasts the shared inference back to the map. The HII
and SF HDUs are then masked views of the same underlying result, so overlapping
pixels cannot acquire different model values. No all-gas model-grid value is
calculated because these grids describe HII regions.

Integrated HII and SF values remain in the run log only. For each region, the
script selects one common raw six-line aperture, sums each line over exactly
those pixels, and sums raw errors in quadrature. pyqz/NebulaBayes apply one
integrated Balmer decrement (with the existing 2.86 floor) and receive the
resulting corrected spectrum. JY22 instead derives its
raw N2/S2/R3 ratios and covariance directly from the raw integrated sums,
without using the integrated Balmer correction.
The script does not average spatial metallicities.

### Locked model settings

pyqz uses its bundled MAPPINGS V 5.0.16 plane-parallel grid with
`[NII]/[SII]+;[OIII]/[SII]+`, `log(P/k)=5.0`, `kappa=inf`, grid sampling 2,
normal input errors, 800 Monte Carlo samples, multivariate KDE on a 101-by-101
grid, seed 20260823, and one process. `[SII]+` is the 6716+6730 flux with its
component errors added in quadrature. The stored abundance and `LogQ` are the
combined KDE estimates. `LOG_U_PYQZ = LOG_Q_PYQZ - 10.47712125472`. The pyqz
pressure is fixed model provenance; it is not an inferred pressure map. Its
bundled range is `8.00 <= 12+log(O/H) <= 8.875` and
`6.5 <= LogQ <= 8.5`, with diagnostic-specific off-grid behaviour retained.

NebulaBayes uses its bundled MAPPINGS 5.1 HII grid with parameter order
`log U`, `log P/k`, `12+log(O/H)`, interpolated shape `[40,20,160]`, linear
interpolation, `grid_error=0.10`, an Hbeta normalisation, all six likelihood
lines, a uniform prior, and `deredden=False` because the script supplies
corrected fluxes. Its grid ranges are `7.061 <= 12+log(O/H) <= 9.304`,
`-4 <= log U <= -2`, and `4.2 <= log(P/k) <= 8.6`.

Before calling NebulaBayes, the script normalises each line to Hbeta and
propagates the Hbeta denominator error under the independent-error
approximation. NebulaBayes cannot represent the covariance shared through that
denominator or through the common Balmer-decrement correction. The reported
best-model reduced chi-square uses observational errors, is evaluated at the
tuple of separate sampled marginal modes, and does not include the separate
10-percent model-grid error term. It is therefore retained as a package
diagnostic rather than treated as the objective that generated the posterior.

The main NebulaBayes O/H, log U, and log(P/k) products are quadrature-normalised
one-dimensional marginal posterior means. These continuous summaries remove
the grid-node banding inherent to the package's sampled marginal modes. The
original sampled modes remain in explicit `*_MODE` HDUs. Parameter-specific
`NB_*_RELIABLE` maps require a finite mean, finite CI68 bounds, and a clear edge
bit; in particular they expose rather than smooth over low-pressure saturation.

JY22 uses the released Peng et al. (2026) interpolated Ji & Yan (2022) grid at
`Peng2026/photoionization_models/photoionization_grid_interpolated.fits`. Its
SHA-256 is
`d7a219b60c9a1ea8b29339988c84b1832028b28b3c417d1c3c58b420831eb38a`.
The released file is 40 by 40 in `log(Z/Zsun)` and log U. Following Peng et
al. (2026), the released line-ratio surface is linearly resampled to a 200 by
200 inference grid over the documented `-4.0 <= log U <= -0.5` range rather
than fitted only at the archive nodes or extended over the file's full
`-4.0 <= log U <= +1.0` extent. This avoids artificial 0.046-dex O/H banding
when high-S/N likelihoods collapse onto individual archive nodes. Oxygen abundance is
`12+log(O/H) = log(Z/Zsun) + 8.69` on the native **pre-depletion** Ji-Yan
scale. Ji & Yan's fixed oxygen depletion gives
`O/H_post = O/H_pre - 0.22 dex`; the main `O_H_JY22` maps use this
post-depletion gas-phase scale, while `O_H_JY22_PREDEP` retains the grid scale.

The JY22 likelihood uses the log ratios N2=`[N II]6583/Halpha`,
S2=`([S II]6716+[S II]6730)/Halpha`, and R3=`[O III]5006/Hbeta`, in that
internal order. The two raw [S II] components are summed before forming S2. A
first-order Jacobian propagates all six independent raw flux
errors through the raw ratios and shared Halpha denominator, and the likelihood
uses the resulting full 3-by-3 covariance. No Balmer-decrement correction enters
N2, S2, or R3. The prior is flat over
the retained grid. Central O/H and log U values are posterior means; uncertainty
maps are marginal equal-tailed 16th and 84th percentiles interpolated through
the ordinary cumulative marginal posterior. The MaNGA-specific
1.25 error inflation used in Ji & Yan validation is not transferred here
(`JY22_ERROR_INFLATION=1.0`).

### New FITS products and QC

The new groups are inserted after the strong-line/C20 products and before the
PyNeb density/temperature products. Every item is an immediate HII/SF pair.

| Method | Products |
| --- | --- |
| pyqz | `O_H_PYQZ_*`, `LOG_Q_PYQZ_*`, `LOG_U_PYQZ_*` and their `*_ERR` maps |
| pyqz QC | `PYQZ_FLAG_*`, `PYQZ_RS_OFFGRID_*`, `PYQZ_VALID_*` |
| NebulaBayes | `O_H_NEBULABAYES_*`, `LOG_U_NEBULABAYES_*`, `LOG_PK_NEBULABAYES_*`, each with `*_CI68_LOW`, `*_CI68_HIGH`, and `*_MODE` |
| NebulaBayes QC | `NB_CHI2_RED_*`, `NB_NLOCALMAX_*`, `NB_FLAG_*`, `NB_OH_RELIABLE_*`, `NB_LOGU_RELIABLE_*`, `NB_LOGPK_RELIABLE_*`, `NB_VALID_*` |
| JY22 | Post-depletion `O_H_JY22_*`, pre-depletion `O_H_JY22_PREDEP_*`, their `*_16`/`*_84` maps, and `LOG_U_JY22_*` with `*_16`/`*_84` |
| JY22 QC | `JY22_CHI2_MIN_*`, `JY22_RESID_N2_*`, `JY22_RESID_S2_*`, `JY22_RESID_R3_*`, `JY22_RESID_NORM_*`, `JY22_FLAG_*`, `JY22_FIT_OK_*`, `JY22_VALID_*` |

Here `*` is `HII` or `SF`. pyqz error maps store the maximum half-extent of
the peak-normalised 0.61 KDE contour; they should not be relabelled as generic
Gaussian 1-sigma errors. NebulaBayes main maps store marginal-posterior means,
with sampled modes and equal-tailed 68-percent bounds retained separately. JY22 maps
store posterior means and separate marginal p16/p84 bounds; `JY22_CHI2_MIN` is
the minimum full-covariance chi-square, not a reduced chi-square.
The primary FITS header records the inference sampling in `JY22SHP` and the
surface interpolation method in `JY22INT`. Finer sampling removes numerical
banding but does not cure a spectrum that lies away from the model surface;
`JY22_FIT_OK` and `JY22_FLAG` remain mandatory downstream QC.

The log-only integrated HII/SF pyqz calculation first runs the unchanged native
pyqz global KDE. If that result is non-finite because a high-S/N Monte Carlo PDF
is narrower than its fixed global mesh, the wrapper replays the same seeded
samples and evaluates the same multivariate KDE and 0.61 contour on a 257-by-257
local adaptive mesh clipped to the MAPPINGS-V q/Z limits. Each integrated log
line records either `native-global-KDE` or `adaptive-local-KDE recovery`.
Central ratios outside the diagnostic grid remain invalid. This fallback is not
enabled for spatial bins and does not change the pyqz FITS maps.

`PYQZ_FLAG` is the raw package flag and `PYQZ_RS_OFFGRID` is the percentage of
Monte Carlo samples outside the grid. `NB_FLAG` is a bitmask: 1 is an O/H edge,
2 a log-U edge, 4 a log(P/k) edge, 8 an open/non-finite 68-percent interval, and
16 a fit exception. `JY22_FLAG` is a bitmask: 1 is an O/H posterior edge, 2 a
log-U posterior edge, 4 an invalid posterior, 8 an invalid/non-positive-definite
covariance, 16 an unexpected fit exception, and 32 empirical model-surface gross
mismatch (`chi2_min > 25`). The value 25 is a local MAUVE QC choice, not a
physically meaningful boundary or a literature-calibrated goodness-of-fit
threshold. For integer model maps, `-99`
means not evaluated; validity is 1 for a finite requested result and 0 for an
attempted invalid fit. In particular, `JY22_VALID=1` requires finite posterior
means, marginal bounds, minimum chi-square, and residuals; it does not assert
that the model fits. `JY22_FIT_OK=1` adds the empirical `chi2_min <= 25`
gross-mismatch screen; continuous chi-square and N2/S2/R3 residuals should be
retained and inspected rather than treating this binary cut as physical.
`JY22_FIT_OK=0` includes both numerically valid poor fits and attempted invalid
fits; inspect `JY22_VALID` and `JY22_FLAG` to distinguish them.
The three metallicity scales must remain separate:
offsets are expected because their grids, assumptions, and inference schemes
differ.

For development smoke tests, mirror one galaxy's gas, spatial-binning, and mask
inputs beneath a temporary `--root` (for example under `/private/tmp`) and run
the same command against that root. Do not use a live science-product directory
for schema testing because `SFR+Z.py` intentionally overwrites its derived
`*_further.fits` output.
