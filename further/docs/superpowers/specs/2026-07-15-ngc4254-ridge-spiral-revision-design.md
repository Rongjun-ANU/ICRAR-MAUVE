# NGC4254 Ridge-Based Spiral-Fit Revision Design

**Date:** 2026-07-15

**Target notebook:** `further/20260713_NGC4254_KTZ_spiral_fit.ipynb`

**Target map:** `further/v3tk_v7.6.8/NGC4254/NGC4254_gas_bin_maps_further.fits`, HDU `LOGSFR_SURFACE_DENSITY_HII`

## 1. Purpose

Revise the NGC4254 notebook so that observed spiral-arm geometry is detected from bright SFR ridges before the KTZ-compatible harmonic source profile is fitted. The present notebook instead selects geometry from a global full-field Fourier/source-rate objective. That objective is sensitive to lopsidedness, faint regions, and inter-arm structure and currently chooses a winding direction that does not follow the visually evident NGC4254 arms.

The revision must:

- preserve every valid HII pixel for morphological detection;
- search base arm multiplicities `m = 1, ..., 6` and both winding signs;
- select `m_arms`, signed pitch, and `Theta0` using a ridge-alignment statistic;
- fit `h_R`, `lambda0_0`, `eta`, and the KTZ harmonic coefficients only after geometry is selected;
- expose all harmonic-profile maxima rather than plotting only the global maximum;
- document the adopted centre, inclination, PA, finite-thickness correction, and coordinate conventions; and
- replace the current three-panel diagnostic with the approved 2-by-2 layout.

The notebook remains a low-dimensional SFR-morphology analysis. It will not infer diffusion, enrichment time, event clustering, pattern speed, or metallicity covariance.

## 2. Geometry provenance and inclination corrections

### 2.1 Adopted geometry

The notebook will obtain or explicitly validate the NGC4254 row in `Brown2021Table1.txt`:

- J2000 optical centre: `12h18m49.68s`, `+14d25m05.52s`;
- inclination: `39 deg`, where zero is face-on; and
- PA: `243 deg` east of north for the kinematically redshifted half.

The undirected major-axis line is equivalently `63 deg`. The directed value fixes the azimuth and `Theta0` convention. The FITS WCS reference coordinate is not the optical centre and must not replace it.

The notebook will print the provenance and adopted values beside the configuration. Missing or inconsistent catalog values are fatal rather than silently replaced by a fallback.

### 2.2 Finite-thickness surface-density correction

Following equation 3 of Huang et al. (2026), use

```text
q0 = 0.2
b_over_a = sqrt((1 - q0**2) * cos(inclination)**2 + q0**2)
```

For NGC4254 this gives `b_over_a approximately 0.787`. The upstream `SFR+Z.py` calculation already multiplies the projected SFR surface density by this factor before writing `LOGSFR_SURFACE_DENSITY_HII`. The notebook will calculate and report the factor as a provenance assertion but will not multiply the FITS values by it again.

This operation is distinct from deprojecting positions. For points in the disc mid-plane, the projected minor coordinate remains

```text
y_disc = y_projected / cos(inclination).
```

Replacing `cos(inclination)` with `b_over_a` in this positional transform would mix an apparent finite-thickness axial ratio with mid-plane foreshortening and is outside this design.

## 3. Data representations

Two representations serve different purposes:

1. **Valid HII pixels:** all finite HII pixels are retained for morphology and ridge detection because image area is the object being traced. No centre hole, outer quantile cut, or shared-completeness cut is allowed.
2. **Independent gas bins:** repeated pixel values are collapsed to one row per `BIN_ID`, retaining centroid and member-pixel area. These rows are used for parameter estimation and spatial resampling so repeated bin values do not create false statistical precision.

The notebook will state this distinction before either representation is used.

## 4. Coordinate transform

The FITS WCS converts pixel coordinates to sky coordinates; no manual array flip or rotation is performed. Tangent-plane offsets are positive east and north. With PA measured east of north,

```text
x_major       = east * sin(PA) + north * cos(PA)
y_projected   = -east * cos(PA) + north * sin(PA)
y_disc        = y_projected / cos(inclination)
R             = hypot(x_major, y_disc)
phi           = atan2(y_disc, x_major)
```

The inverse transform will continue to be verified by a sky-to-disc-to-pixel round trip. A synthetic orientation test must also verify that positive- and negative-pitch spirals retain their intended winding sign after projection and deprojection.

## 5. Axisymmetric background

Fit the robust area-weighted background

```text
Sigma_background(R) = lambda0_0 * exp(-R / h_R).
```

`h_R` is the e-folding length of the radial SFR source field, not an arm width, pitch, optical radius, or necessarily a stellar-disc scale length. The notebook will report both `h_R` and the equivalent logarithmic gradient `-1 / (ln(10) * h_R)` in dex per kpc.

For ridge detection first define the signed logarithmic residual

```text
Delta = log10(Sigma_SFR) - log10(Sigma_background).
```

Do not clip negative residuals: the local flanks on either side of an arm are
part of the ridge evidence.  Remove only broad azimuthal structure with a
mask-normalized circular smoother, giving a local ridge field

```text
C_local = Delta - broad_azimuthal_smooth(Delta).
```

The unclipped measurements also remain available for the later full
source-profile fit.

## 6. Ridge-based geometry detection

### 6.1 Log-polar representation

Map the deprojected HII pixels to `(u, phi)`, where `u = ln(R / R_ref)`. A logarithmic spiral becomes a straight periodic track in this plane. Use quantile edges in `u`, so every radial row contains approximately the same number of valid HII pixels and the very small central area is not discarded by a global-coverage threshold. Construct the log-polar residual map using mask-normalized accumulation and smoothing:

```text
smoothed_residual = smooth(weight * Delta) / smooth(weight).
```

The denominator prevents missing HII regions from being interpreted as zero
SFR.  Cells with no local effective coverage remain masked, but coverage is
tested against a small absolute floor rather than a fraction of the global
maximum.  This keeps all positive-radius HII pixels represented, including the
inner disc.

### 6.2 Candidate family

Search:

- `m_arms = 1, ..., 6`;
- both signed pitch ranges, initially `[-45, -5] deg` and `[5, 45] deg`;
- periodic phase `Theta0`; and
- the same radial range containing all positive-radius valid HII pixels.

The candidate phase is

```text
Theta = (m / tan(pitch)) * ln(R / R_ref) - m * phi + Theta0.
```

Use a narrow periodic ridge kernel around `Theta = 0` and subtract a broader,
unit-normalized kernel at the same location.  This difference-of-Gaussians
matched filter measures a local centre-versus-flanks ridge contrast rather than
brightness relative to a remote inter-arm track.  A configured physical width
`w` in kpc is converted separately at every radius,

```text
sigma_Theta(R) = m * w / (R * abs(sin(pitch))).
```

This radius dependence is required: using one phase width defined at `R_ref`
would make the physical corridor grow in proportion to radius and would score
broad lopsided structure as an arm.

### 6.3 Fair candidate score

The ridge score must compare the unit-normalized narrow arm kernel with its
broader local-flank control. It must be normalized by valid coverage, radial
support, and number of branches so that larger `m` cannot win merely by drawing
more curves. Candidate tables will report at least:

- `m_arms`;
- signed pitch;
- `Theta0`;
- narrow-kernel mean residual;
- broad local-flank mean residual;
- their difference-of-Gaussians response;
- valid coverage; and
- held-out-sector score.

Use spatial-sector holdouts to assess stability and require at least 10 of 12
sectors to return a finite test response. Geometry is selected by the median
held-out local-ridge response multiplied by phase stability and held-out
coverage.  Correct the remaining look-elsewhere bias across `m` with at least
eight row-wise azimuth-scrambled null maps: independently circular-shift each
log-radius row while shifting its mask and coverage with it, record the best
validated score for each `m`, and rank the real candidates by their
`m`-specific null z-score.  This destroys coherent spiral slope while preserving
the radial sampling, one-row residual distribution, and patchy mask.  Geometry
is not selected by the existing full-field intensity objective. Raw full-field
SSE and robust source-model loss may still be reported as secondary
diagnostics.

Adopt `w = 0.25 kpc` as the visible fiducial ridge scale and rerun the real-data
ranking over `w = 0.18, 0.22, 0.25, 0.30, 0.35 kpc`, repeating the same
per-`m` row-scrambled null calibration at every width.  Raw scores from
different widths are not directly comparable.  The fiducial result is the
reported geometry; the width scan is a sensitivity diagnostic, not an
opportunity to choose a scale after seeing the preferred winding.  If `m`, sign,
or pitch changes, report that multimodality explicitly.

### 6.4 Acceptance rule

The NGC4254 result must trace the visibly supported winding direction in both sky and deprojected views. The expected direction under the current disc-coordinate convention is the negative-pitch direction shown by the previous opposite-winding candidate. This expectation is a visual-validation target, not a hard-coded sign prior: both signs remain in the search.

If the highest ridge score still produces an overlay that does not follow the bright arm tracks, the notebook must report that no acceptable single global logarithmic-spiral geometry was recovered. It must not relabel a statistically preferred lopsided mode as an observed arm detection.

## 7. KTZ-compatible arm-profile fit

The harmonic series is inherited from `KTZ_validation.ipynb`:

```text
h(Theta) = sum_n g_n * cos(n * Theta + alpha_n).
```

In the KTZ validation notebook these values are prescribed inputs that shape a synthetic source-rate profile. KTZ does not provide the log-polar Fourier detection algorithm used in the current NGC4254 notebook. The revision will state this explicitly.

After ridge geometry fixes `m_arms`, pitch, and `Theta0`, fit the complete source field to the independent gas-bin catalog:

```text
lambda(R, phi) = lambda0_0 * exp(-R / h_R) * (1 + eta * h(Theta)).
```

Keep the identifiable convention `g_1 = 1` and `alpha_1 = 0`. Fit `lambda0_0`, `h_R`, `eta`, and the remaining `g_n` and `alpha_n` with positivity enforced over a dense phase grid. The full-field fit may not change the ridge-selected winding sign or base arm multiplicity.

`m_arms` will be described as the base periodic template. Harmonic orders can create secondary modulation within that period, so the notebook will also report the number and strength of actual local profile maxima.

## 8. No hidden harmonic maxima

Find every local maximum of the fitted periodic `h(Theta)` on a dense phase grid and refine its phase numerically. For every maximum report:

- phase;
- `h(Theta)`;
- total rate factor `1 + eta * h(Theta)`; and
- whether it is enhanced above the axisymmetric background.

Observed-map panels will show data-derived ridge skeletons, not harmonic maxima. The fitted-model panel will show every enhanced harmonic maximum. Any local maximum below the axisymmetric background will be shown in the phase-profile diagnostic but not called an arm; if plotted on the model map, it must use a thin neutral style and an explicit `below background` label.

## 9. Approved figures

### 9.1 Main 2-by-2 diagnostic

Use the approved Option A layout:

1. observed sky-plane `LOGSFR_SURFACE_DENSITY_HII` plus data-derived arm skeletons;
2. observed deprojected SFR map plus the same skeletons;
3. deprojected fitted KTZ-compatible source model plus all enhanced harmonic maxima; and
4. deprojected observed-minus-model residual plus data-derived skeletons.

Axes, winding sign, colours, and line meanings must be identical across relevant panels. A best alternative may be shown faintly only if it has equal branch treatment and a separate legend entry.

### 9.2 Supporting diagnostics

Add or retain:

- geometry and PA verification on the observed map;
- deprojected HII coverage;
- radial exponential background and `h_R` interpretation;
- ridge-score surfaces or curves for `m = 1, ..., 6` and both signs;
- a ranked candidate table with held-out scores;
- the harmonic arm profile with every maximum marked; and
- a stability summary from spatial-sector resampling.

## 10. Validation strategy

### 10.1 Automated contract tests

Tests must fail before implementation and then verify:

- `Q0 = 0.2` and `b_over_a approximately 0.787` for `i = 39 deg`;
- no second multiplication of the already corrected FITS SFR map;
- every valid HII pixel is retained for ridge detection;
- one independent record is retained per valid gas bin;
- `M_CANDIDATES` contains exactly `1, ..., 6`;
- both pitch signs are evaluated;
- candidate score normalization includes branch count and valid coverage;
- the selected geometry is based on held-out ridge score;
- the comparison across `m` is calibrated against row-scrambled null maps;
- the ridge phase width is converted from kpc separately at every radius;
- no fraction-of-global-maximum coverage cut deletes the inner disc;
- all harmonic local maxima are enumerated;
- the approved 2-by-2 panel contract is present; and
- the saved notebook has no execution errors, unexecuted code cells, or runtime warnings.

### 10.2 Synthetic recovery

Build masked synthetic SFR fields with known geometry and test:

- positive and negative winding;
- representative `m` values spanning 1 through 6;
- pitch and phase recovery within declared tolerances;
- no systematic preference for larger `m` on a null or single-arm field;
- robustness to patchy masks and radial decline; and
- correct sky projection and deprojection of the recovered skeleton.

### 10.3 Real-data acceptance

For NGC4254:

- inspect all four main panels at native saved resolution;
- confirm that the selected skeleton follows the bright SFR arm tracks in both sky and face-on views;
- confirm that the selected winding is not an artefact of axis reversal;
- compare the selected and best alternative ridge scores without unequal branch display; and
- retain a caution if sector holdouts or the declared ridge-width scan show
  multimodal `m`, pitch, or phase.

## 11. Files and scope

Planned implementation changes are limited to:

- `further/20260713_NGC4254_KTZ_spiral_fit.ipynb`; and
- `further/tests/test_20260713_ngc4254_ktz_spiral_fit.py`.

The target FITS product, `KTZ_validation.ipynb`, `SFR+Z.py`, and unrelated local theory or note files will not be modified. Existing user changes in the working tree must remain untouched and unstaged.

## 12. Known limitations

- NGC4254 is disturbed, lopsided, and may contain bifurcating arms that no single global pitch or exact `m` symmetry can describe.
- Patchy HII coverage can make an arm disappear over some radii.
- Ridge-score differences are descriptive stability measures, not formal posterior probabilities.
- The row-scrambled null z-score is an empirical morphology-ranking statistic
  from only eight scrambles, not a Gaussian detection significance or a
  substitute for a posterior model comparison.
- The SFR map is an adaptive-bin product and adjacent pixels are not independent measurements.
- The KTZ-compatible source model is not yet an observed metallicity-fluctuation or diffusion fit.
