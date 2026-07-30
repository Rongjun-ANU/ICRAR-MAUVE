# NGC4254 KTZ Spiral-Geometry Fit Design

**Date:** 2026-07-13

**Output notebook:** `20260713_NGC4254_KTZ_spiral_fit.ipynb`

**Reference notebook:** `KTZ_validation.ipynb`

**Input product:** `v3tk_v7.6.8/NGC4254/NGC4254_gas_bin_maps_further.fits`

**Target HDU:** `LOGSFR_SURFACE_DENSITY_HII`

## Goal

Build a self-contained, tutorial-style notebook that measures a KTZ-compatible logarithmic-spiral source-rate model from the NGC4254 HII-region star-formation-rate surface-density map. The notebook must recover and report the same geometry and arm-profile parameter family used by the first two code cells of `KTZ_validation.ipynb`: `m_arms`, `pitch_angle`, `Theta0`, `h_R`, `eta`, and the harmonic indices, amplitudes, and phases. It must also reproduce the two principal KTZ geometry visualizations: a spiral skeleton and an arm profile.

The notebook is a morphology fit to the SFR source field. It is not yet a metallicity-correlation fit, a KTZ diffusion fit, or a claim that NGC4254 is described perfectly by a stationary global density wave.

## Scientific Model

The observed logarithmic map is converted to a linear surface density before fitting:

$$
\Sigma_{\rm SFR}=10^{\log_{10}\Sigma_{\rm SFR}}.
$$

The KTZ-compatible model is

$$
\lambda(R,\phi)=\lambda_{00}\exp(-R/h_R)\left[1+\eta h(\Theta)\right],
$$

with logarithmic-spiral phase

$$
\Theta(R,\phi)=\frac{m}{\tan p}\ln\left(\frac{R}{R_{\rm ref}}\right)-m\phi+\Theta_0,
$$

and harmonic arm profile

$$
h(\Theta)=\sum_{n=1}^{N_h}g_n\cos(n\Theta+\alpha_n).
$$

Here `m` is an integer arm multiplicity, `p` is the pitch angle, `Theta0` is the phase at `R_ref`, `h_R` is the exponential radial scale, and `eta` controls the fractional non-axisymmetric modulation. The model will use `R_ref = 1 kpc`, matching the reference notebook's interpretation of `R_in = 1 kpc` as the phase reference.

The products `eta * g_n` are identifiable, but `eta` and an arbitrary common scaling of all `g_n` are not. The notebook will therefore impose `g_1 = 1` and `alpha_1 = 0`; `Theta0` absorbs the fundamental phase. Higher-order `g_n` and `alpha_n` remain free. This normalization makes `eta` and the harmonic coefficients reproducible.

## Geometry and Coordinates

The fiducial deprojection will use the current local MAUVE convention already documented for this product:

- distance: 16.5 Mpc;
- centre: RA 12:18:49.68, Dec +14:25:05.52;
- inclination: 39 degrees;
- position angle: 243 degrees, equivalent to 63 degrees modulo 180.

The notebook will state the position-angle convention explicitly and implement the sky-to-disc rotation in one named helper. It will show a coordinate sanity-check plot with the centre, projected major axis, and projected minor axis on the observed map. This check is required because a PA sign or axis-order error directly changes the fitted spiral handedness and pitch.

The WCS in the target HDU supplies the sky coordinates. Angular offsets will be converted to physical offsets using the fixed 16.5 Mpc distance, then rotated and stretched by `1/cos(inclination)` along the projected minor axis. The notebook will report both sky-plane and disc-plane coordinate ranges.

## Sampling and Masking

`LOGSFR_SURFACE_DENSITY_HII` contains repeated pixels from adaptive gas bins. Treating every repeated pixel as an independent datum would create false precision. The notebook will combine the SFR map with `BIN_ID` and create one record per valid gas bin containing:

- bin ID;
- median finite log SFR value;
- sky-coordinate centroid;
- deprojected `x`, `y`, `R`, and `phi`;
- number of member pixels, used as the represented area;
- finite/valid status.

The fit will be area weighted so the centroid representation approximates the original spatial integral without pretending that member pixels are independent measurements. Invalid values will remain masked. The notebook will not fill HII-mask holes or smooth across them.

A configurable radial fitting interval will exclude the unresolved centre and poorly sampled outer edge. Its default will be derived transparently from radial coverage diagnostics rather than hidden inside an optimizer. The chosen limits and the fraction of retained bins and area will be printed.

## Fitting Strategy

### 1. Axisymmetric radial component

Fit `lambda00` and `h_R` to the azimuthally summarized linear SFR field using a robust regression in log space. The radial diagnostic will show individual bin centroids, radial summary points, the fitted exponential, and the retained fitting interval.

Define the normalized residual field

$$
q_i=\frac{\Sigma_{{\rm SFR},i}}{\lambda_{00}\exp(-R_i/h_R)}-1.
$$

This isolates the angular modulation used for the spiral search.

### 2. Discrete arm-number and pitch search

For each integer `m` from 1 through 4 and each pitch on a declared grid, calculate the logarithmic-spiral phase without `Theta0` and measure the weighted complex Fourier coefficient of the residual field. This produces a two-dimensional spiral-power surface over `(m, pitch)` and supplies robust starting values for pitch and phase.

The search must evaluate both winding signs or an equivalent signed-pitch convention. The notebook will report the convention used in the final table and show the leading alternatives, not only the winning value.

### 3. Harmonic arm-profile fit

Fold the normalized residuals by the candidate phase and fit the fundamental plus higher harmonics. The default maximum harmonic order will be three, matching `KTZ_validation.ipynb`. With `g_1 = 1` and `alpha_1 = 0`, the fit returns `eta`, `g_2`, `g_3`, `alpha_2`, and `alpha_3` together with the global phase `Theta0`.

The positivity condition

$$
1+\eta\min_\Theta h(\Theta)>0
$$

will be enforced or checked explicitly so the inferred source-rate field cannot become negative.

### 4. Constrained joint refinement

Starting from the staged solution, jointly refine the continuous parameters with bounded robust least squares while holding integer `m_arms` fixed. The objective will compare the linear SFR field with the full model and retain area weights. Bounds will prevent singular pitch angles, negative radial scales, and negative model intensities.

The notebook will compare the staged and refined objectives. If refinement lands on a bound or violates the positivity constraint, it will retain and clearly label the staged fit rather than silently reporting an invalid solution.

### 5. Uncertainty and model comparison

Uncertainties will come from reproducible bootstrap refits with a fixed random seed. Resampling will operate on spatial sectors or blocks rather than individual neighbouring bins to reduce the effect of spatial dependence. The notebook will report percentile intervals for all continuous parameters and the bootstrap selection frequency of each `m_arms` value.

The discrete candidates `m = 1, 2, 3, 4` will be compared using the same weighted objective plus an information criterion. Because NGC4254 is disturbed and asymmetric, the notebook will distinguish the best global approximation from proof of a literal arm count.

## Notebook Structure

Every code cell will be preceded by a Markdown cell explaining:

1. the scientific purpose;
2. the mathematical operation;
3. the inputs and outputs;
4. the assumptions and potential failure modes;
5. how to interpret the following plot or printed values.

Code cells will also contain step-by-step inline comments. The planned sequence is:

1. title, scope, and interpretation limits;
2. mathematical model and parameter-identifiability explanation;
3. imports, paths, geometry constants, fitting controls, and random seed;
4. FITS loading and HDU/WCS validation;
5. unique-bin table construction;
6. sky-plane map and geometry sanity check;
7. disc deprojection and radial-coverage diagnostics;
8. exponential radial-profile fit;
9. spiral Fourier-power search across `m` and pitch;
10. harmonic arm-profile fit;
11. constrained full-model refinement;
12. bootstrap uncertainty and discrete-model comparison;
13. final KTZ-compatible parameter table;
14. observed-map skeleton overlay and residual map;
15. KTZ-style Cartesian skeleton and fitted arm-profile plot;
16. conclusions, caveats, and how these parameters feed a later KTZ metallicity analysis.

## Required Outputs

The notebook must display, without requiring external saved figures:

- a compact input-data and mask summary;
- the observed HII log-SFR map with centre and orientation axes;
- the deprojected map or centroid representation;
- radial SFR profile and fitted `h_R`;
- spiral-power heatmap or curves for `m = 1..4`;
- ranked arm-number/pitch candidates;
- phase-folded data and fitted harmonic arm profile;
- observed SFR map with the fitted spiral skeleton overlaid;
- a model map and residual map;
- a KTZ-style Cartesian skeleton plot in kpc;
- the KTZ-style `h(Theta)` plot with individual harmonic components;
- a final parameter table containing estimates, uncertainty intervals, units, conventions, and the corresponding reference-notebook variable names.

The final table will include at least:

| Reference variable | NGC4254 result |
|---|---|
| `R_in`, `R_out` | fitted radial domain |
| `m_arms` | selected integer and alternative scores |
| `pitch_angle` | signed degrees |
| `Theta0` | radians modulo `2 pi` |
| `lambda0_0` | fitted linear SFR normalization |
| `h_R` | kpc |
| `eta` | dimensionless modulation |
| `harmonic_n` | `[1, 2, 3]` by default |
| `harmonic_g` | normalized with `g_1 = 1` |
| `harmonic_alpha` | radians with `alpha_1 = 0` |

Diffusion and clustering parameters from the first reference cell (`kappa`, `x0`, `t_star`, and Thomas-process settings) will not be fitted from the SFR morphology map. The notebook will list them separately as out of scope rather than implying that morphology determines them.

## Error Handling and Reproducibility

The notebook will fail early with explicit messages if the FITS file, target HDU, `BIN_ID`, celestial WCS, or required geometry is missing. It will assert matching map shapes and sufficient finite bins. All optimizer bounds, radial cuts, harmonic order, random seeds, and bootstrap counts will be visible in the configuration cell.

The primary notebook should run with `/opt/miniconda3/envs/ICRAR/bin/python` and its `ICRAR` Jupyter kernel. No science product will be overwritten. Plot saving will be disabled by default; displayed outputs will be stored in the executed notebook.

## Verification

Verification will be proportional to the scientific risk:

1. validate notebook JSON and cell ordering;
2. execute the notebook against the real local NGC4254 FITS product;
3. confirm that every code cell completes and every required plot/table is produced;
4. check that all fitted intensities are finite and non-negative;
5. check that the recovered skeleton is geometrically aligned with the high-SFR structure by visual inspection;
6. run a synthetic injection-recovery check using a known KTZ spiral sampled through the real NGC4254 bin geometry;
7. verify that the recovered pitch, phase, radial scale, modulation, and harmonic coefficients agree with injected values within declared tolerances;
8. report convergence warnings, model inadequacy, or bootstrap instability instead of suppressing them.

## Interpretation Limits

The fitted result is the best low-dimensional, global KTZ-compatible summary of the HII SFR morphology under the fixed geometry. NGC4254 is lopsided and disturbed, so one global logarithmic spiral may leave coherent residuals and may not correspond to a unique physical density wave. The HII mask is selection dependent, and the adaptive-bin measurements are spatially correlated.

The fit can provide a source-field template for a later metallicity-correlation forward model. It cannot by itself determine turbulent diffusion parameters, enrichment ages, source clustering, pattern speed, or the dynamical lifetime of the arms.
