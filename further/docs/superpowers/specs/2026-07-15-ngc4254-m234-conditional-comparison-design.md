# NGC4254 Conditional m=2,3,4 Spiral Comparison Design

**Date:** 2026-07-15

**Target notebook:** `further/20260713_NGC4254_KTZ_spiral_fit.ipynb`

**Target map:** `further/v3tk_v7.6.8/NGC4254/NGC4254_gas_bin_maps_further.fits`, HDU `LOGSFR_SURFACE_DENSITY_HII`

## 1. Purpose

Revise the NGC4254 notebook from a single arm-number selector into a
conditional comparison of three negative-winding spiral geometries:

- `m_arms = 2`;
- `m_arms = 3`; and
- `m_arms = 4`.

Each arm number is treated as its own scientific case.  The ridge analysis
independently estimates the best negative pitch angle and `Theta0` for that
fixed `m_arms`, and the KTZ-compatible source-profile fit is then performed at
that frozen geometry.  The notebook will present all three results without
declaring one arm number to be the measured truth.

This change is required because the previous single-winner implementation did
not pass scientific-quality review:

1. Eight row-scrambled null draws produced different winning arm numbers for
   different predeclared random streams.  The winners included `m=3`, `m=4`,
   and `m=6`; pooling 32 draws favored `m=4` rather than the committed `m=3`.
2. The broad 30-bin azimuthal preprocessing was performed before sector
   holdout.  A signal confined to the held sector leaked 49.26% of its
   post-filter L1 amplitude outside that sector and generated a false training
   ridge near the hidden phase.
3. Raw pre-splitting retained a negative-winding `m=3` family near
   `-30` to `-33 deg`, while `m=4` near the `-45 deg` grid boundary remained a
   plausible alternative.  A single exact geometry therefore should not be
   frozen before comparing the source-profile consequences.

The notebook remains a present-day SFR-morphology analysis.  It does not infer
metal diffusion, enrichment time, event clustering, pattern speed, or
metallicity covariance.

## 2. Decisions fixed by the user

The comparison has the following fixed interpretation:

- compare exactly `m = 2, 3, 4`;
- require negative winding for all three cases;
- optimize pitch and phase independently within each fixed arm number;
- fit and report a separate KTZ-compatible profile for each case; and
- do not use an exploratory null score to eliminate any of the three cases.

Constraining the sign makes the comparison test arm number rather than winding
direction.  It also follows the bright-arm direction identified visually in
both the observed and deprojected maps.  The notebook must describe this as an
explicit morphology prior, not as a sign-neutral statistical discovery.

## 3. Unchanged geometry and data contracts

The following validated choices remain unchanged:

- Brown et al. optical centre: `12h18m49.68s +14d25m05.52s`, FK5 J2000;
- inclination: `39 deg`;
- directed receding-side PA: `243 deg` east of north;
- finite-thickness parameter: `q0 = 0.2`;
- equation-3 surface-density factor: `b/a approximately 0.787`, already
  applied upstream and not multiplied into the FITS map again;
- positional minor-axis deprojection: divide by `cos(inclination)`, not `b/a`;
- the FITS WCS handles sky orientation, with no manual north-up/east-left
  array flip; and
- all 681,856 finite HII pixels remain available to morphology detection,
  while the 32,552 unique gas bins remain the independent profile-fit records.

No artificial centre hole or outer radial cut is permitted.  The positive
radius guard exists only because logarithmic-spiral phase is undefined at
exactly zero radius; the live minimum pixel radius is about `0.00664 kpc`.

## 4. Leakage-free ridge validation

### 4.1 Raw log-polar state

The log-polar builder will retain the unsmoothed arrays needed for fold-safe
processing:

- raw weighted residual sum;
- raw coverage;
- fixed full-sample `u` and `phi` edges;
- fixed positional `u` centres; and
- the pixel-to-row and pixel-to-azimuth mapping needed for diagnostics.

The signed logarithmic residual remains

```text
Delta = log10(Sigma_SFR) - log10(Sigma_background).
```

Negative residuals are not clipped.  Missing coverage remains invalid and is
never treated as a zero-valued SFR observation.

### 4.2 Split before smoothing

For each of 12 held-out azimuth sectors, construct training and test fields by
masking the raw weighted sum and raw coverage before Gaussian smoothing.  Both
fields use the same fixed `u` and `phi` edges.  Numerator and denominator are
smoothed identically with `mode=("nearest", "wrap")`, and division occurs only
where local effective coverage exceeds the declared absolute floor.

A four-bin circular guard is removed from both sides of the sector boundary.
This is slightly wider than three times the 1.2-bin azimuthal smoothing sigma.
The guard prevents locally smoothed test values from entering the training
field and makes the train/test supports disjoint at the resolution being
scored.

The held-out selector will score the mask-normalized `radial_residual` field.
It will not apply the 30-bin broad-azimuth filter before sector splitting.  The
later phase-space narrow-minus-broad response already supplies the required
local centre-versus-flanks contrast.  A 30-bin-sigma preprocessing kernel
cannot be validated independently inside a test sector that is itself about
30 bins wide.

The global axisymmetric exponential remains a shared nuisance fit.  Therefore
the notebook will call the result a screened, leakage-controlled sector
robustness score rather than a formally unbiased cross-validation statistic.

### 4.3 Behavioral leakage test

Add a reduced adversarial test in which all injected ridge signal lies inside
one held sector.  The fold-safe training field must contain zero raw signal,
must not recover the hidden phase, and must return no positive validated ridge
score attributable to that signal.  The test must fail if splitting is moved
back after broad full-map preprocessing.

## 5. Arm-number-conditional ridge search

Declare the comparison explicitly:

```text
M_COMPARE = [2, 3, 4]
pitch grid = -45, -44, ..., -5 deg
```

For each fixed `m_arms`:

1. screen every negative pitch using the full leakage-controlled radial
   residual field;
2. retain the configured shortlist of pitches for sector validation;
3. determine phase on the training field and evaluate it only in the disjoint
   held-sector field;
4. require at least 10 of 12 finite held-sector results;
5. calculate median held-sector response, held-sector dispersion, valid-sector
   fraction, and circular phase concentration; and
6. choose the pitch and `Theta0` with the highest positive validated score
   within that fixed arm number.

The physical ridge-width conversion remains

```text
sigma_Theta(R) = m * width_kpc / (R * abs(sin(pitch))).
```

The notebook must report how many radial bands hit the configured minimum or
maximum phase-width clamp.  A pitch at `-45 deg` is a grid-boundary result and
must receive a boundary flag rather than being presented as a well-contained
estimate.

The three selected dictionaries are stored in arm-number order and retain a
clear geometry source label such as
`leakage_controlled_negative_winding_ridge`.  No later source-profile optimizer
may change `m_arms`, pitch, or `Theta0`.

## 6. Null calibration is descriptive

At the fiducial physical ridge width of `0.25 kpc`, use 32 row-scrambled null
draws organized as four predeclared sequential eight-draw blocks with stream
seeds:

```text
4254, 5254, 6254, 7254
```

Every draw circularly shifts residual, coverage, and mask together within each
radial row.  The shift destroys coherent cross-radius spiral slope while
preserving each row's value distribution and sampling pattern.

For each fixed `m=2,3,4`, calculate the best negative-pitch null score and
report:

- pooled 32-draw mean and sample standard deviation;
- descriptive `null_z` for the real conditional geometry;
- the same values in each eight-draw block;
- blockwise pitch rank or candidate rank; and
- whether the descriptive ordering changes across blocks.

The null statistic is the finite leakage-controlled `validated_score` for
every returned fixed-`m` row.  Ridge acceptance remains separate metadata:
below-threshold null rows retain their finite score and must not be replaced
with zero.  A non-finite score is an execution error rather than a synthetic
zero draw.  This keeps the real and null statistics on the same scale.

These null values diagnose how unusual each conditional ridge score is under
row scrambling.  They do not select among `m=2,3,4`, are not Gaussian
detection significances, and are not posterior probabilities.  The notebook
must retain the observed seed sensitivity in its interpretation.

To control memory, the full candidate screen will store scalar rows only.
Sector histograms are recomputed one shortlisted candidate at a time.  The
validated scratch implementation reproduced all 492 prior candidate rows
exactly while reducing histogram payload from 3.8006 GiB to 7.91 MiB.

## 7. Physical-width sensitivity

Repeat the leakage-controlled data search for each fixed arm number at

```text
0.18, 0.22, 0.25, 0.30, 0.35 kpc.
```

Report the selected negative pitch, `Theta0`, validated score, phase
concentration, held-sector count, and pitch-boundary flag for every
`(m_arms, width)` pair.  Raw scores at different widths are not directly
compared because changing the physical kernel changes the score scale.

The 32-draw pooled null calibration is required at the fiducial width only.
Width sensitivity is a conditional morphology diagnostic, not another route
to choose an arm number or an optimized ridge width after inspecting results.

## 8. Three fixed-geometry KTZ-compatible fits

For each of the three ridge geometries, fit the independent gas-bin catalog to

```text
lambda(R, phi) = lambda0_0 * exp(-R / h_R)
                 * (1 + eta * h(Theta)),

h(Theta) = sum_n g_n * cos(n * Theta + alpha_n).
```

The optimizer may fit only:

- `lambda0_0`;
- `h_R`;
- `eta`;
- the non-reference harmonic amplitudes `g_n`; and
- the non-reference phases `alpha_n`.

It may not fit or refine `m_arms`, pitch, or `Theta0`.  Preserve the
identifiability convention `g_1 = 1` and `alpha_1 = 0`, and enforce a positive
source-rate factor on a dense phase grid.

For every case, enumerate and refine every local maximum of the fitted
harmonic profile.  Record each maximum's phase, harmonic value, total rate
factor, and whether it lies above the axisymmetric background.  Observed-map
skeletons remain data-derived `Theta=0` ridges; source-model maxima are shown
only on model/profile diagnostics and are not relabelled as detected arms.

## 9. Figures

### 9.1 Main 2-by-2 comparison

Create one 2-by-2 comparison figure:

1. observed sky-plane SFR map with the `m=2,3,4` data-derived skeletons in
   distinct, consistent colors;
2. observed deprojected SFR map with the same three skeleton sets;
3. the three fitted harmonic arm profiles with all local maxima marked; and
4. a compact comparison of leakage-controlled score, descriptive pooled-null
   result, and ridge-width sensitivity.

The legend must state that all cases use negative winding and that the three
colors identify conditional arm-number hypotheses, not confidence levels.

### 9.2 Detailed model and residual comparison

Create a separate detailed figure with one row per arm number and consistent
columns for:

- the deprojected fitted KTZ-compatible source model; and
- observed-minus-model residual in dex.

Each model panel shows all enhanced source-profile maxima for that conditional
fit.  Each residual panel shows the data-derived skeleton for the same fixed
arm number.  Shared color limits are used where scientifically meaningful so
visual differences are not manufactured by rescaling each row.

Retain supporting figures for geometry, deprojected HII coverage, radial
background, log-polar residuals, synthetic validation, null-block stability,
and sector/width sensitivity.

## 10. Comparative parameter table

The final table contains one row per `m_arms` and no highlighted winner.  It
must include at least:

- fixed `m_arms`;
- negative ridge-selected pitch and `Theta0`;
- geometry-source label;
- validated ridge score and held-sector count;
- phase concentration;
- pooled descriptive null mean, standard deviation, and `null_z`;
- eight-draw block variability summary;
- pitch-boundary flag;
- fitted `lambda0_0`, `h_R`, and gradient in dex/kpc;
- fitted `eta`;
- `harmonic_n`, `harmonic_g`, and `harmonic_alpha`;
- minimum source-rate factor;
- number of enhanced harmonic maxima; and
- weighted log-space SSE or robust cost.

Rows are ordered `m=2`, `m=3`, `m=4`.  The surrounding Markdown must say that
all source-profile parameters are conditional on the chosen arm-number
hypothesis.

## 11. Validation strategy

### 11.1 Source and behavioral contracts

Tests must verify:

- the comparison set is exactly `[2, 3, 4]`;
- the pitch grid is negative only;
- raw accumulators are split before smoothing;
- the four-bin boundary guard is present;
- the broad 30-bin `local_ridge` field does not enter held-out scoring;
- the adversarial held-sector signal gives zero training response;
- the full candidate screen does not retain the 3.80-GiB histogram cache;
- four deterministic null blocks provide exactly 32 draws per conditional
  arm-number family;
- null values are not used to remove an arm-number case;
- all three KTZ fits preserve their supplied geometry exactly;
- all local harmonic maxima are enumerated; and
- the notebook declares the comparative figures and three-row final table.

### 11.2 Numerical checks

Disposable execution must establish:

- all 681,856 finite HII pixels remain represented;
- reduced synthetic negative-winding cases, each run with `m` supplied as a
  one-element comparison set, preserve that fixed `m` and recover injected
  pitch and phase within declared tolerances; this is not an arm-number
  discrimination test;
- the adversarial leakage test passes;
- every real conditional case has at least 10 finite held sectors and a
  positive validated score;
- the null table has exactly `3 * 32` conditional summary draws at the
  fiducial width, plus block identifiers;
- the width table has exactly `3 * 5` rows;
- all three fixed-geometry KTZ fits succeed with positive rate factors;
- every geometry supplied to the KTZ fitter equals the corresponding geometry
  returned unchanged in its fitted dictionary;
- figures contain the required sky, deprojected, profile, model, and residual
  panels; and
- no `RuntimeWarning` is saved in notebook output.

### 11.3 Visual acceptance

Inspect every generated PNG at original resolution.  For each `m` case:

- the negative-winding skeleton must follow the intended winding convention;
- sky and deprojected skeletons must be consistent under the round-trip
  transform;
- crossing portions must not be described as observed arms merely because the
  global logarithmic template passes through them;
- all skeleton branches must be visible; and
- every harmonic maximum shown on a model/profile must match the enumerated
  maxima table.

## 12. Interpretation limits

The notebook must conclude that:

- `m=2`, `m=3`, and `m=4` are conditional morphology hypotheses;
- negative winding is imposed from visual morphology to isolate the arm-number
  comparison;
- a single global logarithmic template is an approximation to NGC4254's
  asymmetric, environmentally perturbed arms;
- pitch and phase may vary with physical ridge width and held sector;
- pooled null values remain descriptive despite using 32 rather than eight
  draws;
- `h_R` describes the axisymmetric SFR source field, not arm width or
  automatically the stellar-disc scale length; and
- KTZ-compatible profile parameters describe present-day SFR morphology, not
  the full Krumholz--Ting stochastic metallicity model.

## 13. Acceptance criteria

The revision is accepted only when:

1. leakage-controlled negative-winding geometries exist for all `m=2,3,4`;
2. no exploratory score removes one of those three cases;
3. three fixed-geometry KTZ fits complete without changing geometry;
4. all valid HII pixels remain represented and no centre hole is introduced;
5. the main 2-by-2 and detailed model/residual figures show all three cases;
6. the parameter table has exactly three conditional rows;
7. source, numerical, and saved-result contracts pass;
8. every saved figure passes visual inspection; and
9. the final interpretation reports instability and conditionality without a
   preferred arm-number claim.
