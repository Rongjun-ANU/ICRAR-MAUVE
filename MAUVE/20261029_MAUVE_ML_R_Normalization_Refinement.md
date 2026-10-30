# Refining the nGIST LIGHT-to-R Stellar Mass-to-Light Conversion

**MAUVE technical report | 29 October 2026 | Rongjun Huang**

**Scope:** a consolidated explanation and implementation specification based on this conversation, inspected local code and products, nGIST GitHub resources, and the executed Setonix LSF test. This report does not implement the refinement or present newly recalculated galaxy masses.

## 1. Main conclusion and status

The current `light_norm_to_r` branch has the appropriate conversion structure: translate fitted template-normalization light weights into R-band luminosity fractions, then average the SSP R-band mass-to-light ratios. The part to refine is the per-SSP normalization factor. Current `Mass.py` measures the mean original template flux inside `LMIN`-`LMAX`; nGIST measures the mean after padding/cropping, resolution convolution, and logarithmic rebinning.

Reconstructing this factor does not require rerunning the galaxy fits. It requires repeating the deterministic template preparation with the same library, run parameters, LSF curves, and numerical routines. The factors do not need to have been stored in the products.

The previously suspected LSF-file mismatch is resolved for the checked Setonix production directory: `lsf_MUSE-WFM-MILES` is byte-for-byte identical to the upstream example named `lsf_MUSE-WFM`. The actual nGIST 7.6.8 loader accepted arbitrary filenames with identical contents and raised `FileNotFoundError` for a nonexistent configured file.

**Current state:** original SSP spectra and BaSTI photometry already exist in `further/`; run configurations and fitted weights are available; the relevant local 7.6.4 functions match the checked upstream 7.6.8 functions; the production LSF files have been located and tested remotely. They have not yet been copied into `further/` during this task. The refined normalization and resulting M/L differences remain uncomputed.

## 2. Notation and the three wavelength definitions

| Symbol | Meaning |
|---|---|
| b | Voronoi-bin index |
| j | SSP-template index |
| T_j(lambda) | Original library SSP spectrum on its original population scale |
| P_j(lambda) | SSP after nGIST crop, convolution, and logarithmic rebinning, before LIGHT normalization |
| n_j | Arithmetic mean of P_j over the prepared logarithmic wavelength samples |
| t_j(lambda) | Mean-one fitting template, P_j/n_j |
| w_bj | Saved normalized nGIST coefficient for SSP j in bin b |
| ell_R,j | Original SSP R-band luminosity, represented by a relative photometric proxy |
| Upsilon_R,j | Tabulated original SSP R-band stellar M/L |
| f_R,bj | Converted, normalized R-band light fraction |
| Upsilon_R,b | Final composite R-band M/L for a bin |

Three wavelength definitions must remain distinct. **Original SSP coverage** is everything in the template FITS file. **The fitting interval** is the galaxy-spectrum region used to constrain the mixture. **The R passband** is a wavelength-dependent filter response used to define luminosity. A template can extend well beyond both the fit and the R passband. It is not an R-only spectrum.

For the checked local E-MILES file, the original grid spans 1680.2-49999.4 Angstrom with 0.9 Angstrom sampling. This is a representative header inspection, not an assertion that every possible MILES/E-MILES distribution has that coverage. The library's broad coverage does not mean all those wavelengths participate in the normalization: nGIST first crops the spectrum.

## 3. Why normalize a template to mean one?

Original SSP spectra have a model population normalization, rather than a requirement that their spectral means equal one. On a common population scale, young and old SSPs can have different luminosities. Cropping or changing the sampling can also change a spectrum's mean even if it had previously been normalized over some other interval.

In LIGHT mode, nGIST defines:

$$
n_j = \frac{1}{N_{\log}}\sum_{i=1}^{N_{\log}}P_j(\lambda_i), \qquad t_j(\lambda_i)=\frac{P_j(\lambda_i)}{n_j}.
$$

Therefore `mean(t_j) = 1`. Dividing all spectral samples by one scalar preserves the template's relative spectral shape, including absorption features and continuum slope. It gives every template the same reference mean brightness, not the same spectrum.

For example, prepared samples `[2, 4, 6]` have mean 4. Dividing by 4 gives `[0.5, 1, 1.5]`, whose mean is one. A second SSP is divided by its own mean, which need not be 4.

This normalization defines a useful coefficient convention. In the simplified stellar mixture `G = sum_j b_j t_j`, each coefficient multiplies a mean-one template. Once normalized across SSPs, these coefficients describe contributions under that template-light convention. They are not automatically luminosity fractions in Johnson V or R. Actual fits also include kinematics, attenuation, multiplicative polynomials and masking; this simplified expression isolates amplitude bookkeeping, not the entire fitting model.

Normalization is not mathematically required to represent a spectrum. Different template scaling implies different coefficient scaling; an equivalent unregularized model can compensate. Regularization and other fitting choices must be treated consistently when comparing alternative fitting parameterizations. This proposal changes post-processing of an existing fit, not the original fitting parameterization.

## 4. Why the conversion contains w times ell_R divided by n

Substitute the mean-one template into its fitted contribution:

$$
b_j t_j(\lambda)=b_j\frac{P_j(\lambda)}{n_j}=\frac{b_j}{n_j}P_j(\lambda).
$$

The amplitude multiplying the unnormalized prepared SSP is therefore `b_j/n_j`. Since its original-population R luminosity is `ell_R,j`, its R-band contribution is proportional to `b_j ell_R,j/n_j`.

Recovering the original spectrum would multiply the normalized spectrum by `n_j`. Keeping the fitted spectral contribution unchanged then requires dividing its coefficient by `n_j`. This is why multiplying the coefficient by `n_j` would go in the wrong direction.

**Illustrative example, not a measurement:** let `n_j = 10`, original R luminosity `ell_R,j = 30`, and fitted coefficient `b_j = 0.4`. The normalized template's corresponding luminosity is `30/10 = 3`, so the fitted luminosity contribution is `0.4 x 3 = 1.2`. Equivalently `(0.4/10) x 30 = 1.2`. The alternative `0.4 x 30 x 10 = 120` does not describe the component fitted.

Saved nGIST weights have a common per-bin scaling removed when their sum is normalized. That omitted scale cancels in both R-band fractions and composite M/L; it need not be recovered for this calculation.

## 5. Final equations and what sums to one

The table supplies an R magnitude and an R-band M/L for each matched SSP:

$$
\ell_{R,j}\propto10^{-0.4R_j}, \qquad \Upsilon_{R,j}=({\rm ML\_R})_j.
$$

The proportionality constant from the magnitude zero point is common to all SSPs and cancels in the ratio. This is a model luminosity proxy, not the observed galaxy flux.

Define the unnormalized R contributions and their fractions:

$$
q_{R,bj}=w_{bj}\frac{\ell_{R,j}}{n_j}, \qquad f_{R,bj}=\frac{q_{R,bj}}{\sum_k q_{R,bk}}.
$$

$$
\sum_j f_{R,bj}=1, \qquad \Upsilon_{R,b}=\sum_j f_{R,bj}\Upsilon_{R,j}.
$$

Equivalently, calculate the final ratio directly:

$$
\boxed{\Upsilon_{R,b}=\frac{\sum_j w_{bj}(\ell_{R,j}/n_j)\Upsilon_{R,j}}{\sum_j w_{bj}(\ell_{R,j}/n_j)}}.
$$

Both the original normalized nGIST weights and the final R fractions sum to one, but describe different light definitions. The intermediate `q_R` values generally do not sum to one. Finite nonzero denominators are required. With nonnegative weights, the result is a convex average of the contributing SSP M/L values.

For numerical work, explicitly validate finite nonnegative input coefficients and positive finite normalization/luminosity terms. Do not silently clip unexpected negative weights or renormalize a bin with a zero denominator.

## 6. What the present scripts do

Both current files inspected in this conversation select `sfh_weight_interpretation = "light_norm_to_r"`: `further/Mass.py` and `extended/Mass.py`. The initial suggestion that the current `extended` script selects `light_r` was corrected after source inspection. Both also provide a `light_v_to_r` fallback.

| Interpretation | Input-weight assumption | Composite M/L calculation |
|---|---|---|
| light_r | Already R-band luminosity fractions | Arithmetic average of SSP R M/L |
| light_norm_to_r | Mean prepared-template light convention | Convert with ell_R/n, normalize, then average |
| light_v_to_r | Approximate V-light convention | Reweight using 10^(0.4(V-R)); approximate for generic fitted-window weights |
| mass | Mass fractions | Harmonic mean: 1 / sum_j(mu_j/Upsilon_R,j) |

The current norm-to-R helper instead estimates `n_j` by `nanmean` of the original linear-grid spectrum between the weight-header `LMIN` and `LMAX`. It omits the preparation padding, LSF convolution, and logarithmic resampling. The final conversion formula is already structured correctly; its normalization-factor input needs refinement.

For a production mode labelled as exact reproduction, failure to reconstruct the factors should be an explicit error or a clearly labelled approximate mode. An automatic approximate fallback must not retain an exact-reconstruction label.

## 7. Required inputs and their provenance

| Quantity | Source | Role |
|---|---|---|
| w_bj | GALID_sfh_weights.fits, WEIGHTS extension, WEIGHTS column | Saved fitted mixture |
| SSP identities | Same file, GRID extension: LOGAGE, METAL, ALPHA | Match all arrays to the same SSP order |
| LMIN, LMAX | Weight-file SFH headers; cross-check CONFIG.SFH | Requested template-preparation interval |
| NORM_TEMP, TEMPLATE_SET, LIBRARY, LSF_TEMP | Weight headers and CONFIG.SFH | Verify preparation mode and library |
| Original spectra | GENERAL.TEMPLATE_DIR + SFH.LIBRARY | Input spectra for reconstructing n_j |
| Redshift | CONFIG.GENERAL.REDSHIFT | Reproduce data-LSF preparation |
| Data LSF | GENERAL.CONFIG_DIR + GENERAL.LSF_DATA | Instrumental FWHM curve |
| Template LSF | GENERAL.CONFIG_DIR + SFH.LSF_TEMP | SSP FWHM curve |
| Galaxy velocity step | Production bin_spectra.hdf5 attribute VELSCALE; standard preparation writes CONFIG.PREPARE_SPECTRA.VELSCALE | Sets logarithmic sampling |
| Template sampling ratio | SFH wrapper, velscale_ratio = 2 | Template step = galaxy step / 2 |
| R magnitude and SSP M/L | BaSTI+Chabrier.dat, R and ML_R | Photometric conversion |
| Numerical behavior | nGIST preparation code and pPXF convolution/rebin utilities | Reproduce original operations |

Use the run's own `CONFIG`, rather than a later master configuration. Relative production paths describe that environment; map them to verified local resources explicitly. If weight-header settings disagree with `CONFIG`, report the conflict rather than choosing silently. For a 7000 product, use that product's own limits, not its directory name.

The user confirms the SSP spectra and photometry are correctly matched; this report takes that as an established input condition and does not reopen the IMF/model matching question.

### 7.1 Checked IC3392 example

The inspected weights array has shape `(4062, 477)`: 4062 bins and 477 SSPs. The grid has `LOGAGE`, `METAL`, `ALPHA`, with age in Gyr recovered as `10**LOGAGE`. The local library contains 477 FITS templates; count agreement is useful but does not replace identity matching.

The product CONFIG gives `GENERAL.REDSHIFT = 0.005566`, `PREPARE_SPECTRA.VELSCALE = 41`, `SFH.LMIN = 4800`, `SFH.LMAX = 8900`, `SFH.NORM_TEMP = LIGHT`, and `SFH.LIBRARY = EMILES_baseFe_safe/`. The verified wrapper therefore prepares templates at 20.5 km/s when the actual galaxy step is 41 km/s. Its LSF names are `lsf_MUSE-WFM-MILES` and `lsf_MILES`.

### 7.2 Inputs outside the galaxy product folder

The existing `further/BaSTI+Chabrier.dat` and `further/EMILES_baseFe_safe/` supply the substantial model resources. The desired LSF placement, following the user's instruction, is directly in `further/`:

```text
further/
  Mass.py
  BaSTI+Chabrier.dat
  EMILES_baseFe_safe/
  lsf_MUSE-WFM-MILES
  lsf_MILES
```

No additional support subdirectory is required. These LSF files are small text resources, not another large model library. The data LSF is 210300 bytes; the checked upstream template LSF is 131 bytes.

The observed cube, distance, stellar EBV map, and inclination are not needed for the bin-level M/L itself. BINID is needed afterward to map the bin results to spatial pixels.

## 8. Detailed reconstruction algorithm

### Step 1: Establish SSP identities

Read weight columns in their saved order. Match each grid entry to the original spectrum and its BaSTI row. Create an explicit mapping keyed by age, metallicity and alpha where applicable. Reproduce nGIST's grid/flattening convention or reorder the recovered factors back into the saved GRID order. Require unique, complete matches; file-list ordering alone is not sufficient.

### Step 2: Prepare the LSF functions

Use the verified `getLSF(config, "SFH")` behavior. The data curve undergoes nGIST's redshift transformation and its implemented resolution-floor adjustment. The template curve is read separately. Both interpolators use linear interpolation with extrapolation enabled. Preserve the actual code, including its floor rule, rather than replace it with a newly interpreted physical formula.

### Step 3: Reconstruct and crop the original wavelength grid

The MILES preparation routine reconstructs the grid from `CRVAL1`, `CDELT1`, and `NAXIS1`. It determines crop indices from the first sorted template and assumes the library shares that grid. A robust reconstruction should check that assumption without changing the resulting selection.

For each side it uses 150 Angstrom overhead when available. Otherwise the boundary is adjusted to remain approximately 5 Angstrom inside the library edge. The selection uses strict inequalities. With sufficient original coverage, 4800-7000 requests samples inside 4650-7150; 4800-8900 requests samples inside 4650-9050. Actual selected endpoints depend on grid sampling.

### Step 4: Convolve the cropped SSP

Evaluate both LSF functions on the retained original wavelengths. Reproduce the resolution check and calculate:

$$
\sigma_{\rm pix}(\lambda)=\frac{\sqrt{{\rm FWHM}_{\rm data}^2(\lambda)-{\rm FWHM}_{\rm SSP}^2(\lambda)}}{2.355\,\Delta\lambda_{\rm SSP}}.
$$

Apply pPXF's same `gaussian_filter1d` routine. If the data resolution is sharper than the template resolution, reproduce/report the original failure condition rather than silently clipping the square root.

### Step 5: Logarithmically rebin

Call the same `log_rebin` routine with the first/last retained wavelengths and template velocity step. A log grid has constant spacing in ln(lambda), corresponding to constant velocity sampling. It does not mean taking the logarithm of the flux. Galaxy and template spectra use compatible sampling, with the verified wrapper using templates sampled twice as finely.

### Step 6: Capture the mean before division

Set `n_j = np.mean(prepared_log_spectrum)`, immediately before nGIST divides that array by its mean. This is neither a mean over only good fitting pixels nor an R-filter integral. Do not replace it with `nanmean`; non-finite prepared input should be diagnosed. Record the identity and all preparation settings with the factor.

### Step 7: Construct R-band fractions

Read `ell_R,j = 10**(-0.4*R_j)` and `Upsilon_R,j = ML_R,j`. Calculate `c_R,j = ell_R,j/n_j`. For each bin form `q_R,bj = w_bj*c_R,j`, then divide by the bin's sum of q. The result is a normalized vector of R-band light fractions.

### Step 8: Calculate and map M/L

Compute `Upsilon_R,b = sum_j(f_R,bj*Upsilon_R,j)`. Map a bin's result to valid BINID pixels afterward. No inclination or distance factor belongs in this M/L conversion. Those parameters enter physical luminosity, mass, or surface density later.

### Process flowchart

```text
RUN CONFIG + SSP GRID + ORIGINAL SSPs + LSF CURVES
                         |
                         v
              Match each grid SSP uniquely
                         |
                         v
       Crop/pad -> LSF convolve -> log rebin
                         |
                         v
              Record prepared mean n_j
                         |
                         v
 BaSTI R magnitude -> ell_R,j/n_j <- saved weights w_bj
                         |
                         v
           q_R,bj = w_bj * ell_R,j/n_j
                         |
                         v
      f_R,bj = q_R,bj / sum_k(q_R,bk); sum f_R = 1
                         |
                         v
       M/L_R,b = sum_j(f_R,bj * BaSTI ML_R,j)
```

### Implementation and cache policy

Do not key a normalization cache only by 7000/8900 or library name. Redshift enters data-LSF processing, and the numerical result can depend on the LSF contents, template contents, velocity step, padding/grid selection and numerical implementation. Cache using these actual inputs or their hashes. Reuse factors across bins of one shared preparation setup; do not assume all galaxies have identical factors.

## 9. Executed GitHub and Setonix verification

### 9.1 Source comparison

The setup README identifies commit `038e9fa4839758e9dae53baeaf93a82bd75d8026` as nGIST 7.6.8. Parsed-function comparisons found the local 7.6.4 `prepareSpectralTemplateLibrary` and `getLSF` functions structurally identical to that commit, including their statements. Its SFH wrapper reads VELSCALE from the production spectra HDF5 and uses a ratio of two for template sampling. This establishes agreement of these checked functions, not all behavior across all nGIST versions.

### 9.2 Production files and identity

The executed SSH check located both production LSF files under:

```text
/software/projects/pawsey1308/ngist_supplementary_public/
    ngistTutorial/configFiles/
```

The production `lsf_MUSE-WFM-MILES` and upstream `tests/gistTutorial/configFiles/lsf_MUSE-WFM` at the checked commit both contain 210300 bytes and have SHA256:

```text
878a2b04299788669a8e97d4912967fbce6966920e533bf2a0ac5b581668e72e
```

They are byte-identical. The production data-file header describes a UDF-10-based curve with minimum FWHM modified to 2.54 Angstrom. The name difference does not encode a special nGIST mode.

### 9.3 Actual loader test on Setonix

The test ran through `ngistenv1308 python`, importing the installed nGIST reporting **7.6.8**. The loader source was `/work/ngist/ngistPipeline/auxiliary/_auxiliary.py` inside that environment.

| Executed test | Result |
|---|---|
| Read production data and template LSFs | Passed |
| Copy identical contents to arbitrary temporary filenames and update CONFIG | Passed; identical interpolated arrays at sampled wavelengths |
| Configure a deliberately nonexistent data-LSF file | FileNotFoundError, as expected |
| Temporary-directory cleanup | Completed automatically |

At wavelengths `[4800, 5500, 6500, 7000, 8500]` Angstrom and redshift 0.005566, the returned data FWHM values were approximately `[2.955903, 2.738083, 2.527179, 2.525941, 2.525941]` Angstrom; template FWHM was 2.51 throughout. These values test the loader and its redshift/floor handling, not the complete normalization reconstruction.

The tests completed with exit code zero. They did not submit a Slurm science job, fit a galaxy, change production inputs or regenerate products. Their log was appended to `Codex_log/2026_10_29.md`.

## 10. MILES/E-MILES range, padding and successful 9100 fits

The production instrumental curve is tabulated over 4000-10000 Angstrom. The template curve has endpoints `(4000, 2.51)` and `(9000, 2.51)`. The last table entry is not an enforced wavelength limit: nGIST allows extrapolation. Equal template endpoints extrapolate to the same 2.51 Angstrom FWHM at 9100 or 9250.

Consequently a 4800-9100 E-MILES fit, whose template preparation can extend to roughly 9250, can pass: the spectra exist, the data curve covers the padded wavelengths, and the template curve extrapolates. The particular historical 9100 fit was not independently audited here; the verified loader behavior explains how such a run can complete.

Numerical acceptance is distinct from the physical resolution prescription. E-MILES has approximately constant 2.51 Angstrom FWHM over its optical region around 3541-8950 Angstrom, with different resolution outside it [R4]. The 7000 preparation interval is comfortably inside that region. A padded 8900 or 9100 run can cross the transition. Exact historical normalization reconstruction must reproduce the LSF used in that fit; changing the physical resolution treatment would be a separate fitting-method change, potentially requiring new fits. A shorter original MILES library cannot provide missing red wavelengths merely because an LSF can be extrapolated.

## 11. Choosing V, R or I for the mass estimate

The measured cube luminosity and the fitted M/L have distinct wavelength requirements. The full cube can supply observed R light outside the fitting interval, while the fitted SSP mixture predicts the mass associated with that light. A 4800-7000 fit does not require using only V for post-processing.

Installed speclite Bessell response tables checked in this conversation have endpoints: V 4700-7000, R 5500-9000 and I 7000-9200 Angstrom. Endpoints include zero-response samples and are not rectangular passbands. Relative to a nominal 4750-9300 cube, V has a small blue-edge shortfall, R has comfortable margins, and I is nominally covered but close to the red edge. Exact comparisons require a consistent rest/observed frame.

R remains a reasonable primary choice, preserving continuity with Paper I, provided the relevant photometry is reliable. It is not intrinsically guaranteed to yield a more accurate mass than V. I is vulnerable to red sky-emission residuals and telluric correction: ESO notes residual sky emission beyond 7000 Angstrom, and the MUSE manual identifies the O2 A band around 7590-7710 Angstrom and red water-vapour absorption [R5, R6]. R also overlaps affected regions. Broadband integration can reduce some narrow-feature effects, but does not guarantee cancellation.

The inspected `Mass.py` photometry integrates cube DATA without subtracting emission lines in that step. If DATA includes nebular lines, R includes Halpha and neighbouring lines whose luminosity is not represented by a purely stellar M/L. This is a separate photometry issue, not solved by correcting n_j.

For the same population and consistent corrected luminosities, masses should satisfy `M = Upsilon_V L_V = Upsilon_R L_R`. A useful later check is a V/R mass comparison from the same fitted mixture. Its difference reflects observed-versus-model colour consistency, subject to edge treatment, attenuation, emission and calibration. No such comparison has been executed here.

## 12. What the two papers do and do not establish

Paper I Section 2.2 explicitly gives `Upsilon_R = sum_j(w_j Upsilon_R,j)` and describes the inputs as SFH light weights. Read literally, this resembles `light_r`; it does not explicitly describe the normalization-window-to-R conversion. That formula is appropriate if the input weights are R-band fractions. The published description alone does not establish the historical preprocessing or script that produced its maps [R2]. No quantitative conclusion about the impact on Paper I follows without reconstructing its actual inputs and comparing results.

Pessa et al. Equation 3 explicitly converts mass weights to V-light weights using `w_LumW = w/Upsilon_V`, then uses weighted means with a denominator [R1]. Their procedure is not equivalent to simply treating generic fitted-window light coefficients as V fractions. The algebra differs because their stated starting convention is mass weights. Band choice and starting weight convention are separate questions.

All consistent routes end with `Upsilon_band = sum_j(f_band,j Upsilon_band,j)`. What matters is establishing the band-light fractions correctly. Our refinement makes that conversion explicit for the current nGIST LIGHT convention.

## 13. Proposed implementation and acceptance checks

This is a specification, not a statement that the change is implemented.

1. Copy the actual production LSF files directly into `further/`; record hashes and source location. Keep the user's existing SSP spectra and photometry unchanged.
2. Load run-specific CONFIG and weight headers. Validate LIGHT mode, library identity and matching settings. Resolve local paths explicitly.
3. Reuse the verified nGIST preparation behavior to recover n_j immediately before each mean division. Record source revision and pPXF version; equivalent algorithms are the target, while bit-level historical identity also depends on numerical versions.
4. Join factors and table properties to the saved GRID order, with no missing or duplicate SSP identities.
5. Compute R contributions, normalized fractions, and M/L. Retain diagnostics for original raw sums and converted-fraction sums.
6. Compare the reconstructed normalized spectra against nGIST's own preparation output under identical inputs. Merely checking that a spectrum has mean one is insufficient to verify the preparation.
7. Validate algebra with synthetic amplitudes: w times the normalized SSP must equal (w/n) times the unnormalized prepared SSP. Confirm direct-ratio and fraction-average M/L calculations agree.
8. Require positive finite factors and denominators, fraction sums near one, and M/L within contributing SSP extrema for nonnegative weights.
9. Run a representative real-product comparison: current approximation versus reconstructed factors, reporting per-SSP factor ratios, per-bin delta log10(M/L_R), percentiles and spatial dependence. Test 7000 and 8900 inputs separately when actual products are available.
10. Only after those checks decide whether to regenerate mass products. Report model-arithmetic comparisons separately from observed mass changes and full scientific validation.

Exact preparation reproduction does not mean zero stellar-mass uncertainty. It fixes the coefficient-convention bookkeeping. It does not remove population degeneracies, attenuation uncertainty, emission contamination, flux-calibration uncertainty or limitations of the fitted populations. These are not reasons to avoid the refinement; they define what it does and does not solve.

## 14. Verification ledger and corrections to earlier statements

| Statement | Consolidated status |
|---|---|
| Current extended/Mass.py selects light_r | Incorrect; inspected source selects light_norm_to_r |
| Current helper exactly mirrors nGIST normalization | Too strong; it omits padding, convolution and log sampling |
| Missing stored factors prevents recovery | Incorrect; deterministic preparation can regenerate them |
| Local 7.6.4 necessarily blocks 7.6.8 reproduction | Not supported; the two checked functions match the identified 7.6.8 commit |
| Differently named production data LSF may differ from upstream | Resolved for checked files: byte-identical |
| nGIST requires a special -MILES filename | Incorrect; the configured filename is read literally |
| Last LSF table row is a hard wavelength limit | Incorrect; interpolation permits extrapolation |
| A completed 9100 run proves the physical LSF is correct | Incorrect; completion and physical applicability are different checks |
| Refined M/L has already been calculated | No; implementation and numerical impact remain pending |

Fresh report-stage checks re-read the current normalization/table code, IC3392 CONFIG and weights/grid structure, and installed preparation/LSF/SFH source. The GitHub and Setonix outcomes above are executed evidence from earlier turns in this same conversation on 29 October 2026, not new production reruns performed while formatting this report. No Codex memory is needed as evidence for this report.

## 15. Sources and reproducibility pointers

- **[R1] Pessa et al. (2023), Resolved stellar population properties of PHANGS-MUSE galaxies.** A&A 673, A147. DOI: https://doi.org/10.1051/0004-6361/202245673. Relevant: fitting interval and Equations 1-3 defining mass-to-light conversion for luminosity-weighted properties.
- **[R2] Huang et al. (2026), MAUVE-MUSE: when metallicity follows or fights star formation - a mass-dependent inversion in Virgo galaxies.** MNRAS 549, stag1019. DOI: https://doi.org/10.1093/mnras/stag1019. Section 2.2: https://academic.oup.com/mnras/article/549/3/stag1019/8698768.
- **[R3] nGIST source, pinned commit 038e9fa4839758e9dae53baeaf93a82bd75d8026.** Repository: https://github.com/geckos-survey/ngist/tree/038e9fa4839758e9dae53baeaf93a82bd75d8026. Relevant paths: `ngistPipeline/prepareTemplates/miles.py`; `ngistPipeline/auxiliary/_auxiliary.py`; `ngistPipeline/starFormationHistories/ppxf_sfh_wrapper.py`; `tests/gistTutorial/configFiles/lsf_MUSE-WFM`; `tests/gistTutorial/configFiles/lsf_MILES`.
- **[R4] Vazdekis et al. (2016), E-MILES stellar population models.** MNRAS 463, 3409-3436. DOI: https://doi.org/10.1093/mnras/stw2231. Relevant: wavelength-dependent spectral resolution of the extended models.
- **[R5] ESO MUSE pipeline note.** https://www.hq.eso.org/sci/publications/announcements/sciann16040.html. Relevant: residual sky emission redward of 7000 Angstrom.
- **[R6] ESO MUSE User Manual, version 12.1, Section 7.1.7.** https://www.hq.eso.org/sci/facilities/paranal/instruments/muse/doc/ESO-261650_MUSE_User_Manual.pdf. Relevant: telluric O2 and water-vapour absorption.

Local inspected sources: `/Users/Igniz/Desktop/ICRAR/further/Mass.py`; `/Users/Igniz/Desktop/ICRAR/extended/Mass.py`; `/Users/Igniz/Desktop/ICRAR/further/v3tk_v7.6.8/IC3392/CONFIG`; its `IC3392_sfh_weights.fits` and `LOGFILE`; local `EMILES_baseFe_safe/`; `/opt/miniconda3/envs/ngistenv/lib/python3.11/site-packages/ngistPipeline/`. The run log for the remote LSF verification is `/Users/Igniz/Desktop/Codex_log/2026_10_29.md`.
