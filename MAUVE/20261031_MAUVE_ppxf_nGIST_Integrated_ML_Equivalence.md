# Integrated stellar M/L from pPXF and nGIST LIGHT weights

Date: 2026-10-31 (Australia/Perth)

Scope: the pPXF stellar-population M/L calculation, the meaning of saved nGIST
LIGHT coefficients, and the `light_r` and `light_norm_to_r` implementations in
`further/Mass.py`.

## 1. Main result and adopted model consistency

The current `light_norm_to_r` calculation is mathematically equivalent to the
pPXF procedure for obtaining the stellar mass-to-light ratio of a fitted SSP
mixture in a specified photometric band. Both evaluate total stellar mass
divided by total band luminosity after accounting for the normalization of
each fitted template. Expressing this calculation as an R-light-weighted mean
of individual SSP M/L values introduces no mathematical approximation.

The project owner confirmed on 2026-10-31 that the MILES photometric table was
downloaded directly from the MILES website with the configuration matching
`EMILES_baseFe_safe`, the template library used for the nGIST products. This
report adopts that matched table/library configuration as an established
project input. It does not treat a suspected model mismatch as an unresolved
objection to the equivalence. The confirmation is user-supplied provenance;
this report does not claim to have independently reconstructed the download.

The three methods can be distinguished as follows:

| Method | Treatment of fitted coefficients | Mixture M/L |
|---|---|---|
| pPXF `sps_lib.mass_to_light` | Uses the normalization factors stored by its SPS helper; computes band luminosities from spectra | Total stellar mass / total band luminosity |
| `light_r` | Treats the saved coefficients directly as R-light fractions | Arithmetic mean of SSP M/L weighted by those coefficients |
| `light_norm_to_r` | Reconstructs nGIST's factors, converts saved LIGHT coefficients into R-light fractions | Exactly the total-mass / total-R-light expression |

For the checked nGIST LIGHT workflow, `light_norm_to_r` is the applicable
interpretation. `light_r` requires a different meaning for the coefficients
and is not generally correct for these saved weights.

Here, "integrated M/L" means the ratio for all SSP components of one fitted
spectrum, normally one spatial bin. A galaxy-wide ratio requires an additional
luminosity-weighted combination of spatial bins (Section 10).

## 2. Definitions and the general mass-to-light identity

For SSP index $j$, define:

| Symbol | Meaning |
|---|---|
| $S_j(\lambda)$ | Original physical SSP spectrum with the library's population normalization |
| $P_j(\lambda)$ | SSP after the preparation used for fitting, before LIGHT normalization |
| $n_j$ | Scalar used to divide that prepared SSP |
| $T_j(\lambda)=P_j(\lambda)/n_j$ | Normalized template supplied to the fit |
| $w_j$ | Coefficient multiplying $T_j$; an overall common scale is immaterial |
| $a_j=w_j/n_j$ | Coefficient of the SSP on its original population-amplitude scale |
| $M_j$ | Stellar mass associated with the original SSP population normalization |
| $L_{R,j}$ | R-band luminosity of that same SSP population |
| $\Upsilon_{R,j}=M_j/L_{R,j}$ | Individual SSP stellar M/L in solar R-band units |

Preparation changes sampling and spectral resolution; the SSP population
amplitude is tracked separately from the scalar LIGHT division. The population
mixture represented by the fit has prepared spectrum

$$
P_{\rm mix}(\lambda)
=\sum_j w_jT_j(\lambda)
=\sum_j\frac{w_j}{n_j}P_j(\lambda).
$$

The same population amplitudes multiply the corresponding stellar masses and
band luminosities. Therefore

$$
\boxed{
\Upsilon_R
=\frac{\sum_j a_jM_j}{\sum_j a_jL_{R,j}}
=\frac{\sum_j(w_j/n_j)M_j}{\sum_j(w_j/n_j)L_{R,j}}.
}
\tag{1}
$$

Multiplying every $w_j$ by the same positive constant changes neither ratio.
The coefficients $a_j$ are population amplitudes; surviving stellar-mass
contributions are $a_jM_j$. Calling $w_j/n_j$ "mass weights" does not remove
the need for the SSP mass array when surviving masses vary between SSPs.

## 3. How pPXF obtains integrated M/L from fitted light weights

### 3.1 The fit and the M/L helper are separate operations

Calling `ppxf(...)` determines coefficients and fitting parameters. The
stellar-population helper `sps_lib.mass_to_light(...)` subsequently combines
those coefficients with the SPS information. M/L is not an automatic output
of every `ppxf(...)` fit.

The pPXF release history records the generalization of this helper to arbitrary
SPS models, filters and redshifts in version 9.1.1, dated 2024-01-18. Its
synthetic-photometry approach removes the requirement for a precomputed
photometric M/L table, while still requiring SSP spectra and stellar masses.
See the [official pPXF changelog](https://pypi.org/project/ppxf/#changelog).

The release-note passage supplied in this discussion states:

```text
V9.1.1: MC, Oxford, 18 January 2024
sps_util.mass_to_light: Now it can calculate the stellar mass-to-light ratio (M*/L) for any stellar population synthesis (SPS) model, any filter, and any redshift, using the output weights from pPXF. No need for pre-computed tables anymore.
```

Thus pPXF provides a function that directly returns the mixture M/L. Its
availability does not mean that another pipeline calls it automatically.

### 3.2 Template normalization in the inspected pPXF 9.4.2 SPS helper

The helper loads the spectral library, wavelength grid, stellar masses, ages,
metallicities and spectral resolutions. When resolution matching is requested,
it broadens the spectra, then logarithmically rebins them.

With `norm_range` specified and `norm_type='mean'`, it measures each template's
mean over the selected log-sampled wavelength interval and divides that
template by its own factor. It retains those factors in `sps.flux`. The fit
then returns coefficients of these individually normalized templates.

With `norm_range=None`, the inspected implementation instead divides all
templates by a common positive median. It retains that common factor in an
array. The same M/L bookkeeping works for both cases.

The SPS helper also retains `templates_full`. Here "full" refers to retaining
the wavelength extent before optional `lam_range` trimming; these spectra have
already undergone the helper's normalization. The distinction matters when
reading the M/L denominator.

The source makes the stored normalization and mass arrays explicit:

```python
self.flux = flux            # factor by which each template was divided
self.mass_no_gas_grid = masses
```

These assignments are in
[sps_util.py, lines 222--223](/opt/miniconda3/envs/ICRAR/lib/python3.13/site-packages/ppxf/sps_util.py:222).
In the notation of this report, `self.flux[j]` is $n_j$ and
`self.mass_no_gas_grid[j]` is $M_j$. The attribute name `flux` should not be
mistaken for the R-band luminosity: it records the normalization divisor.

### 3.3 The M/L calculation

In the inspected [local pPXF source](/opt/miniconda3/envs/ICRAR/lib/python3.13/site-packages/ppxf/sps_util.py:272),
the helper performs synthetic photometry on the normalized SSP spectra and on
a solar reference spectrum in the requested band and redshift convention.
Its conversion through the flux at 10 pc produces luminosities in solar band
units. The stellar mass array represents living stars plus stellar remnants,
excluding gas returned during stellar evolution.

The decisive calculation is reproduced directly from
[sps_util.py, lines 298--300](/opt/miniconda3/envs/ICRAR/lib/python3.13/site-packages/ppxf/sps_util.py:298):

```python
mass_weights = weights/self.flux    # Revert possible templates normalization
lum = p1.flux/p2.flux               # Lum in solar luminosities
mlpop = np.sum(mass_weights*self.mass_no_gas_grid)/np.sum(weights*lum)
```

The code-to-equation correspondence is:

| Code quantity | Mathematical quantity | Interpretation |
|---|---|---|
| `weights` | $w_j$ | Coefficients multiplying the normalized fitting templates |
| `self.flux` | $n_j$ | Per-template normalization divisors |
| `mass_weights` | $w_j/n_j$ | Coefficients restored to the original SSP population scale |
| `self.mass_no_gas_grid` | $M_j$ | Stellar mass for that original SSP population |
| `p1.flux` | Band flux of the normalized SSP, after the code's unit conversion | Obtained by synthetic photometry of `self.templates_full` |
| `p2.flux` | Solar reference band flux | Converts the SSP band flux into solar band-luminosity units |
| `lum` | $\widetilde L_{R,j}=L_{R,j}/n_j$ | Band luminosity of the already normalized SSP |
| `mlpop` | $\Upsilon_R$ | Total mixture stellar mass divided by total mixture band luminosity |

In particular, `lum` is a band-integrated photometric quantity, not the
template's wavelength-by-wavelength flux array. The general helper uses the
requested band; the $R$ subscript in this table describes its R-band use.

Let $\widetilde L_{R,j}$ denote the luminosity measured from the normalized
template. Linearity of the bandpass measurement gives

$$
\widetilde L_{R,j}=\frac{L_{R,j}}{n_j}.
$$

The helper combines mass contributions proportional to $(w_j/n_j)M_j$ and
luminosity contributions proportional to $w_j\widetilde L_{R,j}$. Thus

$$
\Upsilon_R^{\rm pPXF}
=\frac{\sum_j(w_j/n_j)M_j}{\sum_jw_j\widetilde L_{R,j}}
=\frac{\sum_j(w_j/n_j)M_j}{\sum_j(w_j/n_j)L_{R,j}}.
\tag{2}
$$

The denominator does not need another explicit division by $n_j$: the spectra
used for its synthetic photometry were already divided by $n_j$.

The default pPXF band is `SDSS/r`. An R-band comparison must select the same
response and luminosity convention as the adopted R-band table. The physical
formula applies to either band; the names `r` and `R` alone do not identify
identical filters. Model matching confirmed by the project owner is accepted
here; a comparison to a different default pPXF filter would be a different
measurement.

## 4. What nGIST saves, and what its LIGHT weights mean

The checked nGIST source is commit
`038e9fa4839758e9dae53baeaf93a82bd75d8026`. The actual IC3392 log identifies
nGIST 7.6.8. The relevant standard SFH method is `ppxf` and its template
preparation uses `TEMPLATE_SET: miles` and `NORM_TEMP: LIGHT`.

The [template-preparation source](https://github.com/geckos-survey/ngist/blob/038e9fa4839758e9dae53baeaf93a82bd75d8026/ngistPipeline/prepareTemplates/miles.py)
divides each prepared SSP by its mean. Its coefficients consequently describe
contributions in this mean-flux normalization convention. An unweighted mean
on a logarithmic wavelength grid is neither a filter-weighted R luminosity
nor a bolometric luminosity integral.

The [SFH wrapper](https://github.com/geckos-survey/ngist/blob/038e9fa4839758e9dae53baeaf93a82bd75d8026/ngistPipeline/starFormationHistories/ppxf_sfh_wrapper.py)
reshapes the final fitted coefficients and divides them by their total before
returning and saving them. In symbols, if $u_{bj}$ are the coefficients
returned by pPXF for bin $b$, the saved values are

$$
w_{bj}=\frac{u_{bj}}{\sum_k u_{bk}}.
\tag{3}
$$

The corresponding source statements are
([ppxf_sfh_wrapper.py, lines 426--435](https://github.com/geckos-survey/ngist/blob/038e9fa4839758e9dae53baeaf93a82bd75d8026/ngistPipeline/starFormationHistories/ppxf_sfh_wrapper.py#L426-L435)):

```python
normalized_weights = pp.weights / np.sum( pp.weights )

weights = pp.weights.reshape(templates.shape[1:])/pp.weights.sum()
w_row   = np.array([np.reshape(weights, ncomb)])
```

These statements are selected from the same block, with unrelated intervening
lines omitted and trailing comments removed. `normalized_weights` is used to
build the optimal template. `weights` and `w_row` hold the final normalized
coefficient distribution that is returned and saved. They are two uses of the
same sum-to-one convention, not successive mass or R-band conversions.

This operation preserves every relative coefficient. It does not convert the
mixture into R-light fractions or surviving stellar-mass fractions.

Source inspection traced the returned array through serial and parallel
collection to the FITS writer. The writer saves the coefficients in `WEIGHTS`
and the SSP identities in `GRID`. The SFH map writer copies population and
diagnostic columns into maps. The checked pipeline does not call the pPXF M/L
helper or perform an equivalent hidden conversion before saving the weights.

The absence of a conversion is supported by the whole checked path, not only
the name or location of pPXF's helper: the coefficient operations above contain
no $n_j$, $L_{R,j}$ or $M_j$ factors, the pipeline source search found no call
to `mass_to_light`, and the writer stores `w_row` as the `WEIGHTS` column.
The presence of the helper in pPXF's `sps_util` module therefore provides a
capability that the checked nGIST workflow does not invoke.

The inspected IC3392 products corroborate that path:

| Product | Contents |
|---|---|
| `IC3392_sfh_weights.fits` | `WEIGHTS`; `GRID` containing `LOGAGE`, `METAL`, `ALPHA` |
| `IC3392_sfh_maps.fits` | `AGE`, `METAL`, `ALPHA`, `SNR_POSTFIT`, `RED_CHI2`, `EBV` |

Its weights have shape $(4062,477)$, their row sums are unity within floating
point precision, and their header specifies `NORM_TEMP=LIGHT` and
`LIBRARY=EMILES_baseFe_safe/`. No M/L quantity is stored in these two products.
The local intermediate `_sfh.fits` was absent, so its structure was established
from the writer source rather than a direct file read.

## 5. The `light_r` calculation and its required assumption

The `light_r` branch computes

$$
\Upsilon_{R,b}^{\texttt{light\_r}}
=\frac{\sum_jw_{bj}\Upsilon_{R,j}}{\sum_jw_{bj}}.
\tag{4}
$$

An arithmetic mean of SSP M/L values is correct when its weights are the
SSPs' R-band light fractions. However, nGIST's saved LIGHT coefficients are
coefficients of mean-normalized spectral templates. Their associated R-light
contributions are proportional to $w_{bj}L_{R,j}/n_j$, not simply $w_{bj}$.

Consequently `light_r` generally uses the wrong weighting measure for this
workflow. It can agree with the proper result for a single contributing SSP,
identical SSP M/L values, a mixture with constant $L_{R,j}/n_j$ among its
contributors, or an accidental cancellation. Such agreement does not establish
the general interpretation. Its error can be positive or negative; no universal
sign follows from the algebra.

## 6. The `light_norm_to_r` calculation and proof of equivalence

### 6.1 Recover R-band contributions

For each bin, the current
[Mass.py helper](/Users/Igniz/Desktop/ICRAR/further/Mass.py:1086) computes

$$
q_{bj}=w_{bj}\frac{\ell_{R,j}}{n_j},
\qquad
f_{R,bj}=\frac{q_{bj}}{\sum_kq_{bk}},
\tag{5}
$$

where the table-based proxy is

$$
\ell_{R,j}=10^{-0.4\mathcal M_{R,j}}.
\tag{6}
$$

Here $\mathcal M_{R,j}$ is the SSP absolute R magnitude. Physical band
luminosity is $L_{R,j}=C_R\ell_{R,j}$ for one common positive constant $C_R$.
That constant cancels in Equation (5), so no approximation results from using
the magnitude-derived luminosity proxy to obtain fractions.

### 6.2 Combine the SSP M/L values

The final calculation is

$$
\Upsilon_{R,b}^{\texttt{light\_norm\_to\_r}}
=\sum_j f_{R,bj}\Upsilon_{R,j}.
\tag{7}
$$

Substitution gives the complete equivalence:

$$
\begin{aligned}
\Upsilon_{R,b}^{\texttt{light\_norm\_to\_r}}
&=\frac{\sum_j(w_{bj}/n_j)L_{R,j}\Upsilon_{R,j}}
        {\sum_j(w_{bj}/n_j)L_{R,j}}\\
&=\frac{\sum_j(w_{bj}/n_j)L_{R,j}(M_j/L_{R,j})}
        {\sum_j(w_{bj}/n_j)L_{R,j}}\\
&=\frac{\sum_j(w_{bj}/n_j)M_j}
        {\sum_j(w_{bj}/n_j)L_{R,j}}\\
&=\Upsilon_{R,b}^{\rm pPXF}.
\end{aligned}
\tag{8}
$$

The current `Mass.py` implements Equations (5)--(7) directly:

```python
q = weights * (luminosity_r/norm_flux)[None, :]
denominator = q.sum(axis=1)
fractions = q/denominator[:, None]
ml = np.sum(fractions*ml_r[None, :], axis=1)
```

In this excerpt, `weights` is the saved $w_{bj}$ array, `norm_flux` is $n_j$,
`luminosity_r` is the MILES-table proxy $\ell_{R,j}$, `ml_r` is the table's
$\Upsilon_{R,j}$, and `ml` is the resulting $\Upsilon_{R,b}$. The actual
implementation also validates the inputs and resulting fractions. See
[r_light_mass_to_light](/Users/Igniz/Desktop/ICRAR/further/Mass.py:1086).

Replacing the intermediate fractions with their definition makes the source
expression explicit:

$$
\boxed{
\Upsilon_{R,b}
=\frac{\sum_j(w_{bj}/n_j)\ell_{R,j}\Upsilon_{R,j}}
        {\sum_j(w_{bj}/n_j)\ell_{R,j}}
=\frac{\sum_j(w_{bj}/n_j)M_j}
        {\sum_j(w_{bj}/n_j)L_{R,j}}.
}
$$

The last equality uses $L_{R,j}=C_R\ell_{R,j}$ and
$M_j=L_{R,j}\Upsilon_{R,j}$. The common $C_R$ cancels. This connects the
literal `Mass.py` operations to the literal pPXF calculation above.

Thus the MILES-table route and the pPXF spectral-photometry route implement
the same population-mixture identity. Given the adopted matched configuration,
`light_norm_to_r` is a direct implementation of that identity, not an empirical
correction calibrated against the `light_r` result. A numerical comparison of
table photometry with newly integrated spectra can still show tabulation,
rounding or integration differences; these do not change the identity.

### 6.3 Why restoring the population coefficient requires division

The fitted term is $w_j(P_j/n_j)$. Its coefficient of $P_j$ is therefore
$w_j/n_j$. Multiplying the coefficient by $n_j$ would reverse the required
operation. Restoring the spectrum itself and restoring its coefficient are
different bookkeeping operations.

### 6.4 A transparent two-SSP example

Take illustrative SSP values with consistent arbitrary population units:

| SSP | $w_j$ | $n_j$ | $L_{R,j}$ | $\Upsilon_{R,j}$ | $M_j=L_{R,j}\Upsilon_{R,j}$ |
|---|---:|---:|---:|---:|---:|
| 1 | 0.5 | 1 | 1 | 1 | 1 |
| 2 | 0.5 | 2 | 4 | 5 | 20 |

`light_r` gives $(1+5)/2=3$. The actual R contributions are $(0.5,1)$,
so the R fractions are $(1/3,2/3)$ and `light_norm_to_r` gives $11/3$.
The direct ratio gives the identical answer:

$$
\frac{0.5\times1+0.25\times20}{0.5\times1+0.25\times4}
=\frac{5.5}{1.5}=\frac{11}{3}.
$$

These values illustrate the weighting distinction, not a measured galaxy.

## 7. Reconstructing the exact nGIST normalization factors

### 7.1 Read the selected run's metadata

`load_run_setup` reads the selected SFH weights and its own accompanying
`CONFIG`. It cross-checks `LMIN`, `LMAX`, `NORM_TEMP`, `TEMPLATE_SET`, `LIBRARY`
and `LSF_TEMP` between the configuration and saved SFH header. It also checks
that `RUN_ID` matches the selected weights filename.

The reconstruction requires the original SSP library, the data and template
line-spread-function (LSF) tables, systemic redshift, and configured velocity
sampling. Galaxy sampling comes directly from
`CONFIG.PREPARE_SPECTRA.VELSCALE`. The checked nGIST SFH wrapper uses
`velscale_ratio=2`, so template sampling is

$$
\Delta v_{\rm template}=\Delta v_{\rm galaxy}/2.
\tag{9}
$$

For the checked IC3392 run these values are 41 and 20.5 km/s. Intermediate
spectra HDF5 files are not needed to recover this configured sampling.

### 7.2 Reproduce the padded linear-wavelength crop

The original SSP grid is reconstructed from `CRVAL1`, `CDELT1` and `NAXIS1`.
Let its endpoints be $\lambda_{\rm start}$ and $\lambda_{\rm end}$, and let
the SFH limits be $\lambda_{\min}$ and $\lambda_{\max}$. nGIST defines

$$
h_{\rm lo}=\begin{cases}
150 & \lambda_{\min}-\lambda_{\rm start}>150,\\
\lambda_{\min}-\lambda_{\rm start}-5 & \text{otherwise},
\end{cases}
$$

$$
h_{\rm hi}=\begin{cases}
150 & \lambda_{\rm end}-\lambda_{\max}>150,\\
\lambda_{\rm end}-\lambda_{\max}-5 & \text{otherwise}.
\end{cases}
$$

All quantities in these crop rules are in Angstrom. The retained samples obey
strict inequalities:

$$
\lambda_{\min}-h_{\rm lo}<\lambda<\lambda_{\max}+h_{\rm hi}.
\tag{10}
$$

The reconstruction uses the same `np.arange` grid construction and the actual
first and last retained samples as the logarithmic-rebinning bounds. The
normalization is therefore evaluated over the padded prepared spectrum,
not simply over the nominal fitting interval.

### 7.3 Reproduce the run's LSF treatment

`load_lsf` reads wavelength/FWHM tables and follows the checked nGIST
`getLSF` implementation. It divides the data-LSF wavelength and FWHM by
$1+z$. The checked reference then applies its historical floor rule: entries
whose transformed FWHM divided by $1+z$ is below 2.51 Angstrom are assigned
$2.54(1+z)$ Angstrom. The template LSF table is not subjected to that same
data-LSF transformation. Both curves use linear interpolation with
extrapolation.

The second division in the reference floor condition is recorded deliberately.
Reconstruction must reproduce the actual preparation used by the pipeline;
silently replacing that expression with a preferred physical rule would
reconstruct a different template normalization. This is an implementation
description, not an independent endorsement of that floor prescription.

At each retained original SSP wavelength, the Gaussian broadening width in
linear-grid pixels is

$$
\sigma_{\rm pix}(\lambda)
=\frac{\sqrt{{\rm FWHM}_{\rm data}(\lambda)^2
                    -{\rm FWHM}_{\rm template}(\lambda)^2}}
       {2.355\,\Delta\lambda}.
\tag{11}
$$

The reconstruction stops if the template is broader than the data at a
retained wavelength, following the reference's rejection condition. It uses
the same pPXF `gaussian_filter1d` routine for wavelength-dependent smoothing.

### 7.4 Logarithmically rebin, then measure the mean

After cropping and smoothing each SSP, the reconstruction calls pPXF
`log_rebin` with the retained endpoints and the template velocity step from
Equation (9). With $N$ output samples, the required factor is

$$
\boxed{n_j=\frac{1}{N}\sum_{p=1}^{N}P_{pj}.}
\tag{12}
$$

This is the arithmetic mean immediately before nGIST's LIGHT division.
`Mass.py` uses `np.mean`, not `nanmean`. Non-finite prepared spectra and
non-positive normalization factors are rejected.

The nGIST preparation excerpt makes the order of operations visible
([miles.py, lines 248--260](https://github.com/geckos-survey/ngist/blob/038e9fa4839758e9dae53baeaf93a82bd75d8026/ngistPipeline/prepareTemplates/miles.py#L248-L260)):

```python
ssp = gaussian_filter1d(ssp, sigma)
sspNew, logLam2, _ = log_rebin(
    lamRange_spmod, ssp, velscale=velscale
)

templates[:, j, k, i] = sspNew / np.mean(sspNew)
```

This excerpt omits grid assignments and shows the last statement from inside
the `NORM_TEMP == "LIGHT"` branch. The `velscale` passed into this preparation
routine is the template step, already set to the galaxy step divided by two.
The exact divisor that must be recovered is `np.mean(sspNew)`.

The reconstruction in
[Mass.py](/Users/Igniz/Desktop/ICRAR/further/Mass.py:1020) uses the same
convolution and rebinning routines after reconstructing the crop and LSF:

```python
prepared, _, _ = log_rebin(bounds, gaussian_filter1d(cropped, sigma), velscale=template_velscale)
mean = float(np.mean(prepared))
lookup[identity] = mean
```

These statements are selected from the function with validation code omitted.
`prepared` corresponds to nGIST's `sspNew`, and `mean` is the reconstructed
$n_j$. The lookup is then reordered to the saved SSP grid before its values
enter the M/L calculation.

The full sequence is

```text
Original SSP FITS spectrum
    -> strict padded wavelength crop
    -> wavelength-dependent LSF convolution on the linear grid
    -> logarithmic rebinning at CONFIG VELSCALE / 2
    -> arithmetic mean of the resulting spectrum = n_j
    -> nGIST divides the prepared spectrum by n_j before fitting
```

`Mass.py` repeats the steps through measurement of $n_j$; it does not refit the
observed spectrum. It then uses the already saved weights in Equation (5).

### 7.5 Match all quantities to the saved SSP order

Template filenames are parsed into age, metallicity and alpha identities. The
saved `GRID` supplies `LOGAGE`, `METAL` and `ALPHA`; age is recovered as
$10^{\rm LOGAGE}$ in Gyr. The current lookup matches rounded identities in
the saved grid order rather than assuming that filesystem order equals
coefficient order. The matched MILES table supplies R magnitude and R-band
M/L for the adopted baseFe configuration.

Checks reject duplicate template/grid identities, missing required SSPs,
inconsistent original wavelength grids, invalid LSF values, non-finite spectra,
negative weights, and invalid luminosity or normalization factors. The final
fractions must sum to unity and the mixture M/L must lie between the minimum
and maximum M/L of its contributing SSPs. Provenance includes configuration,
LSF and template hashes, numerical-package versions and the Mass.py hash.

### 7.6 Why a raw-window mean was insufficient

An earlier implementation took the mean of the original linearly sampled SSP
within the nominal fitting limits. It omitted the padding, LSF preparation and
logarithmic sampling. Even without smoothing, a mean on an equally spaced log
grid generally differs from a mean on an equally spaced linear grid. In a
continuous approximation they weight wavelength as

$$
\langle S\rangle_{\log\lambda}
\simeq\frac{\int S(\lambda)\,d\lambda/\lambda}
                  {\ln(\lambda_2/\lambda_1)},
\qquad
\langle S\rangle_{\lambda}
\simeq\frac{\int S(\lambda)\,d\lambda}{\lambda_2-\lambda_1}.
$$

The current implementation reproduces the actual discrete operations rather
than substituting either continuous approximation. Neither mean should be
identified automatically with R-band luminosity.

## 8. What is reconstructed and what remains supplied by the table

The normalization factor is measured over the prepared fitting spectrum. The
R-band SSP luminosity and SSP M/L are taken from the matched photometric table.
These serve different purposes: reconstructing a fitted-template scale does
not require deriving the full R luminosity from the truncated fitting window.
In particular, a fit ending at 7000 Angstrom can still use a full-band R
luminosity supplied by its matched SSP table.

pPXF's own `sps_lib` preparation and nGIST's preparation need not use identical
cropping, broadening algorithms or normalization intervals. A fresh default
`sps_lib` object is therefore not automatically the correct normalization
object for previously saved nGIST weights. One must use nGIST's $n_j$, or
transform the coefficients consistently to the fresh object's convention.

If a different normalization uses $n'_j$ for the same underlying prepared
templates, preserving the same mixture requires

$$
w'_j=w_j\frac{n'_j}{n_j},
\qquad \frac{w'_j}{n'_j}=\frac{w_j}{n_j}.
\tag{13}
$$

The M/L is then invariant. This invariance concerns re-expressing a fixed
mixture. Refitting with regularization on a differently normalized weight
distribution can produce a different optimum; the present conversion performs
no such refit.

## 9. Numerical verification and its evidential scope

The checks below were executed on 2026-10-30 and their saved results were
inspected while preparing this report on 2026-10-31. No science pipeline was
rerun for this documentation task.

### 9.1 Reproduction of nGIST template preparation

The existing validator ran the original pinned nGIST preparation functions
alongside the definitions extracted directly from `Mass.py`. For the IC3392
4800--8900 Angstrom run:

| Quantity | Result |
|---|---:|
| Saved spatial bins | 4062 |
| SSPs | 477 |
| Retained prepared-spectrum endpoints | 4650.2--9049.4 Angstrom |
| Galaxy / template velocity step | 41 / 20.5 km/s |
| Maximum relative difference in normalization factors | 0 |
| Maximum absolute difference in normalized template samples | 0 |

These exact numerical matches apply to the checked inputs and software stack:
pPXF 9.4.2, NumPy 2.2.6, SciPy 1.15.2 and Astropy 7.0.2. They do not establish
bitwise identity with every historical production environment.

### 9.2 Direct exercise of the pPXF M/L helper

The separate check used the unmodified installed pPXF helper, a subset of its
bundled E-MILES SPS data, two prescribed mixtures, and three normalization
choices: a common scalar, a mean over 5070--5950 Angstrom, and a mean over
4800--8900 Angstrom. It transformed coefficients consistently for each
normalization and compared pPXF's output, the extracted `Mass.py` helper, and
the direct physical mass/luminosity ratio.

All six cases agreed within relative tolerance $3\times10^{-14}$. For one
mixture under 5070--5950 Angstrom normalization, the results were:

| Calculation | M/L in the test's SDSS/r solar units |
|---|---:|
| pPXF helper | 0.6661734636912048 |
| `light_norm_to_r` algebra | 0.6661734636912047 |
| Direct mass/luminosity ratio | 0.6661734636912048 |
| Unconverted `light_r` assumption | 0.5886175625176435 |

This prescribed-mixture test isolates the M/L arithmetic; it does not perform
a new spectral fit. Its bundled SPS dataset and SDSS/r filter are a separate
algorithm demonstration, not a replacement for the adopted project library
and R-band table.

The installed Accelerate-backed numerical routines emitted floating-point
warnings during initial photometric matrix products. Independent elementwise
quadrature of the solar and SSP fluxes yielded finite values and agreed with
the pPXF results within the same tolerance. The final check includes filter
coverage, finite-value and independent-quadrature assertions.

### 9.3 Direct mass/luminosity ratio for the saved real mixture

Using IC3392's actual saved weights, the reconstructed factors and the current
MILES-table quantities, the direct expression in Equation (1) and the
R-fraction expression in Equation (7) agreed across all 4062 bins. The maximum
relative difference was $6.66\times10^{-16}$.

For the same saved mixtures, the median of
$100(\Upsilon_R^{\texttt{light\_r}}/
\Upsilon_R^{\texttt{light\_norm\_to\_r}}-1)$ was -2.1809 percent.
The sign was negative in 4055 bins and positive in 7 bins. This is a
within-run coefficient-interpretation comparison, not an observed change in
the galaxy or a galaxy-integrated mass correction.

### 9.4 Reproduction commands and evidence

The commands used the ICRAR environment and wrote temporary outputs:

```bash
cd /Users/Igniz/Desktop/ICRAR/further
MPLCONFIGDIR=/private/tmp/mass_norm_refinement/mpl \
  /opt/miniconda3/envs/ICRAR/bin/python validate_mass_light_normalization.py \
  --weights v3tk_v7.6.8/IC3392/IC3392_sfh_weights.fits \
  --ngist-source-dir /private/tmp/mass_norm_refinement \
  --output-dir /private/tmp/ppxf_ml_audit_20261030

MPLCONFIGDIR=/private/tmp/mass_norm_refinement/mpl \
  /opt/miniconda3/envs/ICRAR/bin/python \
  /private/tmp/ppxf_ml_audit_20261030/check_ppxf.py
```

The first command needs the pinned `miles.py` and `_auxiliary.py` in its source
directory. The second check consumes the normalization arrays written by the
first. The scratch script and JSON links below are local evidence and may be
removed by temporary-directory cleanup:

- [pPXF comparison script](/private/tmp/ppxf_ml_audit_20261030/check_ppxf.py)
- [pPXF and real-mixture comparison JSON](/private/tmp/ppxf_ml_audit_20261030/ppxf_comparison.json)
- [Fresh IC3392 reference-comparison JSON](/private/tmp/ppxf_ml_audit_20261030/validation_8900.json)
- [Persistent project validator](/Users/Igniz/Desktop/ICRAR/further/validate_mass_light_normalization.py)
- [Existing project validation report](/Users/Igniz/Desktop/ICRAR/further/mass_light_normalization_validation_20261029/README.md)

## 10. Interpretation boundaries

The equivalence establishes the M/L of the stellar population mixture implied
by the saved fit and adopted SSP models. It does not by itself measure the
uncertainty in the fitted weights, test a different IMF, or validate the
observed luminosity subsequently used to calculate stellar mass. Dust,
continuum/line separation, photometric coverage and distance belong to that
subsequent luminosity calculation. Multiplicative fitting polynomials are not
additional SSP populations with their own physical stellar masses; the
reported population M/L is defined from the inferred SSP mixture.

If $L_{R,b}$ is the consistently defined luminosity of spatial bin $b$, then

$$
M_{*,b}=\Upsilon_{R,b}L_{R,b},
\qquad
\Upsilon_{R,\rm galaxy}
=\frac{\sum_b\Upsilon_{R,b}L_{R,b}}{\sum_bL_{R,b}}.
\tag{14}
$$

An unweighted mean or median of bin M/L values is not generally the
galaxy-integrated M/L. Neither the per-bin error statistics above nor the
equivalence proof provides a measured galaxy-integrated mass change.

For the confirmed matching MILES table and `EMILES_baseFe_safe` configuration,
the methodological conclusion is firm: reconstruct the actual nGIST
normalization factors, convert saved LIGHT coefficients into R-light
fractions, and combine the matched SSP M/L values. This is the same
total-stellar-mass / total-R-luminosity calculation implemented by pPXF.

## 11. Source and provenance register

1. [pPXF official documentation and changelog](https://pypi.org/project/ppxf/).
   Used for package context and the 9.1.1 release-history statement. Detailed
   algorithm descriptions above follow the inspected installed source.
2. [Installed pPXF 9.4.2 SPS helper](/opt/miniconda3/envs/ICRAR/lib/python3.13/site-packages/ppxf/sps_util.py:145)
   and [synthetic photometry implementation](/opt/miniconda3/envs/ICRAR/lib/python3.13/site-packages/ppxf/ppxf_util.py:919).
3. [Pinned nGIST miles preparation](https://github.com/geckos-survey/ngist/blob/038e9fa4839758e9dae53baeaf93a82bd75d8026/ngistPipeline/prepareTemplates/miles.py),
   [LSF preparation](https://github.com/geckos-survey/ngist/blob/038e9fa4839758e9dae53baeaf93a82bd75d8026/ngistPipeline/auxiliary/_auxiliary.py),
   [SFH wrapper](https://github.com/geckos-survey/ngist/blob/038e9fa4839758e9dae53baeaf93a82bd75d8026/ngistPipeline/starFormationHistories/ppxf_sfh_wrapper.py),
   and [map writer](https://github.com/geckos-survey/ngist/blob/038e9fa4839758e9dae53baeaf93a82bd75d8026/ngistPipeline/writeFITS/save_maps_fits.py).
   Local checkout inspected at `/private/tmp/ngist-sfh16-source`.
4. [Current Mass.py](/Users/Igniz/Desktop/ICRAR/further/Mass.py:951),
   [MILES photometric table](/Users/Igniz/Desktop/ICRAR/further/BaSTI+Chabrier.dat),
   and `EMILES_baseFe_safe` under the same `further` directory.
5. [IC3392 run CONFIG](/Users/Igniz/Desktop/ICRAR/further/v3tk_v7.6.8/IC3392/CONFIG),
   [weights](/Users/Igniz/Desktop/ICRAR/further/v3tk_v7.6.8/IC3392/IC3392_sfh_weights.fits),
   and [population maps](/Users/Igniz/Desktop/ICRAR/further/v3tk_v7.6.8/IC3392/IC3392_sfh_maps.fits).
6. Project-owner confirmation, 2026-10-31: the directly downloaded MILES table
   matches the configuration of the nGIST `EMILES_baseFe_safe` template library.
7. [Earlier normalization-refinement report](/Users/Igniz/Desktop/ICRAR/MAUVE/20261029_MAUVE_ML_R_Normalization_Refinement.md).
   This new report provides the explicit pPXF equivalence and adopts the
   subsequently confirmed table/library match.

Documentation prepared with AI assistance from the current source files,
saved numerical verification records and the project owner's model-provenance
confirmation. No production science code or products were changed for this
report.
