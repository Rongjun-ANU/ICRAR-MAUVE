# Publisher photometry findings and mass recommendation

Research and numerical follow-up: 2026-10-31. This supplements the executed SSP
photometry notebook; it does not replace its recorded inputs or outputs. No
production Mass.py, nGIST, SSP, or cube products were modified.

## New primary-source findings

The working official pages distinguish two different products:

- Photometric-library predictions:
  https://research.iac.es/proyecto/miles/pages/predicted-masses-and-photometric-observables-based-on-photometric-libraries.php
  These specify Vega U=B=V=R=0.03 mag and I=0.07 mag. They supply a
  BaSTI/baseFe/Chabrier photometric table and a separate stellar mass table.
- E-MILES SED predictions:
  https://research.iac.es/proyecto/miles/pages/photometric-predictions-based-on-e-miles-seds.php
  These specify Vega zero in all bands, CALSPEC alpha_lyr_stis_008.fits,
  solar R magnitude 4.459, and downloadable responses and SED photometry.

Downloaded resources are saved in publisher_resources/. Base URL:
https://research.iac.es/proyecto/miles/media/tarfiles/SSP_models/v11.0/
Files: out_phot_CH_BASTI, out_mass_CH_BASTI,
UBBVRIJHK_ch_iTp0.00.MAG, filters_UBBVRIJHK.res, HEADER_out_mass.
The CALSPEC reference was retrieved from
https://ssb.stsci.edu/trds/calspec/alpha_lyr_stis_008.fits.

The local BaSTI+Chabrier.dat matches the official out_phot_CH_BASTI in every
data token: 636 rows by 34 columns. This establishes its source product.
The mass table specifies M(*+remn), living stars plus remnants, separately from
M*, Mremn, and returned Mgas, per initial one solar mass. All 477 fitted SSP
age/metallicity keys have corresponding mass entries.

## Numerical cross-checks

1. The official R_Cousins response has 24 samples spanning 5500-9000 A.
   Its tabulated response is exactly identical to the installed speclite
   bessell-R response. The pPXF-selected Johnson-Cousins_R response has 159
   samples spanning 5475-8635 A and is a different curve. Comparing response
   arrays is insufficient by itself: the wavelength weighting must also match.
2. Evaluated all 477 actual safe SSP spectra on Setonix using the official
   R_Cousins response and CALSPEC reference, with the original notebook's
   independently checked integrator. Both photon and energy weighting were
   evaluated. No SFH fits were repeated.
3. The publisher SED table has valid R predictions for 403 of these SSPs;
   74 entries at ages 0.03-0.25 Gyr use -99 sentinels. Sentinels were excluded
   from comparisons, not interpreted as magnitudes.
4. Energy-weighted synthetic magnitudes reproduce those 403 publisher SED R
   values with median difference +0.00051161 mag. After removing that constant,
   maximum absolute residual is 0.00000545 mag. Photon weighting instead gives
   median -0.02143279 mag, maximum absolute 0.02705971 mag. Thus the published
   SED photometry is consistent with energy integration of the supplied response;
   this is a numerical finding, not an inference from its filename. The tiny
   constant has not been assigned to a particular physical constant or flux scale.
5. Adding the documented R=0.03 Vega offset to energy-weighted magnitudes makes
   their median difference from the photometric-library table +0.00001417 mag,
   but the 95th percentile absolute residual is 0.02046175 mag and maximum
   0.03236390 mag. A zero-point shift therefore does not make the full grid
   identical. This agrees qualitatively with the authors' distinction between
   photometric-library and spectral colour predictions:
   https://arxiv.org/html/1205.5496 (Section 3 and Figure 8).
6. Direct publisher masses versus the original table's ML_R and magnitude,
   using solar R=4.459, give median relative reconstruction error -0.00526%,
   maximum absolute 0.59650%. Allowing +/-0.0005 rounding in each tabulated
   ML_R and magnitude encloses 476/477 direct masses. The exception is
   Ech1.30Zm1.49T00.3000_iTp0.00_baseFe.fits: direct 0.7068 versus reconstructed
   0.704618 solar masses. This small published inconsistency is retained.
7. The SED table's mass column differs from out_mass_CH_BASTI for 27 SSPs at
   ages 0.06, 0.08, 0.09 Gyr, by up to about 1.19% when divided by the mass-file
   value. The discrepancy is larger than printed rounding. Do not silently mix
   these two mass sources: pin the adopted file and investigate version provenance.
8. The previous experiment reconstructed masses using adopted solar R=4.432816.
   Replacing that value by 4.459 increases those reconstructed masses by a common
   factor 1.0244095. This qualifies their absolute calibration; it is not a
   correction factor to apply blindly to existing Mass.py galaxy mass maps.

## Recommended calculation

Prefer a spectral calculation of SSP luminosities under the same explicitly
defined operator used for observed stellar light. Retain the proven nGIST LIGHT
normalization reconstruction. The preferred mass source is the direct matching
publisher M(*+remn) grid, with file/version provenance and the young-age
discrepancies documented, rather than recovering mass through an R magnitude.

For SSP j in bin b, let n_j be the actual nGIST template normalization divisor,
w_bj the fitted LIGHT coefficient, s_j(lambda) the original luminosity-density
spectrum per initial model mass, and m_j its surviving stars-plus-remnants mass.
Then a_bj=w_bj/n_j, up to an irrelevant common bin normalization. Define one
linear photometric operator P, with its exact response and wavelength weighting,
and apply it to both the SSP mixture and the observed stellar spectrum in a
consistent wavelength frame. Define L_b,model=P[sum_j a_bj s_j] and
M_b,model=sum_j a_bj m_j. The inferred mass is

    M_b = L_b,observed / L_b,model * M_b,model.

Equivalently compute ML_b,P=M_b,model/L_b,model and multiply by the matched
observed luminosity. If reporting solar-band units, compute the solar-band
luminosity using the same P; use it in both terms so it cancels in final mass.
AB/Vega zero points likewise cancel when used consistently. They do not cancel
if a tabulated M/L is held fixed while the luminosity convention is changed.

This can retain speclite's present bessell-R implementation for the cube and
apply precisely the same photon-weighted operator to the SSPs. It can instead
use another explicitly selected operator on both. It should not multiply the
current native-pPXF Johnson-Cousins_R M/L by speclite bessell-R luminosity.
The package name has no physical priority over consistency of the operator.

The full chain also requires consistent foreground/internal attenuation,
stellar-only flux, observed/rest wavelength treatment, distance, spatial
aperture, and AO gaps/masked spectral samples. Use the same effective operator
for model and data, or a documented model-based correction to a complete band.
The 7000-A fit can infer populations used to predict full R from full SSPs;
do not truncate the R integral merely because the fitting range ends at 7000 A.
Keep unconstrained multiplicative continuum polynomials distinct from physical
attenuation. Sum spatial masses for the total; integrated M/L is sum(M)/sum(L)
over the same footprint, not the mean or median resolved M/L.

## Scope and remaining uncertainty

These findings supersede the earlier lack of table-zero-point and direct-mass
provenance. They do not invalidate the earlier verified nGIST normalization
reproduction, and do not establish a universally more accurate R passband.
The original comparison is valid under its explicit adopted assumptions, but
its mass scale and differing luminosity definitions do not justify choosing a
production estimator by closeness of those M/L numbers alone.

Next implementation should be a separate end-to-end closure check using pinned
mass arrays and one common model/data/solar operator, including a clear decision
on the 27 young-age mass discrepancies. No production conversion was changed
in this research pass.
