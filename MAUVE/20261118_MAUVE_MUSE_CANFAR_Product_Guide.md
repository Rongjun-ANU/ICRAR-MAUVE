---
title: "MAUVE–MUSE products on CANFAR"
subtitle: "Archive locations and FITS data dictionary"
date: "18 November 2026 · MAUVE team reference"
lang: en-AU
---

# 1. Find the products

Open the project product directory at **`/arc/projects/mauve/products`** from a CANFAR terminal, notebook or CARTA session. From outside a session, use **`arc:projects/mauve/products`** with the CANFAR storage client or sign in with a CADC account that has MAUVE project access. The live ARC metadata marks the project product root as private and names the MAUVE group. [Browse the project in CANFAR](https://www.canfar.net/storage/arc/list/projects/mauve/products) · [CANFAR transfer documentation](https://www.opencadc.org/canfar/latest/platform/storage/transfers/)

| Analysis branch | Directory relative to `products/` | What it provides |
|---|---|---|
| Regular 8900 Å (`normal`) | `v3tk_v7.6.8/<GALID>/` | nGIST v7.6.8 and three downstream map files |
| Regular 7000 Å | `v3tk_v7.6.8_7000/<GALID>/` | nGIST v7.6.8 and three downstream map files |
| 7000 Å SFH16 | `v3tk_v7.6.8_7000/<GALID>_SFH16/` | SFH refit, refit provenance and a stellar-mass file. Use the sibling ordinary 7000 Å folder for gas and EW. |

There are **39 cube IDs representing 40 galaxies**. `NGC4567_8` identifies the shared cube for NGC 4567 and NGC 4568. The 18 November inventory includes **NGC4579 in all three branches**, superseding its absence in the earlier SFH16 transfer snapshot.

```text
/arc/projects/mauve/products/
├── v3tk_v7.6.8/
│   ├── <GALID>/                nGIST + mass, gas/SFR/O/H, and EW maps
│   └── {mass,sfr,proxy_ewha}_logs/
├── v3tk_v7.6.8_7000/
│   ├── <GALID>/                nGIST + mass, gas/SFR/O/H, and EW maps
│   ├── <GALID>_SFH16/          SFH refit + SFH16 stellar-mass maps
│   └── {mass,sfr,proxy_ewha}_logs/
├── v0.5/ and v0.6/             earlier product trees
└── configFiles/ and scripts/    earlier configurations and tooling
```

The earlier trees also contain ancillary products. For example, `v0.5/IC3392` and `v0.6/IC3392` contain moment products and quicklooks, and `v0.6/IC3392` contains multithreshold gas maps.

```bash
# In a CANFAR terminal
cd /arc/projects/mauve/products/v3tk_v7.6.8_7000
ls IC3392 IC3392_SFH16
```

For a single download from an authenticated external terminal, use `vcp arc:projects/mauve/products/<RUN>/<GALID>/<FILE> ./` with the literal run, cube ID and filename. CANFAR also documents `canfar data ls` and `canfar data cp` for clients with the data commands. [CANFAR data commands](https://www.opencadc.org/canfar/latest/cli/data/)

# 2. Original nGIST products

Each regular cube directory lists these ten non-empty items. Here `G` stands for the exact cube ID. Both regular runs use nGIST v7.6.8. The run labels refer to their nominal upper fitting limits, with 4800 Å as the configured lower limit in the checked configurations. Consult the per-run `CONFIG` for the settings of an individual module.

| File | Contents |
|---|---|
| `CONFIG` | Configuration saved by the run, including wavelength ranges, templates and module options. |
| `LOGFILE` | nGIST execution log, including module progress and run-completion messages. |
| `G_mask.fits` | Run footprint and exclusion mask used in the fit. |
| `G_spatial_binning_maps.fits` | Voronoi bin assignment (`BINID`), flux and S/N diagnostics, and bin-coordinate maps. |
| `G_kin_maps.fits` | Stellar line-of-sight kinematics and fit diagnostics. |
| `G_sfh_maps.fits` | Stellar-population summaries on the map grid, including age, metallicity and stellar reddening. |
| `G_sfh_weights.fits` | Fitted SSP template weights for the stellar population in each bin. Interpret them through the normalization specified by `CONFIG`. |
| `G_gas_bin_maps.fits` | Emission-line fluxes, errors, velocities and dispersion products from bin-level fits, on the spatial grid. |
| `G_gas_spaxel_maps.fits` | Emission-line maps from individual-spaxel fits. |
| `G_cont_cube.fits` | Fitted stellar-continuum spectral cube used to estimate the local Hα continuum for EW. |

## The 7000 Å SFH16 refit

The SFH16 branch reuses the ordinary 7000 Å binning, saved binned spectra and fixed stellar kinematics and refits the stellar populations with `SFH.MDEG = 16` in place of 12. The refit uses saved spectra with `SFH.SPEC_EMICLEAN = false`. The four nGIST items in every `G_SFH16/` directory are:

```text
CONFIG
LOGFILE
G_SFH16_sfh_maps.fits
G_SFH16_sfh_weights.fits
```

The adjacent `G_SFH16_spatial_binning_maps_further.fits` is the downstream mass product from the SFH16 populations and inherited 7000 Å binning. The mask, KIN maps, regular binned and spaxel gas maps, continuum cube, and downstream gas and EW maps are in the sibling `G/` directory. For a resolved analysis using **7000 Å SFH16 stellar masses**, take the mass file from `G_SFH16/` and the nebular and EW products from `G/`.

The four-item SFH16 archive layout follows the later storage-efficient transfer. Earlier transfer records describe an interim layout with extra copied inputs. [MAUVE Setonix-to-CANFAR transfer report](20261008%20Setonix%20to%20CANFAR%20Uploader%20Report.md)

# 3. Downstream products and logs

| File inside a regular `G/` folder | Use | 8900 Å | 7000 Å |
|---|---|---|---|
| `G_spatial_binning_maps_further.fits` | Stellar mass, mass surface density, M/L and R-band photometry | 39/39 | 39/39 |
| `G_gas_bin_maps_further.fits` | Gas dust correction, BPT, Hα SFR, nebular abundance and plasma diagnostics | 39/39 | 39/39 |
| `G_proxy_EW_maps_further.fits` | Continuum-based proxy EW(Hα) and R-band-based comparison | 39/39 | 39/39 |

All counts refer to non-empty archive objects. All 39 SFH16 directories have the four original SFH16 items and one downstream mass file. Together, the two regular runs and SFH16 hold 273 downstream FITS files.

| Path from each regular run root | Stage and useful provenance |
|---|---|
| `mass_logs/G.log` | `Mass.sh` / `Mass.py`. Input cube and SFH files, M/L normalization, geometry and mass diagnostics. |
| `sfr_logs/G.log` | `SFR.sh` / `SFR+Z.py`. Detection/classification, abundance-method diagnostics and output path. |
| `proxy_ewha_logs/G.log` | `proxy_EWHa.sh` / `proxy_EWHa.py`. Input gas, photometry and continuum paths, EW window and valid-pixel counts. |

Each of these six run-level log directories lists 39 files. All 234 regular wrapper logs contained a success marker on read-back. The per-run nGIST `LOGFILE` is available next to the SFH16 mass product in the shared tree.

## Reading the maps

Select HDUs by `EXTNAME`. Section 7 lists each HDU by name. Check the FITS header for `BUNIT`, comments and any method provenance. A header value of `1` is used for many log, class and flag images in the gas product. Their physical meanings are specified below. A fitted bin can appear on many spatial pixels, so repeated bin-level quantities should be treated as correlated in resolved statistics. For summing extensive images, use the specific image definition and pixel mask.

# 4. Stellar mass and R-band photometry

**File:** `G_spatial_binning_maps_further.fits`, or `G_SFH16_spatial_binning_maps_further.fits` for the SFH16 fit. The inherited nGIST layers in the checked SFH16 file are the primary HDU, `BINID`, `FLUX`, `SNR`, `SNRBIN`, `XBIN` and `YBIN`. Mass appends 11 images.

| Added HDU | Quantity and units |
|---|---|
| `FLUX_R_corr`, `FLUX_R_corr_ERR` | Internally attenuation-corrected synthetic R-band flux and propagated 1σ uncertainty, in nanomaggies. |
| `ML_R` | Stellar-population R-band mass-to-light ratio, in solar masses per solar R-band luminosity. A bin-level result mapped to pixels. |
| `LOGMSTAR`, `LOGMSTAR_ERR` | Base-10 log stellar mass assigned per map pixel (solar masses) and its 1σ error in dex. |
| `LOGMASS_SURFACE_DENSITY`, `LOGMASS_SURFACE_DENSITY_ERR` | Base-10 log inclination-adjusted stellar mass surface density (solar masses per kpc²) and its 1σ error in dex. |
| `MAGNITUDE_R`, `MAGNITUDE_R_ERR` | Stellar-attenuation-corrected synthetic R apparent AB magnitude and 1σ error. |
| `MAGNITUDE_R_UNCORRECTED`, `MAGNITUDE_R_UNCORRECTED_ERR` | Synthetic R apparent AB magnitude before the pipeline's internal stellar attenuation correction, and 1σ error. This is the input to the legacy EW denominator. |

Mass uses synthetic R photometry from the observed cube, assigns bin-averaged flux per pixel and corrects the magnitude for attenuation estimated from **stellar** SFH reddening. Luminosities and physical pixel area use a fixed distance of 16.5 Mpc. The current default `light_norm_to_r` converts the fitted LIGHT coefficients to R-band light fractions before combining template M/L values. Check `MLMODE` and the run log for the adopted convention on an archived file.

$$
M_\star=(M/L)_R L_{R,\mathrm{corr}},\qquad
\Sigma_\star=\frac{M_\star}{A_{\mathrm{pix}}}\,f_i,\qquad
f_i=\sqrt{(1-q_0^2)\cos^2i+q_0^2},\quad q_0=0.2 .
$$

Here $A_{\mathrm{pix}}$ is the projected physical area from the map WCS and $i$ is the catalogue inclination. The pipeline uses $f_i=1$ if no applicable inclination is supplied. The formal mass errors propagate photometry and any available stellar attenuation error at fixed M/L, distance and geometry. For separate component measurements in the combined cube, use a component mask and establish which inclination correction applies to each component.

# 5. Binned gas maps

**File:** `G_gas_bin_maps_further.fits`. Its first HDUs are the inherited nGIST binned gas-map content, including `BIN_ID`, raw fluxes and errors, velocities, dispersions and fit diagnostics. Downstream HDUs are appended. The 7000 Å schema contains 371 HDUs, including 91 inherited HDUs. The 8900 Å schema contains 423 HDUs, including 119 inherited HDUs. Four additional [O II] lines account for the difference, with seven inherited and six attenuation-corrected extensions per line. Section 7 lists the individual names.

## 5.1 Dust correction, BPT and SFR

| Added HDU or pattern | Quantity |
|---|---|
| `GAS_E_BV_BD`, `GAS_E_BV_BD_ERR` | Nebular colour excess from Hα/Hβ and first-order 1σ error, in mag. |
| `<LINE>_FLUX_CORR`, `<LINE>_FLUX_ERR_CORR` | Attenuation-corrected line flux and propagated error, in 10⁻²⁰ erg s⁻¹ cm⁻². Per-line versions with `_HII` and `_SF` apply the corresponding region masks. |
| `NII_BPT`, `SII_BPT` | Independent integer maps for [N II] and [S II] BPT classifications from corrected ratios. |
| `HALPHA_LUMINOSITY_CORR`, `HALPHA_SFR_CORR` | Hα luminosity (erg s⁻¹) and Hα-based equivalent SFR per map pixel (solar masses yr⁻¹). These full maps include substitutions where the Balmer lines fail the detection criteria. |
| `LOGSFR_SURFACE_DENSITY` | Base-10 log of inclination-adjusted Hα SFR surface density (solar masses yr⁻¹ kpc⁻²). |
| `LOGSFR_SURFACE_DENSITY_<REGION>` | Log surface density for `HII`, `SF`, `NONHII`, `NONSF`, `UNCLASSIFIED2`, `UNCLASSIFIED1`, `UPPER`. |
| `HALPHA_SFR_CORR_<REGION>` | Linear surface density in solar masses yr⁻¹ kpc⁻² for the six named regions excluding `UPPER`. |

The gas uses an intrinsic Hα/Hβ ratio of 2.86 and the default Calzetti curve. The decrement is floored at 2.86. With $k(\lambda)=A_\lambda/E(B-V)$, the corrected flux is $F_{\lambda,\mathrm{corr}}=F_\lambda 10^{0.4 k(\lambda) E(B-V)}$. For SFR, $L_{\mathrm{H}\alpha}=4\pi D^2F_{\mathrm{H}\alpha}$ with $D=16.5$ Mpc. The Chabrier conversion is $\mathrm{SFR}=4.98\times10^{-42}L_{\mathrm{H}\alpha}$ in solar masses yr⁻¹ for luminosity in erg s⁻¹. The surface densities use the inclination factor in Section 4.

**BPT codes.** Both maps use −1 for unknown/non-detection or invalid diagnostic and 0 for an uncertain class with usable inputs. Positive classes require stable central and axis-wise ±1σ positions. The [N II] codes 1, 2, 3 mean H II, composite and AGN. The [S II] codes 1, 2, 3 mean H II, LINER and Seyfert.

**Region suffixes.** The default selection uses both diagrams. `_SF` (classification 1) allows H II or composite in [N II] and requires H II in [S II]. `_HII` (classification 2) requires H II in both. The corresponding `NON*` layers hold classified emission outside the respective selection. `UNCLASSIFIED1/2` carry the ambiguous regions for the respective scheme. These pipeline masks have their own detection and uncertainty logic, including treatment of missing non-Balmer lines. A downstream science sample may apply further EW or kinematic cuts.

**Balmer non-detections.** The main SFR pathway tests both Balmer lines against S/N ≥ 3 and a flux floor of 20 in the stored flux units. When either line fails these criteria, it substitutes the larger of the observed Hα and that floor where finite. `_UPPER` carries those proxy values on the upper-limit mask. These are pipeline flux-floor proxies. Use a calibrated censoring model if statistical confidence bounds are needed. Outside the selected regions, continuous masked images generally carry NaN.

## 5.2 Nebular abundances, density and temperature

All oxygen-abundance images below are in $12+\log_{10}(\mathrm{O/H})$. The two region suffixes `_HII` and `_SF` have the meanings in Section 5.1. Errors on log abundances are in dex. Each diagnostic has its own line and validity requirements, so finite coverage can differ by method, galaxy and wavelength range.

### Strong-line abundance calibrations

| HDU stem (append `_HII` or `_SF`) | Calibration |
|---|---|
| `O_H_D16` | Dopita et al. (2016) [N II], [S II] and Hα calibration. |
| `O_H_PG16` | Pilyugin & Grebel (2016) S calibration using oxygen, nitrogen and sulphur line ratios. |
| `O_H_N2S2_N06` | Cubic calibration of log([N II] λ6583 / summed [S II] doublet), labelled N2S2-N06 in the pipeline. |
| `O_H_O3N2_M13`, `O_H_N2_M13` | Marino et al. (2013) O3N2 and N2 calibrations. |
| `O_H_O3N2_PP04`, `O_H_N2_PP04` | Pettini & Pagel (2004) O3N2 and N2 calibrations. |
| `O_H_{O3N2,O3S2,RS32,R3,N2,S2}_C20` | Six Curti et al. (2020) line-ratio branches. |
| `O_H_COMBINED_C20` | Inverse-variance combination of usable C20 estimates. |

The C20 value maps have paired `_HII_ERR` and `_SF_ERR` layers. `COMBINED_C20_METHOD` identifies the dominant weight of the combination (−1 none, 0 O3N2, 1 O3S2, 2 RS32, 3 R3, 4 N2, 5 S2). The index identifies a contributor to the combined estimate. Shared diagnostic lines can induce correlated errors.

### Electron density and temperature-based methods

| HDU stem or family | Quantity and coverage |
|---|---|
| `NE_SII`, `NE_SII_ERR`, `NE_SII_ALL` | PyNeb density from [S II] λ6716/λ6730 and propagated error (cm⁻³). `NE_SII_ALL` fills usable pixels lacking a density estimate with 20 cm⁻³. |
| `NE_SII_FIXED20`, `NE_SII_EXACT_HIDEN` | Density floor and exact high-density solver fallback flags. |
| `TE_NII`, `TE_SIII_BR24` | Temperatures (K) from [N II] and [S III] auroral diagnostics. The [S III] branch requires λ9069 to infer the λ9532 contribution. |
| `TE_OIII_BR24`, `TE_OIII_NII_CHAIN` | High-ionisation temperature inferred from the [S III] or [N II] temperature chain, in K. |
| `O_PLUS_BR24`, `O_DPLUS_BR24`, `O_H_BR24_DIRECT` | Brazzini-style multi-zone O⁺/H⁺ and O²⁺/H⁺ ionic ratios and their summed oxygen abundance. |
| `O_PLUS_NII_OII7325`, `O_DPLUS_NII_OII7325`, `O_H_NII_OII7325` | Semi-direct ionic abundances using [N II] temperature and [O II] λ7325, plus the summed oxygen abundance. |
| `O_H_NII_K25` | Kreckel-style [N II] temperature-based abundance proxy. |
| `BR24_DIRECT_VALID`, `NII_OII7325_VALID`, `NII_K25_VALID` | Method-specific flags for passing the input gates. |

Temperature and ionic/oxygen abundance stems above have paired region and corresponding `_ERR` layers. Validity flags also carry the region suffix. For ionic sums the stored order is $12+\log_{10}(\mathrm{O^+/H^+}+\mathrm{O^{2+}/H^+})$. The 7000 Å fit excludes the red [O II] λλ7320,7330 and [S III] λ9069 diagnostics. Check the original fitted line list and the method's valid map even when its HDU exists. The 8900 Å fit adds the [O II] wavelength region, which can support the semi-direct branch. [S III] λ9069 lies beyond its nominal upper bound.

## 5.3 Photoionisation-grid diagnostics

The current gas writer includes three modelling families. Map availability and finite coverage depend on the archive version, input lines and the fit flags. These grids yield model-dependent gas abundances and ionisation parameters. Use one documented grid and selection consistently when comparing resolved measurements.

| Model and region-paired HDU stems | Interpretation |
|---|---|
| **pyqz** · `O_H_PYQZ`, `LOG_Q_PYQZ`, `LOG_U_PYQZ` | Combined multivariate-KDE estimates of gas abundance and ionisation. $q$ is in cm s⁻¹ and $U=q/c$ is dimensionless. The region+`_ERR` layers report a KDE-contour half-extent. |
| **NebulaBayes** · `O_H_NEBULABAYES`, `LOG_U_NEBULABAYES`, `LOG_PK_NEBULABAYES` | Marginal posterior means for abundance, log ionisation parameter and log pressure $P/k$ in K cm⁻³. Region+`_CI68_LOW/HIGH` give equal-tailed credible bounds and region+`_MODE` gives the sampled marginal mode. |
| **JY22** · `O_H_JY22`, `O_H_JY22_PREDEP`, `LOG_U_JY22` | Posterior means from the Ji & Yan-based grid and the joint N2, S2, R3 log-ratio likelihood with propagated covariance. The first O/H is post-depletion gas-phase abundance. `PREDEP` is the native pre-depletion grid scale. Region+`_16` and `_84` give marginal percentiles. |

### Fit diagnostics

- **pyqz** · `PYQZ_FLAG` is the package flag, `PYQZ_RS_OFFGRID` gives the percentage of Monte Carlo draws off-grid, and `PYQZ_VALID` identifies finite O/H and log $q$ results.
- **NebulaBayes** · `NB_CHI2_RED` records the package fit statistic at marginal modes, `NB_NLOCALMAX` the maximum number of local maxima among marginal PDFs, and `NB_VALID` finite summaries. `NB_FLAG` bits 1, 2 and 4 mark edge solutions in O/H, log $U$ and log $P/k$, respectively. Bits 8 and 16 mark open/non-finite intervals and fit exceptions. `NB_OH_RELIABLE`, `NB_LOGU_RELIABLE` and `NB_LOGPK_RELIABLE` apply per-parameter reliability checks.
- **JY22** · `JY22_CHI2_MIN` gives the minimum fit statistic, and `JY22_RESID_{N2,S2,R3,NORM}` gives observed-minus-model log-ratio residuals and their Euclidean norm at the best grid point. `JY22_VALID` records finite summaries. `JY22_FIT_OK` adds the empirical χ² ≤ 25 fit cut. `JY22_FLAG` bits 1, 2, 4, 8, 16 and 32 mark O/H edge, log $U$ edge, invalid posterior, invalid covariance, exception and gross mismatch, respectively.

The diagnostic stems take `_HII` or `_SF`. In categorical model diagnostics, −99 denotes a location where the method was unevaluated when specified in the header. For a bitmask, multiple bits can be set at once. For a finite abundance map, check the applicable diagnostics alongside the region mask.

# 6. Hα equivalent-width proxies

**File:** `G_proxy_EW_maps_further.fits`. This is a standalone seven-HDU file assembled from the observed binned Hα gas flux, the mass-file photometry and the nGIST continuum model. Its primary HDU carries the Hα image under `OBS_HA6562_FLUX`. The output uses bin-level values on the map grid.

| HDU | Quantity |
|---|---|
| `OBS_HA6562_FLUX` (primary) | Observed Hα flux from the binned gas product, reduced by `BINID` for EW. Check the input flux scaling in `BUNIT`. |
| `CONT_HA_MEAN` | Mean model continuum flux density near Hα in erg s⁻¹ cm⁻² Å⁻¹, in the bin rest-frame convention. |
| `CONT_HA_NPIX` | Count of finite continuum-cube samples contributing to each bin, presented on the map grid. |
| `OBS_R_FLUX` | Observed synthetic R flux derived from `MAGNITUDE_R_UNCORRECTED`, in nanomaggies. |
| `CONT_HA_RBAND` | Legacy continuum estimate obtained by converting that R-band flux density via $f_\lambda=f_\nu c/\lambda^2$ at 6562.8 Å. |
| `PROXY_EWHA` | Observed Hα flux divided by `CONT_HA_MEAN`, in Å. |
| `PSEUDO_EWHA` | Observed Hα flux divided by `CONT_HA_RBAND`, in Å. Legacy comparison. |

The preferred continuum-based proxy uses a passband of 6557.6–6571.6 Å in rest vacuum wavelengths, converted to air when applicable. The current algorithm shifts the window by the systemic redshift and fitted bin Hα velocity, and scales the continuum flux density by $1+z_{\mathrm{bin}}$ to form the rest-frame EW. `CONT_HA_NPIX` and the primary/HDU provenance keywords indicate continuum support and adopted inputs. Check the headers for older products, as implementation details can vary across archived runs.

```python
from astropy.io import fits

root = "/arc/projects/mauve/products/v3tk_v7.6.8_7000"
g = "IC3392"
with fits.open(f"{root}/{g}_SFH16/{g}_SFH16_spatial_binning_maps_further.fits") as f:
    log_sigma_star = f["LOGMASS_SURFACE_DENSITY"].data.copy()
with fits.open(f"{root}/{g}/{g}_gas_bin_maps_further.fits") as f:
    log_sigma_sfr_hii = f["LOGSFR_SURFACE_DENSITY_HII"].data.copy()
with fits.open(f"{root}/{g}/{g}_proxy_EW_maps_further.fits") as f:
    ew_ha = f["PROXY_EWHA"].data.copy()
```

## Provenance and verification

The directory and coverage counts come from an authenticated, read-only CANFAR ARC VOSpace listing on 18 November 2026. The listing establishes canonical object presence and non-zero size. The data dictionary follows `further/Mass.py`, `further/SFR+Z.py`, `further/proxy_EWHa.py` and their shell wrappers. The SFH16 routing follows the dedicated uploader and the live directory listing. The checks included archived IC3392 configurations, mass and EW headers and 117 local 7000 Å FITS files whose lengths match the remote listing. Those local files had 371 HDUs in each further gas file (starting with 91 nGIST HDUs including primary), seven in each EW file, and 18 in each SFH16 mass file. The inventory is a file-level snapshot. The 7000 Å gas index was checked against the complete header sequences of all 39 local gas files. For 8900 Å, all 39 archived SFR logs list the same four additional fitted lines, and the original and downstream file sizes agree with the 28 inherited and 24 corrected [O II] extensions. A direct header read of the archived NGC4607 gas file matched the index through the inherited block and the first corrected-line HDUs. The rest of the 8900 Å index follows the output writer, whose extension names and order were checked against the complete 7000 Å headers.

# 7. FITS HDU index

Names below are printed in uppercase for lookup. FITS name lookup in Astropy is case-insensitive. `PRIMARY` denotes the unnamed HDU at index 0 in the mass and gas files. In the EW file, HDU 0 is named `OBS_HA6562_FLUX`. The gas lists specify names present in the schema. Pixel coverage depends on the diagnostic and its input lines.

## 7.1 Stellar mass, photometry and EW

**Stellar mass and photometry** · 18 HDUs for the regular and SFH16 mass files. The first seven, `PRIMARY` through `YBIN`, are inherited from the nGIST spatial-binning file. Section 4 defines the 11 additions.

```text
PRIMARY                               YBIN                                  LOGMASS_SURFACE_DENSITY
BINID                                 FLUX_R_CORR                           LOGMASS_SURFACE_DENSITY_ERR
FLUX                                  FLUX_R_CORR_ERR                       MAGNITUDE_R
SNR                                   ML_R                                  MAGNITUDE_R_ERR
SNRBIN                                LOGMSTAR                              MAGNITUDE_R_UNCORRECTED
XBIN                                  LOGMSTAR_ERR                          MAGNITUDE_R_UNCORRECTED_ERR
```

**Hα equivalent-width proxies** · Seven HDUs. Section 6 defines the primary flux map and its six companion images.

```text
OBS_HA6562_FLUX                       OBS_R_FLUX                            PSEUDO_EWHA
CONT_HA_MEAN                          CONT_HA_RBAND
CONT_HA_NPIX                          PROXY_EWHA
```

## 7.2 Binned gas maps · Inherited nGIST HDUs

**Both wavelength runs** · The 91 HDUs below carry the original primary, bin assignment, stellar and fitting diagnostics, emission-line maps, and original BPT maps. Each line has fitted flux, velocity and dispersion, their errors, and the dispersion correction term. The inherited BPT maps are `NII_HA_BPT` and `SII_HA_BPT`. Section 5.1 defines the added `NII_BPT` and `SII_BPT` maps.

```text
PRIMARY                               NII5754_SIGMA_ERR                     HA6562_FLUX_ERR
BIN_ID                                NII5754_SIGMA_CORR                    HA6562_VEL
V_STARS2                              OI6300_FLUX                           HA6562_VEL_ERR
SIGMA_STARS2                          OI6300_FLUX_ERR                       HA6562_SIGMA
CHI2_TOT                              OI6300_VEL                            HA6562_SIGMA_ERR
HB4861_FLUX                           OI6300_VEL_ERR                        HA6562_SIGMA_CORR
HB4861_FLUX_ERR                       OI6300_SIGMA                          NII6583_FLUX
HB4861_VEL                            OI6300_SIGMA_ERR                      NII6583_FLUX_ERR
HB4861_VEL_ERR                        OI6300_SIGMA_CORR                     NII6583_VEL
HB4861_SIGMA                          SIII6312_FLUX                         NII6583_VEL_ERR
HB4861_SIGMA_ERR                      SIII6312_FLUX_ERR                     NII6583_SIGMA
HB4861_SIGMA_CORR                     SIII6312_VEL                          NII6583_SIGMA_ERR
OIII4958_FLUX                         SIII6312_VEL_ERR                      NII6583_SIGMA_CORR
OIII4958_FLUX_ERR                     SIII6312_SIGMA                        SII6716_FLUX
OIII4958_VEL                          SIII6312_SIGMA_ERR                    SII6716_FLUX_ERR
OIII4958_VEL_ERR                      SIII6312_SIGMA_CORR                   SII6716_VEL
OIII4958_SIGMA                        OI6363_FLUX                           SII6716_VEL_ERR
OIII4958_SIGMA_ERR                    OI6363_FLUX_ERR                       SII6716_SIGMA
OIII4958_SIGMA_CORR                   OI6363_VEL                            SII6716_SIGMA_ERR
OIII5006_FLUX                         OI6363_VEL_ERR                        SII6716_SIGMA_CORR
OIII5006_FLUX_ERR                     OI6363_SIGMA                          SII6730_FLUX
OIII5006_VEL                          OI6363_SIGMA_ERR                      SII6730_FLUX_ERR
OIII5006_VEL_ERR                      OI6363_SIGMA_CORR                     SII6730_VEL
OIII5006_SIGMA                        NII6548_FLUX                          SII6730_VEL_ERR
OIII5006_SIGMA_ERR                    NII6548_FLUX_ERR                      SII6730_SIGMA
OIII5006_SIGMA_CORR                   NII6548_VEL                           SII6730_SIGMA_ERR
NII5754_FLUX                          NII6548_VEL_ERR                       SII6730_SIGMA_CORR
NII5754_FLUX_ERR                      NII6548_SIGMA                         NII_HA_BPT
NII5754_VEL                           NII6548_SIGMA_ERR                     SII_HA_BPT
NII5754_VEL_ERR                       NII6548_SIGMA_CORR
NII5754_SIGMA                         HA6562_FLUX
```

**8900 Å additions to the inherited block** · 28 HDUs for four [O II] components.

```text
OII7318_FLUX                          OII7319_VEL_ERR                       OII7329_SIGMA_CORR
OII7318_FLUX_ERR                      OII7319_SIGMA                         OII7330_FLUX
OII7318_VEL                           OII7319_SIGMA_ERR                     OII7330_FLUX_ERR
OII7318_VEL_ERR                       OII7319_SIGMA_CORR                    OII7330_VEL
OII7318_SIGMA                         OII7329_FLUX                          OII7330_VEL_ERR
OII7318_SIGMA_ERR                     OII7329_FLUX_ERR                      OII7330_SIGMA
OII7318_SIGMA_CORR                    OII7329_VEL                           OII7330_SIGMA_ERR
OII7319_FLUX                          OII7329_VEL_ERR                       OII7330_SIGMA_CORR
OII7319_FLUX_ERR                      OII7329_SIGMA
OII7319_VEL                           OII7329_SIGMA_ERR
```

## 7.3 Binned gas maps · Attenuation-corrected line HDUs

**Both wavelength runs** · 72 HDUs. The suffixes identify corrected flux, its uncertainty, and the two region selections defined in Section 5.1.

```text
HB4861_FLUX_CORR                      OI6300_FLUX_CORR                      HA6562_FLUX_CORR
HB4861_FLUX_ERR_CORR                  OI6300_FLUX_ERR_CORR                  HA6562_FLUX_ERR_CORR
HB4861_FLUX_CORR_HII                  OI6300_FLUX_CORR_HII                  HA6562_FLUX_CORR_HII
HB4861_FLUX_CORR_SF                   OI6300_FLUX_CORR_SF                   HA6562_FLUX_CORR_SF
HB4861_FLUX_ERR_CORR_HII              OI6300_FLUX_ERR_CORR_HII              HA6562_FLUX_ERR_CORR_HII
HB4861_FLUX_ERR_CORR_SF               OI6300_FLUX_ERR_CORR_SF               HA6562_FLUX_ERR_CORR_SF
OIII4958_FLUX_CORR                    SIII6312_FLUX_CORR                    NII6583_FLUX_CORR
OIII4958_FLUX_ERR_CORR                SIII6312_FLUX_ERR_CORR                NII6583_FLUX_ERR_CORR
OIII4958_FLUX_CORR_HII                SIII6312_FLUX_CORR_HII                NII6583_FLUX_CORR_HII
OIII4958_FLUX_CORR_SF                 SIII6312_FLUX_CORR_SF                 NII6583_FLUX_CORR_SF
OIII4958_FLUX_ERR_CORR_HII            SIII6312_FLUX_ERR_CORR_HII            NII6583_FLUX_ERR_CORR_HII
OIII4958_FLUX_ERR_CORR_SF             SIII6312_FLUX_ERR_CORR_SF             NII6583_FLUX_ERR_CORR_SF
OIII5006_FLUX_CORR                    OI6363_FLUX_CORR                      SII6716_FLUX_CORR
OIII5006_FLUX_ERR_CORR                OI6363_FLUX_ERR_CORR                  SII6716_FLUX_ERR_CORR
OIII5006_FLUX_CORR_HII                OI6363_FLUX_CORR_HII                  SII6716_FLUX_CORR_HII
OIII5006_FLUX_CORR_SF                 OI6363_FLUX_CORR_SF                   SII6716_FLUX_CORR_SF
OIII5006_FLUX_ERR_CORR_HII            OI6363_FLUX_ERR_CORR_HII              SII6716_FLUX_ERR_CORR_HII
OIII5006_FLUX_ERR_CORR_SF             OI6363_FLUX_ERR_CORR_SF               SII6716_FLUX_ERR_CORR_SF
NII5754_FLUX_CORR                     NII6548_FLUX_CORR                     SII6730_FLUX_CORR
NII5754_FLUX_ERR_CORR                 NII6548_FLUX_ERR_CORR                 SII6730_FLUX_ERR_CORR
NII5754_FLUX_CORR_HII                 NII6548_FLUX_CORR_HII                 SII6730_FLUX_CORR_HII
NII5754_FLUX_CORR_SF                  NII6548_FLUX_CORR_SF                  SII6730_FLUX_CORR_SF
NII5754_FLUX_ERR_CORR_HII             NII6548_FLUX_ERR_CORR_HII             SII6730_FLUX_ERR_CORR_HII
NII5754_FLUX_ERR_CORR_SF              NII6548_FLUX_ERR_CORR_SF              SII6730_FLUX_ERR_CORR_SF
```

**8900 Å additions** · 24 HDUs in the same six forms for each [O II] component.

```text
OII7318_FLUX_CORR                     OII7319_FLUX_CORR_HII                 OII7329_FLUX_ERR_CORR_HII
OII7318_FLUX_ERR_CORR                 OII7319_FLUX_CORR_SF                  OII7329_FLUX_ERR_CORR_SF
OII7318_FLUX_CORR_HII                 OII7319_FLUX_ERR_CORR_HII             OII7330_FLUX_CORR
OII7318_FLUX_CORR_SF                  OII7319_FLUX_ERR_CORR_SF              OII7330_FLUX_ERR_CORR
OII7318_FLUX_ERR_CORR_HII             OII7329_FLUX_CORR                     OII7330_FLUX_CORR_HII
OII7318_FLUX_ERR_CORR_SF              OII7329_FLUX_ERR_CORR                 OII7330_FLUX_CORR_SF
OII7319_FLUX_CORR                     OII7329_FLUX_CORR_HII                 OII7330_FLUX_ERR_CORR_HII
OII7319_FLUX_ERR_CORR                 OII7329_FLUX_CORR_SF                  OII7330_FLUX_ERR_CORR_SF
```

## 7.4 Binned gas maps · Extinction, SFR and abundance HDUs

**Gas colour excess, BPT, Hα luminosity and SFR** · 20 HDUs. See Section 5.1.

```text
GAS_E_BV_BD                           LOGSFR_SURFACE_DENSITY_HII            HALPHA_SFR_CORR_HII
GAS_E_BV_BD_ERR                       LOGSFR_SURFACE_DENSITY_SF             HALPHA_SFR_CORR_SF
NII_BPT                               LOGSFR_SURFACE_DENSITY_NONHII         HALPHA_SFR_CORR_NONHII
SII_BPT                               LOGSFR_SURFACE_DENSITY_NONSF          HALPHA_SFR_CORR_NONSF
HALPHA_LUMINOSITY_CORR                LOGSFR_SURFACE_DENSITY_UNCLASSIFIED2  HALPHA_SFR_CORR_UNCLASSIFIED2
HALPHA_SFR_CORR                       LOGSFR_SURFACE_DENSITY_UNCLASSIFIED1  HALPHA_SFR_CORR_UNCLASSIFIED1
LOGSFR_SURFACE_DENSITY                LOGSFR_SURFACE_DENSITY_UPPER
```

**Strong-line oxygen abundance** · 43 HDUs. See Section 5.2.

```text
O_H_D16_HII                           O_H_O3N2_C20_SF                       O_H_O3S2_C20_HII_ERR
O_H_D16_SF                            O_H_O3S2_C20_HII                      O_H_O3S2_C20_SF_ERR
O_H_PG16_HII                          O_H_O3S2_C20_SF                       O_H_RS32_C20_HII_ERR
O_H_PG16_SF                           O_H_RS32_C20_HII                      O_H_RS32_C20_SF_ERR
O_H_N2S2_N06_HII                      O_H_RS32_C20_SF                       O_H_R3_C20_HII_ERR
O_H_N2S2_N06_SF                       O_H_R3_C20_HII                        O_H_R3_C20_SF_ERR
O_H_O3N2_M13_HII                      O_H_R3_C20_SF                         O_H_N2_C20_HII_ERR
O_H_O3N2_M13_SF                       O_H_N2_C20_HII                        O_H_N2_C20_SF_ERR
O_H_N2_M13_HII                        O_H_N2_C20_SF                         O_H_S2_C20_HII_ERR
O_H_N2_M13_SF                         O_H_S2_C20_HII                        O_H_S2_C20_SF_ERR
O_H_O3N2_PP04_HII                     O_H_S2_C20_SF                         O_H_COMBINED_C20_HII_ERR
O_H_O3N2_PP04_SF                      O_H_COMBINED_C20_HII                  O_H_COMBINED_C20_SF_ERR
O_H_N2_PP04_HII                       O_H_COMBINED_C20_SF                   COMBINED_C20_METHOD
O_H_N2_PP04_SF                        O_H_O3N2_C20_HII_ERR
O_H_O3N2_C20_HII                      O_H_O3N2_C20_SF_ERR
```

**Electron density, temperatures, ionic abundances and temperature-based oxygen abundances** · 55 HDUs. See Section 5.2.

```text
NE_SII                                O_PLUS_BR24_HII                       O_DPLUS_NII_OII7325_SF_ERR
NE_SII_ALL                            O_PLUS_BR24_SF                        O_H_NII_OII7325_HII
NE_SII_ERR                            O_PLUS_BR24_HII_ERR                   O_H_NII_OII7325_SF
TE_NII_HII                            O_PLUS_BR24_SF_ERR                    O_H_NII_OII7325_HII_ERR
TE_NII_SF                             O_DPLUS_BR24_HII                      O_H_NII_OII7325_SF_ERR
TE_NII_HII_ERR                        O_DPLUS_BR24_SF                       O_H_NII_K25_HII
TE_NII_SF_ERR                         O_DPLUS_BR24_HII_ERR                  O_H_NII_K25_SF
TE_SIII_BR24_HII                      O_DPLUS_BR24_SF_ERR                   O_H_NII_K25_HII_ERR
TE_SIII_BR24_SF                       O_H_BR24_DIRECT_HII                   O_H_NII_K25_SF_ERR
TE_SIII_BR24_HII_ERR                  O_H_BR24_DIRECT_SF                    NE_SII_FIXED20
TE_SIII_BR24_SF_ERR                   O_H_BR24_DIRECT_HII_ERR               NE_SII_EXACT_HIDEN
TE_OIII_BR24_HII                      O_H_BR24_DIRECT_SF_ERR                BR24_DIRECT_VALID_HII
TE_OIII_BR24_SF                       O_PLUS_NII_OII7325_HII                BR24_DIRECT_VALID_SF
TE_OIII_BR24_HII_ERR                  O_PLUS_NII_OII7325_SF                 NII_OII7325_VALID_HII
TE_OIII_BR24_SF_ERR                   O_PLUS_NII_OII7325_HII_ERR            NII_OII7325_VALID_SF
TE_OIII_NII_CHAIN_HII                 O_PLUS_NII_OII7325_SF_ERR             NII_K25_VALID_HII
TE_OIII_NII_CHAIN_SF                  O_DPLUS_NII_OII7325_HII               NII_K25_VALID_SF
TE_OIII_NII_CHAIN_HII_ERR             O_DPLUS_NII_OII7325_SF
TE_OIII_NII_CHAIN_SF_ERR              O_DPLUS_NII_OII7325_HII_ERR
```

## 7.5 Binned gas maps · Photoionisation-grid HDUs

The pyqz, NebulaBayes and JY22 families add 90 HDUs. Section 5.3 defines their parameter estimates, posterior summaries and diagnostic flags.

```text
O_H_PYQZ_HII                          LOG_U_NEBULABAYES_HII_CI68_HIGH       O_H_JY22_HII_84
O_H_PYQZ_SF                           LOG_U_NEBULABAYES_SF_CI68_HIGH        O_H_JY22_SF_84
O_H_PYQZ_HII_ERR                      LOG_U_NEBULABAYES_HII_MODE            O_H_JY22_PREDEP_HII
O_H_PYQZ_SF_ERR                       LOG_U_NEBULABAYES_SF_MODE             O_H_JY22_PREDEP_SF
LOG_Q_PYQZ_HII                        LOG_PK_NEBULABAYES_HII                O_H_JY22_PREDEP_HII_16
LOG_Q_PYQZ_SF                         LOG_PK_NEBULABAYES_SF                 O_H_JY22_PREDEP_SF_16
LOG_Q_PYQZ_HII_ERR                    LOG_PK_NEBULABAYES_HII_CI68_LOW       O_H_JY22_PREDEP_HII_84
LOG_Q_PYQZ_SF_ERR                     LOG_PK_NEBULABAYES_SF_CI68_LOW        O_H_JY22_PREDEP_SF_84
LOG_U_PYQZ_HII                        LOG_PK_NEBULABAYES_HII_CI68_HIGH      LOG_U_JY22_HII
LOG_U_PYQZ_SF                         LOG_PK_NEBULABAYES_SF_CI68_HIGH       LOG_U_JY22_SF
LOG_U_PYQZ_HII_ERR                    LOG_PK_NEBULABAYES_HII_MODE           LOG_U_JY22_HII_16
LOG_U_PYQZ_SF_ERR                     LOG_PK_NEBULABAYES_SF_MODE            LOG_U_JY22_SF_16
PYQZ_FLAG_HII                         NB_CHI2_RED_HII                       LOG_U_JY22_HII_84
PYQZ_FLAG_SF                          NB_CHI2_RED_SF                        LOG_U_JY22_SF_84
PYQZ_RS_OFFGRID_HII                   NB_NLOCALMAX_HII                      JY22_CHI2_MIN_HII
PYQZ_RS_OFFGRID_SF                    NB_NLOCALMAX_SF                       JY22_CHI2_MIN_SF
PYQZ_VALID_HII                        NB_FLAG_HII                           JY22_RESID_N2_HII
PYQZ_VALID_SF                         NB_FLAG_SF                            JY22_RESID_N2_SF
O_H_NEBULABAYES_HII                   NB_OH_RELIABLE_HII                    JY22_RESID_S2_HII
O_H_NEBULABAYES_SF                    NB_OH_RELIABLE_SF                     JY22_RESID_S2_SF
O_H_NEBULABAYES_HII_CI68_LOW          NB_LOGU_RELIABLE_HII                  JY22_RESID_R3_HII
O_H_NEBULABAYES_SF_CI68_LOW           NB_LOGU_RELIABLE_SF                   JY22_RESID_R3_SF
O_H_NEBULABAYES_HII_CI68_HIGH         NB_LOGPK_RELIABLE_HII                 JY22_RESID_NORM_HII
O_H_NEBULABAYES_SF_CI68_HIGH          NB_LOGPK_RELIABLE_SF                  JY22_RESID_NORM_SF
O_H_NEBULABAYES_HII_MODE              NB_VALID_HII                          JY22_FLAG_HII
O_H_NEBULABAYES_SF_MODE               NB_VALID_SF                           JY22_FLAG_SF
LOG_U_NEBULABAYES_HII                 O_H_JY22_HII                          JY22_FIT_OK_HII
LOG_U_NEBULABAYES_SF                  O_H_JY22_SF                           JY22_FIT_OK_SF
LOG_U_NEBULABAYES_HII_CI68_LOW        O_H_JY22_HII_16                       JY22_VALID_HII
LOG_U_NEBULABAYES_SF_CI68_LOW         O_H_JY22_SF_16                        JY22_VALID_SF
```
