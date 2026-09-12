---
title: "Remaining Star Formation and Rising NSF Occupancy"
subtitle: "An observational audit, literature synthesis and analytical scenario for MAUVE"
author: "Research synthesis prepared for Rongjun / MAUVE"
date: "11 September 2026"
lang: en
---

# 1. Executive answer

**The most economical working scenario is outside-in gas removal combined with fading of the surviving young-star component, while some inner gas remains detectable under a changing mixture of ionizing sources.** Low-column outer regions tend to leave the detected population and enter ND. In the inner, denser regions, retained gas can remain Balmer-detected after its emission ceases to satisfy the SF selection; those regions enter NSF. Older stars, radiation leaking from surviving H II regions, and in some systems shocks, outflows or an AGN can maintain this residual emission. The relative importance of those components can rise even while their absolute luminosities decline.

This scenario connects the observations through **gas retention, young-star luminosity, residual ionization and the actual selection rules**. It also admits a compact analytical model: a gas reservoir feeds a fading H II component; a second emission component mixes with it; and explicit detection, equivalent-width and line-width thresholds convert the luminosity distribution into SF/NSF/ND fractions. Closed-form expressions below predict both occupancy and conditional intensity.

The present work establishes a **testable explanatory model, not a fitted evolutionary history**. The model reproduces the main qualitative combinations, including increasing inner NSF area with decreasing NSF H-alpha intensity. Neither the available stage averages nor the analytical examples uniquely identify a quenching time, ram pressure, star-formation efficiency or ionizing source. Several observations actively constrain how simple the scenario can be:

1. SF loses both occupancy and, over much of the stellar-density range, intensity within the surviving selection. Geometric truncation alone is therefore an incomplete description unless selection and population differences explain the intensity change.
2. NSF is a residual observational class, not a direct measurement of DIG, old-star ionization or quenched gas. Its larger area fraction does not establish a larger non-SF luminosity budget.
3. All forbidden-line ratios do **not** rise everywhere. In the current products, NSF [O III]/H-beta decreases across several high-stellar-density bins from pre-peak to post-peak. A spatially uniform hardening sequence fails this check.
4. Close-to-peak galaxies do not consistently sit between the other stages. The samples are different galaxies, and the stage labels are not a calibrated common clock.
5. Detection and classification can generate part of the pattern. They must be included in a forward model before interpreting a transition density or a broad-line NSF excess as a physical threshold.

The strongest next step is a joint, galaxy-level comparison of **occupancy, surviving SF intensity, absolute NSF emission and excitation**, with molecular gas and stellar-population constraints used to distinguish fuel loss, inefficient star formation and changing ionization.

## 1.1 What was examined and recalculated

The starting document was the explicitly requested eight-page report, *20260909 Remained SF and Increasing NSF along Infall Stages.pdf*. Its accompanying Markdown and related reports from 27 August, 31 August, 1 September and 6 September were read; relevant sections of earlier analytical reports were used for context. The current 9 September notebooks in `further` were then inspected to recover the definitions and estimators behind the plots. Their latest contents include diagnostics beyond the named PDF.

Selected data-loading and summary cells were executed against the live FITS products without modifying the notebooks or products. New scalar summaries and a 10,000-draw whole-system bootstrap were generated for this report. Both target and pre-peak samples were resampled, with the same system draws retained across bins and metrics. The full source notebooks, fitting pipeline and all original plotting cells were not rerun. Provenance, the scope of execution and the independent numerical checks appear in Section 10.

The literature search was a targeted deep review through **11 September 2026**, covering analytical stripping, gas regulation, Virgo molecular gas, DIG/LIER excitation, emission-line fading, shocks and contrary cases of SF enhancement. It was not a formal systematic review or a claim of exhaustive coverage. Section 11 gives 26 primary-source records or author manuscripts, their relevance and access depth.

# 2. What the plots actually measure

## 2.1 Sample, masks and support

The freshly loaded sample contains **26 analysis products: 8 pre-peak, 5 close-to-peak and 13 post-peak**. NGC4567_8 is one combined product and one resampling unit, so these products encompass 27 physical galaxies. The configured post-peak list has 16 candidates; NGC4548, NGC4579 and NGC4689 lack the required products in this run. The stage comparison must therefore use the loaded sample rather than the nominal catalogue count.

The underlying valid domain requires finite stellar surface density, a finite gas-bin identifier, `MASK = 0`, and a finite elliptical radius. The retained domain additionally requires **finite `SNR_POSTFIT > 25`**, with a strict inequality. This is a post-fit continuum-quality selection. It is separate from emission-line detection, which uses raw line S/N at least 3 and raw flux at least 20 in the pipeline units of $10^{-20}\,{\mathrm{erg}\,s^{-1}\,cm^{-2}}$.

The configuration excludes inclinations at or above 80 degrees, checks WCS agreement to 0.05 arcsec, and does not apply an additional inclination correction to the supplied surface-density maps. Elliptical radius comes from the adopted geometry; the combined NGC4567_8 product uses the smaller component-normalized radius. These conventions should remain fixed when fitting the current observables.

Profiles use 0.25-dex stellar-density bins and $0.25R_e$ radial bins. The displayed windows are

$$
7\leq \log_{10}\!\left(\frac{\Sigma_\star}{M_\odot\,{\mathrm{kpc}}^{-2}}\right)<9.5,
\qquad 0\leq R/R_e<2.5. \tag{1}
$$

A galaxy-bin normally requires at least 50 retained pixels. Category-conditional intensities and line diagnostics usually require at least 20 category pixels, with additional common-line support where applicable. **A displayed window is not itself a reliability cut.** The three-galaxy support threshold is diagnostic rather than an enforced mask in all outputs. For example, the highest displayed stellar-density bin has only one close-to-peak contributor, and no post-peak contributor with 20 SF pixels even though a finite occupancy-weighted SF intensity can still be calculated. Such points should not anchor physical claims.

## 2.2 SF, NSF and ND are operational classes

Within the valid domain, the SF selection requires all of the following:

- finite `LOGSFR_SURFACE_DENSITY_HII`, inherited from the pipeline's H II selection using both BPT diagrams;
- finite `PROXY_EWHA` greater than 6 Angstrom;
- finite intrinsic H-alpha dispersion below $45\,{\mathrm{km}\,s^{-1}}$.

Intrinsic dispersion is computed from the square root of observed variance minus the correction variance only when that difference is positive. Unresolved or invalid differences become non-finite and can fail the SF gate. The pipeline configuration uses measured weak non-Balmer fluxes in its BPT classification rather than a fully limit-aware classification.

ND comprises pixels for which **both Balmer lines are not jointly detected**. NSF comprises pixels with both Balmer lines detected that fail the SF selection. Thus NSF can contain low-EW gas, broad gas, non-H II excitation, uncertain BPT classifications, and failures caused by missing EW or dispersion information. It can also contain emission powered partly by young stars. It is not equivalent to zero star formation.

For the retained domain $U$, the intended partition is

$$
U=U_{\mathrm{SF}}\,\dot\cup\,U_{\mathrm{NSF}}\,\dot\cup\,U_{\mathrm{ND}},
\qquad f_{\mathrm{SF}}+f_{\mathrm{NSF}}+f_{\mathrm{ND}}=1. \tag{2}
$$

The partition was checked on the extracted live data. Its physical meaning is **selected pixel occupancy in the MUSE footprint**, not the fraction of an entire optical disc that is forming stars.

Two consequences follow. First, interpreting an NSF line-width excess as independent proof of stronger turbulence is partly circular because width helps define the class. Second, a first-failed-gate diagnostic, evaluated in the order H II classification, EW and dispersion, attributes overlapping failures to the first gate. Those bars are not independent causal contributions.

## 2.3 Separate occupancy, surviving intensity and total contribution

For system $g$ in bin $b$, let $N_{U,g}$ be the usable count and $N_{{\mathrm{SF}},g}$ the SF count. Suppressing the bin index,

$$
f_g=\frac{N_{{\mathrm{SF}},g}}{N_{U,g}},\qquad
T_g=\frac{\sum_{i\in{\mathrm{SF}},g}\Sigma_{{\mathrm{SFR}},i}}{N_{U,g}}.
\tag{3}
$$

The SF-intensity notebook defines the stage statistics

$$
O=\langle f_g\rangle_g,\qquad T=\langle T_g\rangle_g,
\qquad I=\frac{T}{O}. \tag{4}
$$

The decomposition $T=OI$ is exact. Here $O$ answers where SF remains; $I$ measures the intensity of the surviving selection under this particular weighting; and $T$ measures the selected SF contribution per usable area. Crucially,

$$
I=\frac{\sum_g f_g I_g}{\sum_g f_g}, \tag{5}
$$

where $I_g$ is a within-SF mean when SF pixels exist. **This $I$ is not the arithmetic mean of galaxy-level SF intensities.** Systems with more SF occupancy contribute more to it, even though $O$ and $T$ were constructed as equal-system means. Galaxies with zero SF contribute zero to $O$ and $T$; no physical intensity is assigned to their nonexistent SF selection.

The H-alpha category notebook uses a different conditional estimator. For category $c$ it calculates an equal-system mean fraction $F_c$ and an equal-system mean within-category luminosity surface density $I_c$, with the latter evaluated only where category support passes. It also now includes

$$
J_c=\left\langle\frac{\sum_{i\in c,g}\mathcal{L}_{\alpha,i}}{N_{U,g}}\right\rangle_g
=\langle f_{c,g}I_{c,g}\rangle_g. \tag{6}
$$

$J_c$ measures the category's absolute H-alpha contribution per usable area. In general,

$$
J_c\ne F_c I_c. \tag{7}
$$

Even on identical support, covariance between fraction and intensity prevents factorization; here the conditional-support requirements can also differ. Likewise, the mean of within-system luminosity shares is not generally the ratio of mean luminosities. A rising $F_{\mathrm{NSF}}$ and falling $I_{\mathrm{NSF}}$ alone do not determine the sign of $J_{\mathrm{NSF}}$.

## 2.4 H-alpha luminosity, SFR and the ND surrogate

The current product conversion is

$$
\log_{10}\mathcal{L}_{\alpha}
=\log_{10}\Sigma_{\mathrm{SFR}}-\log_{10}C_{\mathrm{SFR}},
\qquad C_{\mathrm{SFR}}\simeq4.98358\times10^{-42}, \tag{8}
$$

with the supplied Chabrier-based calibration and Calzetti attenuation treatment. $\mathcal{L}_{\alpha}$ is in ${\mathrm{erg}\,s^{-1}\,kpc^{-2}}$ when the SFR surface density is in $M_\odot\,{\mathrm{yr}^{-1}\,kpc^{-2}}$. For NSF this inversion is a way to recover **luminosity**, not evidence that the corresponding SFR calibration is valid.

For ND, the all-spaxel map substitutes a quantity involving the larger of noise and measured H-alpha flux when either Balmer line fails. This is not a calibrated probabilistic upper limit and lacks a reliable Balmer-decrement correction. Consequently, ND luminosity shares should be treated as sensitivity surrogates. ND does not mean no gas, no H-alpha or zero SFR. Claims that one of these plotted estimators is a strict upper or lower bound on the true SFR require assumptions about censoring, attenuation, calibration and contamination that have not been established here.

# 3. The empirical pattern after checking the analysis

## 3.1 The stage sequence has two spatial branches

Figure 1 reproduces the central occupancy estimates from the live products. Post-peak systems show strong SF losses at low stellar density and large radius, where ND becomes prominent. At high density and small radius, much of the retained area stays Balmer-detected but enters NSF. This is the basic spatial distinction the physical model must explain.

![Figure 1. Live equal-system SF, NSF and ND fractions under the strict post-fit S/N cut. Bands are pointwise 16th-84th percentiles from 10,000 whole-system resamples. Open markers denote fewer than three contributing systems under the plotted metric's support rule. The three columns share a partition; their uncertainties are correlated.](assets/20260911_sf_nsf_physical_model/figure_01_live_occupancy.png)

At $\log\Sigma_\star=8.625$, the pre-, close- and post-peak SF fractions are approximately **0.499, 0.595 and 0.311**, while NSF fractions are **0.469, 0.378 and 0.563**. The close-to-peak value is not an intermediate point. At $R/R_e=0.125$, SF fractions are approximately **0.515, 0.530 and 0.247**, and NSF fractions are **0.472, 0.451 and 0.694**. At $R/R_e=1.625$, the corresponding SF fractions fall from about **0.429 to 0.241 to 0.027**, with ND rising from about **0.317 to 0.586 to 0.789**.

These examples support an inner-detected/outer-undetected contrast, but they do not identify a universal transition at a particular stellar density. Stellar density and radius are correlated; sensitivity, coverage, morphology and galaxy support change across both profiles.

## 3.2 Remaining SF often fades as well as becoming rarer

The table gives selected **post-peak minus pre-peak logarithmic contrasts**, in dex. Brackets are newly calculated pointwise 16th-84th percentile intervals with both stage samples resampled. They describe the finite sample; they are not simultaneous confidence bands or a many-bin significance test.

| Bin centre | $\Delta\log O$ | $\Delta\log I$ | $\Delta\log T$ |
|---|---:|---:|---:|
| $\log\Sigma_\star=7.625$ | $-1.011\ [-1.403,-0.801]$ | $-0.429\ [-0.522,-0.285]$ | $-1.440\ [-1.783,-1.209]$ |
| $\log\Sigma_\star=8.625$ | $-0.204\ [-0.381,-0.029]$ | $-0.408\ [-0.530,-0.261]$ | $-0.613\ [-0.799,-0.406]$ |
| $R/R_e=1.625$ | $-1.204\ [-1.432,-1.023]$ | $+0.105\ [-0.204,+0.348]$ | $-1.099\ [-1.510,-0.774]$ |

The low-density example loses roughly an order of magnitude in occupancy and another factor of about 2.7 in selected intensity. At intermediate/high density both terms contribute. In the outer radial example, however, the small surviving SF population does not show a secure intensity decrease. Survivor selection can preserve only the bright tail after extensive removal. Thus neither uniform fading nor uniform intensity preservation is an adequate global description.

![Figure 2. Fresh post/pre and close/pre logarithmic contrasts for SF occupancy, selected intensity and their product, followed by NSF fraction, within-category H-alpha intensity and absolute contribution per usable area. Bands resample both stages. NSF fraction and intensity use different support rules, so their plotted shifts need not add to the J shift. Open markers flag fewer than three contributors under the corresponding support rule. Extreme weak-support edge values can fall outside the displayed vertical range.](assets/20260911_sf_nsf_physical_model/figure_02_live_decomposition.png)

The original SF-intensity contrast bootstrap held the pre-peak reference curve fixed while resampling the target. The new bootstrap includes reference uncertainty. This changes the uncertainty interpretation without changing the central estimator. Zero or absent contributions can make logarithmic draws undefined; finite-draw fractions are retained in the output tables, and those tails should not be interpreted through an ordinary symmetric log-error approximation.

## 3.3 Increasing NSF occupancy does not mean increasing NSF power

At $\log\Sigma_\star=8.625$, the post/pre NSF fraction contrast is $+0.080$ dex, its within-category H-alpha intensity contrast is $-0.508$ dex, and $J_{\mathrm{NSF}}$ changes by $-0.351$ dex. The respective intervals are $[-0.038,+0.220]$, $[-0.715,-0.143]$ and $[-0.566,-0.022]$. At $R/R_e=0.125$, the fraction rises by $+0.167$ dex $[+0.056,+0.297]$, whereas intensity decreases by $-0.576$ dex $[-0.941,-0.042]$. The central $J$ change is $-0.206$ dex with a wide interval $[-0.613,+0.351]$.

The physical requirement is therefore more precise than "NSF grows": **an increasing fraction of selected inner area has NSF classifications, and the mean emission within that class is fainter.** The absolute central NSF contribution is less certain. Newly entering faint regions can dilute a smaller bright NSF population, and cross-stage changes in the contributing galaxies can affect the same means. A brightening shock or AGN component in a subset is not excluded merely by a declining class-average intensity.

## 3.4 Excitation needs a qualified, line-specific interpretation

The common-support line profiles require both Balmer lines, [O III], [N II] and both [S II] lines to pass the raw detections. Within each category they use the same pixels and require sufficient common-line support. Stage means remain equal-system means, but the systems and pixels are not matched parcels of gas followed in time. Ratios of these mean line profiles also differ from the notebook's separate within-system median NSF-minus-SF diagnostics.

![Figure 3. Stage contrasts of NSF direct [O III]/H-beta, [N II]/H-alpha and [S II]/H-alpha, from left to right, on the shared detected-line support. The bottom row uses radius instead of stellar density. Bands resample whole systems in both stages. The sign of the [O III] change depends on the binning variable and spatial regime; a universal increase is inconsistent with these profiles.](assets/20260911_sf_nsf_physical_model/figure_03_live_line_ratios.png)

For example, at $\log\Sigma_\star=8.625$, post/pre contrasts are approximately $+0.181$ dex in [N II]/H-alpha, $+0.043$ dex in [S II]/H-alpha, and **$-0.077$ dex in [O III]/H-beta**. The [O III] contrasts become about $-0.151$ and $-0.200$ dex at centres 8.875 and 9.125. At $R/R_e=0.625$, the same three contrasts are approximately $+0.428$, $+0.184$ and $+0.084$ dex. These are different weighted mixtures of the sample, not contradictory measurements of the same material.

A larger NSF-minus-SF contrast in post-peak galaxies can occur because the SF ratio falls more rapidly, even if the absolute NSF ratio also falls. It should not be restated as an absolute stage increase in NSF excitation. Metallicity, ionization parameter, spectral hardness, mixing and the contributing population can all affect these ratios; "hardness" is not a single monotonic interpretation of every BPT axis.

### A specific H-beta denominator issue

The current line-profile cell uses the stage change in [O III]/H-alpha as a proxy for the stage change in [O III]/H-beta, justified by a supposedly fixed corrected Balmer decrement of 2.86. The pipeline clamps observed decrements below 2.86 to 2.86 **when estimating attenuation**. For those pixels it applies zero attenuation; it does not force their corrected H-alpha/H-beta to 2.86. Across retained NSF pixels, observed decrements below 2.86 occur in about **7.2%, 7.2% and 27.5%** of the pre-, close- and post-peak samples.

This justification is therefore incorrect. However, directly using H-beta on the exact common-line, luminosity-weighted support changes the plotted NSF stage contrasts by **less than 0.0096 dex** across the finite bins checked. Bright-pixel weighting and the additional common-line cuts make the numerical effect small here. Figure 3 uses direct H-beta. The separate paired-median diagnostic already uses actual H-beta. No source notebook was edited in this task.

# 4. What the literature contributes

## 4.1 Gas removal and reduced replenishment are a useful starting point

Analytical stripping work distinguishes long pressure pulses, where the restoring-force criterion is useful, from short pulses governed more strongly by impulse. This cautions against mapping an infall label directly to an instantaneous stripping radius. It supplies the physical outer-disc boundary of the model below. [Koppen et al. 2018](https://arxiv.org/abs/1806.05887)

Virgo observations provide a close empirical comparison. VERTICO VII finds that reduced molecular-gas content and reduced star-formation efficiency both contribute to low SFR in environmentally affected galaxies. VERTICO V finds changes inside the cold-gas disc as well as preferential loss of low-surface-density molecular material. These results motivate testing both a gas-quantity term and an efficiency term in MAUVE, without assuming either term is already measured by the present optical diagrams. [Brown et al. 2023](https://arxiv.org/abs/2308.10943), [Watts et al. 2023](https://arxiv.org/abs/2303.07549)

A gas-regulation framework supplies the continuity equation linking inflow, gas mass, star formation and outflow. Applied locally, it is a closure model rather than a direct measurement of cloud-scale transport. NGC4330 demonstrates that spatially resolved quenching histories can constrain a radially varying clock in an individual stripped galaxy; those ages cannot simply be assigned to every MAUVE stage. [Lilly et al. 2013](https://arxiv.org/abs/1303.5059), [Fossati et al. 2018](https://arxiv.org/abs/1801.09685)

## 4.2 More diffuse-looking area can coexist with a small luminosity contribution

GASP XXXII is particularly relevant because it separates DIG-dominated area from integrated DIG emission. It finds a greater DIG-dominated area at fixed SFR surface density in stripped systems without a clear corresponding difference in integrated DIG flux fraction. Its classifications are not identical to MAUVE NSF, but it is a direct precedent for treating area and luminosity as different observables. [Tomicic et al. 2021a](https://arxiv.org/abs/2011.08869)

PHANGS-MUSE work shows how leakage from H II regions and hot evolved stars can contribute differently to diffuse gas. A falling bright H II component can reveal a relatively harder background, particularly in central low-ionization regions. MaNGA studies independently link extended low-EW LIER emission to old stellar populations and show that diffuse gas alters strong-line ratios and inferred abundances. None of these results makes every MAUVE NSF pixel an old-star-powered region. [Belfiore et al. 2022](https://arxiv.org/abs/2111.14876), [Belfiore et al. 2016](https://doi.org/10.1093/mnras/stw1234), [Zhang et al. 2017](https://arxiv.org/abs/1612.02000)

EW offers a useful photon-budget clue, but published WHAN and CALIFA thresholds describe particular data and physical mixtures. They are not an automatic calibration of `PROXY_EWHA`, and a numerical EW boundary is not a universal DIG fraction. Resolution, stellar continuum and tail geometry matter. [Cid Fernandes et al. 2011](https://doi.org/10.1111/j.1365-2966.2011.18244.x), [Lacerda et al. 2018](https://doi.org/10.1093/mnras/stx3022)

## 4.3 Shocks and outflows are plausible additions, not automatic consequences

Stripped tails can require heating beyond ordinary H II-region photoionization, especially when [O I] is enhanced. ESO137-001 and GASP studies illustrate mixtures of photoionization, shocks and other heating rather than one universal component. Their tail diagnostics should not be transferred uncritically to all central MAUVE NSF regions. [Fossati et al. 2016](https://doi.org/10.1093/mnras/stv2400), [Tomicic et al. 2021b](https://arxiv.org/abs/2108.12433)

Projection, rotation and beam smearing can also broaden profiles. Joint excitation and kinematic decomposition is more informative than a width threshold alone. JO201 and spatial decomposition work in NGC1068 show why multiple components and spatial context matter. [Bellhouse et al. 2019](https://arxiv.org/abs/1902.04486), [D'Agostino et al. 2019](https://arxiv.org/abs/1906.07907)

The MAUVE study of NGC4064 is a direct warning against a purely passive interpretation: a SF-driven outflow can coexist with a galaxy undergoing quenching. It supports allowing a physically motivated outflow component in individual systems. Its measured mass-loading estimate is not a justified universal parameter for the other 25 products. [Attwater et al. 2025, accepted manuscript](https://arxiv.org/abs/2512.10574)

## 4.4 A fading stellar population does not make every forbidden line linger

After the ionizing stellar spectrum softens, high-ionization emission can fade rapidly. Photoionization calculations specifically predict decreasing [O III]/H-alpha after star formation shuts down. Therefore "forbidden lines fade more slowly than Balmer lines" cannot serve as a general explanation of the MAUVE sequence. Continued excitation or a changing source mixture is needed where ratios increase; a softer field can help explain where [O III] decreases. [Citro et al. 2017](https://arxiv.org/abs/1704.05462)

H-alpha responds to recent massive-star formation, with a response kernel that depends on the stellar population and radiative transfer. It is not an instantaneous gas-consumption meter. Dust, photon escape and diffuse emission affect its relation to SFR. Small spatial apertures also sample different phases of the gas/star-formation lifecycle. [Tacchella et al. 2022](https://doi.org/10.1093/mnras/stac818), [Kennicutt & Evans 2012](https://arxiv.org/abs/1204.3552), [Kruijssen & Longmore 2014](https://arxiv.org/abs/1401.4459)

## 4.5 Contrary evidence constrains a universal stage story

RPS need not suppress SF immediately. Pressure and radial gas transport can transiently increase dense central gas; observations of other stripped samples find enhanced resolved SF. Conversely, a cosmological jellyfish population need not show population-wide enhancement even when individual objects undergo bursts. Selection, orbit, gas content and the reference population therefore matter. [Fujita & Nagashima 1999](https://arxiv.org/abs/astro-ph/9812378), [Zhu et al. 2024](https://arxiv.org/abs/2309.07037), [Vulcani et al. 2020](https://arxiv.org/abs/2007.04996), [Goller et al. 2023](https://arxiv.org/abs/2304.09199)

Recent FUV/H-alpha work gives additional constraints on knots and excitation, but does not support assuming the same ordering of UV and H-alpha extent in every stripped disc. Recent TNG50 work also describes slow transformations in some environments. Such results broaden the allowable histories rather than calibrating a unique time separation between the MAUVE stages. [George et al. 2025](https://arxiv.org/abs/2505.15066), [Lora et al. 2026](https://arxiv.org/abs/2608.11336)

# 5. A minimal physical scenario

## 5.1 Before strong stripping

The pre-peak sample provides an empirical starting distribution of gas, SF and excitation, but it is already a cluster sample and is not a guaranteed unprocessed field control. A galaxy has a distribution of gas columns and cloud evolutionary states at each radius. Young stars power bright H II regions; diffuse radiation and older stars contribute a background whose importance varies with stellar density, gas covering fraction and source geometry.

## 5.2 During the pressure encounter

External pressure preferentially removes weakly bound, low-column material and alters the supply feeding denser star-forming gas. Some clouds remain, some are compressed, and some gas can move inward. Hence SF occupancy can fall outside while remaining temporarily high or increasing in the centre. This allows close-to-peak departures from a monotonic average sequence without requiring a different global mechanism for every point.

## 5.3 After substantial gas loss

With reduced replenishment, molecular reservoirs can decline and the young-star ionizing luminosity falls. If efficiency also falls, fading can be faster than gas depletion alone predicts. Where the remaining gas is optically detectable, an older stellar continuum and a relatively persistent residual ionizing component push EW and line ratios across the observational selection boundaries. A change in component weights can also increase the measured line width without requiring an increase in the intrinsic width of either component.

This creates an **inner route from SF to NSF**. Farther out, lower residual luminosity and stricter effective Balmer sensitivity create a **route toward ND**. These are statistical routes through observed classes, not tracked trajectories of individual pixels. Some outer material may pass briefly through NSF, while some inner material may already have been NSF before the encounter.

## 5.4 What "persistent residual emission" means

Persistence refers to continued energy input on the relevant interval. Recombination emission is not a long-lived store of the pre-stripping ionization state. Symbolically,

$$
t_{\mathrm{rec}}=\frac{1}{\alpha_B(T)n_e}. \tag{9}
$$

The model therefore requires surviving gas plus continuing ionization/heating. Old stars can supply an ionizing field tied broadly to the stellar population, but the emitted luminosity still depends on gas absorption and geometry. Leaked young-star photons generally fade with their sources unless escape or absorption fractions change. Shock or AGN emission needs its own energy source. A constant additive background is an approximation over a limited interval, not a claim that residual gas radiates forever.

# 6. Analytical model: from gas to the observed classes

## 6.1 Gas retention and an outside-in boundary

The incident ram pressure is

$$
P_{\mathrm{ram}}=\rho_{\mathrm{ICM}}v_{\mathrm{rel}}^{\,2}. \tag{10}
$$

A useful local restoring-pressure approximation is

$$
P_{\mathrm{rest}}\simeq\Sigma_g\max_z\left|\frac{\partial\Phi}{\partial z}\right|
\simeq2\pi G\Sigma_\star\Sigma_g. \tag{11}
$$

The second form assumes an approximately stellar-sheet restoring force. Gas self-gravity, the halo, disc thickness and wind orientation can matter. A long, sustained pressure pulse can be compared with this restoring pressure; for a short pulse, a schematic impulse condition instead involves

$$
\int P_{\mathrm{ram}}(t)\,dt\ \gtrsim\ \Sigma_g v_{\mathrm{esc},z}. \tag{12}
$$

This distinction is central to the analytical treatment of [Koppen et al. (2018)](https://arxiv.org/abs/1806.05887). Neither expression makes projected cluster position a unique pressure measurement.

For exponential stellar and gas discs, using the long-pulse approximation gives

$$
R_{\mathrm{strip}}=\left(R_\star^{-1}+R_g^{-1}\right)^{-1}
\ln\!\left(\frac{2\pi G\Sigma_{\star,0}\Sigma_{g,0}}{P_{\mathrm{peak}}}\right). \tag{13}
$$

Negative values correspond to no radius satisfying this simplified retention condition; large values are limited by the actual disc. A real multiphase medium produces a soft boundary. If $\ln\Sigma_g$ has mean $\mu_g$ and dispersion $s_g$ in a local patch, an illustrative retained-area fraction is

$$
q_{\mathrm{ret}}=1-\Phi\!\left(\frac{\ln\Sigma_{g,\mathrm{crit}}-\mu_g}{s_g}\right),
\qquad \Sigma_{g,\mathrm{crit}}=\frac{P_{\mathrm{peak}}}{2\pi G\Sigma_\star}. \tag{14}
$$

$\Phi$ is the standard normal cumulative distribution. This is a retention model, not yet an SF fraction: retained gas need not be bound, molecular or actively forming stars. It explains why a sharp physical forcing can yield a gradual occupancy change across a broad gas-column distribution.

## 6.2 A two-reservoir model for fading SF

Let $a$ be a diffuse gas reservoir that can feed molecular gas $m$, both expressed per the **same fixed physical area**. Let $t_f$ be the conversion/replenishment time, $t_d$ the molecular depletion time, $k_a$ and $k_m$ effective removal rates, $\mathcal R$ the instantaneous returned fraction and $\eta$ the outflow mass loading. After external supply is cut, a minimal closure is

$$
\psi=\frac{m}{t_d},\qquad
\dot a=-\left(k_a+t_f^{-1}\right)a\equiv-\kappa a, \tag{15}
$$

$$
\dot m=\frac{a}{t_f}-\lambda m,
\qquad \lambda=k_m+\frac{1-\mathcal R+\eta}{t_d}. \tag{16}
$$

Here $\psi$ is the actual young-star SFR surface density in the model. Recycling is returned instantaneously to the tracked molecular reservoir as a simplifying closure; explicit dissociation and reverse phase exchange are omitted. Those assumptions are not statements about detailed ISM chemistry. They permit a transparent analytical baseline, motivated by gas regulation. [Lilly et al. 2013](https://arxiv.org/abs/1303.5059)

For constant coefficients and $\lambda\ne\kappa$,

$$
a(t)=a_0e^{-\kappa t}, \tag{17}
$$

$$
m(t)=m_0e^{-\lambda t}
+\frac{a_0}{t_f}\frac{e^{-\kappa t}-e^{-\lambda t}}{\lambda-\kappa}. \tag{18}
$$

When $\lambda=\kappa$, the second term is $(a_0/t_f)t e^{-\lambda t}$. The total gas budget is

$$
\frac{d(a+m)}{dt}=-k_a a-k_m m-(1-\mathcal R+\eta)\psi. \tag{19}
$$

Phase transfer cancels, as conservation requires. The model can show an initial increase in $m$ if feeding temporarily exceeds its losses, even while total gas declines. Later SF fades as the feed reservoir disappears. Thus a temporary central enhancement and later suppression can be two parts of the same history.

For a sufficiently supply-poor interval, $\psi\propto e^{-t/\tau_q}$. A measured decline then constrains $t/\tau_q$, not $t$ and $\tau_q$ separately. More generally, on matched support,

$$
\Delta\log\psi=\Delta\log\Sigma_{\mathrm{mol}}-\Delta\log t_d. \tag{20}
$$

CO measurements are needed to separate these two terms. Optical SFR alone does not distinguish losing molecular material from forming stars less efficiently. If an area-retention factor is introduced separately, the reservoir densities must be defined per retained area before multiplying by that factor. Applying $q_{\mathrm{ret}}$ to a reservoir already measured per total area would double-count removal.

## 6.3 The young-star H-alpha response

Let $A(t)$ denote the H-alpha luminosity surface density powered by young stars. With a normalized response kernel $K_\alpha$, it is

$$
A(t)=C_\alpha\int_0^\infty\psi(t-u)K_\alpha(u)\,du. \tag{21}
$$

For an illustrative exponential kernel with timescale $\tau_{\mathrm{ion}}$, constant pre-encounter SFR and subsequent exponential fading with timescale $\tau_q$,

$$
\frac{A(t)}{A(0)}=
\frac{\tau_q e^{-t/\tau_q}-\tau_{\mathrm{ion}}e^{-t/\tau_{\mathrm{ion}}}}
{\tau_q-\tau_{\mathrm{ion}}}. \tag{22}
$$

This reduces to the long-timescale SF decline after the short response transient. At equal timescales the limit is $(1+t/\tau_q)e^{-t/\tau_q}$. A few-Myr illustrative kernel does not assign a stellar age to an infall stage; population synthesis, dust and photon transport are required for a calibrated kernel. [Kennicutt & Evans 2012](https://arxiv.org/abs/1204.3552), [Tacchella et al. 2022](https://doi.org/10.1093/mnras/stac818)

## 6.4 Two-component line mixing: dimmer gas can have larger ratios

Let $D(t)$ be H-alpha emission from the residual component and $H=A+D$ the total. For a line $j$ whose component ratios are defined relative to H-alpha,

$$
L_j=r_{j,A}A+r_{j,D}D,\qquad
w=\frac{D}{A+D}, \tag{23}
$$

$$
\frac{L_j}{H}=r_{j,A}+(r_{j,D}-r_{j,A})w. \tag{24}
$$

Fluxes must be added before taking logarithms. If $A\propto e^{-t/\tau_A}$ and $D\propto e^{-t/\tau_D}$,

$$
\dot w=w(1-w)\left(\tau_A^{-1}-\tau_D^{-1}\right). \tag{25}
$$

When $\tau_D>\tau_A$, the residual fraction increases. If also $r_{j,D}>r_{j,A}$, the ratio increases although both component amplitudes, and hence $L_j$, decrease. No increase in absolute residual luminosity is required. If $r_{j,D}<r_{j,A}$, the ratio decreases instead. This latter possibility is essential for the high-density [O III] behaviour.

For [O III]/H-beta, the correct weights use H-beta. Writing component Balmer decrements as $b_A$ and $b_D$,

$$
B_A=A/b_A,\quad B_D=D/b_D,\qquad
\frac{L_{\mathrm{[O\,III]}}}{B_A+B_D}
=r_{3,A}+(r_{3,D}-r_{3,A})\frac{B_D}{B_A+B_D}. \tag{26}
$$

Using $w$ from H-alpha in this equation requires equal effective decrements. A real fit should allow line templates to vary with metallicity, ionization parameter and source spectrum. The algebra supplies the mixing relation; it does not supply the templates or prove which source powers $D$.

As a numerical example, choose $A_0=1$, $D_0=0.2$, $\tau_A=150$ Myr and $\tau_D=600$ Myr, in arbitrary common luminosity units. At 300 Myr, total H-alpha falls from 1.2 to about 0.257, a decline of 0.67 dex. With [N II]/H-alpha component ratios 0.3 and 1.2, [N II] falls from 0.54 to about 0.186, while its ratio to H-alpha rises from 0.45 to about 0.725. These are illustrative parameters, not MAUVE estimates.

## 6.5 EW and line width can change through mixing

For a local stellar continuum surface density $C_\lambda$, measured consistently with the line,

$$
W_\alpha=\frac{A+D}{C_\lambda}. \tag{27}
$$

If the continuum changes slowly while young-star emission fades, EW falls. The actual `PROXY_EWHA` must be calibrated against this expression rather than assumed identical to a standard spectroscopic EW.

For two intrinsic Gaussian line components with widths $\sigma_A$, $\sigma_D$ and centroid separation $\Delta v$, the second central moment is

$$
\sigma_{\mathrm{mix}}^2=(1-w)\sigma_A^2+w\sigma_D^2
+w(1-w)(\Delta v)^2. \tag{28}
$$

Changing component weights can broaden the total line even when neither component becomes more turbulent. The centroid term can peak at intermediate mixing; general monotonic broadening is not guaranteed. Instrumental broadening and spatial beam smearing must be applied before comparing with observed fitted widths. The SF width cut makes an independent test especially necessary.

![Figure 4. A two-component illustrative calculation: absolute lines fade while selected ratios increase, EW falls, and the mixed width can cross the SF cut. The second [O III] case shows that a weaker residual [O III]/H-beta component instead produces a declining ratio. The two [O III] cases are alternative templates, not components added together. All times and amplitudes are illustrative.](assets/20260911_sf_nsf_physical_model/figure_05_mixing_fading_width.png)

## 6.6 Closed-form occupancy and intensity

The gas and emission model becomes directly comparable to the classification plots only after a distribution and an observing rule are specified. For a simple local example, assume

$$
\ln A\sim\mathcal N(\mu,s^2), \tag{29}
$$

with fixed $D$, continuum and detection threshold in one bin. This represents unresolved spatial variation, not independent pixel trials. The fully realistic model needs their joint distribution, and can include an atom at $A=0$ for removed or inactive regions.

In a noiseless approximation with a fixed observed Balmer decrement $b$, joint Balmer detection implies an effective H-alpha threshold

$$
H>L_{\mathrm{det}}\equiv\max(L_{\alpha,\mathrm{lim}},bL_{\beta,\mathrm{lim}}),
\qquad A_d=\max(0,L_{\mathrm{det}}-D). \tag{30}
$$

For this illustrative calculation only, suppose the BPT and width conditions can be represented by a maximum permitted residual fraction $w_c$. The SF threshold becomes

$$
A_s=\max\!\left[A_d,\;6C_\lambda-D,\;
D\frac{1-w_c}{w_c},\;0\right]. \tag{31}
$$

The coefficient 6 carries the Angstrom unit of the EW boundary. This threshold is a simplified projection of the actual selection, not a universal BPT mixing fraction. For example, if $\Delta v=0$, $\sigma_A=25$ and $\sigma_D=70\,{\mathrm{km}\,s^{-1}}$, the $45\,{\mathrm{km}\,s^{-1}}$ cut alone gives

$$
w_{c,\sigma}=\frac{45^2-25^2}{70^2-25^2}\simeq0.3275. \tag{32}
$$

The illustrated occupancy calculation assumes the BPT cut is no more restrictive; a real forward model must evaluate both BPT boundaries using the summed line fluxes. Unequal extinction, line-dependent noise, missing measurements and unresolved widths all require the original observing operator.

Define $F_A(x)=\Phi[(\ln x-\mu)/s]$ for $x>0$, with $F_A(0)=0$. Then

$$
p_{\mathrm{ND}}=F_A(A_d),\qquad p_{\mathrm{SF}}=1-F_A(A_s), \tag{33}
$$

$$
p_{\mathrm{NSF}}=F_A(A_s)-F_A(A_d). \tag{34}
$$

The probabilities sum exactly to one. When $D>L_{\mathrm{det}}$, fading $A$ can transfer probability from SF into NSF while ND remains small. Where $D<L_{\mathrm{det}}$, stronger fading eventually sends the distribution below joint Balmer detection, increasing ND. NSF need not increase forever: it can peak and later disappear as the residual component also fades or is removed.

The unnormalized truncated first moment is

$$
M(l,u)=e^{\mu+s^2/2}
\left[\Phi\!\left(\frac{\ln u-\mu-s^2}{s}\right)
-\Phi\!\left(\frac{\ln l-\mu-s^2}{s}\right)\right]. \tag{35}
$$

The limits $l=0$ and $u=\infty$ are understood continuously. Conditional H-alpha intensities follow directly:

$$
I_{\mathrm{NSF}}=\frac{M(A_d,A_s)}{p_{\mathrm{NSF}}}+D,\qquad
I_{\mathrm{SF}}=\frac{M(A_s,\infty)}{p_{\mathrm{SF}}}+D. \tag{36}
$$

They are defined only when the relevant class probability is nonzero. This makes it possible for $p_{\mathrm{NSF}}$ to increase while $I_{\mathrm{NSF}}$ decreases. For selected SF emission, the measured all-usable-area contribution in the simplified model is

$$
T_{\mathrm{SF}}=C_{\mathrm{SFR}}\,p_{\mathrm{SF}}I_{\mathrm{SF}}. \tag{37}
$$

If the aim is the true young-star SFR, only the young-star component belongs in its calibration; the residual contribution in $I_{\mathrm{SF}}$ is potential contamination. Finally, these local predictions must be aggregated using the actual galaxy-level weights and support rules. Equation (37) does not justify replacing the observed NSF $J$ by the product of two separately averaged profiles.

![Figure 5. Closed-form population example for an inner and an outer environment. A lognormal young-star amplitude fades while the residual component declines more slowly. The inner residual remains detectable, so SF occupancy transfers mainly to NSF and its conditional intensity falls; the faint outer residual produces increasing ND. The time axis, continuum and amplitudes are chosen demonstrations, not fitted stage ages or a calibrated BPT model.](assets/20260911_sf_nsf_physical_model/figure_04_analytic_occupancy.png)

![Figure 6. Analytical gas-reservoir histories and an illustrative soft gas-column retention boundary. Reservoir transfer can briefly raise molecular gas before later fading, even as gas is removed. The right-hand retention calculation concerns gas columns, not the measured SF fraction. These panels test the behaviour of the equations and do not estimate ram pressure or removal rates for MAUVE.](assets/20260911_sf_nsf_physical_model/figure_06_gas_regulator.png)

# 7. Does this scenario account for the observables?

| Observable | Model ingredient that can produce it | What is still not established |
|---|---|---|
| Strong outer/low-density SF occupancy loss | Preferential removal plus a declining young-star amplitude | Gas loss versus fading below sensitivity in each individual region |
| Inner NSF occupancy increase | Detectable residual gas crosses BPT, EW or width gates as the young-star component fades | The source composition of NSF and the role of missing/uncertain diagnostics |
| Outer ND increase | Residual emission falls below joint Balmer detection | Whether undetected gas or low-level SF remains |
| Lower intensity in surviving SF over much of the density range | Less molecular gas, reduced feeding, increased depletion time, or changed selected population | Relative contributions of gas content, efficiency and selection |
| Larger NSF area but lower within-NSF H-alpha intensity | A growing population of faint detected regions and fading component amplitudes | The sign of the total NSF luminosity change in every region |
| Some rising low-ionization ratios despite fainter lines | Increasing relative contribution of a component with larger ratios | Whether that component is old-star ionization, leakage, shocks or AGN |
| Declining high-density [O III]/H-beta | Softer young-star field, a low-[O III] residual template, or population changes | A unique hardness or metallicity interpretation |
| Broader NSF profiles | Mixture weighting, centroid offsets, beam smearing, or actual disturbed motions | Independent evidence for added turbulent/shock energy |
| Non-monotonic close-to-peak profiles | Transient gas redistribution and heterogeneous encounter histories | A common timeline for all systems |

The model explains the **logical compatibility** of the trends with a small set of physical processes. It has not been optimized to the measured curves, and no goodness of fit, posterior parameter constraint or model-selection score is claimed. In particular, fixed line templates cannot capture every observed ratio trend; templates or mixtures must vary across the disc and between systems. That is a testable extension, not a license to tune each bin independently.

## 7.1 Models that are too restrictive

**Pure truncation with unchanged surviving clouds** predicts occupancy loss while the conditional SF population remains unchanged on genuinely matched support. It can explain an outer branch and survivor selection, but needs additional assumptions to account for widespread lower surviving intensity at fixed stellar density.

**Uniform fading of one emission component with fixed line ratios** reduces luminosities and EW but leaves all line ratios unchanged. It cannot explain the excitation changes. Allowing the young-star spectrum to soften changes this conclusion, but does not generically make [O III] stronger relative to Balmer emission.

**A uniformly brightening hard component** is unnecessary and is not favoured by the broadly declining within-NSF intensities. However, its contribution in a subset cannot be excluded through conditional means alone. The relevant test is its absolute luminosity and energy budget, with galaxy-level decomposition.

**One stage label equals one elapsed time** ignores orbital diversity, galaxy structure, prior processing and the non-monotonic close-to-peak profiles. A stage-conditioned distribution of histories is more defensible.

## 7.2 Energy checks that can make the model falsifiable

For an old-star interpretation, a stellar-population model supplies an ionizing-photon rate $Q_{\mathrm{old}}$. An approximate luminosity budget is

$$
L_{\alpha,\mathrm{old}}=\epsilon_\alpha Q_{\mathrm{old}}f_{\mathrm{abs}},
\qquad 0\leq f_{\mathrm{abs}}\leq1, \tag{38}
$$

where $\epsilon_\alpha$ is the case-dependent H-alpha energy emitted per absorbed ionizing photon. If the required $f_{\mathrm{abs}}$ exceeds unity under plausible population models, old stars alone fail. Agreement only establishes sufficiency, because other sources can contribute.

For incident-flow-powered emission, an order-of-magnitude mechanical ceiling is

$$
\dot E_{\mathrm{inc}}\sim\tfrac12\rho_{\mathrm{ICM}}v_{\mathrm{rel}}^{\,3}A_{\mathrm{intercept}},
\qquad L_{\alpha,\mathrm{shock}}\leq\epsilon_{\alpha,\mathrm{shock}}\dot E_{\mathrm{inc}}. \tag{39}
$$

The efficiency includes interception, coupling and the fraction radiated in H-alpha, and cannot exceed unity. Star-formation-driven outflows have a different mechanical budget. A measured optical line width is not the three-dimensional ICM wind speed and cannot be inserted for $v_{\mathrm{rel}}$ without justification.

# 8. A practical fitting and discrimination programme

## 8.1 First fit the observing process

The immediate target is a **joint model of local line emission and its classification**, not a regression through stage-average curves alone. For each system, predict gas retention, young-star amplitude, residual emission and continuum on a physical grid. Sum component line fluxes, apply attenuation, convolve with the spatial PSF and spectral response, add the measured noise model, and run the same line detections, BPT, EW and dispersion gates. Then apply the exact valid-domain and strict post-fit S/N masks, binning and support rules.

This order matters. Thresholding before convolution, assuming an ND surrogate is a flux limit, or fitting ratios without common line support would generate a different observable. The forward model should initially reproduce the pre-peak distributions and within-system correlations before interpreting shifts between stages.

Match or reweight the stage samples in the **joint** distribution of stellar density, radius, galaxy mass, morphology, inclination, resolution and sensitivity where support permits. The stellar-density and radial profiles are two projections of overlapping data, not two independent experiments. Reweighting may also reveal regions in which the samples do not provide a defensible comparison.

## 8.2 Use a small hierarchy of competing models

| Model | Added freedom | Discriminating prediction |
|---|---|---|
| M0: geometric removal | Radius/column-dependent retention only | Surviving intrinsic SF amplitude distribution is preserved on matched support |
| M1: removal plus fading | A small stage/radius-dependent shift in young-star amplitude | Occupancy and conditional intensity change together; fixed-template ratios do not |
| M2: fading plus residual ionization | Residual amplitude constrained by stellar populations and retained gas | Inner NSF growth, low EW and excitation changes can occur without increasing absolute residual power |
| M3: system-specific extra excitation | A constrained shock/outflow/AGN contribution | Excess line luminosity and kinematics beyond M2, localized to relevant systems |
| M4: observational/population null | Sensitivity and sample-composition differences without stage-dependent physics | Apparent stage trends shrink after matched-support forward modelling |

Start with a few shared parameters: a retention strength, a young-star amplitude shift with a low-order radial dependence, and a residual normalization linked to each galaxy's stellar population. Do not fit a different history, mixing fraction and line template in every bin. External gas and stellar-age data should constrain those freedoms before adding more terms.

There are fewer than a dozen independent systems in two stages. A covariance matrix for a large vector of binned summaries is therefore poorly determined; for five close-to-peak systems the sample covariance has rank at most four. Thousands of bootstrap draws do not create more independent galaxies. Use a small set of summaries with regularization checks or a hierarchical model of galaxy-level data, and check held-out systems.

## 8.3 Avoid counting algebraic identities as extra evidence

The three fractions have two independent degrees of freedom. $T=OI$ makes the SF decomposition redundant if all three quantities are treated as independent likelihood terms. Line fluxes and their ratios share information, as do $J$, category fractions and intensities. A likelihood or distance function must preserve these correlations. Millions of spatial pixels are not millions of independent environmental experiments.

Resample whole galaxies/products for stage inference, retaining NGC4567_8 as one unit. Within-galaxy errors can use spatial blocks or an explicit covariance model. Inspect leave-one-system-out changes and individual maps, including the influence of systems such as NGC4383, NGC4064 and NGC4694; do not exclude them merely because they complicate the average.

ND should enter through raw flux/error information or calibrated injection-and-recovery completeness. Its true flux can be positive below detection. Model-based censoring is preferable to assigning a physical upper limit to `max(NOISE, flux)`.

## 8.4 Measurements with the highest discriminating value

1. **CO and H I on matched physical support.** Measure whether lower SF intensity follows lower molecular surface density, a larger depletion time, or both. Check whether low-column gas is preferentially missing and whether the inner reservoir is also depleted. This directly tests Equations (14)-(20).
2. **Stellar populations and UV with explicit response kernels.** Compare recent SF histories with the young-star fading required by the optical lines. Use full histories and dust treatment rather than a universal H-alpha/FUV age conversion. The stellar continuum also constrains the old-star photon budget.
3. **Absolute NSF line emission and component decomposition.** Fit H-alpha, H-beta, [N II], [S II], [O III] and, where supported, [O I] jointly. Require the same component weights to explain multiple lines; freely fitting every ratio loses the test.
4. **Independent kinematics.** Examine line widths under a classification that does not itself use width, then account for rotation and PSF smearing. A BPT-only diagnostic split can be useful as a sensitivity experiment while preserving the original SF/NSF definition for the main result.
5. **Sensitivity and support experiments.** Repeat the observational prediction across measured noise levels and with calibrated completeness. A physical transition should not track an arbitrary detection boundary after sensitivity is accounted for.

For fixed templates, different ratios estimate a common mixing fraction through

$$
w_j=\frac{R_j-r_{j,A}}{r_{j,D}-r_{j,A}}, \tag{40}
$$

using the correct Balmer weighting. Values outside $[0,1]$, or inconsistent values across lines beyond errors and template uncertainty, reject that two-component template pair. This is a more restrictive test than merely showing that a mixing line passes through one BPT projection.

## 8.5 Outcomes that would change the interpretation

The residual-old-star explanation weakens if the stellar photon budget is insufficient, if the required gas absorption fraction is unphysical, or if spatially localized excitation and independently measured kinematics require substantial mechanical power. The fading interpretation weakens if matched-support stellar histories show no corresponding recent decline and the CO/SF relation is unchanged. A strong reduction of the stage trends after sensitivity and composition matching would shift emphasis toward observational selection. Conversely, successful prediction of held-out galaxies' occupancy, absolute lines and excitation with a common parameterization would provide substantially stronger support than fitting the three stage means.

# 9. How this revises the narrative of the recent reports

The recent reports contain a productive distinction between where SF remains and what happens inside the surviving selection. That distinction should remain the structure of the scientific argument. The present inspection adds four refinements.

First, **the all-usable-area selected SF term is a decomposition statistic**. It combines occupancy and intensity and is useful for bookkeeping, but it is not automatically a bound on total true SFR. Second, **NSF is a detected residual class** whose physical composition must be inferred, and ND is a sensitivity class. Neither should be relabelled as an evolutionary endpoint without additional evidence.

Third, **excitation contrasts must preserve their denominator and reference population**. Direct H-beta confirms that the specific denominator shortcut is numerically small for the common-support mean profiles; it does not restore a universal [O III] increase. Absolute NSF stage changes and NSF-minus-SF within-stage contrasts should be discussed separately.

Fourth, **a falling conditional NSF intensity does not by itself rule out a growing hard component in a subset**. The newly added $J$ and luminosity-share diagnostics make the question more concrete, but their averaging and ND treatment must remain explicit. The strongest phrasing supported now is that post-peak inner regions occupy the NSF class more frequently and are generally fainter within that class, with line-dependent excitation changes.

For a paper, the central physical claim can be framed as: the observations are consistent with preferential loss of outer low-column gas and fading of surviving young-star emission, while retained inner gas remains optically detectable under a changing ionization mixture. The data motivate this scenario; they do not yet select a unique mechanism or assign a universal evolutionary timescale.

# 10. Reproducibility, verification and limitations

## 10.1 Local inputs

All paths below identify current local sources inspected for this work. Source notebooks and FITS products were left unchanged.

| ID | Source and inspection scope |
|---|---|
| L1 | `/Users/Igniz/Desktop/ICRAR/MAUVE/20260909 Remained SF and Increasing NSF along Infall Stages.pdf` and its companion `.md`; named report, text and figures |
| L2 | `/Users/Igniz/Desktop/ICRAR/MAUVE/20260906 Where Star Formation Remains along Infall Stages.md`; recent occupancy interpretation |
| L3 | `/Users/Igniz/Desktop/ICRAR/MAUVE/20260901 Close-to-peak and Post-peak rSFMS with SF, NSF and ND.md`; recent resolved-SF interpretation |
| L4 | `/Users/Igniz/Desktop/ICRAR/MAUVE/20260831 Pre-Peak rSFMS with SF, ND and NSF.md`; pre-peak interpretation |
| L5 | `/Users/Igniz/Desktop/ICRAR/MAUVE/20260827 Resolved Star Formation Population through RPS.md`; population framework |
| L6 | `/Users/Igniz/Desktop/ICRAR/MAUVE/20260830_RPS_Local_Star_Formation_Theory_and_Toy_Models.md` and `/Users/Igniz/Desktop/ICRAR/MAUVE/20260813 Analytical Models of Galaxy Evolution in Cluster Environments and Ram Pressure Stripping.md`; selected contextual and analytical sections |
| N1 | `/Users/Igniz/Desktop/ICRAR/further/20260909_check_combined_SF_fraction_categories_by_stage_SNR_postfit.ipynb`; masks, fractions, strict cut, coverage and individual diagnostics |
| N2 | `/Users/Igniz/Desktop/ICRAR/further/20260909_check_Sigma_SFR_intensity_by_stage.ipynb`; definitions, loading, binning, decomposition and uncertainty implementation |
| N3 | `/Users/Igniz/Desktop/ICRAR/further/20260909_check_corrected_Halpha_surface_density_SF_NSF_ND_by_stage.ipynb`; category luminosities, line support, J, shares and paired diagnostics |
| P1 | `/Users/Igniz/Desktop/ICRAR/further/SFR+Z.py`; attenuation, line-detection and luminosity conventions |
| P2 | `/Users/Igniz/Desktop/ICRAR/further/mauve_master_wiki_newclass.fits` and `/Users/Igniz/Desktop/ICRAR/further/MAUVE_effective_radii.csv`; sample and geometry inputs |
| L7 | `/Users/Igniz/Desktop/ICRAR/further/20260905_stage_SF_fraction_analysis.md`; selected sample and definition checks |

The asset manifest records exact input hashes and the scalar extracts used by this report. The loaded sample identities are also retained in `sample_qc.csv`. Cell numbers below are one-based and include Markdown cells.

## 10.2 What was executed

N2 data cells 3, 5, 7, 9 and 11, and N3 data cells 3, 5, 7, 9, 11 and 13 were executed in isolated read-only extraction scripts. Relevant central line-profile and contribution functions were extracted from N3 cells 22 and 25. A direct H-beta calculation was added to the extraction wrapper on the existing common-line support; no production calibration or notebook was changed. N1 definitions and stored outputs were inspected, while its common loading/selection logic was exercised through the dependent calculations.

The report then used a fresh scalar-data calculation for the whole-system bootstrap and figures, with **10,000 draws and seed 20260911**. Every stage was resampled independently, and a stage's same system draw was used across bins and metrics. The stored outputs preserve support counts and non-finite/zero-draw diagnostics. This is not a rerun of every pooled, map or paired-diagnostic notebook branch.

## 10.3 Independent checks

| Check | Fresh result | What the result establishes |
|---|---:|---|
| Live $T-OI$ identity | Maximum absolute residual $8.72\times10^{-17}$ | The reported SF decomposition follows the implemented estimators |
| Analytic reservoir solution against numerical ODE integration | Maximum absolute difference $4.00\times10^{-11}$ in the demonstration units | Correct solution of the stated constant-coefficient equations |
| Total gas conservation residual | Maximum absolute residual $1.73\times10^{-17}$ | Internal mass bookkeeping in the chosen closure |
| Truncated lognormal NSF intensity against numerical quadrature | Absolute difference $1.10\times10^{-14}$ | Correct closed-form conditional moment |
| Occupancy probability sum | Maximum error $1.11\times10^{-16}$ | The simplified model preserves the three-class partition |
| Two-component fading example | All chosen absolute line amplitudes decline while [N II]/H-alpha increases | The proposed coexistence is mathematically possible |
| Direct H-beta versus proxy stage contrasts | Maximum finite difference below 0.0096 dex on the checked support | The specific denominator shortcut is small here despite its invalid general justification |

These checks establish numerical and definitional consistency. They **do not validate a physical mechanism against MAUVE**, determine model parameters or prove causal evolution. The report figures are either freshly calculated sample summaries or explicitly labelled analytical demonstrations.

## 10.4 Reproduction package

The persistent asset directory is

`/Users/Igniz/Desktop/ICRAR/MAUVE/assets/20260911_sf_nsf_physical_model/`.

It contains the input fingerprints; sample QC; per-system SF and H-alpha scalar tables; stage line profiles including H-beta; denominator checks; bootstrap summaries and contrasts; the analytical calculation and six figures in PNG/PDF; numerical check results; a source/claim ledger; and the report rendering tools and QA record. The lightweight figure/calculation entry point is `reproduce_analysis.py`, operating on the saved scalar extracts. The two `inspect_*.py` wrappers document the earlier live extraction and retain their original temporary-workflow destinations; they are provenance scripts rather than a general pipeline command.

Remaining limitations are the small and heterogeneous stage samples, changing support across bins, no complete source-pipeline rerun, no probabilistic ND reconstruction, no fit of physical parameters, no matched CO/stellar-age modelling, and non-unique ionizing-source templates. Prior report/workflow guidance helped organize the audit; the scientific statements and numerical values above were checked against current local inputs or the cited sources.

# 11. References and prioritized reading

The first papers to read for this specific problem are **Brown et al. (gas content and efficiency), Tomicic et al. 2021a (area versus luminosity), Belfiore et al. 2022 (mixed DIG ionization), Citro et al. (the fading counterexample), and Koppen et al. (analytical stripping)**. The list below groups neither NSF nor all stripped systems into a single physical category.

Access labels: **F** means relevant full-text sections were inspected; **A** means the primary abstract/record was inspected; **E** means the abstract plus indexed full-text excerpts were inspected. These labels do not imply every page was read. Metadata for accepted recent manuscripts are reported as such rather than assigned an unverified final journal date.

1. **Koppen, J., Jachym, P., Taylor, R. & Palous, J. (2018).** *Ram Pressure Stripping Made Easy: An Analytical Approach.* MNRAS, 479, 4367. [Author manuscript](https://arxiv.org/abs/1806.05887); [DOI](https://doi.org/10.1093/mnras/sty1610). **F.** Analytical restoring-force and impulse regimes; the basis for distinguishing pressure amplitude from encounter duration.

2. **Lilly, S. J. et al. (2013).** *Gas-regulation of galaxies: the evolution of the cosmic specific star formation rate, the metallicity-mass-star-formation rate relation and the stellar content of halos.* ApJ, 772, 119. [Author manuscript](https://arxiv.org/abs/1303.5059); [DOI](https://doi.org/10.1088/0004-637X/772/2/119). **A.** Conservation framework used as motivation for the simplified local reservoir closure.

3. **Brown, T. et al. (2023).** *VERTICO VII: Environmental quenching caused by suppression of molecular gas content and star formation efficiency in Virgo cluster galaxies.* [Author manuscript](https://arxiv.org/abs/2308.10943). **F.** Closest resolved Virgo comparison for interpreting lower SFR through both molecular content and efficiency.

4. **Watts, A. B. et al. (2023).** *VERTICO V: The environmentally driven evolution of the inner cold gas discs of Virgo cluster galaxies.* PASA. [Author manuscript](https://arxiv.org/abs/2303.07549); [DOI](https://doi.org/10.1017/pasa.2023.14). **A.** Inner gas changes and preferential depletion of low-surface-density molecular material.

5. **Fossati, M. et al. (2018).** *A Virgo Environmental Survey Tracing Ionised Gas Emission (VESTIGE). II. Constraining the quenching time in the stripped galaxy NGC4330.* A&A, 614, A57. [Author manuscript](https://arxiv.org/abs/1801.09685); [DOI](https://doi.org/10.1051/0004-6361/201732373). **A.** An individual-galaxy example of spatially resolved quenching chronology.

6. **Belfiore, F. et al. (2022).** *A tale of two DIGs: The relative role of H II regions and low-mass hot evolved stars in powering the diffuse ionised gas (DIG) in PHANGS-MUSE galaxies.* A&A, 659, A26. [Author manuscript](https://arxiv.org/abs/2111.14876); [DOI](https://doi.org/10.1051/0004-6361/202141859). **E.** Leakage and evolved-star contributions are spatially dependent; useful for residual-component templates.

7. **Belfiore, F. et al. (2016).** *SDSS IV MaNGA - spatially resolved diagnostic diagrams: a proof that many galaxies are LIERs.* MNRAS, 461, 3111. [DOI](https://doi.org/10.1093/mnras/stw1234). **F.** Resolved low-EW emission and its relation to older stellar populations.

8. **Cid Fernandes, R., Stasinska, G., Mateus, A. & Vale Asari, N. (2011).** *A comprehensive classification of galaxies in the Sloan Digital Sky Survey: how to tell true from fake AGN?* MNRAS, 413, 1687. [DOI](https://doi.org/10.1111/j.1365-2966.2011.18244.x). **A.** EW-based classification and the distinction between weak AGN-like ratios and retired stellar populations.

9. **Zhang, K. et al. (2017).** *SDSS-IV MaNGA: the impact of diffuse ionized gas on emission-line ratios, interpretation of diagnostic diagrams and gas metallicity measurements.* MNRAS, 466, 3217. [Author manuscript](https://arxiv.org/abs/1612.02000); [DOI](https://doi.org/10.1093/mnras/stw3308). **A.** DIG changes line diagnostics and abundance inferences.

10. **Lacerda, E. A. D. et al. (2018).** *Diffuse ionized gas in galaxies across the Hubble sequence at the CALIFA resolution.* MNRAS, 474, 3727. [DOI](https://doi.org/10.1093/mnras/stx3022). **F.** Resolved EW regimes and mixed emission; thresholds depend on the observational setting.

11. **Tomicic, N. et al. (2021a).** *GASP. XXXII. Measuring the diffuse ionized gas fraction in ram-pressure stripped galaxies.* [Author manuscript](https://arxiv.org/abs/2011.08869); [DOI](https://doi.org/10.3847/1538-4357/abca93). **F.** Particularly useful evidence that greater DIG-dominated area need not imply a larger integrated DIG flux fraction.

12. **Tomicic, N. et al. (2021b).** *GASP. XXXV. Characteristics of the diffuse ionized gas in gas-stripped galaxies.* ApJ, 922, 131. [Author manuscript](https://arxiv.org/abs/2108.12433); [DOI](https://doi.org/10.3847/1538-4357/ac230e). **F.** Tail excitation and [O I] excess motivate extra heating in appropriate regions.

13. **Fossati, M. et al. (2016).** *MUSE sneaks a peek at extreme ram-pressure stripping events - II. The physical properties of the gas tail of ESO137-001.* MNRAS, 455, 2028. [DOI](https://doi.org/10.1093/mnras/stv2400). **F.** Photoionized knots and additional excitation in a stripped tail.

14. **Bellhouse, C. et al. (2019).** *GASP. XV. A MUSE view of extreme ram-pressure stripping along the line of sight: physical properties of the jellyfish galaxy JO201.* MNRAS, 485, 1157. [Author manuscript](https://arxiv.org/abs/1902.04486); [DOI](https://doi.org/10.1093/mnras/stz460). **F.** Projection, excitation and multiple kinematic components.

15. **D'Agostino, J. J. et al. (2019).** *Separating line emission from star formation, shocks, and AGN ionization in NGC 1068.* MNRAS, 487, 4153. [Author manuscript](https://arxiv.org/abs/1906.07907); [DOI](https://doi.org/10.1093/mnras/stz1611). **F.** Joint spatial, kinematic and line-ratio decomposition; a methodological reference rather than a Virgo analogue.

16. **Citro, A. et al. (2017).** *A methodology to select galaxies just after the quenching of star formation.* MNRAS, 469, 3108. [Author manuscript](https://arxiv.org/abs/1704.05462); [DOI](https://doi.org/10.1093/mnras/stx932). **F.** A direct counterexample to the claim that [O III] generically fades more slowly than Balmer emission.

17. **Tacchella, S. et al. (2022).** *H-alpha emission in local galaxies: star formation, time variability, and the diffuse ionized gas.* MNRAS, 513, 2904. [DOI](https://doi.org/10.1093/mnras/stac818). **F.** Radiative transfer and time variability complicate the H-alpha-to-SFR mapping.

18. **Attwater, A. et al. (2025).** *MAUVE-MUSE: A star formation-driven outflow caught in the act of quenching the stripped Virgo galaxy NGC4064.* Accepted ApJL manuscript. [Author manuscript](https://arxiv.org/abs/2512.10574). **F.** A directly relevant MAUVE case requiring consideration of an outflow during quenching.

19. **Zhu, J., Tonnesen, S. & Bryan, G. L. (2024).** *When and How Ram Pressure Stripping of Low-mass Satellite Galaxies Enhances Star Formation.* ApJ, 960, 54. [Author manuscript](https://arxiv.org/abs/2309.07037); [DOI](https://doi.org/10.3847/1538-4357/acfe6f). **A.** Pressure-driven redistribution can transiently enhance central SF in simulations.

20. **Goller, J. et al. (2023).** *Jellyfish galaxies in the IllustrisTNG simulations: no enhanced population-wide star formation according to TNG50.* MNRAS, 525, 3551. [Author manuscript](https://arxiv.org/abs/2304.09199); [DOI](https://doi.org/10.1093/mnras/stad2551). **A.** Population selection and individual bursts need not produce the same average trend.

21. **George, K. et al. (2025).** *Star formation at different stages of ram-pressure stripping as seen in far-ultraviolet imaging of 13 GASP galaxies.* A&A, 700, A38. [Author manuscript](https://arxiv.org/abs/2505.15066); [DOI](https://doi.org/10.1051/0004-6361/202554945). **A.** UV/H-alpha comparisons provide additional constraints without a universal extent ordering.

22. **Lora, V. et al. (2026).** *From blue to red spirals: Slow galaxy transformation via ram pressure stripping in TNG-50.* Accepted ApJ manuscript. [Author manuscript](https://arxiv.org/abs/2608.11336). **A.** A recent example of gradual transformation; not a calibration of Virgo stage ages.

23. **Kennicutt, R. C., Jr. & Evans, N. J., II (2012).** *Star Formation in the Milky Way and Nearby Galaxies.* ARA&A, 50, 531. [Author manuscript](https://arxiv.org/abs/1204.3552); [DOI](https://doi.org/10.1146/annurev-astro-081811-125610). **E.** Author review of SFR tracers and their response assumptions.

24. **Kruijssen, J. M. D. & Longmore, S. N. (2014).** *An uncertainty principle for star formation - I. Why galactic star formation relations break down below a certain spatial scale.* MNRAS, 439, 3239. [Author manuscript](https://arxiv.org/abs/1401.4459); [DOI](https://doi.org/10.1093/mnras/stu098). **A.** Small apertures sample different gas and young-star evolutionary phases.

25. **Vulcani, B. et al. (2020).** *GASP. XXX. The spatially resolved SFR-mass relation in stripping galaxies in the local universe.* ApJ, 899, 98. [Author manuscript](https://arxiv.org/abs/2007.04996); [DOI](https://doi.org/10.3847/1538-4357/aba4ae). **A.** An observed enhancement case that limits a universal suppression narrative.

26. **Fujita, Y. & Nagashima, M. (1999).** *Effects of Ram Pressure from the Intracluster Medium on the Star Formation Rate of Disk Galaxies in Clusters of Galaxies.* ApJ, 516, 619. [Author manuscript](https://arxiv.org/abs/astro-ph/9812378); [DOI](https://doi.org/10.1086/307139). **A.** Early analytical treatment of pressure effects on star formation as well as gas removal.
