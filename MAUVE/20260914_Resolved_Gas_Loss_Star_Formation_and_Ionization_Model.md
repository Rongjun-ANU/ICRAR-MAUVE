---
title: "Resolved gas loss, star formation, and continuing ionization in Virgo galaxies"
subtitle: "Observational constraints, analytical derivation, and a partial fit to MAUVE-MUSE"
date: "14 September 2026"
lang: en
---

# 1. What the current spatially resolved observations establish

This report develops a physical interpretation of the changing distribution of star-forming (SF), non-SF-classified but Balmer-detected (NSF), and Balmer-nondetected (ND) regions along the adopted Virgo infall stages. The starting point is the measured spatial distribution and emission of these classes. Gas removal, declining young-star ionization, and continued ionization of retained gas are then introduced as a candidate explanation with explicitly stated assumptions. A new fit constrains the attenuation of the surviving-SF intensity profile; it does not independently determine a stripping rate or an infall timescale.

The observational calculations below were refreshed on 14 September 2026 from the live products and selected data cells of the three 9 September stage-analysis notebooks, together with the explicitly requested **20260914_check_SF_gradient_scan_VIVA.ipynb**. Appendix A identifies the files, execution scope, and reproducible exports. The earlier 11 September report is preserved. Throughout, logarithms denoted by $\log_{10}$ are decimal, whereas $\ln$ denotes the natural logarithm.

## 1.1 Sample, class definitions, and the control assumption

The current stage analyses contain 26 independent data products: eight pre-peak, five close-to-peak, and thirteen post-peak. The combined NGC4567/8 product represents two physical galaxies but receives one product weight and one resampling unit; the inventory therefore represents 27 physical galaxies. The missing post-peak candidates NGC4548, NGC4579, and NGC4689 are not silently counted as nondetections.

The analysis domain is the valid mapped disc with finite stellar mass surface density, a valid gas-bin identifier, the nGIST usable-mask value, and a defined elliptical radius. The stage occupancy and intensity calculations additionally require finite post-fit continuum signal-to-noise ratio **strictly greater than 25**. The working bins have widths 0.25 dex in $\log_{10}\Sigma_*$ or $0.25R_e$ in radius, where $R_e$ is the catalogue effective radius; the displayed windows are $7\leq\log_{10}(\Sigma_*/M_\odot\,\mathrm{kpc}^{-2})<9.5$ and $0\leq r/R_e<2.5$. A galaxy contributes to an occupancy bin when at least 50 retained pixels are present. Category-intensity diagnostics impose their additional 20-pixel support criterion where specified.

The SF class requires a finite HII-selected SFR-map value, $\mathrm{EW}(\mathrm{H}\alpha)>6$ Angstrom, and intrinsic H-alpha dispersion below $45\ \mathrm{km\,s^{-1}}$. The HII map combines the [N II] and [S II] BPT criteria: the [N II] HII mask lies below both the empirical [Kauffmann et al. (2003)](#ref-kauffmann) and theoretical [Kewley et al. (2001)](#ref-kewley) boundaries, and the [S II] HII mask uses the latter theoretical boundary, with the pipeline's domain and quality conditions. The intrinsic dispersion is evaluated from the observed and instrumental variances only where their difference is positive. ND means that the joint H-alpha/H-beta detection requirement fails. NSF is the remaining Balmer-detected population that fails the SF definition. Thus NSF includes several possible excitation and quality-selection outcomes; it is not an independently identified population of diffuse ionized gas, evolved-star-powered regions, shocks, or completely quenched regions. A missing EW or usable intrinsic dispersion can also prevent SF classification.

The raw line-detection conditions are signal-to-noise ratio at least three and flux at least 20 in the stored $10^{-20}\ \mathrm{erg\,s^{-1}\,cm^{-2}}$ units. These emission-line thresholds differ from the continuum $\mathrm{S/N}>25$ condition. The analysis retains the surface-density convention already present in the products and applies no additional inclination correction. Physical model surface densities must be projected consistently before comparison with those products.

**Control assumption.** The pre-peak sample is adopted as a reasonable field-like, normally star-forming control at fixed stellar mass surface density and within the common observed support. This is the requested working approximation. Virgo membership does not prevent its use as an internal control. However, field equivalence is an assumption about its baseline gas/SF state, rather than a consequence of the stage label. Residual environmental processing, morphology, stellar populations, and the particular eight-system composition remain possible systematic offsets. The calculations retain NGC4383 in this stage analysis; no exclusion from another scientific branch is imported.

The three stages are cross-sectional populations. They are not three observations of the same galaxies at known time intervals. Differences between them will be called *stage contrasts*; interpreting those contrasts as evolution requires the control assumption and an additional assumption that the stage samples trace a sufficiently comparable family of histories.

## 1.2 Area occupancy and intensity are different observables

For galaxy $g$, spatial bin $b$, and class $c$, let $N_{g,b}$ be the retained usable-pixel count and $N_{c,g,b}$ the class count. Equal native-pixel areas make

$$
f_{c,g,b}=\frac{N_{c,g,b}}{N_{g,b}},
\qquad
F_{c,s,b}=\frac{1}{N_{\mathrm{gal},s,b}}
\sum_{g\in(s,b)}f_{c,g,b}.
\tag{1}
$$

Here $s$ denotes stage, and $N_{\mathrm{gal},s,b}$ counts eligible galaxy products in that stage and bin. The class fractions sum to unity on this common domain. $F_c$ measures the fraction of retained observed area in a class. It does not measure a gas mass fraction or the fraction of the entire unobserved galaxy.

For SF pixels, define the galaxy's category-conditional mean SFR surface density $I_{\mathrm{SF},g,b}$. The stage notebook constructs

$$
O_{s,b}=F_{\mathrm{SF},s,b},
\qquad
T_{s,b}=\left\langle f_{\mathrm{SF},g,b}I_{\mathrm{SF},g,b}\right\rangle_g,
\qquad
I_{\mathrm{SF},s,b}=\frac{T_{s,b}}{O_{s,b}}.
\tag{2}
$$

The brackets denote the arithmetic mean over eligible galaxy products, including a zero SF numerator where no SF pixels survive. Consequently,

$$
T_{s,b}=O_{s,b}I_{\mathrm{SF},s,b},
\qquad
\Delta\log_{10}T=\Delta\log_{10}O+
\Delta\log_{10}I_{\mathrm{SF}}.
\tag{3}
$$

This is an exact estimator identity. The resulting $I_{\mathrm{SF},s,b}$ is an occupancy-weighted mean of the individual-galaxy SF intensities, rather than an arithmetic mean giving each surviving SF population equal weight. $T$ is retained-area-normalized emission attributed to the selected SF component; it is not a direct estimate of all true star formation in ND and NSF pixels.

For NSF emission, the existing notebook separately reports an arithmetic mean of within-category H-alpha luminosity surface densities and an occupancy-weighted contribution:

$$
I_{\alpha,\mathrm{NSF},s,b}
=\left\langle I_{\alpha,\mathrm{NSF},g,b}\right\rangle_{\mathrm{supported}\ g},
\qquad
J_{\alpha,\mathrm{NSF},s,b}
=\left\langle f_{\mathrm{NSF},g,b}
I_{\alpha,\mathrm{NSF},g,b}\right\rangle_g.
\tag{4}
$$

$J_{\alpha,\mathrm{NSF}}$ is defined relative to the same retained-area denominator. In general it is not $F_{\mathrm{NSF}}I_{\alpha,\mathrm{NSF}}$: the galaxy subsets can differ, and fraction and intensity can covary within galaxies. An NSF-area increase therefore cannot establish an NSF-luminosity increase. Nor can the H-alpha luminosity of NSF gas be converted to SFR without identifying its ionizing sources.

## 1.3 The robust pattern is spatial redistribution plus suppression

The observations support two distinct changes: a marked loss of SF occupancy in the outskirts/low-$\Sigma_*$ domain, and a lower intensity in much of the remaining SF population. They also support increasing inner NSF occupancy, with substantial object-to-object variation and non-monotonic behaviour in some close-to-peak bins.

![Figure 1. Refreshed equal-galaxy class-occupancy profiles. The horizontal coordinates are stellar mass surface density and elliptical radius; these are complementary projections, not interchangeable environmental coordinates. Shading denotes 16th-84th percentile whole-galaxy bootstrap intervals. Support changes with position and stage.](assets/20260914_resolved_RPS_academic_model/figure_01_live_occupancy.png)

Several examples quantify the result. Intervals below are 16th-84th percentile whole-galaxy bootstrap ranges, resampling both the pre-peak reference and the comparison stage. They are descriptive uncertainty intervals, not multiple-testing-corrected detections.

| Location | Observed pre-peak to post-peak change | Interpretation within the measured support |
|---|---|---|
| $r/R_e=1.625$ | $F_{\mathrm{SF}}:0.429\rightarrow0.027$; $F_{\mathrm{ND}}:0.317\rightarrow0.789$ | Strong outer loss of detectable SF area and growth of the ND class. |
| $r/R_e=0.125$ | $F_{\mathrm{NSF}}:0.472\rightarrow0.694$; $\Delta\log_{10}F_{\mathrm{NSF}}=0.167\ [0.060,0.296]$ | Growing central area that remains Balmer-detected but fails the SF criteria. |
| $\log_{10}\Sigma_*=7.625$ | $\Delta\log_{10}O=-1.011$; $\Delta\log_{10}I_{\mathrm{SF}}=-0.429\ [-0.523,-0.291]$ | Occupancy loss and weaker remaining SF both contribute. |
| $\log_{10}\Sigma_*=8.625$ | $\Delta\log_{10}I_{\mathrm{SF}}=-0.408\ [-0.529,-0.263]$ | Suppression persists in the denser surviving-SF domain. |
| $r/R_e=1.625$ | $\Delta\log_{10}I_{\mathrm{SF}}=+0.105\ [-0.213,0.345]$ | The sparse outer survivors do not establish uniform intensity suppression. |

Surface-density coordinates in this and subsequent tables use $M_\odot\,\mathrm{kpc}^{-2}$. At $\log_{10}\Sigma_*=7.625$, equation (3) gives $\Delta\log_{10}T=-1.440$, combining the approximately 1.01-dex occupancy decline with the 0.43-dex intensity decline. A model that reduces only the number of SF regions, while leaving the surviving population unchanged, misses this second effect. Conversely, a model that uniformly dims every region without allowing detection and classification transitions misses the spatial occupancy pattern.

![Figure 2. Decomposition of the stage changes into SF occupancy, intensity within surviving SF area, and their exact product, together with NSF H-alpha intensity and contribution diagnostics. A reduced selected-area intensity is distinct from a reduced number of selected regions. The NSF contribution must be computed from galaxy-level products as in equation (4).](assets/20260914_resolved_RPS_academic_model/figure_02_live_decomposition.png)

An additional observational constraint should be included in the proposed narrative: **inner NSF regions become more common without necessarily becoming brighter**. At $\log_{10}\Sigma_*=8.625$, the post/pre contrasts are $+0.080\ [-0.039,0.223]$ dex in NSF fraction, $-0.508\ [-0.721,-0.150]$ dex in NSF conditional H-alpha intensity, and $-0.351\ [-0.566,-0.020]$ dex in its retained-area contribution. The fraction increase in that particular bin has an interval crossing zero, whereas higher-density and central radial bins more clearly support the occupancy increase. At $r/R_e=0.625$, the corresponding NSF contribution interval includes zero. The report therefore does not generalize a luminosity decrease to every inner bin.

## 1.4 Differential fading is an interpretation to test

The third proposed observational statement needs a distinction between a measurement and its physical explanation. The measurements are line-specific changes in luminosity and line ratio, together with changing category membership. Differential decline of emission components is a plausible explanation; it is not yet established as the cause of the NSF-area increase.

For a fixed region and a line $\ell$ with Balmer denominator $B$, the following identity is exact:

$$
\frac{d}{dt}\ln\!\left(\frac{\mathcal L_\ell}{\mathcal L_B}\right)
=\frac{1}{\mathcal L_\ell}\frac{d\mathcal L_\ell}{dt}
-\frac{1}{\mathcal L_B}\frac{d\mathcal L_B}{dt}.
\tag{5}
$$

Here $\mathcal L_\ell$ is luminosity per area. A ratio rises if its numerator declines more slowly in *fractional* terms than its denominator; the numerator can itself become fainter. However, comparing different galaxies, different class-selected pixels, or different subsets of line detections is not a measurement of the time derivatives in equation (5).

The refreshed NSF diagnostic uses a common detected-line mask, with direct corrected H-beta for [O III]/H-beta. It forms ratios from the equal-galaxy mean line intensities, rather than a median of spaxel ratios. Its post/pre contrasts illustrate why the forbidden lines must not be combined into a single claim:

| Location | $\Delta\log_{10}([\mathrm{N\,II}]/\mathrm{H}\alpha)$ | $\Delta\log_{10}([\mathrm{S\,II}]/\mathrm{H}\alpha)$ | $\Delta\log_{10}([\mathrm{O\,III}]/\mathrm{H}\beta)$ |
|---|---|---|---|
| $\log_{10}\Sigma_*=8.625$ | $+0.181\ [-0.053,0.288]$ | $+0.043\ [0.002,0.074]$ | $-0.077\ [-0.133,0.034]$ |
| $r/R_e=0.625$ | $+0.428\ [0.143,0.503]$ | $+0.184\ [0.076,0.219]$ | $+0.084\ [0.022,0.231]$ |

Thus enhanced low-ionization ratios are a useful part of the interpretation, while a universal increase in every BPT ratio is unsupported. A rightward BPT displacement need not also be upward. Moreover, NSF is partly defined by the diagnostic boundaries: finding different line ratios in SF and NSF populations is partly built into the classification. Absolute within-NSF stage contrasts and NSF-minus-SF contrasts within a galaxy answer different questions.

![Figure 3. Refreshed line-intensity contrasts and ratios for NSF regions on common detected-line support. [O III]/H-beta is calculated using H-beta itself. The line-specific results constrain the allowable ionizing-source mixture; they do not establish a common forbidden-line fading time.](assets/20260914_resolved_RPS_academic_model/figure_03_live_line_ratios.png)

The H-beta distinction also matters technically. The attenuation code clips inferred reddening to zero for observed Balmer decrements below 2.86; it does not force their corrected H-alpha/H-beta ratios to 2.86. Consequently, replacing H-beta with a fixed multiple of H-alpha is not an exact identity. The extra direct-H-beta calculation in this report leaves the notebooks unchanged. This correction does not create a universal positive [O III]/H-beta trend.

A more defensible statement of the third finding is: *the evolving NSF population shows line-dependent excitation changes and declining H-alpha intensity in several important inner bins; these motivate a model in which the young-star contribution fades relative to other ionization sources, subject to metallicity, selection, and excitation degeneracies.*

## 1.5 NGC4654: a localized enhancement candidate

The specifically named 14 September notebook was inspected and its data/plane-estimation cells rerun. Its title and some explanatory text retain earlier dates, and it records an unexecuted plotting-layout change. Those textual dates do not alter which artifact was used here. This report relies on fresh numerical extraction and a new report figure rather than treating all stored notebook figures as freshly validated.

The notebook defines a local SF residual relative to each galaxy's own median SFR-versus-$\Sigma_*$ relation, and an environmental residual relative to the median of the pre-peak galaxy medians. It fits a plane in projected east/north coordinates:

$$
\delta_{\mathrm{local}}(E,N)
=\beta_0+\beta_E E+\beta_N N,
\qquad
G_{\mathrm{SF}}=(\beta_E^2+\beta_N^2)^{1/2},
\qquad
\theta_{\mathrm{SF}}=\operatorname{atan2}(\beta_E,\beta_N).
\tag{6}
$$

$E,N$ are projected offsets in kpc, $\beta_0$ is a residual in dex, and the two slopes have units dex per kpc. Position angle is east of north, expressed modulo $360^\circ$. The fitted local plane determines the direction without consuming the VIVA angle. Environmental hemisphere medians then describe whether either side is above the adopted control.

For NGC4654, the local gradient has $G_{\mathrm{SF}}=0.02134\ \mathrm{dex\,kpc^{-1}}$ and $\theta_{\mathrm{SF}}=319.33^\circ$, separated by $4.33^\circ$ from the adopted northwest compression direction, $315^\circ$. The local plane explains $R^2=0.0444$ of the pixel-level local-residual variance. There are 342,926 pixels on the joint SF/residual mask; this count is area support, not the number of independent measurements.

| Diagnostic | At the inferred $319.33^\circ$ direction | At the fixed VIVA $315^\circ$ direction |
|---|---|---|
| Environmental median, facing side | $+0.0354$ dex | $+0.0357$ dex |
| Environmental median, opposite side | $-0.0957$ dex | $-0.0943$ dex |
| Difference between the two medians | $0.1311$ dex | $0.1300$ dex |
| Fraction with standardized residual $z_{\mathrm{env}}>1$, facing side | 22.83% | 22.71% |
| Same fraction, opposite side | 14.06% | 14.25% |

The facing-side median is approximately 8.5% above the control, while the opposite-side median is below it. The 0.131-dex side contrast corresponds to a ratio of approximately 1.35 between the residual-normalized medians; it is not a 35% enhancement relative to the field. The positive tail is more populated on the facing side. Here $z_{\mathrm{env}}$ uses the width of the pre-peak residual distribution, not the standard error on the reference and not a Gaussian detection significance.

![Figure 4. NGC4654 on the freshly reconstructed joint SF mask. The two maps show within-galaxy and environmental residuals, with north upward and east to the left. Cyan denotes the inferred gradient direction and amber the fixed VIVA direction; their lengths are display guides. The cumulative distributions use the fixed VIVA hemispheres. The modest positive facing-side median and its larger positive tail support a candidate local enhancement.](assets/20260914_resolved_RPS_academic_model/figure_04_ngc4654.png)

The notebook calculates no formal gradient significance, position-angle uncertainty, or calibrated spatial-correlation-aware null distribution. Its fixed $\pm45^\circ$ VIVA tolerance is a comparison sector, not a confidence interval. Its gradient mask also does **not** apply the stage notebooks' additional continuum $\mathrm{S/N}>25$ cut. These analyses therefore cannot be treated as a single identical selection.

NGC4654 provides a physically interesting spatial coincidence, especially given the independent gas-compression literature discussed below. It remains a candidate enhancement because footprint structure, spiral arms, stellar-density/radius coupling, the galaxy's tidal history, and uncertainty in the control have not been jointly removed. A coherent plane is a descriptive approximation to a much more structured map.

# 2. What previous studies contribute

The literature review is targeted to the observational constraints and the model ingredients. It is not a systematic census. Primary papers and author manuscripts are cited directly, with the depth of source inspection recorded in the supporting ledgers. Numerical results from other surveys retain their original spatial scale and sample context.

## 2.1 Gas removal and suppression within the surviving disc

[Brown et al. (2023)](#ref-brown) provide a particularly close comparison: VERTICO measurements at 720-pc resolution connect reduced resolved SFR in H I-poor Virgo galaxies with changes in both molecular content and molecular-gas star-formation efficiency. Their early-stage outskirts also show enhancement associated with increased molecular surface density, without a corresponding efficiency increase. These results motivate allowing both $\Sigma_{\mathrm{H_2}}$ and $\tau_{\mathrm{dep}}$ to change; optical SFR alone cannot separate them. [Watts et al. (2023)](#ref-watts) extend the environmental gas-depletion discussion into the inner cold-gas disc. Retention of an inner molecular reservoir therefore need not mean preservation of its original state.

The spatially resolved quenching analysis of NGC4330 by [Fossati et al. (2018)](#ref-fossati) shows how photometry and spectroscopy can provide a chronology beyond an H-alpha snapshot. Its importance here is methodological: a quenching clock requires temporal information and a specified star-formation-history model. It cannot be transferred automatically to the three MAUVE stage labels.

[Koppen et al. (2018)](#ref-koppen) distinguish a long-duration stripping regime controlled by restoring forces from a short-pulse regime controlled by momentum transfer. [Singh et al. (2019)](#ref-singh) derive stripping-radius and stripped-mass estimates for idealized discs. Together they justify a radially selective gas-loss model, but not an uncalibrated conversion of optical ND fraction into stripped gas mass.

## 2.2 Compression can enhance SF, but its sign is conditional

The most directly relevant NGC4654 study is [Lizee et al. (2021)](#ref-lizee). Their deeper CO/H I analysis and analytical/dynamical modelling identify a compressed northwestern region with inferred molecular-gas efficiency about 1.5-2 times higher, conditional on the CO conversion and model assumptions. Their efficiency map uses 12-arcsecond resolution, approximately 1 kpc at their adopted 17-Mpc distance; this is not an individual approximately 100-pc MAUVE element. The result supplies independent physical motivation for examining the northwest optical residuals. [Vollmer (2003)](#ref-vollmer03) also finds that combined gravitational interaction and ram pressure explain NGC4654 more successfully than either process alone in the tested models. Its asymmetry cannot therefore be assigned uniquely to ram pressure from directional agreement.

Compression does not universally shorten the molecular depletion time. [Vollmer et al. (2012)](#ref-vollmer12) find broadly field-like molecular efficiencies and possible windward molecular-fraction increases in several Virgo galaxies, without a significant general windward molecular-efficiency enhancement. [Nehlig et al. (2016)](#ref-nehlig) identify compressed regions that include low molecular efficiency and discuss the competition between compression and turbulent support. A larger molecular fraction, higher SFR per total gas mass, and higher SFR per molecular mass are three distinct statements.

The pressure-dependent framework of [Blitz & Rosolowsky (2006)](#ref-br06) relates molecular-to-atomic ratio to estimated hydrostatic ISM pressure. It is an empirical state relation, not a molecular formation-rate measurement. Its discussion anticipates breakdown around individual-cloud scales. The thermal/dynamical equilibrium model of [Ostriker, McKee & Leroy (2010)](#ref-oml10) similarly describes averages over cloud populations; its bound-cloud reservoir is not identical to all H2. These models inform the physical possibilities without being imposed as exact laws at each approximately 100-pc MAUVE resolution element.

Earlier pressure/SF modelling by [Fujita & Nagashima (1999)](#ref-fujita) and the jellyfish-galaxy model of [Safarzadeh & Loeb (2019)](#ref-safarzadeh) explicitly consider competition between pressure-induced SF and loss of cold gas. The latter uses prescribed cloud efficiency and gas-loss behaviour rather than solving the local class-occupancy problem. Simulations by [Lee et al. (2020)](#ref-lee) produce enhancement or suppression under different wind conditions with strong stellar feedback. These studies argue for a model that admits a localized positive response while its broader gas/SF population declines.

## 2.3 Ionization can outlast local bright HII-region emission

[Belfiore et al. (2022)](#ref-belfiore) model PHANGS-MUSE diffuse emission using leakage from HII regions and contributions from hot evolved stars; their roles vary with location and line species. [Zhang et al. (2017)](#ref-zhang) show that DIG contamination changes low-ionization ratios and diagnostic-diagram positions. These findings support line-specific mixture modelling, including the possibility that an energetically subdominant component matters strongly for particular lines.

[Tomicic et al. (2021a)](#ref-tomicic32) compare stripped and control GASP galaxies and warn that a single EW, surface-brightness, or line-ratio threshold does not universally isolate DIG. Their inferred DIG luminosity fraction is not the same observable as MAUVE NSF area occupancy. [Tomicic et al. (2021b)](#ref-tomicic35) discuss excitation in stripped gas, including evidence for additional heating in tails. Such work motivates shocks as a possible component, not a compulsory explanation for every inner NSF pixel.

An essential counterexample to generic forbidden-line persistence is [Citro et al. (2017)](#ref-citro): in their quenching photoionization models, high-ionization emission can decline rapidly as the ionizing spectrum softens. Thus a model that merely stops young-star formation does not automatically predict that all forbidden/Balmer ratios rise. The relevant hypotheses are changes in source mixtures, spectral shape, ionization parameter, gas structure, and abundance, each with line-dependent consequences.

## 2.4 Relation to the existing MAUVE regulator model

[Huang et al. (2026)](#ref-huang), the user's MAUVE-MUSE paper, defines the local gas reservoir as molecular and writes its SFR, supply, removal, and metallicity budgets in surface-density form. Its equations (16), (17), and (21) supply the notation used here:
$\Sigma_{\mathrm{SFR}}=\Sigma_{\mathrm{H_2}}/\tau_{\mathrm{dep}}$,
$\Sigma_\Phi$ is a supply-rate surface density, and
$\tau_\Phi=\Sigma_{\mathrm{H_2}}/\Sigma_\Phi$.
This report extends that regulator language to environmental loss and an optical observation operator; it does not reinterpret $\Sigma_\Phi$ as a timescale.

The present extension also distinguishes effective supply to a molecular reservoir from gross H I-to-H2 chemistry. Local transport can deliver gas that is already molecular, while formation and dissociation can occur simultaneously. An atomic-reservoir closure will be introduced only as an explicitly restricted realization. The general conservation framework is related to the regulator approach of [Lilly et al. (2013)](#ref-lilly).

# 3. From the measurements to requirements on a physical scenario

The candidate scenario must satisfy several constraints simultaneously. The following table separates a required model behaviour from the mechanism that might supply it.

| Observational constraint | Required behaviour of a successful forward model | Physically motivated possibilities |
|---|---|---|
| Outer ND growth and SF-area loss | More retained observing elements fail joint Balmer detection, with a stronger change in outer/low-density support | Preferential removal of weakly bound gas; depletion of the remaining fuel; declining young-star emission below the local sensitivity limit |
| Lower intensity among selected SF regions | Surviving regions must also have lower apparent H-alpha SFR after applying the selection | Reduced molecular surface density, longer depletion time, or both; possible luminosity-response and selection effects |
| Inner NSF occupancy increases while its H-alpha often fades | Gas and ionization remain detectable after more regions fail the SF cuts; extra NSF area need not add net luminosity | Declining HII dominance, continuing diffuse/old-star/shock ionization, lower EW, changing line width, and selection transitions |
| Low-ionization ratios often increase, while [O III]/H-beta is not uniformly higher | The model must predict each numerator and denominator, including regions that leave the detected-line sample | Evolving mixtures and ionization conditions; metallicity and spectral-shape changes |
| NGC4654's northwest positive residuals | A spatially limited increase or relative preservation of SF must coexist with a declining ensemble | Local compression, increased molecular supply, shorter depletion time, transport; tidal compression remains an alternative |
| Strong galaxy scatter and differing spatial support | The population prediction must preserve galaxy weights and individual footprints | A distribution of gas structure, inclination, orbital history, baseline SF, and disturbance strength |

Two additional requirements follow from this list. First, a *retained-gas* model needs an ionizing-energy budget: retaining atoms alone does not sustain optical emission. Second, an explanation of SF-to-NSF transitions must reproduce the same BPT/EW/width/detection operator that creates the classes. A simplified analytical classifier can clarify the mechanism, but it is not automatically a calibrated implementation of the measured categories.

The proposed physical ordering is therefore conditional: environmental processing removes or disrupts the readily stripped gas reservoir, reducing supply and/or molecular content; surviving young-star emission declines; retained inner gas can continue to absorb ionizing photons or dissipate mechanical energy; its composite spectrum and EW then cross the observational class boundaries. Localized compression can interrupt this ordering in selected regions. The equations below make each link explicit.

# 4. A spatially resolved analytical model

## 4.1 Coordinates, notation, and explicit assumptions

Let $\boldsymbol{x}=(r,\varphi)$ label a location in the galaxy disc, where $r$ is cylindrical radius and $\varphi$ is azimuth. Time $t=0$ denotes the onset of a specified environmental perturbation at that location, rather than the catalogue pre-peak stage. All coefficients may depend on $\boldsymbol{x}$. Spatial arguments are suppressed when solving a local equation. The ordinary derivative $d/dt$ below denotes the evolution of the area-averaged quantity in the stated local aperture, with transport included explicitly or neglected by assumption.

The model is deliberately layered. Conservation equations constrain the gas reservoir. Additional physical closures specify supply and removal. Stellar-population response converts recent SF into young-star-powered luminosity. Continuing ionization supplies other line components. Finally, the observational operator assigns SF, NSF, and ND classes. Only the last layer is directly comparable to the optical class statistics.

| Symbol | Definition and units |
|---|---|
| $\Sigma_*,\Sigma_{\mathrm{HI}},\Sigma_{\mathrm{H_2}}$ | Stellar, atomic, and molecular mass per disc area; $M_\odot\,\mathrm{pc}^{-2}$ or explicitly converted $M_\odot\,\mathrm{kpc}^{-2}$ |
| $\Sigma_{\mathrm{SFR}}$ | True stellar mass formed per area and time; $M_\odot\,\mathrm{yr}^{-1}\,\mathrm{kpc}^{-2}$ |
| $\Sigma_\Phi,\Sigma_{\mathrm{out}}$ | Molecular-reservoir supply and removal rates per area; the same dimensional units as $\Sigma_{\mathrm{SFR}}$ |
| $\tau_{\mathrm{dep}}$ | Molecular depletion time, $\Sigma_{\mathrm{H_2}}/\Sigma_{\mathrm{SFR}}$, after consistent area/time conversion |
| $\tau_\Phi$ | Supply time, $\Sigma_{\mathrm{H_2}}/\Sigma_\Phi$, defined for positive supply |
| $\tau_{\mathrm{conv}}$ | Effective atomic-reservoir transfer time in the restricted two-reservoir closure; not a measured chemical formation time |
| $R,\eta$ | Dimensionless prompt effective return fraction and molecular-reservoir feedback mass-loading factor |
| $k_{\mathrm{HI}},k_{\mathrm{H_2}}$ | Fractional gas-removal coefficients; inverse time |
| $\rho_{\mathrm{ICM}},v_{\mathrm{rel}},P_{\mathrm{ram}}$ | Ambient density, relative speed, and incident dynamical pressure; cgs units in force calculations |
| $\Phi_{\mathrm{grav}},G$ | Gravitational potential per unit mass and Newton's constant; $\Phi_{\mathrm{grav}}$ is distinct from the supply symbol $\Sigma_\Phi$ |
| $\mathcal L_\ell$ | Luminosity in emission line $\ell$ per area; $\mathrm{erg\,s^{-1}\,kpc^{-2}}$ in observational comparisons |
| $\mathcal C_{\lambda,\alpha}$ | Stellar continuum luminosity density per area at H-alpha; $\mathrm{erg\,s^{-1}\,kpc^{-2}\,Angstrom^{-1}}$ |
| $K_\alpha(a),\tau_{\mathrm{ion}}$ | Normalized young-star H-alpha response kernel versus stellar age $a$, and its illustrative exponential response time |
| $C_\alpha$ | H-alpha SFR conversion coefficient, with units $M_\odot\,\mathrm{yr^{-1}}/(\mathrm{erg\,s^{-1}})$ |
| $w_\alpha,\mathcal R_{\ell/B}$ | Continuing-component H-alpha fraction and a dimensionless line ratio |
| $\mu_\alpha,s_\alpha,\mathcal L_{\mathrm{ref}}$ | Mean and dispersion of the assumed log-luminosity distribution, and a reference luminosity surface density making its logarithm dimensionless |
| $\mathcal E$ | Dimensionless integrated logarithmic SF attenuation; it is not an elapsed time |

Gas surface densities in the illustrative calculations share a single associated-helium convention. They can be understood as atomic-associated and molecular-associated total gas masses, rather than counting a helium correction twice. Changing that convention consistently changes inferred depletion times but not the algebraic conservation structure.

The local analytic solution assumes: (i) fixed aperture area; (ii) constant coefficients during each solved interval; (iii) approximately constant $\Sigma_*$ over that interval; (iv) nonnegative reservoirs and removal coefficients; and (v) explicitly specified initial conditions. The instantaneous-recycling approximation treats $R$ as the fraction promptly available to the modelled reservoir. Taking the usual stellar return fraction as that molecular-reservoir return is an additional coarse-grained approximation: actual ejecta need not return directly to H2. A different phase assignment requires extra terms, rather than silently altering the meaning of supply.

Stellar mass can be considered fixed only if $(1-R)\int\Sigma_{\mathrm{SFR}}\,dt\ll\Sigma_*$. The hundred-Myr illustrative trajectories are not an assertion that every approximately 100-pc aperture is a closed, stationary object for that duration. Disc rotation, transport, cloud lifecycles, and evolving apertures limit that interpretation. Population averages can remain useful even when individual apertures exchange material.

The requested notation distinction is dimensional:

$$
[\Sigma_\Phi]=\frac{\mathrm{mass}}{\mathrm{area}\,\mathrm{time}},
\qquad
\tau_\Phi=\frac{\Sigma_{\mathrm{H_2}}}{\Sigma_\Phi},
\qquad
[\tau_\Phi]=\mathrm{time}.
\tag{7}
$$

Equation (7) preserves the definition in [Huang et al. (2026), equation (21)](#ref-huang). The inequality $\tau_\Phi\leq\tau_{\mathrm{dep}}$ printed there is not imposed on the present declining-reservoir extension. Under interrupted supply, $\Sigma_\Phi$ can approach zero and $\tau_\Phi$ can become arbitrarily long. That is a necessary change of regime, not a redefinition of $\Sigma_\Phi$.

## 4.2 Why gas loss is expected to depend on position

**Step 1: incident forcing.** In the usual idealized ram-pressure description,

$$
P_{\mathrm{ram}}(t)=
\rho_{\mathrm{ICM}}(t)\,v_{\mathrm{rel}}^2(t).
\tag{8}
$$

This pressure scale is used in the analytical stripping treatments of [Koppen et al. (2018)](#ref-koppen) and [Singh et al. (2019)](#ref-singh). It is not automatically the isotropic pressure confining a molecular cloud. The relevant wind stress depends on orientation, shielding, and the gas geometry.

**Step 2: restoring force.** For a gas layer with surface density $\Sigma_{\mathrm{gas}}$ in a specified potential, a local force-per-area scale is

$$
P_{\mathrm{rest}}(r)=
\Sigma_{\mathrm{gas}}(r)
\max_z\left|\frac{\partial\Phi_{\mathrm{grav}}(r,z)}{\partial z}\right|
\ \simeq\
2\pi G\,\Sigma_*(r)\Sigma_{\mathrm{gas}}(r).
\tag{9}
$$

The last expression is the stellar-sheet approximation, not a full finite-thickness, bulge-plus-halo calculation. The gas column here is the column accelerated by the wind; it must not be assumed to equal all molecular gas in a resolution element. Equation (9) follows the force-balance approximation used by [Singh et al. (2019), equation (3)](#ref-singh), within the dynamical qualifications emphasized by [Koppen et al. (2018)](#ref-koppen).

**Step 3: radial ordering.** Suppose, only for an analytic illustration, that $\Sigma_*=\Sigma_{*,0}\exp(-r/h_*)$ and $\Sigma_{\mathrm{gas}}=\Sigma_{\mathrm{gas},0}\exp(-r/h_{\mathrm{gas}})$, where $h_*$ and $h_{\mathrm{gas}}$ are positive scale lengths. Substitution into equation (9) gives

$$
P_{\mathrm{rest}}(r)=P_{\mathrm{rest},0}
\exp\!\left[-r\left(\frac{1}{h_*}+\frac{1}{h_{\mathrm{gas}}}\right)\right],
\qquad
P_{\mathrm{rest},0}=2\pi G\Sigma_{*,0}\Sigma_{\mathrm{gas},0}.
\tag{10}
$$

Equating the applied effective normal stress $P_\perp$ to equation (10), taking logarithms, and isolating $r$ yields

$$
r_{\mathrm{strip}}=
\frac{\ln(P_{\mathrm{rest},0}/P_\perp)}
{h_*^{-1}+h_{\mathrm{gas}}^{-1}}.
\tag{11}
$$

This is the exponential-disc construction also used by [Singh et al. (2019), equation (18)](#ref-singh). It applies only within the assumed disc extent and force-balance regime. A negative formal radius means that even the central restoring scale is exceeded in this approximation; an infinite or very large radius does not imply infinite physical disc size.

**Step 4: duration matters.** For a layer accelerated vertically,

$$
\frac{dv_z}{dt}=
\frac{P_\perp(t)}{\Sigma_{\mathrm{gas}}}-g_z,
\qquad
\Delta v_z\simeq
\frac{1}{\Sigma_{\mathrm{gas}}}\int P_\perp(t)\,dt
\quad\text{(short-pulse limit)}.
\tag{12}
$$

$v_z$ is vertical velocity and $g_z$ the opposing gravitational acceleration. The short-pulse expression neglects gravitational impulse during the pulse and assumes a nearly fixed accelerated column. Escape then requires comparison with the appropriate escape speed. It illustrates the impulse regime discussed by [Koppen et al. (2018)](#ref-koppen), rather than replacing orbital dynamics with a universal pressure threshold.

Equations (8)-(12) motivate larger loss coefficients in weakly bound outer gas. They do not determine $k_{\mathrm{HI}}(\boldsymbol{x},t)$ or $k_{\mathrm{H_2}}(\boldsymbol{x},t)$ from the optical data. A truncation radius describes a gas response; the optical ND boundary also depends on subsequent SF fading, ionization, and detectability.

## 4.3 Molecular conservation and the meaning of supply

**Step 1: define the aperture budget.** For a fixed disc-plane aperture of area $\mathcal A$, the vertically integrated continuity equation has the integral form

$$
\frac{d\overline{\Sigma}_{\mathrm{H_2}}}{dt}
=
\frac{1}{\mathcal A}\int_{\mathcal A}S_{\mathrm{H_2}}\,dA
-\frac{1}{\mathcal A}\oint_{\partial\mathcal A}
\Sigma_{\mathrm{H_2}}\,\boldsymbol v_{\mathrm{H_2}}
\cdot\boldsymbol n\,dl.
\tag{13}
$$

The overbar denotes an aperture average; $S_{\mathrm{H_2}}$ is the local net volume-integrated phase/source term per area and time; $\boldsymbol v_{\mathrm{H_2}}$ is the in-plane molecular-gas velocity; and $\boldsymbol n$ is the outward boundary normal. The second term is signed net transport. Equation (13) is mass conservation written for the chosen control area. Dropping it is an assumption about transport, not a consequence of spatial resolution.

**Step 2: group physically distinct gains and losses.** A useful bookkeeping convention is

$$
\begin{aligned}
\Sigma_\Phi&=\Sigma_{\mathrm{form,H_2}}+
\Sigma_{\mathrm{in,H_2}},\\
\Sigma_{\mathrm{out}}&=
\Sigma_{\mathrm{diss,H_2}}+
\Sigma_{\mathrm{boundary,out,H_2}}+
\eta\Sigma_{\mathrm{SFR}}+
\Sigma_{\mathrm{strip,H_2}}.
\end{aligned}
\tag{14}
$$

Every term on the right has units mass per area per time. The first line includes gross formation into H2 and inward delivery of material already molecular. The second separates dissociation, outward transport, feedback removal, and stripping. Boundary gains and losses are the inward and outward parts of equation (13), not extra copies of its net term. One may instead group net phase conversion into a signed supply term, but must state that alternative convention. The restricted positive-supply solution below uses equation (14) with dissociation and additional transport omitted.

**Step 3: recover and extend the local regulator.** With effective prompt return fraction $R$, the reservoir budget becomes

$$
\begin{aligned}
\frac{d\Sigma_{\mathrm{H_2}}}{dt}
&=\Sigma_\Phi-(1-R)\Sigma_{\mathrm{SFR}}-\Sigma_{\mathrm{out}},\\
\Sigma_{\mathrm{SFR}}&=\frac{\Sigma_{\mathrm{H_2}}}{\tau_{\mathrm{dep}}}.
\end{aligned}
\tag{15}
$$

This is the molecular-reservoir form of [Huang et al. (2026), equations (16)-(17)](#ref-huang), related to [Lilly et al. (2013)](#ref-lilly). Splitting $\Sigma_{\mathrm{out}}$ into feedback and environmental terms does not add a second total-removal term. In particular, if $\Sigma_{\mathrm{out}}=\eta\Sigma_{\mathrm{SFR}}+k_{\mathrm{H_2}}\Sigma_{\mathrm{H_2}}$, other phase/transport losses are assumed absent or explicitly absorbed into the stated effective coefficients.

**Step 4: derive the condition for suppression.** Divide the first line of equation (15) by $\Sigma_{\mathrm{H_2}}$, substitute the above removal closure, and differentiate $\ln\Sigma_{\mathrm{SFR}}=\ln\Sigma_{\mathrm{H_2}}-\ln\tau_{\mathrm{dep}}$. For positive reservoirs,

$$
\frac{d\ln\Sigma_{\mathrm{SFR}}}{dt}
=
\frac{\Sigma_\Phi}{\Sigma_{\mathrm{H_2}}}
-\frac{1-R+\eta}{\tau_{\mathrm{dep}}}
-k_{\mathrm{H_2}}
-\frac{1}{\tau_{\mathrm{dep}}}\frac{d\tau_{\mathrm{dep}}}{dt}.
\tag{16}
$$

This derivation exposes four distinct contributions: supply, consumption/feedback, direct molecular loss, and changes in molecular depletion time. Gas supply is not the same parameter as efficiency. A retained molecular reservoir can have declining SFR through a longer $\tau_{\mathrm{dep}}$, and molecular loss can suppress SFR even if $\tau_{\mathrm{dep}}$ stays constant.

Define the integrand $\Gamma_{\mathrm{SFR}}$ as the negative of the right-hand side of equation (16). Integration gives

$$
\mathcal E(\boldsymbol{x},t)=\int_0^t
\Gamma_{\mathrm{SFR}}(\boldsymbol{x},t')\,dt',
\qquad
\frac{\Sigma_{\mathrm{SFR}}(\boldsymbol{x},t)}
{\Sigma_{\mathrm{SFR}}(\boldsymbol{x},0)}
=\exp[-\mathcal E(\boldsymbol{x},t)].
\tag{17}
$$

$\Gamma_{\mathrm{SFR}}$ has units inverse time and $\mathcal E$ is dimensionless. A temporary enhancement corresponds to a negative contribution to $\mathcal E$. Equation (17) is an exact integral identity under equation (16); it does not assign a clock to the infall stages.

## 4.4 A solvable atomic-supply and molecular-loss realization

A two-reservoir closure makes the delayed molecular response explicit. Assume that the local atomic-associated reservoir transfers material into the molecular reservoir at the effective rate $\Sigma_{\mathrm{HI}}/\tau_{\mathrm{conv}}$, with no additional accretion, reverse phase conversion, or radial transport during the solved interval. Prompt return is assigned to the molecular reservoir as stated above. Then

$$
\begin{aligned}
\Sigma_\Phi&=\frac{\Sigma_{\mathrm{HI}}}{\tau_{\mathrm{conv}}},\\
\frac{d\Sigma_{\mathrm{HI}}}{dt}
&=-\frac{\Sigma_{\mathrm{HI}}}{\tau_{\mathrm{conv}}}
-k_{\mathrm{HI}}\Sigma_{\mathrm{HI}},\\
\frac{d\Sigma_{\mathrm{H_2}}}{dt}
&=\frac{\Sigma_{\mathrm{HI}}}{\tau_{\mathrm{conv}}}
-\left(k_{\mathrm{H_2}}+\frac{1-R+\eta}{\tau_{\mathrm{dep}}}\right)
\Sigma_{\mathrm{H_2}}.
\end{aligned}
\tag{18}
$$

Equation (18) is a **new, deliberately restricted closure in this report**, motivated by reservoir regulation and environmental gas loss. It is not a chemical-rate law claimed by the cited papers. In particular, $\tau_{\mathrm{conv}}$ measures transfer from the selected atomic reservoir, whereas $\tau_\Phi$ measures replenishment of the molecular reservoir:

$$
\tau_\Phi=
\tau_{\mathrm{conv}}\frac{\Sigma_{\mathrm{H_2}}}{\Sigma_{\mathrm{HI}}},
\qquad
\kappa_{\mathrm{HI}}=k_{\mathrm{HI}}+\tau_{\mathrm{conv}}^{-1},
\qquad
\lambda_{\mathrm{H_2}}=
k_{\mathrm{H_2}}+\frac{1-R+\eta}{\tau_{\mathrm{dep}}}.
\tag{19}
$$

The two coefficients $\kappa_{\mathrm{HI}}$ and $\lambda_{\mathrm{H_2}}$ are defined inverse-time sums, introduced solely to keep the exact solution legible. They are not separate physical processes or extra independent parameters.

**Step 1: solve the atomic reservoir.** Dividing its equation by $\Sigma_{\mathrm{HI}}$ and integrating gives

$$
\ln\!\left[\frac{\Sigma_{\mathrm{HI}}(t)}{\Sigma_{\mathrm{HI},0}}\right]
=-\kappa_{\mathrm{HI}}t,
\qquad
\Sigma_{\mathrm{HI}}(t)=
\Sigma_{\mathrm{HI},0}e^{-\kappa_{\mathrm{HI}}t}.
\tag{20}
$$

The subscript zero denotes the initial value in the aperture.

**Step 2: use an integrating factor for H2.** Insert equation (20) into the molecular equation and multiply by $e^{\lambda_{\mathrm{H_2}}t}$:

$$
\frac{d}{dt}
\left[e^{\lambda_{\mathrm{H_2}}t}\Sigma_{\mathrm{H_2}}(t)\right]
=
\frac{\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
e^{(\lambda_{\mathrm{H_2}}-\kappa_{\mathrm{HI}})t}.
\tag{21}
$$

**Step 3: integrate from zero to $t$.** For unequal coefficients, evaluation of the exponential integral and division by the integrating factor gives

$$
\Sigma_{\mathrm{H_2}}(t)=
\Sigma_{\mathrm{H_2},0}e^{-\lambda_{\mathrm{H_2}}t}
+\frac{\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
\frac{e^{-\kappa_{\mathrm{HI}}t}-e^{-\lambda_{\mathrm{H_2}}t}}
{\lambda_{\mathrm{H_2}}-\kappa_{\mathrm{HI}}}.
\tag{22}
$$

The second term is the surviving contribution from material supplied after $t=0$. Its numerator and denominator have the same sign, so this supply contribution is nonnegative.

**Step 4: handle the coincident-timescale limit.** Taking the limit $\lambda_{\mathrm{H_2}}\rightarrow\kappa_{\mathrm{HI}}$, or directly integrating equation (21) with a constant right-hand factor, yields

$$
\Sigma_{\mathrm{H_2}}(t)=
e^{-\lambda_{\mathrm{H_2}}t}
\left(\Sigma_{\mathrm{H_2},0}
+\frac{\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}t\right).
\tag{23}
$$

The apparent singularity in equation (22) is therefore removable. Numerical evaluation near equality should use equation (23) or a stable exponential-difference expression.

**Step 5: verify conservation.** Add the two reservoir equations:

$$
\frac{d}{dt}\left(\Sigma_{\mathrm{HI}}+\Sigma_{\mathrm{H_2}}\right)
=
-(1-R+\eta)\Sigma_{\mathrm{SFR}}
-k_{\mathrm{HI}}\Sigma_{\mathrm{HI}}
-k_{\mathrm{H_2}}\Sigma_{\mathrm{H_2}}.
\tag{24}
$$

Internal phase transfer cancels. The total neutral-associated reservoir changes only through net stellar locking/feedback and environmental removal in this closure. This check would fail if the same stripping or phase-transfer term were counted twice.

Rapid atomic loss makes $\Sigma_\Phi$ fall promptly, while the molecular reservoir retains memory of its initial content and recent supply. A larger outer $k_{\mathrm{H_2}}$ makes the outer SF response faster. The model permits finite inner H2 after substantial outer suppression; it does not assume that H2 is immune to stripping.

## 4.5 A more general supply history and identifiable combinations

The atomic closure is not required to obtain a useful response equation. With constant $\lambda_{\mathrm{H_2}}$ but arbitrary supply, equation (15) gives

$$
\Sigma_{\mathrm{H_2}}(t)=
\Sigma_{\mathrm{H_2},0}e^{-\lambda_{\mathrm{H_2}}t}
+\int_0^t
\Sigma_\Phi(t')e^{-\lambda_{\mathrm{H_2}}(t-t')}\,dt'.
\tag{25}
$$

This follows from the same integrating-factor step as equation (21). It is a causal convolution: gas supplied at time $t'$ is reduced by subsequent depletion/removal before being observed at $t$. Equation (25) retains the gas-supply language of the local regulator without identifying supply exclusively with chemistry.

For illustration, suppose an initially equilibrated reservoir has $\Sigma_{\Phi,0}=\lambda_{\mathrm{H_2}}\Sigma_{\mathrm{H_2},0}$ and its supply is reduced to a constant fraction $f_\Phi$ of that value, with $0\leq f_\Phi\leq1$. Direct integration of equation (25) gives

$$
\frac{\Sigma_{\mathrm{SFR}}(t)}{\Sigma_{\mathrm{SFR}}(0)}
=f_\Phi+(1-f_\Phi)e^{-t/\tau_{\mathrm{q}}},
\qquad
\tau_{\mathrm{q}}=\lambda_{\mathrm{H_2}}^{-1}.
\tag{26}
$$

Here $\tau_{\mathrm{dep}}$ is held fixed. Complete interruption gives pure exponential decline; partial supply leaves a nonzero asymptote. A measured suppression factor alone cannot determine both $f_\Phi$ and $t/\tau_{\mathrm{q}}$, much less the separate values of $k_{\mathrm{H_2}}$, $\eta$, $R$, and $\tau_{\mathrm{dep}}$.

Initial enhancement is also possible. From equation (16), it requires

$$
\frac{\Sigma_\Phi}{\Sigma_{\mathrm{H_2}}}
>
\frac{1-R+\eta}{\tau_{\mathrm{dep}}}
+k_{\mathrm{H_2}}
+\frac{1}{\tau_{\mathrm{dep}}}\frac{d\tau_{\mathrm{dep}}}{dt}.
\tag{27}
$$

A decreasing depletion time contributes negatively to the right-hand side. Therefore compression can temporarily raise SF by increasing molecular supply, increasing molecular content through transport, or reducing $\tau_{\mathrm{dep}}$, even while environmental loss grows elsewhere. This condition is more general than assuming that ram pressure always enhances or always suppresses SF.

## 4.6 Spatially limited compression and NGC4654

At a fixed location, define the molecular fraction of neutral gas as $f_{\mathrm{H_2}}=\Sigma_{\mathrm{H_2}}/(\Sigma_{\mathrm{HI}}+\Sigma_{\mathrm{H_2}})$ and the neutral surface density as $\Sigma_{\mathrm{neutral}}=\Sigma_{\mathrm{HI}}+\Sigma_{\mathrm{H_2}}$. Comparing a compressed state with its specified baseline,

$$
\Delta\ln\Sigma_{\mathrm{SFR}}
=
\Delta\ln\Sigma_{\mathrm{neutral}}
+\Delta\ln f_{\mathrm{H_2}}
-\Delta\ln\tau_{\mathrm{dep}}.
\tag{28}
$$

Equation (28) is an identity, obtained by writing $\Sigma_{\mathrm{SFR}}=\Sigma_{\mathrm{neutral}}f_{\mathrm{H_2}}/\tau_{\mathrm{dep}}$. It distinguishes accumulation, phase partition, and efficiency. A 30% increase in molecular surface density at fixed depletion time produces 0.114 dex of enhanced SFR; unchanged molecular content with a 25% shorter depletion time produces 0.125 dex. Optical residuals at that level do not distinguish these possibilities.

The empirical pressure prescription of [Blitz & Rosolowsky (2006), equation (11)](#ref-br06) can be written

$$
\mathcal R_{\mathrm{mol}}=
\frac{\Sigma_{\mathrm{H_2}}}{\Sigma_{\mathrm{HI}}}
=\left(\frac{P_{\mathrm{hyd}}}{P_0}\right)^\alpha,
\qquad
f_{\mathrm{H_2}}=\frac{\mathcal R_{\mathrm{mol}}}{1+\mathcal R_{\mathrm{mol}}}.
\tag{29}
$$

$P_{\mathrm{hyd}}$ is the hydrostatic pressure external to clouds, $P_0$ a fitted normalization, and $\alpha$ a fitted dimensionless exponent. This report does not insert $P_{\mathrm{ram}}$ for $P_{\mathrm{hyd}}$ or apply a fitted kpc-scale relation unmodified at each MAUVE pixel. Such a substitution would require a new stress, vertical-equilibrium, and averaging prescription. Likewise, $\mathcal R_{\mathrm{mol}}$ alone supplies no absolute formation or replenishment clock.

Spatial variation can be represented by $\Sigma_\Phi(r,\varphi,t)$, $k_{\mathrm{H_2}}(r,\varphi,t)$, and $\tau_{\mathrm{dep}}(r,\varphi,t)$ in equations (15)-(27). A northwest sector in NGC4654 can satisfy equation (27) while the opposite hemisphere does not. Its inferred 0.131-dex contrast is of the same descriptive scale as modest changes in equation (28), but the approximately 0.035-dex facing-side offset from the control is the relevant absolute enhancement measure. Neither quantity is fitted here to a pressure or efficiency parameter.

The analytical disc model of [Lizee et al. (2021)](#ref-lizee) is a useful route toward a more physical compression closure because it couples turbulent support, molecular formation, and SF. It also makes nonstationary approximations and requires gas information absent from the present optical fit. This report therefore uses it as a physically specific comparison rather than importing its complete parameter system without the needed observations.

## 4.7 From molecular SF to young-star H-alpha emission

The instantaneous $\Sigma_{\mathrm{SFR}}$ in equation (15) is not exactly the H-alpha-based estimator during rapid change. Young stars formed at earlier times still contribute ionizing photons. Stellar-population calibration therefore introduces a response kernel; the dependence of SFR tracers on recent SF history is discussed by [Kennicutt & Evans (2012)](#ref-ke12). The following exponential kernel is an analytical approximation introduced here, not a population-synthesis prediction.

**Step 1: write the causal luminosity response.** Let $a$ be stellar age and $K_\alpha(a)$ a nonnegative kernel with unit integral. For constant young-photon absorption fraction $f_{\mathrm{abs,HII}}$,

$$
\mathcal L_{\alpha}^{\mathrm{HII}}(t)=
\frac{f_{\mathrm{abs,HII}}}{C_\alpha}
\int_0^\infty K_\alpha(a)\,
\Sigma_{\mathrm{SFR}}(t-a)\,da,
\qquad
\int_0^\infty K_\alpha(a)\,da=1.
\tag{30}
$$

$C_\alpha$ refers to the adopted fully absorbed, dust-corrected calibration; $f_{\mathrm{abs,HII}}$ specifies the fraction assigned to the HII component. The report calculations use $C_\alpha=4.983582\times10^{-42}\ M_\odot\,\mathrm{yr^{-1}}/(\mathrm{erg\,s^{-1}})$, matching the live product coefficient, and $f_{\mathrm{abs,HII}}=1$ in the simplified examples. Their separate continuing component is consequently interpreted as old-star or mechanical emission, not a second allocation of the same young ionizing photons. A leakage model must partition the young-photon budget consistently between HII emission, diffuse absorption, dust absorption, and escape.

**Step 2: choose and solve a simple response kernel.** For

$$
K_\alpha(a)=\frac{1}{\tau_{\mathrm{ion}}}e^{-a/\tau_{\mathrm{ion}}},
\qquad
\tau_{\mathrm{ion}}\frac{d\mathcal L_{\alpha}^{\mathrm{HII}}}{dt}
+\mathcal L_{\alpha}^{\mathrm{HII}}
=\frac{f_{\mathrm{abs,HII}}}{C_\alpha}\Sigma_{\mathrm{SFR}}(t),
\tag{31}
$$

the differential equation follows by differentiating the convolution. $\tau_{\mathrm{ion}}$ is a luminosity-response time for the young stellar population. It is not the hydrogen recombination time.

Assume constant SF before $t=0$ and $\Sigma_{\mathrm{SFR}}(t)=\Sigma_{\mathrm{SFR}}(0)e^{-t/\tau_{\mathrm{q}}}$ afterward. The pre-existing stars supply the initial condition $\mathcal L_\alpha^{\mathrm{HII}}(0)=f_{\mathrm{abs,HII}}\Sigma_{\mathrm{SFR}}(0)/C_\alpha$. Multiplying equation (31) by $e^{t/\tau_{\mathrm{ion}}}$ and integrating the forcing gives

$$
\frac{\mathcal L_{\alpha}^{\mathrm{HII}}(t)}
{\mathcal L_{\alpha}^{\mathrm{HII}}(0)}
=
\frac{\tau_{\mathrm{q}}e^{-t/\tau_{\mathrm{q}}}
-\tau_{\mathrm{ion}}e^{-t/\tau_{\mathrm{ion}}}}
{\tau_{\mathrm{q}}-\tau_{\mathrm{ion}}}.
\tag{32}
$$

For equal timescales, the limit is $(1+t/\tau_{\mathrm{ion}})e^{-t/\tau_{\mathrm{ion}}}$. If gas is removed almost instantaneously, emission still fades over the young-star response. If $\tau_{\mathrm{q}}\gg\tau_{\mathrm{ion}}$, H-alpha closely follows the declining SFR after the short transient. The two-exponential molecular solution in equation (22) is propagated by linear superposition of equation (32), with the initial components summing to the initial reservoir.

At cloud scales, stochastic population sampling and displacement between gas and young stars complicate a smooth kernel. The spatial-scale dependence of gas/SF correlations is explicitly discussed by [Kruijssen & Longmore (2014)](#ref-kl14). The analytic response should therefore be interpreted as a resolved statistical model, not the detailed light curve of every cloud.

## 4.8 Continuing ionization requires both photons or heating and gas

Write the dust-corrected line luminosity as a sum of physically labelled sources:

$$
\mathcal L_\ell=
\mathcal L_\ell^{\mathrm{HII}}
+\mathcal L_\ell^{\mathrm{DIG,young}}
+\mathcal L_\ell^{\mathrm{old}}
+\mathcal L_\ell^{\mathrm{shock}}
+\mathcal L_\ell^{\mathrm{AGN}}.
\tag{33}
$$

The labels denote classical young-star HII emission, diffuse gas powered by young-star photons, evolved-star-powered emission, mechanically excited emission, and AGN-powered emission. They identify ionization/excitation contributions, not mutually exclusive gas phases inferred from the class labels. The relative importance of young leakage and evolved stars is motivated by [Belfiore et al. (2022)](#ref-belfiore); shock-related excitation in stripped gas is motivated by [Tomicic et al. (2021b)](#ref-tomicic35). No component is identified solely by being labelled NSF.

For an evolved-star component, define $q_{\mathrm{H,old}}$ as the ionizing photon production rate per stellar mass, $f_{\mathrm{abs,old}}$ as the fraction absorbed by local gas, $h_{\mathrm P}$ as Planck's constant, $\nu_\alpha$ as H-alpha frequency, and $p_\alpha$ as the number of H-alpha photons per absorbed ionizing photon under the assumed recombination conditions. Then

$$
\mathcal L_\alpha^{\mathrm{old}}
=
h_{\mathrm P}\nu_\alpha p_\alpha
f_{\mathrm{abs,old}}q_{\mathrm{H,old}}\Sigma_*.
\tag{34}
$$

This follows by multiplying absorbed ionizations by H-alpha energy emitted per recombination. Case-B hydrogenic recombination calculations are supplied by [Hummer & Storey (1987)](#ref-hs87); the source fractions and stellar spectrum require their own population model. At fixed stellar population, the slowly changing quantity is the photon supply per stellar mass. The emitted luminosity can still decline because gas covering fraction, absorption, density, and geometry change.

There is also a gas-capacity constraint. In recombination-dominated gas at specified temperature,

$$
\mathcal L_\alpha=
h_{\mathrm P}\nu_\alpha\alpha_\alpha^{\mathrm{eff}}
\int n_e n_p\,dz,
\qquad
\tau_{\mathrm{rec}}=\frac{1}{\alpha_B n_e}.
\tag{35}
$$

Here $\alpha_\alpha^{\mathrm{eff}}$ is the effective H-alpha recombination coefficient, $\alpha_B$ is the case-B total coefficient, and $n_e,n_p$ are electron and proton number densities; the integral is through the emitting column with consistent area units. The second expression follows from recombinations per proton at approximately fixed $n_e$. Sustained emission over a timescale much longer than $\tau_{\mathrm{rec}}$ requires continuing ionization or heating. It is not merely a long-lived afterglow of gas that was ionized once.

Equations (34)-(35) constrain any proposed background normalization. A large stellar surface density can support more evolved-star ionization, but gas removal can reduce the absorbing/emitting column. A shock interpretation must similarly supply mechanical power; an AGN interpretation requires evidence for that source. The illustrative continuing luminosities below are prescribed conditional normalizations, not measurements of $q_{\mathrm{H,old}}$, emission measure, or shock power in MAUVE. Their energy and gas budgets remain to be tested with additional diagnostics.

For compact analytic expressions, define $\mathcal L_\ell^{\mathrm{cont}}$ as the sum of the continuing components in equation (33) other than the local HII term. This is a physically labelled sum, not a newly identified source. A two-component reduction uses $\mathcal L_\ell=\mathcal L_\ell^{\mathrm{HII}}+\mathcal L_\ell^{\mathrm{cont}}$. It is useful only when the effective continuing-component spectrum is specified and its time variation is acknowledged.

## 4.9 Deriving line-ratio changes from component evolution

**Step 1: add fluxes before taking ratios.** Define
$w_\alpha=\mathcal L_\alpha^{\mathrm{cont}}/
(\mathcal L_\alpha^{\mathrm{HII}}+\mathcal L_\alpha^{\mathrm{cont}})$.
For a line with H-alpha denominator,

$$
\mathcal R_{\ell/\alpha}
=(1-w_\alpha)\mathcal R_{\ell/\alpha}^{\mathrm{HII}}
+w_\alpha\mathcal R_{\ell/\alpha}^{\mathrm{cont}}.
\tag{36}
$$

This is an algebraic consequence of luminosity addition, not an interpolation in logarithmic BPT coordinates. Flux-component modelling is motivated by the DIG literature cited above, but equation (36) itself is derived directly from the definition of a ratio.

**Step 2: use the correct Balmer weight.** Let $B_{\mathrm{HII}}$ and $B_{\mathrm{cont}}$ denote the H-alpha/H-beta luminosity ratios of the two components. The continuing fraction for H-beta is

$$
w_\beta=
\frac{w_\alpha/B_{\mathrm{cont}}}
{(1-w_\alpha)/B_{\mathrm{HII}}+w_\alpha/B_{\mathrm{cont}}},
\qquad
\mathcal R_{\ell/\beta}
=(1-w_\beta)\mathcal R_{\ell/\beta}^{\mathrm{HII}}
+w_\beta\mathcal R_{\ell/\beta}^{\mathrm{cont}}.
\tag{37}
$$

Only equal component decrements make $w_\beta=w_\alpha$. Dust and excitation can invalidate that equality, and correcting a mixture with one observed decrement need not recover each intrinsic component separately.

**Step 3: differentiate the mixture weight.** For locally exponential H-alpha components with positive fading times $\tau_{\alpha,\mathrm{HII}}$ and $\tau_{\alpha,\mathrm{cont}}$,

$$
\frac{d\mathcal L_\alpha^j}{dt}
=-\frac{\mathcal L_\alpha^j}{\tau_{\alpha,j}},
\qquad
\frac{dw_\alpha}{dt}
=w_\alpha(1-w_\alpha)
\left(\frac{1}{\tau_{\alpha,\mathrm{HII}}}
-\frac{1}{\tau_{\alpha,\mathrm{cont}}}\right),
\tag{38}
$$

where $j$ labels HII or continuing emission. To obtain the second expression, differentiate the ratio defining $w_\alpha$, substitute the first expression for both components, and collect the product of their fractional contributions. If HII emission fades faster, the continuing fraction grows even though both luminosities decline.

**Step 4: include intrinsic spectral evolution.** Differentiating equation (36) gives

$$
\begin{aligned}
\frac{d\mathcal R_{\ell/\alpha}}{dt}
={}&
\left(\mathcal R_{\ell/\alpha}^{\mathrm{cont}}
-\mathcal R_{\ell/\alpha}^{\mathrm{HII}}\right)
\frac{dw_\alpha}{dt}\\
&+(1-w_\alpha)\frac{d\mathcal R_{\ell/\alpha}^{\mathrm{HII}}}{dt}
+w_\alpha\frac{d\mathcal R_{\ell/\alpha}^{\mathrm{cont}}}{dt}.
\end{aligned}
\tag{39}
$$

For fixed templates, a growing continuing fraction raises a ratio only when that component has the larger intrinsic ratio. [N II]/H-alpha and [S II]/H-alpha may therefore rise while [O III]/H-beta is constant or decreases. Evolving metallicity, ionization parameter, or stellar spectrum enters through the last two terms. The quenching models of [Citro et al. (2017)](#ref-citro) emphasize why those terms cannot always be neglected.

This derivation gives a precise meaning to *differential fading*. The relevant inequality is between the fractional luminosity evolution of physically specified components, together with their line spectra. It does not assert that every forbidden transition has a universally longer fading time than every Balmer transition.

## 4.10 Equivalent width, velocity width, and the observation operator

The model must next generate the quantities used by the actual selection. For positive emission EW,

$$
\mathrm{EW}(\mathrm{H}\alpha)=
\frac{\mathcal L_\alpha^{\mathrm{HII}}+
\mathcal L_\alpha^{\mathrm{cont}}}
{\mathcal C_{\lambda,\alpha}}.
\tag{40}
$$

The continuum can contain old and young stars and can itself evolve. Holding it fixed in a short illustrative interval emphasizes how falling line emission crosses the 6-Angstrom cut. This is a selection transition, not proof that all SF has ceased.

For two line components with intrinsic dispersions $\sigma_{\mathrm{HII}}$, $\sigma_{\mathrm{cont}}$ and centroid difference $\Delta v$, their exact second central moment is

$$
\sigma_\alpha^2=
(1-w_\alpha)\sigma_{\mathrm{HII}}^2+
w_\alpha\sigma_{\mathrm{cont}}^2+
w_\alpha(1-w_\alpha)(\Delta v)^2.
\tag{41}
$$

Equation (41) follows by adding each component's variance and squared offset from the combined centroid. A single-Gaussian fitted width need not equal this moment for strongly non-Gaussian profiles; instrumental convolution and the pipeline fit must be included in a realistic prediction. A broader continuing component can cause the width cut to fail as HII emission fades, but the same cut makes a resulting SF/NSF width difference partly selection-induced.

Define a raw detection indicator $D_\ell$ for each line using its measured flux and error. In the live convention,

$$
D_\ell=
\boldsymbol{1}\!\left[
\frac{F_\ell}{\sigma_{F_\ell}}\geq3
\ \text{and}\ 
F_\ell\geq F_{\mathrm{floor}}\right],
\qquad
F_{\mathrm{floor}}=2.0\times10^{-19}\ 
\mathrm{erg\,s^{-1}\,cm^{-2}}.
\tag{42}
$$

$\boldsymbol{1}$ is an indicator equal to one when its condition is satisfied. $F_\ell$ is observed, attenuated flux per native observing element, not the dust-corrected luminosity density. With physical emitting area $A_{\mathrm{pix}}$, distance $d$, and attenuation $A_\ell^{\mathrm{att}}$ in magnitudes, the idealized flux prediction is

$$
F_{\ell,\mathrm{pred}}=
\frac{A_{\mathrm{pix}}\mathcal L_\ell}
{4\pi d^2}\,10^{-0.4A_\ell^{\mathrm{att}}}.
\tag{43}
$$

The area and luminosity-density conventions must match, and the map must be PSF-convolved and sampled consistently. An actual comparison then adds measurement noise and passes the resulting fluxes through the same attenuation, BPT, EW, and width procedures as the observations.

On the retained valid domain, define $D_B=D_\alpha D_\beta$ and let $C_{\mathrm{SF}}$ indicate satisfaction of the actual HII-map, EW, and intrinsic-width criteria. The exhaustive partition is

$$
P_{\mathrm{ND}}=1-D_B,\qquad
P_{\mathrm{SF}}=D_B C_{\mathrm{SF}},\qquad
P_{\mathrm{NSF}}=D_B(1-C_{\mathrm{SF}}).
\tag{44}
$$

These are indicators for a deterministic noiseless realization and probabilities after averaging over measurement noise. A physical model predicts class counts by summing equation (44) over the same retained observing elements and then applying equations (1)-(4). Applying class fractions directly to a gas mass, or averaging all model pixels across galaxies without their separate weights, changes the observable.

## 4.11 A closed-form occupancy demonstration

A simplified analytical selection layer reveals how different radial class outcomes can arise. It is a **demonstration**, not a calibrated replacement for the pipeline BPT operator.

At a fixed spatial bin and time, suppose the young HII H-alpha luminosity density varies across equal-area elements according to

$$
u_\alpha\equiv
\ln\!\left(\frac{\mathcal L_\alpha^{\mathrm{HII}}}
{\mathcal L_{\mathrm{ref}}}\right)
\sim\mathcal N(\mu_\alpha,s_\alpha^2).
\tag{45}
$$

$\mathcal N$ is a normal distribution. The median young luminosity is $\mathcal L_{\mathrm{ref}}e^{\mu_\alpha}$ and its mean is $\mathcal L_{\mathrm{ref}}e^{\mu_\alpha+s_\alpha^2/2}$. A distribution of initial gas amplitudes with common response coefficients generates the same multiplicative luminosity evolution, so $\mu_\alpha(t)=\mu_\alpha(0)+\ln[\mathcal L_{\alpha,\mathrm{med}}^{\mathrm{HII}}(t)/\mathcal L_{\alpha,\mathrm{med}}^{\mathrm{HII}}(0)]$ under this closure. Scatter $s_\alpha$ is held fixed; neither lognormality nor that invariance is asserted as measured.

Hold the continuing luminosity and continuum fixed across the elements at a given time. For this demonstration only, adopt equal Balmer decrements and collapse the two raw line limits into one effective H-alpha luminosity limit $\mathcal L_{\alpha,\mathrm{det}}$. The minimum young luminosity for joint detection is

$$
\mathcal L_{\mathrm{young,det}}
=\max\!\left(0,\mathcal L_{\alpha,\mathrm{det}}
-\mathcal L_\alpha^{\mathrm{cont}}\right).
\tag{46}
$$

If the continuing component alone exceeds the effective limit, declining young emission need not produce ND. If it lies below the limit, the same young-star fading can produce ND. This is the key distinction between the illustrative inner and outer cases.

For analytic purposes, replace the line-ratio and width selection with an allowed continuing fraction $w_\alpha<w_{\mathrm{lim}}$, where $0<w_{\mathrm{lim}}<1$. This value must eventually be determined by actual line templates and the pipeline boundaries. Combining detection, EW, and this fraction condition gives the young-luminosity SF threshold

$$
\begin{aligned}
\mathcal L_{\mathrm{young,SF}}=\max\bigg\{&
\mathcal L_{\mathrm{young,det}},\
6\,\mathrm{Angstrom}\,\mathcal C_{\lambda,\alpha}
-\mathcal L_\alpha^{\mathrm{cont}},\\
&\frac{1-w_{\mathrm{lim}}}{w_{\mathrm{lim}}}
\mathcal L_\alpha^{\mathrm{cont}},\ 0\bigg\}.
\end{aligned}
\tag{47}
$$

The EW term has luminosity-density units because the continuum is per Angstrom. For coincident centroids and $\sigma_{\mathrm{cont}}>\sigma_{\mathrm{HII}}$, equation (41) supplies a width-based limit
$w_{\sigma}=(45^2-\sigma_{\mathrm{HII}}^2)/
(\sigma_{\mathrm{cont}}^2-\sigma_{\mathrm{HII}}^2)$
when it lies between zero and one. The demonstration uses the more restrictive of this value and an explicitly assumed BPT-proxy limit.

Let $\mathcal F(\ell)$ be the cumulative probability of young luminosity below $\ell$: it is zero for $\ell\leq0$ and otherwise
$\mathcal F(\ell)=\Phi_{\mathrm N}\{[\ln(\ell/\mathcal L_{\mathrm{ref}})-\mu_\alpha]/s_\alpha\}$,
where $\Phi_{\mathrm N}$ is the standard normal cumulative distribution. Direct integration of the allowed intervals gives

$$
\begin{aligned}
F_{\mathrm{ND}}&=\mathcal F(\mathcal L_{\mathrm{young,det}}),\\
F_{\mathrm{SF}}&=1-\mathcal F(\mathcal L_{\mathrm{young,SF}}),\\
F_{\mathrm{NSF}}&=
\mathcal F(\mathcal L_{\mathrm{young,SF}})
-\mathcal F(\mathcal L_{\mathrm{young,det}}).
\end{aligned}
\tag{48}
$$

The fractions add to unity because the thresholds partition the young-luminosity distribution. Inner gas can pass from SF to NSF when the upper threshold is crossed while detection remains satisfied. Outer gas can instead pass below the detection threshold, producing ND growth.

Intensity predictions require luminosity-weighted integrals, not only counts. Completing the square in
$u_\alpha-(u_\alpha-\mu_\alpha)^2/(2s_\alpha^2)$
shifts the Gaussian mean to $\mu_\alpha+s_\alpha^2$ and gives the luminosity moment above threshold:

$$
\begin{aligned}
\mathcal M_\alpha(>\ell)
&=\int_{\ln(\ell/\mathcal L_{\mathrm{ref}})}^\infty
\mathcal L_{\mathrm{ref}}e^u
\frac{e^{-(u-\mu_\alpha)^2/(2s_\alpha^2)}}
{s_\alpha\sqrt{2\pi}}\,du\\
&=\mathcal L_{\mathrm{ref}}e^{\mu_\alpha+s_\alpha^2/2}
\Phi_{\mathrm N}\!\left(
\frac{\mu_\alpha+s_\alpha^2-\ln(\ell/\mathcal L_{\mathrm{ref}})}
{s_\alpha}\right).
\end{aligned}
\tag{49}
$$

For $\ell\leq0$, $\mathcal M_\alpha(>\ell)$ is the full mean. Let $\mathcal M_{\alpha,c}$ denote the young-luminosity moment over a class interval, obtained by subtracting two such moments when needed. Then, for one model population,

$$
I_{\mathrm{SF}}^{\mathrm{app}}
=C_\alpha
\frac{\mathcal M_{\alpha,\mathrm{SF}}
+\mathcal L_\alpha^{\mathrm{cont}}F_{\mathrm{SF}}}
{F_{\mathrm{SF}}},
\qquad
J_{\alpha,\mathrm{NSF}}=
\mathcal M_{\alpha,\mathrm{NSF}}
+\mathcal L_\alpha^{\mathrm{cont}}F_{\mathrm{NSF}}.
\tag{50}
$$

The superscript “app” denotes the apparent H-alpha-based SFR assigned to the selected SF population. It includes any continuing emission still admitted by that selection, just as a mixed spectrum can contribute H-alpha to an observational SFR estimator. It is not the same as the true instantaneous SF in equation (15). For multiple galaxies, the model must construct each galaxy's numerator and occupancy before reproducing the stage estimator.

Equation (49) explains a survivor bias: as the faint end leaves the SF class, the conditional mean remains weighted toward brighter elements. Therefore a flat or even elevated outer surviving-SF intensity does not rule out strong removal of most outer SF area. Conversely, matching a decline in the conditional mean can require a larger change in the underlying luminosity distribution than the observed mean contrast alone suggests.

# 5. Numerical predictions and their relation to the observations

The calculations in this section evaluate the analytical equations with declared parameters. They are illustrative forward predictions, separate from the observational fit in Section 6. Their purpose is to test the internal mathematics and show which observed behaviours follow, which require additional conditions, and which are not quantitatively reproduced.

## 5.1 Gas supply loss, delayed SF fading, and different inner/outer classes

The two example locations begin with a molecular reservoir approximately balanced by supply before environmental removal is switched on. Maintaining a constant pre-perturbation atomic reservoir would require external replenishment before $t=0$; the no-additional-supply assumption in equation (18) applies afterward. The initial gas values below are medians of the multiplicative area-element population, with common HI/H2 ratio and response coefficients. They are not measured MAUVE gas columns.

| Parameter | Inner illustration | Outer illustration |
|---|---:|---:|
| Initial $\Sigma_{\mathrm{HI}}$ | $2.4\ M_\odot\,\mathrm{pc}^{-2}$ | $8.0\ M_\odot\,\mathrm{pc}^{-2}$ |
| Initial $\Sigma_{\mathrm{H_2}}$ | $8.0\ M_\odot\,\mathrm{pc}^{-2}$ | $4.0\ M_\odot\,\mathrm{pc}^{-2}$ |
| $\tau_{\mathrm{dep}}$ | 2000 Myr | 2000 Myr |
| $\tau_{\mathrm{conv}}$ | 1000 Myr | 6666.7 Myr |
| $k_{\mathrm{HI}}^{-1}$ | 160 Myr | 80 Myr |
| $k_{\mathrm{H_2}}^{-1}$ | 650 Myr | 180 Myr |
| Initial $\mathcal L_\alpha^{\mathrm{cont}}/\mathcal L_{\mathrm{ref}}$ | 2.0 | 0.05 |
| Continuing-component exponential fading time | 1200 Myr | 100 Myr |
| $\mathcal C_{\lambda,\alpha}/\mathcal L_{\mathrm{ref}}$ | $1.2\ \mathrm{Angstrom}^{-1}$ | $0.05\ \mathrm{Angstrom}^{-1}$ |

Both use $R=0.4$, $\eta=0$, $\tau_{\mathrm{ion}}=5$ Myr, $\mathcal L_{\mathrm{ref}}=10^{38}\ \mathrm{erg\,s^{-1}\,kpc^{-2}}$, $s_\alpha=1.1$ in natural-log units, and $\mathcal L_{\alpha,\mathrm{det}}=0.8\mathcal L_{\mathrm{ref}}$. The long effective $\tau_{\mathrm{conv}}$ values are reservoir-transfer times, not chemical H2 formation times. They were selected to provide the stated pre-perturbation molecular balance, not inferred from data.

For the selection demonstration, the assumed BPT-proxy limit is 0.45. Intrinsic widths of 20 and $65\ \mathrm{km\,s^{-1}}$, equal centroids, and the observed 45-$\mathrm{km\,s^{-1}}$ threshold give $w_\sigma=0.42484$; hence $w_{\mathrm{lim}}=0.42484$. The numerical BPT proxy is not derived from a full line-ratio grid. A fixed effective detection limit also omits the real pixel-dependent noise and attenuation. These choices make the example analytically transparent while limiting its quantitative interpretation.

![Figure 5. Exact two-reservoir solutions propagated through the exponential young-star response and the analytical class operator. The top row represents the inner example and the bottom row the outer example. Time is illustrative and is not assigned to catalogue stages. Selected-population intensity need not track the decline of the underlying median molecular reservoir. The outer example has no NSF interval under its adopted thresholds, so no outer NSF-intensity curve is defined.](assets/20260914_resolved_RPS_academic_model/figure_05_gas_to_classes.png)

After 600 Myr in this example:

| Predicted quantity | Inner case | Outer case | Connection to the measured pattern |
|---|---|---|---|
| Molecular reservoir relative to initial value | 0.350 | 0.0311 | Retained inner gas coexists with a much larger outer depletion. |
| SF occupancy | $0.653\rightarrow0.248$ | $0.936\rightarrow0.048$ | SF area declines in both cases, especially outside. |
| NSF occupancy | $0.347\rightarrow0.752$ | Remains zero in this simplified outer case | Continuing inner emission keeps regions detected after they fail SF cuts. |
| ND occupancy | Remains zero at this time | $0.064\rightarrow0.952$ | Weak outer continuing emission allows faded regions to fail detection. |
| Apparent selected-SF intensity contrast | $-0.183$ dex | $-0.749$ dex | Survivor selection weakens the mapping from gas depletion to conditional intensity. |
| Conditional NSF H-alpha contrast | $-0.138$ dex | Undefined | Increasing inner NSF area can coexist with declining NSF intensity. |

The example reproduces the *directions* of outer ND growth, inner NSF growth, and lower selected-SF intensity. Its class fractions, times, and amplitudes are not a fit. It does not reproduce the observed approximately $-0.4$ dex inner surviving-SF suppression with these parameters; it predicts only $-0.183$ dex. It also predicts outer NSF occupancy identically zero and no inner ND at 600 Myr, unlike a heterogeneous galaxy population. Those are consequences of the simplified fixed thresholds and backgrounds.

An especially useful failed prediction is the NSF retained-area contribution: it rises from $1.69\times10^{38}$ to $2.66\times10^{38}\ \mathrm{erg\,s^{-1}\,kpc^{-2}}$ in the inner example, despite the lower conditional NSF intensity. Newly admitted NSF elements more than compensate their fading. This differs from the declining contribution in several observed inner density bins. A verbal fading scenario is therefore insufficient to predict the luminosity budget.

For a single common population, writing $J_{\alpha,\mathrm{NSF}}=F_{\mathrm{NSF}}I_{\alpha,\mathrm{NSF}}$ gives

$$
\frac{d\ln J_{\alpha,\mathrm{NSF}}}{dt}<0
\quad\Longleftrightarrow\quad
\frac{d\ln I_{\alpha,\mathrm{NSF}}}{dt}
<
-\frac{d\ln F_{\mathrm{NSF}}}{dt}.
\tag{51}
$$

Thus NSF intensity must fade sufficiently rapidly to outweigh occupancy growth. For observed stage averages, equation (51) must be applied to galaxy-level products before averaging, because equation (4) contains covariance and support differences. A future joint fit must use $F_{\mathrm{NSF}}$, $I_{\alpha,\mathrm{NSF}}$, and $J_{\alpha,\mathrm{NSF}}$ simultaneously. Matching only the class map can select the wrong luminosity evolution.

## 5.2 Increasing low-ionization ratios while all emission fades

A separate line-mixture experiment isolates equation (39). Set the initial young and continuing H-alpha luminosities to $8\mathcal L_{\mathrm{ref}}$ and $2\mathcal L_{\mathrm{ref}}$, with exponential fading times of 200 and 700 Myr, respectively. Adopt equal component Balmer decrements and fixed component line ratios:

| Component | [N II]/H-alpha | [S II]/H-alpha | [O III]/H-beta |
|---|---:|---:|---:|
| HII | 0.25 | 0.18 | 0.60 |
| Continuing | 1.00 | 0.65 | 0.30 |

These are illustrative templates chosen to expose the sign conditions, not calibrated shock or stellar-photoionization models. A real spectrum must be checked against photoionization/shock grids and measured line fluxes.

![Figure 6. A line-mixture calculation in which both components fade. The increasing continuing-component fraction raises [N II]/H-alpha and [S II]/H-alpha, while [O III]/H-beta decreases for the stated templates. The BPT trajectory is rightward and downward. No observed BPT class boundary is imposed in this experiment.](assets/20260914_resolved_RPS_academic_model/figure_06_line_mixing.png)

At 600 Myr, total H-alpha decreases from 10 to approximately $1.247\mathcal L_{\mathrm{ref}}$. The continuing fraction rises from 0.20 to approximately 0.681. The three ratios change from 0.400 to 0.760, from 0.274 to 0.500, and from 0.540 to 0.396, respectively. All three forbidden-line luminosities also decline. Larger low-ionization ratios are therefore compatible with fainter gas emission, and declining [O III]/H-beta is compatible with an increasing contribution from this continuing template.

This experiment addresses the sign pattern, not the observed numerical line ratios or a full excitation mechanism. It does not imply that an old-star spectrum necessarily has the chosen [O III]/H-beta ratio. If the continuing template instead has higher [O III]/H-beta, equation (39) predicts an upward change. This flexibility is physically useful only if constrained by independent source diagnostics.

## 5.3 Predictions that distinguish competing explanations

Gas depletion at nearly fixed $\tau_{\mathrm{dep}}$ predicts lower matched-resolution $\Sigma_{\mathrm{H_2}}$ in suppressed SF regions. A depletion-time change predicts lower $\Sigma_{\mathrm{SFR}}/\Sigma_{\mathrm{H_2}}$ after matching resolution, class support, and SFR timescale. A falling young-star contribution relative to evolved-star ionization predicts coupled changes in continuum-normalized H-alpha, line ratios, and independently inferred recent stellar populations. Shock dominance requires a compatible mechanical-energy budget and excitation/kinematic evidence.

For NGC4654, the terms in equation (28) can be tested using molecular/atomic maps and matched-scale SFR residuals in the fixed northwest sector and the opposite side. Higher molecular content at unchanged depletion time differs from similar molecular content at shorter depletion time. Neither should be inferred simply from a positive optical residual. The published gas results make this a promising individual-galaxy test.

# 6. A quantitative fit to the surviving-SF intensity profiles

## 6.1 What is fitted

The new fit targets the stage contrast in $I_{\mathrm{SF}}$ defined by equation (2), using seven stellar-density bins with centres 7.625, 7.875, 8.125, 8.375, 8.625, 8.875, and 9.125. This is the interval $7.5\leq\log_{10}(\Sigma_*/M_\odot\,\mathrm{kpc}^{-2})<9.25$. It avoids the most weakly supported density edges while retaining a broad resolved range. The interval was selected for this exploratory report; it was not preregistered.

For stage $s$ and bin $b$, define

$$
y_{s,b}=
\log_{10}\!\left[
\frac{I_{\mathrm{SF},s,b}}{I_{\mathrm{SF,pre},b}}\right],
\qquad
I_{\mathrm{SF},s,b}^{\mathrm{model}}
=I_{\mathrm{SF,pre},b}\,10^{a_s}.
\tag{52}
$$

$a_s$ is a constant logarithmic attenuation for one stage. The observed pre-peak profile supplies the spatial dependence of the baseline. The fit retains the resolved stellar-density profile while asking whether a single multiplicative change describes its surviving-SF intensity.

The connection to equation (17) is

$$
\mathcal E_s^{\mathrm{app}}=-a_s\ln 10.
\tag{53}
$$

The superscript “app” is essential. This is an effective attenuation of a category-conditioned H-alpha estimator. Equating it to the true local gas/SF attenuation $\mathcal E(\boldsymbol{x},t)$ additionally requires comparable baseline populations, sufficiently slow luminosity response, controlled contamination, and a model for which elements survive selection. The fitted quantity is identifiable as a profile normalization; the individual physical terms in equation (16) are not.

## 6.2 Derivation of the estimator and uncertainty procedure

Let $s_{s,b}$ be the standard deviation of the bin's contrast under whole-galaxy bootstrap resampling, and define fixed descriptive weights $w_{s,b}=s_{s,b}^{-2}$. The fitting objective and its derivative are

$$
Q(a_s)=\sum_b w_{s,b}(y_{s,b}-a_s)^2,
\qquad
\frac{dQ}{da_s}
=-2\sum_b w_{s,b}(y_{s,b}-a_s).
\tag{54}
$$

Setting the derivative to zero and dividing by the positive weight sum yields

$$
\widehat a_s=
\frac{\sum_b w_{s,b}y_{s,b}}{\sum_b w_{s,b}}.
\tag{55}
$$

This is a weighted least-squares location estimate. The objective is not used as an independent-bin chi-squared likelihood. Bin contrasts share galaxies and the same reference, so their correlations are substantial.

Uncertainty is calculated with 10,000 draws using seed 20260914. In each draw, entire galaxy products are resampled within each stage, retaining each selected product's full set of bins and paired occupancy/intensity numerators. Both the pre-peak and target-stage profiles are recomputed before refitting equation (55), using the full-sample weights. A zero-SF denominator can make a draw undefined; complete finite contrast vectors are retained. There are 9,992 usable close-to-peak draws and 9,998 usable post-peak draws. The excluded fractions are 0.08% and 0.02%, respectively.

This propagates sampled galaxy-to-galaxy variation and reference uncertainty without treating spaxels as independent. It does not include full line-flux systematics, model discrepancy, stage-assignment uncertainty, or uncertainty in field equivalence.

Two sensitivity checks are included. Equal-bin weighting tests dependence on the inverse-variance weights. Removing one entire pre-peak or target-stage galaxy at a time tests sensitivity to individual systems, recomputing the profiles with the same full-fit weights. A linear-density extension,

$$
y_{s,b}=a_{s,0}+a_{s,1}(x_b-8.5),
\qquad
\widehat{\boldsymbol a}_s=
(X^{\mathsf T}WX)^{-1}X^{\mathsf T}W\boldsymbol y_s,
\tag{56}
$$

tests whether a slope is required descriptively. Here $x_b$ is the bin centre, $X$ has columns $1$ and $x_b-8.5$, $W$ is the diagonal weight matrix, and $\boldsymbol a_s=(a_{s,0},a_{s,1})^{\mathsf T}$. Equation (56) follows by differentiating the quadratic objective with respect to both coefficients. The same whole-galaxy draws determine their uncertainty. It is a shape-sensitivity model, not an additional physical gas parameter.

## 6.3 Fit results

![Figure 7. New fits to the resolved surviving-SF intensity contrasts. Points show observed stage contrasts and whole-galaxy bootstrap intervals. Solid lines are fitted multiplicative attenuations of the pre-peak spatial profile; bands show their bootstrap uncertainty. Dashed lines show the linear-density sensitivity model. This fit addresses the intensity branch and is separate from the illustrative gas/class calculation in Figure 5.](assets/20260914_resolved_RPS_academic_model/figure_07_attenuation_fit.png)

| Fitted quantity | Close-to-peak | Post-peak |
|---|---|---|
| $\widehat a_s$ (dex), 16th-84th percentile interval | $-0.141\ [-0.248,-0.059]$ | $-0.390\ [-0.468,-0.298]$ |
| Multiplicative intensity relative to control | $0.722\ [0.565,0.873]$ | $0.407\ [0.341,0.504]$ |
| Effective $\mathcal E_s^{\mathrm{app}}$ | $0.326\ [0.136,0.571]$ | $0.898\ [0.686,1.077]$ |
| RMS residual about constant fit | 0.028 dex | 0.094 dex |
| Equal-bin estimate of $a_s$ | $-0.150$ dex | $-0.405$ dex |
| Leave-one-galaxy range of $a_s$ | $[-0.217,-0.086]$ dex | $[-0.420,-0.340]$ dex |
| Linear-density slope, dex per dex | $-0.048\ [-0.222,0.105]$ | $+0.102\ [-0.132,0.278]$ |

The post-peak selected-SF intensity is approximately 41% of the field-like pre-peak baseline over the fitted interval. Its suppression persists when any one pre-peak or post-peak product is omitted. The close-to-peak profile shows a smaller average suppression despite possible localized enhancement in individual regions.

A constant attenuation summarizes the dominant normalization change but does not explain all spatial structure. The largest post-peak residual is 0.214 dex at $x_b=7.875$, where the observed contrast is about $-0.604$ dex. A linear slope reduces the post-peak RMS only from 0.094 to 0.092 dex, and its bootstrap interval includes zero. It does not account for that local departure. For close-to-peak, the linear RMS decreases to 0.011 dex, but the slope interval also includes zero.

Pairwise bootstrap correlations between bins extend to approximately 0.93 for close-to-peak and 0.83 for post-peak. Treating seven bins as seven independent galaxies would give misleading confidence. The bootstrap and leave-one-galaxy results support a descriptive normalization difference; no chi-squared goodness-of-fit probability, formal model-selection preference, or out-of-sample predictive validation is claimed.

## 6.4 What physical inference is and is not justified

The fitted post-peak factor of 0.407 can be produced by different regulator histories. For complete supply interruption with fixed depletion time, equation (26) gives $t/\tau_{\mathrm{q}}=-\ln(0.407)=0.898$. If an independent clock were 300 Myr, the corresponding conditional $\tau_{\mathrm{q}}$ would be approximately 334 Myr; a 600-Myr clock would give approximately 668 Myr. Neither clock is inferred here.

Partial supply produces a different answer:

$$
\frac{t}{\tau_{\mathrm{q}}}
=-\ln\!\left(
\frac{A_{\mathrm{SFR}}-f_\Phi}{1-f_\Phi}\right),
\qquad
0\leq f_\Phi<A_{\mathrm{SFR}}<1.
\tag{57}
$$

Here $A_{\mathrm{SFR}}$ is a specified true-SFR attenuation factor in equation (26); substituting the apparent fitted factor requires the extra assumptions above. With $A_{\mathrm{SFR}}=0.407$ and $f_\Phi=0.2$, equation (57) gives about 1.35 rather than 0.898. A supply fraction at or above the attenuation cannot reach that state in finite time under this constant-coefficient declining model. Time-dependent depletion and source mixtures introduce further alternatives.

The fit constrains an effective intensity attenuation and rules out an entirely unchanged surviving-SF profile over the fitted range under the adopted control. It does **not** establish whether reduced molecular content, reduced efficiency, direct stripping, or declining supply dominates. It does not fit the NSF fractions, ND fractions, BPT ratios, or NGC4654's directional signal. Those components remain explicit forward-model requirements.

## 6.5 How to progress from this partial fit to a physical fit

A physical fit should first reproduce the actual observation operator using observed continuum, noise, geometry, and gas/SF support for each galaxy. Initial SF/luminosity distributions should be anchored to the pre-peak control, with galaxy-specific nuisance variation. Model line fluxes should be generated before dust correction and class assignment, including nondetections instead of discarding them.

The first joint target should combine SF/NSF/ND occupancy, selected-SF intensity, NSF H-alpha intensity, and its area-normalized contribution on common support. Fitting only line ratios or occupancy would leave the failed prediction in Section 5.1 unconstrained. Gas data can constrain $\Sigma_{\mathrm{H_2}}$ and $\tau_{\mathrm{dep}}$; resolved recent SF histories and dynamical modelling can supply temporal information. NGC4654 merits a separate sector-resolved fit with its tidal history represented.

With only 8/5/13 stage products, many free per-bin timescales would be weakly identified. A parsimonious hierarchical fit should vary a small number of physical histories, retain galaxy-level scatter, and compare supply loss alone, molecular loss with fixed depletion time, and models allowing depletion-time changes. Cross-validation should withhold whole galaxies, not random spaxels from the same galaxy.

# 7. Physical interpretation and use of the model

The model organizes the observations into a conditional evolutionary sequence. In weakly bound outer regions, gas loss and declining supply reduce molecular content and young-star emission. Where little continuing ionization remains above the sensitivity limit, those regions become ND. In inner regions, a retained absorbing gas column and a continuing ionizing/heating source can maintain Balmer emission while the young-star contribution falls. Changes in EW, line ratios, and width then place more area in NSF. The selected SF regions can themselves fade; the new fit quantifies that branch independently of occupancy loss.

This interpretation requires neither universal NSF brightening nor a universal upward BPT motion. Line mixing can increase low-ionization ratios while all line luminosities decline, and [O III]/H-beta can respond differently. The observed NSF luminosity budget is a stronger condition than NSF occupancy alone and should be treated as a central constraint in future modelling.

Localized compression is compatible with this sequence when it transiently increases supply or shortens depletion time enough to satisfy equation (27). NGC4654 supplies a specific candidate: its northwest residuals align descriptively with the published compressed gas region, although the absolute median enhancement is modest and tidal effects remain relevant. Population-level close-to-peak attenuation and a localized positive residual are compatible statements at different spatial and statistical levels.

For current scientific writing, the strongest supported interpretation is that the optical data favour spatial loss of detectable SF area and diminished intensity in much of the surviving SF population, accompanied by a changing excitation/classification mixture in retained inner regions. Outside-in gas loss plus declining molecular regulation and continuing inner ionization is a coherent candidate framework. Its conservation equations are explicit and its optical predictions can be tested. The full observations have not yet selected a unique physical realization.

# 8. Limitations and caveats

**Cross-sectional evolution and field equivalence.** The adopted pre-peak field-like control is retained throughout. Its absolute offset from a separately selected field population is not measured here. Stage samples need not share identical histories, structures, or baseline distributions, and their labels provide no common elapsed time. If pre-peak intensity differs from the field by $\delta_{\mathrm{pre-field}}$ dex, the inferred stage/field contrast shifts by that same amount.

**Class semantics and censoring.** ND means failure of joint Balmer detection, not proven absence of gas or SF. NSF combines excitation, EW, width, and finite-measurement outcomes. The existing ND H-alpha surrogate based on a noise/flux floor is not a formal upper-limit likelihood for true SFR, particularly where the Balmer decrement is uncertain. This report does not use it as such.

**Selection and population composition.** Intensity is conditioned on surviving class membership. Changing support can change the mean without following an individual region. Stellar-density and radial bins overlap physically and share galaxies. Line-detected NSF subsets omit weak-line regions, so their ratio evolution need not describe all NSF area.

**Ionization and luminosity budgets.** Fixed continuing spectra, exponential luminosity decay, constant continuum, and lognormal scatter are illustrative assumptions. Old stars, leakage, shocks, and AGN require separate energetic, spectral, and spatial evidence. Single-screen attenuation and a single-Gaussian width can misrepresent composite spectra. The examples do not yet enforce a measured photon budget and gas emission measure at every location.

**Gas-phase closure and spatial dynamics.** One-way atomic-to-molecular transfer is an effective closure. Real phase cycling, delayed recycling, radial flows, reaccretion, ionized reservoirs, and changing cloud identities can alter it. The stripping threshold omits a full three-dimensional potential, wind geometry, magnetic fields, and shielding. The ionized emitting mass is not separately evolved in the two-neutral-reservoir example.

**Resolution and timescale mismatch.** An approximately 100-pc optical resolution element is not generally a closed regulator over hundreds of Myr. Equilibrium pressure laws calibrated on larger scales require averaging and applicability checks. H-alpha, UV, CO, and H I trace different histories and phases. A gas/SFR ratio on mismatched apertures is not automatically a depletion time for the same material.

**Statistical scope.** Bootstrap intervals describe variation among observed galaxy products; they do not remove sample-selection bias or all shared systematics. The fitting interval and model were selected after exploratory inspection. The partial fit is descriptive, with residual spatial structure, and has not undergone external validation. NGC4654's PA alignment has no calibrated significance, and its gradient mask differs from the stage-intensity mask.

**Incomplete quantitative reproduction.** The representative forward realization reproduces several qualitative trends but fails the NSF contribution decline and understates inner selected-SF suppression. That failure is reported rather than concealed. The intensity fit is a separate constrained component, not proof that the illustrative gas parameters fit the entire dataset. A defensible next result is a joint galaxy-level forward fit that can improve this realization or reject it.

# Appendix A. Data provenance and reproducibility

All source notebooks and FITS products were treated as read-only. The following live files define the observational analyses. Paths are relative to **/Users/Igniz/Desktop/ICRAR/further/**.

::: {.provenance-table}

| ID | Source |
|---|---|
| N1 | 20260909_check_combined_SF_fraction_categories_by_stage_SNR_postfit.ipynb |
| N2 | 20260909_check_Sigma_SFR_intensity_by_stage.ipynb |
| N3 | 20260909_check_corrected_Halpha_surface_density_SF_NSF_ND_by_stage.ipynb |
| N4 | 20260914_check_SF_gradient_scan_VIVA.ipynb |
| P1 | SFR+Z.py |
| C1 | mauve_master_wiki_newclass.fits |
| C2 | MAUVE_effective_radii.csv |
| D1 | Product maps under v3tk_v7.6.8/, using notebook loaders and HDU names |

:::

The six previously recorded fingerprints for N1-N3, P1, C1, and C2 were unchanged from the earlier report. Their current hashes were recorded again, but the observational quantities were nevertheless re-extracted from the live products. N4 was fingerprinted separately.

Fresh execution covered N2 data cells 3, 5, 7, 9, and 11; N3 data cells 3, 5, 7, 9, 11, and 13, plus selected central-estimator functions from cells 22, 25, and 28; and N4 cells 3, 5, 9, 11, 13, 16, 20, 28, and 30. Cell numbers are one-based serialized positions. N1's selection context was retained through the inherited helpers, rather than rerunning every N1 plotting branch. N3 received an extra direct-H-beta diagnostic only in the report extraction namespace. The full notebooks, pooled sensitivities, map galleries, and all historical validation cells were not rerun.

The output folder is **assets/20260914_resolved_RPS_academic_model/** beside this report. It contains refreshed scalar profiles, per-galaxy bin tables, bootstrap contrasts, NGC4654 plane/hemisphere tables, mathematical checks, fitted parameters, and leave-one-galaxy results. Principal scripts are **refresh_observations.py**, **inspect_gradient.py**, **analytical_predictions.py**, and **fit_attenuation.py**; the first records its adapted extraction and observational-bootstrap sources. The science interpreter is **/opt/miniconda3/envs/ICRAR/bin/python**.

The observational bootstrap uses 10,000 complete-galaxy draws per stage with seed 20260914. The fit uses the same seed in its own recorded draw sequence. Slight differences in interval endpoints from conditioning on complete finite fit vectors are expected and are not changes to central measurements. The exact SF/NSF/ND partition and strict continuum cut were checked during extraction. Invalid full-array Balmer divisions outside usable masks produced expected nonfinite values; retained diagnostic selections exclude them.

The fitting script records all seven bin contrasts, fixed weights, constant and linear predictions, bootstrap counts, parameter intervals, and leave-one-galaxy fits. Every reported fit parameter was recomputed in this report run. Mathematical and rendering checks are separate from scientific acceptance: correctly rendered equations and a stable fit do not establish unique causality.

# Appendix B. Equation attribution and numerical verification

The report distinguishes literature equations from algebraic identities and new assumptions:

::: {.equation-ledger}

| Equations | Status and source |
|---|---|
| (1)-(4), (6) | Definitions of live notebook estimators and plane fit, rewritten with explicit notation |
| (5) | Exact logarithmic derivative identity |
| (7), (15) | Molecular regulator definitions adapted from Huang et al. (2026), equations (16), (17), and (21); regulator context from Lilly et al. (2013) |
| (8)-(12) | Pressure, restoring-force, exponential-disc, and impulse approximations as used by Koppen et al. (2018) and Singh et al. (2019), retaining their stated regimes |
| (13)-(14) | Control-area mass conservation and explicit phase/transport bookkeeping |
| (16)-(17), (20)-(27) | Algebraic derivations and integrating-factor solutions under the declared regulator closures |
| (18)-(19) | New restricted two-reservoir closure and definitions of inverse-time combinations |
| (28) | Exact SFR decomposition into neutral content, molecular fraction, and depletion time |
| (29) | Blitz & Rosolowsky (2006), equation (11), with the pressure meaning retained |
| (30)-(32) | Stellar-emission convolution motivated by tracer theory; the exponential kernel and its solution are this report's approximation |
| (33)-(35) | Source addition and photon/recombination accounting; source motivation from Belfiore/Tomicic and recombination framework from Hummer & Storey |
| (36)-(41) | Flux-mixture, derivative, EW, and variance algebra; constant templates and exponential component fading are assumptions |
| (42)-(44) | Observed detection/classification definitions and idealized flux projection |
| (45)-(50) | New lognormal/fixed-threshold demonstration; CDF and truncated moments derived explicitly |
| (51)-(57) | Product condition, partial-fit definitions, least-squares derivation, and regulator inversion algebra |

:::

The reservoir and young-luminosity solutions were checked against an independent adaptive ODE integration: the maximum relative discrepancy was $5.58\times10^{-9}$. The truncated lognormal moment was checked against direct quadrature, with relative discrepancy $5.49\times10^{-16}$. Class fractions sum to unity to numerical precision and remain nonnegative. These checks verify implementation of the equations, not empirical validity of their closures.

# References

<a id="ref-huang"></a>
**Huang, R., et al. (2026).** *MAUVE-MUSE: When Metallicity Follows or Fights Star Formation--A Mass-Dependent Inversion in Virgo Galaxies.* [User-specified published DOI](https://doi.org/10.1093/mnras/stag1019); [full author manuscript, arXiv:2605.31412v1](https://arxiv.org/html/2605.31412v1). Relevant full text inspected; equation numbers refer to that manuscript.

<a id="ref-brown"></a>
**Brown, T., et al. (2023).** *VERTICO VII: Environmental quenching caused by suppression of molecular gas content and star formation efficiency in Virgo Cluster galaxies.* [Full author manuscript, arXiv:2308.10943](https://arxiv.org/html/2308.10943v1). Resolved results and their scale inspected.

<a id="ref-watts"></a>
**Watts, A. B., et al. (2023).** *VERTICO V: The environmentally driven evolution of the inner cold gas discs of Virgo cluster galaxies.* PASA. [DOI:10.1017/pasa.2023.14](https://doi.org/10.1017/pasa.2023.14); [author record](https://arxiv.org/abs/2303.07549). Primary abstract-level result used.

<a id="ref-fossati"></a>
**Fossati, M., et al. (2018).** *A Virgo Environmental Survey Tracing Ionised Gas Emission (VESTIGE). II. Constraining the quenching time in the stripped galaxy NGC4330.* A&A, 614, A57. [DOI:10.1051/0004-6361/201732373](https://doi.org/10.1051/0004-6361/201732373); [author record](https://arxiv.org/abs/1801.09685). Primary abstract used for the resolved-chronology comparison.

<a id="ref-koppen"></a>
**Koppen, J., Jachym, P., Taylor, R., & Palous, J. (2018).** *Ram Pressure Stripping Made Easy: An Analytical Approach.* MNRAS, 479, 4367. [DOI:10.1093/mnras/sty1610](https://doi.org/10.1093/mnras/sty1610); [author manuscript](https://arxiv.org/abs/1806.05887). Analytical regimes inspected for the preceding report; primary record refreshed for this revision.

<a id="ref-singh"></a>
**Singh, A., Gulati, M., & Bagla, J. S. (2019).** *Ram pressure stripping: an analytical approach.* MNRAS, 489, 5582-5593. [Full publisher text](https://academic.oup.com/mnras/article/489/4/5582/5575211); [DOI:10.1093/mnras/stz2523](https://doi.org/10.1093/mnras/stz2523). Equations and assumptions checked.

<a id="ref-lizee"></a>
**Lizee, T., Vollmer, B., Braine, J., & Nehlig, F. (2021).** *Gas compression and stellar feedback in the tidally interacting and ram-pressure stripped Virgo spiral galaxy NGC4654.* A&A, 645, A111. [DOI:10.1051/0004-6361/202038910](https://doi.org/10.1051/0004-6361/202038910); [author manuscript](https://arxiv.org/abs/2011.10531). Relevant full model and observational sections inspected.

<a id="ref-vollmer03"></a>
**Vollmer, B. (2003).** *NGC4654: gravitational interaction or ram pressure stripping?* A&A, 398, 525-540. [DOI:10.1051/0004-6361:20021729](https://doi.org/10.1051/0004-6361:20021729); [author record](https://arxiv.org/abs/astro-ph/0211321). Primary abstract and model context used.

<a id="ref-vollmer12"></a>
**Vollmer, B., Wong, O. I., Braine, J., Chung, A., & Kenney, J. D. P. (2012).** *The influence of the cluster environment on the star formation efficiency of 12 Virgo spiral galaxies.* A&A, 543, A33. [DOI:10.1051/0004-6361/201118690](https://doi.org/10.1051/0004-6361/201118690); [author manuscript](https://arxiv.org/abs/1204.0430). Relevant resolved-efficiency sections inspected.

<a id="ref-nehlig"></a>
**Nehlig, F., Vollmer, B., & Braine, J. (2016).** *Effects of environmental gas compression on the multiphase ISM and star formation.* A&A, 587, A108. [DOI:10.1051/0004-6361/201527021](https://doi.org/10.1051/0004-6361/201527021); [author manuscript](https://arxiv.org/abs/1601.04883). Relevant model and limitations inspected.

<a id="ref-br06"></a>
**Blitz, L., & Rosolowsky, E. (2006).** *The Role of Pressure in GMC Formation II: The H2-Pressure Relation.* ApJ, 650, 933-944. [DOI:10.1086/505417](https://doi.org/10.1086/505417); [author manuscript](https://arxiv.org/abs/astro-ph/0605035). Equation (11), fitted samples, and scale limitations checked.

<a id="ref-oml10"></a>
**Ostriker, E. C., McKee, C. F., & Leroy, A. K. (2010).** *Regulation of Star Formation Rates in Multiphase Galactic Disks: A Thermal/Dynamical Equilibrium Model.* ApJ, 721, 975-994. [DOI:10.1088/0004-637X/721/2/975](https://doi.org/10.1088/0004-637X/721/2/975); [author manuscript](https://arxiv.org/abs/1008.0410). Reservoir definitions and applicability inspected.

<a id="ref-fujita"></a>
**Fujita, Y., & Nagashima, M. (1999).** *Effects of Ram Pressure from the Intracluster Medium on the Star Formation Rate of Disk Galaxies in Clusters of Galaxies.* ApJ, 516, 619. [DOI:10.1086/307139](https://doi.org/10.1086/307139); [author record](https://arxiv.org/abs/astro-ph/9812378). Primary abstract-level conceptual comparison used.

<a id="ref-safarzadeh"></a>
**Safarzadeh, M., & Loeb, A. (2019).** *Explaining the enhanced star formation rate of Jellyfish galaxies in galaxy clusters.* MNRAS Letters, 486, L26-L30. [Full publisher text](https://academic.oup.com/mnrasl/article/486/1/L26/5454771); [DOI:10.1093/mnrasl/slz053](https://doi.org/10.1093/mnrasl/slz053). Pressure, cloud-efficiency, and gas-loss assumptions checked.

<a id="ref-lee"></a>
**Lee, J., Kimm, T., Katz, H., Rosdahl, J., Devriendt, J., & Slyz, A. (2020).** *Dual Effects of Ram Pressure on Star Formation in Multiphase Disk Galaxies with Strong Stellar Feedback.* ApJ, 905, 31. [DOI:10.3847/1538-4357/abc3b8](https://doi.org/10.3847/1538-4357/abc3b8); [author record](https://arxiv.org/abs/2010.11028). Primary abstract and institutional metadata checked.

<a id="ref-belfiore"></a>
**Belfiore, F., et al. (2022).** *A tale of two DIGs: The relative role of HII regions and low-mass hot evolved stars in powering the diffuse ionised gas in PHANGS-MUSE galaxies.* A&A, 659, A26. [DOI:10.1051/0004-6361/202141859](https://doi.org/10.1051/0004-6361/202141859); [author record](https://arxiv.org/abs/2111.14876). Primary abstract and previously retrieved excerpts used for source-mixture motivation.

<a id="ref-zhang"></a>
**Zhang, K., et al. (2017).** *SDSS-IV MaNGA: The Impact of Diffuse Ionized Gas on Emission-line Ratios, Interpretation of Diagnostic Diagrams, and Gas Metallicity Measurements.* [DOI:10.1093/mnras/stw3308](https://doi.org/10.1093/mnras/stw3308); [author record](https://arxiv.org/abs/1612.02000). Primary abstract-level findings used.

<a id="ref-tomicic32"></a>
**Tomicic, N., et al. (2021a).** *GASP XXXII. Measuring the diffuse ionized gas fraction in ram-pressure stripped galaxies.* ApJ. [DOI:10.3847/1538-4357/abca93](https://doi.org/10.3847/1538-4357/abca93); [author record](https://arxiv.org/abs/2011.08869). Primary abstract used for definitions and threshold limitations.

<a id="ref-tomicic35"></a>
**Tomicic, N., et al. (2021b).** *GASP XXXV: Characteristics of the diffuse ionised gas in gas-stripped galaxies.* ApJ. [DOI:10.3847/1538-4357/ac230e](https://doi.org/10.3847/1538-4357/ac230e); [author record](https://arxiv.org/abs/2108.12433). Primary record and abstract-level excitation comparison used.

<a id="ref-citro"></a>
**Citro, A., Pozzetti, L., Quai, S., Moresco, M., Vallini, L., & Cimatti, A. (2017).** *A methodology to select galaxies just after the quenching of star formation.* MNRAS. [DOI:10.1093/mnras/stx932](https://doi.org/10.1093/mnras/stx932); [author record](https://arxiv.org/abs/1704.05462). Primary abstract and stated model trend checked.

<a id="ref-lilly"></a>
**Lilly, S. J., et al. (2013).** *Gas-regulation of galaxies: the evolution of the cosmic sSFR, the metallicity-mass-SFR relation and the stellar content of haloes.* ApJ, 772, 119. [DOI:10.1088/0004-637X/772/2/119](https://doi.org/10.1088/0004-637X/772/2/119); [author record](https://arxiv.org/abs/1303.5059). Primary record used for context; local equations checked against Huang et al.

<a id="ref-ke12"></a>
**Kennicutt, R. C., & Evans, N. J. (2012).** *Star Formation in the Milky Way and Nearby Galaxies.* ARA&A, 50, 531-608. [DOI:10.1146/annurev-astro-081811-125610](https://doi.org/10.1146/annurev-astro-081811-125610); [author manuscript](https://arxiv.org/abs/1204.3552). Tracer review used as background, not as the source of the new exponential kernel.

<a id="ref-kl14"></a>
**Kruijssen, J. M. D., & Longmore, S. N. (2014).** *An uncertainty principle for star formation. I. Why galactic star formation relations break down below a certain spatial scale.* MNRAS, 439, 3239. [DOI:10.1093/mnras/stu098](https://doi.org/10.1093/mnras/stu098); [author record](https://arxiv.org/abs/1401.4459). Primary record used for cloud-scale sampling limitations.

<a id="ref-hs87"></a>
**Hummer, D. G., & Storey, P. J. (1987).** *Recombination-line intensities for hydrogenic ions. I. Case B calculations for H I and He II.* MNRAS, 224, 801-820. [DOI:10.1093/mnras/224.3.801](https://doi.org/10.1093/mnras/224.3.801); [NASA bibliographic record](https://ntrs.nasa.gov/citations/19870046128). Metadata and case-B scope verified; no numerical recombination coefficient is imported.

<a id="ref-kauffmann"></a>
**Kauffmann, G., et al. (2003).** *The host galaxies of active galactic nuclei.* MNRAS, 346, 1055-1077. [Publisher text](https://academic.oup.com/mnras/article/346/4/1055/1062435); [author record](https://arxiv.org/abs/astro-ph/0304239). Citation for the empirical [N II] boundary implemented in the live pipeline.

<a id="ref-kewley"></a>
**Kewley, L. J., et al. (2001).** *Theoretical Modeling of Starburst Galaxies.* ApJ, 556, 121-140. [DOI:10.1086/321545](https://doi.org/10.1086/321545); [author record](https://arxiv.org/abs/astro-ph/0106324). Citation for the theoretical BPT boundaries implemented in the live pipeline.
