# MAUVE: eight findings in the context of the literature

**What previous studies establish, what remains unresolved, and what the stage-resolved MAUVE data can test**

**Date:** 22 October 2026, Australia/Perth.  
**Scope:** resolved star formation, its spatial suppression and possible local enhancement, and the properties of residual non-star-forming ionised emission.  
**Evidence status:** literature review and research interpretation. The eight findings supplied by Rongjun are the observational premises; this report does not independently remeasure their amplitudes or statistical significance. Tentative findings remain tentative.

## Executive assessment

The eight findings can form a coherent paper about **the spatial disappearance of the selected star-forming population and the changing emission from the gas that remains detectable**. The literature already establishes environmental truncation of star formation, examples of compression and local enhancement, extended low-ionisation emission without an active nucleus, and dynamically disturbed stripped gas. MAUVE's strongest opportunity is to connect these phenomena within a common spatial selection and across environmental stages, while separating changes in occupied area from changes within surviving regions.

Three comparisons are especially important for positioning the work. First, [Brown et al. (2023), VERTICO VII](https://arxiv.org/abs/2308.10943), and [Cakir et al. (2026), SAMI](https://arxiv.org/abs/2602.03088), already provide direct resolved-SF evidence for environmental suppression. A broad agreement with a nearby-galaxy rSFMS must therefore be stated alongside stage-dependent departures. Second, [Brown et al. (2026), MAUVE-MUSE](https://arxiv.org/abs/2604.20056), is direct prior work on this survey's ionisation and kinematics. Findings 6 and 7 should extend its cluster-field comparison through stage and local-environment dependence; they are not a first discovery of anomalous MAUVE line ratios. Third, [Baron et al. (2026)](https://arxiv.org/abs/2610.07163), a 5 October preprint, makes hardened radiation leaking from H II regions a particularly timely alternative to an increasing HOLMES contribution. Its publication status is **submitted to ApJ**, not established peer-reviewed consensus.

The most defensible interpretation is conditional: the environmental sequence may reduce both the area and intensity of H II-selected emission; outer regions increasingly fall below Balmer detectability, while inner regions can remain detectable under a different balance of ionising sources. This can raise the NSF fraction while lowering NSF Halpha surface brightness. It does not require increasing absolute heating power. Local compression, shocks, mixing, outflows, and gravitational interactions can coexist with this broad picture, but their individual contributions need separate tests.

The most valuable new result may be **the joint dependence on stage, local stellar density, radius, Halpha brightness, and proximity to remaining SF**, rather than any single marginal trend. The [S II] doublet and line-width hints are potentially useful discriminants, but are currently the weakest links in a causal interpretation.

## 0. Scope, definitions, and how to read this review

### 0.1 Search and evidence boundaries

This is a broad, targeted literature review with an annotated source register, not a claim to have located every relevant paper. Searches were conducted on 22 October 2026 using publisher pages, arXiv author manuscripts, and institutional publication records. Coverage includes foundational work and targeted 2024-2026 updates. Citation trails were followed from PHANGS, GASP, VIVA, VERTICO, SAMI, and MAUVE papers. Reviews are navigation/context sources; specific empirical results are attributed to the primary studies.

The searches covered: resolved SF main sequence and spatial scale; Virgo and RPS stages; NGC 4654 and competing tidal/RPS models; truncated and centrally concentrated SF; DIG fractions and Halpha fading; HOLMES and photon leakage; optical diagnostic diagrams; stripped-gas kinematics; and [S II] density sensitivity. Recent searches explicitly included 2025 and 2026. Peripheral high-redshift density work is included only where it clarifies diagnostic limitations, not as an environmental analogue.

The source register distinguishes **T**, a relevant full-text passage inspected, from **A**, an abstract and bibliographic record inspected. T does not mean every page was audited. Publication years follow the journal record where verified; an earlier arXiv year is not silently used as the publication year. No numerical meta-analysis, exhaustive ADS export, formal completeness estimate, or comprehensive retraction/author-conflict audit was performed. Search-engine crawl dates were not treated as publication dates.

### 0.2 MAUVE premises and live selection context

The user's eight findings determine the order and scientific scope below. The live sources `20261021_plot_paper_2.ipynb` and `20261017_check_NSF_SII_ratio_EW_SF_proximity_sigma_by_stage.ipynb` were read to clarify selection semantics; their numerical cells were not executed. The 11 October paper-planning report was consulted for context, but its old numerical values and candidate status do not supersede the present summary.

In the current plotting notebook, the common usable footprint requires finite stellar surface density, finite gas BIN_ID, nGIST MASK equal to zero, and finite SNR_POSTFIT >25, with galaxy/component inclination strictly <80 degrees. Selected SF additionally requires finite H II SFR, proxy EW(Halpha) >6 Angstrom, and intrinsic Halpha dispersion <45 km/s. ND fails the joint Balmer-detection selection. NSF is the remaining usable population, partitioning the footprint with SF and ND. These are observational categories; NSF is not synonymous with pure DIG, pure HOLMES emission, or physically zero star formation. ND is not a measurement of zero gas mass or zero SFR.

The [S II] stage statistic in the current plotting notebook is a ratio of mean linear intensities on paired valid doublet pixels, with separate equal-galaxy and pooled variants. It is not the arithmetic mean of individual spaxel ratios. This matters when comparing published density estimates. The high-EW reference selection in the 17 October proximity notebook applies EW(Halpha) >20 Angstrom to **SF reference spaxels only**; it does not redefine the whole SF/NSF partition.

**Interpretation labels used below:** agreement means qualitative consistency under stated conditions; tension means a comparison worth testing; discrimination means a proposed measurement could separate explanations. None means an uncomputed significance test has already succeeded.

### 0.3 Keep four quantities separate

| Quantity | Scientific question | Essential boundary |
|---|---|---|
| SF-selected rSFMS | How intense is SF where it still satisfies the SF criteria? | Conditional on survival and selection. |
| SF, NSF, ND area fractions | Where does each observed category occur? | Use the same usable-area denominator within each mass/radius bin. |
| Halpha flux fraction of a component | Which component supplies the line luminosity? | Not interchangeable with its area fraction. |
| SFE = SFR / molecular-gas mass | How efficiently is molecular gas forming stars? | Requires a gas measurement and a conversion-factor assumption. |

A useful accounting identity is `selected-SF SFR / usable area = SF area fraction x mean SFR surface density within selected SF`. It describes the selected contribution, not the true total SFR if young-star emission also lies in NSF or below the detection limit. Setting ND to zero for a specifically defined lower-limit estimator does not establish a physical zero.

## 1. Finding 1: a broadly field-comparable resolved SF main sequence

**MAUVE premise:** Virgo galaxies broadly follow the field rSFMS, providing a baseline for environmental comparisons.

### What previous studies have found

The relation between local stellar and SFR surface density is well established in nearby integral-field surveys. [Cano-Diaz et al. (2016)](https://arxiv.org/abs/1602.02770) demonstrated a spatially resolved SF sequence with CALIFA. [Erroz-Ferrer et al. (2019)](https://arxiv.org/abs/1901.04493) examined SF and ionised-gas properties at much finer scales with the MUSE Atlas of Disks. [Pessa et al. (2021)](https://doi.org/10.1051/0004-6361/202140733) is the particularly relevant PHANGS comparison because it investigates resolved scaling relations from 100 pc to 1 kpc. These studies justify a nearby-disc reference, but do not make a fitted relation universal across spatial scales, galaxy populations, or SF-selection methods.

The direct environmental comparison is more nuanced. [Brown et al. (2023)](https://arxiv.org/abs/2308.10943) studied 33 Virgo satellites at 720 pc resolution and found lower resolved SFRs in H I-poor systems, linked to both molecular-gas content and efficiency. Their SFR estimator uses NUV plus mid-infrared emission and their field control is resolution matched. It is consequently not identical to an Halpha/H II-selected MAUVE relation. [Cakir et al. (2026)](https://arxiv.org/abs/2602.03088) found steeper rSFMS relations in RPS-affected SAMI populations, driven by suppression at low stellar density; the central SF activity was comparable to controls.

### What remains unclear

The unresolved issue is not whether a resolved sequence exists. It is how much of its environmental variation comes from molecular-gas availability, efficiency, the disappearance of SF sites, and the choice of tracer or spatial scale. A sample can look broadly field-like when all surviving SF regions are combined while showing meaningful offsets in particular stages. Conversely, removing low-EW and faint regions can make survivors appear deceptively normal.

PHANGS should be described precisely as a nearby star-forming comparison sample. Its membership and selection should be checked before calling it a strictly isolated, environment-free field control. Pre-peak MAUVE galaxies are an internal reference, not automatically an unprocessed field sample. Their prior group history, orbit, and selection can matter.

### What MAUVE can agree with, challenge, or add

Broad agreement supports using the familiar local SF framework as a baseline. It does **not** establish that Virgo has no effect on the gas-SF cycle. A robust stage-dependent displacement would agree with VERTICO/SAMI at a conceptual level even if the combined locus remains close to PHANGS. A disagreement would become substantive only after harmonising stellar masses, IMF, SFR calibrations, resolution, projected versus deprojected areas, SF masks, and the local-density range.

The useful MAUVE addition is to ask whether approximately normal surviving H II regions coexist with a much smaller SF occupation fraction, and whether later-stage survivors fade as well. This distinguishes two observational routes to a deficit without prematurely attributing either to a unique gas process.

**Priority comparison:** fit or measure offsets over a common stellar-density interval, retain galaxy-to-galaxy scatter, and show the conditional SF relation beside SF occupancy. A pooled spaxel relation is descriptive; whole-galaxy resampling is needed for stage inference. If a resolution-matched VERTICO comparison is undertaken, label the tracer difference explicitly.

**Paper-ready wording:** "The selected SF regions broadly occupy the nearby-disc rSFMS, while the stage-resolved distributions and spatial occupation fractions reveal environmental differences that are not captured by the combined locus."

## 2. Finding 2: resolved SF changes across RPS stages

**MAUVE premise:** pre-peak, close-to-peak, and post-peak systems differ in their resolved SF properties.

### What previous studies have found

[Chung et al. (2009), VIVA](https://doi.org/10.1088/0004-6256/138/6/1741), established the diversity of Virgo H I morphologies, including asymmetry and truncated discs. [Yoon et al. (2017)](https://arxiv.org/abs/1704.00760) connected H I morphology and gas loss with projected cluster phase space to organise stripping histories. Such classifications provide environmental context, rather than exact elapsed times for every object.

[Koopmann & Kenney (2004)](https://arxiv.org/abs/astro-ph/0406243) showed that truncated Halpha discs are a prominent route to reduced SF in Virgo. [Fossati et al. (2018), VESTIGE II](https://doi.org/10.1051/0004-6361/201732373), used the spatially resolved history of NGC 4330 to constrain outside-in quenching. [Crowl & Kenney (2008)](https://arxiv.org/abs/0807.3747) used stellar populations in stripped Virgo discs to examine recent SF truncation. These case-based time constraints complement morphology-based stage assignments.

The newer population results make a stage narrative plausible but not automatic. VERTICO VII finds a relation between H I stripping state and the resolved gas-SF cycle; SAMI's 2026 asymmetric and truncated classes show different spatial SF patterns. The broad synthesis in [Cortese, Catinella & Smith (2021)](https://arxiv.org/abs/2104.02193) emphasises how cold-gas removal and subsequent SF quenching need not occur on the same timescale. [Boselli, Fossati & Sun (2022)](https://arxiv.org/abs/2109.13614) reviews the strong dependence of stripping outcomes on orbit, gas phase, and galaxy structure.

### What remains unclear

A cross-sectional stage sequence cannot uniquely reconstruct one galaxy's trajectory. Stellar mass, bulge fraction, initial gas distribution, encounter geometry, cluster substructure, and prior processing can differ between bins. A close-to-peak sample selected by conspicuous disturbances can also overrepresent unusually active interactions. A difference between stage averages need not mean every individual galaxy evolves monotonically through those averages.

There is also an important distinction between rapid removal of diffuse outer gas and the fate of dense inner gas. The latter may continue forming stars, become inefficient, lose replenishment, flow inward, or be expelled by feedback. Optical SF maps alone do not separate these routes.

### What MAUVE can agree with, challenge, or add

The proposed ordering agrees with a well-supported environmental transformation picture. Its incremental value is strongest if the same galaxy weighting and spatial support connect rSFMS offsets, occupation fractions, and residual emission. MAUVE can test where the late-stage deficit appears first in the observed distributions: area, surviving intensity, or both.

An inner deficit stronger than that reported by Cakir et al. would be a useful conditional tension. It could reflect more advanced stripping, a different galaxy population, feedback, or different definitions. It should not be described as a contradiction until these possibilities are tested.

**Priority comparison:** retain per-galaxy summaries; match or reweight the stellar-mass and structural distributions; perform leave-one-galaxy-out checks; and vary uncertain stage assignments within plausible alternatives. Use age-sensitive stellar diagnostics to test recent quenching where feasible. Do not convert stage labels into a universal clock.

**Paper-ready wording:** "The resolved SF distributions vary systematically with the adopted environmental stages, consistent with progressive transformation, although the stage sequence is a population comparison rather than a calibrated temporal track."

## 3. Finding 3: candidate local SF enhancement, especially in NGC 4654

**MAUVE premise:** local enhancement may coexist with overall suppression; RPS-induced compression is a candidate explanation requiring stronger verification.

### What previous studies have found

NGC 4654 is not a clean single-mechanism experiment. [Vollmer (2003)](https://doi.org/10.1051/0004-6361:20021729) modelled its asymmetric gas distribution using both gravitational interaction and ram pressure. [Lizee et al. (2021)](https://arxiv.org/abs/2011.10531) directly studied gas compression and stellar feedback in this galaxy. They identified a high-H I-surface-density region and inferred enhanced molecular SFE, by roughly a factor of 1.5-2, with conclusions involving gas conversion, pressure, and feedback assumptions. This is the closest physical precedent for the proposed MAUVE feature.

There is broader evidence for enhancement during stripping: [Vulcani et al. (2018)](https://arxiv.org/abs/1810.05164) found elevated SF in their stripping sample, and [Vulcani et al. (2020)](https://arxiv.org/abs/2007.04996) examined the resolved enhancement. However, enhancement is not universal. [Vollmer et al. (2012)](https://doi.org/10.1051/0004-6361/201118690) found no significant systematic windward enhancement of molecular SFE in their Virgo sample. [Tonnesen & Bryan (2012)](https://arxiv.org/abs/1203.0308) showed that stripping and tail SF need not produce a global burst. [Goller et al. (2023)](https://doi.org/10.1093/mnras/stad2551) found no population-wide SF enhancement for simulated jellyfish galaxies, while allowing individual episodes.

### What remains unclear

At least four propositions are different: one side has more SF than the other; one side lies above the matched reference rSFMS; gas accumulation increases SFR at unchanged efficiency; and SFE itself increases. An azimuthal asymmetry establishes only the first unless its absolute amplitude is measured. A high SFR surface density alone cannot distinguish more molecular gas from a shorter depletion time.

NGC 4654 also contains spiral structure and a tidal history, both capable of producing spatial asymmetry. The current opposite-tail or compression direction is an external physical clue with projection uncertainty; it should not define the region used to discover the excess and then be treated as an independent confirmation.

### What MAUVE can agree with, challenge, or add

A significant absolute excess spatially coincident with the compressed gas would support the Lizee et al. interpretation at a different resolution and with ionised-gas/stellar diagnostics. It would not independently prove that RPS is the sole cause. If the excess disappears relative to a matched reference, the correct conclusion is weaker suppression or redistribution, even if the map is strikingly one-sided.

MAUVE can contribute an objective, blind spatial test: find coherent positive residuals relative to a fixed reference; quantify uncertainty using independent spatial elements; then compare their location with H I, CO, dust, and stellar structure. A galaxy-wide null result does not rule out a small compressed sector, just as a local excess does not reverse a total deficit.

**Priority tests:** establish a positive absolute residual, account for searching many sectors or thresholds, compare equivalent radii and stellar densities, and use spatially correlated rather than independent-pixel errors. Only then inspect directional agreement with the external gas morphology. CO is required to distinguish gas accumulation from enhanced SFE. Test whether a stellar asymmetry or tidal model explains the same region.

**Paper-ready wording:** "NGC 4654 contains a candidate local SF excess that is compatible with previously reported compression, but its absolute significance and the relative roles of ram pressure and tides remain to be established."

## 4. Finding 4: different spatial outcomes of declining SF fraction

**MAUVE premise:** later stages have lower SF fraction; outer/low-density locations show more ND, while inner/high-density locations show more NSF.

### What previous studies have found

Virgo H I and Halpha truncation studies provide strong precedent for outer-disc loss, while resolved IFU work has established that detectable non-H II-like emission can persist where recent SF is weak. [Belfiore et al. (2017)](https://doi.org/10.1093/mnras/stw3211) distinguished central and extended LIER systems and linked their emission to different SF configurations. [Lacerda et al. (2018)](https://arxiv.org/abs/1711.07844) found systematic changes in low-EW and SF-dominated emission across morphology and radius. Their categories provide context, not a one-to-one mapping onto MAUVE's NSF mask.

The most important counter-comparison is [Tomicic et al. (2021a), GASP XXXII](https://arxiv.org/abs/2011.08869). In 38 stripped and 33 control galaxies, they found no clear stripped-control difference in the inferred DIG contribution to integrated emission, despite broad variation among galaxies. They also cautioned against a universal single brightness, ratio, or EW boundary for DIG separation. This is not a direct refutation of increasing NSF **area** fraction: the quantity, classification, spatial footprint, and stage comparison are different.

### What remains unclear

The spatial partition alone does not identify a physical transition path. An increased NSF fraction need not imply that NSF gas is being created, accumulating in mass, or becoming more luminous. It can arise because a previously dominant H II contribution fades while other emission remains detectable. An increase in ND can reflect declining emission measure, fewer ionising photons, lower gas covering fraction, or a detection threshold. Even at the same radius, stellar density and galaxy type are coupled.

The central NSF excess may also include pre-existing bulges, AGN, bars, outflows, or projection through diffuse layers. A stage-dependent change in the mixture of these structures can imitate a change caused by environmental processing.

### What MAUVE can agree with, challenge, or add

This finding provides a useful connection between outer truncation and the residual inner ionised medium. Its strongest formulation is that **the loss of selected SF has different observable counterparts in different parts of the disc**. It can add more than a single radial SFR profile because it keeps detectable non-SF-like emission separate from non-detection.

A genuine empirical addition would be a residual stage dependence in the NSF/ND partition at matched local stellar density and radius. A two-dimensional comparison is important: separate one-dimensional profiles cannot by themselves distinguish a density effect from a radial effect. The available common support should be reported rather than extrapolating into bins occupied by only one stage.

**Priority tests:** retain the exact common denominator; tabulate galaxy counts and category counts; test sensitivity to Balmer depth and SF selection; and examine the reasons each NSF spaxel failed the SF criteria. Compare both occupation and luminosity fractions if making a GASP comparison. Where possible, degrade the deeper data to a common detection limit as a check.

**Paper-ready wording:** "The declining selected-SF fraction is accompanied primarily by increasing Balmer non-detection in outer/low-density regions and by increasing detectable NSF emission in inner/high-density regions. These are distinct observational outcomes, not a direct measurement of gas-phase conversion."

## 5. Finding 5: fading Halpha emission in NSF regions

**MAUVE premise:** NSF Halpha becomes fainter from pre-peak to post-peak, consistent with declining contributions from recent SF.

### What previous studies have found

Diffuse Halpha is connected to the radiation field and the structure of the ionised medium. [Haffner et al. (2009)](https://arxiv.org/abs/0901.0941) reviews the observational and theoretical basis. [Belfiore et al. (2022)](https://arxiv.org/abs/2111.14876) modelled the PHANGS DIG using leakage from H II regions together with evolved-star contributions to harder spectral features. [Flores-Fajardo et al. (2011)](https://arxiv.org/abs/1104.0525) explored OB stars plus HOLMES in NGC 891. [Cid Fernandes et al. (2011)](https://arxiv.org/abs/1012.4426) and [Byler et al. (2019)](https://arxiv.org/abs/1904.10978) show why weak emission can persist after vigorous SF ceases, without requiring an AGN.

These papers motivate a fading-young-source interpretation, but they do not establish the specific MAUVE stage ordering. Their primary contribution here is to supply physically plausible sources of the residual emission and observational constraints on their budgets.

The retired-galaxy interpretation also has an important earlier precedent: [Stasinska et al. (2008)](https://arxiv.org/abs/0809.1341) combined stellar-population fits and photoionisation models to show that old populations could explain some apparently LINER-like systems. The fraction explained remained uncertain. Residual emission after SF declines is therefore an established possibility; its spatial and environmental dependence is the question for MAUVE.

### What remains unclear

For optically thin recombination emission at a specified temperature, Halpha surface brightness depends on the line-of-sight integral of electron and proton density times the effective recombination coefficient. Fainter emission can therefore result from lower density, smaller emitting path length or filling factor, fewer absorbed ionising photons, or some combination. Dust and stellar-absorption corrections also matter for weak lines. This is an emission-measure constraint, not by itself a gas-mass measurement.

An NSF-only comparison follows a changing selected population. The regions in the pre-peak NSF sample are not necessarily the ancestors of those in the post-peak NSF sample. Moreover, imposing Balmer detection removes the faintest regions at later stages. A trend among detections may underestimate fading across the full footprint, while its precise magnitude can still be affected by selection.

### What MAUVE can agree with, challenge, or add

Fading NSF emission is compatible with reduced leakage from declining SF and with a larger **relative** contribution from slowly changing old-star radiation. It does not establish increasing absolute HOLMES luminosity. A rising NSF fraction plus lower NSF brightness is physically and mathematically consistent: classification area and emissivity are different measurements.

MAUVE can test whether the fading follows nearby SF activity, stellar surface density, or changing gas conditions. Use Halpha luminosity per unit old stellar mass as a budget diagnostic, EW as a line-to-continuum diagnostic, and spatial association with surviving SF as a leakage diagnostic. Their predictions overlap, so their joint behaviour is more useful than any one axis.

For a stellar-population prediction, the full-absorption HOLMES luminosity is a **maximum available photon-budget estimate under the adopted population model**. A predicted/observed ratio below one can constrain how much Halpha that model could supply; a ratio above one means the budget is sufficient, not that HOLMES dominance is demonstrated. Escaping photons, dust absorption, geometry, and mass normalisation must be accounted for before interpreting this ratio as an actual contribution.

**Priority tests:** compare at common local density and radius; retain upper limits or explicitly model the detection selection; check low-EW stellar-continuum fits; and plot Halpha alongside EW and a young-source radiation proxy. An Halpha profile alone cannot distinguish gas removal from a fading radiation field.

**Paper-ready wording:** "The NSF Halpha emission fades toward later stages, consistent with a declining young-star contribution and/or lower emission measure; the data do not yet uniquely separate changes in the radiation field from changes in the emitting gas."

## 6. Finding 6: increasing excitation-sensitive line ratios in NSF regions

**MAUVE premise:** optical BPT ratios generally increase from pre-peak to post-peak.

### What previous studies have found

[Zhang et al. (2017)](https://arxiv.org/abs/1612.02000) showed that diffuse, low-Halpha-surface-brightness emission can shift common optical ratios away from the H II-region sequence. [Belfiore et al. (2016)](https://doi.org/10.1093/mnras/stw1234) demonstrated the widespread spatial extent of LIER emission and evidence for diffuse stellar ionisation. These results are central reasons not to label every non-SF BPT location as AGN or shock excitation. [Mannucci et al. (2021)](https://arxiv.org/abs/2109.02684) further shows that aperture choices can account for some differences previously attributed to DIG contamination.

For stripping specifically, [Tomicic et al. (2021b), GASP XXXV](https://arxiv.org/abs/2108.12433) found strong [O I] excesses in diffuse tail gas, with additional ionisation/heating implicated. [Fossati et al. (2016)](https://arxiv.org/abs/1510.04283) found that a mixture of mechanisms was needed in ESO 137-001 rather than a pure shock spectrum. Tail results should not be transferred wholesale to retained inner-disc gas.

The closest survey precedent is [Brown et al. (2026)](https://arxiv.org/abs/2604.20056): the early MAUVE sample shows elevated optical ratios and broader gas kinematics relative to PHANGS. The paper favours the greater visibility of diffuse emission as SF is suppressed, with an additional smaller shock-like contribution. Its cluster-field result is closely related to the present hypotheses, but is not itself a controlled pre-/close-/post-peak comparison.

### The active uncertainty: HOLMES versus spectrally hardened leakage

The PHANGS framework of Belfiore et al. (2022) invokes old-star radiation to help reproduce the harder DIG spectrum. The newly submitted [Baron et al. (2026)](https://arxiv.org/html/2610.07163v1) instead argues that propagation-induced hardening of leaked H II radiation can reproduce the joint high- and low-ionisation trends in actively star-forming discs, reducing the need for HOLMES there. It retains a role for HOLMES in sufficiently retired regions. This is a live modelling disagreement, not a settled replacement of one universal explanation by another.

An increasing [N II]/Halpha or [S II]/Halpha ratio is not a direct measurement of increasing ionisation parameter. Harder spectra, lower ionisation parameter, abundance variations, temperature changes, geometry, mixing, and shocks can move diagnostics differently. [O III]/Hbeta is especially important because its response to dilution depends on the spectrum and physical regime. The low-ionisation BPT axis and the [S II] doublet ratio must not be confused.

### What MAUVE can agree with, challenge, or add

The stage trend is compatible with a changing source mixture, but it cannot identify HOLMES dominance from BPT position alone. It can test whether the aggregate change is largely explained by sampling fainter gas or whether later-stage NSF gas remains spectrally different at fixed Halpha brightness, local stellar density, radius, and SF proximity.

The latter comparison is particularly useful for the new leakage debate. Estimate the available radiation from **all relevant neighbouring SF regions**, allowing attenuation and geometry, rather than using only the nearest-region distance. Then ask whether residual line-ratio changes track old stellar mass, the remaining SF radiation field, or kinematic disturbances. The >20 Angstrom proximity-reference cut can help define luminous young-source locations, but does not itself prove the completeness of the leaked-photon budget.

If a pure HOLMES budget is insufficient where ratios are most extreme, that challenges that particular population/absorption model; it does not automatically establish shocks. If leakage is quantitatively insufficient far from surviving SF and weak-EW emission follows old stars, a HOLMES interpretation becomes more credible. A model should reproduce multiple ratios and the absolute Halpha luminosity simultaneously.

**Priority tests:** common-detection-support and censoring checks, because weak Hbeta and [O III] can bias ratio distributions; matched-brightness comparisons; joint density-radius conditioning; and separate central AGN/outflow regions. For ratios of summed fluxes, retain that estimator when comparing models. For additive emitting components, mix linear fluxes before taking logarithms. Mixing incident spectra in one gas cloud requires a photoionisation calculation, not just a straight line in log-BPT space.

**Paper-ready wording:** "The NSF line-ratio distributions shift toward later stages, consistent with a changing mixture of ionising radiation and gas conditions. Joint constraints from Halpha brightness, EW, SF proximity, and kinematics are required to distinguish hardened leakage, evolved-star ionisation, and mechanical excitation."

## 7. Finding 7: a possible increase in NSF velocity dispersion

**MAUVE premise:** NSF emission lines may be slightly broader at later stages.

### What previous studies have found

[Law et al. (2021)](https://arxiv.org/abs/2011.06012) demonstrated a strong connection between optical diagnostic classification and ionised-gas dynamics in MaNGA. This provides a general precedent for broader non-H II-like emission, without making environment its unique cause. The MAUVE cluster-field differences in Brown et al. (2026) are the closest prior comparison.

Stripped tails provide clear examples of disturbed gas: [Fumagalli et al. (2014)](https://arxiv.org/abs/1407.7527) found a transition to broader gas motions downstream in ESO 137-001, and [Luo et al. (2023)](https://arxiv.org/abs/2212.03891) extended the kinematic mapping over much longer tails. Such results establish the plausibility of wake-driven disturbance. Their spatial regime and line widths should not be used as a direct quantitative prediction for a slight shift within MAUVE discs. [Allen et al. (2008)](https://arxiv.org/abs/0805.0204) supplies shock-emission grids; model shock speed is not the same observable as a fitted line-of-sight dispersion.

### What remains unclear

A broader fitted line can reflect turbulence, several unresolved velocity components, beam-smeared rotation, extraplanar gas, stellar feedback, or gravitational disturbance. It does not uniquely measure mechanical heating or shock velocity. An instrumental correction removes the adopted line-spread contribution, but does not remove thermal broadening or spatial/line-of-sight mixing.

There is a specific selection issue in the current analysis: SF requires intrinsic Halpha dispersion <45 km/s. Consequently, some broad Balmer-detected regions are assigned to NSF because of their width. No additional <45 cut is imposed when measuring NSF dispersion, but the classification is already partly width-dependent. An NSF-SF contrast is therefore not independent evidence of a dynamical difference. Even a stage comparison within NSF can be affected if the mix of reasons for failing SF changes with stage.

[Sanchez et al. (2024)](https://arxiv.org/abs/2311.10573) explicitly developed the WHaD classification using Halpha EW and dispersion. This is relevant methodological precedent for combining those observables, but does not validate every numerical boundary or make a subsequent width contrast independent of the classification.

### What MAUVE can agree with, challenge, or add

A small, reproducible later-stage broadening would support changes in the observed ionised-gas dynamics. It would strengthen a mechanical interpretation only if associated with independent signatures such as coherent non-circular velocities, asymmetric profiles, [O I] excess, or spatially identified interfaces. Smooth low-EW emission with modest widths is also compatible with a growing relative contribution from diffuse layers.

Published MAUVE outflows make sample-specific alternatives essential. [Watts et al. (2024)](https://doi.org/10.1093/mnras/stae898) identified a large bipolar outflow in NGC 4383. [Attwater et al. (2026)](https://doi.org/10.3847/2041-8213/ae2c51) found a star-formation-driven outflow in the stripped NGC 4064. [Cortese et al. (2026)](https://arxiv.org/abs/2602.01612) adds multiphase evidence for fountain behaviour in NGC 4383. These demonstrate that external processing and internal feedback can act together; a broad component is not uniquely an RPS signature.

**Priority tests:** repeat the dynamical comparison using an excitation/EW selection independent of sigma, while preserving the original selection as the primary definition; inspect instrumental-resolution uncertainty and positive-radicand selection; control inclination, radial velocity gradients, S/N, and Halpha brightness; and separate outflows/AGN. Use galaxy-level errors and independent resolution elements, not the raw spaxel count.

**Paper-ready wording:** "There is tentative evidence for broader NSF line profiles in later-stage systems, but the small effect, instrumental sensitivity, projection, and the width-dependent SF selection require further checks before interpreting it as enhanced turbulence or shocks."

## 8. Finding 8: possible changes in NSF electron density

**MAUVE premise:** in high-stellar-density regions, pre-peak galaxies have lower [S II] 6717/6731 ratios than later stages, potentially indicating higher electron density.

### What previous studies have found

At fixed electron temperature and atomic data, the [S II] 6717/6731 ratio decreases as density rises through its useful sensitivity range. [Luridiana, Morisset & Shaw (2015)](https://arxiv.org/abs/1410.6662) describes PyNeb's atomic-level treatment and the importance of the input atomic data. [Kewley et al. (2019)](https://arxiv.org/abs/1908.05504) shows why converting diagnostics into density or pressure also depends on the thermal and ionisation structure; calibrations for star-forming nebulae are not universal NSF models.

[Sanders et al. (2016)](https://doi.org/10.3847/0004-637X/816/1/23) discusses the proximity of many local measurements to the low-density limit. [Barnes et al. (2021)](https://doi.org/10.1093/mnras/stab2958), using PHANGS H II regions, explicitly distinguishes the [S II]-weighted density from a volume-averaged density. The commonly quoted limits are approximately 1.4-1.5 at low density and 0.4-0.5 at high density; exact limits and useful sensitivity depend on temperature and atomic data. These papers calibrate interpretation, not the claimed MAUVE stage trend.

There are direct density measurements in stripped systems, but the emitting population matters. [Poggianti et al. (2019), GASP XIII](https://arxiv.org/abs/1811.00823), studied star-forming clumps in 16 galaxies' tails and obtained a median density of about 52 cm^-3 for the measured tail-clump sample. The doublet constrained density in only about half of the clumps. These compact H II regions are different from MAUVE's inner-disc NSF population. Their existence is compatible with diffuse gas dominating elsewhere; it does not predict a universal density increase during stripping.

### What remains unclear

If both stage distributions are consistent with the low-density plateau, a small ratio shift may not resolve different physical densities. Noise, correlated flux errors, sky residuals, line-detection cuts, and the nonlinear inversion can become more important than the inferred difference. Values outside theoretical limits should not simply be clipped into apparently precise densities.

The current stage estimator adds a further boundary: a ratio of intensity averages or sums produces an effective emission-weighted diagnostic. Inverting it does not yield the mean spaxel density, and certainly not the mean density of all retained gas. Unresolved dense structures can contribute disproportionately. The [S II]-emitting low-ionisation zone may sample a different medium in different stages.

### What MAUVE can agree with, challenge, or add

The direction reported by the user is correct for a density interpretation: **lower 6717/6731 in pre-peak implies potentially higher, not lower, electron density**, under matched conditions. It is compatible with later-stage emission being weighted toward a more diffuse phase or with reduced density in the emitting gas. It would be uncomfortable for a simple picture in which the whole later-stage inner NSF medium becomes uniformly denser under compression, but does not rule out compressed clouds or local shocks embedded in diffuse gas.

This distinction is scientifically useful. A lower effective density together with broader profiles would not be paradoxical: velocity dispersion and local electron density measure different properties. Nor can lower density alone explain all BPT shifts. The ionisation parameter depends on both photon flux and gas density, so both have to be constrained.

**Priority tests:** lead with the measured doublet ratio and uncertainty; jointly model the two line fluxes; preserve their shared footprint and covariance; quantify what fraction of the posterior lies in the sensitive regime; and propagate plausible temperature and atomic-data choices. Report density limits when appropriate. Compare paired flux-sum ratios with the distribution of individually constrained ratios as separate statistics. There is no basis here for deriving a dense molecular-gas mass or a unique ISM pressure from this doublet alone.

**Paper-ready wording:** "At high stellar surface density, the NSF [S II] doublet ratio is tentatively higher in later-stage systems, consistent with lower effective electron density if the measurements resolve the diagnostic's sensitive regime."

## 9. Where apparent agreements and disagreements become scientifically useful

### 9.1 Comparison matrix

The judgements below are conditional on the user's findings. They are not new statistical measurements.

| MAUVE finding | Closest established comparison | Present assessment | What would make the comparison decisive? |
|---|---|---|---|
| 1. Broad rSFMS agreement | Pessa 2021; Cano-Diaz 2016; VERTICO III/VII | Agreement with a local sequence can coexist with environmental offsets. | Same scale, selection, calibration, and density interval. |
| 2. Stage-dependent SF | Yoon 2017; Fossati 2018; Brown 2023; Cakir 2026 | Consistent with environmental processing; exact timing is unresolved. | Robust stage contrasts after galaxy matching and uncertain-stage checks. |
| 3. Local enhancement | Lizee 2021; Vulcani 2020; Vollmer 2012; Foster 2025 | Plausible and previously reported in NGC 4654; universality is unsupported. | Absolute positive excess, spatial significance, CO comparison, tidal alternative. |
| 4. Outer ND / inner NSF | Koopmann & Kenney 2004; Belfiore 2017; Tomicic 2021a | Useful spatial extension; NSF area cannot be equated with DIG flux fraction. | Joint density-radius comparison with matched detection depth. |
| 5. NSF Halpha fading | Belfiore 2022; Tomicic 2024; Byler 2019 | Fits reduced young-star contribution or lower emission measure. | Photon budget, EW, UV/Halpha, and gas constraints. |
| 6. Higher BPT ratios | Brown 2026; Zhang 2017; Belfiore 2022; Baron 2026 | Strong precedent; origin remains non-unique. | Residual stage shifts at fixed brightness and source budgets. |
| 7. Slightly broader NSF | Law 2021; Brown 2026; ESO 137-001 studies | Dynamical hint; not yet a shock measurement. | Width-independent selection and resolved kinematic evidence. |
| 8. Higher late-stage [S II] ratio | Barnes 2021; Kewley 2019; Sanders 2016 | Possible lower effective density; sensitive to saturation and weighting. | Flux-space inference and demonstrated sensitivity away from the plateau. |

### 9.2 Five comparisons that should not be presented as simple contradictions

**Field-like relations versus environmental suppression.** [Jimenez-Donaire et al. (2023)](https://arxiv.org/abs/2211.16521) provides a useful example: the ensemble Virgo molecular KS relation resembles that of isolated galaxies while H I-deficient subsets have lower efficiency. The KS relation and rSFMS are distinct, but the statistical lesson applies: ensemble resemblance does not preclude conditional environmental trends. MAUVE should show the strata rather than using the combined locus as a veto on suppression.

**Enhancement versus suppression.** [Roberts & Parker (2020)](https://arxiv.org/abs/2004.12033) and [Vulcani et al. (2024)](https://arxiv.org/abs/2311.13486) strengthen the evidence for enhanced activity in selected stripping populations. [Foster et al. (2025)](https://arxiv.org/abs/2502.13123) found no significant population-wide asymmetric enhancement in a much broader satellite sample. This difference plausibly involves selection, tracer timescale, event duration, and spatial averaging. Those explanations should be tested rather than assumed to resolve all tensions. MAUVE's candidate is a local, individual-galaxy claim.

**Compression versus gas transport.** [Zhu, Tonnesen & Bryan (2024)](https://arxiv.org/abs/2309.07037) attributes enhancement in its simulations to ram-pressure-driven flows that increase central dense-gas supply. Thus, even if environmental enhancement is confirmed, direct windward cloud compression is not the only route. [Steinhauser et al. (2016)](https://arxiv.org/abs/1604.05193) likewise illustrates that different stripping regimes can have different SF outcomes. Spatial geometry and gas dynamics are needed to distinguish these mechanisms.

**More NSF area versus unchanged DIG luminosity fraction.** An area fraction weights faint and bright selected pixels alike; a flux fraction is dominated by luminosity. Moreover, the MAUVE NSF category includes several reasons for failing SF selection. A quantitative GASP comparison should recreate both the measurement and the population, not compare the words "diffuse" or "non-star-forming."

**Higher ratios versus more heating.** A ratio can rise because its Balmer denominator declines faster than the forbidden line. The emission can become fainter overall. Demonstrating more mechanical energy injection requires an energy budget and dynamical context, not only a movement on the BPT diagram. Conversely, an adequate old-star photon budget is necessary for a pure HOLMES interpretation but is not sufficient evidence that the photons are absorbed in the observed gas.

### 9.3 A particularly useful bridge: UV versus Halpha

[Tomicic et al. (2024)](https://doi.org/10.3847/1538-4357/ad7130) compared UV- and Halpha-derived SF in five late-stage GASP systems. DIG-dominated locations can have higher UV-based estimates, which the authors interpret in terms of recent declines in SF; they also identify Halpha-only regions with non-SF-like excitation. This is a more direct temporal comparison for finding 5 than simply interpreting lower Halpha as less recent SF.

For MAUVE, UV/Halpha can help test fading because the tracers respond differently to the recent stellar population. It still needs attenuation, resolution, and source-contamination control. It must not convert NSF Halpha to an SFR as if all of it were powered by young stars. The useful observable can be the relative distribution of UV, Halpha, stellar ages, and excitation, even where a formal SFR ratio is not physically appropriate.

## 10. A connected physical scenario and its alternatives

### 10.1 The economical working scenario

The following is a synthesis proposed for testing, not a unique inversion of the observations.

1. Environmental processing reduces the supply and spatial extent of gas capable of supporting the selected H II population. The observed SF fraction falls, particularly outside the retained inner disc. Surviving regions can also move below the reference rSFMS.
2. Some locations experience local gas accumulation or compression during the interaction, producing short-lived enhancement amid a larger deficit. Gravitational interactions and internal gas flows can contribute.
3. At low stellar density and large radius, declining emission increasingly fails the Balmer detection requirement. In the inner footprint, gas can remain detectable after it ceases to meet the SF selection.
4. Reduced young-star emission changes the relative weights of leaked photons, evolved-star radiation, and any mechanically excited component. The NSF area may rise while its mean Halpha brightness falls.
5. The line-ratio, line-width, and doublet-ratio changes constrain the physical state of this residual emitting population, but do not by themselves tell us whether each change reflects the same parcels of gas or a changing mixture of phases.

This working scenario is compatible with the broad picture reviewed by Cortese et al. (2021) and Boselli et al. (2022). Its distinctive test in this project is the linked spatial accounting. The current evidence does not determine a unique stripping rate, photon escape fraction, gas depletion time, or mechanical-energy injection rate.

### 10.2 Alternative explanations and discriminating observations

| Explanation | What it could reproduce | Observation that would strengthen it | Main ambiguity |
|---|---|---|---|
| Declining SF plus more visible diffuse emission | Less SF area; fainter Halpha; altered ratios | Changes largely follow Halpha/EW and recent-SF indicators | Does not identify the diffuse ionising source. |
| Greater relative HOLMES contribution | Weak-EW emission associated with old stars | Sufficient absorbed photon budget and old-star spatial scaling, including far from SF | Stellar-population and covering-factor uncertainty. |
| Hardened leaking H II radiation | Diffuse emission connected non-locally to SF; elevated ratios | Radiation-field model reproduces brightness and multiple ratios simultaneously | Three-dimensional transport and attenuation are poorly known. |
| Mechanical excitation / mixing | Line-ratio changes with disturbed profiles | Co-spatial kinematic components, interfaces, [O I], and an energy budget | Beam smearing and feedback can mimic parts of the signal. |
| Feedback-assisted late quenching | Central loss of retained gas; outflows | Resolved bipolar geometry and multiphase mass/velocity constraints | Outflow mass and whether gas escapes are uncertain. |
| Gas rarefaction or changing filling factor | Fainter recombination emission and a higher doublet ratio | Density-sensitive measurements away from the low-density limit | Doublet-weighted density does not measure volume or gas mass. |
| Selection and changing galaxy composition | Several apparent stage trends | Trends weaken after matching, forward selection modelling, or leave-one-out tests | Real environmental and selection effects may coexist. |

**A useful falsification target:** if the stage-dependent BPT differences disappear after matching Halpha brightness and the available young-source radiation field, there is less need for a separate stage-dependent excitation mechanism. If they persist, the result justifies additional physics or a different source mixture, but still does not uniquely select shocks or HOLMES. Conditioning on Halpha is a diagnostic decomposition, not a causal estimate: Halpha itself may be affected by the environmental process.

## 11. What remains open, and what MAUVE may contribute

### 11.1 Established results that should be acknowledged as precedent

Outside-in truncation, environmental variation in resolved SF, individual enhancement during stripping, DIG/LIER emission in weakly star-forming systems, and disturbed stripped-gas kinematics all have substantial prior literature. The broad effect of environment on MAUVE excitation and dynamics has also already been published. These are the starting points for the paper.

Historical optical classification papers, including [Kewley et al. (2001)](https://arxiv.org/abs/astro-ph/0106324) and [Kewley et al. (2006)](https://arxiv.org/abs/astro-ph/0605681), remain important for the diagnostic boundaries. Their labels do not uniquely identify the source of spatially extended emission; the later resolved LIER literature is essential when interpreting them.

### 11.2 Strong candidate contributions of this analysis

**A stage-resolved decomposition of SF suppression.** Measure occupancy and survivor intensity with consistent support, and determine whether their relative contributions depend on stellar density and radius. This can connect the rSFMS to the spatial loss of the selected population more directly than either statistic alone.

**The spatial partition of residual detectability.** Establish whether the outer-ND/inner-NSF pattern remains after jointly controlling density, radius, structure, and sensitivity. This can constrain where line emission survives, even when the ionising source is uncertain.

**A controlled test of what changes within NSF.** Determine which spectral changes remain at matched brightness and which arise because later-stage NSF occupies a different part of parameter space. A robust residual would add to Brown et al. (2026); a null residual would also be informative because it would support redistribution among existing emission regimes.

**A testbed between active discs and retired systems.** The remaining SF in stripped galaxies provides a useful range of local radiation fields. Combined photon-budget, spatial, and EW measurements can test the relative roles of leakage and evolved stars in this intermediate regime. A claim of decisively solving the debate would require forward models and uncertainties beyond the present eight statements.

**A focused NGC 4654 case study.** The contribution would be an independently defined, quantified optical residual and its relation to previous gas results. It should explicitly acknowledge the earlier compression/SFE study and mixed tidal/RPS history.

**A density/kinematics constraint on the residual medium.** This is potentially novel in the precise stage-and-location form, but should remain a secondary result until sensitivity and selection are demonstrated. The literature search alone cannot certify a global "first" claim.

### 11.3 Prioritised analyses, tied to questions rather than checklist volume

| Priority | Analysis | Question answered | Interpretation if the result is null |
|---|---|---|---|
| 1 | Per-galaxy occupancy and survivor-intensity contrasts on common support | Where does the observed SF deficit reside? | A combined rSFMS offset may be dominated by population composition. |
| 2 | Joint local-density/radius and galaxy-structure controls | Does stage retain information beyond location and sample makeup? | The marginal stage trend need not reflect a local environmental effect. |
| 3 | NSF diagnostics matched in Halpha, EW, and observational depth | Does NSF change internally or through selection into different regimes? | Fading/reweighting may explain much of the aggregate spectral change. |
| 4 | HOLMES budget plus spatially propagated young-source budget | Which sources can energetically supply the observed emission? | A plausible spectrum without enough photons is insufficient. |
| 5 | Absolute, blind NGC 4654 residual test | Is there enhancement rather than weaker suppression? | Retain asymmetry as descriptive; drop the stronger enhancement claim. |
| 6 | Width-independent excitation sample and kinematic component checks | Is broader NSF independent of its definition and instrumental effects? | Do not use width to support environmental turbulence. |
| 7 | Doublet-flux likelihood and temperature/atomic-data sensitivity | Is the density contrast actually resolved? | Report ratios and limits rather than a measured density evolution. |
| 8 | Matched H I/CO and UV/stellar-population comparisons | Is gas loss, efficiency change, or fading responsible? | Optical categories alone cannot complete the causal account. |

These are proposed analyses, not checks executed for this review. They should preserve the existing primary masks; alternative definitions should be explicitly labelled sensitivity tests. A galaxy with many spaxels is not many independent galaxies, and a combined product must not be silently counted as several independent measurements.

## 12. Suggested discussion sequence

Begin with the field-comparable surviving-SF baseline, then show how the distributions and occupied area differ across stages. Place local enhancement immediately after the broader suppression result to show that the two can coexist. Use the outer-ND/inner-NSF partition as the bridge to the residual ionised medium. Discuss Halpha fading before line ratios, since brightness changes provide an essential baseline for interpreting those ratios. End the observational sequence with the tentative kinematic and doublet-ratio results.

The physical discussion can then ask, in order: how much is explained by loss/fading of the selected H II population; how much is explained by reweighting existing emission sources; what requires changes in the gas itself; and what external data can distinguish gas removal, suppressed efficiency, leakage, HOLMES, and mechanical excitation. This order allows the stronger results to carry the narrative while the candidate enhancement and weak density/dynamical trends remain useful tests.

**Suggested central statement:** "MAUVE links the declining spatial occupation and changing intensity of selected star formation to the detectability and spectral properties of the residual ionised medium across Virgo environmental stages. The combined trends favour spatially selective suppression and a changing balance of emission sources, while leaving the roles of photon leakage, evolved stars, and mechanical excitation to be constrained jointly."

## Appendix A. Annotated reference register

The register includes core and supporting papers so that the discussion can be expanded without treating every reference as equally direct. Each entry gives its principal use and a boundary. **F1-F8** refer to the eight findings; **T/A** describe access as defined in Section 0. Journal details are supplied where checked; the stable primary link is the retrieval reference. Author lists are abbreviated for readability. The annotations are evidence summaries, not quotations.

### A1. Resolved SF baselines and environmental dependence

**R01. Cano-Diaz et al. (2016).** *Spatially-Resolved Star Formation Main Sequence of Galaxies in the CALIFA Survey.* ApJL, 821, L26. [Primary manuscript](https://arxiv.org/abs/1602.02770); DOI: 10.3847/2041-8205/821/2/L26. **F1; A; observed population.** Establishes a local SF sequence at roughly kiloparsec scales. Use as a baseline precedent; its slope is selection- and resolution-specific.

**R02. Erroz-Ferrer et al. (2019).** *The MUSE Atlas of Disks (MAD): resolving star formation rates and gas metallicities on <100 pc scales.* MNRAS, 484, 5009-5027. [Primary manuscript](https://arxiv.org/abs/1901.04493). **F1, F6; A; observed population.** Provides nearby-disc resolved SF and emission-line context. Fine spatial sampling does not eliminate diffuse-emission or population-selection issues.

**R03. Pessa et al. (2021).** *Star formation scaling relations at ~100 pc from PHANGS: Impact of completeness and spatial scale.* A&A, 650, A134. [Publisher](https://doi.org/10.1051/0004-6361/202140733); [manuscript](https://arxiv.org/abs/2104.09536). **F1, F4; A; observed population.** Direct reference for scale/completeness dependence; all three resolved scaling relations are recovered at 100 pc. Non-detection treatment affects fitted slopes.

**R04. Vulcani et al. (2019).** *GASP - XX. From the loose spatially resolved to the tight global SFR-mass relation in local spiral galaxies.* MNRAS, 488, 1597-1617. [Publisher](https://doi.org/10.1093/mnras/stz1829). **F1; T; observed control sample.** Demonstrates substantial internal and between-galaxy variation. The 30 undisturbed galaxies include both field and cluster systems; "undisturbed" is not equivalent to "field."

**R05. Ellison et al. (2018).** *Star formation is boosted (and quenched) from the inside-out: radial star formation profiles from MaNGA.* MNRAS, 474, 2039-2054. [Primary manuscript](https://arxiv.org/abs/1711.00915); DOI: 10.1093/mnras/stx2882. **F1, F2, F4; A; observed population.** Shows how global SF offsets map into radial structure. It is an alternative internal/population context, not evidence that RPS must quench inside-out.

**R06. Medling et al. (2018).** *The SAMI Galaxy Survey: Spatially Resolving the Main Sequence of Star Formation.* [Primary manuscript](https://arxiv.org/abs/1801.04283); DOI: 10.1093/mnras/sty127. **F1, F2, F4; A; observed population.** Links resolved SF patterns to environment and galaxy properties. SF-selected samples are clean but incomplete representations of all SF-related emission.

**R07. Schaefer et al. (2017).** *The SAMI Galaxy Survey: spatially resolving the environmental quenching of star formation in GAMA galaxies.* [Primary manuscript](https://arxiv.org/abs/1609.02635); DOI: 10.1093/mnras/stw2289. **F2, F4; A; observed population.** Denser environments show more concentrated SF and steeper profiles, with substantial scatter. Broad environment does not isolate a single gas-removal mechanism.

**R08. Brown et al. (2023).** *VERTICO VII: Environmental quenching caused by suppression of molecular gas content and star formation efficiency in Virgo Cluster galaxies.* ApJ, 956, 37. [Primary manuscript](https://arxiv.org/abs/2308.10943); DOI: 10.3847/1538-4357/acf195. **F1-F4; T; direct Virgo comparison.** Finds suppression associated with both molecular content and SFE; early-stage outer enhancement can instead reflect more molecular gas. NUV+WISE3 at 720 pc and a HERACLES control differ from H II-selected optical measurements.

**R09. Jimenez-Donaire et al. (2023).** *VERTICO III: The Kennicutt-Schmidt relation in Virgo cluster galaxies.* A&A, 671, A3. [Primary manuscript](https://arxiv.org/abs/2211.16521); DOI: 10.1051/0004-6361/202244718. **F1-F3; A; direct Virgo comparison.** The ensemble molecular KS relation resembles the field, with lower efficiency in H I-deficient subsets. A KS result is not itself an rSFMS measurement.

**R10. Cakir et al. (2026).** *The SAMI Galaxy Survey: Quenching of Star Formation in Clusters III. Ram-Pressure-Affected Galaxy Populations.* PASA, 43, e029. [Primary manuscript](https://arxiv.org/abs/2602.03088); DOI: 10.1017/pasa.2026.10157. **F1-F4; T; observed population.** RPS candidates show low-density/outer suppression and steeper rSFMS; central activity is comparable to controls. Morphological classes are not interchangeable with the MAUVE stage labels.

### A2. Virgo stripping, spatial truncation, and temporal context

**R11. Koopmann & Kenney (2004a).** *Massive Star Formation Rates and Radial Distributions from Halpha Imaging of 84 Virgo Cluster and Isolated Spiral Galaxies.* ApJ, 613, 851-865. [Primary manuscript](https://arxiv.org/abs/astro-ph/0209547); DOI: 10.1086/423190. **F2, F4; A; direct Virgo comparison.** Outer truncation accounts for much of the integrated SF deficit, while inner SF can persist. The initial preprint is from 2002; the publication is 2004.

**R12. Koopmann & Kenney (2004b).** *Halpha Morphologies and Environmental Effects in Virgo Cluster Spiral Galaxies.* ApJ, 613, 866-885. [Primary manuscript](https://arxiv.org/abs/astro-ph/0406243); DOI: 10.1086/423191. **F2-F4; A; direct Virgo comparison.** Establishes the prevalence and diversity of truncated SF morphologies. Mixed tidal and hydrodynamic effects remain possible.

**R13. Chung et al. (2009).** *VLA Imaging of Virgo Spirals in Atomic Gas (VIVA). I. The Atlas and the H I Properties.* AJ, 138, 1741-1816. [Publisher](https://doi.org/10.1088/0004-6256/138/6/1741); [manuscript](https://arxiv.org/abs/0909.0781). **F2-F4; A; observed atlas.** The key gas-morphology and kinematics reference. The original article DOI is given here; the arXiv related-link metadata can instead point to its later erratum.

**R14. Yoon et al. (2017).** *A History of H I Stripping in Virgo: A Phase-space View of VIVA Galaxies.* [Primary manuscript](https://arxiv.org/abs/1704.00760); DOI: 10.3847/1538-4357/aa6579. **F2, F4; A; observed population plus orbital interpretation.** Associates stripping morphologies with infall/backsplash context. Projection and multiple histories limit individual time assignments.

**R15. Vollmer (2009).** *A holistic view on ram pressure stripping in the Virgo cluster - The first complete model-based time sequence.* A&A, 502, 427-435. [Publisher](https://doi.org/10.1051/0004-6361/200911892); [manuscript](https://arxiv.org/abs/0905.1770). **F2; T; dynamical modelling.** Provides a model-based sequence for six Virgo galaxies. Useful for timing methodology, not a universal empirical clock.

**R16. Crowl & Kenney (2008).** *The Stellar Populations of Stripped Spiral Galaxies in the Virgo Cluster.* AJ, 136, 1623-1644. [Primary manuscript](https://arxiv.org/abs/0807.3747); DOI: 10.1088/0004-6256/136/4/1623. **F2-F4; A; stellar-population constraints.** Recent outer-disc truncation with surviving central SF; strong outer bursts are not generally required. Does not exclude a small local enhancement in another part of a galaxy.

**R17. Fossati et al. (2018).** *A Virgo Environmental Survey Tracing Ionised Gas Emission (VESTIGE). II. Constraining the quenching time in the stripped galaxy NGC 4330.* A&A, 614, A57. [Publisher](https://doi.org/10.1051/0004-6361/201732373); [manuscript](https://arxiv.org/abs/1801.09685). **F2, F4; A; individual case.** Resolves an outside-in history with retained inner SF. Times depend on the fitted history and geometry.

**R18. Cortese, Catinella & Smith (2021).** *The Dawes Review 9: The role of cold gas stripping on the star formation quenching of satellite galaxies.* PASA, 38, e035. [Primary review](https://arxiv.org/abs/2104.02193). **F1-F5; A; review.** Framework for relating gas loss to delayed or incomplete SF suppression. Use to locate primary evidence, not as a substitute for the direct comparisons above.

**R19. Boselli, Fossati & Sun (2022).** *Ram pressure stripping in high-density environments.* A&A Review, 30, 3. [Primary review](https://arxiv.org/abs/2109.13614); DOI: 10.1007/s00159-022-00140-3. **F2-F8; A; review.** Broad multiphase environmental context and remaining uncertainties. Disc and tail regimes must be distinguished.

### A3. Enhancement and its counterexamples

**R20. Vollmer (2003).** *NGC 4654: gravitational interaction or ram pressure stripping?* A&A, 398, 525-540. [Primary manuscript](https://arxiv.org/abs/astro-ph/0211321); DOI: 10.1051/0004-6361:20021729. **F3; A; individual dynamical model.** Combined tides and RPS reproduce the observed system better than either alone in the tested models. The causal conclusion is model-dependent.

**R21. Lizee et al. (2021).** *Gas compression and stellar feedback in the tidally interacting and ram-pressure stripped Virgo spiral galaxy NGC 4654.* A&A, 645, A111. [Primary manuscript](https://arxiv.org/abs/2011.10531); DOI: 10.1051/0004-6361/202038910. **F3; A; direct same-galaxy precedent.** Reports compressed gas and enhanced molecular efficiency in the high-H I-density region. SFE and an rSFMS residual have different denominators, and alpha_CO matters.

**R22. Vollmer et al. (2012).** *The influence of the cluster environment on the star formation efficiency of 12 Virgo spiral galaxies.* A&A, 543, A33. [Primary manuscript](https://arxiv.org/abs/1204.0430); DOI: 10.1051/0004-6361/201118690. **F3; T; observed population.** No significant general windward molecular-SFE enhancement, despite gas-fraction and total-gas-efficiency changes in some regions. Retain the gas denominator when citing SFE.

**R23. Vulcani et al. (2018).** *Enhanced star formation in both disks and ram pressure stripped tails of GASP jellyfish galaxies.* [Primary manuscript](https://arxiv.org/abs/1810.05164); DOI: 10.3847/2041-8213/aae68b. **F2, F3; A; selected stripping population.** Reports elevated integrated disc SF and additional tail SF. Selected jellyfish are not a representative census of all satellites.

**R24. Vulcani et al. (2020).** *GASP XXX. The spatially resolved SFR-Mass relation in stripping galaxies in the local universe.* [Primary manuscript](https://arxiv.org/abs/2007.04996); DOI: 10.3847/1538-4357/aba4ae. **F1-F3; T; observed population.** Resolved enhancement relative to controls at approximately kiloparsec scales. Truncated systems require qualification; do not infer equal enhancement in every stripping stage.

**R25. Roberts & Parker (2020).** *Ram pressure stripping candidates in the Coma Cluster: Evidence for enhanced star formation.* MNRAS, 495, 554-569. [Primary manuscript](https://arxiv.org/abs/2004.12033); DOI: 10.1093/mnras/staa1213. **F3; A; observed candidates.** Adds evidence for enhanced SF in selected disturbed systems. Coma and candidate-selection effects limit direct comparison to Virgo.

**R26. Vulcani et al. (2024).** *Evidence for enhanced star formation rates in z ~0.35 cluster galaxies undergoing ram pressure stripping.* A&A, 682, A117. [Primary manuscript](https://arxiv.org/abs/2311.13486); DOI: 10.1051/0004-6361/202348135. **F3; A; observed population.** Extends enhancement evidence to a different epoch and cluster population. Integrated disc-plus-tail and local measurements should be distinguished.

**R27. Foster et al. (2025).** *No Evidence of Asymmetrically Enhanced Star Formation in Infalling Galaxies in UNIONS.* [Primary manuscript](https://arxiv.org/abs/2502.13123), marked accepted by ApJ; journal volume/DOI not verified here. **F3; A; broad satellite population.** Finds no significant population-wide asymmetric enhancement with u-band indicators. Does not exclude uncommon, short-lived, or smaller-scale Halpha excesses.

**R28. Goller et al. (2023).** *Jellyfish galaxies with the IllustrisTNG simulations - No enhanced population-wide star formation according to TNG50.* MNRAS, 525, 3551-3570. [Publisher](https://doi.org/10.1093/mnras/stad2551); [manuscript](https://arxiv.org/abs/2304.09199). **F3; T; cosmological simulation.** Individual bursts coexist with no ensemble enhancement relative to matched analogues. Subgrid physics and the comparison population matter.

**R29. Tonnesen & Bryan (2012).** *Star formation in ram pressure stripped galactic tails.* MNRAS, 422, 1609-1624. [Primary manuscript](https://arxiv.org/abs/1203.0308); DOI: 10.1111/j.1365-2966.2012.20737.x. **F2, F3; T; simulation.** Moderate stripping can truncate disc SF without a global burst. This is one regime, not a general prohibition on enhancement.

**R30. Zhu, Tonnesen & Bryan (2024).** *When and How Ram Pressure Stripping in Low-mass Satellite Galaxies Enhances Star Formation.* ApJ, 960, 54. [Primary manuscript](https://arxiv.org/abs/2309.07037); DOI: 10.3847/1538-4357/acfe6f. **F2, F3; A; simulation.** Central gas supply can rise through RPS-driven flows, producing enhancement. Distinguishes transport from a simple direct-compression interpretation.

**R31. Steinhauser, Schindler & Springel (2016).** *Simulations of ram-pressure stripping in galaxy-cluster interactions.* [Primary manuscript](https://arxiv.org/abs/1604.05193); DOI: 10.1051/0004-6361/201527705. **F2, F3; A; simulation.** Different wind strengths produce enhancement or rapid quenching in the studied configurations. Numerical and initial-condition dependence limits literal transfer to MAUVE stages.

### A4. Residual emission, DIG, and ionising sources

**R32. Haffner et al. (2009).** *The warm ionized medium in spiral galaxies.* Reviews of Modern Physics, 81, 969. [Primary review](https://arxiv.org/abs/0901.0941). **F5-F8; T; review.** Physical framework for diffuse gas, line ratios, radiation transport, and emission measure. Background source; not direct evidence of a Virgo stage sequence.

**R33. Cid Fernandes et al. (2011).** *A comprehensive classification of galaxies in the SDSS: How to tell true from fake AGN?* [Primary manuscript](https://arxiv.org/abs/1012.4426); DOI: 10.1111/j.1365-2966.2011.18244.x. **F4-F6; A; classification and population modelling.** Uses EW plus line ratio to distinguish weak nuclear activity from retired systems. Integrated thresholds are useful diagnostics, not universal physical purity criteria for individual resolved pixels.

**R34. Flores-Fajardo et al. (2011).** *Ionization of the diffuse gas in galaxies: Hot low-mass evolved stars at work.* MNRAS, 415, 2182-2192. [Primary manuscript](https://arxiv.org/abs/1104.0525); DOI: 10.1111/j.1365-2966.2011.18848.x. **F5, F6, F8; A; photoionisation models of NGC 891.** Demonstrates a viable mixed OB/HOLMES scenario. Geometry and abundance assumptions are specific to the modelled extraplanar gas.

**R35. Belfiore et al. (2016).** *SDSS IV MaNGA - spatially resolved diagnostic diagrams: a proof that many galaxies are LIERs.* MNRAS, 461, 3111-3134. [Publisher](https://doi.org/10.1093/mnras/stw1234). **F4-F7; T; observed population.** Extended weak emission and stellar-population/spatial evidence favour distributed sources in many systems. The nuclear term LINER should not imply an AGN origin for extended LIER emission.

**R36. Belfiore et al. (2017).** *SDSS-IV MaNGA - the spatially resolved transition from star formation to quiescence.* MNRAS, 466, 2570. [Publisher](https://doi.org/10.1093/mnras/stw3211); [manuscript](https://arxiv.org/abs/1609.01737). **F2, F4-F6; T; observed population.** Relates central/extended LIER patterns to SF and structural properties. This is not an RPS-selected sequence.

**R37. Zhang et al. (2017).** *SDSS-IV MaNGA: The Impact of Diffuse Ionized Gas on Emission-line Ratios, Interpretation of Diagnostic Diagrams, and Gas Metallicity Measurements.* MNRAS, 466, 3217-3243. [Primary manuscript](https://arxiv.org/abs/1612.02000); DOI: 10.1093/mnras/stw3308. **F5, F6; A; observed population.** Shows systematic DIG-related ratio shifts at low Halpha brightness. Important for selection and excitation interpretation; metallicity science itself is outside this report's scope.

**R38. Lacerda et al. (2018).** *Diffuse ionized gas in galaxies across the Hubble sequence at the CALIFA resolution.* MNRAS, 474, 3727-3739. [Primary manuscript](https://arxiv.org/abs/1711.07844); DOI: 10.1093/mnras/stx3022. **F4-F6; A; observed population.** Identifies low-EW, mixed, and SF-complex regimes with systematic radial/morphological variation. Its EW categories differ from MAUVE's multi-condition SF mask.

**R39. Byler et al. (2019).** *Self-consistent predictions for LIER-like emission lines from post-AGB stars.* AJ, 158, 2. [Primary manuscript](https://arxiv.org/abs/1904.10978); DOI: 10.3847/1538-3881/ab1b70. **F5, F6; A; stellar and photoionisation models.** Predicts LIER-like spectra and weak Halpha EW for old populations. Demonstrates plausibility rather than proving that any particular MAUVE region is powered this way.

**R40. Tomicic et al. (2021a).** *GASP. XXXII. Measuring the Diffuse Ionized Gas Fraction in Ram-pressure-stripped Galaxies.* ApJ, 907, 22. [Primary manuscript](https://arxiv.org/abs/2011.08869); DOI: 10.3847/1538-4357/abca93. **F4-F6; A; observed stripped/control population.** No clear stripped-control difference in inferred DIG flux fraction; a universal single diagnostic threshold is inadequate. Flux fractions cannot be directly compared with NSF area fractions.

**R41. Tomicic et al. (2021b).** *GASP XXXV: Characteristics of the diffuse ionised gas in gas-stripped galaxies.* ApJ, 922, 131. [Primary manuscript](https://arxiv.org/abs/2108.12433); DOI: 10.3847/1538-4357/ac230e. **F5-F7; A; observed population.** Diffuse tail gas shows strong [O I] excess and evidence for additional processes. Tail-specific heating is not a prediction that retained discs are shock-dominated.

**R42. Mannucci et al. (2021).** *The diffuse ionized gas (DIG) in star-forming galaxies: the influence of aperture effects on local H II regions.* [Primary manuscript](https://arxiv.org/abs/2109.02684). **F1, F6; A; observational comparison.** Aperture differences explain part of some line-ratio discrepancies attributed to DIG. Useful counterweight to automatic contamination corrections.

**R43. Belfiore et al. (2022).** *A tale of two DIGs: The relative role of H II regions and low-mass hot evolved stars in powering the diffuse ionised gas (DIG) in PHANGS-MUSE galaxies.* A&A, 659, A26. [Full manuscript](https://arxiv.org/html/2111.14876v3); DOI: 10.1051/0004-6361/202141859. **F5, F6; T, especially Sections 3.2 and 5.3-5.5; observations plus models.** Key leakage/HOLMES comparison. Its galaxy-integrated budget must not be assigned unchanged to every local NSF region.

**R44. Tomicic et al. (2024).** *Spatially resolved comparison of SFRs from UV and Halpha in GASP gas-stripped galaxies.* ApJ, 976, 90. [Primary manuscript](https://arxiv.org/abs/2409.01351); DOI: 10.3847/1538-4357/ad7130. **F2, F5, F6; A; five late-stage systems.** Tracer differences support recent SF fading in DIG-dominated locations. Non-SF-powered Halpha is not a reliable young-star SFR measurement.

**R45. Baron et al. (2026).** *DIGging into H II regions: the role of ionizing photon leakage and spectral hardening in coupling H II regions and diffuse ionized gas.* [Full manuscript](https://arxiv.org/html/2610.07163v1). **F5, F6; T, especially Section VI.3; submitted preprint, 5 October 2026.** Important current alternative to the two-source interpretation in active discs. Model-dependent; its conclusions must not be generalised automatically to retired or strongly stripped regions.

### A5. Gas dynamics, outflows, and direct MAUVE precedents

**R46. Allen et al. (2008).** *The MAPPINGS III Library of Fast Radiative Shock Models.* ApJS, 178, 20-55. [Primary manuscript](https://arxiv.org/abs/0805.0204); DOI: 10.1086/589652. **F6, F7; A; model grid.** Shock and precursor spectra span a range of physical parameters. Observed sigma is not a grid shock velocity.

**R47. Fumagalli et al. (2014).** *MUSE sneaks a peek at extreme ram-pressure stripping events. I. A kinematic study of the archetypal galaxy ESO137-001.* [Primary manuscript](https://arxiv.org/abs/1407.7527). **F4, F7; A; individual tail.** Retained rotation and broader downstream motions show the dynamical complexity of stripping. A tail is not directly analogous to an inner-disc NSF bin.

**R48. Fossati et al. (2016).** *MUSE sneaks a peek at extreme ram-pressure stripping events. II. The physical properties of the gas tail of ESO137-001.* [Primary manuscript](https://arxiv.org/abs/1510.04283). **F5-F8; A; individual tail.** Line ratios favour multiple emission mechanisms, with ordinary H II-like knots embedded in diffuse gas. The arXiv submission is 2015; the journal publication is 2016.

**R49. Law et al. (2021).** *SDSS-IV MaNGA: Refining Strong Line Diagnostic Classifications Using Spatially Resolved Gas Dynamics.* [Primary manuscript](https://arxiv.org/abs/2011.06012). **F6, F7; A; observed population.** Connects gas line width with optical diagnostic location and other stellar/galaxy properties. The correlation is not unique to environmental disturbance.

**R50. Luo et al. (2023).** *Tracing the kinematics of the whole ram pressure stripped tails in ESO 137-001.* [Primary manuscript](https://arxiv.org/abs/2212.03891). **F7; A; individual extended tail.** Maps the persistence and loss of rotation and increased dispersion along the wake. Geometrical and hydrodynamic interpretations remain model-dependent.

**R51. Watts et al. (2024).** *MAUVE: a 6 kpc bipolar outflow launched from NGC 4383, one of the most H I-rich galaxies in the Virgo cluster.* MNRAS, 530, 1968-1983. [Publisher](https://doi.org/10.1093/mnras/stae898). **F3, F6, F7; A; individual MAUVE galaxy.** Establishes a feedback-driven outflow that can affect excitation and line profiles. Important sample-specific alternative to direct RPS excitation.

**R52. Attwater et al. (2026).** *MAUVE-MUSE: A Star-formation-driven Outflow Caught in the Act of Quenching the Stripped Virgo Galaxy NGC 4064.* ApJL, 996, L40. [Publisher](https://doi.org/10.3847/2041-8213/ae2c51); [manuscript](https://arxiv.org/abs/2512.10574). **F2, F4, F7; A; individual MAUVE galaxy.** Demonstrates feedback operating after substantial environmental stripping. The published year is 2026, despite the 2025 preprint.

**R53. Cortese et al. (2026).** *MAUVE: Cold neutral gas in the outflow of NGC 4383 and evidence for a fountain flow.* PASA, 43, e034. [Primary manuscript](https://arxiv.org/abs/2602.01612); DOI: 10.1017/pasa.2026.10159. **F4, F7; A; multiphase case.** Constrains which gas phases participate in the outflow and the possible fountain outcome. Ionised outflow visibility is not equivalent to permanent removal of the cold reservoir.

**R54. Brown et al. (2026).** *MAUVE-MUSE: Ionization and Kinematic Signatures of Environmental Effects on Virgo Cluster Disks.* ApJ, 1004, 73. [Full manuscript](https://arxiv.org/html/2604.20056v1); DOI: 10.3847/1538-4357/ae632e. **F4-F7; T; direct survey precedent.** Closest prior study for excitation/kinematics. Its binned-map, equal-galaxy summaries and selection must be matched before comparing amplitudes with the current analysis.

### A6. Density and diagnostic interpretation

**R55. Luridiana, Morisset & Shaw (2015).** *PyNeb: a new tool for analyzing emission lines. I. Code description and validation of results.* A&A, 573, A42. [Primary manuscript](https://arxiv.org/abs/1410.6662); DOI: 10.1051/0004-6361/201323152. **F8; A; atomic diagnostic method.** Provides density/temperature inference and configurable atomic data. A successful inversion does not guarantee that a noisy ratio constrains density.

**R56. Sanders et al. (2016).** *The MOSDEF Survey: Electron Density and Ionization Parameter at z ~2.3.* ApJ, 816, 23. [Publisher](https://doi.org/10.3847/0004-637X/816/1/23); [manuscript](https://arxiv.org/abs/1509.03636). **F8; T; diagnostic comparison.** Relevant here for low-density-limit behaviour and local comparison samples. High-redshift SF nebulae are not direct analogues of Virgo NSF gas.

**R57. Kewley et al. (2019).** *Theoretical ISM Pressure and Electron Density Diagnostics for Local and High-redshift Galaxies.* ApJ, 880, 16. [Primary manuscript](https://arxiv.org/abs/1908.05504); DOI: 10.3847/1538-4357/ab16ed. **F6, F8; A; photoionisation diagnostics.** Shows the importance of temperature, abundance, and different emitting zones in converting ratios to density/pressure. Do not apply a star-forming-nebula calibration uncritically to mixed NSF gas.

**R58. Barnes et al. (2021).** *Comparing the pre-SNe feedback and environmental pressures for 6000 H II regions across 19 nearby spiral galaxies.* MNRAS, 508, 5362-5389. [Publisher](https://doi.org/10.1093/mnras/stab2958); [manuscript](https://arxiv.org/abs/2110.05492). **F8; T; observed population plus physical inference.** Explicit discussion of doublet sensitivity and emission-weighted/clumpy densities. Its H II-region model should not be equated with a volume average of NSF gas.

**R59. Kewley, Nicholls & Sutherland (2019).** *Understanding Galaxy Evolution Through Emission Lines.* ARA&A, 57, 511-570. [Primary review](https://doi.org/10.1146/annurev-astro-081817-051832). **F6-F8; T; review.** Summarises source, temperature, abundance, density, and selection degeneracies. Provides the diagnostic framework rather than an environmental finding.

**R60. Kewley et al. (2001).** *Theoretical Modeling of Starburst Galaxies.* ApJ, 556, 121. [Primary manuscript](https://arxiv.org/abs/astro-ph/0106324); DOI: 10.1086/321545. **F6; A; theoretical classification.** Basis for commonly used model demarcations. Crossing a boundary does not identify which non-young-star source supplies the emission.

**R61. Kewley et al. (2006).** *The Host Galaxies and Classification of Active Galactic Nuclei.* MNRAS, 372, 961-976. [Primary manuscript](https://arxiv.org/abs/astro-ph/0605681); DOI: 10.1111/j.1365-2966.2006.10859.x. **F6; A; empirical classification.** Important Seyfert/LINER diagnostic context. Later spatially resolved work changes how extended LIER emission should be interpreted physically.

### A7. Additional historical and selection-specific precedents

**R62. Stasinska et al. (2008).** *Can retired galaxies mimic active galaxies? Clues from the Sloan Digital Sky Survey.* MNRAS, 391, L29-L33. [Primary manuscript](https://arxiv.org/abs/0809.1341); DOI: 10.1111/j.1745-3933.2008.00550.x. **F4-F6; A; stellar-population and photoionisation modelling.** Early demonstration that old-star-powered emission can occupy the apparent AGN wing. The fraction of systems explained was uncertain; the result does not identify the source in an individual MAUVE spaxel.

**R63. Sanchez et al. (2024).** *WHaD diagram: Classifying the ionizing source with one single emission line.* A&A, 682, A71. [Primary manuscript](https://arxiv.org/abs/2311.10573); DOI: 10.1051/0004-6361/202347711. **F4, F6, F7; A; diagnostic classification.** Uses Halpha EW and line width to classify emission. Relevant selection context, with the consequent need to avoid treating a width-selected contrast as independent physical evidence.

**R64. Poggianti et al. (2019).** *GASP XIII. Star formation in gas outside galaxies.* MNRAS, 482, 4466-4502. [Primary manuscript](https://arxiv.org/abs/1811.00823); DOI: 10.1093/mnras/sty2999. **F3, F5-F8; T, abstract and density conclusions; observed tail population.** Measures physical properties of star-forming tail clumps and distinguishes diffuse tail emission. Density-constrained clumps are not a representative measurement of all stripped gas, and are not equivalent to inner-disc NSF regions.

## Appendix B. Source coverage, provenance, and limitations

### B1. Search clusters and retained evidence

| Search cluster | Representative queries | Main source families retained |
|---|---|---|
| Resolved SF baseline | "PHANGS rSFMS completeness spatial scale"; "CALIFA resolved star formation main sequence" | Pessa, Cano-Diaz, MAD, SAMI, GASP control |
| Environmental SF and stages | "Virgo resolved star formation VERTICO"; "SAMI ram pressure resolved star formation 2026" | VERTICO III/VII, SAMI 2026, VIVA, Yoon, VESTIGE |
| Enhancement and null results | "NGC4654 gas compression star formation"; "ram pressure no enhanced star formation" | Lizee, Vollmer, GASP, Coma, UNIONS, TNG50, wind simulations |
| Residual ionised emission | "GASP diffuse ionized gas"; "PHANGS HOLMES leakage" | Tomicic, Belfiore, Zhang, Lacerda, Byler, Baron |
| Kinematics and feedback | "stripped ionized gas velocity dispersion"; "MAUVE ionization kinematic" | Law, ESO 137-001, Brown 2026, MAUVE outflows |
| Electron density | "S II low density limit"; "electron density diagnostics Kewley" | PyNeb, Sanders, Kewley, Barnes |

The register contains **64 distinct papers**. It includes sources that qualify or oppose simple enhancement/HOLMES/shock narratives. It deliberately excludes metallicity-science extensions, unrelated high-redshift environmental studies, and papers identified only through unverifiable secondary claims. A comprehensive census would require a reproducible bibliographic database export, citation-network screening, and wider full-text review. This report is a substantial evidence map, not such a census.

### B2. Local evidence inspected

- `/Users/Igniz/Desktop/ICRAR/AGENTS.md`: project interpretation and preservation rules.
- `/Users/Igniz/Desktop/ICRAR/further/20261021_plot_paper_2.ipynb`: source definitions of common support, SF/NSF/ND, and diagnostic estimators.
- `/Users/Igniz/Desktop/ICRAR/further/20261017_check_NSF_SII_ratio_EW_SF_proximity_sigma_by_stage.ipynb`: source definitions of the high-EW SF proximity reference and NSF diagnostics.
- `/Users/Igniz/Desktop/ICRAR/MAUVE/20261011_MAUVE_Paper_II_Findings_Synthesis_and_Paper_Plan.md`: contextual reading only; old numbers were not adopted as fresh results.
- `/Users/Igniz/Desktop/ICRAR/MAUVE/20261001_Connected_HI_Stripping_and_HOLMES_Line_Ratio_Model.md`: targeted reference pointers; the relevant Belfiore primary text was checked afresh.

Prior report-workflow memory was used to preserve the evidence/literature/interpretation distinction and the Markdown/PDF delivery convention. It was not used as current scientific evidence. No FITS data, notebooks, palettes, masks, or pipeline products were edited. No scientific pipeline or significance calculation was run for this review.

### B3. Claim boundaries to retain in a manuscript

Do not replace "candidate local enhancement" with a confirmed RPS-triggered burst. Do not translate increasing NSF area into increasing gas mass or heating power. Do not equate declining Halpha with a unique gas-loss timescale. Do not equate non-H II BPT classification with pure HOLMES, shocks, or an AGN. Do not call a stage sequence a longitudinal observation of the same gas. Do not claim density evolution when the doublet remains consistent with its low-density limit. Do not claim a first MAUVE ionisation/kinematics result without addressing Brown et al. (2026).

The report was prepared with AI-assisted source discovery, synthesis, and document production, including one bounded bibliography/references review worker. The main agent retained scientific synthesis and final acceptance. Literature access and uncertainty are recorded above so that human collaborators can trace and refine the argument before manuscript use.
