# MAUVE Paper II: what we have found, how it fits together, and how to write the paper

**Research synthesis and paper plan | 11 October 2026**

**Scope:** spatially resolved star formation and ionized emission across the pre-peak, close-to-peak, and post-peak environmental stages. Stellar mass, star-forming area, surviving-SF intensity, spatial structure, and excitation diagnostics are included. Gas-phase metallicity science is reserved for Paper III.

**Status:** a source-based synthesis of the current 8900 Angstrom analysis, with a focused live reproduction of the central SF results. This is a planning report, not a completed manuscript or a validation of the forthcoming 7000 Angstrom reduction. Notebook identifiers N001-N128 refer to the complete inventory in Appendix A; R identifiers refer to the dated reports in Appendix B. Cell numbers are one-based in the files inspected on this date.

## 1. The paper already present in the results

The strongest coherent paper is about **where star formation survives, how strongly those surviving regions form stars, and what emission remains outside the selected star-forming population as Virgo galaxies progress through environmental stages**. An rSFMS comparison alone misses much of this transformation because it measures only the locations that still satisfy the star-forming selection.

The current evidence supports three connected empirical statements.

1. **Selected star-forming area declines strongly toward post-peak.** In the strict-continuum-S/N sample, the equal-product mean SF fraction across each product's usable footprint is 0.494, 0.332, and 0.101 for pre-, close-, and post-peak respectively. The resolved profiles show that this is spatially structured, with especially strong outer-disc/low-density loss and substantial suppression even in the inner footprint. These central values were reproduced from the live FITS files for this report.
2. **The remaining selected SF is also weaker, particularly at fixed stellar surface density.** The later analyses support an additional intensity shift beyond disappearance from the SF selection. However, its strength depends on the estimator and location: sparse outer survivors do not show uniform fading. A September fit summarizes the post-peak intensity offset as about -0.39 dex over its fitted density interval; that fit is a descriptive normalization change, not a measured gas-consumption timescale.
3. **The emission left outside the selected SF population changes in both spatial occupancy and spectral character.** Outer regions increasingly fail joint Balmer detection, while central/high-density regions often remain Balmer-detected but fail the SF criteria. NSF Halpha emission can become fainter while its local area fraction increases, and its BPT distributions shift. This does not mean NSF grows everywhere, or that the harder-looking emission is uniquely powered by shocks, AGN, or old stars.

The most useful conceptual link is an exact observational accounting:

> HII-traced SFR per usable area = SF occupancy x intensity within the surviving SF area.

The paper can therefore show two ways that an environment-associated SFR deficit appears in resolved data, then ask what the residual ionized emission tells us about the changing population. A gas-loss plus fading-young-source picture can connect the observations, but the present work does **not** establish one complete physical mechanism. The simplest fixed-spectrum HOLMES model has quantitative tensions, and the strongest local enhancement candidate, NGC4654, has not passed its newest absolute-enhancement test.

**Recommended central claim:** Across the adopted Virgo environmental stages, the detectable HII-region population becomes spatially depleted and, in the common stellar-density range, fainter among its survivors. The remaining Balmer-detected population develops different excitation distributions. These changes are consistent with spatially selective environmental suppression and a changing balance of ionizing sources, without yet identifying a unique gas-loss or excitation mechanism.

**Recommended working title:** *MAUVE-MUSE: Where star formation survives and fades across Virgo environmental stages.* A more method-led alternative is *MAUVE-MUSE: Separating the loss of star-forming area from the fading of surviving regions in Virgo galaxies.* Neither title depends on proving universal enhancement or HOLMES dominance.

## 2. What was reviewed and what was verified

### 2.1 Coverage of the workspace

The recursive inventory found **128 notebooks**, including two undated notebooks, spanning 15 March to 10 October 2026. All 128 were parsed, their cell sources and text outputs extracted, their execution/output/error status inventoried, and their file hashes recorded. The review concentrated on the stage-resolved mass/SFR, category, emission-line, spatial, and gas-response branches. Metallicity-only notebooks were inventoried but their science was deliberately excluded. Earlier reduction/QC and computing notebooks were used as methodological provenance, not as extra environmental findings.

This is a broad workspace review with targeted scientific depth; it is **not** a line-by-line code audit of every notebook or an execution of all 128. A saved execution count or a PASS string is evidence of a stored run, not proof that the current source and current inputs have all been rerun together. The dated filenames identify the notebook lineage, not an immutable snapshot: several files have later edits and copied headings.

The supporting survey also indexed **82 dated MAUVE Markdown reports from April through October**. The most informative are the August occupancy/distribution sequence, the September 9-14 measurement and interpretation reports, and the October excitation/model updates. The October 9 report named `20261009 .md` is empty. Image-only notes were not treated as independently quantified results. The NGC4254 spiral work is a separate exploratory branch; the newest environmental result is the October 10 NGC4654 notebook.

The broad folder also contains the active mass, SFR, and EW pipelines; master/geometry catalogues; normal and 7000 product locations; QC figures/reports; model experiments; and older research exports. Their roles matter more than the number of saved figures: many notebook copies reuse the same underlying galaxies and are not independent confirmations.

### 2.2 Evidence levels used here

| Label | What it means in this report |
|---|---|
| Live central check | Selected loading and central-estimator cells were run against current FITS products; exact checks are described below. |
| Saved notebook result | Read from the current notebook's stored numerical output and checked against its source definition; the full analysis was not rerun here. |
| Dated-report result | Read from a report and, where relevant, its exported tables; it retains that report's estimator and uncertainty convention. |
| Physical interpretation | A proposed explanation of the measurements, with alternatives and limitations stated. |
| Proposed analysis | Work needed for the paper or the 7000 Angstrom transition; it has not already been completed. |

For the live check, source cells **3, 5, 7, 9, and 11 of N103** were executed in isolation. They load the actual maps, reconstruct the existing masks, and calculate per-product measurements and central stage estimators. A read-only wrapper also records the unchanged SF/ND/NSF counts after the same S/N cut. No bootstrap, spectral fitting, post-processing pipeline, or notebook-output rewrite was run. The fraction partition is exact. The largest difference from the saved four-decimal category fractions is **0.00004043**; the largest difference from the saved five-decimal log-decomposition values is **0.0000050 dex**. These are rounding-level agreements, not independent confirmation of the calibration.

The PDF and report figures are generated from these checked central tables and explicitly identified saved results. No source notebook or FITS product was changed. The only stored notebook exception in the inventory is a `SyntaxError` in N055, the July 16 NGC4254 spiral experiment; that notebook is not a basis for the main paper. N106 also contains an explicit warning that a plotting-layout edit had not been verified when its accompanying prose was written. Neither an old warning nor later execution counts alone settle source/output consistency.

### 2.3 Relationship to Paper I

Paper I studied the resolved mass-metallicity relation and its SFR dependence in 14 Virgo spirals, including the mass-dependent sign change. Paper II can be a distinct observational question using the expanded workspace: **how does the spatial distribution of SF, its conditional intensity, and the observable non-SF-selected emission vary with environmental stage?** The metallicity calibrations, rMZR, and secondary abundance relations should not form a parallel results thread here. The published first-paper identity and scope were checked through its [author preprint](https://arxiv.org/abs/2605.31412); the supplied ADS page was unavailable to the web tool.

## 3. Definitions that must remain fixed through the story

### 3.1 Sample and independent units

The live primary analysis contains **26 product groups representing 27 physical galaxies**. NGC4567 and NGC4568 share `NGC4567_8`; that product is currently one unit in weighting and resampling. There are 8 pre-peak, 5 close-to-peak, and 13 post-peak product groups. The catalogue contains 40 physical galaxies, while the local normal product tree has 35 matching mass and gas further products. Those three counts describe different populations and must not be used interchangeably.

After the inclination selection, 29 catalogue product groups are candidates. NGC4548, NGC4579, and NGC4689 lack the required products in the current local tree and are excluded; they are **not non-detections**. The 7000 product directory currently has no FITS/FITS.gz files, although separate QC documents exist elsewhere. This is a local availability statement, not a statement about remote production progress.

Environmental stages come from the catalogue `class` mapping: 1 pre-peak, 2 close-to-peak, 3 post-peak. They are a cross-sectional ordering of different galaxies. They are not elapsed times, and the current calculations do not measure how a particular observed spaxel moves between stages. Pre-peak is an internal, relatively less-processed reference; it is not automatically an unperturbed field control. The environmental-history reports and pressure experiments further show why interactions, projection, and heterogeneous histories must remain visible [R15; N067-N072].

### 3.2 The usable footprint and the three observational categories

The mature September/October analyses generally use inclination **i < 80 degrees**, finite stellar surface density and gas-bin geometry, the nGIST spatial mask equal to zero, finite deprojected radius, and **SNR_POSTFIT > 25**. The initial valid footprint is constructed before line-detection requirements, so the denominator can contain genuine line non-detections. Mass/gas WCS agreement is checked, and the combined NGC4567/8 product uses the minimum elliptical R/Re of its components. The notebooks apply no additional inclination correction to the supplied surface-density maps.

**Selected SF** requires a finite `LOGSFR_SURFACE_DENSITY_HII` value, finite proxy EW(Halpha) > 6 Angstrom, and finite intrinsic Halpha dispersion below 45 km/s. The intrinsic width is sqrt(observed width squared minus correction squared); non-positive differences are unavailable rather than valid zero-width SF. Thus the implemented selection effectively requires a positive resolved width as well as the upper cut.

**ND** is usable area without both Halpha and Hbeta passing the product's raw detection criteria. **NSF** is usable area with both Balmer lines detected, but outside the selected SF mask. For the inspected representative product the thresholds are S/N >= 3 and flux >= 20 in the pipeline's stored flux units; the loader reads those values from each product. The same header records Calzetti dust correction, the Chabrier SFR coefficient 4.98358 x 10^-42, and `BPTMODE=both`.

| Category | Operational meaning | Interpretation to avoid |
|---|---|---|
| SF | HII-map, EW, and intrinsic-width selection within usable area | All true star formation, or a complete census of individual HII regions |
| ND | Failure of joint Balmer detection | Zero gas, zero SFR, or a formal SFR upper limit from the stored ordinate |
| NSF | Balmer-detected but fails one or more SF gates | A pure DIG, shock, AGN, or old-star component |
| TNSF | NSF and NII_BPT=3 and SII_BPT in {2,3} | A physically certified absence of young-star ionization |

SF + ND + NSF exactly partition the retained footprint. The NSF reason breakdown uses the first failed gate in a prescribed order; its components are bookkeeping categories, not mutually exclusive physical causes. In particular, low EW, unavailable measurements, ambiguous BPT classification, and high width can all remove a location from SF. The broad pipeline `_SF` maps used in some early notebooks include a different classification branch from the `_HII` plus EW/width selection used in the mature paper analysis. A suffix alone does not establish comparable selections.

The pipeline's all-spaxel Halpha-equivalent SFR values in ND/NSF must not be promoted to true-SFR measurements or upper/lower confidence bounds. The August 31 category notebook explicitly warns about substituted/floored ND values. September 11 corrected the earlier language of upper/lower SFR boundaries [N092; R06; R09].

### 3.3 TNSF and the NSF complement

The current dual-BPT mask is stricter than the earlier shorthand 'remove Composite'. NII code 3 is the AGN-side class; SII codes 2 and 3 are LINER and Seyfert. Consequently, **NSF minus TNSF is not identical to Composite**: it also includes other detected regions failing the joint selection. The October 8 report's shorthand must be updated in any manuscript. An EW(Halpha) < 3 Angstrom subsample is a useful low-EW diagnostic, but it does not by itself identify the ionizing source in every resolved element. The [WHAN classification paper](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1365-2966.2011.18244.x) motivates the diagnostic context, not a guarantee of component purity.

### 3.4 Four different averages answer four different questions

For a product g and a chosen spatial/density bin, define f_g = N_SF/N_usable, I_g as the **linear** mean Sigma_SFR in selected SF pixels, and T_g as their summed Sigma_SFR divided by N_usable. Then T_g = f_g I_g. A product with no selected SF has T_g = 0 and an undefined conditional intensity.

The N103 stage estimator is O = mean_g(f_g), T = mean_g(T_g), and **I_eff = T/O** when O > 0. Therefore T = O I_eff exactly and the positive pre/stage ratios obey:

> Delta log T = Delta log O + Delta log I_eff.

I_eff is an **occupancy-weighted mean of the product-level conditional intensities**. It is not simply mean_g(I_g), a median of spaxel log SFR, or the Q3 median residual. This distinction is essential when combining the August, September, and model-fit numbers. Averaging f and I separately and multiplying them generally breaks the identity.

Equal-product category fractions describe a typical contributing product. Pooled-pixel fractions describe the observed pixel population and can be dominated by large footprints. Pixels also share fitted bins and the PSF, so neither millions of pixels nor 25 distinct gas bins are automatically that many independent physical observations. At differing distances, pixel pooling is not exactly physical-area weighting across galaxies.

The usual comparison windows are log10(Sigma_star/[M_sun kpc^-2]) = 7.0-9.5 and R/Re = 0-2.5, using 0.25-wide bins. They are **display/common-range choices, not guarantees of adequate galaxy support in every bin**. The saved SFR table contains a close-to-peak 9.25-9.5 bin with only one usable product. Some bins have finite central intensities even when few products meet the separate 20-SF-pixel support flag. The final paper must enforce, not merely print, its support requirements.

## 4. How the discoveries developed

| Period and notebook family | What the work added | Role in Paper II |
|---|---|---|
| March-June: N001-N052 | Reduction, tied-line, wavelength, mass, EW, and product comparisons; mixed early science | Calibration provenance and justification for the new reduction |
| July: N053-N057 | NGC4254 spiral experiments; first infall-coloured rSFMS/rMZR view | Stage-rSFMS starting point; spiral branch remains separate |
| Early/mid August: N058-N078 | Individual-stage relations, pooled versus equal weighting, density distributions, Re restrictions, pressure criteria | Identifies sample weighting, selection, and spatial support as central issues |
| Late August: N085-N090 | SF occupancy, surviving-SF distribution, and first spatial-asymmetry tests | Establishes the two-component SFR question and tests enhancement/broadening |
| August 31-September 9: N091-N105 | SF/ND/NSF bookkeeping, continuum-S/N control, galaxy bootstrap, exact SFR decomposition, absolute lines | Main empirical backbone |
| September 14-October 4: N106, N108-N110 | Directional candidate, line fading, HII invariance/mixing tests, gas/HOLMES models | Mechanism tests and limits; not one validated evolutionary model |
| October 6-9: N111-N125 | Stronger dispersion support, dual-BPT TNSF, signed-distance PDFs, conditional excitation tests | Refines what can be said about residual emission |
| October 10: N126 | NGC4654 spatially clustered inference and held-out external comparison | Separates detected asymmetry from unproven absolute enhancement |

The dates do not represent 128 independent findings. The scientific progression is from a survivor-only relation, to selection-aware population accounting, to tests of what powers the remaining emission. Several negative or qualifying results are useful: they prevent the narrative from resting on selection-induced enhancement, broadening that is not detected, or an over-constrained ionization model.

## 5. The observational findings, with their limits

### 5.1 The clearest signal is loss of selected SF occupancy

The live full-footprint equal-product fractions are:

| Stage | Product groups | SF | ND | NSF |
|---|---:|---:|---:|---:|
| Pre-peak | 8 | 0.494 | 0.137 | 0.369 |
| Close-to-peak | 5 | 0.332 | 0.386 | 0.282 |
| Post-peak | 13 | 0.101 | 0.628 | 0.272 |

These are means of each product's fraction after SNR_POSTFIT > 25, over its entire usable footprint. They are not mass-matched stage effects and have no newly computed uncertainty intervals here. They nevertheless make an important wording correction unavoidable: **the global mean NSF footprint fraction does not increase monotonically**. The increase discussed in the notes concerns particular central/high-density locations, and the NSF share of the remaining detected population can rise even while its share of the whole footprint falls. Ratios formed from these mean fractions would define another estimator and should be labelled accordingly.

The resolved profiles are more informative than the footprint averages. The September 14 report gives SF fraction 0.429 -> 0.027 and ND fraction 0.317 -> 0.789 at R/Re = 1.625. At R/Re = 0.125, NSF fraction changes from 0.472 to 0.694 [R10, section 1.3]. The mature N104 plots show post-peak SF fractions of order 0.1-0.2 in much of the inner disc and near-zero fractions beyond about 1.5 Re, with changing support. These are representative profile features, not universal physical thresholds or measurements of each galaxy's truncation radius.

The earlier N085 analysis already checked the full footprint, R <= Re, inclination restrictions, a common global-mass subset, HII versus broader SF classification, and individual-product diversity. Its saved output shows the post-peak occupancy deficit survives those checks. At the same time, these checks used an earlier analysis definition and cannot replace rerunning the final common selection. In particular, an effect inside Re means it is not solely an outer-edge phenomenon; it does not establish a causal local stellar-density dependence independent of radius.

**Paper-ready direction:** lead with the observed spatial loss of the selected SF population. Call f_SF an occupancy or selected-SF area fraction. 'Survival' can be a useful shorthand, provided the text does not imply tracked regions or a directly measured temporal survival probability.

![Figure 1. Current central stage accounting, reconstructed from live products. Panel A shows the equal-product mean full-footprint category fraction; panels B and C show selected examples of the post/pre SFR decomposition. These are central estimates, not new significance tests.](assets/20261011_paper_II_synthesis/figure_01_central_findings.png)

### 5.2 Survivors shift downward; there is no supported close-to-peak broadening signal

N087 asks whether the surviving SF distribution is shifted, broadened, or develops stronger enhanced/suppressed tails relative to an internal pre-peak relation. It uses product-level statistics and label permutations. The following are its **saved August 29 values**, before the later strict-continuum-S/N framework; they should be treated as a useful independent-estimator check, not numerically merged with N103.

| Statistic | Pre-peak | Close-to-peak | Post-peak |
|---|---:|---:|---:|
| Median product residual M, dex | -0.001 | -0.138 | -0.243 |
| Median product width W, dex | 0.395 | 0.384 | 0.375 |
| Median enhanced fraction | 0.175 | 0.102 | 0.055 |
| Median suppressed fraction | 0.097 | 0.149 | 0.238 |

The notebook defines its tails using the pre-peak scatter-normalized residual; they are not percentages of the entire usable disc. The post/pre median-residual difference is -0.241 dex with an unadjusted two-sided galaxy-label permutation p = 0.0372. The enhanced-fraction difference is -0.120 with p = 0.0020. The post/pre width change is -0.019 dex with p = 0.559, and the suppressed-fraction difference has p = 0.135. Multiple metrics were examined, so the individual p-values should not be presented as a set of independent discoveries. The close-to-peak sample has only five products and none of its main changes is secure in that permutation comparison.

This rules against building the paper around a claimed **broadening or enhancement-rich phase near peak** on the strength of these results. The positive result is that the post-peak survivor distribution loses its enhanced tail and shifts toward lower residuals, without evidence for increased robust width or total tail fraction. The low absolute p-value for the enhanced-fraction change is promising, but the exact estimator, reference construction, and multiplicity choices should be frozen for the 7000 rerun.

The September 14 fit [R10, section 6] addresses a different statistic: the linear-intensity decomposition profile. It reports attenuation -0.141 [-0.248,-0.059] dex for close-to-peak and -0.390 [-0.468,-0.298] dex for post-peak, with **16th-84th percentile** whole-galaxy bootstrap intervals. The latter corresponds to a factor 0.407 [0.341,0.504]. Post-peak leave-one-galaxy estimates range from -0.420 to -0.340 dex. The density-slope intervals include zero, and the bins are highly correlated. These are descriptive fits, not 95% detection statements or a model-selection result.

There is no contradiction between N085's weak per-bin survivor evidence, N087's residual-distribution shift, and the later attenuation fit. They differ in selection, support, statistic, and how information is combined across density bins. The report should preserve that progression rather than declaring the earliest cautious conclusion wrong.

### 5.3 The SFR deficit changes its composition across the disc

The strongest bridge between the rSFMS and footprint results is N103's exact decomposition. Representative **live-reproduced central post/pre contrasts** are:

| Coordinate/bin centre | Delta log occupancy | Delta log survivor intensity | Delta log HII-traced total |
|---|---:|---:|---:|
| log Sigma_star = 7.625 | -1.011 | -0.429 | -1.440 |
| log Sigma_star = 8.125 | -0.803 | -0.336 | -1.139 |
| log Sigma_star = 8.625 | -0.204 | -0.408 | -0.613 |
| R/Re = 0.625 | -0.384 | -0.438 | -0.821 |
| R/Re = 1.625 | -1.204 | +0.105 | -1.099 |

At log Sigma_star = 7.625, most of the logarithmic decline is in occupancy, with additional dimming among selected SF. At 8.625, the conditional-intensity term is larger than the occupancy term. At 1.625 Re, the loss of area dominates while the few survivors do not have a negative central intensity offset. Thus **'disappearance dominates everywhere' is too strong**, even though area loss is the clearest broad sample transformation.

The September 14 resampling of both stages gives the outer survivor offset +0.105 [-0.213,+0.345] dex and the 8.625 intensity offset -0.408 [-0.529,-0.263] dex, again 16th-84th percentiles. The original N103 comparison holds the pre-peak reference fixed in its bootstrap; its plotted bands and the later report's two-stage intervals are not the same uncertainty calculation. The final paper should resample the reference and comparison consistently.

The accounting refers to **selected HII-traced SFR over the observed usable area**. It is not total multiwavelength SFR, a reconstruction of hidden SF, or a decomposition of a galaxy's true quenching rate. The main advance is observational: it measures what the survivor-only rSFMS omits.

### 5.4 Inner detected NSF becomes relatively prominent without universal brightening

N105 and the September reports compare NSF occupancy with its absolute corrected Halpha intensity. This was the correct next question: a larger NSF area alone cannot tell us whether a luminous new excitation component was created or whether previously subdominant emission became more visible as young-star emission faded.

The saved N105 summary gives median post/pre NSF Halpha-intensity contrasts of -0.674 dex across its fixed stellar-density bins and -0.545 dex across its radial bins. These are medians of bin contrasts, not population fits or independent-bin averages with calibrated errors. The corresponding median NSF-fraction contrasts are only +0.018 and -0.031 dex, hiding substantial spatial structure. At density 8.625, the September 14 report gives Delta log f_NSF = +0.080 [-0.039,+0.223] and Delta log I_NSF = -0.508 [-0.721,-0.150]. The fraction increase in that bin is uncertain; central radial bins make the growing NSF occupancy clearer.

The defensible conclusion is that **NSF's increased local prominence is not a uniform increase in its absolute Halpha brightness**. A declining class-average intensity does not rule out new or brighter shocks/AGN/HOLMES-associated emission in a subset, because membership changes and multiple sources can contribute. The earlier statement that this 'rules out creation' is stronger than the data warrant [R08 versus R09-R10].

For NSF contribution to usable-area luminosity, retain the product-level f_g I_g accounting. A product of separately averaged NSF fractions and conditional intensities is generally not the mean contribution. This is the same covariance issue that motivated the exact SF decomposition.

### 5.5 Excitation changes are real descriptive results, but depend on line and weighting

The absolute-line work asks whether ratios increase because the numerator brightens or because the denominator fades more. At density 8.625, the September common-detected-line export has NSF mean Halpha surface luminosities 3.558, 2.077, and 1.111 x 10^39 erg s^-1 kpc^-2 for pre/close/post, with 6/5/13 contributing products. The ratios of equal-product linear line means are:

| Stage | N2 = [NII]/Halpha | S2 = [SII]sum/Halpha | O3 = [OIII]/Hbeta |
|---|---:|---:|---:|
| Pre-peak | 0.3731 | 0.3773 | 0.8978 |
| Close-to-peak | 0.4367 | 0.3857 | 0.4678 |
| Post-peak | 0.5663 | 0.4168 | 0.7528 |

N2 and S2 rise in this example, while O3 is neither monotonically rising nor larger post-peak than pre-peak. At a different radial coordinate, all three contrasts can be positive. Therefore 'all forbidden lines fade more slowly than Balmer lines' and 'all BPT ratios increase along infall' should be removed from the central claim. The [OIII] numerator and actual Hbeta denominator must be retained separately: clipping negative extinction does not force every corrected Balmer decrement to 2.86 [R10-R12].

The October BPT-distance notebooks add a compact distributional view. They calculate the signed global shortest distance to the relevant KE01 curve, with equal metric weight on the two log-ratio axes. This is a geometric coordinate in dex, **not a physical mixing fraction or a statistical distance including flux errors**. The classification curves originate in [Kewley et al. (2001)](https://arxiv.org/abs/astro-ph/0106324).

| NSF statistic | Pre-peak | Close-to-peak | Post-peak |
|---|---:|---:|---:|
| NII signed-distance median, pooled pixels | -0.282 | -0.267 | +0.027 |
| NII signed-distance median, equal-product PDF | -0.347 | -0.277 | -0.053 |
| SII signed-distance median, pooled pixels | -0.123 | -0.050 | +0.082 |
| SII signed-distance median, equal-product PDF | -0.051 | -0.022 | +0.003 |

These saved N120/N121 summaries support a post-peak shift in NSF distributions, especially in NII, but show that the sign of the median and the strength of the shift depend on weighting. A pooled NII median crossing zero is not the same claim as a typical galaxy crossing a boundary. SF remains on the negative side; much of the SF/NSF difference is partly built into their definitions and should not be counted as independent evidence for environmental excitation changes.

N120's additive NSF partition finds TNSF fractions of approximately 0.96%, 1.28%, and 11.69% of the finite NII-distance NSF sample. The SII counterpart is very similar but has its own denominator. These are **pooled finite-distance fractions**, not the fraction of total luminosity or the fraction of galaxies hosting a new source. TNSF and NSF-TNSF share the NSF normalization, so their densities sum to the NSF density. In N121 the catalogue PDFs are individually normalized; equal-product TNSF PDFs also condition on the subset of products containing TNSF. Those two plot families cannot be interpreted using the same area-under-curve logic.

### 5.6 The October controls strengthen the diagnostic question, not a unique mechanism

N125 repeats the additive distance distributions in separate stellar-density and radial intervals. These are separate one-dimensional controls, not joint matching in both coordinates. The plotted 8.5-9.0 density interval still shows the post-peak NSF shift. A final population claim needs per-product effect sizes and common joint support; a binned visual persistence alone does not remove all structural confounding.

N124 goes further by comparing NSF ratios at fixed corrected Halpha surface-brightness intervals. In the saved 39.0 <= log Sigma_Halpha < 39.5 interval, all three stages contribute 8/5/13 products. Post/pre pooled median shifts are approximately +0.121 dex in N2, +0.081 in S2, and +0.174 in O3. The notebook's own proposed direction, a possible higher pre-peak ratio at fixed brightness, is therefore not what this example shows. The result is descriptive and conditioning is only within a finite 0.5-dex brightness bin. Positive finite forbidden fluxes are not a common robust forbidden-line detection gate; shared denominators, extinction, within-bin brightness distributions, and galaxy dominance still matter. This test is promising for mechanism discrimination but is not already an sSFR-controlled result.

N123's high-intrinsic-width (>45 km/s) catalogue overlaps the NSF/TNSF populations; it is not a fourth mutually exclusive component of the footprint. Its distributions show why high width and high line ratio should be compared directly rather than assumed interchangeable. The SF width cut also means the selected SF catalogue cannot test a high-width SF tail.

The dedicated N111 median-dispersion analysis is the strongest warning against a simple turbulence narrative. It requires at least **25 distinct fitted gas bins per supporting product/spatial bin** and at least **5 NSF-supporting products per stage bin**, with a shared parent-NSF support convention for the restricted catalogues. It uses medians of within-product medians for the equal-product branch. Its saved plots and the October 6 note do not establish a clean monotonic stage trend. Restricted TNSF/EW subsets can have fewer effective contributors than the parent support count, and linewidths can also reflect unresolved rotation, fitting limitations, and source mixtures.

Earlier line-profile notebooks using other estimators/support do show some positive post/pre width contrasts. The appropriate synthesis is **estimator- and selection-sensitive linewidth evidence**, not either 'there is universally more turbulence' or 'all width differences vanish'. N115's [NII] width variant is another consistency check. High ratios do not independently prove shocks, and an EW < 3 cut does not independently prove HOLMES.

### 5.7 NGC4654: spatial asymmetry is detected; absolute enhancement remains unestablished

The August/September direction scans and plane fits identified a northwest-facing candidate near a coarse VIVA compression direction. The September 14 result had PA about 319 degrees versus the adopted 315-degree compass direction, but a small explained variance (R-squared about 0.044), no formal PA uncertainty, and no correlation-aware significance. Its preferred side was only +0.035 dex above the internal pre-peak reference, while the opposite side was about -0.096 dex. A side contrast is not the same thing as enhancement above a control [N106; R10].

N126 is a material advance. It defines coarse cells and gas-bin-aware spatial clusters, fits an m=1 term with radial and m=2 controls, compares pre-peak systems, then separates discovery and buffered confirmation regions. Its current saved result is:

| Measurement | Saved result | Meaning |
|---|---|---|
| EW coherent m=1 amplitude | 0.113 dex; Holm p = 0.004 | Conditional asymmetry detection at the adopted scale |
| SFR coherent m=1 amplitude | 0.163 dex; Holm p = 0.004 | Conditional asymmetry detection at the adopted scale |
| Held-out positive-side residual | +0.100 dex; 95% CI [-0.120,+0.249] | Absolute enhancement not established |
| Held-out positive-minus-negative | +0.103 dex; 95% CI [-0.020,+0.223] | External-residual hemisphere contrast also uncertain |
| One-sided absolute-enhancement test | p = 0.223 | Does not support a detection under this test |

![Figure 2. Selected saved diagnostics summarized without rerunning their inference. Panel A: August Q3 median survivor residuals with its stored 95% bootstrap intervals. Panel B: October 10 NGC4654 external-reference estimates with stored 95% intervals. The analyses use different samples and estimators and are not combined into one test.](assets/20261011_paper_II_synthesis/figure_02_survivors_and_asymmetry.png)

The Fourier fits use 1,803 coarse cells and 249 clusters, not native-spaxel sample sizes. Doubling the cell/block scales leaves similar amplitudes. The stellar mass map also has a smaller but detected asymmetry, which prevents an exclusive attribution to environmental gas physics. The external-control transfer diagnostics show substantial heterogeneity and some controls without supported predictions; the empirical six-control exceedance ranks are 0.571 for EW and 0.143 for SFR, not matched-population significance tests. The target and historical direction were already known, so this remains an exploratory retrospective analysis despite its data separation.

**Recommended placement:** a qualified case study or appendix, after the population result. Its most useful message is that directional structure can coexist with population suppression. It should not anchor a claim of triggered star formation. VIVA supplies contextual HI morphology; a coarse compass direction without a measured uncertainty is not a wind vector or an alignment probability. The [VIVA atlas](https://arxiv.org/abs/0909.0781) is the primary source for that context.

### 5.8 What the pressure and spiral branches contribute

N067-N072 compare different expected-HI, optical-size, initial-HI-size, and rotation-speed prescriptions in a Koppen-style pressure calculation. They help assess environmental histories and systematic dependence of the adopted stage context. Their 40-physical-galaxy tables are not the optical 26-product sample. Some models leave pressure undefined for inclination-zero catalogue entries; common scatter and projected geometry also limit population inference. Treat them as context or stage-sensitivity support, not a new measured elapsed-time axis. Do not replace the optical stage labels with whichever model yields the clearest ordering.

The NGC4254 KTZ/spiral experiments are substantial but answer another question. Some fits condition on winding, compare different arm geometries, or explore harmonic orders. N053 explicitly declines to claim an arm-number winner or a fit to the full diffusion/enrichment model. They should remain outside the main Paper II narrative unless a tightly bounded structural-control use becomes necessary. Their existence does not require a spiral-arm section in this paper.

## 6. What physical picture is consistent, and where it fails

### 6.1 A coherent interpretation at the level the data support

A consistent sequence is: environmental processing changes the gas supply/distribution; fewer locations satisfy the detectable HII-region selection; surviving selected regions in the common density domain tend to have lower Halpha-traced intensity; low-density outer emission increasingly enters ND, while retained central emission can fail the SF gates and display a different mix of excitation. A decreasing young-star contribution can make other ionization channels more prominent without making every forbidden line brighter.

This explains why the results fit together without requiring each stage to be intermediate in every observable. Different gas phases, local structure, source mixtures, and galaxy histories can respond differently. The measured relationship is stage-associated and cross-sectional. The [ram-pressure review by Boselli, Fossati and Sun](https://arxiv.org/abs/2109.13614) supplies the broader physical context for multiphase gas removal and SF changes; it does not validate a unique mechanism for this particular sample.

The empirical sequence is more secure than its causal arrows. The footprint categories are outputs of detection and classification, not directly measured gas phases. A region can leave SF by fading below a limit, changing EW or ratios, having an unresolved width, or losing gas. The paper should show how much each observable changes before interpreting these paths.

### 6.2 Three useful equations, with separate evidential status

**Observed accounting:** T = O I_eff. This is exact for the implemented positive-occupancy estimator and is a strong organizing principle.

**Molecular-reservoir interpretation:** at fixed conversion between molecular gas and SFR, Sigma_SFR = Sigma_H2 / tau_dep. A lower SFR can reflect less molecular gas, a longer depletion time, or both. The optical result alone does not choose between gas amount and efficiency. A later CO comparison would need common resolution and support; the [VERTICO survey](https://arxiv.org/abs/2111.00937) provides the relevant observational resource, but no new resolved gas-consumption measurement is made here.

**Line-mixture interpretation:** for two components with line-to-Balmer ratios R_y and R_h, R_obs = (1-w) R_y + w R_h, where w is the appropriate Balmer-luminosity fraction. Ratios must mix in linear space. A common Halpha/Hbeta mixing weight requires compatible decrements and attenuation assumptions. If all components tied to young stars fade in the same proportion, their relative weights do not change; 'less SF' alone is not a complete explanation for changing ratios.

The PHANGS analysis [A tale of two DIGs](https://doi.org/10.1051/0004-6361/202141859) provides primary evidence that leaked young-star photons and hot low-mass evolved stars can contribute differently to diffuse emission. It motivates a multi-source model, but it does not make the MAUVE operational NSF or TNSF labels physical component measurements.

### 6.3 The model sequence is productive, but it has not solved the data

The September model reports evolve from a general gas-loss plus line-response framework toward stricter HI-only stripping and young-plus-HOLMES models. That progression matters: adding an explicit old-star photon budget exposes constraints that flexible differential fading can hide.

The September 30 HII-invariance/common-weight notebook N109 directly tests prerequisites. In its strict HII-like selection E50_L39.5_S35, only three matched stellar-density/radius bins remain. The post/pre O3 contrast is -0.437 dex with a saved interval [-0.573,-0.290], labelled 'offset beyond tolerance'; N2 and S2 do not jointly establish the required practical invariance. The empirical two-endpoint weights move toward the diffuse-like endpoint post-peak, but depend strongly on the chosen template selection and often approach the boundary. Those fitted weights are conditional diagnostics, not calibrated physical DIG fractions. Held-out-ratio residuals and endpoint sensitivity prevent treating a single changing weight as a demonstrated explanation.

The October 2 Model 0 report [R12] exposes a second issue. At the density-8.625 example, its generous specified local old-star normalization gives Halpha_HOLMES = 4.050 x 10^37 erg s^-1 kpc^-2, compared with total NSF Halpha 3.558 x 10^39 pre-peak and 1.111 x 10^39 post-peak. The corresponding adopted-model weights are only 1.14% and 3.65%. With the chosen fixed spectra, the N2 increase is too small and the O3 evolution has the wrong direction. The inverse fixed-spectrum exercise requires a negative O3/Hbeta component ratio (-4.82), which is unphysical under those assumptions. This is a conditional central-value inconsistency, not a statistical exclusion of all HOLMES models.

The same report's no-feedback, fixed-depletion example needs roughly four Gyr to reach the specified fading factor. Its comparison with one Gyr is explicitly hypothetical; the stages have not supplied a measured clock. This is a test of assumed rates, not a measured quenching time or proof that molecular gas must be directly stripped.

N110 then uses **different illustrative assumptions**: feedback loading eta = 3, an imposed initial 30% HOLMES Balmer fraction, and user-selected young/old spectra. It prints the difference between its model normalization and observed comparison rather than fitting them. Its saved two-Gyr N2/S2/O3 changes are +0.415/+0.429/+0.380 dex, compared with its observed NSF median contrasts +0.042/+0.025/-0.071 in the chosen interval. There is no assigned two-Gyr separation between the observed stages. Consequently, those curves demonstrate possible temporal and directional responses; they do not repair the quantitative comparison by themselves.

### 6.4 The appropriate physical conclusion

The report should retain **fading young-star emission with a changing relative contribution from other sources** as a testable interpretation. It should not conclude that old stars uniquely explain the BPT shift, that shocks are required by linewidths, or that a constant local two-component model already fits all observables. Low-EW spectra, AGN/bulge contamination, transport of ionizing photons, changing nebular conditions, detection censoring, and population mixing are still relevant alternatives or modifiers. Abundance-related degeneracies can be acknowledged as limits without opening a metallicity results programme in Paper II.

A successful eventual physical model must reproduce **absolute line intensities, ratios, category occupancy, and their spatial support together**. Matching a BPT track alone is insufficient. It should generate observables before applying the same line-noise, dust, EW, width, and classification procedures. This is a longer-term interpretation goal; a complete radiation-transfer model is not necessary for a strong empirical Paper II if the current physical claims remain appropriately bounded.

## 7. A single paper narrative, with clear priorities

### 7.1 The main research question

**At fixed local stellar mass surface density and with explicit control of usable spatial coverage, how do the occupancy and conditional intensity of selected star formation change across Virgo environmental stages, and what happens to the ionized emission outside that population?**

This gives one primary question with a diagnostic extension. It avoids making the paper simultaneously about rSFMS slopes, spiral modes, complete gas-regulator histories, metallicity, shocks, and local enhancement.

### 7.2 The argument in reading order

First show the survivor-only rSFMS and explain why it cannot characterize the full observed disc. Then reveal where SF disappears from the selection, separating ND from detected NSF. Next quantify occupancy and intensity in one exact decomposition, so the reader sees what the conventional relation misses. Finally show the absolute lines and excitation of the residual detected population and assess whether fading/relative-source-change is plausible. Local asymmetry and analytical models then become tests and qualifications of this picture, rather than separate stories competing for the title.

The conclusion can be stronger scientifically by preserving two non-results: there is no current evidence for a general enhancement-rich close-to-peak distribution, and the simplest universal excitation explanation does not satisfy all current tests. These constrain the physical story instead of weakening the central empirical result.

### 7.3 Claim hierarchy

| Priority | Claim | Current assessment |
|---|---|---|
| Primary | Post-peak selected-SF occupancy is strongly reduced and spatially structured | Best supported; central values reproduced live; final 7000 validation needed |
| Primary | Surviving SF intensity also changes at fixed density | Supported by several estimators, with different selections/uncertainties; not uniform at all radii |
| Primary | An rSFMS alone misses the occupancy term | Exact consequence of the measured selection and decomposition |
| Supporting | Central NSF prominence can rise while conditional Halpha intensity falls | Supported in specified bins; not a global monotonic NSF increase |
| Supporting | NSF excitation distributions shift post-peak | Strong descriptive evidence; quantify galaxy-level uncertainty and selection effects |
| Interpretation | Fading young sources increases the relative importance of other excitation | Plausible; simple fixed-spectrum/budget models are insufficient |
| Optional | NGC4654 has directional SF structure | Conditional asymmetry detected; absolute enhancement and RPS causation unestablished |
| Do not claim | Stages yield a quenching clock, all NSF is HOLMES, or all ratios/widths rise | Not established |

### 7.4 What makes this distinct and worthwhile

The contribution is a spatially resolved, selection-aware decomposition in cluster galaxies: **the environment-associated loss of SF cannot be described completely by the locus of its surviving regions**. The residual-emission analysis then tests whether changing class fractions imply newly luminous gas or changing visibility of different excitation sources. This connects mass-SFR parameter space to where activity remains in a disc.

The [Pessa et al. completeness and spatial-scale study](https://arxiv.org/abs/2104.09536) is useful framing: excluding non-detections can alter resolved scaling relations. The MAUVE advance should be stated in terms of its cluster-stage comparison and explicit occupancy/intensity accounting. It should not claim that previous IFU work ignored environments or that all pre-peak Virgo systems are equivalent to the PHANGS field sample. The August 3 report already showed that changing to the comparison SF selection weakens the initial impression of universal enhancement above PHANGS [R02].

## 8. Recommended manuscript structure

The following is a concrete structure to draft once the 7000 products and final selections are frozen. Approximate lengths are editorial targets, not journal limits. A main text of roughly 7,000-9,000 words plus appendices is reasonable for this scope.

### 1. Introduction: why a surviving-SF relation is incomplete

About 800-1,000 words. Introduce cluster environmental processing, spatially selective SF changes, and the distinction between an active-region relation and a whole observed-disc census. State the primary question and the reason for using both stellar-density and radial coordinates. Briefly distinguish this project from Paper I; reserve metallicity for future work. End with the three empirical objectives: occupancy, survivor intensity, residual emission.

### 2. Sample, spectral products, and measurement definitions

About 1,300-1,700 words. Describe the physical-galaxy sample, product groups, environmental-stage provenance and uncertainties, 7000 fitting choice, spatial resolution/coverage, mass and Halpha calibration, geometry, usable footprint, masks, and classification. Explicitly distinguish the broader pipeline SF branch from the paper's selected HII population. Give per-stage sample attrition and the treatment of NGC4567/8. State the independent sampling unit and show continuum/line support. This is the place to explain what ND and NSF mean.

### 3. Where star formation remains

About 1,200-1,500 words. Start with a compact rSFMS comparison to orient the reader, then the SF/ND/NSF maps and equal-product occupancy profiles versus density and radius. Give per-galaxy scatter and supported stage contrasts. Show whether the pattern persists within Re and on overlapping global-mass/radius support. End with the need to measure intensity among survivors separately.

### 4. How strongly the surviving population forms stars

About 1,200-1,500 words. Present the exact occupancy-intensity-total decomposition, then the conditional residual distribution and tail balance. Emphasize the density dependence of the decomposition and the different radial result for outer survivors. Keep mean-linear and median-log results distinct. A constant attenuation fit can be a concise summary if it survives the revised products; do not translate it directly into a time.

### 5. The ionized emission outside the selected SF population

About 1,200-1,600 words. Present absolute NSF Halpha/line intensities alongside occupancy, then NII/SII/OIII ratios and BPT-distance distributions. Use common line support, equal-product summaries, and a transparent NSF/TNSF partition. Include one high-value control, such as joint structural matching or fixed-Halpha comparison. Keep low-EW and linewidth checks concise; they constrain interpretations but do not label a unique source.

### 6. Discussion: spatially selective suppression and the changing emission mix

About 1,000-1,400 words. Connect the two SF terms to possible gas loss/supply/efficiency changes without claiming the optical data separate them. Explain how ratios can rise during absolute fading. Present the successes and failures of the simple source-mixture picture. NGC4654 can be a short case study showing why asymmetry differs from triggered enhancement. Discuss stage heterogeneity, selection, resolution, and the 7000 comparison. Resolved CO or more detailed source modelling belongs to the next discriminating tests, not to an unsupported conclusion.

### 7. Conclusions

About 350-500 words. Four or five numbered findings: selected area loss; surviving-intensity change; their relative contributions by location; line-specific residual-emission evolution; limits on enhancement and unique excitation/clock inference. Include only claims that have passed the final common analysis. Use one clearly labelled forward-looking sentence for gas/ionization follow-up.

### Proposed main figures and tables

| Item | Scientific job | Existing starting point |
|---|---|---|
| Figure 1 | Sample/stage overview, coverage, representative category maps | N104 and geometry/master tables |
| Figure 2 | Survivor-only rSFMS and the selection problem | N058-N078, N092-N096; final mask rerun |
| Figure 3 | SF/ND/NSF occupancy versus density and radius, with support | N104 |
| Figure 4 | Exact occupancy/intensity/total decomposition | N103; report Figure 1 |
| Figure 5 | Survivor residual shift and tail balance with galaxy-level summaries | N087, harmonized to final mask |
| Figure 6 | Absolute NSF line/Halpha intensity versus local occupancy | N105, N108, N112 |
| Figure 7 | BPT distributions and the explicit NSF partition on controlled support | N120-N125 |
| Figure 8, optional | Mechanism test or NGC4654 asymmetry versus absolute enhancement | N109/N126; choose one if space is tight |
| Table 1 | Physical sample, product grouping, stage, geometry, coverage, exclusions | Live inventory plus stage references |
| Table 2 | Prespecified effect sizes, intervals, support, and robustness | Final galaxy-level contrasts |

Eight figures are a ceiling, not a requirement. If the main text becomes long, combine Figures 2 and 5, and move NGC4654 to an appendix. Do not let a large gallery obscure the principal measurements. Appendices should hold wavelength/pointing comparisons, per-galaxy maps, mask/detection sensitivity, weighting/support alternatives, linewidth/TNSF details, and conditional model assumptions/failure tests. Pressure-model variants and spiral fits need only a short context reference unless they are used directly.

## 9. The 7000 Angstrom transition is part of the scientific argument

### 9.1 Why the transition matters

The requested plan is to adopt fits ending at 7000 Angstrom because they are expected to reduce pointing bias. That is the user's stated motivation. The local historical comparisons establish wavelength/configuration sensitivity, but this review did not find a completed quantitative demonstration of the pointing-bias improvement across the final sample. The report therefore treats the improvement as a testable expectation.

The main observables survive the wavelength change: Hbeta, [OIII], Halpha, [NII], and both [SII] lines lie below 7000 Angstrom. The broader fit can still affect them through continuum/absorption modelling, tied-line choices, fitted kinematics, stellar population inference, and binning. The migration is therefore more than changing an input directory. Stellar mass, EW, dust corrections, line ratios, and category membership can all move.

The June IC3392 comparison reports an 8900-minus-7000 median pixel age shift of +2.94 Gyr, stellar-width shift -6.97 km/s, and Halpha-velocity shift -2.73 km/s [R16; N050]. NGC4064 has corresponding shifts about +1.06 Gyr, -5.36 km/s, and -2.20 km/s [R17; N051]. These examples are sensitivity evidence, not a sample-wide correction. The May same-bin mass tests also change velocity scale and use a 9100 upper limit; their mass shifts cannot be attributed solely to the proposed 8900-to-7000 switch.

### 9.2 A controlled comparison that can demonstrate reduced pointing bias

**First, preserve provenance.** Freeze the existing 8900 product and notebook manifests. For the 7000 branch record run ID, actual fit interval, templates, polynomial choices, line ties, binning, velocity scale, LSF prescription, masks, and software version. Regenerate mass, SFR, and proxy-EW products consistently. Compare matching configurations where feasible, and label any remaining coupled changes.

**Second, compare matching sky locations.** Use a matched set of galaxies and WCS-aligned sky footprints. Report both a common-footprint comparison and each run's native usable footprint. If binning differs, do not equate BIN_ID numbers between reductions; compare physical sky cells at a common effective resolution or explicitly map the overlap.

**Third, measure pointing structure directly.** Use the known pointing/exposure membership or mosaic boundaries. In overlaps or comparable regions on each side of seams, measure residual offsets in log Sigma_star, continuum/EW, Balmer flux/decrement, line ratios, and intrinsic widths after controlling for the smooth radial/structural trend. Summarize robust per-pointing offsets and the extra between-pointing variance with spatial-block uncertainty. A visual smoother map alone is not sufficient; require the pointing-associated component to decrease without erasing physical gradients or worsening fit residuals. No numerical tolerance is pre-imposed here: choose it relative to the environmental effect sizes and calibration uncertainties before looking at the final stage contrasts.

**Fourth, separate value shifts from selection shifts.** Produce a joint SF/ND/NSF transition matrix on common usable sky pixels and separately count pixels entering/leaving the footprint. Repeat for the EW, width, continuum-S/N, Balmer, and dual-BPT gates. Check whether pointing seams coincide with class boundaries. Report changes by stage, product, and pointing so that a stage-dependent reduction bias cannot cancel in a global average.

**Fifth, rerun the frozen core questions.** Recompute occupancy, survivor intensity, their exact total, and the principal line/BPT comparisons with the same estimator and support policy. Compare effect sizes, galaxy rankings, and leave-one-product results with the 8900 branch. A change in the p-value alone is not the acceptance test. The final science values should come from 7000; the broader-range result becomes a documented sensitivity experiment.

### 9.3 Acceptance conditions before manuscript claims are frozen

| Requirement | Concrete evidence needed |
|---|---|
| Products are internally consistent | Matching runs, readable products/HDUs, WCS/shape checks, correct bin expansion, complete per-product inventory |
| Pointing bias is reduced | Quantified overlap/seam or pointing-offset statistics, with spatial uncertainty and residual-quality checks |
| Category changes are understood | Common-footprint transition tables and first-failed-gate changes, separated from missing area |
| Core environmental effects persist or are revised honestly | Final supported stage effect sizes and uncertainty from the same frozen estimator |
| No false stability from sample changes | Matched-galaxy comparison plus native-sample comparison and explicit attrition |
| Source diagnostics remain interpretable | Corrected direct denominators, line-support audit, width/LSF sensitivity, and final dual-BPT mask |

If the change mainly removes seam structure and preserves the occupancy/intensity findings, it strengthens the paper. If the magnitude changes, replace the numbers and explain why. If one excitation or local-enhancement result weakens, the empirical occupancy/decomposition paper can still stand. If the core stage effect disappears on controlled support, the paper narrative must be revised rather than protected by selecting favourable versions.

## 10. The shortest defensible route to a complete paper

### 10.1 Essential work, in dependency order

1. **Freeze the question, final estimators, and exclusions.** Keep mass/SFR and residual excitation in scope; exclude abundance results. Select a primary galaxy/product weighting, exact support gates, and a small set of headline contrasts before inspecting the 7000 stage outcome. Preserve the current 8/5/13 grouping only as a snapshot, not a hard-coded target.
2. **Complete the 7000 provenance and pointing comparison.** Resolve product availability and regenerate all derived quantities consistently. Quantify the bias improvement described above.
3. **Harmonize the core notebook branches.** N104 and N103 should share one final sample/geometry/footprint. Recompute the Q3 survivor statistic with that same selection. Keep the existing notebooks as provenance and build a clearly dated final analysis entry point or shared loader, rather than silently editing historical results.
4. **Complete population-level uncertainty.** Resample whole product groups in both reference and comparison, carry the reference through residual construction, and assess leave-one-product sensitivity. Choose a small number of primary contrasts or a prespecified joint summary rather than interpreting every pointwise interval as an independent detection. Resolved bins and gas-bin support address coverage/correlation but do not replace galaxy-level inference.
5. **Control structure and observational selection.** Compare stages on common local-density/radius support and overlapping global mass/inclination/resolution. Treat the small close-to-peak sample as a power limitation. Recheck key claims after excluding combined/interacting or AGN/outflow-dominated cases one at a time, without assuming every such galaxy should be removed from the primary sample.
6. **Select the minimal residual-emission diagnostic set.** Keep absolute Halpha plus NII/SII/OIII changes and one controlled BPT distribution. Report detection attrition and noisy positive-flux tails. Decide whether nondetections require a forward/noise model for a particular claim; do not turn substituted ND values into upper limits. If proposing quantitative source fractions, require a photon budget and held-out line checks.
7. **Draft the empirical results before extending the physical model.** Build Figures 1-7 and the two sample/effect tables from the frozen products, then write the discussion to their actual strength. The detailed gas-model and NGC4654 branches should not delay a sound core paper unless they are made central claims.

### 10.2 High-value discriminating checks

| Question | Minimal useful test | What a contrary result would change |
|---|---|---|
| Is the SF deficit only loss of selected area? | Final T/O decomposition and common-mask survivor residuals | If intensity is unchanged, narrow the title/story to occupancy |
| Is the density signal just radial structure? | Joint density-radius matching with product counts | If it disappears, describe radial/structural dependence without an independent density claim |
| Is the NSF shift only a changed brightness distribution? | Per-product fixed-Halpha contrasts on matched support, narrower bins, actual brightness balance | If removed, emphasize population mixing rather than source-spectrum evolution |
| Can a fixed young-plus-old mixture explain all ratios? | Photon-normalized templates, common denominators, held-out lines, raw-flux forward comparison | If it fails, keep it as a limited explanatory model rather than a measured source decomposition |
| Is higher dispersion a population result? | Robust medians, valid LSF subtraction, distinct gas-bin support, beam-smearing/nuclear sensitivity | If not, do not invoke universal turbulent/shock heating |
| Is NGC4654 absolutely enhanced? | Stable external reference and held-out spatial comparison under the final reduction | If still inconclusive, retain asymmetry only or move the case to an appendix |

### 10.3 Claims to retire or rewrite now

Do not carry forward the early claims that the entire pre-peak sample is enhanced above PHANGS, that NSF grows globally with stage, that every forbidden/Balmer ratio rises, that all post-peak NSF gas has larger widths, or that declining NSF mean brightness rules out new excitation. Do not equate TNSF's complement with Composite, low EW with a pure HOLMES spectrum, a BPT distance with source fraction, a direction fit with an independent wind measurement, or an infall class with a calibrated time.

These revisions do not remove the paper's central result. They make it precise: **the observed transformation has a strong occupancy component, a location-dependent surviving-intensity component, and a residual-emission population whose excitation changes require more than one simple universal explanation.**

## 11. Deliverable and verification record

The report is accompanied by a supporting asset directory containing the notebook inventory and SHA-256 snapshot, report index, selected extracted tables, the read-only central-check script and log, live central measurements, and the report-figure/PDF build scripts. The notebook and dated-report inventories provide traceability beyond the narrative. The PDF contains the same report body and appendices as the Markdown; no scientific results exist only in the PDF.

**Checks actually performed:** recursive notebook inventory and JSON parsing; source/markdown/text-output extraction; saved-table extraction; targeted mask/estimator inspection; visual inspection of selected occupancy, decomposition, line, dispersion, BPT, mixture, and NGC4654 saved figures; the five-cell live central reproduction; representative FITS header inspection and current product availability; source-link/reference checks; report figure creation; PDF text/layout checks and rendered-page inspection.

**Not performed:** a new spectral fit, complete post-processing run, all-notebook execution, new full bootstrap/permutation/model fitting, an external field-control analysis, a remote 7000-product audit, or a new measurement of pointing bias. Saved uncertainty intervals and old model fits retain their original provenance. A brief memory lookup helped locate the analysis families, but every scientific statement retained here was checked against current local files or cited primary sources; memory is not used as observational evidence.

The scientific synthesis and final acceptance were performed by the main agent. One bounded read-only helper indexed dated report claims and line locators; it did not execute science or make the physical interpretation. Original notebook SHA-256 values are checked again before delivery. Routine font-cache warnings during the isolated live loading did not affect the numerical check.

## Appendix A. Complete notebook inventory

The following inventory includes all 128 notebooks discovered outside checkpoint directories. IDs are stable within this report. E/C records code cells with a non-null saved execution count / total code cells; it is an output-status indicator, not a reproducibility certificate. Family assignments indicate review scope, not scientific validity. In mixed mass/SFR/metallicity notebooks only the mass/SFR part enters this synthesis. The source of every substantive claim is identified by notebook ID and relevant cells in the main text or the source ledger.

### Reduction, fitting and calibration provenance

- **N001** [20260315_check_IC3392_v3tk_runs.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260315_check_IC3392_v3tk_runs.ipynb) - E/C 9/9; 46 saved PNG outputs.
- **N002** [20260427_check_NGC4383_v3tk_runs.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260427_check_NGC4383_v3tk_runs.ipynb) - E/C 9/9; 46 saved PNG outputs.
- **N003** [20260427_check_NGC4419_v3tk_runs.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260427_check_NGC4419_v3tk_runs.ipynb) - E/C 9/250; 46 saved PNG outputs.
- **N004** [20260430_check_KIN.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260430_check_KIN.ipynb) - E/C 13/13; 13 saved PNG outputs.
- **N005** [20260503_check_KIN_IC3392_3runs.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260503_check_KIN_IC3392_3runs.ipynb) - E/C 13/13; 28 saved PNG outputs.
- **N006** [20260503_check_KIN_IC3392_4runs.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260503_check_KIN_IC3392_4runs.ipynb) - E/C 9/9; 19 saved PNG outputs.
- **N007** [20260503_check_KIN_subctraction_IC3392_3runs.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260503_check_KIN_subctraction_IC3392_3runs.ipynb) - E/C 8/8; 3 saved PNG outputs.
- **N008** [20260505_check_GAS_kin_IC3392_5runs.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260505_check_GAS_kin_IC3392_5runs.ipynb) - E/C 11/11; 12 saved PNG outputs.
- **N009** [20260506_check_GAS_kin_IC3392_VELSCALE40_20.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260506_check_GAS_kin_IC3392_VELSCALE40_20.ipynb) - E/C 11/11; 31 saved PNG outputs.
- **N010** [20260506_check_kinematics.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260506_check_kinematics.ipynb) - E/C 18/19; 22 saved PNG outputs.
- **N011** [20260507_check_GAS_pairwise_fluxerr_IC3392.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260507_check_GAS_pairwise_fluxerr_IC3392.ipynb) - E/C 16/16; 25 saved PNG outputs.
- **N012** [20260509_check_GAS_pairwise_fluxerr_IC3392_7000_samebinning.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260509_check_GAS_pairwise_fluxerr_IC3392_7000_samebinning.ipynb) - E/C 19/19; 27 saved PNG outputs.
- **N013** [20260511_check_OIII_IC3392_highSN_sigma_velocity_outliers.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260511_check_OIII_IC3392_highSN_sigma_velocity_outliers.ipynb) - E/C 10/10; 5 saved PNG outputs.
- **N014** [20260511_check_OIII_IC3392_highSN_sigma_velocity_outliers_uptoSII.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260511_check_OIII_IC3392_highSN_sigma_velocity_outliers_uptoSII.ipynb) - E/C 14/14; 6 saved PNG outputs.
- **N015** [20260512_check_OIII_IC3392_highSN_sigma_velocity_outliers_tied_lines.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260512_check_OIII_IC3392_highSN_sigma_velocity_outliers_tied_lines.ipynb) - E/C 10/10; 7 saved PNG outputs.
- **N016** [20260513_check_OIII_IC3392_used_lines_highSN_sigma_velocity_outliers.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260513_check_OIII_IC3392_used_lines_highSN_sigma_velocity_outliers.ipynb) - E/C 15/15; 10 saved PNG outputs.
- **N017** [20260514_check_IC3392_SIII_driver_linecube_diagnostics.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260514_check_IC3392_SIII_driver_linecube_diagnostics.ipynb) - E/C 13/13; 11 saved PNG outputs.
- **N018** [20260514_check_OIII_IC3392_use_ArIII_noSIII_highSN_sigma_velocity_outliers.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260514_check_OIII_IC3392_use_ArIII_noSIII_highSN_sigma_velocity_outliers.ipynb) - E/C 15/15; 8 saved PNG outputs.
- **N019** [20260514_check_OIII_IC3392_use_trick_SIII_tied_highSN_sigma_velocity_outliers.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260514_check_OIII_IC3392_use_trick_SIII_tied_highSN_sigma_velocity_outliers.ipynb) - E/C 15/15; 10 saved PNG outputs.
- **N020** [20260515_check_OIII_IC3392_use_SIII6312_only_highSN_sigma_velocity_outliers.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260515_check_OIII_IC3392_use_SIII6312_only_highSN_sigma_velocity_outliers.ipynb) - E/C 13/13; 8 saved PNG outputs.
- **N021** [20260515_check_OIII_IC3392_use_SIII9068_only_highSN_sigma_velocity_outliers.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260515_check_OIII_IC3392_use_SIII9068_only_highSN_sigma_velocity_outliers.ipynb) - E/C 13/13; 8 saved PNG outputs.
- **N022** [20260515_check_SIII9068_IC3392_onlySIII9068_previous_bins_fit.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260515_check_SIII9068_IC3392_onlySIII9068_previous_bins_fit.ipynb) - E/C 9/9; 3 saved PNG outputs.
- **N023** [20260517_check_OIII_IC3392_use_red_highSN_sigma_velocity_outliers.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260517_check_OIII_IC3392_use_red_highSN_sigma_velocity_outliers.ipynb) - E/C 13/13; 8 saved PNG outputs.
- **N024** [20260518_check_IC3392_background_galaxy.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260518_check_IC3392_background_galaxy.ipynb) - E/C 1/1; 0 saved PNG outputs.
- **N025** [20260518_check_SFH_mass_IC3392_9100_40_use_vs_7000_53_samebin.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260518_check_SFH_mass_IC3392_9100_40_use_vs_7000_53_samebin.ipynb) - E/C 15/15; 6 saved PNG outputs.
- **N026** [20260518_check_SFH_mass_NGC4383_9100_40_use_vs_7000_53_samebin.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260518_check_SFH_mass_NGC4383_9100_40_use_vs_7000_53_samebin.ipynb) - E/C 15/15; 6 saved PNG outputs.
- **N027** [20260518_check_SFH_mass_NGC4419_9100_40_use_vs_7000_53_samebin.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260518_check_SFH_mass_NGC4419_9100_40_use_vs_7000_53_samebin.ipynb) - E/C 15/15; 6 saved PNG outputs.
- **N028** [20260525_initial_further_check.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260525_initial_further_check.ipynb) - E/C 28/28; 17 saved PNG outputs.
- **N030** [20260528_check_ngist_memory_runtime.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260528_check_ngist_memory_runtime.ipynb) - E/C 11/11; 2 saved PNG outputs.
- **N033** [20260529_test_galaxy_color_palette.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260529_test_galaxy_color_palette.ipynb) - E/C 11/11; 5 saved PNG outputs.
- **N034** [20260601_investigate_ngist_sfh_sky_residuals.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260601_investigate_ngist_sfh_sky_residuals.ipynb) - E/C 13/13; 5 saved PNG outputs.
- **N038** [20260605_check_binsize_extended_further.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260605_check_binsize_extended_further.ipynb) - E/C 11/11; 10 saved PNG outputs.
- **N043** [20260610_NGC4254_KIN_ESPEC_failure_chain.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260610_NGC4254_KIN_ESPEC_failure_chain.ipynb) - E/C 14/14; 9 saved PNG outputs.
- **N050** [20260623_check_IC3392_7000_8900.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260623_check_IC3392_7000_8900.ipynb) - E/C 6/6; 2 saved PNG outputs.
- **N051** [20260624_check_NGC4064_7000_8900.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260624_check_NGC4064_7000_8900.ipynb) - E/C 6/6; 2 saved PNG outputs.

### Metallicity branch (science excluded from this report)

- **N029** [20260526_check_te_metallicity_maps.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260526_check_te_metallicity_maps.ipynb) - E/C 10/10; 46 saved PNG outputs.
- **N031** [20260529_check_SII_NII_metallicity_NGC4189.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260529_check_SII_NII_metallicity_NGC4189.ipynb) - E/C 6/9; 5 saved PNG outputs.
- **N032** [20260529_initial_further_check_fake.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260529_initial_further_check_fake.ipynb) - E/C 19/34; 21 saved PNG outputs.
- **N035** [20260602_initial_check_OH_SFR_relation.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260602_initial_check_OH_SFR_relation.ipynb) - E/C 9/9; 61 saved PNG outputs.
- **N036** [20260602_initial_check_rMZR_fake.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260602_initial_check_rMZR_fake.ipynb) - E/C 20/21; 15 saved PNG outputs.
- **N037** [20260604_initial_check_rMZR.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260604_initial_check_rMZR.ipynb) - E/C 22/23; 17 saved PNG outputs.
- **N039** [20260607_initial_check_OH_SFR_relation_NGC4294_NGC4298.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260607_initial_check_OH_SFR_relation_NGC4294_NGC4298.ipynb) - E/C 18/18; 17 saved PNG outputs.
- **N041** [20260608_Moran_MAUVE_EW.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260608_Moran_MAUVE_EW.ipynb) - E/C 2/2; 35 saved PNG outputs.
- **N042** [20260608_Moran_MAUVE_SF.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260608_Moran_MAUVE_SF.ipynb) - E/C 2/2; 28 saved PNG outputs.
- **N044** [20260616_initial_check_rMZR.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260616_initial_check_rMZR.ipynb) - E/C 20/20; 15 saved PNG outputs.
- **N046** [20260617_initial_check_rMZR_infall.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260617_initial_check_rMZR_infall.ipynb) - E/C 16/16; 11 saved PNG outputs.
- **N047** [20260618_initial_check_rMZR_binned_by_EWHa.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260618_initial_check_rMZR_binned_by_EWHa.ipynb) - E/C 16/16; 11 saved PNG outputs.
- **N048** [20260618_initial_check_rMZR_binned_by_sSFR.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260618_initial_check_rMZR_binned_by_sSFR.ipynb) - E/C 16/16; 11 saved PNG outputs.
- **N049** [20260621_initial_check_rMZR_all_metallicity.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260621_initial_check_rMZR_all_metallicity.ipynb) - E/C 7/7; 5 saved PNG outputs.
- **N052** [20260624_initial_check_rMZR_binned_by_gas.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260624_initial_check_rMZR_binned_by_gas.ipynb) - E/C 6/6; 10 saved PNG outputs.
- **N065** [20260812_initial_check_rMZR_binned_by_mass_to_gas.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260812_initial_check_rMZR_binned_by_mass_to_gas.ipynb) - E/C 6/6; 8 saved PNG outputs.
- **N079** [20260824_QC_IC3392_updated_model_grid_rerun.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260824_QC_IC3392_updated_model_grid_rerun.ipynb) - E/C 17/17; 9 saved PNG outputs.
- **N080** [20260824_check_IC3392_new_metallicity_diagnostics.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260824_check_IC3392_new_metallicity_diagnostics.ipynb) - E/C 14/14; 9 saved PNG outputs.
- **N081** [20260824_diagnose_fix_IC3392_NebulaBayes_JY22.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260824_diagnose_fix_IC3392_NebulaBayes_JY22.ipynb) - E/C 11/11; 5 saved PNG outputs.
- **N082** [20260825_check_photonization_model_metallicty.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260825_check_photonization_model_metallicty.ipynb) - E/C 9/9; 7 saved PNG outputs.
- **N083** [20260825_check_photonization_model_rMZR.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260825_check_photonization_model_rMZR.ipynb) - E/C 9/9; 29 saved PNG outputs.
- **N084** [20260826_check_photonization_model_rMZR_extended.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260826_check_photonization_model_rMZR_extended.ipynb) - E/C 9/9; 29 saved PNG outputs.
- **N100** [20260906_check_rMZR.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260906_check_rMZR.ipynb) - E/C 7/7; 40 saved PNG outputs.

### Early structural/WHAN diagnostics

- **N040** [20260607_initial_check_radial_profiles.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260607_initial_check_radial_profiles.ipynb) - E/C 10/10; 31 saved PNG outputs.
- **N045** [20260617_Moran_MAUVE_WHAN.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260617_Moran_MAUVE_WHAN.ipynb) - E/C 4/4; 105 saved PNG outputs.

### Spiral/KTZ branch (separate exploration)

- **N053** [20260713_NGC4254_KTZ_spiral_fit.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260713_NGC4254_KTZ_spiral_fit.ipynb) - E/C 15/15; 9 saved PNG outputs.
- **N054** [20260715_NGC4254_spiral.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260715_NGC4254_spiral.ipynb) - E/C 14/14; 7 saved PNG outputs.
- **N055** [20260716_NGC4254_spiral.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260716_NGC4254_spiral.ipynb) - E/C 3/17; 17 saved PNG outputs; stored exception: SyntaxError.
- **N056** [20260717_NGC4254_spiral.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260717_NGC4254_spiral.ipynb) - E/C 17/17; 17 saved PNG outputs.
- **N127** [KTZ_validation.ipynb](/Users/Igniz/Desktop/ICRAR/further/KTZ_validation.ipynb) - E/C 9/9; 7 saved PNG outputs.

### Stage rSFMS and structural controls

- **N057** [20260730_initial_check_rSFMS_rMZR_infall.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260730_initial_check_rSFMS_rMZR_infall.ipynb) - E/C 10/10; 5 saved PNG outputs.
- **N058** [20260803_initial_check_rSFMS_pre_peak_infall.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260803_initial_check_rSFMS_pre_peak_infall.ipynb) - E/C 18/18; 14 saved PNG outputs.
- **N059** [20260803_initial_check_rSFMS_pre_peak_infall_SF.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260803_initial_check_rSFMS_pre_peak_infall_SF.ipynb) - E/C 18/18; 14 saved PNG outputs.
- **N060** [20260804_initial_check_rSFMS_infall_stage_histograms.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260804_initial_check_rSFMS_infall_stage_histograms.ipynb) - E/C 13/13; 10 saved PNG outputs.
- **N061** [20260804_initial_check_rSFMS_near_peak_infall.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260804_initial_check_rSFMS_near_peak_infall.ipynb) - E/C 17/17; 14 saved PNG outputs.
- **N062** [20260804_initial_check_rSFMS_post_peak_infall.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260804_initial_check_rSFMS_post_peak_infall.ipynb) - E/C 25/25; 18 saved PNG outputs.
- **N063** [20260804_initial_check_rSFMS_post_peak_infall_SF.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260804_initial_check_rSFMS_post_peak_infall_SF.ipynb) - E/C 25/25; 18 saved PNG outputs.
- **N064** [20260805_initial_check_rSFMS_infall_stage_weighted_histograms.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260805_initial_check_rSFMS_infall_stage_weighted_histograms.ipynb) - E/C 13/13; 10 saved PNG outputs.
- **N066** [20260813_initial_check_rSFMS_rMZR_galaxy.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260813_initial_check_rSFMS_rMZR_galaxy.ipynb) - E/C 9/9; 38 saved PNG outputs.
- **N073** [20260819_initial_check_rSFMS_infall_stage_histograms.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260819_initial_check_rSFMS_infall_stage_histograms.ipynb) - E/C 10/10; 7 saved PNG outputs.
- **N074** [20260819_initial_check_rSFMS_rMZR_galaxy_Re.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260819_initial_check_rSFMS_rMZR_galaxy_Re.ipynb) - E/C 11/11; 38 saved PNG outputs.
- **N075** [20260820_initial_check_rSFMS_rMZR_galaxy_weighted_by_galaxy.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260820_initial_check_rSFMS_rMZR_galaxy_weighted_by_galaxy.ipynb) - E/C 9/9; 38 saved PNG outputs.
- **N076** [20260820_initial_check_rSFMS_rMZR_galaxy_weighted_by_galaxy_Re.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260820_initial_check_rSFMS_rMZR_galaxy_weighted_by_galaxy_Re.ipynb) - E/C 11/11; 38 saved PNG outputs.
- **N077** [20260821_initial_check_rSFMS_rMZR_galaxy.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260821_initial_check_rSFMS_rMZR_galaxy.ipynb) - E/C 10/10; 39 saved PNG outputs.
- **N078** [20260821_initial_check_rSFMS_rMZR_galaxy_Re.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260821_initial_check_rSFMS_rMZR_galaxy_Re.ipynb) - E/C 12/12; 39 saved PNG outputs.

### HI pressure and stage-context experiments

- **N067** [20260815_recompute_Koppen_criteria_MAUVE.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260815_recompute_Koppen_criteria_MAUVE.ipynb) - E/C 15/15; 4 saved PNG outputs.
- **N068** [20260817_recompute_Koppen_AMIGA_Wang_MAUVE.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260817_recompute_Koppen_AMIGA_Wang_MAUVE.ipynb) - E/C 17/17; 5 saved PNG outputs.
- **N069** [20260817_recompute_Koppen_criteria_MAUVE_ModelA.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260817_recompute_Koppen_criteria_MAUVE_ModelA.ipynb) - E/C 15/15; 4 saved PNG outputs.
- **N070** [20260817_recompute_Koppen_criteria_MAUVE_ModelC.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260817_recompute_Koppen_criteria_MAUVE_ModelC.ipynb) - E/C 15/15; 4 saved PNG outputs.
- **N071** [20260818_recompute_Koppen_criteria_MAUVE_ModelD.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260818_recompute_Koppen_criteria_MAUVE_ModelD.ipynb) - E/C 15/15; 4 saved PNG outputs.
- **N072** [20260818_recompute_Koppen_criteria_MAUVE_ModelE.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260818_recompute_Koppen_criteria_MAUVE_ModelE.ipynb) - E/C 15/15; 4 saved PNG outputs.

### SF occupancy and survivor distributions

- **N085** [20260827_check_SF_fraction.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260827_check_SF_fraction.ipynb) - E/C 15/15; 16 saved PNG outputs.
- **N086** [20260828_check_SF_fraction_pre_peak.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260828_check_SF_fraction_pre_peak.ipynb) - E/C 6/6; 17 saved PNG outputs.
- **N087** [20260829_check_surviving_SF_distribution.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260829_check_surviving_SF_distribution.ipynb) - E/C 13/13; 7 saved PNG outputs.
- **N091** [20260831_check_SF_fraction_pre_peak.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260831_check_SF_fraction_pre_peak.ipynb) - E/C 9/9; 26 saved PNG outputs.
- **N092** [20260831_check_pre_peak_rSFMS_categories.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260831_check_pre_peak_rSFMS_categories.ipynb) - E/C 8/8; 16 saved PNG outputs.
- **N094** [20260901_check_close_to_peak_SF_fraction_and_rSFMS_categories.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260901_check_close_to_peak_SF_fraction_and_rSFMS_categories.ipynb) - E/C 14/14; 27 saved PNG outputs.
- **N095** [20260901_check_post_peak_SF_fraction_and_rSFMS_categories.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260901_check_post_peak_SF_fraction_and_rSFMS_categories.ipynb) - E/C 14/14; 67 saved PNG outputs.
- **N096** [20260901_check_pre_peak_SF_fraction_and_rSFMS_categories.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260901_check_pre_peak_SF_fraction_and_rSFMS_categories.ipynb) - E/C 14/14; 42 saved PNG outputs.
- **N097** [20260902_check_combined_SF_fraction_categories_by_stage.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260902_check_combined_SF_fraction_categories_by_stage.ipynb) - E/C 12/12; 30 saved PNG outputs.
- **N098** [20260903_check_combined_SF_fraction_categories_by_stage.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260903_check_combined_SF_fraction_categories_by_stage.ipynb) - E/C 12/12; 78 saved PNG outputs.
- **N099** [20260904_check_combined_SF_fraction_categories_by_stage_SNR_postfit.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260904_check_combined_SF_fraction_categories_by_stage_SNR_postfit.ipynb) - E/C 12/12; 80 saved PNG outputs.
- **N102** [20260908_check_EW_Halpha_line_scaling_relations.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260908_check_EW_Halpha_line_scaling_relations.ipynb) - E/C 8/8; 2 saved PNG outputs.
- **N103** [20260909_check_Sigma_SFR_intensity_by_stage.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260909_check_Sigma_SFR_intensity_by_stage.ipynb) - E/C 9/9; 5 saved PNG outputs.
- **N104** [20260909_check_combined_SF_fraction_categories_by_stage_SNR_postfit.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260909_check_combined_SF_fraction_categories_by_stage_SNR_postfit.ipynb) - E/C 15/15; 81 saved PNG outputs.
- **N105** [20260909_check_corrected_Halpha_surface_density_SF_NSF_ND_by_stage.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260909_check_corrected_Halpha_surface_density_SF_NSF_ND_by_stage.ipynb) - E/C 17/17; 23 saved PNG outputs.

### Spatial asymmetry and directional comparisons

- **N088** [20260830_check_SF_direction_scan_VIVA.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260830_check_SF_direction_scan_VIVA.ipynb) - E/C 19/19; 63 saved PNG outputs.
- **N089** [20260830_check_SF_gradient_scan_VIVA.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260830_check_SF_gradient_scan_VIVA.ipynb) - E/C 20/20; 93 saved PNG outputs.
- **N090** [20260830_check_SF_spatial_asymmetry.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260830_check_SF_spatial_asymmetry.ipynb) - E/C 14/14; 7 saved PNG outputs.
- **N093** [20260901_check_SF_gradient_scan_VIVA_all.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260901_check_SF_gradient_scan_VIVA_all.ipynb) - E/C 21/21; 154 saved PNG outputs.
- **N106** [20260914_check_SF_gradient_scan_VIVA.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260914_check_SF_gradient_scan_VIVA.ipynb) - E/C 22/22; 103 saved PNG outputs.
- **N126** [20261010_check_NGC4654_blind_SF_asymmetry_enhancement.ipynb](/Users/Igniz/Desktop/ICRAR/further/20261010_check_NGC4654_blind_SF_asymmetry_enhancement.ipynb) - E/C 13/13; 4 saved PNG outputs.

### Line fading, source mixtures and excitation

- **N108** [20260930_check_emission_line_fading_BPT_ratios_Halpha_sigma_SF_NSF_by_stage.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260930_check_emission_line_fading_BPT_ratios_Halpha_sigma_SF_NSF_by_stage.ipynb) - E/C 20/20; 20 saved PNG outputs.
- **N109** [20260930_test2_HII_stage_invariance_test3_common_mixing_weight.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260930_test2_HII_stage_invariance_test3_common_mixing_weight.ipynb) - E/C 21/21; 12 saved PNG outputs.
- **N110** [20261004_Numerical_Simple_Gas_Fading_and_BPT_Model.ipynb](/Users/Igniz/Desktop/ICRAR/further/20261004_Numerical_Simple_Gas_Fading_and_BPT_Model.ipynb) - E/C 14/14; 10 saved PNG outputs.
- **N111** [20261006_check_Halpha_sigma_NSF_TNSF_EWHa_lt3_by_stage.ipynb](/Users/Igniz/Desktop/ICRAR/further/20261006_check_Halpha_sigma_NSF_TNSF_EWHa_lt3_by_stage.ipynb) - E/C 13/13; 15 saved PNG outputs.
- **N112** [20261006_check_emission_line_fading_BPT_ratios_Halpha_sigma_NSF_by_stage.ipynb](/Users/Igniz/Desktop/ICRAR/further/20261006_check_emission_line_fading_BPT_ratios_Halpha_sigma_NSF_by_stage.ipynb) - E/C 25/25; 24 saved PNG outputs.
- **N113** [20261006_check_emission_line_fading_BPT_ratios_Halpha_sigma_TNSF_by_stage.ipynb](/Users/Igniz/Desktop/ICRAR/further/20261006_check_emission_line_fading_BPT_ratios_Halpha_sigma_TNSF_by_stage.ipynb) - E/C 25/25; 24 saved PNG outputs.
- **N114** [20261006_check_the_BPT_line_ratios_NSF_TNSF_EWHa_lt3_by_stage.ipynb](/Users/Igniz/Desktop/ICRAR/further/20261006_check_the_BPT_line_ratios_NSF_TNSF_EWHa_lt3_by_stage.ipynb) - E/C 11/11; 19 saved PNG outputs.
- **N115** [20261007_check_NII6583_sigma_NSF_TNSF_EWHa_lt3_by_stage.ipynb](/Users/Igniz/Desktop/ICRAR/further/20261007_check_NII6583_sigma_NSF_TNSF_EWHa_lt3_by_stage.ipynb) - E/C 13/13; 15 saved PNG outputs.
- **N116** [20261007_check_NSF_BPT_NII_SII_by_stage_and_galaxy.ipynb](/Users/Igniz/Desktop/ICRAR/further/20261007_check_NSF_BPT_NII_SII_by_stage_and_galaxy.ipynb) - E/C 8/8; 27 saved PNG outputs.
- **N117** [20261007_check_corrected_Halpha_surface_density_BPT_PDFs_by_stage.ipynb](/Users/Igniz/Desktop/ICRAR/further/20261007_check_corrected_Halpha_surface_density_BPT_PDFs_by_stage.ipynb) - E/C 27/27; 21 saved PNG outputs.
- **N118** [20261007_check_orthogonal_distance_to_KE01_BPT_PDFs_by_stage.ipynb](/Users/Igniz/Desktop/ICRAR/further/20261007_check_orthogonal_distance_to_KE01_BPT_PDFs_by_stage.ipynb) - E/C 14/14; 7 saved PNG outputs.
- **N119** [20261008_check_corrected_emission_line_BPT_PDFs_SF_NSF_TNSF_partition_by_stage.ipynb](/Users/Igniz/Desktop/ICRAR/further/20261008_check_corrected_emission_line_BPT_PDFs_SF_NSF_TNSF_partition_by_stage.ipynb) - E/C 27/27; 21 saved PNG outputs.
- **N120** [20261008_check_orthogonal_distance_to_KE01_BPT_PDFs_SF_NSF_TNSF_partition_by_stage.ipynb](/Users/Igniz/Desktop/ICRAR/further/20261008_check_orthogonal_distance_to_KE01_BPT_PDFs_SF_NSF_TNSF_partition_by_stage.ipynb) - E/C 14/14; 7 saved PNG outputs.
- **N121** [20261008_check_orthogonal_distance_to_KE01_BPT_equal_galaxy_weighted_PDFs_by_stage.ipynb](/Users/Igniz/Desktop/ICRAR/further/20261008_check_orthogonal_distance_to_KE01_BPT_equal_galaxy_weighted_PDFs_by_stage.ipynb) - E/C 14/14; 7 saved PNG outputs.
- **N122** [20261008_check_six_catalogues_BPT_NII_SII_by_stage.ipynb](/Users/Igniz/Desktop/ICRAR/further/20261008_check_six_catalogues_BPT_NII_SII_by_stage.ipynb) - E/C 13/13; 6 saved PNG outputs.
- **N123** [20261009_check_BPT_line_ratio_density_histograms_SF_NSF_TNSF_intrinsic_Halpha_sigma_gt45_by_stage.ipynb](/Users/Igniz/Desktop/ICRAR/further/20261009_check_BPT_line_ratio_density_histograms_SF_NSF_TNSF_intrinsic_Halpha_sigma_gt45_by_stage.ipynb) - E/C 17/17; 11 saved PNG outputs.
- **N124** [20261009_check_NSF_BPT_line_ratio_histograms_at_fixed_corrected_Sigma_Ha_by_stage.ipynb](/Users/Igniz/Desktop/ICRAR/further/20261009_check_NSF_BPT_line_ratio_histograms_at_fixed_corrected_Sigma_Ha_by_stage.ipynb) - E/C 7/7; 12 saved PNG outputs.
- **N125** [20261009_check_orthogonal_distance_to_KE01_BPT_PDFs_SF_NSF_TNSF_partition_by_stage.ipynb](/Users/Igniz/Desktop/ICRAR/further/20261009_check_orthogonal_distance_to_KE01_BPT_PDFs_SF_NSF_TNSF_partition_by_stage.ipynb) - E/C 18/18; 11 saved PNG outputs.

### Presentation, computing and scratch

- **N101** [20260907_MAUVE_seminar_poster_image_options.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260907_MAUVE_seminar_poster_image_options.ipynb) - E/C 9/9; 14 saved PNG outputs.
- **N107** [20260929_Pawsey_proposal_budget_accounting.ipynb](/Users/Igniz/Desktop/ICRAR/further/20260929_Pawsey_proposal_budget_accounting.ipynb) - E/C 12/12; 2 saved PNG outputs.
- **N128** [tmp_debug.ipynb](/Users/Igniz/Desktop/ICRAR/further/tmp_debug.ipynb) - E/C 12/12; 6 saved PNG outputs.


## Appendix B. Dated-report source register

- **R01** [20260731 rSFMS and rMZR by infall stages.md](/Users/Igniz/Desktop/ICRAR/MAUVE/20260731%20rSFMS%20and%20rMZR%20by%20infall%20stages.md) - First stage-rSFMS impression; line 5.

- **R02** [20260803 Pre peak rSFMS.md](/Users/Igniz/Desktop/ICRAR/MAUVE/20260803%20Pre%20peak%20rSFMS.md) - Selection dependence of the PHANGS comparison; lines 3-15.

- **R03** [20260827 Resolved Star Formation Population through RPS.md](/Users/Igniz/Desktop/ICRAR/MAUVE/20260827%20Resolved%20Star%20Formation%20Population%20through%20RPS.md) - Occupancy and survivor-intensity distinction; lines 3-35.

- **R04** [20260828 SF fraction in Pre-peak.md](/Users/Igniz/Desktop/ICRAR/MAUVE/20260828%20SF%20fraction%20in%20Pre-peak.md) - Individual-galaxy diversity; lines 3-15.

- **R05** [20260906 Where Star Formation Remains along Infall Stages.md](/Users/Igniz/Desktop/ICRAR/MAUVE/20260906%20Where%20Star%20Formation%20Remains%20along%20Infall%20Stages.md) - Strict-S/N spatial category pattern; lines 3-21.

- **R06** [20260831 Pre-Peak rSFMS with SF, ND and NSF.md](/Users/Igniz/Desktop/ICRAR/MAUVE/20260831%20Pre-Peak%20rSFMS%20with%20SF,%20ND%20and%20NSF.md) - ND proxy and category interpretation; lines 3-9; compare current N092.

- **R07** [20260901 Close-to-peak and Post-peak rSFMS with SF, NSF and ND.md](/Users/Igniz/Desktop/ICRAR/MAUVE/20260901%20Close-to-peak%20and%20Post-peak%20rSFMS%20with%20SF,%20NSF%20and%20ND.md) - Historical galaxy-by-galaxy interpretation; later quantitative work qualifies it.

- **R08** [20260909 Remained SF and Increasing NSF along Infall Stages.md](/Users/Igniz/Desktop/ICRAR/MAUVE/20260909%20Remained%20SF%20and%20Increasing%20NSF%20along%20Infall%20Stages.md) - Early synthesis, including wording later corrected; lines 21-59.

- **R09** [20260911_Physical_Scenario_for_Remaining_SF_and_Increasing_NSF.md](/Users/Igniz/Desktop/ICRAR/MAUVE/20260911_Physical_Scenario_for_Remaining_SF_and_Increasing_NSF.md) - Source audit and two-stage resampling; sections near lines 115-177 and 568-578.

- **R10** [20260914_Resolved_Gas_Loss_Star_Formation_and_Ionization_Model.md](/Users/Igniz/Desktop/ICRAR/MAUVE/20260914_Resolved_Gas_Loss_Star_Formation_and_Ionization_Model.md) - Key empirical extraction and attenuation fit; sections 1.3-1.5 and 6.3-6.4; lines 87-165 and 1043-1076.

- **R11** [20260930_HI_Stripping_and_HII_DIG_Line_Ratio_Model.md](/Users/Igniz/Desktop/ICRAR/MAUVE/20260930_HI_Stripping_and_HII_DIG_Line_Ratio_Model.md) - Why proportionally fading young components cannot change their mixture; lines 472-498.

- **R12** [20261002_Model0_Gas_Fading_and_HOLMES_Derivation.md](/Users/Igniz/Desktop/ICRAR/MAUVE/20261002_Model0_Gas_Fading_and_HOLMES_Derivation.md) - Photon normalization and quantitative model tensions; section 8, especially lines 624-853.

- **R13** [20261003 The Origin of Fading Emission Lines and Increasing BPT Line Ratios along Infall Stages.md](/Users/Igniz/Desktop/ICRAR/MAUVE/20261003%20The%20Origin%20of%20Fading%20Emission%20Lines%20and%20Increasing%20BPT%20Line%20Ratios%20along%20Infall%20Stages.md) - Simplified analytical narrative; model time is not the observed stage clock.

- **R14** [20261008 Distribution of SF NSF TNSF BPT Distances Across Infall Stages.md](/Users/Igniz/Desktop/ICRAR/MAUVE/20261008%20Distribution%20of%20SF%20NSF%20TNSF%20BPT%20Distances%20Across%20Infall%20Stages.md) - Additive PDF interpretation; complement shorthand must be checked against the current dual-BPT mask.

- **R15** [20260812 Comprehensive Literature Atlas and Environmental Histories of the MAUVE Virgo Galaxy Sample.md](/Users/Igniz/Desktop/ICRAR/MAUVE/20260812%20Comprehensive%20Literature%20Atlas%20and%20Environmental%20Histories%20of%20the%20MAUVE%20Virgo%20Galaxy%20Sample.md) - Environmental-history context; 40-galaxy scope, not the current optical sample.

- **R16** [20260623 IC3392 7000 vs 8900.md](/Users/Igniz/Desktop/ICRAR/MAUVE/20260623%20IC3392%207000%20vs%208900.md) - Wavelength-range sensitivity; numerical table at lines 16-74.

- **R17** [20260624 NGC4064 7000 vs 8900.md](/Users/Igniz/Desktop/ICRAR/MAUVE/20260624%20NGC4064%207000%20vs%208900.md) - Second range-comparison example; numerical table at lines 16-74.

- **R18** [20261006 Halpha Dispersion along Infall Stages.md](/Users/Igniz/Desktop/ICRAR/MAUVE/20261006%20Halpha%20Dispersion%20along%20Infall%20Stages.md) - Support-sensitive linewidth interpretation; historical TNSF shorthand superseded by current source.

- **R19** [20261001_Connected_HI_Stripping_and_HOLMES_Line_Ratio_Model.md](/Users/Igniz/Desktop/ICRAR/MAUVE/20261001_Connected_HI_Stripping_and_HOLMES_Line_Ratio_Model.md) - Intermediate population-normalized model and retained limitations.

- **R20** [20260830_RPS_Local_Star_Formation_Theory_and_Toy_Models.md](/Users/Igniz/Desktop/ICRAR/MAUVE/20260830_RPS_Local_Star_Formation_Theory_and_Toy_Models.md) - Local enhancement theoretical possibilities; not an observed triggering measurement.

- **R21** [20260831_Molecular_Gas_Supply_Definitions_and_Observational_Inference.md](/Users/Igniz/Desktop/ICRAR/MAUVE/20260831_Molecular_Gas_Supply_Definitions_and_Observational_Inference.md) - Distinguishes reservoir supply, HI-to-H2 conversion, and cloud accretion.

- **R22** [20260518 SFH and Stellar Mass Check.md](/Users/Igniz/Desktop/ICRAR/MAUVE/20260518%20SFH%20and%20Stellar%20Mass%20Check.md) - Mass sensitivity in tests that also vary range and velocity scale; not a pure 7000/8900 correction.

## Appendix C. Verified external context

These are primary papers used for specific context rather than a new systematic literature review. Local measurement claims are supported by the notebook/report records above.

- Huang et al. (2026), *MAUVE-MUSE: When Metallicity Follows or Fights Star Formation -- A Mass-Dependent Inversion in Virgo Galaxies*. [Author preprint](https://arxiv.org/abs/2605.31412); published identifier supplied by the user: [MNRAS stag1019](https://doi.org/10.1093/mnras/stag1019). Used only to delimit Paper I.
- Pessa et al. (2021), *Star formation scaling relations at about 100 pc from PHANGS: Impact of completeness and spatial scale*. [Primary paper](https://arxiv.org/abs/2104.09536). Relevant to survivor selection, resolution, and comparison-sample definitions.
- Boselli, Fossati & Sun (2022), *Ram Pressure Stripping in High-Density Environments*. [Author review](https://arxiv.org/abs/2109.13614). Physical context; not a calibration of the MAUVE stage clock.
- Belfiore et al. (2022), *A tale of two DIGs*. [Publisher paper](https://doi.org/10.1051/0004-6361/202141859). Context for young-star leakage and evolved-star contributions; MAUVE fractions remain to be constrained.
- Kewley et al. (2001), *Theoretical Modeling of Starburst Galaxies*. [Primary paper](https://arxiv.org/abs/astro-ph/0106324). Basis for the adopted KE01 boundaries, not a source-fraction calibration.
- Cid Fernandes et al. (2011), *A comprehensive classification of galaxies in the Sloan Digital Sky Survey: how to tell true from fake AGN?* [Publisher paper](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1365-2966.2011.18244.x). WHAN/weak-line context; resolved purity is not assumed.
- Chung et al. (2009), *VIVA, VLA Imaging of Virgo spirals in Atomic gas: I. The Atlas and the HI Properties*. [Primary paper](https://arxiv.org/abs/0909.0781). HI morphology context for the directional work.
- Brown et al. (2021), *VERTICO: The Virgo Environment Traced In CO Survey*. [Primary paper](https://arxiv.org/abs/2111.00937). Resource for a future gas-amount/efficiency test, not an analysis completed in this report.
