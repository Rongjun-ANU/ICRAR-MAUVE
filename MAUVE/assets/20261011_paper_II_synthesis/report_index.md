# MAUVE dated-report evidence index for second-paper review

Read-only source indexing by bounded subagent on 2026-10-11. All report paths below are relative to `/Users/Igniz/Desktop/ICRAR/MAUVE/` unless specified. Values are **claims or calculations reported in these documents**, not newly verified measurements. No science notebooks, FITS loading, or science pipelines were executed. Local global/ICRAR AGENTS were read; no MAUVE-local AGENTS was found by the targeted search. No memory was used.

## High-priority report sources and exact locators

### Initial stage contrasts and selection dependence

- `20260731 rSFMS and rMZR by infall stages.md:5`: reports pre-peak rSFMS median about 0.3 dex higher than post-peak across log stellar density 7-9. Its metallicity statement at line 7 is out of scope.
- `20260803 Pre peak rSFMS.md:3-15`: initially describes pre-peak as above PHANGS; after performing Pessa21 SF selection without EW/dispersion cuts, the text says the relation becomes more consistent and there is not significant enhancement across the sample. This is an explicit within-report qualification of the enhancement impression.
- `20260804 Pre peak and Post peak rSFMS and their histograms.md:19-23`: pre-peak-reference offsets look systematically lower toward later stages; post-peak scatter and narrower mass support are acknowledged.
- Same report `:31-51`: reports missing low-density (<7) selected SF pixels in post-peak; SFR distributions differ most at log stellar density 7.5-8, less at 8-8.5, and little at 8.5-9. These are visual report claims rather than tests given in the text.
- `20260820 rSFMS and rMZR with Infall Stages within Effective Radius.md:3-7,21,29` and the 20260821 report `:3-7,21,32`: selected HII pixel distribution shifts toward higher stellar density in post-peak; restricting to Re leaves rSFMS slope/normalization visually similar. Re hierarchy is NGVS g, GOLDMine H, with NGC4064 half-mass Re 4.27 +/-1.05 kpc; all selected HII pixels of NGC4598 stated to lie outside Re. Actual catalogue/notebook should govern current use.

### Occupancy and remaining-SF intensity

- `20260827 Resolved Star Formation Population through RPS.md:3-9`: defines SF survival/occupancy fraction in density or radial bins; SF uses HII, EW(Ha)>6 A, intrinsic sigma(Ha)<45 km/s. Says inclination<80 selection gives 26 products with stage counts 8/5/13.
- Same `:13-35`: post-peak occupancy is lower across supported density/radial bins including inside Re. Remaining SF is below pre-peak at intermediate/high density, but surviving-SF radial offsets fluctuate and outer suppression is unclear. This qualification predates the September numerical review.
- `20260828 SF fraction in Pre-peak.md:3-15`: object-level diversity in density versus radial occupancy; NGC4254/4294 maintain 0.6-0.8 SF fraction radially but decline below log stellar density 7.5. NGC4383 and NGC4567_8 omitted from that particular discussion because of inclined starburst/double-system character, not a general current sample exclusion.
- `20260831 Pre-Peak rSFMS with SF, ND and NSF.md:3-9`: opening language calls ND an upper limit, but the category table explicitly says ND is a stored pipeline proxy/substitution, **not a formal upper limit**, and NSF ordinate is Ha-equivalent SFR, not automatically a star-formation measurement.
- `20260901 Close-to-peak and Post-peak rSFMS with SF, NSF and ND.md:3-7`: initial impression is similar/slightly lower ND fraction and increasing NSF; notes optical footprint concern. Many subsequent object comments say ND is non-SF; these are interpretive visual statements and were qualified in later reports.
- Same `:17-59`: five close-to-peak case notes. NGC4654 northwest is called best local-enhancement candidate at line 49; NGC4694 is considered peculiar/post-like at 59. Same `:69-191` records all 13 post-peak cases, with candidate undetected SF explicitly allowed for NGC4380/4394 at lines 101/111. No formal class-independent measurement is given for those diagnoses.
- `20260906 Where Star Formation Remains along Infall Stages.md:3-21`: switches to POSTFIT_SNR>25 and adds equal-galaxy plots. Reports outside-in occupancy loss pre->close, then broad loss close->post; outer/low-density area tends to ND, inner/high-density area to NSF. Approximate visual thresholds: post SF fraction 0.1-0.2 inside Re; nearly zero outside 1.5Re; ND/NSF crossover log density ~8.25 and ~0.6Re. Text itself stresses heterogeneity and only five close-to-peak systems. These approximate crossover values are explicitly not treated as universal thresholds in 20260911.
- `20260909 Remained SF and Increasing NSF along Infall Stages.md:3-13`: display windows log density 7-9.5 and radius 0-2.5Re; continuum S/N>25; 10,000 whole-galaxy bootstrap; equal-galaxy chosen as primary and pooled as supplement.
- Same `:21-27`: separates selected-SF conditional intensity from selected-SF luminosity divided by all usable area; calls these upper/lower boundaries, a claim later explicitly corrected. Reports close-to-peak ~0.2 dex lower than pre over density; post another ~0.2 at 8.2-9.2 and larger decrease at lower density. Radial surviving-SF trend weaker outside Re.
- Same `:39-59`: increasing NSF occupancy with declining conditional Ha intensity; line 45 claims this rules out brighter NSF-component creation, and lines 49-59 broadly describe higher line ratios/width and slower forbidden-line fading. Later reports explicitly narrow these claims (below).

### September 11 source audit and numerical extraction

`20260911_Physical_Scenario_for_Remaining_SF_and_Increasing_NSF.md` is especially useful because it records the numerical extraction and caveats:

- `:15-21`: explanatory model, not fitted evolutionary history; NSF residual class; not all ratios rise; close stage is not consistently intermediate; stage labels are not a calibrated clock.
- `:27-29,603-621`: selected loading/summary cells were executed then, plus fresh 10,000-draw whole-system bootstrap resampling reference and comparison. Not a rerun of full notebooks/pipeline. Recorded identity and analytical checks establish numerical consistency, not physical validation.
- `:37-50`: strict finite SNR_POSTFIT>25, raw line S/N>=3 and flux>=20; 0.25 dex/0.25Re bins; normal occupancy support>=50 retained pixels; conditional categories usually>=20; three-galaxy support not enforced in all original outputs and high-density edge has weak contribution.
- `:54-73`: SF/NSF/ND operational masks; unresolved/nonfinite intrinsic dispersion may fail SF; ND is joint Balmer detection failure; NSF is detected but SF gate fails.
- `:115-127`: pipeline calibration and ND surrogate; no validated true-SFR upper/lower bounds.
- `:137`: at log density 8.625, SF fractions pre/close/post=0.499/0.595/0.311 and NSF=0.469/0.378/0.563. At r/Re=.125 SF=.515/.530/.247 and NSF=.472/.451/.694. At r/Re=1.625 SF=.429/.241/.027 and ND=.317/.586/.789.
- `:143-155`: post/pre log contrasts with 16-84 percentile intervals: at density7.625 occupancy -1.011[-1.403,-.801], surviving intensity -.429[-.522,-.285], product -1.440[-1.783,-1.209]; at8.625 occupancy -.204[-.381,-.029], intensity -.408[-.530,-.261], product -.613[-.799,-.406]. At1.625Re occupancy -1.204[-1.432,-1.023], intensity +.105[-.204,.348], product -1.099[-1.510,-.774]. Outer selected survivors do not show secure fading. Original notebook intensity bootstrap held pre reference fixed; report resamples it.
- `:159-161`: at density8.625 NSF fraction +.080[-.038,.220] dex, NSF intensity -.508[-.715,-.143], contribution J -.351[-.566,-.022]. At central .125Re NSF fraction +.167[.056,.297], intensity -.576[-.941,-.042], J -.206[-.613,.351]. Thus mean class fading does not rule out a brightening hard component in a subset.
- `:165-177`: common detected-line support; line profiles use ratios of equal-system linear line means. At density8.625 post/pre N2 +.181, S2 +.043, O3 -.077 dex (O3 -.151 and -.200 at8.875/9.125); at .625Re N2+.428/S2+.184/O3+.084. O3/Ha proxy is not algebraically equal to O3/Hb because dust clipping does not force corrected decrement2.86. Direct-Hb impact in this specific support is <.0096 dex. Below2.86 observed decrements in NSF are 7.2%,7.2%,27.5%.
- `:568-578`: explicit revisions to September9: all-area SF is decomposition not SFR bound; NSF not unique physical component; ND sensitivity class; distinguish within-stage NSF-SF contrast from absolute NSF stage change; declining conditional mean cannot eliminate subset brightening.

Notebook links recorded by this report (`:594-596`; all under `/Users/Igniz/Desktop/ICRAR/further/`):

1. `20260909_check_combined_SF_fraction_categories_by_stage_SNR_postfit.ipynb`
2. `20260909_check_Sigma_SFR_intensity_by_stage.ipynb`
3. `20260909_check_corrected_Halpha_surface_density_SF_NSF_ND_by_stage.ipynb`

Assets/reproduction package: `assets/20260911_sf_nsf_physical_model/` (`:623-631`).

### September 14 refresh, attenuation fit, and directional candidate

`20260914_Resolved_Gas_Loss_Star_Formation_and_Ionization_Model.md` reports a subsequent fresh extraction:

- `:16-26`: 26 products=27 physical galaxies because NGC4567/8 weighted/resampled as one; stage counts8/5/13. Missing NGC4548/4579/4689 not counted as ND. Operational pre-peak reference is a working field-like assumption, not measured field equivalence. Cross-sectional stage labels do not determine elapsed time.
- `:18-22`: exact domain and thresholds; no added inclination correction; SF HII+EW>6+intrinsic sigma<45; missing EW/sigma can cause failure. ND fails joint Balmer detection.
- `:87-99`: refresh of September11 occupancy/intensity table, same central values and slightly different bootstrap endpoints. In density8.625 NSF fraction contrast interval includes0; at some central bins NSF contribution interval includes0. No blanket luminosity decrease in every inner bin.
- `:116-129`: N2/S2/O3 line-specific contrast table; O3 at density8.625 -.077[-.133,.034] and at .625Re +.084[.022,.231]. Boundary-defined SF/NSF ratio difference is partly built into classification.
- `:133-165`: `20260914_check_SF_gradient_scan_VIVA.ipynb` inspected with selected data/plane cells rerun. NGC4654 G=.02134 dex/kpc, PA319.33deg, 4.33deg from adopted NW315deg; R2=.0444. 342,926 pixels are area support, not independent samples. Facing environmental median +.0354 dex, opposite -.0957, contrast .1311; at fixed315deg +.0357/-.0943/.1300. Facing standardized tail>1:22.83%, opposite14.06%. +.0354 is ~8.5% above control; .1311 side contrast~1.35 is NOT 35% above field. No formal significance/PA uncertainty/spatial-correlated null; this mask lacks continuum>25.
- `:974-1037`: describes fit estimator (see full report for exact support/weights; not independently reproduced here).
- `:1043-1057`: constant intensity attenuation close/pre -.141[-.248,-.059]dex and post/pre -.390[-.468,-.298]; factors .722[.565,.873] and .407[.341,.504]. RMS .028/.094dex; equal-bin estimates -.150/-.405; leave-one-galaxy ranges [-.217,-.086]/[-.420,-.340]. Linear slopes -.048[-.222,.105]/+.102[-.132,.278]dex/dex include0; bins correlated up to.93/.83. No formal model preference or predictive validation.
- `:1061-1076`: attenuation alone cannot identify gas content, efficiency, stripping/supply, or clock; conditional example 300Myr clock ->334Myr decay time,600->668 are not inferred measurements. Fit does not fit occupancy/NSF/ND/BPT/NGC4654.
- `:1096-1108`: cross-sectional, censoring, survivor selection, line-support, energy-budget, phase-closure, and resolution limitations.
- `:1122-1141`: same three Sep9 notebooks plus Sep14 gradient; old fingerprints unchanged then but quantities freshly extracted. Assets `assets/20260914_resolved_RPS_academic_model/` include per-galaxy tables, bootstrap contrasts, fits, leave-one-galaxy files, and scripts `refresh_observations.py`, `inspect_gradient.py`, `analytical_predictions.py`, `fit_attenuation.py`.

### October excitation/dispersion update

- `20261006 Halpha Dispersion along Infall Stages.md:3-9`: reports unreliable high dispersions up to500km/s at low-density/outer post NSF caused by a few spatial bins with many pixels. Adds support requiring25 distinct bins and an unspecified galaxy count (literal sentence omits number). With support cuts, equal-galaxy median shows **no stage trend**; pooled case hard to tell. This supersedes the uncomplicated higher post-peak width impression in Sep9:49.
- Same `:5-7`: historical TNSF described by excluding NII COMP; EW<3 described as HOLMES selection. These are report labels; parent should check current masks and avoid treating an EW cut as unique source proof.
- Same `:15-23`: post/pre line-ratio differences persist but diminish under stricter selection/EW cut.
- `20261007 Distribution of NSF BPT Distances Across Infall Stages.md:3-16`: signed shortest distance to KE01; close/pre PDFs almost identical, post NSF moves positive.
- `20261008 Distribution of SF NSF TNSF BPT Distances Across Infall Stages.md:14-27`: additive NSF partition. SF/NSF each unit-area; TNSF/remainder share NSF denominator so sum exactly to NSF. Report calls remainder COMP, and reports pre/close similar with post positive shift/increased TNSF. It explicitly says current evidence cannot distinguish unveiled old component versus new/harder environmentally affected component. The remainder=COMP statement may not survive a later dual-BPT TNSF mask; current notebook is authoritative.
- `20261009 .md`: empty report; no extractable finding.

## Theory/report chain: distinguish demonstrations, data anchors and fitted results

- `20260812 Environmental Histories of the MAUVE Virgo Galaxy Sample.md:5-21`: literature synthesis says stages mix different mechanisms and are not unique clock. Its40-physical-galaxy master scope is not the 26-product science sample. Uses legacy unresolved citation markers and sandbox image paths; primary references need checking if reused. The companion `20260812 Comprehensive Literature Atlas and Environmental Histories of the MAUVE Virgo Galaxy Sample.md` is a2450-line literature atlas, not a measurement of the second-paper optical statistics.
- `20260813 Analytical Models of Galaxy Evolution in Cluster Environments and Ram Pressure Stripping.md:3-16`: systematic-style scoping review; reports80 DOI papers plus one foundational source; model families and proposed roadmap, not newly measured stripping histories.
- `20260823 Finite-Thickness Arbitrary-Orientation Ram Pressure Stripping - Full Derivation.md:40-66,1619-1641`: derivation/validity-audit report. Its exact claim status should be read from executive/version1 sections; not an observed optical stage fit.
- `20260830_RPS_Local_Star_Formation_Theory_and_Toy_Models.md:13-34`: separates gas amount, efficiency, loss; E=g f e sqrt(d p)/s is explicitly a closure, not universal validated RPS law. `:1161-1209` conclusions/reproducibility cover analytical/toy status.
- `20260831_Molecular_Gas_Supply_Definitions_and_Observational_Inference.md:11-29`: explicitly literature/derivation, no MAUVE rate measurement. Distinguishes reservoir replenishment, HI->H2 conversion and GMC accretion. Instrument inventory is user-supplied; exact per-galaxy beams not established.
- Sep11 (above): analytical model plus separately extracted empirical summaries; not fitted evolutionary history.
- Sep14 (above): quantitative fit only to surviving-SF intensity; gas/class/ionization models remain demonstrations.
- `20260928_Derivation_and_MAUVE_Anchored_Gas_Line_Predictions.md` and `20260929_Detailed_Resolved_Gas_and_Emission_Line_Derivation.md`: data-scaled two-reservoir/line derivation, with explicit sections separating what reproduced versus not (Sep28:747-769; Sep29:1201-1223). These lead into more restrictive Sep30-Oct2 models rather than independent stage-data discoveries.
- `20260930_HI_Stripping_and_HII_DIG_Line_Ratio_Model.md:472-498`: OB fading alone with fixed photon absorption fractions leaves HII/leaked-DIG weights unchanged. `:895-907` explicitly states incomplete self-contained quantitative explanation. This is a change from merely choosing differential component fading to requiring a physical source of changing weights.
- `20261001_Connected_HI_Stripping_and_HOLMES_Line_Ratio_Model.md:9-24,949-959`: adds population-normalized HOLMES term; numerical budget/amplitude and consumption-timescale limits retained.
- `20261002_Model0_Gas_Fading_and_HOLMES_Derivation.md:9-17`: restricted Model0, constant local coefficients, young HII+leaked photons grouped; constant HOLMES and both decrements2.86; explicitly not full explanation. Numerical failures retained, not statistical rejection/clock.
- Oct2 `:624-638`: reuses Sep14 exported line profiles and recorded fingerprints, **does not rerun current FITS**. At density8.625 NSF common-support galaxy counts6/5/13; Ha=3.558e39/2.077e39/1.111e39 erg/s/kpc2; Hb=1.244e39/7.266e38/3.892e38; N2=.3731/.4367/.5663; S2=.3773/.3857/.4168; O3=.8978/.4678/.7528. These ratios of means are distinct from medians/PDFs.
- Oct2 `:640-655`: generous specified local old-mass/photon model at bin center yields Ha_HOLMES=4.050e37, weights1.138% pre/3.647% post. Not a universal ceiling.
- Oct2 `:730-740`: illustrative HOLMES spectra N2=1.5,S2=1,O3=3, young component calibrated at pre-point; post predicted/observed N2=.4017/.5663, S2=.3931/.4168, O3=.9511/.7528. Spectra are assumed sensitivity choices, not measured/grid values.
- Oct2 `:778-792`: inverse fixed-spectrum test demands N2_HOLMES7.99,S2_HOLMES1.94,O3_HOLMES-4.82 under fiducial budget; conditional mean-value failure, not statistical rejection of HOLMES.
- Oct2 `:815-826`: required fading fraction.30426 takes no-supply minimum3.966Gyr, full balanced instantaneous4.223Gyr (response4.226) for gamma_H2=.3/Gyr. Hypothetical1Gyr constraint would require depletion<=.504Gyr under R=.4,no feedback; no observed stage duration inferred.
- Oct2 `:830-853`: data-anchor post/pre luminosity fractions Ha.3122,Hb.3128,NII.4739,SII.3449,OIII.2623; OIII fades faster here. Model cannot reproduce all selected line means with fixed local spectrum/budget; no fit of whole data/occupancy.
- Oct2 `:1418`: explicit scope: reused extraction record, no full map/bootstrap/stellar-population/CO-HI/radiation-transfer/photoionization rerun.
- `20261003 The Origin of Fading Emission Lines and Increasing BPT Line Ratios along Infall Stages.md:7,36,318`: concise analytical narrative; t=0 not observed infall stage, HI-only sink hypothesis, constant efficiency not constant SFR, temporal peak not spatial enhancement measurement. This narrative should not replace Oct2's quantitative failure tests.

## Wavelength-range, mass, EW and kinematic provenance

- `20260401 Updated Proxy EW(Ha).md:5-83,105-110`: continuum-cube Ha-window proxy replaces R-band approximation; per-pixel velocity-shifted air window6555.789-6569.785A and reference6562.795; old pseudo EW retained separately. Visual similarity claimed, not a new quantitative calibration validation.
- `20260420 Halpha dispersion.md:3-15,34`: old instrumental subtraction48.69km/s and final sigma_int<77.5, with14-product fractions31.61% HII,29.89% HII+EW,29.28% final. This threshold/sample is historical; do not mix with current45km/s gate.
- `20260427 v3tk_v7.6.8 and v3_v7.2.1 comparison.md:9-11` and Apr28:9-13: old/new comparisons changed fitting range4800-7000 to4800-9100 and binning/config/version; not pure wavelength-only controls.
- `20260504 Check IC3392 Stellar and Gas Kinematic.md:7-15,53,76`: old/new7000 similar; broader9100 causes more changes; version versus config/range proxies explicitly distinguished.
- `20260505 IC3392 GAS kinematics including 8900 test.md:3-23`: changing end9100->8900 at fixed v768/VELSCALE41 shifts gas velocity maps+10/+11km/s; offset to old7000 reduced from-13/-14 to-3/-4km/s.
- `20260506 IC3392 GAS kinematics VELSCALE40 and VELSCALE20 tests.md:3-25`: tests wavelength/grid sensitivity and nonmonotonic oversampling; OIII tied-SIII concern.
- `20260517 Original spectral masking for spectra higher than 8250A.md:3-16`: red masking relaxation did not remove OIII biased kinematics in selected10 bins.
- `20260518 SFH and Stellar Mass Check.md:79-89,172-182,258-266`: compares7000/53 samebin versus9100/40, not the eventual pure7000/8900 switch. Reported total stellar-mass differences IC3392+.058dex,NGC4383+.178dex,NGC4419-.055dex. Median density shifts8.2787->8.3433,7.0749->7.2259,7.8453->7.8497. Age/EBV/M-L also change.
- `20260623 IC3392 7000 vs 8900.md:16-74`: same-version compared medians and median pixel residuals: age residual+2.940Gyr,stellar velocity-3.278km/s,stellar dispersion-6.965km/s,Ha velocity-2.732km/s; Ha flux medians29.421 vs26.827 (pipeline units not defined in table).
- `20260624 NGC4064 7000 vs 8900.md:16-74`: age residual+1.059Gyr,stellar velocity-2.091,dispersion-5.357,Ha velocity-2.198km/s; Ha flux22.148 vs21.872.
- Targeted case-insensitive search of all2026 MAUVE markdown for pointing/mosaic bias found **no explicit measurement demonstrating less pointing bias at7000**. That motivation comes from the user's latest instruction; local dated reports above establish broader sensitivity to range but do not quantify the requested pointing improvement.

## Out-of-scope and non-result inventory

Metallicity-only or predominantly metallicity reports intentionally excluded from scientific synthesis: April4/7/10 Secondary Relation Is Real; June2 Turnover Galaxies in OH SFR Space; June6 three Carton/Local Metallicity models; June7 Total property in rMZR/rFMR; June18 Bin by EWHa or sSFR (mixed selection context only); June22 Spearman calibration trends; Aug13 FMR/rFMR/gFMR analytical models; Aug25 photonization-model metallicities; Aug30 Nebular Metallicity atlas; Aug31 FMR derivation/inversion. Mixed April13/14 integrated-SFMS reports and July31/Aug20/Aug21 rSFMS+rMZR reports are usable only for their SF/mass content.

JWST proposal/review documents (Sep6,12,13) are planning/context, not additional measured second-paper findings. Computing/storage/installer/uploader documents (Sep29,Oct8) are infrastructure. Oct9 blank report adds nothing. Image-only `20260915 GM Maybe Local SF Enhancement.md` has no textual values; its five images require direct visual inspection to extract any additional claim. May/Junes QC notes concerning error thresholds, masks and foreground positions are methods provenance, not separate environmental results.

## Supplemental precise locators and source URLs

- Inventory contains82 dated Markdown reports from20260401 to20261009; the October9 file is empty. Full path/line-count inventory is `/private/tmp/mauve_second_paper_20261011/report_inventory.tsv`.
- June23 IC3392 report contains comparison-map asset at line5, `assets/image-20260623163127223.png`, and pairwise/residual plot at line81, `assets/image-20260623164313802.png`. June24 NGC4064 assets are at line5 (`assets/image-20260624100957024.png`) and82 (`assets/image-20260624101218978.png`). Their text says ~3km/s and~2km/s offsets, respectively, at line7. No explicit pointing-bias result appears in their prose; images were not scientifically interpreted by this indexing worker.
- October3 conclusion at line550 is stronger than its own conditional setup: it states a linked gas-loss -> SFR-fading -> increasing HOLMES weight/ratios narrative. Preserve the scope at lines7/36/318 and October2 quantitative tests when quoting it.
- Primary-context URLs already cited by local sources, **not freshly web-verified by this indexing worker**: PaperI [DOI](https://doi.org/10.1093/mnras/stag1019) and [author manuscript](https://arxiv.org/html/2605.31412v1), Sep14:1174; Brown2023 VERTICO VII [author manuscript](https://arxiv.org/html/2308.10943v1), Sep14:1177, and [DOI](https://doi.org/10.3847/1538-4357/acf195), Oct3:567; Watts2023 VERTICO V [DOI](https://doi.org/10.1017/pasa.2023.14), Sep14:1180; NGC4654 Lizee2021 [DOI](https://doi.org/10.1051/0004-6361/202038910), [manuscript](https://arxiv.org/abs/2011.10531), Sep14:1192; Vollmer2003 NGC4654 tides/RPS [DOI](https://doi.org/10.1051/0004-6361:20021729), Sep14:1195; Belfiore2022 two DIGs [DOI](https://doi.org/10.1051/0004-6361/202141859), [manuscript](https://arxiv.org/abs/2111.14876), Sep14:1219; Cid Fernandes2011 EW/WHAN [DOI](https://doi.org/10.1111/j.1365-2966.2011.18244.x), Oct3:573; Byler2019 post-AGB spectra [DOI](https://doi.org/10.3847/1538-3881/ab1b70), Oct3:570; Koppen2018 stripping regime [DOI](https://doi.org/10.1093/mnras/sty1610), Sep14:1186; Kruijssen/Longmore2014 resolution/sampling [DOI](https://doi.org/10.1093/mnras/stu098), Sep14:1240.

## Verification of this index

Used `rg --files`, dated-name/line-count inventory, line-numbered reads, targeted `rg` for notebook links, model status, wavelength range, and pointing/mosaic bias. Quoted/extracted report content only; no independent validation of reported science numbers or literature assertions. Parent owns scientific interpretation, current notebook mask/estimator verification, final source selection, and final report acceptance.
