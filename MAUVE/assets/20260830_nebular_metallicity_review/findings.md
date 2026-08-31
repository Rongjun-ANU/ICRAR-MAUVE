# Findings and gap matrix

## Research question slots

| Slot | Evidence needed | Status |
|---|---|---|
| Atomic states, forbidden emission, collision strengths and recombination | Atomic data/calculations and nebular physics sources | Addressed in Sections 2-6; original PyNeb checks |
| Ionization, heating/cooling, transfer, flux formation | Primary photoionization physics and atomic references | Addressed in Sections 2-4 |
| Direct and semi-direct O/H | Temperature diagnostics, ionic sums, temperature relations, biases | Addressed in Sections 7-8; temperature/ICF limits explicit |
| Strong-line atlas | Exact formulas, calibration samples/domains, independent assessments | Addressed in Section 9; unavailable historical coefficients not reproduced |
| Model diagnostics | Grid parameters/priors, likelihood, identifiability, validation | Addressed in Section 10; grid/engine distinction explicit |
| U versus hardness; NII/SII | Ionization thresholds, controlled-grid evidence, abundance-ratio effects | Addressed in Section 11; conditional sensitivity rankings |
| MUSE capabilities | Official coverage/notch and red auroral implementations | Addressed in Section 13 and per-method notes; all three user windows |
| Non-oxygen and non-optical | Primary N/O, S/H, He/H, C/H and UV/IR papers | Addressed separately in Sections 15-16 |
| Systematics and practical recommendations | Temperature structures, ADF, DIG, dust, averaging and selection | Addressed in Sections 12 and 14; no live-data performance claim |

## Final evidence status

The gap matrix above reflects the completed synthesis. The discovery notes below are a chronological working record: words such as "pending" and "still to read" describe the earlier search stage, not current unresolved tasks. Final evidence/access decisions are in source_claim_ledger.md, source_registry.json, and Appendix C of the report. There are 80 bibliography entries and 40 claim clusters. The review is a broad, source-traceable atlas, not an exhaustive PRISMA review or a transcription of every published calibration polynomial. Original Pagel79, Z94, and PT05 full-text access remained incomplete; no unchecked historical coefficients are presented as verified. Current preprints and conditional accuracy claims are labeled.

The added density-versus-temperature section derives the competing radiative/collisional rates, explains low/high-density saturation and Boltzmann sensitivity, and ties both diagnostics to ionic-abundance emissivities. The numerical demonstration uses fixed-state PyNeb calculations, not a self-consistent photoionization grid or observed MAUVE measurements. Scientific reasoning and final acceptance remain with the main agent.

## Current local context

- Destination exists, contains prior dated detailed Markdown/PDF reviews and assets directory.
- User asks for a literature/physics report, not an audit or rerun of existing FITS products. Discuss 7000/8900 truncation conditionally and do not claim live products were inspected.
- Prior review memory used only to recover report/evidence/QA conventions; all scientific claims will use fresh primary-source checks.
- Pandoc, Poppler and Node are installed; no TeX engine on PATH. Bundled Python dependencies available. Decide equation-safe PDF route during preflight.

## Searches and key evidence

Initial discovery (2026-08-30):

- Marino+2013 original A&A full text/PDF located: https://www.aanda.org/articles/aa/full_html/2013/11/aa21956-13/aa21956-13.html ; DOI 10.1051/0004-6361/201321956. Abstract distinguishes direct-reference O3N2/N2 random scatters 0.18/0.16 dex from smaller ONS-based precision. Must not mix the two equations or accuracy claims.
- Pilyugin & Grebel 2016 original full text read: https://academic.oup.com/mnras/article/457/4/3678/2589035 ; DOI 10.1093/mnras/stw238. R/S use SUMMED [OIII]4959+5007 and [NII]6548+6584 divided by Hbeta; N2 here is not log([NII]6584/Halpha). 313 reference objects selected for consistency with C-method estimates; agreement is relative to Te reference, not proof of absolute truth. Equations 4-7 need transcription check.
- Dopita+2016 original identified: https://arxiv.org/abs/1601.01337 ; exact polynomial and domain still to read. D13/pyqz source https://arxiv.org/abs/1307.5950 explicitly fits q/U and O/H in model grids and discusses kappa distributions.
- Curti+2017 original full text located: https://academic.oup.com/mnras/article/465/2/1384/2417485 ; stacked Te calibration of multiple strong-line indices. Coefficient table still to read.
- Kewley & Ellison 2008 primary study https://arxiv.org/abs/0801.1849 establishes up to 0.7 dex scale variation for its SDSS calibration set; conversions remove inter-calibration mismatch, not absolute errors.
- Luridiana+2015 PyNeb DOI 10.1051/0004-6361/201323152: solves n-level statistical equilibrium for CELs and interpolates recombination emissivities. PyNeb is not a photoionization-forward-model code.
- Nicholls+2012 kappa hypothesis https://arxiv.org/abs/1204.3880 identified; must pair with Maxwellian relaxation counter-evidence before synthesis.

## Outstanding high-impact gaps after first pass

Exact PP04/M13/D16/PG16/Curti/KD02 formula conventions; independent U/N/O sensitivity tests; current direct-method temperature/ADF debates; source-preserving renderer for many equations; MUSE red-auroral papers and official mode windows (worker).

User clarification: default MUSE 4800-9300 A, with 7000 and 8900 A upper fit limits explicitly retained. No alternate baseline.

PG16 equations 6-7 source-read: upper S coefficients [8.424,0.030,0.751,-0.349,0.182,0.508], lower [8.072,0.789,0.726,1.069,-0.170,0.022]; branch log(sum NII/Hbeta) >= -0.6. Calibration-sample rms about 0.048 dex is not universal accuracy.

Curti17 Table 2 source-read: log(R)=sum c_n x^n with x=12+log(O/H)-8.69. For N2 c=[-0.489,1.513,-2.554,-5.293,-2.867]; O3N2 c=[0.281,-4.765,-2.268]. Table distinguishes line-ratio rms from abundance-direction sigma, and valid monotonic branch for each diagnostic. Stacks reduce intrinsic dispersion, so do not equate table sigma with error for individual regions. DOI/full-text and formulas retained in native read_records.

Dopita16 PDF pp. 3-5, equations 1-3 verified: y=log(NII/SII)+0.264log(NII/Halpha); A_O=8.77+y, optional +0.45(y+0.3)^5. Linear up to ~9.05; reference abundance 8.77, log U -3.5 to -2.0, log P/k 5.2-6.7. Explicit N/O(O/H) dependence and uncertain EUV/dust. Use 6584 (paper text contains 6484 typo), summed SII doublet. Fifth-order optional not silently substituted.

Ferland+2016 primary repository https://uknowledge.uky.edu/physastron_facpub/475/ (RevMexAA52,261): electron thermalization fast compared with heating/cooling; challenges kappa as nebular ADF explanation. NIST ASD v5.12 page https://www.nist.gov/pml/atomic-spectra-database documents official levels, wavelengths, A-values, and thresholds.

Worker found Brazzini+2024 (A&A691,A173, DOI 10.1051/0004-6361/202451007; arXiv2410.00106) MUSE red-auroral path. Root verification pending. Worker found Zamora & Diaz2023 (MNRAS525,5767; DOI10.1093/mnras/stad2090), equations 7-8 explicitly substitute SIII9532/9069=2.44; Brazzini instead uses 9069 directly in PyNeb, not that explicit substitution.

Access errors: A&A M13 HTML/PDF and arXiv C20 PDF tool errors; use one alternate primary endpoint rather than repeated retries. No local TeX installed; installed Chrome + Pandoc MathML is available as a potential math-preserving PDF route.

Independent reliability sources discovered:
- Mendez-Delgado+2023 Nature618,249, DOI10.1038/s41586-023-05956-2; https://arxiv.org/abs/2305.11578. Observational support for high-ionization-zone temperature inhomogeneities explaining ADF; do not present a universal settled correction.
- Juan de Dios & Rodriguez2017 https://arxiv.org/abs/1704.06009: 52 atomic-data combinations cause 0.1-0.2 dex differences at low density and >0.6-0.8 dex in some ratios at ne>1e4; not a modern universal error floor.
- Bresolin+2016 M83 https://arxiv.org/abs/1607.06840: direct abundances broadly consistent with young stars after dust correction, empirical strong-line methods can underpredict high-Z stellar values; independent checks have population/element assumptions.
- Yates+2020 A&A634,A107 source identified: https://www.aanda.org/articles/aa/pdf/2020/02/aa36506-19.pdf ; semi-direct T2-T3 prescriptions can underpredict O/H by up to 0.6 dex in low-ionization systems in their sample. Exact scope to verify.
- KD02 https://arxiv.org/abs/astro-ph/0206495: explicit abundance+q solving; N2O2 low-U-sensitivity claim above roughly half solar applies under model N/O relation, not general independence from hardness/N/O.
- Perez-Montero2014 https://academic.oup.com/mnras/article/441/3/2663/1133246 DOI10.1093/mnras/stu753: HCm matches Te well with aurorals; without aurorals needs empirical O/H-U/N/O restrictions for tight precision. Its chi-square weighting is not automatically a fully specified Bayesian posterior.
- Lopez-Sanchez+2012 https://academic.oup.com/mnras/article/426/4/2630/1010699 is a useful forward-model recovery benchmark (not independent proof that a grid's abundance zero point is physically true).

Model-method sources now located: BOND Vale Asari+2016 https://arxiv.org/abs/1605.01057 (DOI10.1093/mnras/stw971) frees N/O, stellar age and geometry; uses ArIII/NeIII to break branch and HeI/Hbeta for hardness, hence not full-information nearby-MUSE-only. NebulaBayes Thomas+2018 https://arxiv.org/abs/1803.00740 DOI10.3847/1538-4357/aab3db is a generic Bayesian grid interface originally demonstrated with Seyferts and also HII; not intrinsically one metallicity calibration. Ji & Yan2022 A&A659,A112 primary PDF located (aa42312-21); optical 3D surface comparison yet to read.

Parent read Brazzini24 arXiv HTML: 95/31497 simultaneous reliable auroral detections (0.30%), 969 NII5756,173 SIII6312,427 redOII. Its sample is selected auroral detections, not a completeness result for MAUVE. Methods/diagnostic-U details need targeted section reads. M13 arXiv PDF now accessible at https://arxiv.org/pdf/1307.5316. C20 and JiYan publisher tool errors require alternate arXiv endpoints.

Correction to preliminary worker note: Brazzini24 does mention SIII9531; worker verified 2.469 for SIII9531/9069 in Section3.3 U proxy. Do not conflate that convention with Zamora&Diaz2023's 2.44 thermometer substitution. Root will verify these exact sections.

Atomic source discovery: NIST selection rules and emissivity https://www.nist.gov/pml/atomic-spectroscopy-compendium-basic-ideas-notation-data-and-formulas/atomic-spectroscopy ; Storey&Hummer1995 original https://adsabs.harvard.edu/pdf/1995mnras.272...41s DOI10.1093/mnras/272.1.41; Storey,Sochi&Badnell2014 OIII collision strengths https://strathprints.strath.ac.uk/51191/ ; OII recombination coefficients Storey+2017 https://arxiv.org/abs/1703.09982. MIT Press official book record confirms Osterbrock&Ferland second edition (publisher lists Nov2005; conventional citation2006); book not fully source-read, so list as further reading rather than pretend equation-level access.

Primary kinetic counter-evidence: Draine & Kreisch2018 https://arxiv.org/abs/1803.10003 confirms local Maxwellian through relevant excitation energies; temperature structure distinct from nonthermal velocity distribution. Cloudy C17 paper https://arxiv.org/abs/1705.10877; not claimed latest version. JiYan22 available via https://arxiv.org/abs/2110.00612 ; IZI via https://arxiv.org/abs/1410.8146 and original author software page https://users.obs.carnegiescience.edu/gblancm/izi/.

Installed ICRAR Python has pyneb, numpy, scipy, matplotlib, fitz. Plan optional reproducible fixed-state forward/inverse emissivity calculation, explicitly not a full thermal/ionization-equilibrium grid or observed data.

Recent follow-up found Easeman,Schady,Wuyts&Yates2024 (online2023), MNRAS527,5484 DOI10.1093/mnras/stad3464: https://academic.oup.com/mnras/article/527/3/5484/7420514. Primary study of 671 MAD regions uses sulfur-Te-based O/H reference, finds D16 least clear U trend among tested diagnostics but ~0.2 dex offset; N2/O3N2 compress range and U-correlate, sub-HII artifacts. Crucial caveat: reference is sulfur abundance converted via S/O, not a universal independent oxygen truth. Root will read exact S/O and temperature assumptions.

Contemporary high-z calibrations located: Nakajima+2022 EMPRESS https://arxiv.org/abs/2206.02824 (6.9-8.9); Sanders+2024 paper first arXiv2023 https://arxiv.org/abs/2303.08149 (46 galaxies,z1.4-8.7,7.0-8.4); AURORA later sample https://arxiv.org/abs/2508.10099 (139 galaxies,41 new; publication/date verified at Caltech https://authors.library.caltech.edu/records/4nwj8-0w395). Distinguish rest-optical observed-NIR from truly rest-UV/IR abundance physics.
