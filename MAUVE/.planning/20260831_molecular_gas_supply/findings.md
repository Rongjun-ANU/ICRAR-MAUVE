# Evidence and gap matrix

## Initial recovery
- User addition: available MUSE 4800-9300 A, ALMA CO(2-1), MeerKAT HI, UVIT and HST. Report must include present-data feasibility and conditional additional observations; do not assume HST cluster-age filters or high-resolution CO product exists.
- User ADS bibcode: 2026MNRAS.549g1019H. Title: MAUVE-MUSE: When Metallicity Follows or Fights Star Formation -- A Mass-Dependent Inversion in Virgo Galaxies.
- Primary arXiv record: https://arxiv.org/abs/2605.31412 (v1, submitted 2026-05-29; 20 pages; accepted MNRAS).
- Full source: https://arxiv.org/html/2605.31412v1 ; model in Section 4. DOI discovered: 10.1093/mnras/stag1019, publisher access still to verify.
- No MAUVE-specific AGENTS.md found. Global, Desktop, and ICRAR instructions read. No unrelated files changed.
- Memory was used only to recover report/QA conventions; all scientific claims must be checked live.

## Gap matrix
| Claim family | Evidence | Confidence | Gap / next action |
|---|---|---|---|
| Exact Paper I supply definition | Full HTML located | Pending detailed read | Read Section 4 and numbered equations |
| Gross vs net phase conversion | Conservation-law derivation planned | Not yet finalized | Primary chemistry and cycling sources |
| GMC boundary accretion | None yet | Open | Resolved observational measurements + limitations |
| Observable clocks | None yet | Open | Cloud lifecycles, continuity, chemistry, chemical evolution |
| Survey resolution | None yet | Open | Published MUSE/ALMA/MeerKAT/VLA products, brightness sensitivity |
| HI alternatives | None yet | Open | Dust, extinction, absorption, multi-resolution inference |

## Search history
1. Exact ADS bibcode and Huang MAUVE gas supply. Broad bathtub term produced irrelevant results; use exact astrophysical phrases subsequently.
2. Direct arXiv abstract/full HTML and DOI. arXiv succeeded, ADS/DOI failed in web tool.

## Paper I exact model read (arXiv full Section 4)
- Eq16 defines Sigma_g as molecular gas; Eq17: dSigma_g/dt = Sigma_Phi -(1-R)Sigma_SFR - Sigma_out.
- Eq18 carries same-composition sinks; Eq19 defines Z=Sigma_Z/Sigma_g; Eq20: Sigma_g dZ/dt=y(1-R)Sigma_SFR+(Z_Phi-Z)Sigma_Phi.
- Eq21 writes tau_Phi=Sigma_g/Sigma_Phi <= tau_dep. The inequality is an extra condition (Phi>=SFR), not a conservation-law identity; important for any new estimator.
- Table3 calls Sigma_Phi local accretion/replenishment from surroundings. This is broader than chemical H2 formation. Return fraction R is immediately returned to ISM, not automatically the molecular phase.
- Eq26 uses base-10 log differentials with equality: quantitative first-order relation needs ln(10) factor. Discuss only insofar as needed for inference, not a general paper critique.
- Time coordinate organizes spaxel variations; paper does not measure actual time derivatives. Cannot directly treat adjacent spaxels as time samples without a dynamical/lifecycle model.
- Need avoid assuming instantaneous well-mixed recycling or metal-poor local HI; new formulation must separate phase, boundary, and composition.

## First literature wave
- Chevance et al. 2020 primary PHANGS cloud lifecycle: https://academic.oup.com/mnras/article/493/2/2872/5681410 ; 10-30 Myr CO lifetimes and 100-300 pc separations; CO lifetime is not chemical H2 formation time.
- Kim et al. 2022 54-galaxy lifecycle: https://academic.oup.com/mnras/article/516/2/3006/6673429 . Full text follow-up needed.
- Sternberg et al. 2014 primary analytic HI-H2 theory identified, arXiv:1404.5042; follow-up needed.
- OUP journal issue confirms volume549 issue3 July2026, stag1019. Publisher article link still to resolve.

## Consequential new literature (second wave)
- Bialy et al. 2025, ApJ 982:24, DOI 10.3847/1538-4357/adb3a6, The Molecular Cloud Life Cycle I: Constraining H2 Formation and Dissociation Rates with Observations. Original 2024 preprint https://arxiv.org/abs/2408.06416 ; v2 HTML https://arxiv.org/html/2408.06416v2 ; published PDF https://www.mso.anu.edu.au/~krumholz/publications/2025/bialy25a.pdf . Directly addresses target: columns plus fluorescent H2 lines estimate formation/dissociation separately. MUST read published formulae and transfer limits, not just abstract.
- Burkhart et al. 2024, ApJ975:269, DOI10.3847/1538-4357/ad75f8, companion lifecycle simulation; arXiv:2402.01587.
- Burkhart et al. 2025 Eos observational application: https://www.nature.com/articles/s41550-025-02541-7 . Reports formation0.02 and dissociation0.32 Msun pc^-2 Myr^-1, nearby diffuse cloud and coarsely smoothed fluorescent spectroscopy; not a Virgo validation.
- Bialy et al. 2026 cosmic-ray-excited H2 detection: https://www.nature.com/articles/s41550-025-02771-9 . Potential critical contamination/control for IR H2; follow-up needed.
- Fukui et al.2009 LMC CO-HI envelopes, https://arxiv.org/abs/0909.0382 DOI10.1088/0004-637X/705/1/144: inferred 0.05 Msun/yr over10Myr from densities/linewidths/evolution; not direct boundary flux measurement.
- Braine et al.2012 outer M33 cloud, https://arxiv.org/abs/1210.6470 :10pc CO, conditional 2e-4 Msun/yr inflow given orientation. Important example of geometric degeneracy.
- Schneider et al.2023 ionized carbon assembly tracer found at https://www.nature.com/articles/s41550-023-01901-5 ; verify before inclusion.

## Search record additions
3. Resolved cloud accretion + LMC/M33 queries; cloud lifecycle and Sternberg HI-H2 theory.
4. Exact Bialy formation/dissociation title, published DOI, Eos application, 2026 H2 diagnostics. This revealed a direct rate-estimation path missing from a simple regulator-only treatment.

## Runtime
- Direct PDF download failed sandbox DNS; approved escalated curl succeeded for Paper I. Temporary source PDF retained only in task scratch.

## Primary-method reconciliation
- Publisher Paper I recovered via OUP issue link: https://academic.oup.com/mnras/article/549/3/stag1019/8698768 redirects to article-minimal; Eq16-24/Table3 match arXiv content read so far.
- Bialy published PDF read through web tool; arXiv v2 also recovered and PDF downloaded. Mean molecular mass includes helium factor1.4, mbar=2mH*1.4. Adopt report convention1.36 and explicitly rescale any published numerical coefficient.
- Bialy Eq7 true formation=mbar integral R*n*nHI ds; Eq10-11 uses HI-weighted mean Rn. Eq12 closure k0=2e-16 s^-1, alpha1.3: formation=0.14 fHI N21^2.3 Msun pc^-2 Myr^-1 (1.4 helium). NOT a universal observational law; simulated column-density to volume-density calibration.
- Bialy Eq9 FUV photon estimator: diss=mbar*4pi*p/(1-p)*I_tot*tau/(1-exp(-tau)); p~0.15; tau=1.9N21; published numerical0.30*I5*N21/(1-exp(-1.9N21)). IR Eq14 divides by ~3.6 photons/cascade and requires removing shock/CR/formation-pumped emission.
- Published tests: idealized SILCC-Zoom cloud snapshots, four times x3 views; mean relative errors27% formation/31% dissociation, no instrument degradation. Do NOT quote30% as Virgo accuracy. Nonlinear beam averaging requires explicit testing.
- Jeffreson et al.2024 Clouds of Theseus: https://academic.oup.com/mnras/article/527/3/7093/7424987 ; distinguishes cloud identity from molecular survival and high concurrent supply/ejection. Simulations, not observations.
- Glover & Mac Low2007 turbulent H2 formation: https://arxiv.org/abs/astro-ph/0605121 DOI10.1086/512227. Fast dense structures invalidate a single mean-density chemistry time.

## Completed instrument lane and final targeted discovery
- Worker resolution_evidence returned primary-source records for MAUVE/MHONGOOSE/VLA/VIVA/THINGS/HeViCS/UVIT/HST and HI alternatives; no files changed. Root verified consequential MAUVE and MHONGOOSE full-text passages independently.
- MAUVE programme paper Catinella et al2025 DOI10.18727/0722-6691/5393: newer MAUVE-ALMA ~50pc (12m+ACA), MUSE~100-200pc; VERTICO 600-800pc is a different data product. Actual user cubes not inspected.
- MHONGOOSE de Blok et al2024 A&A688A109 DOI10.1051/0004-6361/202348297: Table5 highest-resolution beam8.2x7.1arcsec, rms0.219mJy/beam per1.4km/s, 3sigma logNHI19.77 over16km/s. These are benchmark data, not MAUVE measurements. Maximum MeerKAT baseline7.7km does not give1arcsec HI.
- NRAO official resolution table at1.5GHz gives VLA A1.3arcsec/B4.3arcsec under uniform full-synthesis assumptions; natural weighting coarser and brightness sensitivity/short spacings remain limiting.
- HST filters/depths unknown. UVIT PSF~1.3-1.5arcsec FUV, NUV channel availability date dependent. NUV/U/B/V/I coverage, not arbitrary HST imaging, can support cluster-age fitting.
- Dust FIR beams ~9-37arcsec HeViCS, not high-res HI. Balmer-extinction EDGE-CALIFA calibration assumes constant6Msun/pc2 HI, so cannot validate HI reconstruction independently. NaD existing MUSE neutral-gas kinematics is foreground/abundance/ionization/covering-factor limited; not a total HI mass map.
- Di Teodoro & Peek2021 arXiv:2110.01618: radial HI mass-flux method, generally small flows; boundary transport not phase conversion.
- den Brok et al2021 arXiv:2103.10442 and2025 arXiv:2506.09125: R21 environmental variation;2025 ALMA study is1.7kpc, not100pc universal calibration. CO2-1 alone cannot separate mass/excitation variations.
- JWST official docs: NIRSpec IFU0.6-5.3um,3x3arcsec,R up to~2700; MIRI4.9-27.9um. Neither measures local158um CII; NIRSpec spectral resolution cannot resolve fewkm/s cloud accretion but can constrain excitation with line ratios.
- Bialy et al2026 NatAstron10,540-547 DOI10.1038/s41550-025-02771-9: abstract-level verified detection of CR-excited H2 in B68 with JWST, supports excitation-contamination warning; no detailed unpublished line strengths used.
- Kennicutt & Evans2012 primary review Table1: logC_Halpha41.27, response age0-3-10Myr; FUV0-10-100Myr. These are population kernels, not individual cloud-age clocks. Use5.37e-42 if quoting exact table calibration, not silently mix5.5e-42.
- Incorrect guessed arXiv1301.3500 was unrelated and excluded. Correct Bolatto et al2013 arXiv1301.3498/DOI10.1146/annurev-astro-082812-140944 verified.

## Root synthesis choices (new derivations, not observational findings)
- Recommend explicitly defined effective molecular-reservoir replenishment, with gross chemical formation, dissociation, and molecular boundary flux separate. PaperI sink/recycling convention determines whether Phi is gross or net.
- Existing five datasets enable gas-state/clock/transport constraints and conditional coarse-beam chemistry proxies, but not a unique fine-resolution gross formation/destruction field.
- Two-tier analysis: high-resolution CO+optical+stellar lifecycle at measured common beam, HI-containing quantities at native/coarse HI beam; a hierarchical forward model is an alternative with explicitly prior-dependent fine structure.
- Chemical reaction physics supplies a rate coefficient (time scale) without secular temporal monitoring. FUV or IR fluorescence counts photoexcitation/dissociation only after excitation separation; current broadband UVIT/HST do not provide that measurement.
