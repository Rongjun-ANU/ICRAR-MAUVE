# Claim/source ledger

Cutoff/access date: 2026-08-30. Source IDs refer to source_registry.json and Appendix D.
This is a scoping review, not a PRISMA database export. Original derivations are
identified as such and checked numerically where possible. Native selected
search/read batches are archived separately. Source availability is not evidence
that every paper was read cover to cover.

| Claim group | Primary evidence and locator | Evidence/access | Interpretation or limit retained |
|---|---|---|---|
| Forbidden transitions, levels and selection rules | B01-B02 NIST; B05 collision calculation | Official atomic documentation / primary data | Quantum derivation is explanatory; no new scattering calculation claimed |
| Statistical equilibrium and emissivities | B03 PyNeb, B04 hydrogen, B05 O III, B06 O II | Primary descriptions and installed runtime | Original fixed-state calculations are not photoionization equilibrium |
| Thermal/ionization forward chain | B07 Cloudy C17; B48 CL01; B54 Gutkin16 | Primary implementation descriptions | C17 explicitly not labeled latest release |
| Kappa versus local Maxwellian | B08 proposal; B09-B10 kinetic counter-evidence | Opposing primary studies | Spatial temperature variation is distinct from nonthermal velocities |
| Density versus temperature ratios | B01-B06; B62 atomic-data experiment | Atomic model plus original calculation | Low/high-density saturation, zone weighting, and quenching retained |
| Direct O/H ionic sum and ICFs | B11 Izotov06; B15/B16 MUSE examples | Primary method examples | O III temperature is transferred in red-only MUSE route |
| Red thermal relation | B15, inspected arXiv v1 Appendix B | Full source relation read | t(O III)=0.80t(S III)+0.20; 1270 +/- 170 K scatter propagated conceptually |
| Sulfur-derived oxygen | B17; B14 Section methods | Primary text read | S/H plus fixed S/O is not independent direct O/H |
| Semi-direct low-ionization bias | B12 Yates20 | Primary full-text discussion | Up to 0.6 dex in tested regime, not global correction |
| Recombination-line abundances | B06 and B13 | Primary data/measurements | O II RL traces recombining O++; faintness/ADF unresolved |
| PP04 exact fits | B25, equations 1-3 | Full source equations checked by root and worker | Single N II/Halpha; log O3N2; domains and empirical dispersion |
| M13 exact fits | B36, equations 2/4 | Full source equations checked | Te-reference coefficients not confused with CALIFA-ONS fit |
| R23/branch lineage | B20-B22,B24,B26-B27 | Mixed access | Original Pagel79/Z94/PT05 full text not all accessible; no unchecked polynomial reproduction |
| KD02/KK04 | B24, B26; KK04 eq13/17 | Primary definitions and method read | Joint q/OH inference, missing blue O II in nearby MUSE |
| N2O2 | B24; B67 | Model and abundance-pattern evidence | Relatively weak U sensitivity is conditional; N/O/attenuation remain |
| D16 | B37, eq1-3 and model ranges | Full text checked | Correct Ap&SS DOI; optional fifth-order distinct from linear |
| PG16 R/S | B38 eq4-7 | Full text checked by root and worker | Summed O III/N II over Hbeta; branch logN2=-0.6; conditional rms |
| ON/ONS/NS/C | B33-B35; NS11 eq8 | Primary equation/definition check | Branch rules differ; MUSE NS/red-C possible, full ON/ONS unavailable |
| Curti suite | B39/B40; C20 Table2 | All nine coefficient rows source-checked | Polynomial direction, scatter axis, shared-flux covariance retained |
| Hybrid M08/M10 | B31 Table4/eq1; B32 Section2.1 | Primary method read | Hybrid scale and selection/averaging rules, not labels alone |
| S/Ar/Ne indices | B29-B30, B68 | Primary diagnostic papers | Relative-element, temperature and excitation leverage distinguished |
| Local analog and high-z fits | B41-B46 | Primary metadata/methods; exact Bian eq11-12 checked | Population scope distinguished from observed bandpass |
| Genesis status | B45 DOI and arXiv v3 | Publication metadata freshly checked | Published Jan2026, not preprint-only; conditional validation not absolute zero point |
| AURORA status | B46 DOI | Publication metadata freshly checked | 2025 first preprint, 2026 journal record |
| DESIRED | B47 arXiv v1 Tables2-4 | Submitted preprint, full HTML tables | 27 definitions/coverage checked; t2 scales not mixed; no settled-consensus claim |
| T04 / CL01 | B48-B49; T04 Section3.1 and eq1 | Primary text read | 0.03 dex formal uncertainty excludes systematics; R23 compression not full inference |
| Inference engines | B50-B56 originals | Primary methods and versioned HCm page | Grid choice/priors, N/O and U identifiability distinguished |
| Scale offsets / model recovery | B58, B60, B52, B56 | Independent comparisons and synthetic test | 0.7 dex is spread in tested scales, not a universal error floor |
| Stellar benchmark | B61 M83 | Primary independent comparison | Dust, population and elemental-mixture caveats |
| DIG disagreement | B63 global mixing versus B64 paired regions | Opposing/qualified primary tests | Different estimands and sample conditions, no single global correction |
| Temperature structure / ADF | B57, B65, original two-zone calculation | Model test and primary observational evidence | No universal t2 prescription or fully settled physical cause asserted |
| Recent O+ temperature reversal | B66 | Submitted preprint | Unresolved and not generalized to N+/S+ or all metal-rich gas |
| MUSE coverage and AO notch | B18-B19 official ESO; user windows | Official docs plus wavelength arithmetic | 4800 baseline preserved; mode/redshift/detection caveats |
| Other elemental abundances | B67-B71 plus B11/B13 | Primary N/O, S/H, He and multi-element work | Ionic abundance, ICF, depletion and O/H proxy kept separate |
| UV methods | B72-B73 | Primary UV/optical examples | Rest-optical observed-NIR distinguished from rest-UV diagnostics |
| IR methods | B74-B79 | Primary observations/model diagnostics | Weak temperature dependence, density and hydrogen normalization still matter |
| IR OIII/NIII | B75/B77 | Primary ratio definitions | Weighted/unweighted versions distinct; N/O-O/H assumption explicit |
| Phase context | B80 review plus basic definitions | Review only for context | Not used as sole support for numerical calibration claims |

## Closure

All user-requested scientific families have supporting primary evidence and an
explicit applicability/accuracy discussion. Remaining limits: not every historic
polynomial was recovered; preprints remain provisional; no real MAUVE data or
production-calibration performance was assessed. Original numerical checks are
documented in physics_checks.json and verify_physics.py.
