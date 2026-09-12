# Claim and inference ledger - 11 September 2026

Reference numbers refer to the dated report and source_registry.json. Local identifiers refer to report Section 10.1. Scientific interpretation and final acceptance were retained by the main agent.

| Claim | Evidence | Status / limit |
|---|---|---|
| 26 products, stage counts 8/5/13; combined pair one unit | N1-N3, live extraction, sample_qc.csv | Fresh local evidence; not nominal catalogue counts |
| Strict SNR_POSTFIT > 25; SF/NSF/ND masks and support | N1 cells 3/5/15 and N2/N3 shared loading code | Fresh source and array checks; NSF is not pure DIG |
| T = O I and occupancy-weighted I | N2 cells 9/11; stage_sfr_profiles.csv | Fresh central data; numerical identity checked |
| J generally differs from F times conditional I | N3 cells 11/13/25 | Algebra plus source inspection; covariance and support differ |
| Inner NSF rises while outer ND rises | Fresh per-system extracts and bootstrap, Figures 1-2 | Descriptive cross-sectional sample; not tracked pixels |
| Remaining SF intensity often declines, but outer survivor intensity need not | fresh_two_stage_contrasts.csv | Pointwise whole-system bootstrap, both stages; selection matters |
| NSF intensity decline does not fix absolute J | Same CSV; direct J estimator | Central J has broad uncertainty; subset brightening not excluded |
| [OIII]/Hb does not rise uniformly | Direct common-support Hb profiles and fresh_true_bpt_contrasts.csv | Actual denominator used; projections and populations differ |
| Fixed corrected Balmer decrement justification is invalid but small here | P1 calculate_BD; N3 cell22; balmer/denominator audit CSVs | Less than 0.0096 dex on checked common support only |
| ND map substitution is not a calibrated upper limit | P1 plus N3 loading | No censored flux likelihood or completeness reconstruction performed |
| Long-pulse pressure and short-pulse impulse differ | R1 relevant analytical sections | Simplified restoring geometry; stage is not pressure/time |
| Gas content and efficiency can both contribute in Virgo | R3 full-text sections; R4 primary abstract | Motivation for MAUVE CO test; not measured by this optical audit |
| DIG area and flux fraction differ | R11 sections 4.3-5.3 | Different classification from MAUVE NSF |
| Residual ionization may mix leakage and old stars | R6 indexed excerpts; R7 full-text sections; R9-R10 | Templates and photon budget must be tested; not all NSF old-star powered |
| EW criteria are resolution/definition dependent | R8, R10, R11 | PROXY_EWHA not assumed calibrated to published WHAN cuts |
| Extra heating occurs in some stripped tails | R12-R13 relevant full-text sections | Not transferred to every inner-disc NSF region |
| Width alone does not establish shocks | R14-R15; width is an SF gate locally | Independent PSF/rotation-corrected decomposition required |
| NGC4064 supports allowing an outflow during quenching | R18 accepted manuscript sections | Individual case; no universal mass loading adopted |
| Stellar fading can decrease [OIII]/Balmer | R16 model sections | Direct counterexample to generic forbidden-line persistence |
| Halpha has a time response and transfer effects | R17, R23; R24 spatial-scale context | Exponential kernel illustrative, not calibrated ages |
| Enhancement and suppression both occur across stripped samples | R19, R20, R25, R26 primary abstracts | Counterevidence to universal monotonic sequence |
| Recent UV and slow-transformation studies do not set MAUVE stage ages | R21-R22 primary abstracts | No universal UV extent ordering or imported clock |
| Gas reservoir solution, mixing relations, lognormal class probabilities | Report Equations 15-37; reproduce_analysis.py | Main-agent derivation under stated assumptions; numerical checks are not an observational fit |
| One residual component can yield fainter emission with larger ratios | Direct toy calculation and linear flux addition | Mathematical possibility; line templates are not identified |
| Old-star-only and shock hypotheses can fail energy constraints | Report Equations 38-39 | Proposed tests; budgets not measured in this run |

## Search and review record

Targeted searches covered ram-pressure analytic models and impulse, VERTICO gas/SFE, VESTIGE quenching, DIG area versus luminosity, PHANGS/MaNGA/WHAN, stripped-tail excitation, line fading, Halpha response, SF enhancement counterexamples, and 2025-2026 relevant work. Primary arXiv, publisher and author/institution records were used. No formal database-completeness or PRISMA claim is made.

A bounded read-only agent checked 12 excitation references, four additional abstract/metadata records, and draft citation/path/figure consistency. It flagged a wrong local filename and two caption/plot mismatches; those were corrected. Main-agent checks covered science, masks, estimators, calibration semantics, uncertainty and final acceptance.

## Adversarial checks retained in the report

- NSF residual classification, ND censoring, width-selection circularity and first-failed-gate dependence.
- Common line support and ratio-of-means versus paired median contrasts.
- Both-stage bootstrap; sparse edge support; correlated marginals; finite system count and covariance rank.
- Non-monotonic close-to-peak values and high-density [OIII] declines.
- Ongoing energy input required; conditional intensity cannot rule out a bright subset.
- Explicit no-fit status; no unique timescale, ionization source or efficiency attribution.
