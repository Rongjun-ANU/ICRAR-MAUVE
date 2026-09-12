# MAUVE SF/NSF physical scenario - 11 September 2026

The canonical report is `../../20260911_Physical_Scenario_for_Remaining_SF_and_Increasing_NSF.md`, with a matching PDF. The report uses current local evidence, primary literature, and explicitly illustrative analytical calculations. No physical parameter fit was performed.

## Evidence and calculation files

- `input_fingerprints.json`: SHA-256 fingerprints of the three source notebooks, SFR+Z.py, the master catalogue, and effective-radius input.
- `sample_qc.csv`: loaded system inventory and quality checks.
- `sfr_by_galaxy_bin.csv`, `halpha_by_galaxy_bin.csv`: scalar extracts from selected live notebook data calculations.
- `stage_sfr_profiles.csv`, `stage_halpha_profiles.csv`, `stage_j_profiles.csv`: central stage estimators.
- `stage_bpt_line_profiles_with_hbeta.csv`: common-line-support profiles with H-beta included directly.
- `balmer_live_check.json`, `balmer_common_support_check.csv`, `oiii_denominator_audit.csv`: Balmer-decrement and denominator checks.
- `fresh_summary_bootstrap.csv`, `fresh_two_stage_contrasts.csv`, `fresh_true_bpt_contrasts.csv`: new 10,000-draw whole-system resampling, including reference-stage uncertainty. Inspect support and finite-draw columns before interpreting logarithmic tails.
- `toy_occupancy.csv`, `numerical_checks.json`: model demonstrations and numerical validation of the stated equations.
- `source_registry.json`, `source_claim_ledger.md`: verified source records, access depth, claim support and inference limits.
- `dom_audit.json`, `pdf_audit.json`, `verification_summary.json`: final conversion, structural and acceptance evidence. Contact sheets are retained alongside the final QA record.

## Figure mapping

| Report figure | File stem | Type |
|---|---|---|
| 1 | figure_01_live_occupancy | Fresh observational summary |
| 2 | figure_02_live_decomposition | Fresh observational summary |
| 3 | figure_03_live_line_ratios | Fresh observational summary with direct H-beta |
| 4 | figure_05_mixing_fading_width | Analytical demonstration |
| 5 | figure_04_analytic_occupancy | Analytical demonstration |
| 6 | figure_06_gas_regulator | Analytical demonstration |

Every figure has PNG and PDF versions. Figure filenames 4/5 reflect production order; report captions reflect reading order.

## Demonstration inputs

The Figure 4 mixing example uses A0=1, D0=0.2, tau_A=150 Myr, tau_D=600 Myr, a common component Balmer decrement of 2.86, sigma_A=25 and sigma_D=70 km/s, a 20 km/s centroid separation, and continuum C=0.1 in the matching relative units per Angstrom. Component [NII]/Halpha ratios are 0.3/1.2, [SII]/Halpha 0.2/0.6, and [OIII]/Hbeta 0.5/1.8. The alternative low-[OIII] residual uses 0.2 instead of 1.8. All are illustrative templates.

The Figure 5 inner/outer inputs (median A, log scatter, D0, continuum C, effective Halpha detection threshold, tau_A, tau_D) are respectively `(0.4, 0.85, 0.18, 0.07, 0.08, 300, 600)` and `(0.2, 0.7, 0.01, 0.015, 0.08, 120, 350)`, with times in Myr and amplitudes in common arbitrary surface-luminosity units. The allowed residual fraction is 0.327485 from the width-only example with zero centroid separation; BPT is assumed no more restrictive in this illustration. This differs deliberately from Figure 4's 20 km/s offset.

Figure 6 uses initial diffuse/molecular reservoirs 20/10 in consistent arbitrary surface-mass units, transfer time 300 Myr, molecular removal rate 0.0005 per Myr, depletion time 1500 Myr, returned fraction 0.4 and mass loading 0.3. The diffuse removal rates are 0, 0.004 and 0.012 per Myr. Retention curves use relative peak pressures 0.02/0.08/0.25, effective radial scale 0.65 Re and log column scatter 0.75. These parameters are not measured MAUVE constraints.

## Reproduce the lightweight analysis

Run `reproduce_analysis.py` with the ICRAR Python environment. It reads only the scalar snapshots in this directory, writes the six figures/tables, and checks the gas ODE solution, mass budget, truncated lognormal moment and probability partition. Seed: 20260911. It does not write source FITS products or notebooks.

`inspect_live.py` and `inspect_halpha.py` preserve the isolated extraction procedure used during the audit. They refer to the original temporary workflow destinations and live absolute input paths; they are provenance wrappers, not a general replacement for the source notebooks. The report states the exact source cells exercised. Full map, pooled and paired-diagnostic notebook branches were not all rerun.

## PDF production

The report is converted with Pandoc's native MathML, the supplied `equation_numbers.lua` filter, and `report.css`. The isolated `render_report.cjs` blocks network resources, uses a temporary Chrome profile, and prints the self-contained HTML. `assemble_pdf.py` adds the dated cover, headers and page numbers. `qa_report.py` checks page bounds, glyphs and links and creates contact sheets from rendered page PNGs.

Conversion must pass with no Pandoc warnings. Check and enforce every error list in `dom_audit.json`; the renderer writes those lists but does not itself fail on their contents. Structural checks are supplemented by all-page contact-sheet inspection and readable views of complex pages. The initial forced-chapter pagination was revised to remove nearly empty continuation pages.

No source notebook correction, pipeline rerun, physical fit or observational validation of the proposed mechanism is implied by these numerical/layout checks.
