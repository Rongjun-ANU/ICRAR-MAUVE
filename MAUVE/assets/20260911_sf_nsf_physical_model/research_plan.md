# Research plan and evidence boundary

Requested output: an English report dated 11 September 2026, Markdown and PDF,
in the MAUVE directory, explaining the 9 September report and recent results.
The user additionally requested inspection of the relevant further notebooks.

1. Read the named PDF and companion Markdown, relevant recent reports, and
   the three 20260909 notebooks. Trace critical selection and correction code.
2. Recompute selected data cells and central estimators against live FITS,
   without editing notebooks or rerunning the science pipeline. Distinguish
   these checks from stored 10,000-draw bootstrap outputs.
3. Search primary literature for stripping, gas supply, remaining SF, DIG,
   evolved stars, shocks, emission mixing, and counterexamples to passive fading.
4. Derive a minimal gas-loss plus emission-mixture forward model, including
   the actual observation/selection operator. Check analytical limits and
   numerical examples. Do not present an illustrative model as a data fit.
5. Deliver a clear report with local provenance, verified bibliography,
   falsifiable predictions, identifiability limits, and staged fitting procedure.
6. Render mathematical expressions and figures; check all pages, links and
   clipping. Record completed work in the daily job log.

Critical review at scoping: population stage differences are cross-sectional,
NSF is a selection complement, and a fixed plotting window does not ensure
adequate support. A unique physical history is not identifiable from these
profiles alone. The scope is therefore explanatory model construction and a
practical fitting prescription, not a claimed inference of ram pressure or age.

Delegation: one bounded read-only worker checked excitation-paper metadata and
source-supported findings, then identified an existing PDF build workflow.
All scientific interpretation, calculations and final acceptance remain with
the main agent.

## Completion and critical review

All six steps are complete. In addition to selected live data-cell execution,
the report uses a fresh 10,000-draw whole-system bootstrap of scalar extracts,
with both stages resampled. No physical parameters were fitted and source
notebooks/FITS products were not edited.

The evidence review retained the non-monotonic close-to-peak sample, the
high-density [OIII]/Hbeta decreases, changing support, NSF selection circularity,
ND censoring, and the invalid fixed-Balmer-decrement justification. Direct
Hbeta changes the checked common-support contrasts by less than 0.0096 dex.

The final model combines gas retention, a two-reservoir supply calculation,
young-star fading, residual ionization, line mixing and analytic selection
probabilities. Pure truncation, fixed-template fading, universal brightening
and a single stage clock are explicitly tested as inadequate restrictions.

Final acceptance: 27-page PDF and canonical Markdown; 40 numbered equations,
6 figures, 26 source records; numerical checks passed; zero DOM conversion or
overflow errors, broken PDF links, out-of-page words or replacement glyphs.
All pages were inspected in contact sheets and 15 pages at readable resolution;
the final changed pages were re-inspected. Six source fingerprints are unchanged.
See verification_summary.json for hashes and the precise execution boundary.
