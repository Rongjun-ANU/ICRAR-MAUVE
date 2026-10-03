# Atomic spatial-response revision: verification

Verified on 3 October 2026. The October 2 report filenames are retained.

- Section 3.5 now has one heading and five numbered equations, with spatial positions x, x_1, and x_2. The detailed temporal and physical discussion is in Appendix I. Lee and Cramer appear in the final bibliography.
- The spatial difference multiplies the newly supplied, surviving molecular component, excluding the common surviving initial H2 column. This preserves equality at the initial time.
- `spatial_predictions.py` recomputed the atomic-only cases with HI factors 0.5, 1, and 1.5, common initial H2, and matched coefficients. Independent DOP853 integration, supply quadrature, equal/reversed-rate cases, initial conditions, ordering, and peak checks passed. Maximum error was 2.55e-12 in initial-H2 units.
- The leading example peaks at 0.09839 Gyr with SFR/initial SFR = 1.00686. At 1 Gyr, the leading/trailing offsets from the contemporaneous matched reference are +0.01545/-0.01602 dex. These are illustrative predictions, not fitted MAUVE results.
- `check_report.py` passed after final installation. It confirms 84 sequential equations, references and figure paths, numerical-table agreement, unchanged sections apart from equation renumbering, and unchanged pre-existing asset hashes. Its input snapshots are stored in this directory.
- Pandoc/MathML conversion and offline browser rendering passed. `dom_audit.json` records no math errors, broken anchors/images, unconverted math, or overflow.
- `qa_report.py` and Poppler rendering passed. All 35 final PDF pages were visually reviewed; changed sections were additionally inspected at full-page scale. The PDF has no invalid links, replacement glyphs, out-of-page words, or near-empty pages. See `pdf_audit.json` and `revision_audit.json` for final hashes.
- The emission-line calculations and figures were retained; no new observational extraction, full science-pipeline execution, or model fit was performed.
- Delegation: `/root/source_check` performed bounded read-only editorial consistency checks. The main agent performed the scientific reasoning, calculations, edits, and final acceptance.
