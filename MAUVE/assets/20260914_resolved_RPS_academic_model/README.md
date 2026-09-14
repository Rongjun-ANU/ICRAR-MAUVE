# Reproduction record: resolved gas-loss and ionization report

Report date: 14 September 2026. The canonical Markdown and PDF are in the parent MAUVE directory under the basename `20260914_Resolved_Gas_Loss_Star_Formation_and_Ionization_Model`.

## Scope

The report refreshes the current stage statistics and the explicitly requested 14 September NGC4654 gradient calculation, develops a spatial gas-reservoir and emission model, and fits the normalization of the surviving-SF intensity profile. The pre-peak sample is treated as the field-like control under the user's working assumption.

The source notebooks, pipeline and catalogues were not edited. Seven recorded source fingerprints were checked again at completion. FITS maps were read through the notebook loaders; the individual large map files were not all fingerprinted. Selected data cells were executed, not the complete notebooks or their full plotting galleries. Exact cells and masks are documented in Appendix A of the report.

## Numerical products

- `input_fingerprints.json`, `ngc4654_audit.json`, and `final_source_invariance.json`: input identities and final invariance checks.
- `sfr_by_galaxy_bin.csv`, `halpha_by_galaxy_bin.csv`, and `fresh_*`: refreshed galaxy-level estimates and 10,000-draw whole-galaxy bootstrap contrasts.
- `gradient_summary.csv`, `preferred_sides.csv`, `alignment.csv`, and `ngc4654_audit.json`: the exact requested gradient estimator and hemisphere comparison, including the fixed VIVA direction.
- `analytical_predictions.csv`, `analytical_prediction_landmarks.csv`, and `line_mixing_predictions.csv`: illustrative predictions, not fitted physical parameters.
- `analytical_checks.json`: independent ODE and quadrature comparisons plus category-partition checks.
- `attenuation_fit_results.json`, `attenuation_fit_profiles.csv`, and `attenuation_leave_one_galaxy.csv`: the new partial fit, whole-galaxy bootstrap intervals, mass-dependent sensitivity, and leave-one-galaxy results.
- `figure_01` through `figure_07`: standalone PNG and PDF figures used in the report.
- `primary_source_ledger.md`: the bounded literature worker's source audit. The report bibliography records additional sources inspected by the main agent and the depth of inspection.

The fitted post-peak surviving-SF intensity factor is 0.407, with a 16th-84th percentile interval of 0.341-0.504. This fit constrains a descriptive attenuation component. It does not fit a stripping timescale, supply rate, source decomposition, or complete SF/NSF/ND forward model. The illustrative realization fails to reproduce the declining NSF luminosity contribution in some measured bins; that failure is retained in the report.

## Recompute the science products

Run these commands sequentially from this asset directory. They use the installed ICRAR environment and read the existing local products. The extraction scripts write report assets only.

```sh
SCI_PY=/opt/miniconda3/envs/ICRAR/bin/python
"$SCI_PY" inspect_live.py
"$SCI_PY" inspect_halpha.py
"$SCI_PY" observational_bootstrap.py
"$SCI_PY" inspect_gradient.py
"$SCI_PY" analytical_predictions.py
"$SCI_PY" fit_attenuation.py
```

These are the exact adapted extraction scripts retained for this report. Alternatively, `refresh_observations.py` reconstructs the first three scripts from the preserved 11 September asset sources, records the input comparison, and executes them. That wrapper therefore additionally depends on `../20260911_sf_nsf_physical_model/`. The new numerical model and fit do not use the earlier report's illustrative-model branch.

Invalid divisions in full-map Balmer arrays can produce warnings outside usable selections during extraction. The retained diagnostic masks exclude those nonfinite values. All calculations reported here completed; warnings were not treated as evidence that excluded pixels had detections. The fit avoids platform-dependent matrix-multiplication warnings using explicit finite-array sums.

## Render the report

Pandoc converts the canonical Markdown to standalone embedded-resource HTML with MathML. The Lua filter numbers the 57 display equations. The isolated Chrome renderer blocks remote resources and uses no personal browser profile. The assembly step supplies a title page, page furniture and navigable PDF outlines.

```sh
TASK_ROOT=/Users/Igniz/Desktop/ICRAR/MAUVE
TASK_ASSETS="$TASK_ROOT/assets/20260914_resolved_RPS_academic_model"
TASK_REPORT="$TASK_ROOT/20260914_Resolved_Gas_Loss_Star_Formation_and_Ionization_Model"
TASK_NODE=/Users/Igniz/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node
TASK_PDF_PY=/Users/Igniz/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3
mkdir -p /private/tmp/mauve_20260914
/opt/homebrew/bin/pandoc "$TASK_REPORT.md" --standalone --toc --toc-depth=2 --mathml --embed-resources --resource-path="$TASK_ROOT" --css="$TASK_ASSETS/report.css" --lua-filter="$TASK_ASSETS/equation_numbers.lua" --fail-if-warnings -o "$TASK_ASSETS/report.html"
"$TASK_NODE" "$TASK_ASSETS/render_report.cjs" "$TASK_ASSETS/report.html" /private/tmp/mauve_20260914/body.pdf "$TASK_ASSETS/dom_audit.json"
"$TASK_PDF_PY" "$TASK_ASSETS/assemble_pdf.py" /private/tmp/mauve_20260914/body.pdf "$TASK_REPORT.pdf"
```

Headless Chrome needs a permitted local process launch in this environment. The science interpreter supplies PyMuPDF and Pillow; the separate bundled PDF interpreter supplies ReportLab and pypdf.

## Verify the output

For a fresh render, use an empty QA directory so old page images cannot be mistaken for current pages.

```sh
mkdir -p /private/tmp/mauve_20260914/qa_review
/opt/homebrew/bin/pdftoppm -r 96 -png "$TASK_REPORT.pdf" /private/tmp/mauve_20260914/qa_review/page
"$SCI_PY" "$TASK_ASSETS/qa_report.py" "$TASK_REPORT.pdf" /private/tmp/mauve_20260914/qa_review
```

`dom_audit.json` checks equations, anchors, images and layout widths. `pdf_audit.json` records final pagination, SHA-256, outlines, links and page-bound checks. The retained contact sheets document the complete visual inspection; individual full-page PNGs are scratch outputs. `verification_summary.json` records final acceptance and its scientific limits. Numerical and rendering acceptance do not establish the physical uniqueness of the proposed explanation.
