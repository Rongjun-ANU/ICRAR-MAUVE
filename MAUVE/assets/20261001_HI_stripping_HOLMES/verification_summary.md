# Verification - 1 October 2026

Final deliverables:

- `../../20261001_Connected_HI_Stripping_and_HOLMES_Line_Ratio_Model.md`
- `../../20261001_Connected_HI_Stripping_and_HOLMES_Line_Ratio_Model.pdf`

The 29-page report contains 65 consecutively numbered equations, 9 tables,
3 scientific figures, 11 literature entries, and a complete symbol glossary.
Earlier dated reports, notebooks, and FITS products were not modified.

## Mathematical and data checks

Executed `MPLCONFIGDIR=/private/tmp/mauve_20261001/mpl
/opt/miniconda3/envs/ICRAR/bin/python model_predictions.py` in the equivalent
absolute-path invocation recorded in the report.

- Gas analytical solutions versus independent `solve_ivp`, including balanced,
  closed-no-RPS, and enhanced-conversion branches: maximum relative difference
  7.684524070439668e-14.
- Filtered young-Halpha solution versus independent response ODE: maximum
  normalized absolute difference 3.326183772855984e-11.
- Non-balanced equal-rate solution: absolute difference 2.4424906541753444e-15.
- Peak condition, positive-supply lower bound, and initial total-mass upper
  bound passed in the implemented examples.
- Unequal Balmer-decrement weight identity agreed to printed machine precision.
- Three-component examples have rising ratios while the forbidden lines fade.
- Six extraction-provenance fingerprints matched. Anchors were regenerated
  from the existing five-line scalar export, whose hash is in the JSON audit.
- Read current `Mass.py` and `SFR+Z.py` to check their common inclination/area
  correction. Supporting-source hashes are saved separately. This does not
  independently validate all historical FITS headers or stellar M/L conventions.

No full sample fitting, bootstrap, new map extraction, stellar-population fit,
photoionization grid, or large science pipeline was run. These checks validate
the algebra and arithmetic, not the physical closures. The selected bright NSF
anchors are not reproduced by the fixed-spectrum local HOLMES benchmark.

## Document checks

- Pandoc with `--mathml --fail-if-warnings` passed.
- Offline isolated Chrome/Playwright rendering passed: 65 display equations,
  293 inline mathematical expressions, no MathML errors, unconverted math,
  overflowing elements, missing figures, or broken anchors.
- Assembled PDF has 29 pages, 43 outline entries, 63 internal links, and
  22 external links; no invalid links, replacement glyphs, out-of-page words,
  or nearly empty pages.
- `pdftoppm -r 95 -png` rendered all pages. Every page was inspected on contact
  sheets, with detailed equation/figure checks. After final changes, image hashes
  showed that only pages 9, 10, 28, and 29 changed; each was inspected again at
  full rendered resolution. Glossary tables and figure captions remain intact.
- Final PDF SHA-256:
  `1af447a2db69b1241b38ab6d4c4b03ca13653675c8b8d62886ae3ed4b366bd4c`.
- Final QA JSON/contact sheets are saved in `pdf_qa/`.

The bundled Python lacked PyMuPDF; the PDF audit was successfully executed with
the existing ICRAR Python. Initial font-cache warnings did not prevent figure
generation; the resulting figures were visually checked.

## Source and editorial review

One bounded read-only worker, `/root/primary_sources`, checked primary-source
normalizations, mass conventions, spectral-model distinctions, and bibliography,
then checked the final Markdown against `numerical_audit.json`. It reported no
concrete numbering, copied-value, notation, or attribution errors. The main agent
owned physical interpretation, derivation, numerical choices, and final acceptance.

The report explicitly distinguishes source-verified equations, model closures,
direct derivations, exported observational means, and conditional numerical
tests. Primary-source limitations, including the inaccessible Lilly PDF and the
indirect check of the Case-B numerical coefficient, are stated in Appendix C.
