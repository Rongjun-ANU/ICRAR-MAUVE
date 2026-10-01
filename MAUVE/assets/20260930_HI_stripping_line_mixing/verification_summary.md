# Verification of the 30 September 2026 report

The report retains the date of the request; final rendering checks completed
after the local date changed to 1 October. No earlier reports, science
notebooks, or FITS products were modified.

## Delivered artifacts

- `../../20260930_HI_Stripping_and_HII_DIG_Line_Ratio_Model.md`
- `../../20260930_HI_Stripping_and_HII_DIG_Line_Ratio_Model.pdf`
- PDF: 34 pages; 68 sequential numbered equations; 69 display-math blocks;
  370 inline-math blocks; 11 tables; 3 figures; 12 primary references.
- Markdown SHA-256:
  `1c4a6272d45b71bfdc111e5d75bd8542e98993fa70bc1fc206a0f5cfa66c6b0b`
- PDF SHA-256:
  `81dd65c517e53d6a3fce8190fa2c70afea69e69383a8a290579bc07b575df667`

## Executed numerical checks

`/opt/miniconda3/envs/ICRAR/bin/python model_predictions.py` completed with all
assertions satisfied. Maximum relative analytic/ODE difference was
6.73303e-12; the equal-rate branch differed by 1.52945e-15. Initial molecular
balance, monotonic balanced decline, positivity, and the zero-supply bound
passed. Synthetic common-weight recovery had zero printed error; the
unequal-Balmer identity differed by 1.11022e-16.

The finite-interval line interpolation is consistent with a shared,
gas-shaped young source only after imposing the inferred photon-allocation
history. All allocation fractions remained in [0,1] and summed to unity.
This verifies budget feasibility, not a measured or RPS-predicted allocation.

An additional independent ODE check of the maintained-external-supply solution
(equation 57) differed by 2.40705e-12; its asymptotic SFR factor is 0.620954.
A numerical evaluation of the exact midpoint decomposition (equation 41)
differed by 5.55112e-17.

Six recorded notebook/pipeline/context-file fingerprints matched the scalar
extraction record. Anchors were regenerated from the unchanged scalar line
export. Large emission-map FITS inputs were not individually rehashed and the
science notebooks were not rerun.

## Document and PDF checks

- Pandoc conversion with `--mathml --fail-if-warnings`: exit 0.
- Offline headless Chrome DOM audit: no failed math conversion, `merror`,
  missing images, horizontal overflow, broken internal anchors, or unintended
  local Markdown links.
- `assemble_pdf.py`: successful 34-page final file and page furniture.
- Poppler `pdftoppm -r 96 -png`: exit 0 for all 34 pages.
- `qa_report.py`: exit 0; no out-of-page words, replacement glyphs, invalid
  links, or near-empty pages. 56 outline entries; 57 internal and 50 external
  PDF links.
- All 34 final pages inspected through three contact sheets; detailed
  inspection included the integrating-factor derivation, explicit SFR
  solution, Balmer-weight algebra, conditional numerical rates, flux-model
  matrix, and multi-line comparison figure/table pages.
- Equation labels are exactly 1--68; 72 explicit singular equation references
  were checked against existing labels. All referenced assets exist.

The initial two TeX conversion errors were corrected before final generation.
Small comparison tables and their captions are kept together; long notation
tables retain repeating headers. Final audits and contact sheets are saved
under `pdf_qa/`.

## Source review and limits

The main agent retrieved both requested conversation texts using the
connected app, inspected live previous reports and the author's November
derivation note, verified primary-source claims, and owns the scientific
interpretation and numerical acceptance. Linked conversation image
attachments were unavailable in textual retrieval and were not claimed as
inspected figures.

One bounded read-only worker, `/root/primary_sources`, checked literature and
editorial provenance. It changed no files. Its actionable notes on the Brown
publication metadata, equation-18 TeX typo, and notation/acronym definitions
were checked and incorporated by the main agent.

The constant-spectrum baseline is a hypothesis. The HI-only gas example
predicts a one-Gyr SFR factor of 0.769991, not the descriptive surviving-SF
factor 0.407. The NSF Halpha/NII interpolation agrees by construction and
fails some other ratio means with the chosen templates. No full MAUVE gas,
source-fraction, selection, or covariance-aware spectral fit was performed.
