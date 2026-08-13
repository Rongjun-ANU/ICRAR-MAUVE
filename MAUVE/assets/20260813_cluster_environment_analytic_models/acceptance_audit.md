# Acceptance audit: analytical cluster-environment review

**Audit date:** 2026-08-13 (Australia/Perth)

## Evidence and bibliography

- Seeded DOI records: 80.
- Crossref resolutions: 80 passed, 0 failed.
- Foundational DOI-less records: 1 (Cavaliere & Fusco-Femiano 1976, checked through ADS).
- Primary-source publisher/arXiv/ADS spot checks: 64 of 81 works (79.0 percent).
- DOI coverage in final Markdown: 80 of 80 DOI URLs present.
- Local Markdown link coverage: 13 local targets checked, 0 missing.
- The review explicitly distinguishes analytic physics, semi-analytic population models, simulation-calibrated fits, empirical clocks, validation studies, and tailored individual-galaxy inversions.

## Artifact checks

- Final Markdown: 1,290 lines, 15,838 words.
- Final PDF: A4, 41 pages, 1.1 MB.
- PDF reopened with pypdf: passed.
- Extracted PDF text: 125,747 characters; 0 Unicode replacement characters.
- PDF link annotations: 486; 0 annotation rectangles outside the page boxes.
- Figures: four PNG and four PDF files, all non-empty.
- Reproducible Python scripts: all four passed `py_compile`.

## Visual PDF inspection

- All 41 PDF pages were rendered to PNG with `pdftoppm` at 105 dpi.
- Pages 1-41 were inspected in 11 contact sheets.
- Pages 14, 17, 35, 40, and 41 were additionally inspected at full rendered resolution to test four-column tables, the tailored Virgo table, the full-reference section, final-page flow, and the manifest.
- Result: no observed clipping, overlapping text, missing figures, truncated tables, or blank content pages. Page 5 is intentionally a short continuation of the table of contents.

## Known scope limits

- This is a systematic-style scoping review, not a formal PRISMA review or a proof of universal bibliographic completeness.
- Bibliographic identity was verified for every DOI, while central scientific claims received primary-source spot checks; metadata resolution alone is not treated as scientific validation.
- Individual-galaxy histories remain prior- and model-dependent. The report recommends posterior distributions and competing mechanisms rather than deterministic stage labels.
- PDF output is not tagged for accessibility, as confirmed by `pdfinfo`.
