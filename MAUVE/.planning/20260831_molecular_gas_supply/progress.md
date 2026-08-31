# Progress

## 2026-08-31
- Read relevant workspace instructions and selected Deep Research, academic source-verification, file-planning, delegation, and PDF instructions.
- Established scope, plan, source classes, success criteria, and initial evidence gaps.
- Located full Paper I on arXiv; exact model reading is next.
- Build environment: pandoc, pdftoppm, pdfinfo, node available. Standard TeX paths unavailable; PDF runtime preflight pending.
- No raw science data processed; no final report yet.

## Source-recovery checkpoint
- Paper I full model recovered from arXiv HTML/PDF and verified against publisher HTML Eq16-24/Table3.
- User's added data inventory incorporated. Typical survey resolutions will be explicitly distinguished from actual unreported data product properties.
- One read-only worker /root/resolution_evidence checking primary instrument/survey sources; no scientific decisions delegated.
- Key new route identified: Bialy2025 formation and dissociation estimators; Burkhart2025 Eos application. Published PDF read through web, arXivv2 PDF downloaded for formula QA.
- Source/definition gaps now include density-column calibration transport to Virgo, gross-vs-net bookkeeping, molecular return fractions, phase metallicity, and rate identifiability from snapshots.
- Bialy journal download did not finish within first yield; avoided relying on missing local PDF, used accessible web PDF and successful arXiv alternative.

## Synthesis and writing
- Reconciled Paper I's source term with gross formation, destruction, molecular boundary transport, phase-specific stellar return and the chosen sink definition.
- Derived chemical and continuity equations, snapshot identifiability, lifetime/throughput distinctions, observable conversion factors, metallicity inversion, nonlinear beam averaging and explicit uncertainty propagation.
- Incorporated all five user-confirmed facilities and the exact MUSE wavelength window. Added sufficiency matrix, ranked additions, a chemical-spectroscopy pilot and an H I upgrade decision procedure.
- Verified core published chemical coefficients and the Paper I timescale/base-10 differential details. Supporting abstract-only sources are labelled as such in Appendix C.
- Read-only worker /root/resolution_evidence completed instrument/survey checks, six bibliographic metadata checks and an editorial/provenance review. No file writes or scientific decisions delegated.

## Build and final verification
- Wrote the complete canonical source and requested Markdown; generated 2 analytical figures and 27 passing mathematical/unit checks with the ICRAR Python environment.
- Rendered with Pandoc MathML, an equation-number Lua filter, offline headless Chromium, and ReportLab/pypdf. No global package installation or configuration changes.
- Diagnosed and corrected legacy TeX commands/control characters, and preserved all 60 visible equation labels. Final Pandoc build passed --fail-if-warnings.
- DOM inspection found no math errors, unconverted math, broken internal anchors, missing images or overflow. The final PDF contains 35 pages, 10 tables and 160 web link annotations.
- Rendered all 35 pages with pdftoppm; inspected four contact sheets and high-resolution pages 1, 4, 7, 10, 13, 16, 17, 18, 19, 21, 22, 27, 29, 30 and 35. Both scientific figures were also viewed directly.
- No visual defects requiring changes were found. Long tables continue with repeated headers and intact rows. Main text and mathematics remain readable at A4.
- Persisted build instructions, claim-to-source ledger and a hash-bound visual QA record. Structural QA verifies source equality, numbered equations, headings, figures and PDF links; it does not validate empirical closures or every URL's future reachability.
- Scientific result: the current inventory supports conditional, scale-defined supply/formation constraints, not a uniquely measured fine-scale gross conversion rate. Actual delivered beams, noise, filters, completeness and masks remain prerequisites for applying the report.
