# Task plan: nebular metallicity from atomic physics to MUSE

## Objective and scope

Create detailed English Markdown and PDF reports in /Users/Igniz/Desktop/ICRAR/MAUVE. Explain the forward chain from elemental abundance and radiation to observed fluxes, then the inverse abundance problem. Default to gas-phase oxygen abundance in H II regions and optical MUSE spectroscopy; separately cover other elements and non-optical methods. User clarification: adopt 4800-9300 A as the MUSE baseline; compare fits restricted to 4800-7000 and 4800-8900 A. Search cutoff 2026-08-30. Do not claim exhaustive enumeration of every published coefficient variant.

## Phases

1. Scope, instructions, and evidence architecture: complete.
2. Primary-source discovery and verification: complete.
3. Follow-up on disagreements, U/hardness, semi-direct, and wavelength gaps: complete, residual access gaps explicitly bounded.
4. Physical synthesis and canonical report-source.md: complete.
5. PDF production, equation/source/applicability QA, delivery records, job log: complete.

## Skills and responsibilities

- Explicit Deep Research plugin: overall workflow, gap matrix, canonical source, provenance, delivery.
- Academic Research Suite: literature discovery, source verification, and synthesis phases inline; physics-appropriate evidence grading, not clinical RCT grading. The explicit report brief supplies scope; no extra scope-confirmation pause.
- Planning with files: preserve plan and evidence in this task-specific assets directory.
- PDF: professionally render and visually verify final PDF; operation marker before authoring.
- Parallel dispatch: one read-only source worker for official MUSE and alternative-wavelength/element methods. Root owns all science, formulas, recommendations, and acceptance.
- Verification before completion: fresh numeric, source-format, and artifact checks before acceptance.
- Systematic debugging: reproduce and fix Pandoc legacy-math and Chromium equation-pagination issues; preserve all scientific expressions.

## Evidence protocol

Use original calibration papers and independent validation papers, official atomic/instrument documentation, and primary physical calculations. Reviews are discovery maps, not sole support for consequential claims. Record exact index definition, sum versus single line, valid range, sample, reference abundance scale, conditional scatter versus absolute accuracy, U/hardness/N/O/density/geometry sensitivity, and MUSE applicability. Verify coefficients in full text. Mark unavailable evidence, model assumptions, and unresolved debates.

## Deliverables

- 20260830 Nebular Metallicity - Atomic Physics Calibration Atlas and MUSE Guide.md
- Matching PDF.
- Minimal reproducibility/evidence assets in this directory; intermediate renders under /private/tmp.

## Errors and decisions

- User added a dedicated explanation of density-sensitive versus temperature-sensitive ratios and their connection to abundance. Added a standalone section with two-level derivations, diagnostic/coverage tables, and reproduced SII density ratios.
- PyNeb 1.1.30 numerical bridge completed: isothermal abundance recovered to 0.1%; equal-EM 8000/12000 K mixture gives 11149 K and -0.159 dex ionic bias; coverage checks passed. No observed MAUVE data used.
- PDF operation marker ran successfully exactly once before canonical-report authoring.
- Primary-source gaps retained for original full texts of Pagel79, Z94, PT05; use verified metadata/method-level evidence only and do not transcribe unchecked polynomials.

- Some batched instruction outputs were truncated; re-read the missing content in smaller calls before acting.
- No MAUVE-local AGENTS.md exists; global, Desktop, and ICRAR instructions read.
- Do not modify existing science products, pipelines, notebooks, or older reports; no Git operation needed.
- Local N2 vs S2 papers.md is only a discovery lead: contains mismatched citation links and contradictory wording, so verify independently and do not reproduce its claims uncritically.

## Final acceptance

Both requested reports are saved at the exact MAUVE destination. Final PDF: 55 pages, 80 source entries, 91 display and 180 inline math expressions, 20 tables. Fresh PyNeb numeric checks, warning-free Pandoc conversion, DOM checks, PDF structural/link checks, canonical/final Markdown byte identity, and visual coverage passed. All-page contact-sheet review and 20 detailed page inspections were followed by a two-page wavelength-table correction; the corrected pages were re-inspected and the other 53 final page renders were pixel-identical to the reviewed version. Final hashes and limitations are recorded in verification_summary.json. Existing MAUVE science products and pipelines were not changed.
