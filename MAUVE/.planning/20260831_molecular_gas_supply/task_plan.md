# Molecular gas supply: research plan

## Goal
Deliver a detailed, source-traceable English Markdown and PDF report in MAUVE defining the gas-supply surface-density rate in Huang et al. (2026), deriving observational estimators, and treating HI/CO/optical resolution honestly.

## Scope and assumptions
- Follow the user-selected Deep Research workflow; academic source verification supplements it, not an additional approval gate.
- Literature through 2026-08-31; primary papers, original survey descriptions, official instrument documentation.
- Source recovery before interpretation. Distinguish gross chemistry, net phase exchange, reservoir supply, and cloud-boundary accretion.
- No observational data reduction, source-product changes, or claim of measured MAUVE rates.
- User-confirmed current inventory: MUSE 4800-9300 Angstrom, ALMA CO(2-1), MeerKAT 21 cm, UVIT, HST. Explicitly assess sufficiency and rank additional observations; actual beams, sensitivities, coverage, UV/HST filters are not specified.
- Main agent owns scientific definitions, derivations, uncertainty propagation, reconciliation, and final acceptance. One bounded read-only worker may check instrument/survey evidence.
- Scratch downloads/builds in /private/tmp/mauve-gas-supply-YFFLQ3; final requested formats only at MAUVE top level.

## Phases
1. Recover Paper I and relevant local context: complete.
2. Literature discovery and targeted verification: complete.
3. Derive estimators, identify degeneracies, design observing tiers: complete.
4. Write report; verify math/citations; render and inspect PDF: complete.

## Critical questions
1. Which reservoir and sink terms does Paper I actually use?
2. When can supply mean HI-to-H2 conversion or GMC accretion?
3. Which rate combinations are identifiable from snapshots?
4. What can MUSE + existing CO do, and what does HI add at its native beam?
5. What independent clock or flux is needed to infer a rate?

## Errors and alternatives
- ADS abstract and DOI resolver returned internal errors; arXiv abstract and full HTML are accessible. Locate publisher full text separately.
- First combined output was truncated; selected instructions re-read in bounded chunks.
- xelatex/lualatex not on PATH or /Library/TeX/texbin. Preflight bundled PDF runtime; do not install or change global environment without need.
- Used Pandoc native MathML plus an isolated offline Chromium renderer and a ReportLab/pypdf cover assembly; no TeX installation needed.
- Fixed source-level legacy math commands and control characters after a minimal Pandoc reproduction. Added a Lua equation-label filter because native MathML conversion did not retain visible equation tags.
- Isolated Chromium needed a scoped sandbox escalation. It used no user profile and blocked HTTP(S) resources.

## Acceptance
- Exact Paper I notation, numbered equations, reservoir definition and assumptions reconciled.
- Primary support for consequential literature/instrument claims; explicit inaccessible-source limits.
- Mathematical/unit/sign checks, transparent illustrative examples and no invented observations.
- MD/PDF content agreement; readable equations/tables, links, complete text, PDF visual QA.
- Concise final links and job log.

## Final acceptance record
- Deliverables: 20260831_Molecular_Gas_Supply_Definitions_and_Observational_Inference.md and .pdf in the requested MAUVE directory.
- 35 PDF pages, 60 numbered display equations, 10 tables and 2 original analytical figures.
- 27 exact/numerical/dimensional checks passed. Canonical and delivered Markdown agree byte-for-byte.
- All PDF pages rendered and inspected in contact sheets; 15 selected pages inspected at 1600-pixel resolution, including chemical equations, ranked observations, conclusions and references. No clipping, missing symbols, broken equations or overlapping page furniture found.
- Structural/DOM checks passed; 81 headings and all 60 equation labels found, relative figures resolved, web link annotations retained.
- Source-access limits, model dependence and unknown actual MAUVE product properties explicitly recorded. No empirical MAUVE rates or exposure-time promises made.
- Build scripts, math/DOM/artifact/visual verification records and a claim-to-source ledger retained with the companion assets.
