# Molecular gas supply report: reproducibility and verification

Prepared 2026-08-31 for the requested MAUVE research report. The two user-facing deliverables are at the MAUVE directory root:

- `20260831_Molecular_Gas_Supply_Definitions_and_Observational_Inference.md`
- `20260831_Molecular_Gas_Supply_Definitions_and_Observational_Inference.pdf`

## Contents

- `report-source.md`: complete canonical source; byte-identical to the delivered Markdown. Figure links are relative to the MAUVE directory, not this asset directory.
- `source_ledger.md`: claim-to-source crosswalk, access limitations and native research provenance where retained.
- `verify_equations.py`, `verification_math.json`: 27 exact/numerical/dimensional checks and two constructed analytical figures. These do not fit or measure MAUVE data.
- `report.css`, `equation_numbers.lua`, `render_report.cjs`, `assemble_pdf.py`: the PDF build chain. The renderer uses an isolated browser without a user profile and blocks HTTP(S) resources.
- `verification_dom.json`: MathML, equation-label, image, link and horizontal-overflow checks before PDF assembly.
- `audit_artifacts.py`, `verification_artifacts.json`: delivered-source equality, PDF structure, equation labels, headings, figure paths, web-link annotations and SHA-256 checks.
- `verification_visual.json`: primary-agent visual acceptance tied to the exact PDF SHA-256. A changed PDF requires a new visual review.

The persistent research plan is in `../../.planning/20260831_molecular_gas_supply/`. Build scratch, downloaded papers and page PNGs were kept in `/private/tmp/mauve-gas-supply-YFFLQ3`, not mixed with science data. Scratch may expire; the report and its two figure assets are self-contained for rebuilding.

## Build environment and commands

The successful build used local Pandoc, Node, Google Chrome, Poppler, the ICRAR Python environment, and the bundled Codex Python/Playwright runtime. No dependencies were installed and no global configuration was changed. `verification_math.json` records Python, NumPy, Astropy and Matplotlib versions.

Run from `/Users/Igniz/Desktop/ICRAR/MAUVE`. The paths below describe this verified host, not a portable environment lock. A browser-launch sandbox approval may be needed; it does not authorize browser-session access or remote resources.

```sh
task_build_dir=$(mktemp -d /private/tmp/mauve-gas-supply-rebuild-XXXXXX)
MPLCONFIGDIR="$task_build_dir/mpl-cache" /opt/miniconda3/envs/ICRAR/bin/python assets/20260831_molecular_gas_supply/verify_equations.py
cp assets/20260831_molecular_gas_supply/report-source.md 20260831_Molecular_Gas_Supply_Definitions_and_Observational_Inference.md
pandoc 20260831_Molecular_Gas_Supply_Definitions_and_Observational_Inference.md --from=markdown --to=html5 --standalone --mathml --lua-filter=assets/20260831_molecular_gas_supply/equation_numbers.lua --toc --toc-depth=2 --metadata=toc-title:Contents --embed-resources --css=assets/20260831_molecular_gas_supply/report.css --resource-path=. --fail-if-warnings --output="$task_build_dir/report.html"
node assets/20260831_molecular_gas_supply/render_report.cjs "$task_build_dir/report.html" "$task_build_dir/body.pdf" assets/20260831_molecular_gas_supply/verification_dom.json
/Users/Igniz/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 assets/20260831_molecular_gas_supply/assemble_pdf.py "$task_build_dir/body.pdf" 20260831_Molecular_Gas_Supply_Definitions_and_Observational_Inference.pdf
pdftoppm -scale-to 850 -png 20260831_Molecular_Gas_Supply_Definitions_and_Observational_Inference.pdf "$task_build_dir/page"
/Users/Igniz/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 assets/20260831_molecular_gas_supply/audit_artifacts.py "$task_build_dir"
```

Do not copy the old passing visual status onto a changed PDF. Inspect all new contact sheets and selected full-size pages, then write a new hash-bound review. The artifact audit marks visual QA pending if its record does not match the PDF.

## Verified result and limitations

- 35 PDF pages; 60 visible numbered display equations; 10 tables; 2 original analytical figures.
- 27 mathematical/unit checks passed, including exact rational budget identities, chemistry factor-of-two conventions, metallicity algebra, unit conversions, numerical coefficients and constructed examples.
- Final Pandoc build passed `--fail-if-warnings`; DOM checks found no math errors, unconverted math, broken internal anchors, missing images or overflow.
- All 35 pages rendered and viewed in contact sheets; pages 1, 4, 7, 10, 13, 16, 17, 18, 19, 21, 22, 27, 29, 30 and 35 inspected at 1600-pixel maximum dimension. Both figures also viewed directly.
- All 60 labels and 81 headings found in extracted PDF text; 160 web-link annotations retained. Heading checks include the contents page and are not a line-by-line semantic proof of all PDF prose.
- Canonical/delivered Markdown equality and final hashes are recorded in `verification_artifacts.json`.
- No actual MAUVE cubes, photometry, filters, beams, masks or noise products were reduced or inferred. Numerical examples are not measurements, fitted constraints or exposure forecasts.
- Source verification is targeted, not a claim of exhaustive literature coverage or continued reachability of every URL. Abstract-only limitations are explicit in the report and ledger.
- The assembled PDF is not accessibility-tagged. It retains readable text, rendered mathematics and clickable web citations, but no stronger accessibility claim is made.

## Scientific handoff

Use separate names for effective molecular-reservoir replenishment, gross chemical formation, net chemical production, CO-visible lifecycle throughput and GMC boundary accretion. Applying the report requires a data-product audit and an explicit closure for unobserved rates. The ranked observing additions depend on which of these quantities is the scientific target; they are not a request to obtain every listed dataset.
