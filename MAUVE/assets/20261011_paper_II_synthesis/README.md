# Evidence accompanying the 11 October 2026 Paper II synthesis

The report body and its PDF are two directories above this directory. Source notebooks and FITS products were not edited. The report distinguishes live central measurements, saved inference, historical model demonstrations, and proposed analysis.

- `notebook_inventory.json`: all 128 notebooks, saved execution/output status, family assignment, and source SHA-256.
- `report_inventory.tsv`: the 82 dated reports indexed; `report_index.md` contains the read-only claim index; `report_sources.json` records the 22 principal report hashes.
- `live_*.csv`, `live_check.log`, `central_result_check.json`: central results from source cells 3, 5, 7, 9, and 11 of N103; no bootstrap, fitting, or notebook execution-output rewrite.
- `Nxxx_cY_oZ_tW.csv`: selected tables extracted from stored notebook HTML outputs. `c` is the one-based cell; `o` and `t` identify the output and HTML table. These are rendered-output extracts and can inherit pandas display truncation. They do not replace full source tables or a fresh inference run.
- `claim_evidence_ledger.csv`: main claim-to-source mapping.
- `figure_01_central_findings.png`: live central category fractions and selected SFR-decomposition values.
- `figure_02_survivors_and_asymmetry.png`: saved Q3 and NGC4654 intervals; the two analyses are not combined statistically.
- `verification.json`: source hashes, links, PDF text/bounds, and render inventory from the final verification run.

Rebuild figures with the local research Python:

```sh
MPLCONFIGDIR=/private/tmp/mauve_second_paper_20261011/mpl /opt/miniconda3/envs/ICRAR/bin/python /Users/Igniz/Desktop/ICRAR/MAUVE/assets/20261011_paper_II_synthesis/plot_report_figures.py
```

Rebuild the PDF with the bundled PDF runtime:

```sh
/Users/Igniz/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 /Users/Igniz/Desktop/ICRAR/MAUVE/assets/20261011_paper_II_synthesis/build_pdf.py
```

The same bundled Python can run `verify_report.py`. For visual checks, first render the PDF using Poppler `pdftoppm -r 60 -png`, with output prefix `/private/tmp/mauve_second_paper_20261011/report_page`. The verification script creates contact sheets from those images; visual acceptance still requires inspecting them.

`check_central_results.py` is the isolated live science check. It requires the current local further FITS products and the ICRAR Python environment, runs only the five source cells listed above, and writes its results here. Rerunning it changes the evidence snapshot, so preserve this directory first if comparing a later reduction. The originally executed version wrote the same files to temporary storage before copying them here; only the output directory was adjusted for reuse.

`prepare_report_assets.py` archives the original inventory/assembly process and requires the temporary extracted evidence from this review. Use `plot_report_figures.py` for a standalone rebuild of the final figures; the assembly script contains an earlier figure layout.

The 7000 Angstrom products were not available locally for a corresponding science rerun on this date. That remains a prerequisite for the final manuscript claims.
