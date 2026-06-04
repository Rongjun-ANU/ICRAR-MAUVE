# ICRAR And MAUVE Project Instructions For Rongjun

These instructions apply under `/Users/Igniz/Desktop/ICRAR`. Inherit the global
and Desktop-level instructions first.

## Core Project Rules

- Treat the live local tree as the source of truth for science outputs,
  notebooks, FITS products, logs, configs, and repo state.
- Stay in `quick_check` for vague small tasks and instruction-loading questions.
  For a request that only asks whether AGENTS would load, inspect only AGENTS
  files and do not read memory.
- Load deeper ICRAR, MAUVE, nGIST, or science-specific guidance only when the
  task actually concerns those files, products, workflows, or scientific rules.
- Use memory only when the user asks about prior context or current local files
  cannot determine the requested scope.
- Preserve the user's latest narrowed scope literally. If they narrow repos,
  galaxies, notebook cells, HDUs, lines, or plots, remove stale broader machinery
  instead of layering optional branches on top.
- Prefer `further` for active second-project work and use `extended` mainly as
  older first-project/reference material unless the user says otherwise.
- Do not run full science pipelines that overwrite products unless explicitly
  requested. Prefer temporary output under `/private/tmp` for behavior checks.

## Workspace Map

The main research root is `/Users/Igniz/Desktop/ICRAR`. It is a parent of nested
project repos and should not be treated as one simple monorepo.

Important areas:

- `/Users/Igniz/Desktop/ICRAR/further`: active MAUVE second-project scripts,
  notebooks, nGIST products, QC work, and current science exploration.
- `/Users/Igniz/Desktop/ICRAR/extended`: older/first-project MAUVE pipeline
  scripts and products; often the source for copied files.
- `/Users/Igniz/Desktop/ICRAR/MAUVE`: notes, papers, slides, manuscript context,
  and project planning.
- `/Users/Igniz/Desktop/ICRAR/nGIST_v7.6.8_setup`: Setonix config generation,
  scripts, status tooling, and README workflow docs.
- `/Users/Igniz/Desktop/ICRAR/nGIST_dev_7.6.4`: older nGIST development
  experiments.
- `/Users/Igniz/Desktop/ICRAR/v3tk_to_VRI`: VRI cube conversion and image
  products.
- `/Users/Igniz/Desktop/ICRAR/v3tk_masking_VRI`: VRI-based spatial-mask repo.
- `/Users/Igniz/Desktop/ICRAR/pPXF/ppxf_examples` and
  `/Users/Igniz/Desktop/ICRAR/pPXF/ppxf_examples_try`: nested git examples.

Always check the repo boundary before committing, diffing, staging, or cleaning.

## Python And Notebooks

- Prefer `/opt/miniconda3/envs/ICRAR/bin/python` for ICRAR Python checks.
- If Jupyter or `nbconvert` fails with `PermissionError: [Errno 1] Operation not
  permitted` from `find_available_port`, retry with writable temporary
  `MPLCONFIGDIR`, `IPYTHONDIR`, `JUPYTER_CONFIG_DIR`, and
  `JUPYTER_RUNTIME_DIR`; if port binding remains blocked, use direct
  in-process notebook execution against the real local files.
- For notebooks, inspect the live structure before editing: config cells, helper
  cells, plotting cells, markdown interpretation, and current cell numbering.
- Respect exact cell/layout requests. If the user says "after cell 19", insert
  after the actual current cell 19.
- Clear stale notebook outputs when code or interpretation changes and old
  outputs could mislead.

## `further` Rules

- Active products are typically nested under `v3tk_v7.6.8/<GALID>/...`; do not
  assume older flat `*_extended.fits` layouts in copied notebooks or scripts.
- Keep galaxy counts dynamic. Do not preserve fixed `<=14` assumptions when the
  live sample has changed.
- If EW or map shapes differ by an edge row/column, consider explicit overlap
  alignment and NaN padding instead of crashing or silently trimming science
  pixels.
- Distinguish QC-ready from science-ready samples. Science-ready usually requires
  `_further` mass/SFR/metallicity products and the requested HDUs.
- Handle combined cube IDs carefully, especially `NGC4567_8`; reuse existing
  local combined-ID logic.
- Product discovery should accept uppercase/lowercase gas-bin variants and
  `.fits.gz` where relevant.
- For `SFR+Z.py`, preserve Calzetti as the default dust law unless the user
  explicitly changes it.
- Use corrected fluxes for BPT and line-ratio diagnostics unless raw flux is
  explicitly requested.
- Keep non-detections distinct from unclassified when BPT semantics matter.
- For Te-metallicity work, separate exploratory/raw maps from paper-quality
  masks; broad low-S/N tails are not automatically formula bugs.
- For recent `SII/NII` notebook work, use linear `SII/NII` when requested and
  treat it as a mixed abundance plus DIG/shock-sensitive stress-test axis.
- When asked to comment out `savefig`, comment out only `savefig(...)` calls and
  multiline continuations, not the whole analysis cell.

## `extended` Rules

- `extended` is older source material for scripts later copied into `further`;
  reset from it only when that is the intended scope.
- For careful reviews, syntax checks are not enough. Check math, astrophysical
  assumptions, masks, calibration validity, uncertainty propagation, FITS
  headers, logs, and actual output behavior.
- For `Mass.py`, the `light_norm_to_r` branch is expected when the nGIST config
  uses `SFH ... NORM_TEMP: LIGHT`.
- BPT maps should use corrected flux ratios. Current semantics are
  `-1 = unknown/non-detection`, `0 = unclassified`, and positive values for
  stable classes.
- For C20 metallicity branches, do not collapse a combined result to the method
  with the smallest formal scatter; combine finite methods transparently.
- When changing scripts with a beginning changelog, update that changelog in the
  same pass.

## nGIST And Setonix

For `/Users/Igniz/Desktop/ICRAR/nGIST_v7.6.8_setup`:

- The active generator is usually `config_setup/make_gist_config_try.py`; treat
  older generators as provenance unless the user says otherwise.
- Read `config_setup/MAUVE_MasterConfig_v7.6.8_setonix.yaml` before stating
  current defaults.
- Generated Setonix products use `/scratch/pawsey1308/mauve/products/v3tk_v7.6.8/`.
- `CONFIG` beside each product `LOGFILE` is authoritative for enabled modules
  and should drive remaining-work estimates.
- For status or memory diagnosis, check per-galaxy `LOGFILE`s, not only summary
  tables.
- In Setonix `finished_efficiency` output, verify which column is peak memory;
  prior examples used `MaxRSS` as the second-last memory column.
- Treat queue splits and galaxy lists as snapshots. Re-check current scripts,
  logs, and scheduler output before reporting exact current membership.
- Official nGIST restart behavior for gas level `BOTH` is limited: completed BIN
  gas can be reused, but SPAXEL gas pPXF is not checkpointed per spaxel/chunk.
- When testing `27_status.sh` locally on macOS, account for BSD `date -j -f`
  versus GNU `date -d`.

## VRI Repositories

For `/Users/Igniz/Desktop/ICRAR/v3tk_to_VRI` and
`/Users/Igniz/Desktop/ICRAR/v3tk_masking_VRI`:

- Treat them as separate target repos and do not carry edits into adjacent VRI
  folders unless requested.
- Keep README files synchronized with behavior changes, especially input formats,
  output naming, and reproduction commands.
- Preserve direct `.fits.gz` support when modifying path helpers or wrappers.
- `v3tk_to_VRI` converts v3tk cubes into VRI products and downstream images.
- `v3tk_masking_VRI` builds nGIST-compatible spatial masks from VRI products and
  diagnostic overlays; mask convention is `0 = unmasked`, `1 = masked`.
- The mosaic arranger uses OR-Tools CP-SAT when available and can fall back to a
  heuristic mode. Preserve proof reports when layout optimality matters.
- Re-check current file counts before documenting sample size or release status.
- If products exceed GitHub's normal 100 MB limit, use Git LFS intentionally.

## Scientific Review Style

- For "is this really true?" questions, verify source code, notebook outputs or
  FITS products, and cited papers/methods together.
- Separate what the code does, what the data show, and what the literature
  supports.
- State applicability ranges, masking choices, uncertainty propagation, and
  assumptions explicitly.
- Translate caveats into concrete QC steps or paper-readiness criteria.
- For MAUVE next-project planning, prefer paper-shaped, physically interpretable
  100 pc local-ISM questions before high-cost extensions.

## Documentation And Verification

- If the user asks for a detailed README, update docs in the same pass as code
  behavior.
- Write reproducible docs: inputs, environment, exact commands, outputs,
  verification, and caveats.
- Avoid stale claims such as fixed galaxy counts, fixed queue lists, or exact
  package counts unless verified in the current task or clearly marked as a
  snapshot.
- Match verification to risk: syntax/import checks, one-galaxy smoke tests,
  one-HDU/header inspection, one representative notebook block, or safe dry-runs
  before larger reruns.
