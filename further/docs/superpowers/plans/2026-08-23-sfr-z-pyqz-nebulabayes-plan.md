# SFR+Z pyqz and NebulaBayes Integration Plan

**Goal:** Extend `further/SFR+Z.py` with two independent, model-grid metallicity diagnostics: pyqz gas-phase oxygen abundance and NebulaBayes HII-grid oxygen abundance. Also retain the physical parameters actually constrained by each model: pyqz `LogQ` and derived `log U`, and NebulaBayes `log U` and `log(P/k)`.

**Architecture:** Keep the existing Calzetti-corrected flux/error maps, BPT classifications, spatial BINID collapse, integrated-flux convention, legacy metallicities, PyNeb products, and output ordering intact. Put the testable model adapters in a small helper module, evaluate each unique spatial bin once, broadcast the result back to pixels, then make immediate HII/SF output pairs. Initialise the NebulaBayes grid once per process. Never infer either HII-region model outside the existing HII or SF selections.

**Execution note:** Implement task by task with scientific-review checkpoints. Stop at the pyqz compatibility gate if the port cannot reproduce the official 0.8.4 fixtures; do not emit an `O_H_PYQZ` product from an approximate replacement.

---

## 1. Locked scientific contract

### Shared inputs and selections

- Use the existing dust-corrected fluxes and propagated errors for `HB4861`, `OIII5006`, `HA6562`, `NII6583`, `SII6716`, and `SII6730`.
- Map the nGIST names to the model-grid names explicitly: `OIII5006 -> OIII5007` and `SII6730 -> SII6731`.
- Preserve `extinction_law = "calzetti"` as the default. NebulaBayes must receive already corrected measurements with `deredden=False`; it must not perform a second extinction correction.
- Evaluate only where the existing `mask_HII` or `mask_SF` is true and every required corrected flux and error is finite and strictly positive. This is an input-domain gate, not a new S/N threshold. Do not add an implicit `S/N >= 3` model cut.
- Use one common six-line validity aperture for both methods so their map coverage and integrated apertures are directly comparable.
- Run inference once on the union `(mask_HII | mask_SF) & model_input_valid`. Apply `mask_HII` and `mask_SF` only when producing the paired maps. A bin that belongs to both selections must have identical underlying inference values in both products.
- Leave all model outputs outside the requested selection as NaN, and retain an explicit validity map for downstream QC.
- Do not calculate an all-gas/total-region model-grid value. These are HII-region grids; only the HII and SF integrated measurements are in scope.

### pyqz 0.8.4.0

Use the final-conversation configuration literally:

```python
PYQZ_DIAGNOSTIC = "[NII]/[SII]+;[OIII]/[SII]+"
PYQZ_QZS = ["LogQ", "gas[O]+12"]
PYQZ_PK = 5.0
PYQZ_STRUCT = "pp"
PYQZ_KAPPA = np.inf
PYQZ_SAMPLING = 2
PYQZ_ERROR_PDF = "normal"
PYQZ_SRS = 800
PYQZ_KDE_METHOD = "multiv"
PYQZ_KDE_QZ_SAMPLING = 101j
PYQZ_KDE_DO_SINGLES = True
PYQZ_FLAG_LEVEL = 2.0
PYQZ_RANDOM_SEED = 20260823
```

- Feed `[NII] = NII6583`, `[OIII] = OIII5006`, and `[SII]+ = SII6716 + SII6730`; propagate the two [S II] component errors in quadrature. Use the explicit column sequence `["[NII]", "std[NII]", "[SII]+", "std[SII]+", "[OIII]", "std[OIII]"]`.
- Use the MAPPINGS V 5.0.16, plane-parallel, `log(P/k)=5.0`, `kappa=inf` grid bundled with pyqz. The pressure is a fixed model assumption, not a measured pressure product.
- Treat `8.00 <= gas[O]+12 <= 8.875` and `6.5 <= LogQ <= 8.5` as the package-level MAPPINGS V bounds; retain pyqz's diagnostic-specific off-grid behaviour inside those limits.
- Read the combined KDE columns `<gas[O]+12{KDE}>`, `err(gas[O]+12{KDE})`, `<LogQ{KDE}>`, and `err(LogQ{KDE})`. With one diagnostic, do not mislabel the direct across-diagnostic standard deviation as a flux-propagated uncertainty.
- Derive `log U = LogQ - log10(c / (cm s^-1)) = LogQ - 10.4771212547`; the numerical uncertainty is unchanged by this constant shift.
- Preserve the raw pyqz `flag` integer and `rs_offgrid` percentage. The pyqz uncertainty is the maximum half-extent of its peak-normalised 0.61 KDE contour around the contour-vertex mean; describe it that way in FITS comments rather than calling it a Gaussian 1-sigma error.
- Keep `KDE_method="multiv"` only if the compatibility test confirms that current statsmodels reproduces the archived result within tolerance. Never switch silently to `gauss`.
- Process bins in stable ascending BINID order with one deterministic random stream. Start with `nproc=1`; enable parallel execution only after repeated-run identity and seed partitioning are tested.

### NebulaBayes 0.9.9

Use the built-in MAPPINGS 5.1 HII grid with the exact six-line list:

```python
NB_LINE_LIST = [
    "Hbeta", "OIII5007", "Halpha",
    "NII6583", "SII6716", "SII6731",
]
NB_GRID_PARAMS = ["log U", "log P/k", "12 + log O/H"]
NB_INTERPD_GRID_SHAPE = [40, 20, 160]
NB_INTERP_ORDER = 1
NB_GRID_ERROR = 0.10
NB_NORM_LINE = "Hbeta"
NB_DEREDDEN = False
NB_PRIOR = "Uniform"
NB_LIKELIHOOD_LINES = None
```

- The installed grid spans `7.061 <= 12+log(O/H) <= 9.304`, `-4 <= log U <= -2`, and `4.2 <= log(P/k) <= 8.6`. Record those bounds and flag posterior peaks at grid edges.
- Initialise `NB_Model` once, with only the six required grid columns, and reuse it for every bin and integrated spectrum.
- Store the marginal-posterior peak (`Estimate`) and the package's equal-tailed `CI68_low` and `CI68_high` for all three parameters. Do not call the peak a median and do not turn an asymmetric interval into a symmetric error.
- Store the package's reduced chi-square, maximum number of local maxima across the three marginal PDFs, and an integer QC bitmask for edge/open-interval/exception states. Do not invent an arbitrary "broad posterior" cutoff; the CI limits preserve the posterior width directly.
- NebulaBayes divides errors by Hbeta but does not propagate uncertainty in the normalising Hbeta flux. Before the call, form `R_i = F_i/F_Hbeta` and `sigma(R_i) = R_i * sqrt((sigma_i/F_i)^2 + (sigma_Hbeta/F_Hbeta)^2)` under the independent-error approximation. Pass these normalised values so the API divides by a normalisation of one. State that covariance shared by the Hbeta denominator and by the common Balmer-decrement correction remains unmodelled because NebulaBayes accepts independent per-line errors only.
- Treat `grid_error=0.10` as a model-grid likelihood term. Note separately that NebulaBayes' reported best-model reduced chi-square is calculated from observational errors and does not include that grid-error term.

### Explicit non-goals

- Do not add JY20/JY22, IZI, HCm, or any third metallicity method.
- Do not add [O II], [O I], He I, or auroral-line constraints to these two fits.
- Do not modify BPT boundaries, HII/SF semantics, dust correction, existing metallicity equations, PyNeb calculations, or current integrated legacy values.
- Do not combine or average pyqz and NebulaBayes metallicities; they remain separate model-dependent scales.
- Do not write a pyqz pressure map. Its `log(P/k)=5` is fixed and belongs in provenance metadata only.
- Do not write a redundant linear `P/k` map; `LOG_PK_NEBULABAYES` is the native inferred parameter.

---

## 2. Current compatibility gate

The packages are installed in `/opt/miniconda3/envs/ICRAR`, but installation alone is not a runnable preflight:

- `pyqz==0.8.4.0` fails on import because its package uses Python-2 absolute imports; its computational source also contains Python-2 syntax, removed NumPy aliases, and the removed private `matplotlib._cntr` API.
- `NebulaBayes==0.9.9` fails on import because SciPy 1.15 no longer exports `cumtrapz` and `simps`. After narrow API aliases, its built-in shorthand loader also needs an Astropy-HDU compatibility workaround and float64 grid parameters.
- `statsmodels==0.14.6` is installed. A deterministic one-spectrum smoke test now reaches its `KDEMultivariate` path, but official pyqz fixture parity is still required before science-product integration.

Recommended implementation:

1. Do not edit `site-packages`.
2. Put the NebulaBayes compatibility bootstrap and explicit HDU-1 grid loader in the project helper module. Limit the shims to the exact tested version and fail fast on an unsupported version.
3. Keep the official pyqz installation unchanged. At process startup, verify `pyqz==0.8.4.0` and exact source hashes, copy the installed package and bundled grids to a disposable directory under the platform temporary root, apply the narrow Python-3 compatibility transformations there, compile-check the transformed files, and import that runtime copy. Port only package imports, Python-2 syntax/iterator behaviour, removed NumPy aliases, expected zero-dispersion warning handling, and contour extraction via `contourpy`; do not change the interpolation, Monte Carlo, KDE, flags, or grid values.
4. Fail closed if the version, hashes, patch-match counts, compile step, or runtime import differs from the tested installation. A pinned external legacy runner remains the fallback only if official pyqz parity tests reject the runtime overlay.

Compatibility status on 2026-08-23: `tests/test_model_grid_compat.py` passes a finite bundled-grid NebulaBayes HII fit, pyqz's upstream grid-node/off-grid checks, the archived pyqz example tolerance, and a deterministic `multiv` KDE repeat while confirming identical before/after hashes for the installed pyqz sources. This establishes package capability without editing `site-packages`; the production `srs=800` convergence/repeat gate and representative MAUVE-bin validation remain outstanding before `SFR+Z.py` integration.

The pyqz product is authorised only after all parity tests in Task 1 pass.

---

## 3. File map

- Modify: `further/SFR+Z.py` - configuration, orchestration, HII/SF maps, integrated summaries, FITS metadata, and changelog.
- Create: `further/model_grid_compat.py` - process-local, exact-version compatibility bootstrap for NebulaBayes and pyqz; no installed-package writes.
- Create: `further/model_grid_diagnostics.py` - pure input preparation, version-gated adapters, result extraction, bin broadcasting, flags, and integrated helpers.
- Create: `further/tests/test_model_grid_compat.py` - package-integrity, import, bundled-grid, end-to-end smoke, and deterministic-repeat tests.
- Create: `further/tests/test_model_grid_diagnostics.py` - fast compatibility, unit, masking, uncertainty, and schema-contract tests.
- Modify: `further/README.md` - dependencies, model assumptions, run behaviour, new products, and QC interpretation.
- Modify during execution: this plan - check off steps and record any approved deviation.

No FITS data product is edited in place during development; representative runs use a mirrored tree under `/private/tmp`.

Commands below assume the repository root `/Users/Igniz/Desktop/ICRAR` as the working directory.

---

## 4. Output contract

Use `{REGION}` below for immediate `HII` then `SF` pairs in the ordered output builder.

### pyqz HDUs

- `O_H_PYQZ_{REGION}` and `O_H_PYQZ_{REGION}_ERR`
- `LOG_Q_PYQZ_{REGION}` and `LOG_Q_PYQZ_{REGION}_ERR`
- `LOG_U_PYQZ_{REGION}` and `LOG_U_PYQZ_{REGION}_ERR`
- `PYQZ_FLAG_{REGION}` - raw pyqz flag; `-99` means not evaluated
- `PYQZ_RS_OFFGRID_{REGION}` - percentage of Monte Carlo samples outside the grid
- `PYQZ_VALID_{REGION}` - 1 only when the input and finite KDE O/H/LogQ result are valid

### NebulaBayes HDUs

- `O_H_NEBULABAYES_{REGION}`, `O_H_NEBULABAYES_{REGION}_CI68_LOW`, and `O_H_NEBULABAYES_{REGION}_CI68_HIGH`
- `LOG_U_NEBULABAYES_{REGION}`, `LOG_U_NEBULABAYES_{REGION}_CI68_LOW`, and `LOG_U_NEBULABAYES_{REGION}_CI68_HIGH`
- `LOG_PK_NEBULABAYES_{REGION}`, `LOG_PK_NEBULABAYES_{REGION}_CI68_LOW`, and `LOG_PK_NEBULABAYES_{REGION}_CI68_HIGH`
- `NB_CHI2_RED_{REGION}` - package-reported observational-error reduced chi-square
- `NB_NLOCALMAX_{REGION}` - maximum marginal-PDF local-maximum count across O/H, log U, and log(P/k)
- `NB_FLAG_{REGION}` - `-99` not evaluated; otherwise value 1 O/H edge, value 2 log U edge, value 4 log(P/k) edge, value 8 any non-finite/open 68% interval, and value 16 fit exception/zero posterior, combined as a bitmask
- `NB_VALID_{REGION}` - 1 only when all three estimates are finite

All abundance and logarithmic-parameter HDUs use dimensionless `BUNIT=1` plus explicit comments. `LOG_PK` comments must define `log10(P/k / (K cm^-3))`; `PYQZ_RS_OFFGRID` uses `%`; flags and validity maps use integer types.

Add primary-header provenance with compact FITS-safe keys for package version, grid version, diagnostic, pressure, geometry, kappa, sampling, KDE, random seed, NebulaBayes shape/order/grid error/prior, and `deredden=False`. Add `HISTORY` records for the compatibility layer and the unmodelled inter-line covariance. Extend the output helper to accept method-specific references instead of reusing the current PyNeb-only `refs=True` header fields.

Insert these products in `build_ordered_output_hdul()` after the strong-line/C20 block and before the PyNeb block. Do not add them to the older append block that is discarded when the ordered HDU list is rebuilt.

Integrated results remain log output only, matching the current script convention. Print, for HII and SF separately, the common valid pixel/bin count, pyqz O/H/LogQ/log U with KDE uncertainties and flags, and NebulaBayes O/H/log U/log(P/k) with 68% intervals, reduced chi-square, and QC flags.

---

## 5. Implementation tasks

### Task 1: Establish executable, numerically checked adapters

**Files:** create `model_grid_compat.py`, `model_grid_diagnostics.py`, `tests/test_model_grid_compat.py`, and `tests/test_model_grid_diagnostics.py`.

- [x] Add failing-then-passing tests that assert exact installed versions/hashes and exercise both compatibility routes without changing installed files.
- [x] Implement the version-gated NebulaBayes bootstrap: alias removed SciPy integration functions to their current equivalents, supply equivalent NumPy names required by 0.9.9, read `NB_HII_grid.fits.gz` from HDU 1, cast the three parameter columns to float64, and initialise the model with the six-line list.
- [x] Prove the NebulaBayes adapter with a bundled HII-grid example: finite estimates, expected parameter names/order, finite 68% intervals, and finite reduced chi-square.
- [x] Implement an exact-hash-guarded, disposable pyqz 0.8.4 runtime overlay for the `get_global_qz` path. Replace `_cntr` with `contourpy` without storing a modified upstream package in the repository or editing `site-packages`.
- [x] Run the upstream grid-node and off-grid tests. Recover grid-node LogQ/O/H at the upstream rounding criteria, return NaN outside the interpolation grid, and preserve `get_global_qz` off-grid flag 8.
- [x] Run the archived `example_input.csv` values with a fixed seed. Direct outputs match the archived five-decimal values; `multiv` KDE centres/errors agree within 0.02 dex; raw flag and off-grid percentage remain zero.
- [x] Run an identical one-spectrum pyqz smoke call twice with `srs=80` and require bitwise-identical outputs. Repeat at the production `srs=800` setting with the official fixture before `SFR+Z.py` integration.

Verification command:

```bash
/opt/miniconda3/envs/ICRAR/bin/python -m pytest \
  further/tests/test_model_grid_compat.py \
  further/tests/test_model_grid_diagnostics.py -q
```

### Task 2: Implement pure model-input and result helpers

**Files:** modify `model_grid_diagnostics.py` and its tests.

- [x] Implement the common finite-positive six-line validity mask without a hard S/N threshold.
- [x] Implement [S II] summed flux and quadrature error construction.
- [x] Implement Hbeta-normalised NebulaBayes inputs with denominator-error propagation and tests for the analytic formula.
- [x] Implement explicit result extractors for pyqz KDE columns and NebulaBayes `DF_estimates`; fail on missing or renamed fields.
- [x] Implement `LogQ -> log U` with `log10(c_cm_s)` and a regression test for the constant `10.4771212547`.
- [x] Implement the NebulaBayes QC bitmask, local-maximum summary, pyqz raw flag preservation, and `-99` not-evaluated sentinel.
- [x] Implement stable unique-BINID extraction and broadcasting. Assert that all required binned flux/error values are constant within each bin before selecting a representative pixel.

### Task 3: Integrate the map calculations into SFR+Z.py

**Files:** modify `SFR+Z.py` and tests.

- [x] Add a dated changelog entry and visible method configuration beside `ENABLE_TE_METALLICITY_PRODUCTS`.
- [x] Add fail-fast imports only when each model product is enabled. Report interpreter, package versions, grid identity, and settings to the redirected log.
- [x] Insert the model-grid section immediately after the complete HII/SF corrected-line views are created and before the PyNeb block.
- [x] Build the common validity mask, collect stable ascending unique BINIDs from the HII/SF union, run pyqz and NebulaBayes once per bin, and broadcast results.
- [x] Create every HII/SF map by applying the existing masks to the common underlying result arrays. Preserve NaN outside region masks and `-99` in integer QC maps.
- [x] Catch and count per-bin model exceptions without aborting the other bins, but print the BINID and exception summary. Abort the run if every eligible bin fails or if the model/grid preflight fails.
- [x] Add progress and elapsed-time reporting. Do not lower `srs`, interpolation shape, or grid-error settings in response to runtime without a documented convergence decision.

### Task 4: Add HII- and SF-integrated model inference

**Files:** modify `SFR+Z.py`, `model_grid_diagnostics.py`, and tests.

- [x] Build one common six-line raw-flux aperture for each of `mask_HII` and `mask_SF`.
- [x] Sum raw fluxes over exactly that shared aperture, sum raw errors in quadrature, derive one integrated Balmer decrement with the existing 2.86 floor, and apply one integrated Calzetti correction using the existing helpers.
- [x] Feed the corrected integrated spectrum and errors to each model exactly once. Do not average spaxel/bin metallicities to obtain an integrated value.
- [x] Print values, uncertainty/CI semantics, valid pixel/bin counts, and QC. Do not add integrated scalar HDUs and do not run an all-region fit.
- [x] Leave the existing legacy integrated calculations unchanged; the new common-aperture helper is model-specific so old logged values do not drift.

### Task 5: Extend the ordered FITS schema and provenance

**Files:** modify `SFR+Z.py` and tests.

- [x] Add an output-contract test for every HDU name, HII/SF adjacency, dtype, BUNIT, shape, WCS, and comment.
- [x] Append the pyqz and NebulaBayes groups only in `build_ordered_output_hdul()` after C20 and before PyNeb.
- [x] Add package/grid/settings metadata to the copied primary header and method-specific references to the new image headers.
- [x] Assert unique EXTNAME values and that the original gas HDUs remain first and unchanged.
- [x] Keep current products in their present relative order and preserve all legacy data arrays when model products are disabled.

### Task 6: Document the reproducible workflow

**Files:** modify `README.md` and the `SFR+Z.py` changelog.

- [x] Document the six input lines, nGIST-to-grid wavelength names, Calzetti-first flow, exact HII/SF selections, no-extra-S/N rule, and integrated-aperture rule.
- [x] Document all pyqz/NebulaBayes settings, grid ranges, fixed versus inferred quantities, package versions, compatibility approach, random seed, and uncertainty semantics.
- [x] Add the new HDU table and QC flag definitions.
- [x] State the independent-error/covariance limitation and that cross-method offsets are expected because the grids and inference schemes differ.
- [x] Add the normal run command and a safe representative smoke-test recipe.

### Task 7: Scientific convergence, regression, and visual QA

**Files:** tests plus temporary outputs under `/private/tmp`; no committed science product.

- [x] Reconstruct a raw representative input from the original-HDU prefix of a current `_further` product, mirror its spatial-binning and mask files under `/private/tmp`, and run there so no live product is overwritten. Execution used a four-BINID IC3392 cutout because it is the smallest current product while exercising the same schema.
- [x] With both new methods disabled, compare every legacy output HDU against the pre-change baseline using exact equality where possible and `allclose(equal_nan=True)` for floating products; explain any header-only provenance differences.
- [x] With methods enabled, verify shape/WCS, HII/SF masking, common coverage, constant values within each BINID, identical overlap values, valid-map consistency, grid-bound compliance, and correct not-evaluated sentinels.
- [ ] For pyqz, compare `srs=400`, `800`, and `1600` on a representative low/high-S/N and near-grid-edge subset. Accept 800 only if 800-to-1600 centre shifts are below 0.01 dex and uncertainty changes are below 10%; otherwise increase the production setting explicitly.
- [ ] For NebulaBayes, compare `[40,20,160]` with `[60,30,240]`, and `[80,40,320]` if needed. Require median estimate shifts below 0.01 dex, CI-bound shifts below 0.02 dex or 10% of the baseline interval width, and no systematic edge-flag changes.
- [x] Benchmark 10 and 100 unique bins. Optimise batching or deterministic parallelism if necessary; do not weaken the scientific settings silently.
- [ ] Produce temporary PNG/PDF QC panels: pyqz O/H, NebulaBayes O/H, their difference, both log U maps and their difference, NebulaBayes log(P/k), pyqz off-grid fraction/flag, and NebulaBayes chi-square/flag. Inspect coverage and grid-edge pile-ups; do not use cross-method agreement as a pass/fail criterion.
- [x] Run the complete target tests and a syntax check:

```bash
/opt/miniconda3/envs/ICRAR/bin/python -m py_compile \
  further/SFR+Z.py further/model_grid_compat.py further/model_grid_diagnostics.py
/opt/miniconda3/envs/ICRAR/bin/python -m pytest \
  further/tests/test_model_grid_compat.py \
  further/tests/test_model_grid_diagnostics.py -q
```

- [x] Inspect the final diff, confirm no generated FITS/plots/logs are staged, and append the required Codex job log with checks, results, runtime, compatibility route, and residual risks.

Execution notes (2026-08-23):

- The four-bin IC3392 smoke run completed twice with bitwise-identical model
  arrays. It wrote 44 model HDUs between `COMBINED_C20_METHOD` and `NE_SII`,
  passed FITS verification, preserved all 118 original image arrays exactly,
  and produced finite HII/SF map and integrated results for both methods.
- The disabled-method regression contained no model HDUs and matched every
  enabled-run legacy HDU name and data array exactly after removing the new
  model block. Only the intentional primary-header model provenance differed.
- The production `srs=800` pyqz setting is deterministic. On four adjacent
  high-metallicity bins, 800-to-1600 O/H shifts were at most 0.00384 dex, while
  LogQ shifts reached 0.0110 dex and individual uncertainty changes reached
  about 12 percent. The locked final-conversation setting remains 800; the
  stricter convergence checkbox stays open rather than silently changing it.
- NebulaBayes O/H was stable across `[40,20,160]`, `[60,30,240]`, and
  `[80,40,320]`, but the discrete marginal peaks for log U and especially
  log(P/k) moved by more than the aspirational threshold as grid resolution
  changed. The 68-percent intervals and edge-flag pattern were much more
  stable. The locked `[40,20,160]` setting remains unchanged and this residual
  discretisation sensitivity is explicitly retained for follow-up.
- IC3392 benchmarks were 1.93/18.62 seconds for pyqz and 0.21/1.70 seconds for
  NebulaBayes over 10/100 bins. pyqz returned 9/10 and 99/100 finite KDE results;
  the non-finite spectrum was retained through the native QC fields rather than
  treated as a Python exception.
- Full-galaxy visual QA remains open. A four-pixel cutout is not scientifically
  useful for map panels, and no live full-galaxy product was overwritten merely
  to manufacture a figure during implementation.

---

## 6. Acceptance criteria

The implementation is complete only when:

1. The official pyqz compatibility fixtures and deterministic-repeat tests pass; otherwise no pyqz-labelled HDUs are written.
2. Both methods use the existing corrected fluxes and HII/SF masks, with no hidden S/N or classification change.
3. Every requested map and integrated HII/SF value has the correct uncertainty/credible-interval output and QC metadata.
4. pyqz pressure is visibly identified as fixed, while NebulaBayes pressure is stored as inferred `log(P/k)`.
5. Existing HDUs and legacy numerical products are unchanged apart from explicit new provenance/header additions.
6. The representative one-galaxy smoke run, convergence tests, FITS-schema tests, and visual QC all pass.
7. README and script changelog describe the executable settings, provenance, limitations, and output contract exactly.

## Sources fixed for implementation

- [pyqz 0.8.4 documentation](https://fpavogt.github.io/pyqz/)
- [pyqz installation and statsmodels KDE guidance](https://fpavogt.github.io/pyqz/installation.html)
- [pyqz 0.8.4 PyPI release](https://pypi.org/project/pyqz/)
- [Official pyqz repository](https://github.com/fpavogt/pyqz)
- [Official NebulaBayes repository and HII-grid description](https://github.com/ADThomas-astro/NebulaBayes)
