# Integrated pyqz Adaptive-KDE Recovery Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Recover numerically lost integrated pyqz LogQ and gas-phase O/H values for spectra whose central diagnostic ratios are inside the locked MAPPINGS-V grid, while retaining invalid results for genuinely off-grid spectra and leaving all spatial-bin results on the existing pyqz path.

**Architecture:** `run_pyqz_spectra` continues to call pyqz 0.8.4 first. An explicit `adaptive_kde_fallback=True` option may replay the same Monte Carlo sequence only after a non-finite native result or native exception, interpolate those samples through the same pyqz diagnostic, and evaluate the same multivariate KDE and 0.61 contour on a local adaptive mesh. `SFR+Z.py` opts in only for the two one-spectrum integrated HII/SF calls and logs whether recovery was used.

**Tech Stack:** Python 3.13 in the ICRAR Conda environment, NumPy, SciPy truncated-normal sampling, statsmodels `KDEMultivariate`, contourpy, pyqz 0.8.4 through `model_grid_compat.py`, pytest.

---

### Task 1: Lock the fallback contract with failing unit tests

**Files:**
- Modify: `tests/test_model_grid_diagnostics.py` near the existing pyqz batch-adapter test

- [x] **Step 1: Add a reusable synthetic pyqz double**

Add a fake whose native `get_global_qz` returns the observed non-finite combined-KDE shape, while `interp_qz` provides a deterministic on-grid transform and exposes the same metadata attributes used by the project fallback:

```python
class NarrowPdfFakePyqz:
    def __init__(self, *, central_on_grid=True):
        self.central_on_grid = central_on_grid
        self.interp_calls = 0
        self.pyqzm = SimpleNamespace(
            M_version="MV",
            PDF_cont_level=0.61,
            QZs_lim={
                "MV": {
                    "LogQ": np.array([6.5, 8.5]),
                    "gas[O]+12": np.array([8.0, 8.875]),
                }
            },
            diagnostics={
                "[NII]/[SII]+;[OIII]/[SII]+": {"coeffs": [[1, 0], [0, 1]]}
            },
        )

    def get_global_qz(self, data, data_cols, which_grids, **kwargs):
        del data, data_cols, which_grids, kwargs
        columns = [
            "<LogQ{KDE}>",
            "err(LogQ{KDE})",
            "<gas[O]+12{KDE}>",
            "err(gas[O]+12{KDE})",
            "flag",
            "rs_offgrid",
        ]
        return np.array([[np.nan, np.nan, np.nan, np.nan, 1234.0, 0.0]]), columns

    def interp_qz(self, qz, ratio_values, diagnostic, **kwargs):
        del diagnostic, kwargs
        self.interp_calls += 1
        first, second = (np.asarray(value, dtype=float) for value in ratio_values)
        if qz == "LogQ":
            output = 7.2 + 0.08 * first + 0.03 * second
        else:
            output = 8.5 + 0.02 * first - 0.04 * second
        if not self.central_on_grid:
            output = np.full_like(output, np.nan)
        return output
```

- [x] **Step 2: Add tests for the required behavior**

Add four focused tests:

```python
def test_pyqz_adaptive_fallback_recovers_narrow_on_grid_pdf():
    from model_grid_diagnostics import BinSpectra, run_pyqz_spectra

    spectra = BinSpectra(
        bin_ids=np.array([4]),
        fluxes=np.array([[10.0, 5.0, 28.6, 4.0, 3.0, 2.0]]),
        errors=np.full((1, 6), 0.01),
        pixel_counts=np.ones(1, dtype=int),
    )
    run = run_pyqz_spectra(
        NarrowPdfFakePyqz(), spectra, random_seed=17, srs=800,
        adaptive_kde_fallback=True,
    )

    assert run.results["valid"].tolist() == [1]
    assert np.isfinite(run.results["o_h"][0])
    assert np.isfinite(run.results["o_h_err"][0])
    assert np.isfinite(run.results["log_q"][0])
    assert np.isfinite(run.results["log_q_err"][0])
    assert run.failures == ()
    assert run.recoveries and run.recoveries[0][0] == 4


def test_pyqz_adaptive_fallback_is_deterministic_and_restores_numpy_state():
    from model_grid_diagnostics import BinSpectra, run_pyqz_spectra

    spectra = BinSpectra(
        bin_ids=np.array([4]),
        fluxes=np.array([[10.0, 5.0, 28.6, 4.0, 3.0, 2.0]]),
        errors=np.full((1, 6), 0.01),
        pixel_counts=np.ones(1, dtype=int),
    )
    state_before = np.random.get_state()
    first = run_pyqz_spectra(
        NarrowPdfFakePyqz(), spectra, random_seed=17, srs=800,
        adaptive_kde_fallback=True,
    )
    state_after = np.random.get_state()
    second = run_pyqz_spectra(
        NarrowPdfFakePyqz(), spectra, random_seed=17, srs=800,
        adaptive_kde_fallback=True,
    )

    for name in first.results:
        np.testing.assert_array_equal(first.results[name], second.results[name])
    for before, after in zip(state_before, state_after):
        np.testing.assert_array_equal(before, after)


def test_pyqz_adaptive_fallback_keeps_central_off_grid_spectrum_invalid():
    from model_grid_diagnostics import BinSpectra, run_pyqz_spectra

    spectra = BinSpectra(
        bin_ids=np.array([8]),
        fluxes=np.array([[10.0, 5.0, 28.6, 4.0, 3.0, 2.0]]),
        errors=np.full((1, 6), 0.01),
        pixel_counts=np.ones(1, dtype=int),
    )
    run = run_pyqz_spectra(
        NarrowPdfFakePyqz(central_on_grid=False), spectra,
        adaptive_kde_fallback=True,
    )

    assert run.results["valid"].tolist() == [0]
    assert run.recoveries == ()


def test_pyqz_adaptive_fallback_is_opt_in():
    from model_grid_diagnostics import BinSpectra, run_pyqz_spectra

    fake = NarrowPdfFakePyqz()
    spectra = BinSpectra(
        bin_ids=np.array([4]),
        fluxes=np.array([[10.0, 5.0, 28.6, 4.0, 3.0, 2.0]]),
        errors=np.full((1, 6), 0.01),
        pixel_counts=np.ones(1, dtype=int),
    )
    run = run_pyqz_spectra(fake, spectra)

    assert run.results["valid"].tolist() == [0]
    assert run.recoveries == ()
    assert fake.interp_calls == 0
```

- [x] **Step 3: Run the new tests and confirm the expected RED state**

Run:

```bash
/opt/miniconda3/envs/ICRAR/bin/python -m pytest \
  tests/test_model_grid_diagnostics.py -k 'pyqz_adaptive_fallback' -q
```

Expected: failures because `run_pyqz_spectra` has no `adaptive_kde_fallback` parameter and `ModelBatchRun` has no `recoveries` field.

### Task 2: Implement the fail-closed adaptive local KDE

**Files:**
- Modify: `model_grid_diagnostics.py:147-152`
- Modify: `model_grid_diagnostics.py:1104-1198`
- Test: `tests/test_model_grid_diagnostics.py`

- [x] **Step 1: Extend batch provenance without changing existing callers**

```python
@dataclass(frozen=True)
class ModelBatchRun:
    """One complete model pass, including isolated per-spectrum outcomes."""

    results: dict[str, np.ndarray]
    failures: tuple[tuple[int, str], ...]
    recoveries: tuple[tuple[int, str], ...] = ()
```

- [x] **Step 2: Add private helpers for identical MC replay and local KDE extraction**

The helpers must:

```python
PYQZ_ADAPTIVE_KDE_GRID_SIZE = 257
PYQZ_ADAPTIVE_KDE_PADDING_BW = 6.0


def _adaptive_pyqz_kde_result(
    pyqz,
    data,
    data_columns,
    *,
    diagnostic,
    qzs,
    pk,
    kappa,
    struct,
    sampling,
    error_pdf,
    srs,
    flag_level,
):
    """Replay pyqz MC samples and evaluate its KDE on a local mesh."""
```

Implementation requirements:

1. Reject configurations other than one diagnostic, two q/z fields, normal errors, positive `srs`, and the existing multivariate fallback call.
2. Reproduce pyqz's per-line `scipy.stats.truncnorm` sampling, including the central row and exact-zero-error behavior.
3. Build the requested log ratios and call `pyqz.interp_qz` with the same diagnostic coefficients, pressure, kappa, geometry, and interpolation sampling.
4. Return no recovery when either central q/z interpolation is non-finite.
5. Keep only jointly finite MC q/z pairs; compute `rs_offgrid = round(100 * rejected / srs, 1)`.
6. Use pyqz's bandwidth formula `1.06 * std * n**(-1/6)` for both axes and `statsmodels.nonparametric.KDEMultivariate(var_type="cc")`.
7. Construct 257-point local axes from sample extrema plus six bandwidths, clipped to `pyqz.pyqzm.QZs_lim[pyqz.pyqzm.M_version]`.
8. Evaluate and peak-normalise the density, extract `pyqz.pyqzm.PDF_cont_level` with contourpy, and require a contour containing the KDE peak.
9. Match pyqz's contour-centroid and maximum-half-extent summaries and its one-diagnostic flag comparison.
10. Return the existing `PYQZ_FLOAT_FIELDS`/`PYQZ_INTEGER_FIELDS` contract with `valid=1`; otherwise raise a descriptive exception so the caller can retain the native invalid record.

- [x] **Step 3: Replay only after a native failure and preserve the process RNG**

Extend the wrapper signature with `adaptive_kde_fallback: bool = False`. Immediately before the native call, capture the NumPy state. Immediately after it returns or raises, capture the advanced state. If and only if the native record is invalid and fallback is enabled, restore the pre-call state, run `_adaptive_pyqz_kde_result`, then restore the post-call state in `finally`. Record a successful recovery in `recoveries`; retain the native record and native exception in `failures` if recovery does not succeed.

- [x] **Step 4: Run the focused tests and existing pyqz adapter test**

Run:

```bash
/opt/miniconda3/envs/ICRAR/bin/python -m pytest \
  tests/test_model_grid_diagnostics.py \
  -k 'pyqz_adaptive_fallback or pyqz_batch_adapter' -q
```

Expected: all selected tests pass; the original adapter still reports `[1, 0, 1]` validity and restores NumPy state.

### Task 3: Wire recovery only into integrated HII/SF calculations

**Files:**
- Modify: `SFR+Z.py:328-337`
- Modify: `SFR+Z.py:3729-3745`
- Modify: `SFR+Z.py:4168-4240`
- Modify: `SFR+Z.py:4315-4340`
- Modify: `tests/test_model_grid_diagnostics.py:968-995`

- [x] **Step 1: Add a source-contract test before production wiring**

Extend `test_sfr_wires_fixed_model_grid_fields_and_provenance` with assertions that the source contains exactly two semantically relevant markers:

```python
assert "adaptive_kde_fallback=True" in source
assert "adaptive-local-KDE recovery" in source
```

Also assert the spatial `pyqz_bin_run` call does not contain the opt-in argument by inspecting the text between `pyqz_bin_run = run_pyqz_spectra(` and its closing `)`.

- [x] **Step 2: Confirm the source-contract test fails**

Run:

```bash
/opt/miniconda3/envs/ICRAR/bin/python -m pytest \
  tests/test_model_grid_diagnostics.py::test_sfr_wires_fixed_model_grid_fields_and_provenance -q
```

Expected: failure because integrated recovery has not yet been wired.

- [x] **Step 3: Enable and report the fallback only in the integrated call**

Pass `adaptive_kde_fallback=True` in the one-row integrated call only. Add `print_model_recoveries(label, recoveries)` alongside `print_model_failures`, store a per-region boolean such as `pyqz_recovered`, and append either `native-global-KDE` or `adaptive-local-KDE recovery` to each stored integrated summary line. Do not pass the option to `pyqz_bin_run`; do not alter FITS map fields or headers.

- [x] **Step 4: Update the `SFR+Z.py` changelog and README**

Document that the native pyqz call remains first, recovery is integrated-only, the same MAPPINGS-V grid/diagnostic/MC seed/sampling and 0.61 contour estimator are retained, and central off-grid spectra remain invalid. Preserve the user's existing JY22 2026-08-29 changelog entry.

- [x] **Step 5: Run the source-contract test**

Run the exact command from Step 2. Expected: pass.

### Task 4: Validate against all stored integrated apertures

**Files:**
- Modify: `.planning/integrated_pyqz_diagnosis_20260829/audit_all_integrated_centres.py`
- Modify: `.planning/integrated_pyqz_diagnosis_20260829/findings.md`
- Modify: `.planning/integrated_pyqz_diagnosis_20260829/progress.md`

- [x] **Step 1: Extend the read-only audit to call the opt-in wrapper**

For each reconstructed HII/SF one-row spectrum, run the updated wrapper with the production seed, `srs=800`, the locked diagnostic, and `adaptive_kde_fallback=True`. Emit galaxy, region, validity, central direct interpolation, recovered estimates/errors, raw/adaptive flag, off-grid fraction, and recovery provenance.

- [x] **Step 2: Run the all-aperture audit without writing products**

Run:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/private/tmp/mpl-pyqz-fix \
/opt/miniconda3/envs/ICRAR/bin/python \
  .planning/integrated_pyqz_diagnosis_20260829/audit_all_integrated_centres.py
```

Expected: 66/70 finite integrated pyqz results; only NGC4294 HII/SF and NGC4396 HII/SF remain invalid because their central ratios are outside the model grid.

- [x] **Step 3: Check deterministic repeat and estimator convergence**

Repeat the audit with the same seed and compare outputs exactly. On IC3392 SF and NGC4254 HII, compare 129-, 257-, and 513-point local meshes; require central estimates and contour extents at 257 versus 513 to agree within `5e-4` dex before retaining 257.

### Task 5: Full verification and review

**Files:**
- Review: `model_grid_diagnostics.py`
- Review: `SFR+Z.py`
- Review: `README.md`
- Review: `tests/test_model_grid_diagnostics.py`

- [x] **Step 1: Run syntax and complete focused test coverage**

```bash
/opt/miniconda3/envs/ICRAR/bin/python -m py_compile \
  model_grid_diagnostics.py model_grid_compat.py SFR+Z.py
/opt/miniconda3/envs/ICRAR/bin/python -m pytest \
  tests/test_model_grid_diagnostics.py -q
```

Expected: compilation succeeds and the complete test file passes.

- [x] **Step 2: Inspect only the task diff and preserve unrelated edits**

Use a zero-context diff for the touched blocks and confirm the existing JY22 empirical-QC changes are unmodified. Confirm no FITS products or stored SFR logs changed.

- [x] **Step 3: Request bounded read-only review**

Ask one reviewer to inspect the fallback, opt-in boundary, tests, and RNG handling. The primary agent verifies every finding against the live code before accepting it.

- [x] **Step 4: Record the audit job**

Append a concise entry to `/Users/Igniz/Desktop/Codex_log/2026_08_29.md` listing files changed, commands run, 66/70 recovery outcome, the four genuine off-grid apertures, and the fact that no production product was overwritten.

No commit is included: the working tree already contains unrelated user-owned JY22 changes, and the user requested a correct local fix rather than repository publication.
