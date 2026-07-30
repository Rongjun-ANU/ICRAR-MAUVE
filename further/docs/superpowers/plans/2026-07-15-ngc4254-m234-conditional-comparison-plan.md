# NGC4254 Conditional m=2,3,4 Spiral Comparison Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Convert the NGC4254 notebook from a seed-sensitive single-arm-number selector into a leakage-controlled, negative-winding comparison of fixed `m=2`, `m=3`, and `m=4` ridge geometries and their conditional KTZ-compatible source-profile fits.

**Architecture:** Retain raw log-polar numerator and coverage arrays, split train/test azimuth sectors before smoothing, and use a four-bin guard plus the phase-space narrow-minus-broad response. Search negative pitch independently within each fixed arm number, attach 32-draw blocked null diagnostics without selecting among arm numbers, then fit three geometry-frozen KTZ profiles and present comparative skeleton, model, residual, profile, and parameter outputs.

**Tech Stack:** Python 3.13 in the ICRAR conda environment, Jupyter/nbformat/nbclient, NumPy, pandas, SciPy, Astropy FITS/WCS/coordinates, Matplotlib, and `unittest`.

---

## File structure and constraints

**Modify only:**

- `further/20260713_NGC4254_KTZ_spiral_fit.ipynb` — the complete scientific workflow and saved outputs.
- `further/tests/test_20260713_ngc4254_ktz_spiral_fit.py` — source, behavioral, and saved-result contracts.
- `further/docs/superpowers/plans/2026-07-15-ngc4254-m234-conditional-comparison-plan.md` — final execution record only in Task 7.

**Do not touch or stage:** any other modified or untracked workspace path.  The worktree contains unrelated user-owned changes.

**Python:** always use `/opt/miniconda3/envs/ICRAR/bin/python`.

**Notebook editing:** inspect live cell IDs first.  Use a temporary `nbformat` revision script under `/private/tmp`, created with `apply_patch`.  Preserve all existing cell IDs and order.  Clear outputs and execution counts from the first changed code cell onward until the final full execution.

**Reference design:** `further/docs/superpowers/specs/2026-07-15-ngc4254-m234-conditional-comparison-design.md`.

**Known live invariants:**

- finite HII pixels: `681856`;
- independent HII gas bins: `32552`;
- minimum positive morphology radius: about `0.00664195 kpc`;
- Brown geometry: centre `12h18m49.68s +14d25m05.52s`, inclination `39 deg`, directed PA `243 deg`;
- upstream equation-3 `b/a`: about `0.787` and already applied.

### Notebook revision-script pattern

Every task-specific revision script must follow this structure:

```python
from pathlib import Path
import nbformat

ROOT = Path("/Users/Igniz/Desktop/ICRAR/further")
NOTEBOOK = ROOT / "20260713_NGC4254_KTZ_spiral_fit.ipynb"


def find_cell(nb, cell_id):
    matches = [cell for cell in nb.cells if cell.id == cell_id]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one cell id={cell_id!r}; found {len(matches)}")
    return matches[0]


def clear_from(nb, first_index):
    for cell in nb.cells[first_index:]:
        if cell.cell_type == "code":
            cell.execution_count = None
            cell.outputs = []


nb = nbformat.read(NOTEBOOK, as_version=4)
# Task-specific source replacements go here.
clear_from(nb, FIRST_CHANGED_INDEX)
nbformat.validate(nb)
nbformat.write(nb, NOTEBOOK)
```

---

### Task 1: Replace the single-winner contracts with comparative m=2,3,4 contracts

**Files:**

- Modify: `further/tests/test_20260713_ngc4254_ktz_spiral_fit.py`
- Modify: `further/20260713_NGC4254_KTZ_spiral_fit.ipynb`

- [ ] **Step 1: Add code-source extraction helpers to the test module**

Add imports and helpers:

```python
import ast

import numpy as np
from scipy.ndimage import gaussian_filter, gaussian_filter1d


def notebook_code_source():
    return "\n\n".join(
        "".join(cell["source"])
        for cell in load_notebook()["cells"]
        if cell["cell_type"] == "code"
    )


def extracted_function_namespace(function_names):
    tree = ast.parse(notebook_code_source())
    wanted = set(function_names)
    selected = [
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in wanted
    ]
    found = {node.name for node in selected}
    if found != wanted:
        raise AssertionError(f"Missing functions: {sorted(wanted - found)}")
    namespace = {
        "np": np,
        "gaussian_filter": gaussian_filter,
        "gaussian_filter1d": gaussian_filter1d,
        "LOGPOLAR_SMOOTH_SIGMA": (0.8, 1.2),
        "LOGPOLAR_AZIMUTH_BROAD_SIGMA_BINS": 30.0,
        "RIDGE_GUARD_BINS": 4,
        "RIDGE_N_SECTORS": 12,
    }
    module = ast.Module(body=selected, type_ignores=[])
    exec(compile(module, "notebook-functions", "exec"), namespace)
    return namespace
```

- [ ] **Step 2: Replace stale source contracts with the approved comparison contracts**

Add these methods and remove the obsolete tests that require `m=1..6`, both signs, a single `ridge_geometry`, one deterministic eight-draw stream, or a single KTZ fit:

```python
def test_conditional_m234_negative_winding_contract(self):
    text = notebook_source()
    for required in [
        "M_COMPARE = np.array([2, 3, 4]",
        "RIDGE_PITCH_GRID_DEG = np.arange(-45.0, -4.9, 1.0)",
        "conditional_ridge_search",
        "ridge_geometry_table",
        "leakage_controlled_negative_winding_ridge",
        "RIDGE_M234_REAL_COMPLETE",
    ]:
        self.assertIn(required, text)
    self.assertNotIn("np.arange(1, 7", text)
    self.assertNotIn("np.arange(5.0, 45.1", text)


def test_leakage_controlled_fold_contract(self):
    text = notebook_source()
    for required in [
        "raw_weighted_sum",
        "raw_coverage",
        "preprocess_logpolar_raw",
        "build_sector_fold_maps",
        "RIDGE_GUARD_BINS = 4",
        "circular_dilate",
        "circular_erode",
        'field="radial_residual"',
        "RIDGE_LEAKAGE_GUARD_PASS",
        "histogram_cache_payload_bytes = 0",
    ]:
        self.assertIn(required, text)
    self.assertNotIn('field="local_ridge"', text[text.index("def conditional_ridge_search"):])
    self.assertNotIn("histogram_cache = {}", text)


def test_blocked_nulls_are_descriptive_not_selective(self):
    text = notebook_source()
    for required in [
        "RIDGE_NULL_BLOCK_SEEDS = np.array([4254, 5254, 6254, 7254]",
        "RIDGE_NULL_DRAWS_PER_BLOCK = 8",
        "build_blocked_null_diagnostics",
        "null_block_summary",
        "pooled_null_summary",
        "RIDGE_NULL_BLOCKS_COMPLETE",
        "descriptive_not_model_selection",
    ]:
        self.assertIn(required, text)
    self.assertNotIn("accepted.iloc[0].to_dict()", text)


def test_three_fixed_geometry_ktz_fits_are_declared(self):
    text = notebook_source()
    for required in [
        "fit_ktz_profile_fixed_geometry",
        "conditional_fits",
        "conditional_fit_tables",
        "profile_maxima_by_m",
        "KTZ_M234_PROFILE_FITS_COMPLETE",
        "HARMONIC_M234_MAXIMA_COMPLETE",
        '"geometry_source": geometry["geometry_source"]',
    ]:
        self.assertIn(required, text)


def test_m234_comparative_figures_and_table_are_declared(self):
    text = notebook_source()
    for required in [
        "plt.subplots(2, 2",
        "observed sky-plane",
        "observed deprojected",
        "conditional KTZ harmonic profiles",
        "ridge-width and descriptive-null comparison",
        "plt.subplots(3, 2",
        "conditional source model",
        "observed - model residual",
        "M234_COMPARISON_LAYOUT_COMPLETE",
        "M234_MODEL_RESIDUAL_LAYOUT_COMPLETE",
        "FINAL_M234_PARAMETER_TABLE_COMPLETE",
    ]:
        self.assertIn(required, text)


def test_adversarial_sector_signal_has_zero_training_field(self):
    namespace = extracted_function_namespace([
        "preprocess_logpolar_raw",
        "sector_columns",
        "circular_dilate",
        "circular_erode",
        "logpolar_fold_masks",
        "run_adversarial_leakage_check",
    ])
    result = namespace["run_adversarial_leakage_check"](
        n_u=12, n_phi=48, n_sectors=4, held_sector=1, guard_bins=2)
    self.assertEqual(result["raw_training_signal_l1"], 0.0)
    self.assertLessEqual(result["max_abs_training_residual"], 1e-14)
```

- [ ] **Step 3: Replace the saved-result contract markers**

Keep the geometry/data provenance markers and require these revised workflow markers:

```python
for marker in [
    "BROWN_GEOMETRY_PASS",
    "BA_FACTOR_PASS",
    "UPSTREAM_BA_PASS",
    "WCS_REFERENCE_NOT_CENTER_PASS",
    "FIT_DOMAIN_ALL_VALID_HII_PASS",
    "RIDGE_PIXEL_DOMAIN_PASS",
    "RIDGE_LEAKAGE_GUARD_PASS",
    "RIDGE_M234_SYNTHETIC_PASS",
    "RIDGE_BRANCH_NORMALIZATION_PASS",
    "RIDGE_M234_REAL_COMPLETE",
    "RIDGE_NULL_BLOCKS_COMPLETE",
    "RIDGE_WIDTH_M234_COMPLETE",
    "RIDGE_SECTOR_M234_COMPLETE",
    "KTZ_M234_PROFILE_FITS_COMPLETE",
    "HARMONIC_M234_MAXIMA_COMPLETE",
    "DEPROJECTION_ROUNDTRIP_PASS",
    "SYNTHETIC_WINDING_SIGN_PASS",
    "M234_COMPARISON_LAYOUT_COMPLETE",
    "M234_MODEL_RESIDUAL_LAYOUT_COMPLETE",
    "FINAL_M234_PARAMETER_TABLE_COMPLETE",
    "POSITIVITY_CHECK_PASS",
]:
    self.assertIn(marker, joined)
self.assertGreaterEqual(image_count, 9)
self.assertNotIn("RuntimeWarning", joined)
```

- [ ] **Step 4: Run the new tests and verify RED**

Run:

```bash
/opt/miniconda3/envs/ICRAR/bin/python -m unittest \
  tests.test_20260713_ngc4254_ktz_spiral_fit.NotebookContractTests.test_conditional_m234_negative_winding_contract \
  tests.test_20260713_ngc4254_ktz_spiral_fit.NotebookContractTests.test_leakage_controlled_fold_contract \
  tests.test_20260713_ngc4254_ktz_spiral_fit.NotebookContractTests.test_blocked_nulls_are_descriptive_not_selective \
  tests.test_20260713_ngc4254_ktz_spiral_fit.NotebookContractTests.test_three_fixed_geometry_ktz_fits_are_declared \
  tests.test_20260713_ngc4254_ktz_spiral_fit.NotebookContractTests.test_m234_comparative_figures_and_table_are_declared \
  tests.test_20260713_ngc4254_ktz_spiral_fit.NotebookContractTests.test_adversarial_sector_signal_has_zero_training_field -v
```

Expected: six failures because the live notebook still contains the prior single-winner implementation and placeholders.

- [ ] **Step 5: Update visible configuration and top-level Markdown**

In cell `212df4a7`, replace the old arm/sign/null declarations with:

```python
M_COMPARE = np.array([2, 3, 4], dtype=int)
RIDGE_PITCH_GRID_DEG = np.arange(-45.0, -4.9, 1.0)
RIDGE_GUARD_BINS = 4
RIDGE_NULL_BLOCK_SEEDS = np.array([4254, 5254, 6254, 7254], dtype=int)
RIDGE_NULL_DRAWS_PER_BLOCK = 8
RIDGE_NULL_TOTAL_DRAWS = int(
    len(RIDGE_NULL_BLOCK_SEEDS) * RIDGE_NULL_DRAWS_PER_BLOCK)
```

Remove `M_CANDIDATES`, the positive pitch grid, and `RIDGE_N_NULL`.  Update cells `6bbac1c3`, `3de897d1`, and `23f5d07e` to state that the notebook compares exactly `m=2,3,4`, imposes negative winding as a visual morphology prior, and does not declare an arm-number winner.

- [ ] **Step 6: Verify the configuration contract and commit**

Run the first source contract and static AST parsing.  Stage exactly the notebook and test file, check the staged diff, and commit:

```bash
git -C /Users/Igniz/Desktop/ICRAR add -- \
  further/20260713_NGC4254_KTZ_spiral_fit.ipynb \
  further/tests/test_20260713_ngc4254_ktz_spiral_fit.py
git -C /Users/Igniz/Desktop/ICRAR diff --cached --check
git -C /Users/Igniz/Desktop/ICRAR commit -m "test: define conditional NGC4254 m234 workflow"
```

---

### Task 2: Retain raw log-polar state and build leakage-free folds

**Files:**

- Modify: `further/20260713_NGC4254_KTZ_spiral_fit.ipynb`
- Test: `further/tests/test_20260713_ngc4254_ktz_spiral_fit.py`

- [ ] **Step 1: Replace the log-polar builder in cell `0d0ce897`**

Keep the radial exponential fitter and its diagnostics.  Replace `build_log_polar_contrast` with these complete helpers:

```python
def preprocess_logpolar_raw(
        raw_state, allowed_phi=None, include_local=True,
        smooth_sigma=LOGPOLAR_SMOOTH_SIGMA,
        broad_sigma_phi=LOGPOLAR_AZIMUTH_BROAD_SIGMA_BINS):
    weighted_sum = raw_state["raw_weighted_sum"].copy()
    coverage = raw_state["raw_coverage"].copy()
    if allowed_phi is not None:
        allowed_phi = np.asarray(allowed_phi, dtype=bool)
        weighted_sum[:, ~allowed_phi] = 0.0
        coverage[:, ~allowed_phi] = 0.0
    numerator = gaussian_filter(
        weighted_sum, smooth_sigma, mode=("nearest", "wrap"))
    denominator = gaussian_filter(
        coverage, smooth_sigma, mode=("nearest", "wrap"))
    valid = denominator > 0.05
    radial_residual = np.full_like(numerator, np.nan, dtype=float)
    radial_residual[valid] = numerator[valid] / denominator[valid]

    local_ridge = np.full_like(radial_residual, np.nan)
    broad_azimuthal = np.full_like(radial_residual, np.nan)
    if include_local:
        broad_numerator = gaussian_filter1d(
            np.where(valid, radial_residual * denominator, 0.0),
            broad_sigma_phi, axis=1, mode="wrap")
        broad_denominator = gaussian_filter1d(
            np.where(valid, denominator, 0.0),
            broad_sigma_phi, axis=1, mode="wrap")
        broad_azimuthal = np.divide(
            broad_numerator, broad_denominator,
            out=np.zeros_like(broad_numerator), where=broad_denominator > 0)
        local_ridge[valid] = radial_residual[valid] - broad_azimuthal[valid]

    if allowed_phi is not None:
        valid[:, ~allowed_phi] = False
        denominator[:, ~allowed_phi] = 0.0
        radial_residual[:, ~allowed_phi] = np.nan
        local_ridge[:, ~allowed_phi] = np.nan
        broad_azimuthal[:, ~allowed_phi] = np.nan
    return {
        "radial_residual": radial_residual,
        "local_ridge": local_ridge,
        "coverage": denominator,
        "valid": valid,
        "u": raw_state["u"],
        "phi": raw_state["phi"],
        "u_edges": raw_state["u_edges"],
        "phi_edges": raw_state["phi_edges"],
        "broad_azimuthal": broad_azimuthal,
    }


def build_log_polar_contrast(
        pixel_table, radial_parameters,
        n_u=LOGPOLAR_N_U, n_phi=LOGPOLAR_N_PHI):
    radius = pixel_table["radius_kpc"].to_numpy(float)
    azimuth = pixel_table["azimuth_rad"].to_numpy(float)
    observed_log = pixel_table["log_sfr"].to_numpy(float)
    background_log = (np.log10(radial_parameters["lambda0_0"])
                      - radius / (np.log(10.0) * radial_parameters["h_R"]))
    signed_residual = observed_log - background_log
    u = np.log(radius / R_REF_KPC)
    u_edges = np.quantile(u, np.linspace(0.0, 1.0, n_u + 1))
    if np.any(np.diff(u_edges) <= 0):
        raise RuntimeError("Quantile log-radius edges are not strictly increasing")
    u_edges[0] = np.nextafter(u_edges[0], -np.inf)
    u_edges[-1] = np.nextafter(u_edges[-1], np.inf)
    phi_edges = np.linspace(-np.pi, np.pi, n_phi + 1)
    raw_weighted_sum, _, _ = np.histogram2d(
        u, azimuth, bins=(u_edges, phi_edges), weights=signed_residual)
    raw_coverage, _, _ = np.histogram2d(
        u, azimuth, bins=(u_edges, phi_edges))
    row_index = np.clip(
        np.searchsorted(u_edges, u, side="right") - 1, 0, n_u - 1)
    row_count = np.bincount(row_index, minlength=n_u)
    u_centres = np.divide(
        np.bincount(row_index, weights=u, minlength=n_u), row_count,
        out=np.full(n_u, np.nan), where=row_count > 0)
    raw_state = {
        "raw_weighted_sum": raw_weighted_sum,
        "raw_coverage": raw_coverage,
        "u": u_centres,
        "phi": 0.5 * (phi_edges[:-1] + phi_edges[1:]),
        "u_edges": u_edges,
        "phi_edges": phi_edges,
    }
    return raw_state, preprocess_logpolar_raw(raw_state, include_local=True)


logpolar_raw, logpolar = build_log_polar_contrast(hii_pixels, radial_fit)
assert np.isfinite(logpolar["radial_residual"][logpolar["valid"]]).all()
assert len(hii_pixels) == int(valid_hii_pixels.sum())
```

Keep the existing 1-by-2 radial-residual/local-ridge diagnostic, now reading from `logpolar`.

- [ ] **Step 2: Add raw fold helpers to cell `16219bba` before the search**

```python
def sector_columns(phi, n_sectors=RIDGE_N_SECTORS):
    sector = np.floor(
        (np.asarray(phi) + np.pi) / (2.0*np.pi) * n_sectors).astype(int)
    return np.clip(sector, 0, n_sectors - 1)


def circular_dilate(mask, guard_bins):
    mask = np.asarray(mask, dtype=bool)
    return np.logical_or.reduce([
        np.roll(mask, shift)
        for shift in range(-guard_bins, guard_bins + 1)
    ])


def circular_erode(mask, guard_bins):
    mask = np.asarray(mask, dtype=bool)
    return np.logical_and.reduce([
        np.roll(mask, shift)
        for shift in range(-guard_bins, guard_bins + 1)
    ])


def logpolar_fold_masks(phi, held_sector, n_sectors=RIDGE_N_SECTORS,
                        guard_bins=RIDGE_GUARD_BINS):
    held = sector_columns(phi, n_sectors) == int(held_sector)
    train_allowed = ~circular_dilate(held, guard_bins)
    test_allowed = circular_erode(held, guard_bins)
    if np.any(train_allowed & test_allowed):
        raise RuntimeError("Training and test log-polar supports overlap")
    return train_allowed, test_allowed


def build_sector_fold_maps(raw_state, n_sectors=RIDGE_N_SECTORS,
                           guard_bins=RIDGE_GUARD_BINS):
    folds = []
    for held_sector in range(n_sectors):
        train_allowed, test_allowed = logpolar_fold_masks(
            raw_state["phi"], held_sector, n_sectors, guard_bins)
        folds.append({
            "held_sector": held_sector,
            "train": preprocess_logpolar_raw(
                raw_state, train_allowed, include_local=False),
            "test": preprocess_logpolar_raw(
                raw_state, test_allowed, include_local=False),
        })
    return folds
```

- [ ] **Step 3: Add the adversarial behavioral helper**

```python
def run_adversarial_leakage_check(
        n_u=12, n_phi=48, n_sectors=4,
        held_sector=1, guard_bins=2):
    u_edges = np.linspace(np.log(0.5), np.log(8.0), n_u + 1)
    phi_edges = np.linspace(-np.pi, np.pi, n_phi + 1)
    u = 0.5 * (u_edges[:-1] + u_edges[1:])
    phi = 0.5 * (phi_edges[:-1] + phi_edges[1:])
    held = sector_columns(phi, n_sectors) == held_sector
    raw_coverage = np.full((n_u, n_phi), 20.0)
    signal = np.zeros((n_u, n_phi), dtype=float)
    signal[:, held] = 0.4
    raw_state = {
        "raw_weighted_sum": raw_coverage * signal,
        "raw_coverage": raw_coverage,
        "u": u,
        "phi": phi,
        "u_edges": u_edges,
        "phi_edges": phi_edges,
    }
    train_allowed, _ = logpolar_fold_masks(
        phi, held_sector, n_sectors, guard_bins)
    raw_training_signal_l1 = float(np.sum(np.abs(
        raw_state["raw_weighted_sum"][:, train_allowed])))
    train = preprocess_logpolar_raw(
        raw_state, train_allowed, include_local=False)
    finite_train = train["radial_residual"][train["valid"]]
    max_abs_training_residual = (
        float(np.max(np.abs(finite_train))) if len(finite_train) else 0.0)
    return {
        "raw_training_signal_l1": raw_training_signal_l1,
        "max_abs_training_residual": max_abs_training_residual,
    }


adversarial_leakage = run_adversarial_leakage_check()
assert adversarial_leakage["raw_training_signal_l1"] == 0.0
assert adversarial_leakage["max_abs_training_residual"] <= 1e-14
print("RIDGE_LEAKAGE_GUARD_PASS")
```

- [ ] **Step 4: Run behavioral and source tests, then commit**

Run the leakage source contract, adversarial behavioral test, all-pixel test, and static parse.  Expected: all pass.  Commit:

```bash
git -C /Users/Igniz/Desktop/ICRAR add -- \
  further/20260713_NGC4254_KTZ_spiral_fit.ipynb \
  further/tests/test_20260713_ngc4254_ktz_spiral_fit.py
git -C /Users/Igniz/Desktop/ICRAR diff --cached --check
git -C /Users/Igniz/Desktop/ICRAR commit -m "feat: split NGC4254 ridge folds before smoothing"
```

---

### Task 3: Implement the memory-safe conditional ridge search

**Files:**

- Modify: `further/20260713_NGC4254_KTZ_spiral_fit.ipynb`
- Test: `further/tests/test_20260713_ngc4254_ktz_spiral_fit.py`

- [ ] **Step 1: Replace the Task 3 phase-search functions in cell `16219bba`**

Keep `wrap_angle`, `spiral_phase`, and `arm_profile` unchanged.  Replace the prior sector-histogram search with:

```python
def phase_row_histograms(logpolar_map, m_arms, pitch_angle,
                         field="radial_residual", n_phase=RIDGE_N_PHASE):
    uu, pp = np.meshgrid(
        logpolar_map["u"], logpolar_map["phi"], indexing="ij")
    valid = logpolar_map["valid"] & np.isfinite(logpolar_map[field])
    row = np.broadcast_to(
        np.arange(len(logpolar_map["u"]))[:, None], valid.shape)[valid]
    phase = np.mod(
        (m_arms / np.tan(np.deg2rad(pitch_angle))) * uu[valid]
        - m_arms * pp[valid], 2.0*np.pi)
    phase_index = (
        np.floor(phase / (2.0*np.pi) * n_phase).astype(int) % n_phase)
    joint = row * n_phase + phase_index
    size = len(logpolar_map["u"]) * n_phase
    coverage = logpolar_map["coverage"][valid]
    values = logpolar_map[field][valid]
    weighted = np.bincount(
        joint, weights=coverage * values, minlength=size).reshape(
            len(logpolar_map["u"]), n_phase)
    support = np.bincount(
        joint, weights=coverage, minlength=size).reshape(
            len(logpolar_map["u"]), n_phase)
    return weighted, support


def variable_width_ridge_response(
        weighted_hist, support_hist, u_centres, m_arms, pitch_angle,
        core_width_kpc=RIDGE_CORE_WIDTH_KPC,
        broad_ratio=RIDGE_BROAD_RATIO,
        n_radial_bands=LOGPOLAR_N_RADIAL_BANDS):
    n_u, n_phase = weighted_hist.shape
    if n_u % n_radial_bands:
        raise ValueError("LOGPOLAR_N_U must be divisible by radial-band count")
    rows_per_band = n_u // n_radial_bands
    weighted_band = weighted_hist.reshape(
        n_radial_bands, rows_per_band, n_phase).sum(axis=1)
    support_band = support_hist.reshape(
        n_radial_bands, rows_per_band, n_phase).sum(axis=1)
    u_band = np.asarray(u_centres).reshape(
        n_radial_bands, rows_per_band).mean(axis=1)
    radius_band = R_REF_KPC * np.exp(u_band)
    sigma_phase = (m_arms * core_width_kpc
                   / (radius_band * abs(np.sin(np.deg2rad(pitch_angle)))))
    raw_sigma_bins = sigma_phase / (2.0*np.pi) * n_phase
    sigma_bins = np.clip(raw_sigma_bins, 0.65, n_phase / 10.0)
    clamp_low = raw_sigma_bins < 0.65
    clamp_high = raw_sigma_bins > n_phase / 10.0

    narrow_mean = np.full_like(weighted_band, np.nan, dtype=float)
    broad_mean = np.full_like(weighted_band, np.nan, dtype=float)
    narrow_support = np.zeros_like(support_band, dtype=float)
    for radial_index, sigma_bin in enumerate(sigma_bins):
        narrow_weighted = gaussian_filter1d(
            weighted_band[radial_index], sigma_bin, mode="wrap")
        narrow_denominator = gaussian_filter1d(
            support_band[radial_index], sigma_bin, mode="wrap")
        broad_weighted = gaussian_filter1d(
            weighted_band[radial_index], broad_ratio * sigma_bin, mode="wrap")
        broad_denominator = gaussian_filter1d(
            support_band[radial_index], broad_ratio * sigma_bin, mode="wrap")
        narrow_mean[radial_index] = np.divide(
            narrow_weighted, narrow_denominator,
            out=np.full(n_phase, np.nan), where=narrow_denominator > 0.05)
        broad_mean[radial_index] = np.divide(
            broad_weighted, broad_denominator,
            out=np.full(n_phase, np.nan), where=broad_denominator > 0.05)
        narrow_support[radial_index] = narrow_denominator
    return {
        "response": narrow_mean - broad_mean,
        "narrow_mean": narrow_mean,
        "broad_mean": broad_mean,
        "support": narrow_support,
        "clamp_low_count": int(clamp_low.sum()),
        "clamp_high_count": int(clamp_high.sum()),
    }


def aggregate_radial_response(response, support, min_radial_fraction=0.35):
    good = np.isfinite(response) & (support > 0.05)
    count = good.sum(axis=0)
    mean_response = np.divide(
        np.where(good, response, 0.0).sum(axis=0), count,
        out=np.full(response.shape[1], np.nan), where=count > 0)
    positive_fraction = np.divide(
        (good & (response > 0)).sum(axis=0), count,
        out=np.full(response.shape[1], np.nan), where=count > 0)
    median_response = np.array([
        np.median(response[good[:, index], index]) if count[index] else np.nan
        for index in range(response.shape[1])
    ])
    radial_fraction = count / response.shape[0]
    score = mean_response * positive_fraction + 0.5 * median_response
    score[radial_fraction < min_radial_fraction] = np.nan
    return {
        "score": score,
        "mean": mean_response,
        "median": median_response,
        "positive_fraction": positive_fraction,
        "radial_fraction": radial_fraction,
    }


def ridge_curve_for_map(logpolar_map, m_arms, pitch_angle,
                        core_width_kpc, min_radial_fraction=0.35):
    weighted, support = phase_row_histograms(
        logpolar_map, m_arms, pitch_angle, field="radial_residual")
    response = variable_width_ridge_response(
        weighted, support, logpolar_map["u"], m_arms, pitch_angle,
        core_width_kpc=core_width_kpc)
    aggregate = aggregate_radial_response(
        response["response"], response["support"], min_radial_fraction)
    return response, aggregate


def conditional_ridge_search(
        raw_state, m_compare=M_COMPARE,
        pitch_grid=RIDGE_PITCH_GRID_DEG,
        core_width_kpc=RIDGE_CORE_WIDTH_KPC,
        guard_bins=RIDGE_GUARD_BINS,
        allow_no_acceptance=False):
    phase_values = np.linspace(
        0.0, 2.0*np.pi, RIDGE_N_PHASE, endpoint=False)
    full_map = preprocess_logpolar_raw(
        raw_state, allowed_phi=None, include_local=False)
    records = []
    branch_normalization = "mean_not_sum"
    for m_arms in np.asarray(m_compare, dtype=int):
        for pitch_angle in np.asarray(pitch_grid, dtype=float):
            response, aggregate = ridge_curve_for_map(
                full_map, int(m_arms), float(pitch_angle), core_width_kpc)
            if not np.isfinite(aggregate["score"]).any():
                continue
            best_index = int(np.nanargmax(aggregate["score"]))
            records.append({
                "m_arms": int(m_arms),
                "pitch_angle": float(pitch_angle),
                "winding_sign": -1,
                "Theta0": float((-phase_values[best_index]) % (2.0*np.pi)),
                "core_width_kpc": float(core_width_kpc),
                "full_ridge_score": float(aggregate["score"][best_index]),
                "narrow_mean": float(np.nanmean(
                    response["narrow_mean"][:, best_index])),
                "broad_flank_mean": float(np.nanmean(
                    response["broad_mean"][:, best_index])),
                "positive_radial_fraction": float(
                    aggregate["positive_fraction"][best_index]),
                "radial_coverage": float(
                    aggregate["radial_fraction"][best_index]),
                "phase_width_clamp_low_bands": response["clamp_low_count"],
                "phase_width_clamp_high_bands": response["clamp_high_count"],
                "pitch_boundary": bool(np.isclose(pitch_angle, pitch_grid[0])
                                       or np.isclose(pitch_angle, pitch_grid[-1])),
                "branch_normalization": branch_normalization,
                "held_out_score": np.nan,
                "held_out_score_std": np.nan,
                "phase_stability": np.nan,
                "valid_held_out": 0,
                "validated_score": np.nan,
            })
    table = pd.DataFrame(records)
    if table.empty:
        raise RuntimeError("No conditional ridge candidate has finite support")

    shortlist = (table.sort_values("full_ridge_score", ascending=False)
                 .groupby("m_arms", group_keys=False)
                 .head(RIDGE_SHORTLIST_PER_FAMILY).index)
    folds = build_sector_fold_maps(
        raw_state, RIDGE_N_SECTORS, guard_bins)
    fold_records = []
    for record_index in shortlist:
        m_arms = int(table.at[record_index, "m_arms"])
        pitch_angle = float(table.at[record_index, "pitch_angle"])
        held_scores = []
        train_phases = []
        for fold in folds:
            _, train = ridge_curve_for_map(
                fold["train"], m_arms, pitch_angle, core_width_kpc)
            if not np.isfinite(train["score"]).any():
                continue
            train_index = int(np.nanargmax(train["score"]))
            _, test = ridge_curve_for_map(
                fold["test"], m_arms, pitch_angle, core_width_kpc,
                min_radial_fraction=0.08)
            test_score = float(test["score"][train_index])
            if not np.isfinite(test_score):
                continue
            train_theta0 = float(
                (-phase_values[train_index]) % (2.0*np.pi))
            held_scores.append(test_score)
            train_phases.append(train_theta0)
            fold_records.append({
                "m_arms": m_arms,
                "pitch_angle": pitch_angle,
                "held_sector": int(fold["held_sector"]),
                "test_score": test_score,
                "train_Theta0": train_theta0,
            })
        valid_held_out = len(held_scores)
        if valid_held_out:
            held_out_score = float(np.median(held_scores))
            phase_stability = float(np.abs(np.mean(
                np.exp(1j * np.asarray(train_phases)))))
            table.at[record_index, "held_out_score"] = held_out_score
            table.at[record_index, "held_out_score_std"] = float(
                np.std(held_scores))
            table.at[record_index, "phase_stability"] = phase_stability
            table.at[record_index, "valid_held_out"] = valid_held_out
            table.at[record_index, "validated_score"] = (
                max(held_out_score, 0.0)
                * valid_held_out / RIDGE_N_SECTORS
                * phase_stability)

    geometry_rows = []
    for m_arms in np.asarray(m_compare, dtype=int):
        accepted = table.loc[
            (table["m_arms"] == m_arms)
            & (table["valid_held_out"] >= RIDGE_MIN_HELD_OUT_SECTORS)
            & (table["validated_score"] > 0)].sort_values(
                ["validated_score", "full_ridge_score"], ascending=False)
        if accepted.empty:
            if allow_no_acceptance:
                continue
            raise RuntimeError(
                f"No leakage-controlled negative ridge accepted for m={m_arms}")
        geometry = accepted.iloc[0].to_dict()
        geometry["geometry_source"] = (
            "leakage_controlled_negative_winding_ridge")
        geometry_rows.append(geometry)
    histogram_cache_payload_bytes = 0
    return (pd.DataFrame(geometry_rows).sort_values("m_arms").reset_index(drop=True),
            table.sort_values(
                ["m_arms", "validated_score", "full_ridge_score"],
                ascending=[True, False, False], na_position="last").reset_index(drop=True),
            pd.DataFrame(fold_records), histogram_cache_payload_bytes)
```

- [ ] **Step 2: Replace the synthetic Section 10**

Retain the existing physical-ridge generator but use exactly these negative-winding cases:

```python
synthetic_cases = [(2, -22.0), (3, -30.0), (4, -38.0)]
theta0_true = 0.7
synthetic_recoveries = []
for case_index, (m_true, pitch_true) in enumerate(synthetic_cases):
    pixels, radial_parameters = synthetic_sfr_pixels(
        m_true, pitch_true, theta0_true, RNG_SEED + case_index)
    synthetic_raw, _ = build_log_polar_contrast(pixels, radial_parameters)
    recovered, _, _, _ = conditional_ridge_search(
        synthetic_raw, m_compare=np.array([m_true]),
        pitch_grid=RIDGE_PITCH_GRID_DEG)
    geometry = recovered.iloc[0].to_dict()
    synthetic_recoveries.append(geometry)
    assert int(geometry["m_arms"]) == m_true
    assert geometry["pitch_angle"] < 0
    assert abs(geometry["pitch_angle"] - pitch_true) <= 2.0
    assert abs(wrap_angle(geometry["Theta0"] - theta0_true)) <= np.deg2rad(5.0)
print("RIDGE_M234_SYNTHETIC_PASS")

synthetic_full = preprocess_logpolar_raw(
    synthetic_raw, allowed_phi=None, include_local=False)
test_weighted, test_support = phase_row_histograms(
    synthetic_full, 4, -38.0)
response_1 = variable_width_ridge_response(
    test_weighted, test_support, synthetic_full["u"], 4, -38.0)
response_6 = variable_width_ridge_response(
    6.0 * test_weighted, 6.0 * test_support,
    synthetic_full["u"], 4, -38.0)
assert np.allclose(
    response_1["response"], response_6["response"],
    equal_nan=True, atol=1e-12)
print("RIDGE_BRANCH_NORMALIZATION_PASS")
```

- [ ] **Step 3: Run the synthetic block in memory**

Execute notebook code only through cell `c9d26c41` under Matplotlib `Agg`.  Require the leakage marker, the two synthetic markers, preservation of each supplied fixed arm number, recovery of negative pitch and phase within the declared tolerances, explicit wording that this is not arm-number discrimination, and zero `RuntimeWarning`.

- [ ] **Step 4: Verify memory-safe source and commit**

Run the two search contracts and static parse.  Confirm `histogram_cache_payload_bytes == 0` in the disposable execution.  Commit:

```bash
git -C /Users/Igniz/Desktop/ICRAR add -- \
  further/20260713_NGC4254_KTZ_spiral_fit.ipynb \
  further/tests/test_20260713_ngc4254_ktz_spiral_fit.py
git -C /Users/Igniz/Desktop/ICRAR diff --cached --check
git -C /Users/Igniz/Desktop/ICRAR commit -m "feat: compare leakage-controlled m234 ridges"
```

---

### Task 4: Add descriptive blocked nulls, width sensitivity, and sector diagnostics

**Files:**

- Modify: `further/20260713_NGC4254_KTZ_spiral_fit.ipynb`
- Test: `further/tests/test_20260713_ngc4254_ktz_spiral_fit.py`

- [ ] **Step 1: Replace Section 11 with the three-case real search and blocked null helpers**

Use this complete code:

```python
def scramble_logpolar_raw(raw_state, rng):
    result = {
        key: (value.copy() if isinstance(value, np.ndarray) else value)
        for key, value in raw_state.items()
    }
    for row_index in range(len(result["u"])):
        shift = int(rng.integers(0, len(result["phi"])))
        result["raw_weighted_sum"][row_index] = np.roll(
            result["raw_weighted_sum"][row_index], shift)
        result["raw_coverage"][row_index] = np.roll(
            result["raw_coverage"][row_index], shift)
    return result


def build_blocked_null_diagnostics(raw_state):
    rows = []
    for block_index, block_seed in enumerate(RIDGE_NULL_BLOCK_SEEDS):
        block_rng = np.random.default_rng(int(block_seed))
        for draw_index in range(RIDGE_NULL_DRAWS_PER_BLOCK):
            scrambled = scramble_logpolar_raw(raw_state, block_rng)
            geometries, _, _, _ = conditional_ridge_search(
                scrambled, allow_no_acceptance=True)
            for m_arms in M_COMPARE:
                match = geometries.loc[geometries["m_arms"] == m_arms]
                if len(match) != 1:
                    raise RuntimeError(
                        f"Expected one blocked-null row for m={m_arms}")
                validated_score = float(match.iloc[0]["validated_score"])
                if not np.isfinite(validated_score):
                    raise RuntimeError(
                        f"Non-finite blocked-null score for m={m_arms}")
                rows.append({
                    "block_index": block_index,
                    "block_seed": int(block_seed),
                    "draw_index": draw_index,
                    "m_arms": int(m_arms),
                    "validated_score": validated_score,
                    "best_null_score": validated_score,
                })
    draws = pd.DataFrame(rows)
    expected = (len(RIDGE_NULL_BLOCK_SEEDS)
                * RIDGE_NULL_DRAWS_PER_BLOCK * len(M_COMPARE))
    if len(draws) != expected:
        raise RuntimeError(f"Expected {expected} blocked null rows; got {len(draws)}")
    pooled = (draws.groupby("m_arms")["best_null_score"]
              .agg(null_mean="mean", null_std="std", null_count="count")
              .reset_index())
    pooled["null_std_floor"] = np.maximum(pooled["null_std"], 1e-6)
    blocks = (draws.groupby(["block_index", "block_seed", "m_arms"])
              ["best_null_score"]
              .agg(block_null_mean="mean", block_null_std="std",
                   block_null_count="count")
              .reset_index())
    return draws, pooled, blocks


ridge_geometry_table, ridge_candidate_table, ridge_fold_table, cache_bytes = (
    conditional_ridge_search(logpolar_raw))
assert cache_bytes == 0
assert ridge_geometry_table["m_arms"].tolist() == M_COMPARE.tolist()
assert (ridge_geometry_table["pitch_angle"] < 0).all()

null_draws, pooled_null_summary, null_block_summary = (
    build_blocked_null_diagnostics(logpolar_raw))
ridge_geometry_table = ridge_geometry_table.merge(
    pooled_null_summary, on="m_arms", how="left", validate="one_to_one")
ridge_geometry_table["null_z"] = (
    (ridge_geometry_table["validated_score"]
     - ridge_geometry_table["null_mean"])
    / ridge_geometry_table["null_std_floor"])
ridge_geometry_table["null_role"] = "descriptive_not_model_selection"
display(ridge_geometry_table)
display(null_block_summary)
print("RIDGE_NULL_BLOCKS_COMPLETE")
```

- [ ] **Step 2: Append the five-width conditional search**

```python
width_rows = []
for width_kpc in RIDGE_WIDTH_SENSITIVITY_KPC:
    width_geometries, _, _, _ = conditional_ridge_search(
        logpolar_raw, core_width_kpc=float(width_kpc))
    for row in width_geometries.itertuples(index=False):
        width_rows.append({
            "core_width_kpc": float(width_kpc),
            "m_arms": int(row.m_arms),
            "pitch_angle": float(row.pitch_angle),
            "Theta0": float(row.Theta0),
            "validated_score": float(row.validated_score),
            "phase_stability": float(row.phase_stability),
            "valid_held_out": int(row.valid_held_out),
            "pitch_boundary": bool(row.pitch_boundary),
        })
ridge_width_sensitivity = pd.DataFrame(width_rows)
assert len(ridge_width_sensitivity) == len(M_COMPARE) * len(
    RIDGE_WIDTH_SENSITIVITY_KPC)
display(ridge_width_sensitivity)
print("RIDGE_WIDTH_M234_COMPLETE")
```

- [ ] **Step 3: Replace Section 12 with selected-candidate fold diagnostics**

```python
selected_fold_rows = []
for geometry in ridge_geometry_table.to_dict("records"):
    use = ridge_fold_table.loc[
        (ridge_fold_table["m_arms"] == int(geometry["m_arms"]))
        & np.isclose(ridge_fold_table["pitch_angle"],
                     float(geometry["pitch_angle"]))].copy()
    use["selected_Theta0"] = float(geometry["Theta0"])
    selected_fold_rows.append(use)
ridge_sector_results = pd.concat(selected_fold_rows, ignore_index=True)
assert len(ridge_sector_results) >= len(M_COMPARE) * RIDGE_MIN_HELD_OUT_SECTORS

fig, ax = plt.subplots(figsize=(9.2, 5.2), constrained_layout=True)
for m_arms, group in ridge_sector_results.groupby("m_arms"):
    ax.plot(group["held_sector"], group["test_score"], "o-",
            label=f"m={m_arms}")
ax.axhline(0.0, color="0.5", lw=0.8)
ax.set(xlabel="held azimuth-sector index",
       ylabel="leakage-controlled test response",
       title="Conditional m=2,3,4 held-sector responses")
ax.legend()
plt.show()
print("RIDGE_SECTOR_M234_COMPLETE")
print("RIDGE_M234_REAL_COMPLETE")
```

- [ ] **Step 4: Execute the real conditional search and validate**

Run through Section 12 against the live FITS file under `Agg`.  Require:

- three negative-winding geometries in arm-number order;
- at least 10 finite held sectors each;
- positive validated score each;
- exactly 96 null rows (`3 * 32`);
- exactly 15 width rows (`3 * 5`);
- zero histogram cache payload; and
- zero `RuntimeWarning`.

Save scratch evidence under `/private/tmp/ngc4254_m234_*` for parent verification.

- [ ] **Step 5: Verify and commit**

Run the null/source contracts and static parse.  Commit:

```bash
git -C /Users/Igniz/Desktop/ICRAR add -- \
  further/20260713_NGC4254_KTZ_spiral_fit.ipynb \
  further/tests/test_20260713_ngc4254_ktz_spiral_fit.py
git -C /Users/Igniz/Desktop/ICRAR diff --cached --check
git -C /Users/Igniz/Desktop/ICRAR commit -m "feat: add descriptive m234 null and width diagnostics"
```

---

### Task 5: Fit three fixed-geometry KTZ-compatible profiles

**Files:**

- Modify: `further/20260713_NGC4254_KTZ_spiral_fit.ipynb`
- Test: `further/tests/test_20260713_ngc4254_ktz_spiral_fit.py`

- [ ] **Step 1: Replace Section 9 with the fixed-geometry fitter and maxima enumerator**

Use the complete `fit_ktz_profile_fixed_geometry` and `enumerate_profile_maxima` implementations below:

```python
def fit_ktz_profile_fixed_geometry(table, geometry, radial_guess):
    radius = table["radius_kpc"].to_numpy(float)
    azimuth = table["azimuth_rad"].to_numpy(float)
    observed_log = np.log(table["sfr_linear"].to_numpy(float))
    sqrt_weight = np.sqrt(table["area_pix"].to_numpy(float))
    sqrt_weight /= np.nanmedian(sqrt_weight)
    theta = spiral_phase(
        radius, azimuth, int(geometry["m_arms"]),
        float(geometry["pitch_angle"]), float(geometry["Theta0"]))
    initial = np.array([
        np.log(radial_guess["lambda0_0"]), np.log(radial_guess["h_R"]),
        np.log(0.2), 0.3, 0.15, 0.0, 0.0])
    lower = np.array([
        -50.0, np.log(0.2), np.log(1e-4), 0.0, 0.0, -np.pi, -np.pi])
    upper = np.array([
        50.0, np.log(30.0), np.log(2.0), 1.5, 1.5, np.pi, np.pi])
    dense_theta = np.linspace(0.0, 2.0*np.pi, 8192, endpoint=False)

    def unpack(values):
        return {
            "lambda0_0": float(np.exp(values[0])),
            "h_R": float(np.exp(values[1])),
            "eta": float(np.exp(values[2])),
            "harmonic_n": HARMONIC_N.copy(),
            "harmonic_g": np.array([1.0, values[3], values[4]]),
            "harmonic_alpha": np.array([0.0, values[5], values[6]]),
        }

    def residual(values):
        parameters = unpack(values)
        factor = 1.0 + parameters["eta"] * arm_profile(
            theta, parameters["harmonic_g"], parameters["harmonic_alpha"])
        dense_factor = 1.0 + parameters["eta"] * arm_profile(
            dense_theta, parameters["harmonic_g"],
            parameters["harmonic_alpha"])
        model_log = (np.log(parameters["lambda0_0"])
                     - radius / parameters["h_R"]
                     + np.log(np.clip(factor, 1e-10, None)))
        positivity_penalty = min(
            0.0, float(np.min(dense_factor)) - 1e-3) * 1e5
        return np.concatenate([
            sqrt_weight * (observed_log - model_log),
            [positivity_penalty]])

    result = least_squares(
        residual, initial, bounds=(lower, upper),
        loss="soft_l1", f_scale=0.25, max_nfev=3000)
    parameters = unpack(result.x)
    final_residual = residual(result.x)
    minimum_factor = float(np.min(
        1.0 + parameters["eta"] * arm_profile(
            dense_theta, parameters["harmonic_g"],
            parameters["harmonic_alpha"])))
    parameters.update({
        "m_arms": int(geometry["m_arms"]),
        "pitch_angle": float(geometry["pitch_angle"]),
        "Theta0": float(geometry["Theta0"]),
        "geometry_source": geometry["geometry_source"],
        "weighted_log_sse": float(np.sum(final_residual[:-1]**2)),
        "robust_cost": float(2.0 * result.cost),
        "minimum_rate_factor": minimum_factor,
        "success": bool(result.success and minimum_factor > 0.0),
        "message": result.message,
    })
    if not parameters["success"]:
        raise RuntimeError(f"Fixed-geometry KTZ fit failed: {parameters}")
    return parameters


def enumerate_profile_maxima(model, n_grid=65536):
    theta = np.linspace(0.0, 2.0*np.pi, n_grid, endpoint=False)
    profile = arm_profile(
        theta, model["harmonic_g"], model["harmonic_alpha"])
    maxima = (profile > np.roll(profile, 1)) & (
        profile >= np.roll(profile, -1))
    indices = np.flatnonzero(maxima)
    if len(indices) == 0:
        raise RuntimeError("No harmonic-profile maximum found")
    grid_step = 2.0*np.pi / n_grid

    def profile_scalar(value):
        wrapped_value = float(value % (2.0*np.pi))
        return float(arm_profile(
            np.array([wrapped_value]), model["harmonic_g"],
            model["harmonic_alpha"])[0])

    rows = []
    for index in indices:
        centre = float(theta[index])
        refined = minimize_scalar(
            lambda offset: -profile_scalar(centre + offset),
            bounds=(-grid_step, grid_step), method="bounded",
            options={"xatol": 1e-12})
        theta_peak = float((centre + refined.x) % (2.0*np.pi))
        h_peak = profile_scalar(theta_peak)
        rate_factor = float(1.0 + model["eta"] * h_peak)
        rows.append({
            "theta_peak": theta_peak,
            "h_peak": h_peak,
            "rate_factor": rate_factor,
            "enhanced_above_background": bool(rate_factor > 1.0),
        })
    return pd.DataFrame(rows).sort_values(
        "rate_factor", ascending=False).reset_index(drop=True)
```

- [ ] **Step 2: Append three profile fits to the real-data cell after ridge completion**

```python
conditional_fits = {}
conditional_fit_tables = {}
profile_maxima_by_m = {}
for geometry in ridge_geometry_table.to_dict("records"):
    m_arms = int(geometry["m_arms"])
    model = fit_ktz_profile_fixed_geometry(fit_table, geometry, radial_fit)
    assert model["m_arms"] == m_arms
    assert np.isclose(model["pitch_angle"], geometry["pitch_angle"])
    assert np.isclose(model["Theta0"], geometry["Theta0"])
    model["gradient_dex_per_kpc"] = float(
        -1.0 / (np.log(10.0) * model["h_R"]))
    table_m = fit_table.copy()
    theta_fit = spiral_phase(
        table_m["radius_kpc"].to_numpy(float),
        table_m["azimuth_rad"].to_numpy(float),
        model["m_arms"], model["pitch_angle"], model["Theta0"])
    rate_factor = 1.0 + model["eta"] * arm_profile(
        theta_fit, model["harmonic_g"], model["harmonic_alpha"])
    model_sfr = (model["lambda0_0"]
                 * np.exp(-table_m["radius_kpc"].to_numpy(float)
                          / model["h_R"])
                 * rate_factor)
    table_m["theta_fit"] = np.mod(theta_fit, 2.0*np.pi)
    table_m["rate_factor"] = rate_factor
    table_m["model_sfr"] = model_sfr
    table_m["residual_dex"] = np.log10(
        table_m["sfr_linear"].to_numpy(float)
        / np.clip(model_sfr, 1e-30, None))
    maxima = enumerate_profile_maxima(model)
    maxima.insert(0, "m_arms", m_arms)
    conditional_fits[m_arms] = model
    conditional_fit_tables[m_arms] = table_m
    profile_maxima_by_m[m_arms] = maxima

assert sorted(conditional_fits) == M_COMPARE.tolist()
assert all(model["minimum_rate_factor"] > 0
           for model in conditional_fits.values())
print("KTZ_M234_PROFILE_FITS_COMPLETE")
print("HARMONIC_M234_MAXIMA_COMPLETE")
print("POSITIVITY_CHECK_PASS")
```

- [ ] **Step 3: Run the three fits in a disposable partial execution**

Load the live ridge geometry table from the Task 4 scratch evidence or execute through Section 12, then run Section 9 definitions and the appended invocation.  Verify all three fits succeed, geometry is unchanged, minima are positive, and maxima tables are nonempty.

- [ ] **Step 4: Verify and commit**

Run the three-fit contract and static parse.  Commit:

```bash
git -C /Users/Igniz/Desktop/ICRAR add -- \
  further/20260713_NGC4254_KTZ_spiral_fit.ipynb \
  further/tests/test_20260713_ngc4254_ktz_spiral_fit.py
git -C /Users/Igniz/Desktop/ICRAR diff --cached --check
git -C /Users/Igniz/Desktop/ICRAR commit -m "feat: fit conditional KTZ profiles for m234"
```

---

### Task 6: Build comparative skeleton, model, residual, profile, and table outputs

**Files:**

- Modify: `further/20260713_NGC4254_KTZ_spiral_fit.ipynb`
- Test: `further/tests/test_20260713_ngc4254_ktz_spiral_fit.py`

- [ ] **Step 1: Replace Section 13 with projection helpers and the 2-by-2 comparison**

Retain the validated `skeleton_xy_for_phase` and `project_disc_to_pixel` formulas from the previous design.  Use this complete comparison setup after the round-trip and synthetic sign checks:

```python
M_COLORS = {2: "cyan", 3: "lime", 4: "dodgerblue"}
radius_grid = np.geomspace(RIDGE_R_IN_KPC, RIDGE_R_OUT_KPC, 1000)
data_ridge_curves_by_m = {
    int(row.m_arms): skeleton_xy_for_phase(
        row._asdict(), radius_grid, theta_target=0.0)
    for row in ridge_geometry_table.itertuples(index=False)
}

fig, axes = plt.subplots(2, 2, figsize=(14.5, 12.0), constrained_layout=True)
ax_sky, ax_face, ax_profile, ax_width = axes.ravel()

sky_image = ax_sky.imshow(
    log_sfr_map, origin="lower", cmap="magma", vmin=vmin, vmax=vmax)
for m_arms, curves in data_ridge_curves_by_m.items():
    for branch_index, (x_curve, y_curve) in enumerate(curves):
        xp, yp = project_disc_to_pixel(x_curve, y_curve)
        ax_sky.plot(
            xp, yp, color=M_COLORS[m_arms], lw=1.5,
            label=f"m={m_arms}, negative winding" if branch_index == 0
            else "_nolegend_")
ax_sky.set(xlabel="image x (pixel)", ylabel="image y (pixel)",
           title="observed sky-plane log SFR + conditional data ridges")
ax_sky.set_aspect("equal")
ax_sky.legend(fontsize=8)
fig.colorbar(sky_image, ax=ax_sky, shrink=0.82,
             label=r"log $\Sigma_{\rm SFR}$")

face = ax_face.scatter(
    hii_pixels["x_disc_kpc"], hii_pixels["y_disc_kpc"],
    c=hii_pixels["log_sfr"], s=0.35, cmap="magma",
    vmin=vmin, vmax=vmax, linewidth=0, rasterized=True)
for m_arms, curves in data_ridge_curves_by_m.items():
    for x_curve, y_curve in curves:
        ax_face.plot(x_curve, y_curve, color=M_COLORS[m_arms], lw=1.5)
ax_face.set(xlabel="disc major-axis x (kpc)",
            ylabel="disc minor-axis y (kpc)",
            title="observed deprojected log SFR + conditional data ridges")
ax_face.set_aspect("equal")
fig.colorbar(face, ax=ax_face, shrink=0.82,
             label=r"log $\Sigma_{\rm SFR}$")

theta_grid = np.linspace(0.0, 2.0*np.pi, 1600)
for m_arms in M_COMPARE:
    model = conditional_fits[int(m_arms)]
    profile = arm_profile(
        theta_grid, model["harmonic_g"], model["harmonic_alpha"])
    ax_profile.plot(theta_grid, profile, color=M_COLORS[int(m_arms)],
                    lw=2.0, label=f"m={m_arms}")
    maxima = profile_maxima_by_m[int(m_arms)]
    ax_profile.scatter(maxima["theta_peak"], maxima["h_peak"],
                       color=M_COLORS[int(m_arms)], s=34)
ax_profile.set(xlabel=r"$\Theta$", ylabel=r"$h(\Theta)$",
               title="conditional KTZ harmonic profiles + all maxima",
               xlim=(0.0, 2.0*np.pi))
ax_profile.legend()

for m_arms, group in ridge_width_sensitivity.groupby("m_arms"):
    ax_width.plot(group["core_width_kpc"], group["pitch_angle"], "o-",
                  color=M_COLORS[int(m_arms)], label=f"m={m_arms}")
    fiducial = ridge_geometry_table.loc[
        ridge_geometry_table["m_arms"] == m_arms].iloc[0]
    ax_width.annotate(
        f"pooled null z={fiducial['null_z']:.2f}",
        (RIDGE_CORE_WIDTH_KPC, fiducial["pitch_angle"]),
        xytext=(5, 7), textcoords="offset points", fontsize=7,
        color=M_COLORS[int(m_arms)])
ax_width.axvline(RIDGE_CORE_WIDTH_KPC, color="black", ls="--", lw=1)
ax_width.set(xlabel="physical ridge width (kpc)",
             ylabel="selected negative pitch (deg)",
             title="ridge-width and descriptive-null comparison")
ax_width.legend()
plt.show()
print("M234_COMPARISON_LAYOUT_COMPLETE")
```

- [ ] **Step 2: Replace Section 14 with the 3-by-2 model/residual figure**

```python
fig, axes = plt.subplots(3, 2, figsize=(13.5, 18.0), constrained_layout=True)
all_residuals = np.concatenate([
    conditional_fit_tables[int(m)]["residual_dex"].to_numpy(float)
    for m in M_COMPARE])
residual_limit = float(np.nanpercentile(np.abs(all_residuals), 95))
for row_index, m_arms in enumerate(M_COMPARE):
    m_arms = int(m_arms)
    model = conditional_fits[m_arms]
    table_m = conditional_fit_tables[m_arms]
    ax_model, ax_residual = axes[row_index]
    model_scatter = ax_model.scatter(
        table_m["x_disc_kpc"], table_m["y_disc_kpc"],
        c=np.log10(table_m["model_sfr"]), s=2.0,
        cmap="magma", vmin=vmin, vmax=vmax, linewidth=0)
    for peak in profile_maxima_by_m[m_arms].itertuples(index=False):
        if not peak.enhanced_above_background:
            continue
        for x_curve, y_curve in skeleton_xy_for_phase(
                model, radius_grid, theta_target=peak.theta_peak):
            ax_model.plot(x_curve, y_curve, color=M_COLORS[m_arms], lw=1.1)
    ax_model.set(xlabel="disc major-axis x (kpc)",
                 ylabel="disc minor-axis y (kpc)",
                 title=f"m={m_arms} conditional source model")
    ax_model.set_aspect("equal")
    fig.colorbar(model_scatter, ax=ax_model,
                 label=r"model log $\Sigma_{\rm SFR}$")

    residual_scatter = ax_residual.scatter(
        table_m["x_disc_kpc"], table_m["y_disc_kpc"],
        c=table_m["residual_dex"], s=2.0, cmap="coolwarm",
        vmin=-residual_limit, vmax=residual_limit, linewidth=0)
    for x_curve, y_curve in data_ridge_curves_by_m[m_arms]:
        ax_residual.plot(x_curve, y_curve,
                         color=M_COLORS[m_arms], lw=1.3)
    ax_residual.set(xlabel="disc major-axis x (kpc)",
                    ylabel="disc minor-axis y (kpc)",
                    title=f"m={m_arms} observed - model residual")
    ax_residual.set_aspect("equal")
    fig.colorbar(residual_scatter, ax=ax_residual, label="residual (dex)")
plt.show()
print("M234_MODEL_RESIDUAL_LAYOUT_COMPLETE")
```

- [ ] **Step 3: Replace Section 15 with the three-row table**

```python
parameter_rows = []
for geometry in ridge_geometry_table.to_dict("records"):
    m_arms = int(geometry["m_arms"])
    model = conditional_fits[m_arms]
    maxima = profile_maxima_by_m[m_arms]
    blocks = null_block_summary.loc[
        null_block_summary["m_arms"] == m_arms,
        "block_null_mean"].to_numpy(float)
    parameter_rows.append({
        "m_arms": m_arms,
        "pitch_angle_deg": geometry["pitch_angle"],
        "Theta0_rad": geometry["Theta0"],
        "geometry_source": geometry["geometry_source"],
        "validated_ridge_score": geometry["validated_score"],
        "valid_held_sectors": int(geometry["valid_held_out"]),
        "phase_stability": geometry["phase_stability"],
        "pooled_null_mean": geometry["null_mean"],
        "pooled_null_std": geometry["null_std"],
        "descriptive_null_z": geometry["null_z"],
        "block_null_mean_range": f"{blocks.min():.5f}..{blocks.max():.5f}",
        "pitch_boundary": bool(geometry["pitch_boundary"]),
        "lambda0_0": model["lambda0_0"],
        "h_R_kpc": model["h_R"],
        "gradient_dex_per_kpc": model["gradient_dex_per_kpc"],
        "eta": model["eta"],
        "harmonic_n": model["harmonic_n"].tolist(),
        "harmonic_g": model["harmonic_g"].tolist(),
        "harmonic_alpha": model["harmonic_alpha"].tolist(),
        "minimum_rate_factor": model["minimum_rate_factor"],
        "enhanced_maxima": int(maxima["enhanced_above_background"].sum()),
        "weighted_log_sse": model["weighted_log_sse"],
    })
parameter_table = pd.DataFrame(parameter_rows).sort_values("m_arms")
assert parameter_table["m_arms"].tolist() == [2, 3, 4]
display(parameter_table)
for m_arms in M_COMPARE:
    display(profile_maxima_by_m[int(m_arms)])
print("FINAL_M234_PARAMETER_TABLE_COMPLETE")
```

Update Section 16 Markdown to state all three results are conditional, negative winding is imposed to isolate arm number, null values are descriptive, pitch/phase depend on ridge width and held sector, `h_R` is an SFR e-folding scale, and none of these fits estimates KTZ diffusion or metallicity-covariance parameters.

- [ ] **Step 4: Verify plots/table source and commit**

Run the comparative layout contract, explanation contract, and static parse.  Commit:

```bash
git -C /Users/Igniz/Desktop/ICRAR add -- \
  further/20260713_NGC4254_KTZ_spiral_fit.ipynb \
  further/tests/test_20260713_ngc4254_ktz_spiral_fit.py
git -C /Users/Igniz/Desktop/ICRAR diff --cached --check
git -C /Users/Igniz/Desktop/ICRAR commit -m "feat: compare m234 spiral skeletons and KTZ models"
```

---

### Task 7: Execute, inspect, accept, and document the notebook

**Files:**

- Modify: `further/20260713_NGC4254_KTZ_spiral_fit.ipynb`
- Modify: `further/docs/superpowers/plans/2026-07-15-ngc4254-m234-conditional-comparison-plan.md`
- Append job log: `/Users/Igniz/Desktop/Codex_log/2026_07_16.md`

- [x] **Step 1: Run all source and behavioral contracts before execution**

Run:

```bash
/opt/miniconda3/envs/ICRAR/bin/python -m unittest \
  tests.test_20260713_ngc4254_ktz_spiral_fit -v
```

Expected before execution: every source/behavioral test passes; only the saved-result contract fails because outputs are cleared.

- [x] **Step 2: Execute the notebook with writable Jupyter directories**

Run:

```bash
MPLCONFIGDIR=/private/tmp/ngc4254-m234-mpl \
IPYTHONDIR=/private/tmp/ngc4254-m234-ipython \
JUPYTER_CONFIG_DIR=/private/tmp/ngc4254-m234-jupyter-config \
JUPYTER_RUNTIME_DIR=/private/tmp/ngc4254-m234-jupyter-runtime \
/opt/miniconda3/envs/ICRAR/bin/python - <<'PY'
from pathlib import Path
import nbformat
from nbclient import NotebookClient

path = Path("20260713_NGC4254_KTZ_spiral_fit.ipynb")
nb = nbformat.read(path, as_version=4)
client = NotebookClient(
    nb, timeout=7200, kernel_name="python3",
    resources={"metadata": {"path": str(path.parent.resolve())}})
client.execute()
nbformat.write(nb, path)
print("NOTEBOOK_EXECUTION_PASS")
PY
```

If local port binding is sandbox-blocked, rerun the exact command with the minimum required escalation.  Do not reduce null draws, arm cases, widths, or validation thresholds to shorten execution.

- [x] **Step 3: Run the complete test suite**

Run the same `unittest` command again.  Expected: all tests, including saved outputs, pass; at least nine PNG outputs; no `RuntimeWarning`.

- [x] **Step 4: Extract and inspect every PNG**

Extract images to `/private/tmp/ngc4254_m234_figures` with a read-only `nbformat` script.  Inspect every PNG using `view_image` at original detail.  Confirm:

- all sky/deprojected skeletons wind negatively;
- all branches for `m=2`, `m=3`, and `m=4` are visible;
- the same color identifies the same arm-number case across panels;
- the 2-by-2 comparison and 3-by-2 model/residual layouts are legible;
- model maxima match their maxima tables;
- no skeleton is described as an observed arm where it visibly crosses unrelated emission; and
- axes, color bars, legends, and annotations are not clipped.

If a visual defect is found, return to the relevant task, patch via TDD, re-execute, and inspect again.

- [x] **Step 5: Record exact live results in this plan**

Append an `## Execution record` section containing:

- three selected negative pitches and `Theta0` values;
- held-sector counts and phase concentrations;
- pooled and blockwise descriptive null summaries;
- the 15-row width table summary;
- three KTZ `h_R`, gradient, `eta`, coefficient arrays, positivity minima, and costs;
- maxima counts;
- notebook runtime and warning count; and
- figure-inspection outcome.

- [x] **Step 6: Append the required job log**

Using `apply_patch` on a temporary log fragment and the minimum approved copy/write operation, append a concise entry to `/Users/Igniz/Desktop/Codex_log/2026_07_16.md`.  Include request, files changed, commits, tests, numerical checks, visual checks, unresolved conditionality, subagent roles, and confirmation that unrelated dirty files were untouched.

- [x] **Step 7: Stage only the scoped notebook, tests, design, and plan; verify and commit**

```bash
git -C /Users/Igniz/Desktop/ICRAR add -- \
  further/20260713_NGC4254_KTZ_spiral_fit.ipynb \
  further/tests/test_20260713_ngc4254_ktz_spiral_fit.py \
  further/docs/superpowers/specs/2026-07-15-ngc4254-m234-conditional-comparison-design.md \
  further/docs/superpowers/plans/2026-07-15-ngc4254-m234-conditional-comparison-plan.md
git -C /Users/Igniz/Desktop/ICRAR diff --cached --check
git -C /Users/Igniz/Desktop/ICRAR diff --cached --name-only
git -C /Users/Igniz/Desktop/ICRAR commit -m "fix: finalize conditional NGC4254 m234 comparison"
```

Expected staged names: exactly the notebook, its test file, the corrected design, and this plan. The final review required the design/test additions to keep the null statistic, fixed-m synthetic claim, and saved-output freshness contract synchronized. Never stage the unrelated dirty workspace files.

## Execution record

Execution completed on 2026-07-16 against the live `NGC4254_gas_bin_maps_further.fits` product. The final corrected in-place Jupyter run took 118.554 s from the first execute request to the final idle event. All 15 code cells have non-null execution counts; the notebook has 33 unique cell IDs, zero error outputs, zero `RuntimeWarning` or `ResourceWarning` messages, and nine embedded PNG outputs. The saved `NOTEBOOK_SOURCE_SHA256` fingerprint matches the current cell IDs, types, and sources, so outputs are not stale. The complete contract suite passed 43/43 tests.

### Conditional ridge geometries

Negative winding was imposed for every arm-number case. These rows are parallel conditional results, not an arm-number ranking.

| m | pitch (deg) | Theta0 (rad) | held sectors | phase stability | validated score | status |
|---:|---:|---:|---:|---:|---:|:---|
| 2 | -45.0 | 1.169371 | 9/12 | 0.691157 | 0.036330 | below threshold |
| 3 | -31.0 | 6.230825 | 12/12 | 0.996595 | 0.038074 | accepted |
| 4 | -45.0 | 2.356194 | 12/12 | 0.999667 | 0.022002 | accepted |

The m=2 and m=4 pitches touch the declared -45 degree search boundary. The m=2 row remains visible as a sensitivity case but fails the unchanged requirement of at least 10 valid held sectors.

### Descriptive blocked-null evidence

Each pooled row contains 32 draws. The z values are descriptive only and were not used to select an arm number.

| m | pooled mean | pooled std | descriptive z |
|---:|---:|---:|---:|
| 2 | 0.017069 | 0.012251 | 1.572242 |
| 3 | 0.015611 | 0.010654 | 2.108400 |
| 4 | 0.012512 | 0.009187 | 1.032962 |

Block means and within-block standard deviations, eight draws per `(seed, m)` pair:

| seed | m=2 | m=3 | m=4 |
|---:|:---|:---|:---|
| 4254 | 0.019911 +/- 0.014339 | 0.018183 +/- 0.009351 | 0.011127 +/- 0.008405 |
| 5254 | 0.016662 +/- 0.012329 | 0.010578 +/- 0.009148 | 0.012772 +/- 0.008435 |
| 6254 | 0.019625 +/- 0.008448 | 0.012169 +/- 0.011772 | 0.014059 +/- 0.010161 |
| 7254 | 0.012076 +/- 0.013794 | 0.021513 +/- 0.010154 | 0.012091 +/- 0.011114 |

### Ridge-width sensitivity

The five widths produce exactly 15 conditional rows. Status counts are 3/5 accepted for m=2 and 5/5 accepted for both m=3 and m=4.

| width (kpc) | m | pitch (deg) | Theta0 (rad) | held | score | accepted |
|---:|---:|---:|---:|---:|---:|:---:|
| 0.18 | 2 | -40.0 | 1.047198 | 10 | 0.003853 | yes |
| 0.18 | 3 | -32.0 | 4.084070 | 12 | 0.038255 | yes |
| 0.18 | 4 | -45.0 | 2.338741 | 12 | 0.013535 | yes |
| 0.22 | 2 | -45.0 | 1.169371 | 8 | 0.053338 | no |
| 0.22 | 3 | -30.0 | 0.366519 | 12 | 0.055132 | yes |
| 0.22 | 4 | -45.0 | 2.356194 | 12 | 0.027840 | yes |
| 0.25 | 2 | -45.0 | 1.169371 | 9 | 0.036330 | no |
| 0.25 | 3 | -31.0 | 6.230825 | 12 | 0.038074 | yes |
| 0.25 | 4 | -45.0 | 2.356194 | 12 | 0.022002 | yes |
| 0.30 | 2 | -29.0 | 0.558505 | 12 | 0.009282 | yes |
| 0.30 | 3 | -31.0 | 6.265732 | 12 | 0.047169 | yes |
| 0.30 | 4 | -44.0 | 2.670354 | 12 | 0.015318 | yes |
| 0.35 | 2 | -29.0 | 0.593412 | 12 | 0.010596 | yes |
| 0.35 | 3 | -36.0 | 4.468043 | 12 | 0.032605 | yes |
| 0.35 | 4 | -45.0 | 2.286381 | 12 | 0.007522 | yes |

### Fixed-geometry KTZ-compatible source-profile fits

| m | lambda0,0 | h_R (kpc) | gradient (dex/kpc) | eta | harmonic g | harmonic alpha (rad) | min factor | weighted log SSE | robust cost | maxima (all/enhanced) | bound note |
|---:|---:|---:|---:|---:|:---|:---|---:|---:|---:|:---:|:---|
| 2 | 0.091046 | 3.090290 | -0.140535 | 0.058503 | [1.000000, 1.500000, 1.500000] | [0.000000, -2.414377, 2.236639] | 0.794064 | 116471.586894 | 15936.011107 | 3/3 | g2 and g3 at upper bounds |
| 3 | 0.089783 | 3.119706 | -0.139210 | 0.208699 | [1.000000, 0.967117, 0.617365] | [0.000000, 0.751508, 0.543879] | 0.775121 | 115026.578819 | 15701.681331 | 3/2 | no active bounds |
| 4 | 0.088905 | 3.134492 | -0.138553 | 0.028582 | [1.000000, 1.500000, 1.500000] | [0.000000, 0.484076, 0.249337] | 0.934987 | 117348.146057 | 16008.125890 | 3/2 | g2 and g3 at upper bounds |

All fitted rate factors remain positive. The m=2 and m=4 harmonic shapes are explicitly reported as constraint-limited because both higher-order amplitudes reached their upper bounds. The notebook now displays the compact fit summary, all nine `(m, n, g_n, alpha_n)` coefficient rows, the complete 37-column audit table, and all maxima tables without hiding the requested quantities behind column elision.

### Figure inspection and acceptance

All nine saved PNGs were extracted and inspected at original detail. The two skeleton panels show the same negative-winding conditional geometries on the observed sky-plane and deprojected SFR maps; the 2-by-2 comparison and 3-by-2 model/residual layouts are legible; source-profile maxima and data-ridge provenance are distinguished; m=2 is consistently dashed/open and labeled below threshold; axes, legends, color bars, and residual-saturation disclosures are not clipped. After the null-statistic correction, eight unaffected PNGs remained byte-identical and the one changed blocked-null distribution figure was re-inspected at original detail. No arm-number winner is claimed.
