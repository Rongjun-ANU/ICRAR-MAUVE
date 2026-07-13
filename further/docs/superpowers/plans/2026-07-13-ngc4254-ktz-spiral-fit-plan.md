# NGC4254 KTZ Spiral-Geometry Fit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create and execute `20260713_NGC4254_KTZ_spiral_fit.ipynb`, a tutorial-style notebook that fits the KTZ logarithmic-spiral source-rate parameters to NGC4254's HII SFR surface-density map and produces the requested spiral-skeleton and arm-profile plots.

**Architecture:** Keep the scientific workflow self-contained in one notebook. Separate it into small, named code cells for data loading, unique-bin aggregation, deprojection, radial fitting, spiral search, harmonic fitting, constrained refinement, uncertainty estimation, and plotting; place an explanatory Markdown cell before every code cell. Use synthetic injection-recovery assertions as the regression test for the fitting mathematics, then execute the identical functions on the live NGC4254 product.

**Tech Stack:** Python 3 in the ICRAR conda environment, NumPy, SciPy, Astropy FITS/WCS/coordinates, Matplotlib, pandas, nbformat, and nbclient.

---

## File Map

- Create: `further/20260713_NGC4254_KTZ_spiral_fit.ipynb` — complete tutorial, fit, saved tables, and saved plots.
- Create: `further/tests/test_20260713_ngc4254_ktz_spiral_fit.py` — fast notebook-structure and saved-result contract checks; does not load the 3.1 GB FITS product.
- Modify: `further/docs/superpowers/plans/2026-07-13-ngc4254-ktz-spiral-fit-plan.md` — check off completed steps during execution.

No existing notebook, FITS product, or KTZ source file is modified.

### Task 1: Create the notebook contract and tutorial scaffold

**Files:**
- Create: `further/tests/test_20260713_ngc4254_ktz_spiral_fit.py`
- Create: `further/20260713_NGC4254_KTZ_spiral_fit.ipynb`

- [ ] **Step 1: Write the failing notebook-contract test**

Create the test with exact required sections, code-cell ordering, and output contracts:

```python
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "20260713_NGC4254_KTZ_spiral_fit.ipynb"


def load_notebook():
    return json.loads(NOTEBOOK.read_text(encoding="utf-8"))


def test_notebook_has_explanation_before_every_code_cell():
    nb = load_notebook()
    cells = nb["cells"]
    assert cells[0]["cell_type"] == "markdown"
    for index, cell in enumerate(cells):
        if cell["cell_type"] == "code":
            assert index > 0
            assert cells[index - 1]["cell_type"] == "markdown"
            explanation = "".join(cells[index - 1]["source"])
            assert len(explanation.split()) >= 35


def test_notebook_declares_required_sections_and_parameters():
    nb = load_notebook()
    text = "\n".join("".join(cell["source"]) for cell in nb["cells"])
    required = [
        "LOGSFR_SURFACE_DENSITY_HII",
        "BIN_ID",
        "m_arms",
        "pitch_angle",
        "Theta0",
        "h_R",
        "eta",
        "harmonic_n",
        "harmonic_g",
        "harmonic_alpha",
        "Synthetic injection-recovery",
        "Spiral skeleton",
        "Arm profile",
        "Interpretation limits",
    ]
    for item in required:
        assert item in text


def test_executed_notebook_has_saved_result_contract():
    nb = load_notebook()
    code_cells = [cell for cell in nb["cells"] if cell["cell_type"] == "code"]
    assert code_cells
    assert all(cell["execution_count"] is not None for cell in code_cells)
    text_outputs = []
    image_count = 0
    for cell in code_cells:
        for output in cell.get("outputs", []):
            if output.get("output_type") == "stream":
                text_outputs.extend(output.get("text", []))
            data = output.get("data", {})
            image_count += int("image/png" in data)
    joined = "".join(text_outputs)
    assert "SYNTHETIC_RECOVERY_PASS" in joined
    assert "REAL_FIT_COMPLETE" in joined
    assert "POSITIVITY_CHECK_PASS" in joined
    assert image_count >= 6
```

- [ ] **Step 2: Run the contract test and verify it fails**

Run:

```bash
/opt/miniconda3/envs/ICRAR/bin/python -m pytest further/tests/test_20260713_ngc4254_ktz_spiral_fit.py -q
```

Expected: fail because `20260713_NGC4254_KTZ_spiral_fit.ipynb` does not exist.

- [ ] **Step 3: Add the notebook JSON scaffold**

Use `apply_patch` to add a valid nbformat 4.5 notebook. Start with these Markdown sections, each followed by one code cell:

1. title, goal, and interpretation limits;
2. KTZ source-rate equations and parameter-identifiability convention;
3. imports and visible configuration;
4. FITS/WCS validation;
5. unique-bin construction;
6. sky-to-disc deprojection;
7. radial exponential fit;
8. spiral search;
9. harmonic fit and constrained refinement;
10. synthetic injection-recovery;
11. bootstrap uncertainty;
12. real-data fit;
13. diagnostic plots;
14. final parameter table;
15. scientific interpretation limits.

The configuration cell must contain these explicit values:

```python
from pathlib import Path
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from scipy.optimize import least_squares

ROOT = Path.cwd().resolve()
FITS_PATH = ROOT / "v3tk_v7.6.8/NGC4254/NGC4254_gas_bin_maps_further.fits"
SFR_HDU = "LOGSFR_SURFACE_DENSITY_HII"
BIN_HDU = "BIN_ID"

DISTANCE_MPC = 16.5
CENTER = SkyCoord("12h18m49.68s", "+14d25m05.52s", frame="icrs")
INCLINATION_DEG = 39.0
POSITION_ANGLE_DEG = 243.0
R_REF_KPC = 1.0
M_CANDIDATES = np.arange(1, 5, dtype=int)
PITCH_GRID_DEG = np.concatenate((np.linspace(-45.0, -5.0, 161),
                                 np.linspace(5.0, 45.0, 161)))
HARMONIC_N = np.array([1, 2, 3], dtype=int)
INNER_RADIUS_KPC = 0.5
OUTER_COVERAGE_QUANTILE = 0.98
N_BOOTSTRAP = 100
N_BOOTSTRAP_SECTORS = 12
RNG_SEED = 4254
rng = np.random.default_rng(RNG_SEED)

plt.style.use("default")
warnings.filterwarnings("default")
print(f"Input product: {FITS_PATH}")
```

- [ ] **Step 4: Validate the notebook format**

Run:

```bash
/opt/miniconda3/envs/ICRAR/bin/python -m jupyter nbconvert --validate further/20260713_NGC4254_KTZ_spiral_fit.ipynb
```

Expected: exit 0 with a valid notebook; the result-contract test may still fail because the notebook is not executed.

- [ ] **Step 5: Commit the scaffold and test**

```bash
git add further/20260713_NGC4254_KTZ_spiral_fit.ipynb further/tests/test_20260713_ngc4254_ktz_spiral_fit.py
git commit -m "test: define NGC4254 KTZ notebook contract"
```

### Task 2: Load the real FITS maps and construct independent bin records

**Files:**
- Modify: `further/20260713_NGC4254_KTZ_spiral_fit.ipynb`

- [ ] **Step 1: Add failing data-contract assertions**

At the end of the FITS-loading cell, add assertions before defining the loader:

```python
assert FITS_PATH.exists(), f"Missing input FITS file: {FITS_PATH}"
```

In the following unique-bin cell, initially call `build_bin_catalog(...)` before it is defined so an isolated execution stops with `NameError`. This establishes the red phase for the data transformation.

- [ ] **Step 2: Execute to confirm the undefined loader fails**

Run the notebook in a temporary output path with `nbclient`; expected failure: `NameError: name 'build_bin_catalog' is not defined`.

```bash
/opt/miniconda3/envs/ICRAR/bin/python - <<'PY'
import nbformat
from nbclient import NotebookClient
p = "further/20260713_NGC4254_KTZ_spiral_fit.ipynb"
nb = nbformat.read(p, as_version=4)
NotebookClient(nb, timeout=600, kernel_name="python3").execute()
PY
```

- [ ] **Step 3: Implement map validation and bin aggregation**

Replace the failing call with these helpers and then call them:

```python
def load_maps(path, sfr_hdu=SFR_HDU, bin_hdu=BIN_HDU):
    if not path.exists():
        raise FileNotFoundError(path)
    with fits.open(path, memmap=True) as hdul:
        missing = [name for name in (sfr_hdu, bin_hdu) if name not in hdul]
        if missing:
            raise KeyError(f"Missing required HDUs: {missing}")
        log_sfr = np.asarray(hdul[sfr_hdu].data, dtype=float).copy()
        bin_id = np.asarray(hdul[bin_hdu].data, dtype=float).copy()
        header = hdul[sfr_hdu].header.copy()
    if log_sfr.shape != bin_id.shape:
        raise ValueError(f"Map shape mismatch: {log_sfr.shape} versus {bin_id.shape}")
    wcs = WCS(header).celestial
    if not wcs.has_celestial:
        raise ValueError("Target HDU does not contain a celestial WCS")
    return log_sfr, bin_id, header, wcs


def build_bin_catalog(log_sfr, bin_id, wcs):
    valid = np.isfinite(log_sfr) & np.isfinite(bin_id) & (bin_id >= 0)
    yy, xx = np.nonzero(valid)
    ids = bin_id[valid].astype(np.int64)
    values = log_sfr[valid]
    unique_ids, inverse = np.unique(ids, return_inverse=True)
    area_pix = np.bincount(inverse).astype(float)
    x_centroid = np.bincount(inverse, weights=xx) / area_pix
    y_centroid = np.bincount(inverse, weights=yy) / area_pix
    sum_log_sfr = np.bincount(inverse, weights=values)
    mean_log_sfr = sum_log_sfr / area_pix
    ra_deg, dec_deg = wcs.pixel_to_world_values(x_centroid, y_centroid)
    table = pd.DataFrame({
        "bin_id": unique_ids,
        "x_pix": x_centroid,
        "y_pix": y_centroid,
        "ra_deg": ra_deg,
        "dec_deg": dec_deg,
        "area_pix": area_pix,
        "log_sfr": mean_log_sfr,
        "sfr_linear": np.power(10.0, mean_log_sfr),
    })
    if len(table) < 100:
        raise ValueError(f"Only {len(table)} valid HII bins; spiral fitting is not supported")
    return table, valid


log_sfr_map, bin_id_map, sfr_header, celestial_wcs = load_maps(FITS_PATH)
bin_table, valid_hii_pixels = build_bin_catalog(log_sfr_map, bin_id_map, celestial_wcs)
assert bin_table["bin_id"].is_unique
assert np.all(bin_table["area_pix"] > 0)
assert np.all(np.isfinite(bin_table["sfr_linear"]))
print(f"Valid HII pixels: {valid_hii_pixels.sum():,}")
print(f"Independent gas bins: {len(bin_table):,}")
```

- [ ] **Step 4: Execute through the aggregation cell**

Expected current-product values: map shape `(1190, 1012)`, roughly 681,856 finite HII pixels, and 32,552 independent valid bins. Treat a material mismatch as a live-data change to investigate, not a reason to hard-code the old counts.

- [ ] **Step 5: Commit the loader**

```bash
git add further/20260713_NGC4254_KTZ_spiral_fit.ipynb
git commit -m "feat: load NGC4254 HII SFR bins"
```

### Task 3: Implement and test deprojection and the radial exponential fit

**Files:**
- Modify: `further/20260713_NGC4254_KTZ_spiral_fit.ipynb`

- [ ] **Step 1: Add synthetic deprojection and radial-fit assertions**

Add a test cell that will fail until the helpers exist:

```python
test_ra = np.array([CENTER.ra.deg])
test_dec = np.array([CENTER.dec.deg])
test_xy = deproject_sky(test_ra, test_dec)
assert np.allclose(test_xy[0], 0.0, atol=1e-10)
assert np.allclose(test_xy[1], 0.0, atol=1e-10)

r_test = np.linspace(0.5, 8.0, 200)
lambda_true, h_true = 2.5e-2, 3.4
sfr_test = lambda_true * np.exp(-r_test / h_true)
radial_test = fit_radial_exponential(r_test, sfr_test, np.ones_like(r_test))
assert abs(radial_test["h_R"] - h_true) < 0.02
assert abs(radial_test["lambda0_0"] / lambda_true - 1.0) < 0.01
```

- [ ] **Step 2: Run and verify the test cell fails with undefined helpers**

Expected: `NameError` for `deproject_sky`.

- [ ] **Step 3: Implement the coordinate and radial helpers**

```python
ARCSEC_TO_KPC = DISTANCE_MPC * 1.0e3 / 206265.0


def deproject_sky(ra_deg, dec_deg):
    coords = SkyCoord(np.asarray(ra_deg) * u.deg, np.asarray(dec_deg) * u.deg)
    dlon, dlat = CENTER.spherical_offsets_to(coords)
    east_kpc = dlon.to_value(u.arcsec) * ARCSEC_TO_KPC
    north_kpc = dlat.to_value(u.arcsec) * ARCSEC_TO_KPC
    pa = np.deg2rad(POSITION_ANGLE_DEG)
    inc = np.deg2rad(INCLINATION_DEG)
    major = east_kpc * np.sin(pa) + north_kpc * np.cos(pa)
    minor_projected = -east_kpc * np.cos(pa) + north_kpc * np.sin(pa)
    minor = minor_projected / np.cos(inc)
    radius = np.hypot(major, minor)
    azimuth = np.arctan2(minor, major)
    return major, minor, radius, azimuth


def fit_radial_exponential(radius, sfr_linear, area_weight):
    radius = np.asarray(radius, dtype=float)
    sfr_linear = np.asarray(sfr_linear, dtype=float)
    weights = np.sqrt(np.asarray(area_weight, dtype=float))
    valid = (np.isfinite(radius) & np.isfinite(sfr_linear) &
             np.isfinite(weights) & (sfr_linear > 0) & (weights > 0))
    r = radius[valid]
    y = np.log(sfr_linear[valid])
    w = weights[valid] / np.nanmedian(weights[valid])
    initial = np.array([np.nanmedian(y), np.log(3.0)])

    def residual(theta):
        log_lambda0, log_h = theta
        return w * (y - (log_lambda0 - r / np.exp(log_h)))

    result = least_squares(residual, initial, loss="soft_l1", f_scale=0.3,
                           bounds=([-50.0, np.log(0.2)], [50.0, np.log(30.0)]))
    if not result.success:
        raise RuntimeError(result.message)
    return {
        "lambda0_0": float(np.exp(result.x[0])),
        "h_R": float(np.exp(result.x[1])),
        "cost": float(2.0 * result.cost),
        "success": bool(result.success),
    }
```

Apply `deproject_sky` to `bin_table`, define `R_out` from the configured coverage quantile, retain `INNER_RADIUS_KPC <= R <= R_out`, fit the radial model, and calculate `q = sfr_linear / radial_model - 1`.

- [ ] **Step 4: Run the synthetic checks and plot the real radial profile**

Expected: both synthetic assertions pass, all deprojected coordinates are finite, `h_R` is positive and interior to its bounds, and the real radial plot visibly follows the central tendency rather than isolated bright bins.

- [ ] **Step 5: Commit the geometry and radial fit**

```bash
git add further/20260713_NGC4254_KTZ_spiral_fit.ipynb
git commit -m "feat: deproject NGC4254 and fit radial SFR"
```

### Task 4: Implement spiral search, harmonics, and synthetic recovery

**Files:**
- Modify: `further/20260713_NGC4254_KTZ_spiral_fit.ipynb`

- [ ] **Step 1: Add a failing known-spiral injection-recovery cell**

Generate a noiseless KTZ field on the real bin coordinates with:

```python
INJECTED = {
    "m_arms": 2,
    "pitch_angle": 22.0,
    "Theta0": 0.45,
    "eta": 0.32,
    "harmonic_g": np.array([1.0, 0.35, 0.15]),
    "harmonic_alpha": np.array([0.0, 0.55, -0.30]),
}
```

Call the not-yet-defined search and harmonic fitter, then assert:

```python
assert recovered["m_arms"] == INJECTED["m_arms"]
assert abs(recovered["pitch_angle"] - INJECTED["pitch_angle"]) < 0.6
assert abs(np.angle(np.exp(1j * (recovered["Theta0"] - INJECTED["Theta0"])))) < 0.08
assert abs(recovered["eta"] - INJECTED["eta"]) < 0.03
assert np.allclose(recovered["harmonic_g"][1:], INJECTED["harmonic_g"][1:], atol=0.06)
print("SYNTHETIC_RECOVERY_PASS")
```

- [ ] **Step 2: Execute and verify failure for undefined spiral functions**

Expected: `NameError` for `search_spiral_geometry`.

- [ ] **Step 3: Implement the KTZ phase, profile, search, and fit**

```python
def wrap_angle(angle):
    return np.angle(np.exp(1j * np.asarray(angle)))


def spiral_phase(radius, azimuth, m_arms, pitch_angle, theta0=0.0):
    pitch = np.deg2rad(pitch_angle)
    return ((m_arms / np.tan(pitch)) * np.log(radius / R_REF_KPC)
            - m_arms * azimuth + theta0)


def arm_profile(theta, harmonic_g, harmonic_alpha):
    theta = np.asarray(theta, dtype=float)
    value = np.zeros_like(theta)
    for n, g, alpha in zip(HARMONIC_N, harmonic_g, harmonic_alpha):
        value += g * np.cos(n * theta + alpha)
    return value


def search_spiral_geometry(radius, azimuth, q, area_weight):
    w = np.asarray(area_weight, dtype=float)
    w = w / w.sum()
    rows = []
    for m in M_CANDIDATES:
        for pitch in PITCH_GRID_DEG:
            base = spiral_phase(radius, azimuth, int(m), float(pitch), theta0=0.0)
            coefficient = np.sum(w * q * np.exp(-1j * base))
            rows.append((int(m), float(pitch), float(np.abs(coefficient)),
                         float(np.angle(coefficient))))
    search = pd.DataFrame(rows, columns=["m_arms", "pitch_angle", "power", "Theta0"])
    best = search.loc[search["power"].idxmax()].to_dict()
    best["m_arms"] = int(best["m_arms"])
    return best, search


def fit_harmonic_profile(radius, azimuth, q, area_weight, m_arms,
                         pitch_angle, theta0_initial):
    sqrt_w = np.sqrt(np.asarray(area_weight, dtype=float))
    sqrt_w /= np.nanmedian(sqrt_w)

    def unpack(parameters):
        theta0, log_eta, g2, g3, alpha2, alpha3 = parameters
        return (theta0, np.exp(log_eta), np.array([1.0, g2, g3]),
                np.array([0.0, alpha2, alpha3]))

    def residual(parameters):
        theta0, eta, harmonic_g, harmonic_alpha = unpack(parameters)
        theta = spiral_phase(radius, azimuth, m_arms, pitch_angle, theta0)
        model_q = eta * arm_profile(theta, harmonic_g, harmonic_alpha)
        profile_grid = arm_profile(np.linspace(0.0, 2.0 * np.pi, 2048,
                                               endpoint=False),
                                   harmonic_g, harmonic_alpha)
        floor = 1.0 + eta * np.min(profile_grid)
        penalty = np.array([min(0.0, floor - 1.0e-3) * 1.0e4])
        return np.concatenate((sqrt_w * (q - model_q), penalty))

    initial = np.array([theta0_initial, np.log(0.2), 0.2, 0.1, 0.0, 0.0])
    lower = np.array([-4.0 * np.pi, np.log(1.0e-4), -1.5, -1.5, -np.pi, -np.pi])
    upper = np.array([4.0 * np.pi, np.log(2.0), 1.5, 1.5, np.pi, np.pi])
    result = least_squares(residual, initial, bounds=(lower, upper),
                           loss="soft_l1", f_scale=0.25, max_nfev=3000)
    theta0, eta, harmonic_g, harmonic_alpha = unpack(result.x)
    return {
        "Theta0": float(theta0 % (2.0 * np.pi)),
        "eta": float(eta),
        "harmonic_g": harmonic_g,
        "harmonic_alpha": harmonic_alpha,
        "cost": float(2.0 * result.cost),
        "success": bool(result.success),
        "message": result.message,
    }
```

Add this wrapper so the synthetic and real fits use the identical path:

```python
def recover_spiral_parameters(radius, azimuth, q, area_weight):
    discrete, search_table = search_spiral_geometry(radius, azimuth, q, area_weight)
    harmonic = fit_harmonic_profile(
        radius=radius,
        azimuth=azimuth,
        q=q,
        area_weight=area_weight,
        m_arms=discrete["m_arms"],
        pitch_angle=discrete["pitch_angle"],
        theta0_initial=discrete["Theta0"],
    )
    if not harmonic["success"]:
        raise RuntimeError(f"Harmonic fit failed: {harmonic['message']}")
    recovered = {
        "m_arms": int(discrete["m_arms"]),
        "pitch_angle": float(discrete["pitch_angle"]),
        "Theta0": harmonic["Theta0"],
        "eta": harmonic["eta"],
        "harmonic_n": HARMONIC_N.copy(),
        "harmonic_g": harmonic["harmonic_g"],
        "harmonic_alpha": harmonic["harmonic_alpha"],
        "search_power": float(discrete["power"]),
        "harmonic_cost": harmonic["cost"],
    }
    return recovered, search_table
```

- [ ] **Step 4: Run synthetic recovery and verify red-green behaviour**

Expected: the inserted known spiral is recovered within the declared tolerances and prints `SYNTHETIC_RECOVERY_PASS`. If the signed-pitch convention returns the geometrically equivalent conjugate solution, make the convention canonical and test that canonical form; do not simply widen tolerances.

- [ ] **Step 5: Commit the recovered spiral model**

```bash
git add further/20260713_NGC4254_KTZ_spiral_fit.ipynb
git commit -m "feat: recover KTZ spiral geometry and harmonics"
```

### Task 5: Add constrained full-field refinement and spatial bootstrap

**Files:**
- Modify: `further/20260713_NGC4254_KTZ_spiral_fit.ipynb`

- [ ] **Step 1: Add failing refinement invariants**

For the synthetic field, assert that joint refinement does not worsen the weighted objective and that the model intensity stays positive:

```python
refined_synthetic = refine_full_model(synthetic_table, recovered, radial_test)
assert refined_synthetic["objective"] <= refined_synthetic["initial_objective"] + 1.0e-8
assert refined_synthetic["minimum_rate_factor"] > 0.0
```

- [ ] **Step 2: Run and verify failure for undefined refinement**

Expected: `NameError: name 'refine_full_model' is not defined`.

- [ ] **Step 3: Implement bounded full-model refinement**

The parameter vector is

```python
[log_lambda0, log_h_R, pitch_angle, Theta0, log_eta,
 g2, g3, alpha2, alpha3]
```

Keep `m_arms` fixed per candidate. Evaluate

```python
radial = np.exp(log_lambda0 - radius / np.exp(log_h_R))
theta = spiral_phase(radius, azimuth, m_arms, pitch_angle, Theta0)
rate_factor = 1.0 + np.exp(log_eta) * arm_profile(theta, g, alpha)
model = radial * rate_factor
```

Implement the complete refinement as:

```python
def refine_full_model(table, spiral_guess, radial_guess):
    radius = table["radius_kpc"].to_numpy(float)
    azimuth = table["azimuth_rad"].to_numpy(float)
    observed = table["sfr_linear"].to_numpy(float)
    sqrt_w = np.sqrt(table["area_pix"].to_numpy(float))
    sqrt_w /= np.nanmedian(sqrt_w)
    scale = np.nanmedian(observed)
    m_arms = int(spiral_guess["m_arms"])
    g0 = np.asarray(spiral_guess["harmonic_g"], dtype=float)
    a0 = np.asarray(spiral_guess["harmonic_alpha"], dtype=float)

    initial = np.array([
        np.log(radial_guess["lambda0_0"]),
        np.log(radial_guess["h_R"]),
        spiral_guess["pitch_angle"],
        spiral_guess["Theta0"],
        np.log(spiral_guess["eta"]),
        g0[1], g0[2], a0[1], a0[2],
    ])
    if initial[2] < 0:
        pitch_lower, pitch_upper = -60.0, -2.0
    else:
        pitch_lower, pitch_upper = 2.0, 60.0
    lower = np.array([-50.0, np.log(0.2), pitch_lower, -4.0*np.pi,
                      np.log(1.0e-4), -1.5, -1.5, -np.pi, -np.pi])
    upper = np.array([50.0, np.log(30.0), pitch_upper, 4.0*np.pi,
                      np.log(2.0), 1.5, 1.5, np.pi, np.pi])
    dense_theta = np.linspace(0.0, 2.0*np.pi, 4096, endpoint=False)

    def unpack(parameters):
        log_lambda0, log_h, pitch, theta0, log_eta, g2, g3, alpha2, alpha3 = parameters
        return {
            "lambda0_0": np.exp(log_lambda0),
            "h_R": np.exp(log_h),
            "pitch_angle": pitch,
            "Theta0": theta0,
            "eta": np.exp(log_eta),
            "harmonic_g": np.array([1.0, g2, g3]),
            "harmonic_alpha": np.array([0.0, alpha2, alpha3]),
        }

    def evaluate(parameters):
        p = unpack(parameters)
        theta = spiral_phase(radius, azimuth, m_arms,
                             p["pitch_angle"], p["Theta0"])
        rate_factor = 1.0 + p["eta"] * arm_profile(
            theta, p["harmonic_g"], p["harmonic_alpha"])
        radial = p["lambda0_0"] * np.exp(-radius / p["h_R"])
        dense_factor = 1.0 + p["eta"] * arm_profile(
            dense_theta, p["harmonic_g"], p["harmonic_alpha"])
        return radial * rate_factor, float(np.min(dense_factor))

    def residual(parameters):
        model, minimum_factor = evaluate(parameters)
        data_residual = sqrt_w * (observed - model) / scale
        positivity_penalty = np.array([min(0.0, minimum_factor - 1.0e-3) * 1.0e4])
        return np.concatenate((data_residual, positivity_penalty))

    initial_objective = float(np.sum(residual(initial)**2))
    result = least_squares(residual, initial, bounds=(lower, upper),
                           loss="soft_l1", f_scale=0.25, max_nfev=5000)
    model, minimum_factor = evaluate(result.x)
    refined = unpack(result.x)
    refined.update({
        "m_arms": m_arms,
        "harmonic_n": HARMONIC_N.copy(),
        "objective": float(np.sum(residual(result.x)**2)),
        "initial_objective": initial_objective,
        "minimum_rate_factor": minimum_factor,
        "success": bool(result.success and np.isfinite(model).all() and minimum_factor > 0.0),
        "message": result.message,
        "on_boundary": bool(np.any(np.isclose(result.x, lower, rtol=0.0, atol=1.0e-5)) or
                            np.any(np.isclose(result.x, upper, rtol=0.0, atol=1.0e-5))),
    })
    if not refined["success"]:
        raise RuntimeError(f"Full refinement failed: {refined['message']}; "
                           f"minimum rate factor={minimum_factor}")
    return refined
```

- [ ] **Step 4: Implement sector bootstrap**

Assign each bin to one of `N_BOOTSTRAP_SECTORS` equal azimuth sectors. For each bootstrap realization, sample sectors with replacement, refit the radial and spiral models, and store a tidy result row. Use this complete implementation:

```python
def fit_one_table(table):
    radial = fit_radial_exponential(
        table["radius_kpc"], table["sfr_linear"], table["area_pix"])
    radial_model = radial["lambda0_0"] * np.exp(-table["radius_kpc"].to_numpy() /
                                                radial["h_R"])
    q = table["sfr_linear"].to_numpy() / radial_model - 1.0
    spiral, search = recover_spiral_parameters(
        table["radius_kpc"].to_numpy(), table["azimuth_rad"].to_numpy(),
        q, table["area_pix"].to_numpy())
    refined = refine_full_model(table, spiral, radial)
    return radial, spiral, refined, search


def run_sector_bootstrap(table, n_bootstrap=N_BOOTSTRAP):
    work = table.copy()
    wrapped_phi = np.mod(work["azimuth_rad"].to_numpy(), 2.0*np.pi)
    sector = np.floor(wrapped_phi / (2.0*np.pi) * N_BOOTSTRAP_SECTORS).astype(int)
    sector = np.clip(sector, 0, N_BOOTSTRAP_SECTORS - 1)
    work["bootstrap_sector"] = sector
    seed_sequence = np.random.SeedSequence(RNG_SEED)
    child_seeds = seed_sequence.spawn(n_bootstrap)
    rows = []
    for iteration, child_seed in enumerate(child_seeds):
        local_rng = np.random.default_rng(child_seed)
        selected = local_rng.integers(0, N_BOOTSTRAP_SECTORS,
                                      size=N_BOOTSTRAP_SECTORS)
        pieces = [work.loc[work["bootstrap_sector"] == value].copy()
                  for value in selected]
        sample = pd.concat(pieces, ignore_index=True)
        try:
            radial, spiral, refined, _ = fit_one_table(sample)
            rows.append({
                "iteration": iteration,
                "success": True,
                "error": "",
                "m_arms": refined["m_arms"],
                "pitch_angle": refined["pitch_angle"],
                "Theta0": refined["Theta0"] % (2.0*np.pi),
                "lambda0_0": refined["lambda0_0"],
                "h_R": refined["h_R"],
                "eta": refined["eta"],
                "g2": refined["harmonic_g"][1],
                "g3": refined["harmonic_g"][2],
                "alpha2": refined["harmonic_alpha"][1],
                "alpha3": refined["harmonic_alpha"][2],
                "objective": refined["objective"],
            })
        except Exception as exc:
            rows.append({"iteration": iteration, "success": False,
                         "error": f"{type(exc).__name__}: {exc}"})
    results = pd.DataFrame(rows)
    valid_fraction = float(results["success"].mean())
    if valid_fraction < 0.8:
        warnings.warn(f"Only {valid_fraction:.1%} of bootstrap fits succeeded")
    return results


def summarize_bootstrap(results):
    valid = results.loc[results["success"]].copy()
    names = ["pitch_angle", "Theta0", "lambda0_0", "h_R", "eta",
             "g2", "g3", "alpha2", "alpha3"]
    summary = pd.DataFrame({
        name: {
            "p16": valid[name].quantile(0.16),
            "median": valid[name].median(),
            "p84": valid[name].quantile(0.84),
            "n_valid": valid[name].notna().sum(),
        }
        for name in names
    }).T
    m_frequency = valid["m_arms"].value_counts(normalize=True).sort_index()
    return summary, m_frequency
```

Failed fits remain in the returned table with their error messages instead of being dropped silently.

Summarize median, 16th percentile, 84th percentile, and valid-fit count for every continuous parameter. Summarize selection frequency for each `m_arms` if the bootstrap reruns the discrete search. Print a warning if fewer than 80 percent of bootstrap fits are valid.

- [ ] **Step 5: Run synthetic refinement and a 10-bootstrap smoke test**

Temporarily set `N_BOOTSTRAP = 10`; expected: synthetic refinement passes, all intensities remain positive, and at least eight bootstrap fits succeed. Restore `N_BOOTSTRAP = 100` before the final execution.

- [ ] **Step 6: Commit refinement and uncertainty support**

```bash
git add further/20260713_NGC4254_KTZ_spiral_fit.ipynb
git commit -m "feat: refine NGC4254 spiral fit with bootstrap"
```

### Task 6: Fit NGC4254 and create the requested plots and parameter table

**Files:**
- Modify: `further/20260713_NGC4254_KTZ_spiral_fit.ipynb`

- [ ] **Step 1: Run the real-data fit**

Apply the radial mask, fit all four `m` candidates, rank them with the same weighted objective and AIC-like score, refine the winning candidate, and run the spatial bootstrap. Print `REAL_FIT_COMPLETE` only after checking finite parameters and successful optimization.

- [ ] **Step 2: Add the observed/deprojected data panels**

Create these displayed figures:

1. observed `LOGSFR_SURFACE_DENSITY_HII` with WCS axes, centre, major axis, and minor axis;
2. deprojected bin-centroid map colored by log SFR;
3. radial profile plus the fitted exponential;
4. spiral-power curves or heatmap for `m = 1..4`.

Use percentile-based color limits and label all units and conventions.

- [ ] **Step 3: Add the spiral skeleton overlays**

Solve the fitted ridge condition `Theta = 2*pi*k` for each arm:

```python
phi_skeleton = ((m_arms / np.tan(np.deg2rad(pitch_angle)))
                * np.log(radius_grid / R_REF_KPC)
                + Theta0 - 2.0 * np.pi * k) / m_arms
```

Display both:

- the fitted skeleton in deprojected Cartesian kpc coordinates, matching the reference notebook style;
- the skeleton projected back into the observed WCS map, overlaid on the HII log-SFR field.

The projected overlay is the visual acceptance check. If the skeleton systematically misses the main high-SFR ridges, report model inadequacy rather than tuning plot limits to hide it.

- [ ] **Step 4: Add the KTZ arm-profile plot**

For `Theta` from 0 to `2*pi`, plot every fitted harmonic as a dashed curve, the total `h(Theta)` as a black curve, and area-weighted phase-folded SFR residual summaries with uncertainty bars. Add a second axis or panel showing `1 + eta*h(Theta)` and its zero line. Print `POSITIVITY_CHECK_PASS` only if its sampled minimum is strictly positive.

- [ ] **Step 5: Add model and residual maps**

Evaluate the fitted rate at every valid bin centroid and display observed, model, and residual values with shared limits where meaningful. Include residual versus radius and residual versus phase summaries so coherent failures are visible.

- [ ] **Step 6: Add the final parameter table**

Build a pandas table with columns:

```python
["reference_variable", "estimate", "p16", "p84", "unit_or_convention", "status"]
```

Include `R_in`, `R_out`, `m_arms`, `pitch_angle`, `Theta0`, `lambda0_0`, `h_R`, `eta`, `harmonic_n`, `harmonic_g`, and `harmonic_alpha`. Add separate rows marking `kappa`, `x0`, `t_star`, and clustering parameters as `not fitted from SFR morphology`.

- [ ] **Step 7: Write the concluding Markdown cell**

State what the best global KTZ model captures, where coherent residuals remain, how stable `m` and pitch are across bootstrap sectors, and why these parameters can be used as a fixed source-field template but do not determine diffusion, enrichment time, clustering, or pattern speed.

- [ ] **Step 8: Commit the real-data result cells**

```bash
git add further/20260713_NGC4254_KTZ_spiral_fit.ipynb
git commit -m "feat: visualize NGC4254 KTZ spiral fit"
```

### Task 7: Execute, inspect, and verify the final notebook

**Files:**
- Modify: `further/20260713_NGC4254_KTZ_spiral_fit.ipynb`
- Modify: `further/docs/superpowers/plans/2026-07-13-ngc4254-ktz-spiral-fit-plan.md`

- [ ] **Step 1: Execute in process against the real FITS file**

Use the ICRAR kernel and writable runtime directories:

```bash
MPLCONFIGDIR=/private/tmp/mpl-ngc4254-ktz \
IPYTHONDIR=/private/tmp/ipython-ngc4254-ktz \
JUPYTER_CONFIG_DIR=/private/tmp/jupyter-config-ngc4254-ktz \
JUPYTER_RUNTIME_DIR=/private/tmp/jupyter-runtime-ngc4254-ktz \
/opt/miniconda3/envs/ICRAR/bin/python - <<'PY'
from pathlib import Path
import nbformat
from nbclient import NotebookClient

path = Path("further/20260713_NGC4254_KTZ_spiral_fit.ipynb")
notebook = nbformat.read(path, as_version=4)
client = NotebookClient(
    notebook,
    timeout=3600,
    kernel_name="python3",
    resources={"metadata": {"path": str(path.parent)}},
)
client.execute()
nbformat.write(notebook, path)
PY
```

Expected: exit 0 and saved outputs in every code cell.

- [ ] **Step 2: Run the contract tests**

```bash
/opt/miniconda3/envs/ICRAR/bin/python -m pytest further/tests/test_20260713_ngc4254_ktz_spiral_fit.py -q
```

Expected: all tests pass, including `SYNTHETIC_RECOVERY_PASS`, `REAL_FIT_COMPLETE`, `POSITIVITY_CHECK_PASS`, and at least six saved PNG outputs.

- [ ] **Step 3: Validate notebook JSON and scan runtime outputs**

```bash
/opt/miniconda3/envs/ICRAR/bin/python -m jupyter nbconvert --validate further/20260713_NGC4254_KTZ_spiral_fit.ipynb
/opt/miniconda3/envs/ICRAR/bin/python - <<'PY'
import json
from pathlib import Path
p = Path("further/20260713_NGC4254_KTZ_spiral_fit.ipynb")
nb = json.loads(p.read_text())
errors = [o for c in nb["cells"] for o in c.get("outputs", [])
          if o.get("output_type") == "error"]
print({"cells": len(nb["cells"]), "errors": len(errors)})
assert not errors
PY
```

Expected: valid notebook and zero saved error outputs.

- [ ] **Step 4: Render the notebook for visual QA**

```bash
/opt/miniconda3/envs/ICRAR/bin/python -m jupyter nbconvert \
  --to html \
  --output-dir /private/tmp/ngc4254-ktz-render \
  further/20260713_NGC4254_KTZ_spiral_fit.ipynb
```

Inspect the rendered output for clipped labels, blank panels, misleading color scales, geometry-axis errors, skeleton/ridge alignment, and readable Markdown equations. If HTML inspection is insufficient, extract the saved PNG outputs to `/private/tmp` and inspect them individually.

- [ ] **Step 5: Re-read the approved specification as a requirements checklist**

Confirm every required parameter, plot, caveat, error check, and synthetic recovery gate in `further/docs/superpowers/specs/2026-07-13-ngc4254-ktz-spiral-fit-design.md` is present. Record any model inadequacy or unstable parameter in the notebook conclusion rather than calling the fit definitive.

- [ ] **Step 6: Mark this plan complete and inspect the final diff**

Run:

```bash
git status --short
git diff --stat HEAD
git diff --check
```

Ensure only the new notebook, its focused contract test, and this plan's checkbox updates are in scope. Do not stage the existing unrelated `KTZ_validation.ipynb`, `20260713_KTZ.md`, or `KTZ_Theory_z_fluctuation_III/` changes.

- [ ] **Step 7: Commit the executed notebook and verified plan**

```bash
git add further/20260713_NGC4254_KTZ_spiral_fit.ipynb \
        further/tests/test_20260713_ngc4254_ktz_spiral_fit.py \
        further/docs/superpowers/plans/2026-07-13-ngc4254-ktz-spiral-fit-plan.md
git commit -m "test: verify NGC4254 KTZ spiral notebook"
```

Report the exact fit parameters, bootstrap stability, execution result, test result, and any visual/scientific limitation in the final handoff.
