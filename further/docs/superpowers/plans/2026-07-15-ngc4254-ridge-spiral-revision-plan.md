# NGC4254 Ridge-Based Spiral-Fit Revision Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the global-intensity geometry selector in `20260713_NGC4254_KTZ_spiral_fit.ipynb` with a mask-aware, ridge-based `m=1..6` signed-pitch search, then fit and display the KTZ-compatible harmonic source profile without hiding secondary maxima.

**Architecture:** Preserve the self-contained tutorial notebook but separate its inference into two stages. Stage 1 uses every valid HII pixel in deprojected log-polar space to select `m_arms`, signed pitch, and `Theta0` from a radius-corrected, local centre-versus-flanks ridge response with azimuth-sector validation; Stage 2 freezes that geometry and fits the radial and harmonic source-rate parameters to independent gas bins. Fast JSON contract tests protect notebook structure and saved outputs, while in-notebook synthetic assertions test the numerical functions on both winding signs and all arm multiplicities. A declared physical-width scan exposes the real-data `m`/pitch sensitivity instead of hiding it.

**Tech Stack:** Python 3 from `/opt/miniconda3/envs/ICRAR`, NumPy, SciPy (`optimize`, `ndimage`), Astropy FITS/WCS/coordinates, pandas, Matplotlib, nbformat, nbclient, and built-in `unittest`.

---

## File Map

- Modify: `further/20260713_NGC4254_KTZ_spiral_fit.ipynb` — geometry provenance, finite-thickness validation, ridge search, fixed-geometry KTZ fit, saved tables, and figures.
- Modify: `further/tests/test_20260713_ngc4254_ktz_spiral_fit.py` — fast source/output contracts for the revised workflow.
- Modify: `further/docs/superpowers/plans/2026-07-15-ngc4254-ridge-spiral-revision-plan.md` — check off completed steps and append the final execution record.

Do not modify `KTZ_validation.ipynb`, `SFR+Z.py`, the target FITS product, `20260713_KTZ.md`, or `KTZ_Theory_z_fluctuation_III/`. The live target notebook is already modified in the working tree; it is the source of truth and must not be replaced by `HEAD`.

## Editing and verification conventions

- Use `/opt/miniconda3/envs/ICRAR/bin/python` for every Python check.
- Before each edit, inspect `git status --short` and stage only the exact files named by that task.
- Patch notebook cells by stable Markdown heading or cell id, never by assumed cell number.
- Use `apply_patch` to create or update a temporary mechanical notebook-revision script at `/private/tmp/revise_ngc4254_ridge_notebook.py`. The script must load the live notebook, replace only named cells, preserve every unrelated cell, clear outputs from each changed cell onward, and write the notebook with nbformat.
- Use this helper in that temporary script:

```python
from pathlib import Path
import nbformat

NOTEBOOK = Path("/Users/Igniz/Desktop/ICRAR/further/20260713_NGC4254_KTZ_spiral_fit.ipynb")


def markdown_index(nb, heading):
    matches = [i for i, cell in enumerate(nb.cells)
               if cell.cell_type == "markdown" and heading in cell.source]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one heading {heading!r}; found {matches}")
    return matches[0]


def replace_code_after(nb, heading, source):
    md_index = markdown_index(nb, heading)
    code_index = md_index + 1
    if code_index >= len(nb.cells) or nb.cells[code_index].cell_type != "code":
        raise RuntimeError(f"Heading {heading!r} is not followed by a code cell")
    nb.cells[code_index].source = source.rstrip() + "\n"
    for cell in nb.cells[code_index:]:
        if cell.cell_type == "code":
            cell.outputs = []
            cell.execution_count = None
    return code_index


nb = nbformat.read(NOTEBOOK, as_version=4)
# Task-specific replacements or insertions go here.
nbformat.write(nb, NOTEBOOK)
```

- Source-only commits may contain cleared outputs. The saved-result contract is intentionally run only after the final full execution.

### Task 1: Lock the revised contract and geometry provenance

**Files:**
- Modify: `further/tests/test_20260713_ngc4254_ktz_spiral_fit.py`
- Modify: `further/20260713_NGC4254_KTZ_spiral_fit.ipynb`

- [ ] **Step 1: Add failing geometry and workflow contract tests**

Add these helpers below `load_notebook()`:

```python
def notebook_source():
    nb = load_notebook()
    return "\n".join("".join(cell["source"]) for cell in nb["cells"])


def notebook_output_text():
    pieces = []
    for cell in load_notebook()["cells"]:
        if cell["cell_type"] != "code":
            continue
        for output in cell.get("outputs", []):
            if output.get("output_type") == "stream":
                value = output.get("text", "")
                pieces.extend([value] if isinstance(value, str) else value)
    return "".join(pieces)
```

Add these methods to `NotebookContractTests`:

```python
def test_brown_geometry_and_finite_thickness_contract(self):
    text = notebook_source()
    for required in [
        "Brown2021Table1.txt",
        "finite_thickness_axis_ratio",
        "Q0 = 0.2",
        "SURFACE_DENSITY_ALREADY_CORRECTED = True",
        "APPLY_BA_CORRECTION_IN_NOTEBOOK = False",
        "UPSTREAM_SFR_LOG",
        "UPSTREAM_BA_PASS",
        'frame=FK5(equinox="J2000")',
        "WCS_REFERENCE_NOT_CENTER_PASS",
        "BROWN_GEOMETRY_PASS",
        "BA_FACTOR_PASS",
    ]:
        self.assertIn(required, text)
    self.assertNotIn("log_sfr_map + np.log10(B_OVER_A)", text)
    self.assertNotIn("log_sfr_map * B_OVER_A", text)

def test_ridge_geometry_precedes_ktz_profile(self):
    text = notebook_source()
    for required in [
        "M_CANDIDATES = np.arange(1, 7",
        "build_log_polar_contrast",
        "held_out_ridge_search",
        "RIDGE_SYNTHETIC_SIGN_PASS",
        "RIDGE_SYNTHETIC_M_PASS",
        "RIDGE_SYNTHETIC_PHASE_PASS",
        "RIDGE_M_NORMALIZATION_PASS",
        "RIDGE_BRANCH_NORMALIZATION_PASS",
        "RIDGE_NULL_NORMALIZATION_PASS",
        "RIDGE_WIDTH_SENSITIVITY_COMPLETE",
        "RIDGE_REAL_FIT_COMPLETE",
        "RIDGE_SECTOR_STABILITY_COMPLETE",
        "KTZ_PROFILE_FIT_COMPLETE",
    ]:
        self.assertIn(required, text)
    self.assertLess(text.index("RIDGE_REAL_FIT_COMPLETE"),
                    text.index("KTZ_PROFILE_FIT_COMPLETE"))

def test_all_harmonic_maxima_and_approved_layout_are_declared(self):
    text = notebook_source()
    for required in [
        "enumerate_profile_maxima",
        "HARMONIC_MAXIMA_COMPLETE",
        "plt.subplots(2, 2",
        "observed sky-plane",
        "observed deprojected",
        "fitted KTZ-compatible source model",
        "observed - model residual",
        "SYNTHETIC_WINDING_SIGN_PASS",
        "FINAL_PARAMETER_TABLE_COMPLETE",
        "DEPROJECTED_SFR_SKELETON_PASS",
    ]:
        self.assertIn(required, text)

def test_legacy_global_selector_and_bootstrap_are_removed(self):
    text = notebook_source()
    for forbidden in [
        "search_spiral_geometry",
        "fit_harmonic_profile",
        "recover_spiral_parameters",
        "refine_full_model",
        "compare_candidate_modes",
        "opposite_winding_fit",
        "run_sector_bootstrap",
        "bootstrap_summary",
        "valid_bootstrap",
    ]:
        self.assertNotIn(forbidden, text)
```

In `test_notebook_declares_required_sections_and_parameters`, replace the stale required strings `Synthetic injection-recovery` and `Spiral skeleton` with `Synthetic ridge-recovery` and `data-derived skeleton` respectively.

Replace `test_deprojection_is_explained_and_ridge_uses_harmonic_peak` completely with:

```python
def test_deprojection_and_data_model_ridge_semantics(self):
    text = notebook_source()
    for required in [
        "east = major*sin(PA) - projected_minor*cos(PA)",
        "minor = projected_minor/cos(inclination)",
        "DEPROJECTION_ROUNDTRIP_PASS",
        "SYNTHETIC_WINDING_SIGN_PASS",
        "data_ridge_curves",
        "enhanced_model_curves",
        "data-derived skeleton",
        "source-model maxima",
    ]:
        self.assertIn(required, text)
    self.assertNotIn("opposite_winding_fit", text)
    self.assertNotIn("RIDGE_PHASE_CHECK_PASS", text)
```

These edits are required because the revised observed panels deliberately use the independently selected `Theta=0` data ridge, not the global maximum of the later harmonic source profile.

Update `test_executed_notebook_has_saved_result_contract` so its required markers are:

```python
for marker in [
    "BROWN_GEOMETRY_PASS",
    "BA_FACTOR_PASS",
    "UPSTREAM_BA_PASS",
    "WCS_REFERENCE_NOT_CENTER_PASS",
    "FIT_DOMAIN_ALL_VALID_HII_PASS",
    "RIDGE_PIXEL_DOMAIN_PASS",
    "RIDGE_SYNTHETIC_SIGN_PASS",
    "RIDGE_SYNTHETIC_M_PASS",
    "RIDGE_SYNTHETIC_PHASE_PASS",
    "RIDGE_M_NORMALIZATION_PASS",
    "RIDGE_BRANCH_NORMALIZATION_PASS",
    "RIDGE_NULL_NORMALIZATION_PASS",
    "RIDGE_WIDTH_SENSITIVITY_COMPLETE",
    "RIDGE_REAL_FIT_COMPLETE",
    "RIDGE_SECTOR_STABILITY_COMPLETE",
    "KTZ_PROFILE_FIT_COMPLETE",
    "HARMONIC_MAXIMA_COMPLETE",
    "DEPROJECTION_ROUNDTRIP_PASS",
    "SYNTHETIC_WINDING_SIGN_PASS",
    "DEPROJECTED_SFR_SKELETON_PASS",
    "FINAL_PARAMETER_TABLE_COMPLETE",
    "POSITIVITY_CHECK_PASS",
]:
    self.assertIn(marker, joined)
self.assertGreaterEqual(image_count, 7)
```

- [ ] **Step 2: Run only the new source contracts and verify RED**

Run:

```bash
/opt/miniconda3/envs/ICRAR/bin/python -m unittest \
  tests.test_20260713_ngc4254_ktz_spiral_fit.NotebookContractTests.test_brown_geometry_and_finite_thickness_contract \
  tests.test_20260713_ngc4254_ktz_spiral_fit.NotebookContractTests.test_ridge_geometry_precedes_ktz_profile \
  tests.test_20260713_ngc4254_ktz_spiral_fit.NotebookContractTests.test_all_harmonic_maxima_and_approved_layout_are_declared -v
```

Expected: three assertion failures because the live notebook still declares `m=1..4`, has no ridge-search functions, and has no equation-3 provenance markers.

- [ ] **Step 3: Replace the configuration and geometry explanation**

In the temporary notebook-revision script, update the Markdown before the configuration cell to say that Brown et al. supplies the optical centre, inclination, and directed receding-side PA; equation 3 of [Huang et al. (2026)](https://academic.oup.com/mnras/article/549/3/stag1019/8698768) supplies the already-applied surface-density factor; and coordinate deprojection still uses `cos(i)`. State explicitly that the FITS WCS handles sky orientation, so the array is not manually flipped north-up/east-left before the tangent-plane rotation.

In that same Markdown, replace the stale `m=1..4` statement with the exact
`m=1..6, both signed-pitch branches` search contract. In the bin-catalog
Markdown, replace “bootstrap resampling” with “independent-bin KTZ fitting and
leave-one-sector-out morphology stability”; the old full-model bootstrap will
no longer exist. In the radial-background cell, delete the obsolete `q_search`
initializer, its clipping rule, and its comments rather than leaving dead
global-Fourier setup beside the new ridge workflow.

Replace the configuration code with the existing imports/configuration plus these exact additions and substitutions:

```python
import re
from astropy.coordinates import FK5, SkyCoord
from scipy.ndimage import gaussian_filter, gaussian_filter1d
from scipy.optimize import least_squares, minimize_scalar

BROWN_TABLE = ROOT / "Brown2021Table1.txt"
UPSTREAM_SFR_LOG = ROOT / "sfr_logs/NGC4254.log"
Q0 = 0.2
SURFACE_DENSITY_ALREADY_CORRECTED = True
APPLY_BA_CORRECTION_IN_NOTEBOOK = False
M_CANDIDATES = np.arange(1, 7, dtype=int)
RIDGE_PITCH_GRID_DEG = np.concatenate((np.arange(-45.0, -4.9, 1.0),
                                       np.arange(5.0, 45.1, 1.0)))
RIDGE_CORE_WIDTH_KPC = 0.25
RIDGE_WIDTH_SENSITIVITY_KPC = np.array([0.18, 0.22, 0.25, 0.30, 0.35])
RIDGE_BROAD_RATIO = 3.0
RIDGE_N_PHASE = 360
RIDGE_N_SECTORS = 12
RIDGE_MIN_HELD_OUT_SECTORS = 10
RIDGE_SHORTLIST_PER_FAMILY = 5
RIDGE_N_NULL = 8
LOGPOLAR_N_U = 120
LOGPOLAR_N_PHI = 360
LOGPOLAR_N_RADIAL_BANDS = 24
LOGPOLAR_SMOOTH_SIGMA = (0.8, 1.2)
LOGPOLAR_AZIMUTH_BROAD_SIGMA_BINS = 30.0
```

Add this complete geometry loader in the same cell, then replace the hard-coded `CENTER`, inclination, and PA assignments with its result:

```python
def load_brown_geometry(path, galaxy="NGC 4254"):
    if not path.exists():
        raise FileNotFoundError(path)
    rows = [line.rstrip("\n") for line in path.read_text(encoding="utf-8").splitlines()
            if line.startswith(galaxy)]
    if len(rows) != 1:
        raise RuntimeError(f"Expected one {galaxy} row in {path}; found {len(rows)}")
    fields = rows[0].split("\t")
    if len(fields) < 6:
        raise ValueError(f"Malformed Brown table row: {rows[0]}")
    ra_match = re.fullmatch(r"(\d+)\^h(\d+)\^m(\d+)\.s(\d+)", fields[1])
    dec_match = re.fullmatch(r"([+-]?\d+)deg(\d+)'(\d+)\.\"(\d+)", fields[2])
    if ra_match is None or dec_match is None:
        raise ValueError(f"Unrecognized Brown coordinates: {fields[1:3]}")
    hh, mm, ss, sf = ra_match.groups()
    dd, dm, ds, df = dec_match.groups()
    sign = "-" if dd.startswith("-") else "+"
    dd_abs = dd.lstrip("+-")
    coordinate = SkyCoord(
        f"{hh}h{mm}m{ss}.{sf}s",
        f"{sign}{dd_abs}d{dm}m{ds}.{df}s",
        frame=FK5(equinox="J2000"),
    )
    return {
        "center": coordinate,
        "inclination_deg": float(fields[4]),
        "position_angle_deg": float(fields[5]),
        "source_row": rows[0],
    }


def finite_thickness_axis_ratio(inclination_deg, q0=Q0):
    cosine = np.cos(np.deg2rad(float(inclination_deg)))
    return float(np.sqrt((1.0 - q0**2) * cosine**2 + q0**2))


BROWN_GEOMETRY = load_brown_geometry(BROWN_TABLE)
CENTER = BROWN_GEOMETRY["center"]
INCLINATION_DEG = BROWN_GEOMETRY["inclination_deg"]
POSITION_ANGLE_DEG = BROWN_GEOMETRY["position_angle_deg"]
B_OVER_A = finite_thickness_axis_ratio(INCLINATION_DEG)
EXPECTED_CENTER = SkyCoord("12h18m49.68s", "+14d25m05.52s",
                           frame=FK5(equinox="J2000"))
assert CENTER.separation(EXPECTED_CENTER).to_value(u.arcsec) < 1e-6
assert np.isclose(INCLINATION_DEG, 39.0)
assert np.isclose(POSITION_ANGLE_DEG % 360.0, 243.0)
assert np.isclose(B_OVER_A, 0.787, atol=5e-4)
assert SURFACE_DENSITY_ALREADY_CORRECTED
assert not APPLY_BA_CORRECTION_IN_NOTEBOOK
upstream_text = UPSTREAM_SFR_LOG.read_text(encoding="utf-8")
upstream_matches = re.findall(
    r"Inclination correction ENABLED: applying b/a = ([0-9.]+)",
    upstream_text,
)
if not upstream_matches:
    raise RuntimeError(f"No enabled b/a correction record in {UPSTREAM_SFR_LOG}")
assert np.isclose(float(upstream_matches[-1]), B_OVER_A, atol=5e-4)
print(f"Brown catalog row: {BROWN_GEOMETRY['source_row']}")
print(f"Adopted FK5 J2000 centre: {CENTER.to_string('hmsdms')}")
print(f"Inclination={INCLINATION_DEG:.1f} deg; directed PA={POSITION_ANGLE_DEG:.1f} deg east of north")
print("BROWN_GEOMETRY_PASS")
print(f"Equation-3 finite-thickness factor b/a={B_OVER_A:.4f}; already applied upstream")
print("BA_FACTOR_PASS")
print(f"Confirmed upstream correction record: b/a={float(upstream_matches[-1]):.3f}")
print("UPSTREAM_BA_PASS")
```

In the FITS-loading code cell, immediately after `celestial_wcs` is created, verify and explain that the projection reference coordinate is not the Brown optical centre:

```python
wcs_reference = SkyCoord(
    celestial_wcs.wcs.crval[0] * u.deg,
    celestial_wcs.wcs.crval[1] * u.deg,
    frame=FK5(equinox="J2000"),
)
wcs_center_separation_arcsec = float(
    wcs_reference.separation(CENTER).to_value(u.arcsec))
assert wcs_center_separation_arcsec > 1.0
print(f"FITS CRVAL is a WCS projection reference, not the Brown optical centre; "
      f"separation={wcs_center_separation_arcsec:.3f} arcsec")
print("WCS_REFERENCE_NOT_CENTER_PASS")
```

The expected live value is about `26.54 arcsec`; the inequality avoids treating that measured value as an immutable FITS constant while still making a silent centre substitution fatal.

- [ ] **Step 4: Run the geometry contract and static notebook validation**

Run:

```bash
/opt/miniconda3/envs/ICRAR/bin/python -m unittest \
  tests.test_20260713_ngc4254_ktz_spiral_fit.NotebookContractTests.test_brown_geometry_and_finite_thickness_contract -v
/opt/miniconda3/envs/ICRAR/bin/python - <<'PY'
import ast, nbformat
nb = nbformat.read("20260713_NGC4254_KTZ_spiral_fit.ipynb", as_version=4)
for cell in nb.cells:
    if cell.cell_type == "code":
        ast.parse(cell.source)
print("STATIC_NOTEBOOK_PASS")
PY
```

Expected: the geometry contract and static parse pass. Do not run the saved-result contract because outputs are intentionally cleared.

- [ ] **Step 5: Commit the geometry contract**

```bash
git -C /Users/Igniz/Desktop/ICRAR add -- \
  further/20260713_NGC4254_KTZ_spiral_fit.ipynb \
  further/tests/test_20260713_ngc4254_ktz_spiral_fit.py
git -C /Users/Igniz/Desktop/ICRAR diff --cached --check
git -C /Users/Igniz/Desktop/ICRAR commit -m "feat: validate NGC4254 geometry provenance"
```

### Task 2: Build the all-pixel deprojected log-polar local-ridge map

**Files:**
- Modify: `further/20260713_NGC4254_KTZ_spiral_fit.ipynb`
- Modify: `further/tests/test_20260713_ngc4254_ktz_spiral_fit.py`

- [ ] **Step 1: Add a failing all-pixel morphology contract**

Add this test method:

```python
def test_ridge_map_uses_all_valid_hii_pixels_and_mask_normalization(self):
    text = notebook_source()
    for required in [
        "build_hii_pixel_geometry",
        "build_log_polar_contrast",
        "gaussian_filter(weighted_sum",
        "gaussian_filter(coverage",
        "mode=(\"nearest\", \"wrap\")",
        "np.quantile(u",
        "signed_residual",
        "local_ridge",
        "LOGPOLAR_AZIMUTH_BROAD_SIGMA_BINS",
        "gradient_dex_per_kpc",
        "RIDGE_PIXEL_DOMAIN_PASS",
    ]:
        self.assertIn(required, text)
    self.assertNotIn("0.05 * np.nanmax(denominator)", text)
    self.assertNotIn("np.maximum(observed_log - background_log, 0.0)", text)
```

Run it and expect failure because the two functions and marker are absent:

```bash
/opt/miniconda3/envs/ICRAR/bin/python -m unittest \
  tests.test_20260713_ngc4254_ktz_spiral_fit.NotebookContractTests.test_ridge_map_uses_all_valid_hii_pixels_and_mask_normalization -v
```

- [ ] **Step 2: Make the deprojection frame-consistent and add the all-pixel helper**

In the existing `deproject_sky` function, update the docstring to say `FK5 J2000` and replace its coordinate construction with the exact frame used by `CENTER` and the FITS WCS:

```python
coords = SkyCoord(np.asarray(ra_deg) * u.deg,
                  np.asarray(dec_deg) * u.deg,
                  frame=FK5(equinox="J2000"))
```

Do not leave the current default-ICRS `SkyCoord(...)` call in place: Astropy's `spherical_offsets_to` rejects non-matching FK5/ICRS frames. Keep the tangent-plane rotation equations unchanged.

Then insert a new Markdown/code pair after the existing bin-centroid deprojection cell. The Markdown must explain that pixels represent spatial support for morphology while unique bins remain the independent fit records. Use this code:

```python
def build_hii_pixel_geometry(log_sfr, valid_mask, wcs):
    yy, xx = np.nonzero(valid_mask)
    ra_deg, dec_deg = wcs.pixel_to_world_values(xx, yy)
    major, minor, radius, azimuth = deproject_sky(ra_deg, dec_deg)
    table = pd.DataFrame({
        "x_pix": xx.astype(float),
        "y_pix": yy.astype(float),
        "x_disc_kpc": major,
        "y_disc_kpc": minor,
        "radius_kpc": radius,
        "azimuth_rad": azimuth,
        "log_sfr": log_sfr[valid_mask],
    })
    keep = np.isfinite(table).all(axis=1) & (table["radius_kpc"] > 0)
    result = table.loc[keep].reset_index(drop=True)
    if len(result) != int(valid_mask.sum()):
        raise RuntimeError(f"Lost {int(valid_mask.sum()) - len(result)} valid HII pixels")
    return result


hii_pixels = build_hii_pixel_geometry(log_sfr_map, valid_hii_pixels, celestial_wcs)
assert len(hii_pixels) == int(valid_hii_pixels.sum())
RIDGE_R_IN_KPC = float(hii_pixels["radius_kpc"].min())
RIDGE_R_OUT_KPC = float(hii_pixels["radius_kpc"].max())
RIDGE_PIVOT_RADIUS_KPC = float(np.exp(
    np.median(np.log(hii_pixels["radius_kpc"].to_numpy(float)))))
print(f"Ridge morphology pixels: {len(hii_pixels):,} / {int(valid_hii_pixels.sum()):,}")
print(f"Ridge pixel radius: {RIDGE_R_IN_KPC:.4f} to {RIDGE_R_OUT_KPC:.4f} kpc; "
      f"phase pivot={RIDGE_PIVOT_RADIUS_KPC:.3f} kpc")
print("RIDGE_PIXEL_DOMAIN_PASS")
```

- [ ] **Step 3: Explain the radial background and add the mask-normalized log-polar builder**

Update the radial-background Markdown to say that `h_R` is the e-folding length of the axisymmetric SFR source field: after increasing radius by `h_R`, the background falls by a factor of `e`. It is not an arm width, pitch, optical radius, or automatically a stellar-disc scale length. Immediately after `radial_fit`, report its equivalent dex-per-kpc slope, then use this complete helper and invocation:

```python
radial_fit["gradient_dex_per_kpc"] = float(
    -1.0 / (np.log(10.0) * radial_fit["h_R"]))
print(f"Axisymmetric SFR background h_R={radial_fit['h_R']:.3f} kpc; "
      f"gradient={radial_fit['gradient_dex_per_kpc']:.4f} dex/kpc")


def build_log_polar_contrast(
        pixel_table, radial_parameters,
        n_u=LOGPOLAR_N_U, n_phi=LOGPOLAR_N_PHI,
        smooth_sigma=LOGPOLAR_SMOOTH_SIGMA,
        broad_sigma_phi=LOGPOLAR_AZIMUTH_BROAD_SIGMA_BINS):
    radius = pixel_table["radius_kpc"].to_numpy(float)
    azimuth = pixel_table["azimuth_rad"].to_numpy(float)
    observed_log = pixel_table["log_sfr"].to_numpy(float)
    background_log = np.log10(radial_parameters["lambda0_0"]) \
        - radius / (np.log(10.0) * radial_parameters["h_R"])
    signed_residual = observed_log - background_log
    u = np.log(radius / R_REF_KPC)

    # Equal-count radial rows retain the central HII pixels and prevent the
    # much larger outer annuli from setting a global coverage threshold.
    u_edges = np.quantile(u, np.linspace(0.0, 1.0, n_u + 1))
    if np.any(np.diff(u_edges) <= 0):
        raise RuntimeError("Quantile log-radius edges are not strictly increasing")
    u_edges[0] = np.nextafter(u_edges[0], -np.inf)
    u_edges[-1] = np.nextafter(u_edges[-1], np.inf)
    phi_edges = np.linspace(-np.pi, np.pi, n_phi + 1)
    weighted_sum, _, _ = np.histogram2d(u, azimuth, bins=(u_edges, phi_edges),
                                        weights=signed_residual)
    coverage, _, _ = np.histogram2d(u, azimuth, bins=(u_edges, phi_edges))
    numerator = gaussian_filter(weighted_sum, smooth_sigma,
                                mode=("nearest", "wrap"))
    denominator = gaussian_filter(coverage, smooth_sigma,
                                  mode=("nearest", "wrap"))
    valid = denominator > 0.05
    radial_residual = np.full_like(numerator, np.nan, dtype=float)
    radial_residual[valid] = numerator[valid] / denominator[valid]

    # Remove diffuse/lopsided azimuthal structure without treating masks as
    # zero. Negative flanks remain informative and are deliberately not clipped.
    broad_numerator = gaussian_filter1d(
        np.where(valid, radial_residual * denominator, 0.0),
        broad_sigma_phi, axis=1, mode="wrap")
    broad_denominator = gaussian_filter1d(
        np.where(valid, denominator, 0.0),
        broad_sigma_phi, axis=1, mode="wrap")
    broad_azimuthal = np.divide(
        broad_numerator, broad_denominator,
        out=np.zeros_like(broad_numerator), where=broad_denominator > 0)
    local_ridge = np.full_like(radial_residual, np.nan)
    local_ridge[valid] = radial_residual[valid] - broad_azimuthal[valid]

    row_index = np.clip(np.searchsorted(u_edges, u, side="right") - 1,
                        0, n_u - 1)
    row_count = np.bincount(row_index, minlength=n_u)
    u_centres = np.divide(
        np.bincount(row_index, weights=u, minlength=n_u), row_count,
        out=np.full(n_u, np.nan), where=row_count > 0)
    phi_centres = 0.5 * (phi_edges[:-1] + phi_edges[1:])
    return {
        "radial_residual": radial_residual,
        "local_ridge": local_ridge,
        "coverage": denominator,
        "valid": valid,
        "u": u_centres,
        "phi": phi_centres,
        "u_edges": u_edges,
        "phi_edges": phi_edges,
    }


logpolar = build_log_polar_contrast(hii_pixels, radial_fit)
assert np.isfinite(logpolar["local_ridge"][logpolar["valid"]]).all()
assert len(hii_pixels) == int(valid_hii_pixels.sum())

fig, axes = plt.subplots(1, 2, figsize=(14.0, 5.2), constrained_layout=True,
                         sharex=True, sharey=True)
for ax, field, title in zip(
        axes, ["radial_residual", "local_ridge"],
        ["signed residual after radial background",
         "local ridge field after broad-azimuth removal"]):
    image = ax.pcolormesh(
        logpolar["phi_edges"], logpolar["u_edges"], logpolar[field],
        shading="auto", cmap="coolwarm")
    ax.set(xlabel=r"disc azimuth $\phi$ (rad)",
           ylabel=r"$\ln(R/R_{\rm ref})$", title=title)
    fig.colorbar(image, ax=ax, label="residual (dex)")
plt.show()
```

The Markdown above this cell must explicitly record the live falsification that
motivated the implementation: a `0.05 * max(coverage)` cut retained only about
14% of log-polar cells and removed every row below roughly `2.5 kpc`, despite
the pixels having entered the histogram.  The quantile rows and absolute local
coverage floor are the fix; they do not create a centre hole.

- [ ] **Step 4: Verify GREEN and commit**

Run the new contract, all existing source contracts except the saved-result test, and the static parse. Expected: pass.

```bash
/opt/miniconda3/envs/ICRAR/bin/python -m unittest \
  tests.test_20260713_ngc4254_ktz_spiral_fit.NotebookContractTests.test_ridge_map_uses_all_valid_hii_pixels_and_mask_normalization \
  tests.test_20260713_ngc4254_ktz_spiral_fit.NotebookContractTests.test_fit_domain_uses_every_valid_hii_bin -v
git -C /Users/Igniz/Desktop/ICRAR add -- \
  further/20260713_NGC4254_KTZ_spiral_fit.ipynb \
  further/tests/test_20260713_ngc4254_ktz_spiral_fit.py
git -C /Users/Igniz/Desktop/ICRAR diff --cached --check
git -C /Users/Igniz/Desktop/ICRAR commit -m "feat: build mask-aware log-polar SFR ridge field"
```

### Task 3: Implement and validate the held-out ridge search

**Files:**
- Modify: `further/20260713_NGC4254_KTZ_spiral_fit.ipynb`
- Modify: `further/tests/test_20260713_ngc4254_ktz_spiral_fit.py`

- [ ] **Step 1: Add a failing numerical-search contract**

Add:

```python
def test_ridge_search_has_local_flank_control_and_sector_holdout(self):
    text = notebook_source()
    for required in [
        "phase_sector_row_histograms",
        "variable_width_ridge_response",
        "aggregate_radial_response",
        "held_out_ridge_search",
        "sigma_phase = m_arms * core_width_kpc /",
        "narrow_mean - broad_mean",
        "held_out_score",
        'table["valid_held_out"] >= RIDGE_MIN_HELD_OUT_SECTORS',
        "phase_stability",
        "positive_radial_fraction",
        "scramble_logpolar_rows",
        "build_m_null_calibration",
        "null_z",
        "ridge_width_null_draws",
        "branch_normalization = \"mean_not_sum\"",
    ]:
        self.assertIn(required, text)
    self.assertNotIn("RIDGE_N_PHASE // 2", text)
```

Run it and expect failure because the scoring implementation is absent.

- [ ] **Step 2: Preserve the phase primitives and replace only the global Fourier selector**

Replace the Section 8 Markdown first. Explain the log-polar phase, the
radius-dependent conversion from a constant physical ridge width to phase
width, the narrow-minus-broad local-flank control, mean-not-sum branch
normalization, and azimuth-sector holdouts. State that `m_arms`, signed pitch,
and `Theta0` are frozen by the validated local-ridge score. Explicitly show

```text
d_perp approximately R * abs(sin(pitch)) * abs(Theta) / m
sigma_Theta(R) = m * width_kpc / (R * abs(sin(pitch))).
```

Also state that the rejected implementation used one phase width evaluated at
`R_ref`; at `R=5` and `10 kpc` its nominal `0.35 kpc` corridor had broadened to
about `1.75` and `3.5 kpc`, respectively, and therefore rewarded global
lopsidedness. Remove the current claims about a weighted complex Fourier
coefficient and a later optimizer refining pitch or phase.

In the code cell under `## 8. Define the KTZ spiral phase and search arm number and signed pitch` (currently id `16219bba`), retain `wrap_angle`, `spiral_phase`, and `arm_profile` unchanged because every later source-model calculation uses them. Delete `search_spiral_geometry`, its real-data invocation, and the Fourier-power figure. Append these functions after the three retained primitives. The old Fourier method may be reintroduced only as a clearly named secondary diagnostic; it must not select geometry.

```python
def phase_sector_row_histograms(
        logpolar_map, m_arms, pitch_angle, field="local_ridge",
        n_phase=RIDGE_N_PHASE, n_sectors=RIDGE_N_SECTORS):
    uu, pp = np.meshgrid(logpolar_map["u"], logpolar_map["phi"], indexing="ij")
    valid = logpolar_map["valid"] & np.isfinite(logpolar_map[field])
    n_u = len(logpolar_map["u"])
    u = uu[valid]
    phi = pp[valid]
    row = np.broadcast_to(np.arange(n_u)[:, None], valid.shape)[valid]
    ridge_value = logpolar_map[field][valid]
    coverage = logpolar_map["coverage"][valid]
    base_phase = np.mod((m_arms / np.tan(np.deg2rad(pitch_angle))) * u
                        - m_arms * phi, 2.0 * np.pi)
    phase_index = np.floor(base_phase / (2.0 * np.pi) * n_phase).astype(int) % n_phase
    sector_index = np.floor((phi + np.pi) / (2.0 * np.pi) * n_sectors).astype(int)
    sector_index = np.clip(sector_index, 0, n_sectors - 1)
    joint_index = (sector_index * n_u + row) * n_phase + phase_index
    size = n_sectors * n_u * n_phase
    weighted = np.bincount(joint_index, weights=coverage * ridge_value,
                           minlength=size).reshape(n_sectors, n_u, n_phase)
    support = np.bincount(joint_index, weights=coverage,
                          minlength=size).reshape(n_sectors, n_u, n_phase)
    return weighted, support


def variable_width_ridge_response(
        weighted_hist, support_hist, u_centres, m_arms, pitch_angle,
        core_width_kpc=RIDGE_CORE_WIDTH_KPC,
        broad_ratio=RIDGE_BROAD_RATIO,
        n_radial_bands=LOGPOLAR_N_RADIAL_BANDS):
    """Local centre-minus-flanks response at a constant physical width."""
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
    sine_pitch = abs(np.sin(np.deg2rad(pitch_angle)))
    sigma_phase = m_arms * core_width_kpc / (radius_band * sine_pitch)
    sigma_bins = np.clip(
        sigma_phase / (2.0 * np.pi) * n_phase, 0.65, n_phase / 10.0)

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
    return narrow_mean - broad_mean, narrow_mean, broad_mean, narrow_support


def aggregate_radial_response(response, support, min_radial_fraction=0.35):
    """Give each approximately equal-count radial band one vote."""
    good = np.isfinite(response) & (support > 0.05)
    count = good.sum(axis=0)
    total = np.where(good, response, 0.0).sum(axis=0)
    mean_response = np.divide(
        total, count, out=np.full(response.shape[1], np.nan), where=count > 0)
    positive_count = (good & (response > 0)).sum(axis=0)
    positive_fraction = np.divide(
        positive_count, count,
        out=np.full(response.shape[1], np.nan), where=count > 0)
    median_response = np.array([
        np.median(response[good[:, index], index]) if count[index] else np.nan
        for index in range(response.shape[1])
    ])
    radial_fraction = count / response.shape[0]
    combined = mean_response * positive_fraction + 0.5 * median_response
    combined[radial_fraction < min_radial_fraction] = np.nan
    return {
        "score": combined,
        "mean": mean_response,
        "median": median_response,
        "positive_fraction": positive_fraction,
        "radial_fraction": radial_fraction,
    }


def held_out_ridge_search(logpolar_map, m_candidates=M_CANDIDATES,
                          pitch_grid=RIDGE_PITCH_GRID_DEG,
                          core_width_kpc=RIDGE_CORE_WIDTH_KPC,
                          allow_no_acceptance=False):
    records = []
    histogram_cache = {}
    branch_normalization = "mean_not_sum"
    phase_values = np.linspace(0.0, 2.0 * np.pi, RIDGE_N_PHASE, endpoint=False)
    for m_arms in m_candidates:
        for pitch_angle in pitch_grid:
            weighted, support = phase_sector_row_histograms(
                logpolar_map, int(m_arms), float(pitch_angle))
            histogram_cache[(int(m_arms), float(pitch_angle))] = (weighted, support)
            total_weighted = weighted.sum(axis=0)
            total_support = support.sum(axis=0)
            response, narrow, broad, response_support = variable_width_ridge_response(
                total_weighted, total_support, logpolar_map["u"],
                int(m_arms), float(pitch_angle), core_width_kpc)
            full = aggregate_radial_response(response, response_support)
            if not np.isfinite(full["score"]).any():
                continue
            full_index = int(np.nanargmax(full["score"]))
            records.append({
                "m_arms": int(m_arms),
                "pitch_angle": float(pitch_angle),
                "winding_sign": int(np.sign(pitch_angle)),
                "Theta0": float((-phase_values[full_index]) % (2.0*np.pi)),
                "core_width_kpc": float(core_width_kpc),
                "full_ridge_score": float(full["score"][full_index]),
                "narrow_mean": float(np.nanmean(narrow[:, full_index])),
                "broad_flank_mean": float(np.nanmean(broad[:, full_index])),
                "positive_radial_fraction": float(
                    full["positive_fraction"][full_index]),
                "radial_coverage": float(full["radial_fraction"][full_index]),
                "branch_normalization": branch_normalization,
                "held_out_score": np.nan,
                "held_out_score_std": np.nan,
                "phase_stability": np.nan,
                "valid_held_out": 0,
                "validated_score": np.nan,
            })
    table = pd.DataFrame(records)
    if table.empty:
        raise RuntimeError("No ridge candidate had finite full-map support")

    # Screen five pitches per m/sign on the full map, then determine phase on
    # the other 11 sectors and score it only in the omitted sector. The notebook
    # calls this screened sector validation, not an unbiased posterior score.
    shortlist = (table.sort_values("full_ridge_score", ascending=False)
                 .groupby(["m_arms", "winding_sign"], group_keys=False)
                 .head(RIDGE_SHORTLIST_PER_FAMILY).index)
    for record_index in shortlist:
        m_arms = int(table.at[record_index, "m_arms"])
        pitch_angle = float(table.at[record_index, "pitch_angle"])
        weighted, support = histogram_cache[(m_arms, pitch_angle)]
        total_weighted = weighted.sum(axis=0)
        total_support = support.sum(axis=0)
        held_out_scores = []
        held_out_phases = []
        for sector in range(RIDGE_N_SECTORS):
            train_response, _, _, train_support = variable_width_ridge_response(
                total_weighted - weighted[sector],
                total_support - support[sector], logpolar_map["u"],
                m_arms, pitch_angle, core_width_kpc)
            train = aggregate_radial_response(train_response, train_support)
            if not np.isfinite(train["score"]).any():
                continue
            train_index = int(np.nanargmax(train["score"]))
            test_response, _, _, test_support = variable_width_ridge_response(
                weighted[sector], support[sector], logpolar_map["u"],
                m_arms, pitch_angle, core_width_kpc)
            test = aggregate_radial_response(
                test_response, test_support, min_radial_fraction=0.08)
            if np.isfinite(test["score"][train_index]):
                held_out_scores.append(float(test["score"][train_index]))
                held_out_phases.append(float(
                    (-phase_values[train_index]) % (2.0 * np.pi)))
        valid_held_out = len(held_out_scores)
        if valid_held_out:
            phase_stability = float(np.abs(np.mean(
                np.exp(1j * np.asarray(held_out_phases)))))
            held_out_score = float(np.median(held_out_scores))
            table.at[record_index, "held_out_score"] = held_out_score
            table.at[record_index, "held_out_score_std"] = float(
                np.std(held_out_scores))
            table.at[record_index, "phase_stability"] = phase_stability
            table.at[record_index, "valid_held_out"] = valid_held_out
            table.at[record_index, "validated_score"] = (
                max(held_out_score, 0.0)
                * (valid_held_out / RIDGE_N_SECTORS)
                * phase_stability)

    accepted = table.loc[
        (table["valid_held_out"] >= RIDGE_MIN_HELD_OUT_SECTORS)
        & (table["validated_score"] > 0)].copy()
    if accepted.empty and not allow_no_acceptance:
        raise RuntimeError("No ridge candidate passed sector-support acceptance")
    ranked = table.sort_values(
        ["validated_score", "full_ridge_score"], ascending=False,
        na_position="last").reset_index(drop=True)
    if accepted.empty:
        return None, ranked
    accepted = accepted.sort_values(
        ["validated_score", "full_ridge_score"], ascending=False)
    return accepted.iloc[0].to_dict(), ranked
```

Add `RIDGE_SHORTLIST_PER_FAMILY = 5` beside the other configuration constants.
The shortlist is an explicit runtime compromise: every `m`, sign, and pitch is
screened by the full-map local-ridge response, while sector validation is run on
the five best pitches within each `m`/sign family.  The final Markdown must not
call `validated_score` a posterior probability or a fully nested
cross-validation statistic.

- [ ] **Step 3: Replace the legacy synthetic cell before touching the real-data selector**

Replace the existing `## 10. Synthetic injection-recovery test` Markdown and code cell (currently id `c9d26c41`) with a section named `Synthetic ridge-recovery for arm number, winding, and phase`, followed by this generator and assertions. This removes every call to the deleted `recover_spiral_parameters` and `refine_full_model` functions:

```python
def synthetic_sfr_pixels(m_arms, pitch_angle, theta0, seed,
                         n_pixels=80000, mask_fraction=0.25,
                         arm_amplitude_dex=0.48,
                         ridge_sigma_kpc=0.22, radial_h_kpc=3.5):
    """Masked SFR pixels with a known radial decline and logarithmic ridges."""
    local_rng = np.random.default_rng(seed)
    u = local_rng.uniform(np.log(0.35), np.log(11.0), n_pixels)
    azimuth = local_rng.uniform(-np.pi, np.pi, n_pixels)
    radius = R_REF_KPC * np.exp(u)
    theta = ((m_arms / np.tan(np.deg2rad(pitch_angle))) * u
             - m_arms * azimuth + theta0)
    wrapped = wrap_angle(theta)
    perpendicular_kpc = (np.abs(wrapped) * radius
                         * abs(np.sin(np.deg2rad(pitch_angle))) / m_arms)
    ridge_dex = arm_amplitude_dex * np.exp(
        -0.5 * (perpendicular_kpc / ridge_sigma_kpc)**2)
    lambda0_0 = 0.05
    background_log = (np.log10(lambda0_0)
                      - radius / (np.log(10.0) * radial_h_kpc))
    observed_log = background_log + ridge_dex
    observed_log += local_rng.normal(0.0, 0.035, size=n_pixels)
    keep = local_rng.random(n_pixels) >= mask_fraction
    keep &= ~((azimuth > 0.35) & (azimuth < 0.70) & (radius > 5.5))
    keep &= ~((azimuth > -2.20) & (azimuth < -1.95) & (radius < 2.2))
    pixels = pd.DataFrame({
        "radius_kpc": radius[keep],
        "azimuth_rad": azimuth[keep],
        "log_sfr": observed_log[keep],
    })
    radial_parameters = {"lambda0_0": lambda0_0, "h_R": radial_h_kpc}
    return pixels, radial_parameters


def synthetic_logpolar(m_arms, pitch_angle, theta0, seed, **kwargs):
    pixels, radial_parameters = synthetic_sfr_pixels(
        m_arms, pitch_angle, theta0, seed, **kwargs)
    return build_log_polar_contrast(pixels, radial_parameters)


synthetic_cases = [(1, -18.0), (2, 21.0), (3, -24.0),
                   (4, 27.0), (5, -30.0), (6, 33.0)]
theta0_true = 0.7
synthetic_recoveries = []
for case_index, (m_true, pitch_true) in enumerate(synthetic_cases):
    synthetic = synthetic_logpolar(
        m_true, pitch_true, theta0_true, RNG_SEED + case_index)
    recovered, _ = held_out_ridge_search(
        synthetic, m_candidates=np.arange(1, 7),
        pitch_grid=RIDGE_PITCH_GRID_DEG,
    )
    synthetic_recoveries.append(recovered)
    assert recovered["m_arms"] == m_true, (m_true, recovered)
    assert np.sign(recovered["pitch_angle"]) == np.sign(pitch_true), recovered
    assert abs(recovered["pitch_angle"] - pitch_true) <= 2.0, recovered
    phase_error = abs(wrap_angle(recovered["Theta0"] - theta0_true))
    assert phase_error <= np.deg2rad(5.0), (m_true, phase_error, recovered)
print("RIDGE_SYNTHETIC_SIGN_PASS")
print("RIDGE_SYNTHETIC_M_PASS")
print("RIDGE_SYNTHETIC_PHASE_PASS")

assert synthetic_recoveries[0]["m_arms"] == 1
print("RIDGE_M_NORMALIZATION_PASS")

# Scaling both numerator and support must leave a support-normalized response
# unchanged; this is the direct branch-count normalization invariant.
test_weighted, test_support = phase_sector_row_histograms(
    synthetic, 3, -24.0)
response_1, _, _, _ = variable_width_ridge_response(
    test_weighted.sum(axis=0), test_support.sum(axis=0),
    synthetic["u"], 3, -24.0)
response_6, _, _, _ = variable_width_ridge_response(
    6.0 * test_weighted.sum(axis=0), 6.0 * test_support.sum(axis=0),
    synthetic["u"], 3, -24.0)
assert np.allclose(response_1, response_6, equal_nan=True, atol=1e-12)
print("RIDGE_BRANCH_NORMALIZATION_PASS")
```

Run the notebook only through this synthetic cell in a disposable in-memory
execution. Expected: all five synthetic markers print. These cases exercise
both signs, every `m`, phase recovery, a 25% random mask plus two coherent
missing regions, an exponential radial decline, and the branch-count
normalization invariant. The real-map row-scrambled null calibration is added
in the next step. If any injected case fails, fix the score normalization or
phase convention; do not relax the assertions around a wrong result.

- [ ] **Step 4: Replace the legacy real-data comparison with the ridge search**

Replace the Section 11 Markdown and title with `## 11. Select the real NGC4254 ridge geometry, then fit its source profile`. Explain the two-stage inference: all valid HII pixels select geometry from the signed, local ridge field and screened sector validation; independent bins then fit the full un-clipped source field without changing that geometry. Explain that row-wise circular scrambling preserves each radial row's residual distribution, mask, and coverage while destroying a coherent spiral slope, and that the resulting per-`m` null distribution corrects the remaining look-elsewhere bias. Remove the current statements about clipped fractional-residual Fourier initialization, eight `m=1..4` refinements, full-model objective ranking, and base-`m=1` winning as a lopsided mode.

Call `null_z` an empirical morphology-ranking z-score. With only eight
scrambles it is not a Gaussian detection significance, posterior odds, or a
formal uncertainty on `m` or pitch. The same eight-scramble procedure must be
repeated independently at every declared physical ridge width; raw scores from
different widths must not be compared as though they shared one scale.

Replace the entire code cell below it (currently id `da64c8d1`) with the block below. At this task boundary it selects geometry only; Task 4 will append the fixed-geometry KTZ profile call. This removal is deliberate: no references to `compare_candidate_modes`, `candidate_fit_table`, `candidate_models`, or `opposite_winding_fit` may remain.

```python
def scramble_logpolar_rows(logpolar_map, seed):
    """Destroy cross-radius spiral coherence but preserve every row and mask."""
    local_rng = np.random.default_rng(seed)
    result = {
        key: (value.copy() if isinstance(value, np.ndarray) else value)
        for key, value in logpolar_map.items()
    }
    for row_index in range(len(result["u"])):
        shift = int(local_rng.integers(0, len(result["phi"])))
        for key in ["radial_residual", "local_ridge", "coverage", "valid"]:
            result[key][row_index] = np.roll(
                result[key][row_index], shift)
    return result


def build_m_null_calibration(
        logpolar_map, core_width_kpc=RIDGE_CORE_WIDTH_KPC,
        n_null=RIDGE_N_NULL):
    rows = []
    for null_index in range(n_null):
        scrambled = scramble_logpolar_rows(
            logpolar_map, RNG_SEED + 1000 + null_index)
        _, null_table = held_out_ridge_search(
            scrambled, core_width_kpc=core_width_kpc,
            allow_no_acceptance=True)
        for m_arms in M_CANDIDATES:
            family = null_table.loc[
                (null_table["m_arms"] == m_arms)
                & (null_table["valid_held_out"] >= RIDGE_MIN_HELD_OUT_SECTORS)
                & np.isfinite(null_table["validated_score"])]
            best_null_score = (float(family["validated_score"].max())
                               if len(family) else 0.0)
            rows.append({
                "null_index": null_index,
                "m_arms": int(m_arms),
                "core_width_kpc": float(core_width_kpc),
                "best_null_score": best_null_score,
            })
    draws = pd.DataFrame(rows)
    summary = (draws.groupby("m_arms")["best_null_score"]
               .agg(null_mean="mean", null_std="std", null_count="count")
               .reset_index())
    if not (summary["null_count"] == n_null).all():
        raise RuntimeError("Incomplete per-m null calibration")
    summary["null_std_floor"] = np.maximum(summary["null_std"], 1e-6)
    return draws, summary


def apply_m_null_calibration(candidate_table, null_summary):
    calibrated = candidate_table.merge(
        null_summary, on="m_arms", how="left", validate="many_to_one")
    calibrated["null_z"] = (
        (calibrated["validated_score"] - calibrated["null_mean"])
        / calibrated["null_std_floor"])
    accepted = calibrated.loc[
        (calibrated["valid_held_out"] >= RIDGE_MIN_HELD_OUT_SECTORS)
        & (calibrated["validated_score"] > 0)
        & np.isfinite(calibrated["null_z"])].copy()
    if accepted.empty:
        raise RuntimeError(
            "No real ridge candidate survived support and null calibration")
    accepted = accepted.sort_values(
        ["null_z", "validated_score"], ascending=False)
    return accepted.iloc[0].to_dict(), calibrated, accepted


_, ridge_candidate_table_raw = held_out_ridge_search(logpolar)
ridge_null_draws, ridge_null_summary = build_m_null_calibration(logpolar)
ridge_geometry, ridge_candidate_table, accepted_real = apply_m_null_calibration(
    ridge_candidate_table_raw, ridge_null_summary)
display(accepted_real.groupby(["m_arms", "winding_sign"], as_index=False)
        .first().sort_values("null_z", ascending=False))
display(ridge_null_summary)
print("RIDGE_NULL_NORMALIZATION_PASS")

# The fiducial 0.25-kpc result is reported. Every sensitivity width receives
# the same per-m null calibration; raw scores across widths are not comparable.
ridge_width_rows = []
ridge_width_null_draws = []
for width_kpc in RIDGE_WIDTH_SENSITIVITY_KPC:
    if np.isclose(width_kpc, RIDGE_CORE_WIDTH_KPC):
        width_geometry = ridge_geometry
        width_draws = ridge_null_draws.copy()
    else:
        _, width_table_raw = held_out_ridge_search(
            logpolar, core_width_kpc=float(width_kpc))
        width_draws, width_summary = build_m_null_calibration(
            logpolar, core_width_kpc=float(width_kpc))
        width_geometry, _, _ = apply_m_null_calibration(
            width_table_raw, width_summary)
    ridge_width_null_draws.append(width_draws)
    ridge_width_rows.append({
        "core_width_kpc": float(width_kpc),
        "m_arms": int(width_geometry["m_arms"]),
        "pitch_angle": float(width_geometry["pitch_angle"]),
        "Theta0": float(width_geometry["Theta0"]),
        "validated_score": float(width_geometry["validated_score"]),
        "null_z": float(width_geometry["null_z"]),
    })
ridge_width_sensitivity = pd.DataFrame(ridge_width_rows)
ridge_width_null_draws = pd.concat(
    ridge_width_null_draws, ignore_index=True)
assert len(ridge_width_null_draws) == (
    len(RIDGE_WIDTH_SENSITIVITY_KPC) * RIDGE_N_NULL * len(M_CANDIDATES))
display(ridge_width_sensitivity)
print("RIDGE_WIDTH_SENSITIVITY_COMPLETE")

if ridge_geometry["pitch_angle"] >= 0:
    raise RuntimeError(
        "The sign-neutral local-ridge/null-calibrated search did not recover "
        "the visually supported negative winding; report no acceptable global "
        "logarithmic geometry instead of relabelling a crossing curve as an arm."
    )
if ridge_geometry["validated_score"] <= 0 or ridge_geometry["null_z"] <= 0:
    raise RuntimeError("Best ridge geometry is not above its per-m null baseline")
print("RIDGE_REAL_FIT_COMPLETE")
print(ridge_geometry)
```

Do not assert a particular `m_arms`; the expanded, null-calibrated search must
determine it. Do assert the sign and positive validated/null-calibrated score
because they are real-data acceptance conditions from the approved design.
The independently prototyped live expectation is a fiducial `m=3`, pitch near
`-33 deg`; treat that value as a regression clue, not a hard-coded selector.

- [ ] **Step 5: Replace the legacy full-model bootstrap with ridge-sector stability**

Replace the entire `## 12. Spatial-sector bootstrap uncertainty` Markdown and code cell (currently id `a55074f1`). Explain that the previous bootstrap refitted the now-deleted full-field Fourier selector, while the revised diagnostic leaves out one azimuth sector at a time and reruns the ridge selector. Use:

```python
def omit_logpolar_sector(logpolar_map, omitted_sector,
                         n_sectors=RIDGE_N_SECTORS):
    phi = logpolar_map["phi"]
    sector_index = np.floor(
        (phi + np.pi) / (2.0*np.pi) * n_sectors).astype(int)
    sector_index = np.clip(sector_index, 0, n_sectors - 1)
    result = {
        key: (value.copy() if isinstance(value, np.ndarray) else value)
        for key, value in logpolar_map.items()
    }
    result["valid"][:, sector_index == omitted_sector] = False
    result["coverage"][~result["valid"]] = 0.0
    result["radial_residual"][~result["valid"]] = np.nan
    result["local_ridge"][~result["valid"]] = np.nan
    return result


sector_rows = []
for omitted_sector in range(RIDGE_N_SECTORS):
    try:
        _, sector_table = held_out_ridge_search(
            omit_logpolar_sector(logpolar, omitted_sector))
        sector_table = sector_table.merge(
            ridge_null_summary, on="m_arms", how="left", validate="many_to_one")
        sector_table["null_z"] = (
            (sector_table["validated_score"] - sector_table["null_mean"])
            / sector_table["null_std_floor"])
        sector_candidates = sector_table.loc[
            (sector_table["valid_held_out"] >= RIDGE_MIN_HELD_OUT_SECTORS)
            & (sector_table["validated_score"] > 0)
            & np.isfinite(sector_table["null_z"])].sort_values(
                ["null_z", "validated_score"], ascending=False)
        if sector_candidates.empty:
            raise RuntimeError("No null-calibrated sector-omission candidate")
        sector_geometry = sector_candidates.iloc[0].to_dict()
        sector_rows.append({
            "omitted_sector": omitted_sector,
            "success": True,
            "error": "",
            "m_arms": sector_geometry["m_arms"],
            "pitch_angle": sector_geometry["pitch_angle"],
            "Theta0": sector_geometry["Theta0"],
            "held_out_score": sector_geometry["held_out_score"],
            "null_z": sector_geometry["null_z"],
        })
    except Exception as exc:
        sector_rows.append({
            "omitted_sector": omitted_sector,
            "success": False,
            "error": f"{type(exc).__name__}: {exc}",
        })

ridge_sector_results = pd.DataFrame(sector_rows)
valid_sector_results = ridge_sector_results.loc[
    ridge_sector_results["success"]].copy()
sector_valid_fraction = len(valid_sector_results) / RIDGE_N_SECTORS
if sector_valid_fraction < 0.8:
    print("MODEL_STABILITY_WARNING: fewer than 80% of sector omissions succeeded")
ridge_m_frequency = (valid_sector_results["m_arms"].value_counts(normalize=True)
                     .sort_index())
ridge_negative_fraction = float(
    (valid_sector_results["pitch_angle"] < 0).mean())
same_m = valid_sector_results[
    (valid_sector_results["m_arms"] == ridge_geometry["m_arms"])
    & (np.sign(valid_sector_results["pitch_angle"])
       == np.sign(ridge_geometry["pitch_angle"]))].copy()
if len(same_m):
    pivot_u = np.log(RIDGE_PIVOT_RADIUS_KPC / R_REF_KPC)
    same_m["pivot_phase"] = (
        same_m["Theta0"]
        + same_m["m_arms"]
        / np.tan(np.deg2rad(same_m["pitch_angle"])) * pivot_u)
    ridge_phase_stability = float(np.abs(np.mean(
        np.exp(1j * same_m["pivot_phase"]))))
else:
    ridge_phase_stability = np.nan
display(ridge_sector_results)
display(ridge_m_frequency.rename("selection_fraction").to_frame())
print(f"Sector omission valid fraction={sector_valid_fraction:.1%}; "
      f"negative-winding fraction={ridge_negative_fraction:.1%}; "
      f"same-m-and-winding pivot-phase concentration={ridge_phase_stability:.3f} "
      f"at R={RIDGE_PIVOT_RADIUS_KPC:.3f} kpc")
if (ridge_m_frequency.max() < 0.8
        or 0.1 < ridge_negative_fraction < 0.9
        or (np.isfinite(ridge_phase_stability) and ridge_phase_stability < 0.8)):
    print("MODEL_STABILITY_WARNING: m, winding, or phase is sector-sensitive; "
          "retain this caution in the final interpretation.")

fig, ax = plt.subplots(figsize=(8.8, 4.8), constrained_layout=True)
scatter = ax.scatter(
    valid_sector_results["omitted_sector"],
    valid_sector_results["pitch_angle"],
    c=valid_sector_results["m_arms"], cmap="viridis",
    vmin=0.5, vmax=6.5, s=55)
ax.axhline(ridge_geometry["pitch_angle"], color="black", ls="--",
           label="all-sector selection")
ax.set(xlabel="omitted azimuth-sector index", ylabel="selected pitch (deg)",
       title="Leave-one-sector-out ridge-geometry stability")
ax.legend(fontsize=8)
fig.colorbar(scatter, ax=ax, ticks=np.arange(1, 7), label="selected m")
plt.show()
print("RIDGE_SECTOR_STABILITY_COMPLETE")
```

This cell must contain no reference to `search_spiral_geometry`, `compare_candidate_modes`, `fit_one_table`, or `run_sector_bootstrap`.

- [ ] **Step 6: Verify contracts, synthetic recovery, and commit**

Run the two ridge-search source contracts, static parse, and disposable
execution through the synthetic block. Then execute through the real selector
against the live FITS file and require: all synthetic markers, all six injected
cases correct, eight null draws for each `m` at each of the five declared ridge
widths, a finite negative-pitch real winner, and no runtime warnings. The fresh
prototype expectation at the fiducial width is `m=3`,
pitch about `-33 deg`, while the width table must make its `m=3`/`m=4` and sign
sensitivity visible rather than converting it into a false certainty. Then:

```bash
git -C /Users/Igniz/Desktop/ICRAR add -- \
  further/20260713_NGC4254_KTZ_spiral_fit.ipynb \
  further/tests/test_20260713_ngc4254_ktz_spiral_fit.py
git -C /Users/Igniz/Desktop/ICRAR diff --cached --check
git -C /Users/Igniz/Desktop/ICRAR commit -m "feat: select spiral geometry from SFR ridges"
```

### Task 4: Fit the KTZ harmonic profile at fixed ridge geometry

**Files:**
- Modify: `further/20260713_NGC4254_KTZ_spiral_fit.ipynb`
- Modify: `further/tests/test_20260713_ngc4254_ktz_spiral_fit.py`

- [ ] **Step 1: Add a failing fixed-geometry and maxima contract**

Add:

```python
def test_ktz_fit_cannot_replace_ridge_geometry(self):
    text = notebook_source()
    for required in [
        "fit_ktz_profile_fixed_geometry",
        '"geometry_source": "null_calibrated_local_ridge"',
        "enumerate_profile_maxima",
        "minimize_scalar",
        "enhanced_above_background",
    ]:
        self.assertIn(required, text)
```

Run it and expect failure because the current nonlinear fit still refines pitch and phase.

- [ ] **Step 2: Explain KTZ scope and implement the fixed-geometry profile fit**

Replace the `## 9. Fit the harmonic arm profile and refine the complete source-rate field` Markdown. Define “KTZ-compatible” here as the spiral-modulated source-rate component used by the working KTZ26 validation notebook, which extends the Krumholz--Ting stochastic-metallicity framework. State that this notebook is fitting only present-day SFR morphology: it does not fit diffusion, injection ages, clustering, pattern speed, or metallicity covariance. State explicitly that `KTZ_validation.ipynb` prescribes the harmonic source profile and does **not** contain the ridge detector introduced here.

Replace the entire code cell currently id `ab860a7f` with the function below. This removal must eliminate `fit_harmonic_profile`, `recover_spiral_parameters`, `refine_full_model`, and `compare_candidate_modes`. The new optimizer deliberately excludes `m_arms`, pitch, and `Theta0` from its parameter vector. Define the function here but do not call it yet; the synthetic ridge tests in Section 10 must execute before the real-data call in Section 11.

```python
def fit_ktz_profile_fixed_geometry(table, geometry, radial_guess):
    radius = table["radius_kpc"].to_numpy(float)
    azimuth = table["azimuth_rad"].to_numpy(float)
    observed_log = np.log(table["sfr_linear"].to_numpy(float))
    sqrt_weight = np.sqrt(table["area_pix"].to_numpy(float))
    sqrt_weight /= np.nanmedian(sqrt_weight)
    theta = spiral_phase(radius, azimuth, int(geometry["m_arms"]),
                         float(geometry["pitch_angle"]), float(geometry["Theta0"]))
    initial = np.array([
        np.log(radial_guess["lambda0_0"]), np.log(radial_guess["h_R"]),
        np.log(0.2), 0.3, 0.15, 0.0, 0.0,
    ])
    lower = np.array([-50.0, np.log(0.2), np.log(1e-4),
                      0.0, 0.0, -np.pi, -np.pi])
    upper = np.array([50.0, np.log(30.0), np.log(2.0),
                      1.5, 1.5, np.pi, np.pi])
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
        profile = arm_profile(theta, parameters["harmonic_g"],
                              parameters["harmonic_alpha"])
        factor = 1.0 + parameters["eta"] * profile
        dense_factor = 1.0 + parameters["eta"] * arm_profile(
            dense_theta, parameters["harmonic_g"], parameters["harmonic_alpha"])
        model_log = (np.log(parameters["lambda0_0"])
                     - radius / parameters["h_R"]
                     + np.log(np.clip(factor, 1e-10, None)))
        positivity_penalty = min(0.0, float(np.min(dense_factor)) - 1e-3) * 1e5
        return np.concatenate([sqrt_weight * (observed_log - model_log),
                               [positivity_penalty]])

    result = least_squares(residual, initial, bounds=(lower, upper),
                           loss="soft_l1", f_scale=0.25, max_nfev=3000)
    parameters = unpack(result.x)
    final_residual = residual(result.x)
    minimum_factor = float(np.min(1.0 + parameters["eta"] * arm_profile(
        dense_theta, parameters["harmonic_g"], parameters["harmonic_alpha"])))
    parameters.update({
        "m_arms": int(geometry["m_arms"]),
        "pitch_angle": float(geometry["pitch_angle"]),
        "Theta0": float(geometry["Theta0"]),
        "geometry_source": "null_calibrated_local_ridge",
        "weighted_log_sse": float(np.sum(final_residual[:-1]**2)),
        "positivity_penalty_residual": float(final_residual[-1]),
        "robust_cost": float(2.0 * result.cost),
        "minimum_rate_factor": minimum_factor,
        "success": bool(result.success and minimum_factor > 0.0),
        "message": result.message,
    })
    if not parameters["success"]:
        raise RuntimeError(f"Fixed-geometry KTZ fit failed: {parameters}")
    return parameters
```

- [ ] **Step 3: Enumerate every harmonic maximum**

Add:

```python
def enumerate_profile_maxima(model, n_grid=65536):
    theta = np.linspace(0.0, 2.0*np.pi, n_grid, endpoint=False)
    profile = arm_profile(theta, model["harmonic_g"], model["harmonic_alpha"])
    maxima = (profile > np.roll(profile, 1)) & (profile >= np.roll(profile, -1))
    indices = np.flatnonzero(maxima)
    if len(indices) == 0:
        raise RuntimeError("No harmonic-profile maximum found")
    grid_step = 2.0 * np.pi / n_grid

    def profile_scalar(value):
        wrapped_value = float(value % (2.0*np.pi))
        return float(arm_profile(
            np.array([wrapped_value]),
            model["harmonic_g"], model["harmonic_alpha"])[0])

    rows = []
    for index in indices:
        centre = float(theta[index])
        refined = minimize_scalar(
            lambda offset: -profile_scalar(centre + offset),
            bounds=(-grid_step, grid_step), method="bounded",
            options={"xatol": 1e-12})
        theta_peak = float((centre + refined.x) % (2.0*np.pi))
        h_peak = profile_scalar(theta_peak)
        factor = float(1.0 + model["eta"] * h_peak)
        rows.append({
            "theta_peak": theta_peak,
            "h_peak": h_peak,
            "rate_factor": factor,
            "enhanced_above_background": bool(factor > 1.0),
        })
    return pd.DataFrame(rows).sort_values("rate_factor", ascending=False).reset_index(drop=True)


```

Append `enumerate_profile_maxima` to the same Section 9 code cell after `fit_ktz_profile_fixed_geometry`. Then append the following invocation block to the Section 11 real-data cell, after Task 3's `RIDGE_REAL_FIT_COMPLETE` print. This order makes the geometry freeze auditable:

```python
real_fit = fit_ktz_profile_fixed_geometry(fit_table, ridge_geometry, radial_fit)
real_fit["gradient_dex_per_kpc"] = float(
    -1.0 / (np.log(10.0) * real_fit["h_R"]))
theta_fit = spiral_phase(
    fit_table["radius_kpc"].to_numpy(float),
    fit_table["azimuth_rad"].to_numpy(float),
    real_fit["m_arms"], real_fit["pitch_angle"], real_fit["Theta0"],
)
rate_factor = 1.0 + real_fit["eta"] * arm_profile(
    theta_fit, real_fit["harmonic_g"], real_fit["harmonic_alpha"])
model_sfr = (real_fit["lambda0_0"]
             * np.exp(-fit_table["radius_kpc"].to_numpy(float) / real_fit["h_R"])
             * rate_factor)
fit_table = fit_table.copy()
fit_table["theta_fit"] = np.mod(theta_fit, 2.0*np.pi)
fit_table["rate_factor"] = rate_factor
fit_table["model_sfr"] = model_sfr
fit_table["residual_dex"] = np.log10(
    fit_table["sfr_linear"].to_numpy(float)
    / np.clip(model_sfr, 1e-30, None))
profile_maxima = enumerate_profile_maxima(real_fit)
display(profile_maxima)
assert profile_maxima["enhanced_above_background"].any()
print("KTZ_PROFILE_FIT_COMPLETE")
print("HARMONIC_MAXIMA_COMPLETE")
print("POSITIVITY_CHECK_PASS")
print(real_fit)
```

Task 5 will replace the final parameter table so `m_arms`, pitch, and `Theta0` are labelled `ridge-selected`, while `h_R`, `lambda0_0`, `eta`, `g`, and `alpha` are labelled `fixed-geometry KTZ fit`.

- [ ] **Step 4: Verify GREEN and commit**

Run the fixed-geometry contract and static parse. In a disposable partial execution, verify that `RIDGE_REAL_FIT_COMPLETE` occurs before `KTZ_PROFILE_FIT_COMPLETE` and that the two geometry dictionaries have identical `m_arms`, pitch, and `Theta0`.

```bash
git -C /Users/Igniz/Desktop/ICRAR add -- \
  further/20260713_NGC4254_KTZ_spiral_fit.ipynb \
  further/tests/test_20260713_ngc4254_ktz_spiral_fit.py
git -C /Users/Igniz/Desktop/ICRAR diff --cached --check
git -C /Users/Igniz/Desktop/ICRAR commit -m "feat: fit KTZ profile at ridge-selected geometry"
```

### Task 5: Build the approved 2-by-2 diagnostic and expose all ridges

**Files:**
- Modify: `further/20260713_NGC4254_KTZ_spiral_fit.ipynb`
- Modify: `further/tests/test_20260713_ngc4254_ktz_spiral_fit.py`

- [ ] **Step 1: Run the layout contract to verify RED**

Run `test_all_harmonic_maxima_and_approved_layout_are_declared`. Expected: failure until the existing 1-by-3 plot is replaced.

- [ ] **Step 2: Replace the legacy overlay cell and verify projection handedness**

Replace the `## 13. What the “KTZ model” is, and how its fitted ridge is drawn` Markdown so it distinguishes the `Theta=0` data ridge from maxima of the fitted harmonic profile, reiterates that FITS WCS handles sky orientation without a manual flip, and removes every claim about an opposite-winding overlay. Replace the entire following code cell (currently id `1e4ee86a`) with the helpers and checks below. This deletes `profile_peak_phase`, `skeleton_xy`, `opposite_curves`, `opposite_winding_fit`, and the old 1-by-3 figure.

```python
def skeleton_xy_for_phase(geometry, radius_grid, theta_target=0.0):
    curves = []
    m_arms = int(geometry["m_arms"])
    pitch = np.deg2rad(float(geometry["pitch_angle"]))
    theta0 = float(geometry["Theta0"])
    for branch in range(m_arms):
        phi = ((m_arms / np.tan(pitch)) * np.log(radius_grid / R_REF_KPC)
               + theta0 - theta_target - 2.0*np.pi*branch) / m_arms
        curves.append((radius_grid * np.cos(phi), radius_grid * np.sin(phi)))
    return curves


def project_disc_to_pixel(x_disc, y_disc):
    """Invert the fixed deprojection and return image-pixel coordinates."""
    pa_rad = np.deg2rad(POSITION_ANGLE_DEG)
    inc_rad = np.deg2rad(INCLINATION_DEG)
    major = np.asarray(x_disc)
    minor_projected = np.asarray(y_disc) * np.cos(inc_rad)
    east_kpc = major * np.sin(pa_rad) - minor_projected * np.cos(pa_rad)
    north_kpc = major * np.cos(pa_rad) + minor_projected * np.sin(pa_rad)
    sky = CENTER.spherical_offsets_by(
        east_kpc / ARCSEC_TO_KPC * u.arcsec,
        north_kpc / ARCSEC_TO_KPC * u.arcsec,
    )
    return celestial_wcs.world_to_pixel(sky)


# Real-coordinate round trip.
roundtrip_sample = bin_table.iloc[
    np.linspace(0, len(bin_table) - 1, 64, dtype=int)]
roundtrip_x, roundtrip_y = project_disc_to_pixel(
    roundtrip_sample["x_disc_kpc"], roundtrip_sample["y_disc_kpc"])
roundtrip_error_pix = np.hypot(
    roundtrip_x - roundtrip_sample["x_pix"],
    roundtrip_y - roundtrip_sample["y_pix"])
assert np.nanmax(roundtrip_error_pix) < 0.1, np.nanmax(roundtrip_error_pix)
print(f"Maximum sky-disc-pixel round-trip error: "
      f"{np.nanmax(roundtrip_error_pix):.4f} pixel")
print("DEPROJECTION_ROUNDTRIP_PASS")

# Projection/deprojection must preserve the sign of d ln(R) / d phi = tan(pitch).
orientation_radius = np.geomspace(max(0.5, R_IN_KPC), min(8.0, R_OUT_KPC), 600)
for orientation_pitch in (-25.0, 25.0):
    orientation_geometry = {
        "m_arms": 2, "pitch_angle": orientation_pitch, "Theta0": 0.4}
    orientation_curve = skeleton_xy_for_phase(
        orientation_geometry, orientation_radius, theta_target=0.0)[0]
    xp, yp = project_disc_to_pixel(*orientation_curve)
    ra_deg, dec_deg = celestial_wcs.pixel_to_world_values(xp, yp)
    _, _, radius_back, phi_back = deproject_sky(ra_deg, dec_deg)
    recovered_slope = float(np.polyfit(
        np.unwrap(phi_back), np.log(radius_back / R_REF_KPC), 1)[0])
    expected_slope = float(np.tan(np.deg2rad(orientation_pitch)))
    assert np.sign(recovered_slope) == np.sign(expected_slope)
    assert np.isclose(recovered_slope, expected_slope, rtol=2e-3, atol=2e-4)
print("SYNTHETIC_WINDING_SIGN_PASS")


radius_grid = np.geomspace(RIDGE_R_IN_KPC, RIDGE_R_OUT_KPC, 1000)
data_ridge_curves = skeleton_xy_for_phase(ridge_geometry, radius_grid, theta_target=0.0)
enhanced_model_curves = []
for row in profile_maxima.itertuples(index=False):
    if row.enhanced_above_background:
        enhanced_model_curves.append({
            "theta_peak": row.theta_peak,
            "rate_factor": row.rate_factor,
            "curves": skeleton_xy_for_phase(real_fit, radius_grid, row.theta_peak),
        })
```

The observed panels must use only `data_ridge_curves`. The model panel must show every entry in `enhanced_model_curves`; line width or alpha may scale monotonically with `rate_factor`.

- [ ] **Step 3: Replace the main diagnostic with Layout A**

Append this complete layout structure to the same Section 13 code cell:

```python
fig, axes = plt.subplots(2, 2, figsize=(13.5, 12.0), constrained_layout=True)
ax_sky, ax_face, ax_model, ax_residual = axes.ravel()
disc_limit = 1.03 * RIDGE_R_OUT_KPC

sky_image = ax_sky.imshow(log_sfr_map, origin="lower", cmap="magma",
                          vmin=vmin, vmax=vmax)
for index, (x_curve, y_curve) in enumerate(data_ridge_curves):
    xp, yp = project_disc_to_pixel(x_curve, y_curve)
    ax_sky.plot(xp, yp, color="cyan", lw=1.8,
                label="ridge-selected skeleton" if index == 0 else "_nolegend_")
ax_sky.set(xlim=(0, log_sfr_map.shape[1]-1),
           ylim=(0, log_sfr_map.shape[0]-1),
           xlabel="image x (pixel)", ylabel="image y (pixel)",
           title="observed sky-plane log SFR + data-derived skeleton")
ax_sky.set_aspect("equal")
ax_sky.legend(fontsize=8)
fig.colorbar(sky_image, ax=ax_sky, shrink=0.82, label=r"log $\Sigma_{\rm SFR}$")

face_scatter = ax_face.scatter(hii_pixels["x_disc_kpc"], hii_pixels["y_disc_kpc"],
                               c=hii_pixels["log_sfr"], s=0.35, cmap="magma",
                               vmin=vmin, vmax=vmax, linewidth=0, rasterized=True)
for x_curve, y_curve in data_ridge_curves:
    ax_face.plot(x_curve, y_curve, color="cyan", lw=1.8)
ax_face.set(xlabel="disc major-axis x (kpc)", ylabel="disc minor-axis y (kpc)",
            title="observed deprojected log SFR + data-derived skeleton",
            xlim=(-disc_limit, disc_limit), ylim=(-disc_limit, disc_limit))
ax_face.set_aspect("equal")
fig.colorbar(face_scatter, ax=ax_face, shrink=0.82, label=r"log $\Sigma_{\rm SFR}$")

model_scatter = ax_model.scatter(fit_table["x_disc_kpc"], fit_table["y_disc_kpc"],
                                 c=np.log10(fit_table["model_sfr"]), s=2.0,
                                 cmap="magma", vmin=vmin, vmax=vmax, linewidth=0)
for peak_index, peak in enumerate(enhanced_model_curves):
    width = 0.8 + 1.5 * (peak["rate_factor"] - 1.0)
    colour = f"C{peak_index % 10}"
    for branch_index, (x_curve, y_curve) in enumerate(peak["curves"]):
        label = (fr"model maximum {peak_index+1}: factor={peak['rate_factor']:.2f}"
                 if branch_index == 0 else "_nolegend_")
        ax_model.plot(x_curve, y_curve, color=colour, lw=width, label=label)
ax_model.set(xlabel="disc major-axis x (kpc)", ylabel="disc minor-axis y (kpc)",
             title="fitted KTZ-compatible source model + all enhanced maxima",
             xlim=(-disc_limit, disc_limit), ylim=(-disc_limit, disc_limit))
ax_model.set_aspect("equal")
ax_model.legend(fontsize=7)
fig.colorbar(model_scatter, ax=ax_model, shrink=0.82,
             label=r"model log $\Sigma_{\rm SFR}$")

limit = np.nanpercentile(np.abs(fit_table["residual_dex"]), 95)
residual_scatter = ax_residual.scatter(
    fit_table["x_disc_kpc"], fit_table["y_disc_kpc"],
    c=fit_table["residual_dex"], s=2.0, cmap="coolwarm",
    vmin=-limit, vmax=limit, linewidth=0)
for x_curve, y_curve in data_ridge_curves:
    ax_residual.plot(x_curve, y_curve, color="cyan", lw=1.8, alpha=0.9)
ax_residual.set(xlabel="disc major-axis x (kpc)", ylabel="disc minor-axis y (kpc)",
                title="observed - model residual + data-derived skeleton",
                xlim=(-disc_limit, disc_limit), ylim=(-disc_limit, disc_limit))
ax_residual.set_aspect("equal")
fig.colorbar(residual_scatter, ax=ax_residual, shrink=0.82, label="residual (dex)")

plt.show()
print("DEPROJECTED_SFR_SKELETON_PASS")
```

- [ ] **Step 4: Update the arm-profile plot**

Replace the `## 14. KTZ-style Spiral skeleton and Arm profile` Markdown and the entire following code cell (currently id `c7c05d44`). The Markdown must state that these are source-model maxima and are not the data-derived geometry selector. The code below removes stale references to `curves`, `opposite_curves`, and `theta_ridge`, retains the harmonic components and phase-folded measurements, and marks every refined maximum. Solid coloured lines denote enhanced maxima; thin grey dotted lines denote below-background maxima.

```python
theta_grid = np.linspace(0.0, 2.0*np.pi, 1600)
profile_total = arm_profile(
    theta_grid, real_fit["harmonic_g"], real_fit["harmonic_alpha"])
rate_factor_grid = 1.0 + real_fit["eta"] * profile_total
if np.min(rate_factor_grid) <= 0:
    raise ValueError("Fitted source-rate profile becomes non-positive")

final_background = (real_fit["lambda0_0"]
    * np.exp(-fit_table["radius_kpc"].to_numpy(float) / real_fit["h_R"]))
fit_table["q_final_background"] = (
    fit_table["sfr_linear"].to_numpy(float) / final_background - 1.0)
phase_edges = np.linspace(0.0, 2.0*np.pi, 25)
phase_index = np.digitize(fit_table["theta_fit"], phase_edges) - 1
phase_rows = []
for phase_bin in range(len(phase_edges) - 1):
    use = phase_index == phase_bin
    if use.sum() >= 10:
        values = fit_table.loc[use, "q_final_background"].to_numpy(float)
        median = float(np.median(values))
        robust_sigma = float(1.4826 * np.median(np.abs(values - median)))
        phase_rows.append({
            "theta": 0.5 * (phase_edges[phase_bin] + phase_edges[phase_bin + 1]),
            "q_over_eta": median / real_fit["eta"],
            "q_over_eta_sem": robust_sigma / np.sqrt(use.sum()) / real_fit["eta"],
        })
phase_summary = pd.DataFrame(phase_rows)

fig, ax = plt.subplots(figsize=(10.5, 5.8), constrained_layout=True)
for n, g, alpha in zip(
        HARMONIC_N, real_fit["harmonic_g"], real_fit["harmonic_alpha"]):
    ax.plot(theta_grid, g * np.cos(n*theta_grid + alpha), "--", lw=1,
            label=fr"component $n={n}$")
ax.plot(theta_grid, profile_total, "k-", lw=2.2, label=r"total $h(\Theta)$")
ax.errorbar(phase_summary["theta"], phase_summary["q_over_eta"],
            yerr=phase_summary["q_over_eta_sem"], fmt="o", ms=3,
            color="tab:red", alpha=0.75, label="phase-folded data / eta")
for peak_index, row in enumerate(profile_maxima.itertuples(index=False)):
    if row.enhanced_above_background:
        colour = f"C{(peak_index + 3) % 10}"
        linestyle = "-"
        width = 1.4
        label = f"enhanced maximum {peak_index + 1}"
    else:
        colour = "0.55"
        linestyle = ":"
        width = 0.9
        label = f"below-background maximum {peak_index + 1}"
    ax.axvline(row.theta_peak, color=colour, ls=linestyle, lw=width,
               label=label)
    ax.annotate(f"factor={row.rate_factor:.2f}",
                (row.theta_peak, row.h_peak), xytext=(4, 7),
                textcoords="offset points", rotation=90, fontsize=7,
                color=colour)
ax.axhline(0.0, color="0.5", lw=0.7)
ax.set(xlim=(0.0, 2.0*np.pi), xlabel=r"$\Theta$", ylabel=r"$h(\Theta)$",
       title="Fixed-geometry KTZ-compatible harmonic source profile")
ax.set_xticks([0.0, np.pi, 2.0*np.pi], [r"$0$", r"$\pi$", r"$2\pi$"])
ax.legend(fontsize=7, ncol=2)
plt.show()
print(f"Minimum fitted rate factor: {np.min(rate_factor_grid):.4f}")
print("POSITIVITY_CHECK_PASS")
```

- [ ] **Step 5: Add the `m=1..6` candidate diagnostic**

Append this separate candidate figure to the Section 11 real-data cell. It is the fair comparison; do not overlay an alternative with more branches on the observed panels.

```python
fig, (ax, ax_width) = plt.subplots(
    1, 2, figsize=(14.5, 5.5), constrained_layout=True)
for m_arms, group in ridge_candidate_table.groupby("m_arms"):
    group = group.loc[np.isfinite(group["null_z"])].sort_values("pitch_angle")
    ax.plot(group["pitch_angle"], group["null_z"],
            marker=".", ms=3, lw=1.0, label=f"m={m_arms}")
ax.scatter([ridge_geometry["pitch_angle"]], [ridge_geometry["null_z"]],
           marker="*", s=180, color="black", zorder=5, label="selected ridge")
ax.axvline(0.0, color="0.6", lw=0.8)
ax.set(xlabel="signed pitch angle (deg)",
       ylabel="m-specific row-scrambled null z-score",
       title="Local-ridge candidates: m=1..6 and both winding signs")
ax.legend(ncol=4, fontsize=8)

width_scatter = ax_width.scatter(
    ridge_width_sensitivity["core_width_kpc"],
    ridge_width_sensitivity["pitch_angle"],
    c=ridge_width_sensitivity["m_arms"], cmap="viridis",
    vmin=0.5, vmax=6.5, s=75)
ax_width.axhline(0.0, color="0.6", lw=0.8)
ax_width.axvline(RIDGE_CORE_WIDTH_KPC, color="black", ls="--", lw=1,
                 label="fiducial width")
ax_width.set(xlabel="physical ridge-core width (kpc)",
             ylabel="selected signed pitch (deg)",
             title="Declared ridge-width sensitivity")
ax_width.legend(fontsize=8)
fig.colorbar(width_scatter, ax=ax_width, ticks=M_CANDIDATES,
             label="selected m")
plt.show()
```

- [ ] **Step 6: Replace the stale bootstrap parameter table and interpretation**

Replace the `## 15. Final KTZ-compatible parameter table` Markdown and code cell (currently id `95fdc1a2`) completely. Explain that ridge-sector resampling quantifies geometry stability, while uncertainty on the fixed-geometry KTZ amplitude/profile parameters is not estimated by this notebook. Use:

```python
selected_m_frequency = float(
    ridge_m_frequency.get(ridge_geometry["m_arms"], 0.0))
parameter_rows = [
    ("R_in", R_IN_KPC, "kpc", "all positive-radius valid HII bins; no hole"),
    ("R_out", R_OUT_KPC, "kpc", "all valid HII bins; no quantile cut"),
    ("ridge_R_in", RIDGE_R_IN_KPC, "kpc", "all positive-radius valid HII pixels"),
    ("ridge_R_out", RIDGE_R_OUT_KPC, "kpc", "all valid HII morphology pixels"),
    ("m_arms", real_fit["m_arms"], "integer base mode",
     f"ridge-selected; sector frequency={selected_m_frequency:.1%}"),
    ("pitch_angle", real_fit["pitch_angle"], "signed deg",
     f"ridge-selected; negative sector fraction={ridge_negative_fraction:.1%}"),
    ("Theta0", real_fit["Theta0"], "rad modulo 2 pi",
     f"ridge-selected; same-m/winding pivot-phase concentration="
     f"{ridge_phase_stability:.3f} at {RIDGE_PIVOT_RADIUS_KPC:.3f} kpc"),
    ("ridge_core_width", RIDGE_CORE_WIDTH_KPC, "kpc",
     "declared fiducial local-ridge scale"),
    ("held_out_score", ridge_geometry["held_out_score"], "local response",
     "median screened sector-validation statistic"),
    ("validated_score", ridge_geometry["validated_score"], "local response",
     "held-out score times sector coverage and phase stability"),
    ("null_z", ridge_geometry["null_z"], "m-specific z-score",
     f"row-scrambled calibration; N_null={RIDGE_N_NULL}"),
    ("null_mean", ridge_geometry["null_mean"], "local response",
     "best-score null mean for selected m"),
    ("null_std", ridge_geometry["null_std"], "local response",
     "sample standard deviation of best-score null for selected m"),
    ("lambda0_0", real_fit["lambda0_0"], "SFR surface-density normalization",
     "fixed-geometry KTZ fit; uncertainty not estimated"),
    ("h_R", real_fit["h_R"], "kpc",
     "axisymmetric SFR e-folding length; fixed-geometry KTZ fit"),
    ("radial_gradient", real_fit["gradient_dex_per_kpc"], "dex/kpc",
     "equals -1 / (ln(10) h_R)"),
    ("eta", real_fit["eta"], "dimensionless",
     "fixed-geometry KTZ fit; uncertainty not estimated"),
    ("harmonic_n", str(real_fit["harmonic_n"].tolist()), "integer orders", "fixed"),
    ("harmonic_g", np.array2string(real_fit["harmonic_g"], precision=4),
     "g1=1 normalization", "fixed-geometry KTZ fit"),
    ("harmonic_alpha", np.array2string(real_fit["harmonic_alpha"], precision=4),
     "rad; alpha1=0", "fixed-geometry KTZ fit"),
    ("kappa, x0, t_star, clustering", np.nan, "KTZ mixing parameters",
     "not fitted from SFR morphology"),
]
parameter_table = pd.DataFrame(
    parameter_rows, columns=["parameter", "estimate", "unit_or_convention", "status"])
display(parameter_table)
display(profile_maxima)
display(ridge_m_frequency.rename("sector_selection_fraction").to_frame())
print("FINAL_PARAMETER_TABLE_COMPLETE")
```

Replace the `## 16. Interpretation limits and the next KTZ step` Markdown so
it describes local-flank ridge selection, per-`m` row-scrambled null
calibration, the declared physical-width scan, and leave-one-sector-out
stability, not the deleted global Fourier/full-model bootstrap. State plainly
that the fiducial live prototype selects `m=3`, pitch near `-33 deg`, while
nearby negative-pitch `m=4` solutions and one positive-pitch width-scale result
show that a single global logarithmic template is not unique. Retain the
cautions about NGC4254's lopsided/bifurcating morphology and about this not
being a metallicity-diffusion fit.

- [ ] **Step 7: Verify layout source contract and commit**

Run the layout contract and static parse. Expected: pass.

```bash
git -C /Users/Igniz/Desktop/ICRAR add -- \
  further/20260713_NGC4254_KTZ_spiral_fit.ipynb \
  further/tests/test_20260713_ngc4254_ktz_spiral_fit.py
git -C /Users/Igniz/Desktop/ICRAR diff --cached --check
git -C /Users/Igniz/Desktop/ICRAR commit -m "feat: show sky and deprojected spiral diagnostics"
```

### Task 6: Execute, inspect, and accept the revised notebook

**Files:**
- Modify: `further/20260713_NGC4254_KTZ_spiral_fit.ipynb`
- Modify: `further/docs/superpowers/plans/2026-07-15-ngc4254-ridge-spiral-revision-plan.md`

- [ ] **Step 1: Run the full test suite before execution and verify the expected RED state**

Run:

```bash
/opt/miniconda3/envs/ICRAR/bin/python -m unittest \
  tests.test_20260713_ngc4254_ktz_spiral_fit -v
```

Expected: source contracts pass; only the saved-result contract fails because edited cells have no execution counts or saved markers.

- [ ] **Step 2: Execute the notebook against the live FITS product**

Use writable temporary Jupyter directories and preserve a failed attempt for diagnosis:

```bash
MPLCONFIGDIR=/private/tmp/mpl-ngc4254-ridge \
IPYTHONDIR=/private/tmp/ipython-ngc4254-ridge \
JUPYTER_CONFIG_DIR=/private/tmp/jupyter-config-ngc4254-ridge \
JUPYTER_RUNTIME_DIR=/private/tmp/jupyter-runtime-ngc4254-ridge \
/opt/miniconda3/envs/ICRAR/bin/python - <<'PY'
from pathlib import Path
import nbformat
from nbclient import NotebookClient

path = Path("20260713_NGC4254_KTZ_spiral_fit.ipynb")
notebook = nbformat.read(path, as_version=4)
client = NotebookClient(notebook, timeout=3600, kernel_name="python3",
                        resources={"metadata": {"path": str(path.parent)}})
try:
    client.execute()
finally:
    nbformat.write(notebook, "/private/tmp/20260713_NGC4254_ridge_attempt.ipynb")
nbformat.write(notebook, path)
print("EXECUTION_COMPLETE")
PY
```

If local kernel binding is sandbox-blocked, rerun this exact command with the required approval. Expected: `EXECUTION_COMPLETE`, no traceback, and all numerical markers.

- [ ] **Step 3: Run the complete fresh verification gate**

Run:

```bash
/opt/miniconda3/envs/ICRAR/bin/python -m unittest \
  tests.test_20260713_ngc4254_ktz_spiral_fit -v
/opt/miniconda3/envs/ICRAR/bin/python - <<'PY'
import ast
import nbformat

nb = nbformat.read("20260713_NGC4254_KTZ_spiral_fit.ipynb", as_version=4)
errors = []
unexecuted = []
pngs = 0
streams = []
for index, cell in enumerate(nb.cells):
    if cell.cell_type != "code":
        continue
    ast.parse(cell.source)
    if cell.execution_count is None:
        unexecuted.append(index)
    for output in cell.get("outputs", []):
        if output.output_type == "error":
            errors.append((index, output.ename, output.evalue))
        if output.output_type == "stream":
            streams.append(output.text)
        if output.output_type in ("display_data", "execute_result"):
            pngs += int("image/png" in output.get("data", {}))
joined = "\n".join(streams)
assert not errors, errors
assert not unexecuted, unexecuted
assert "RuntimeWarning" not in joined
assert pngs >= 7, pngs
assert "RIDGE_REAL_FIT_COMPLETE" in joined
assert "KTZ_PROFILE_FIT_COMPLETE" in joined
assert "DEPROJECTED_SFR_SKELETON_PASS" in joined
print(f"NOTEBOOK_VERIFY_PASS cells={len(nb.cells)} pngs={pngs}")
PY
```

Expected: every unittest passes; notebook verification prints `NOTEBOOK_VERIFY_PASS`; zero errors, zero unexecuted code cells, and no runtime warnings.

- [ ] **Step 4: Extract and inspect every saved figure**

Mechanically extract PNG outputs to `/private/tmp/ngc4254_ridge_figures/`:

```bash
/opt/miniconda3/envs/ICRAR/bin/python - <<'PY'
import base64
from pathlib import Path
import nbformat

notebook = nbformat.read("20260713_NGC4254_KTZ_spiral_fit.ipynb", as_version=4)
output_dir = Path("/private/tmp/ngc4254_ridge_figures")
output_dir.mkdir(parents=True, exist_ok=True)
paths = []
for cell_index, cell in enumerate(notebook.cells):
    if cell.cell_type != "code":
        continue
    for output_index, output in enumerate(cell.get("outputs", [])):
        data = output.get("data", {})
        if "image/png" not in data:
            continue
        path = output_dir / f"cell_{cell_index:02d}_output_{output_index:02d}.png"
        path.write_bytes(base64.b64decode(data["image/png"]))
        paths.append(path)
print("\n".join(str(path) for path in paths))
assert len(paths) >= 7, len(paths)
PY
```

Open every printed path with the local image viewer at original detail. At minimum verify:

- the optical centre and PA axes are sensible in the sky panel;
- the observed sky and face-on panels show identical physical skeletons;
- the selected skeleton winds in the visually supported negative direction;
- the curves pass through the bright NGC4254 SFR tracks rather than merely covering more area;
- the model panel shows every enhanced harmonic maximum listed in the table;
- the residual panel does not reuse stale geometry; and
- all axes, legends, colour bars, and titles are readable without clipping.

If the selected ridge is visually wrong, stop at the ridge score and diagnose its contrast construction, masking, phase sign, or branch normalization. Do not force completion by changing the plot sign alone.

- [ ] **Step 5: Record final parameters and stability honestly**

Append a compact execution record to this plan containing:

- valid HII pixel and independent-bin counts;
- `b/a` validation value;
- selected `m_arms`, signed pitch, `Theta0`, held-out/validated ridge scores,
  and the per-`m` null z-score;
- best opposite-sign candidate, null-z difference, and the complete physical
  ridge-width sensitivity table;
- fitted `h_R`, its dex-per-kpc gradient, `eta`, `g`, and `alpha`;
- every harmonic maximum and whether it is enhanced;
- sector stability; and
- remaining limitations or rejection of a single global spiral.

- [ ] **Step 6: Commit only the executed notebook and plan record**

Run `git status --short`, confirm that unrelated notebooks, notes, theory files, and `.superpowers/` remain unstaged, then:

```bash
git -C /Users/Igniz/Desktop/ICRAR add -- \
  further/20260713_NGC4254_KTZ_spiral_fit.ipynb \
  further/docs/superpowers/plans/2026-07-15-ngc4254-ridge-spiral-revision-plan.md
git -C /Users/Igniz/Desktop/ICRAR diff --cached --check
git -C /Users/Igniz/Desktop/ICRAR diff --cached --stat
git -C /Users/Igniz/Desktop/ICRAR commit -m "chore: execute ridge-based NGC4254 spiral fit"
```

- [ ] **Step 7: Append the required job log**

Append a concise entry to `/Users/Igniz/Desktop/Codex_log/2026_07_15.md` with the request, files inspected/changed, exact verification, final result, unresolved risks, commit ids, and delegated-worker metadata. Do not include secrets or the 3.1 GB FITS product in Git.
