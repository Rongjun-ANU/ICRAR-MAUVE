import ast
import json
from pathlib import Path
import unittest
import warnings

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter, gaussian_filter1d


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "20260713_NGC4254_KTZ_spiral_fit.ipynb"


def load_notebook():
    assert NOTEBOOK.exists(), f"Notebook has not been created: {NOTEBOOK}"
    return json.loads(NOTEBOOK.read_text(encoding="utf-8"))


def notebook_source():
    nb = load_notebook()
    return "\n".join("".join(cell["source"]) for cell in nb["cells"])


def notebook_code_source():
    nb = load_notebook()
    return "\n".join(
        "".join(cell["source"])
        for cell in nb["cells"]
        if cell["cell_type"] == "code"
    )


def notebook_function_nodes():
    code = notebook_code_source()
    tree = ast.parse(code)
    nodes = [
        node for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    ]
    return code, nodes


def notebook_function_names():
    _, nodes = notebook_function_nodes()
    return {node.name for node in nodes}


def function_source(function_name):
    code, nodes = notebook_function_nodes()
    matches = [node for node in nodes if node.name == function_name]
    assert len(matches) == 1, (
        f"Expected one notebook FunctionDef named {function_name!r}; "
        f"found {len(matches)}"
    )
    source = ast.get_source_segment(code, matches[0])
    assert source is not None, f"Could not recover source for {function_name!r}"
    return source


def extracted_function_namespace(function_names, extra_namespace=None):
    _, nodes = notebook_function_nodes()
    requested = set(function_names)
    definitions = [
        node for node in nodes if node.name in requested
    ]
    found = {node.name for node in definitions}
    missing = sorted(requested - found)
    assert not missing, f"Notebook functions not found: {missing}"
    namespace = {
        "np": np,
        "gaussian_filter": gaussian_filter,
        "gaussian_filter1d": gaussian_filter1d,
        "LOGPOLAR_SMOOTH_SIGMA": (0.8, 1.2),
        "LOGPOLAR_AZIMUTH_BROAD_SIGMA_BINS": 30.0,
        "LOGPOLAR_N_U": 4,
        "LOGPOLAR_N_PHI": 8,
        "LOGPOLAR_N_RADIAL_BANDS": 4,
        "RIDGE_GUARD_BINS": 1,
        "RIDGE_N_SECTORS": 4,
        "RIDGE_N_PHASE": 72,
        "RIDGE_CORE_WIDTH_KPC": 0.25,
        "RIDGE_BROAD_RATIO": 3.0,
        "RIDGE_SHORTLIST_PER_FAMILY": 2,
        "RIDGE_MIN_HELD_OUT_SECTORS": 3,
        "M_COMPARE": np.array([2, 3, 4], dtype=int),
        "RIDGE_PITCH_GRID_DEG": np.array([-30.0, -20.0]),
        "R_REF_KPC": 1.0,
        "pd": pd,
    }
    if extra_namespace is not None:
        namespace.update(extra_namespace)
    module = ast.Module(body=definitions, type_ignores=[])
    exec(compile(ast.fix_missing_locations(module), str(NOTEBOOK), "exec"), namespace)
    return namespace


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


class NotebookContractTests(unittest.TestCase):
    def test_notebook_has_explanation_before_every_code_cell(self):
        cells = load_notebook()["cells"]
        self.assertEqual(cells[0]["cell_type"], "markdown")
        for index, cell in enumerate(cells):
            if cell["cell_type"] == "code":
                self.assertGreater(index, 0)
                self.assertEqual(cells[index - 1]["cell_type"], "markdown")
                explanation = "".join(cells[index - 1]["source"])
                self.assertGreaterEqual(len(explanation.split()), 35)

    def test_notebook_declares_required_sections_and_parameters(self):
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
            "Synthetic ridge-recovery",
            "data-derived skeleton",
            "Arm profile",
            "Interpretation limits",
        ]
        for item in required:
            self.assertIn(item, text)

    def test_executed_notebook_has_saved_result_contract(self):
        nb = load_notebook()
        code_cells = [cell for cell in nb["cells"] if cell["cell_type"] == "code"]
        self.assertTrue(code_cells)
        self.assertTrue(all(cell["execution_count"] is not None for cell in code_cells))
        text_outputs = []
        image_count = 0
        for cell in code_cells:
            for output in cell.get("outputs", []):
                if output.get("output_type") == "stream":
                    text = output.get("text", [])
                    text_outputs.extend([text] if isinstance(text, str) else text)
                data = output.get("data", {})
                image_count += int("image/png" in data)
        joined = "".join(text_outputs)
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

    def test_fit_domain_uses_every_valid_hii_bin(self):
        nb = load_notebook()
        text = "\n".join("".join(cell["source"]) for cell in nb["cells"])
        self.assertIn("USE_ALL_VALID_HII_BINS = True", text)
        self.assertIn("FIT_DOMAIN_ALL_VALID_HII_PASS", text)
        self.assertNotIn("OUTER_COVERAGE_QUANTILE", text)

    def test_ridge_map_uses_all_valid_hii_pixels_and_mask_normalization(self):
        function_names = notebook_function_names()
        for required_function in [
            "build_hii_pixel_geometry",
            "build_log_polar_contrast",
            "preprocess_logpolar_raw",
        ]:
            self.assertIn(required_function, function_names)
        build_source = function_source("build_log_polar_contrast")
        for required in [
            "np.quantile(u",
            "signed_residual",
            "raw_weighted_sum",
            "raw_coverage",
            "preprocess_logpolar_raw(raw_state, include_local=True)",
        ]:
            self.assertIn(required, build_source)
        preprocess_source = function_source("preprocess_logpolar_raw")
        for required in [
            'mode=("nearest", "wrap")',
            "local_ridge",
            "LOGPOLAR_AZIMUTH_BROAD_SIGMA_BINS",
        ]:
            self.assertIn(required, preprocess_source)
        code = notebook_code_source()
        for required in [
            "gradient_dex_per_kpc",
            "RIDGE_PIXEL_DOMAIN_PASS",
        ]:
            self.assertIn(required, code)
        self.assertNotIn("0.05 * np.nanmax(denominator)", code)
        self.assertNotIn("np.maximum(observed_log - background_log, 0.0)", code)

    def test_m234_negative_winding_configuration_contract(self):
        code = notebook_code_source()
        for required in [
            "M_COMPARE = np.array([2, 3, 4], dtype=int)",
            "RIDGE_PITCH_GRID_DEG = np.arange(-45.0, -4.9, 1.0)",
            "RIDGE_GUARD_BINS = 4",
            "RIDGE_NULL_BLOCK_SEEDS = np.array([4254, 5254, 6254, 7254], dtype=int)",
            "RIDGE_NULL_DRAWS_PER_BLOCK = 8",
            "RIDGE_NULL_TOTAL_DRAWS = int(len(RIDGE_NULL_BLOCK_SEEDS) * RIDGE_NULL_DRAWS_PER_BLOCK)",
        ]:
            self.assertIn(required, code)
        self.assertNotIn("M_CANDIDATES =", code)
        self.assertNotIn("np.arange(1, 7", code)
        self.assertNotIn("np.arange(5.0, 45.1", code)
        self.assertNotIn("RIDGE_N_NULL =", code)

    def test_conditional_m234_negative_winding_implementation_contract(self):
        self.assertIn("conditional_ridge_search", notebook_function_names())
        conditional_source = function_source("conditional_ridge_search")
        self.assertIn("leakage_controlled_negative_winding_ridge", conditional_source)
        self.assertIn("ridge_geometry_table", conditional_source)
        self.assertIn("ranked_candidates", conditional_source)

    def test_raw_logpolar_sector_folds_are_declared(self):
        function_names = notebook_function_names()
        for required_function in [
            "preprocess_logpolar_raw",
            "build_sector_fold_maps",
            "circular_dilate",
            "circular_erode",
            "sector_columns",
            "logpolar_fold_masks",
            "run_adversarial_leakage_check",
        ]:
            self.assertIn(required_function, function_names)
        preprocess_source = function_source("preprocess_logpolar_raw")
        for required in [
            "raw_weighted_sum",
            "raw_coverage",
            "allowed_phi",
            "include_local",
            'mode=("nearest", "wrap")',
        ]:
            self.assertIn(required, preprocess_source)
        fold_source = function_source("build_sector_fold_maps")
        self.assertIn("preprocess_logpolar_raw", fold_source)
        self.assertIn("include_local=False", fold_source)
        masks_source = function_source("logpolar_fold_masks")
        self.assertIn("circular_dilate", masks_source)
        self.assertIn("circular_erode", masks_source)
        code = notebook_code_source()
        self.assertIn("RIDGE_LEAKAGE_GUARD_PASS", code)

    def test_synthetic_logpolar_returns_raw_and_processed_maps(self):
        sentinel_pixels = object()
        sentinel_radial = object()
        sentinel_raw = {"raw": "state"}
        sentinel_processed = {"u": np.array([0.0]), "valid": np.ones((1, 1), bool)}

        def fake_synthetic_sfr_pixels(*args, **kwargs):
            return sentinel_pixels, sentinel_radial

        def fake_build_log_polar_contrast(pixels, radial_parameters):
            self.assertIs(pixels, sentinel_pixels)
            self.assertIs(radial_parameters, sentinel_radial)
            return sentinel_raw, sentinel_processed

        namespace = extracted_function_namespace(
            ["synthetic_logpolar"],
            extra_namespace={
                "synthetic_sfr_pixels": fake_synthetic_sfr_pixels,
                "build_log_polar_contrast": fake_build_log_polar_contrast,
            },
        )
        raw_result, processed_result = namespace["synthetic_logpolar"](
            3, -30.0, 0.7, 4254, noise_sigma=0.0)
        self.assertIs(raw_result, sentinel_raw)
        self.assertIs(processed_result, sentinel_processed)

    def test_three_synthetic_m234_cases_recover_with_fixed_tolerances(self):
        namespace = extracted_function_namespace(
            [
                "preprocess_logpolar_raw",
                "build_log_polar_contrast",
                "sector_columns",
                "circular_dilate",
                "circular_erode",
                "logpolar_fold_masks",
                "build_sector_fold_maps",
                "wrap_angle",
                "phase_row_histograms",
                "variable_width_ridge_response",
                "aggregate_radial_response",
                "ridge_curve_for_map",
                "conditional_ridge_search",
                "synthetic_sfr_pixels",
                "synthetic_logpolar",
            ],
            extra_namespace={
                "LOGPOLAR_N_U": 120,
                "LOGPOLAR_N_PHI": 360,
                "LOGPOLAR_N_RADIAL_BANDS": 24,
                "RIDGE_GUARD_BINS": 4,
                "RIDGE_N_SECTORS": 12,
                "RIDGE_N_PHASE": 360,
                "RIDGE_CORE_WIDTH_KPC": 0.25,
                "RIDGE_BROAD_RATIO": 3.0,
                "RIDGE_SHORTLIST_PER_FAMILY": 5,
                "RIDGE_MIN_HELD_OUT_SECTORS": 10,
                "RIDGE_PITCH_GRID_DEG": np.arange(-45.0, -4.9, 1.0),
                "R_REF_KPC": 1.0,
            },
        )
        cases = [(2, -22.0), (3, -30.0), (4, -38.0)]
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", RuntimeWarning)
            for case_index, (m_true, pitch_true) in enumerate(cases):
                raw_state, _ = namespace["synthetic_logpolar"](
                    m_true, pitch_true, 0.7, 4254 + case_index)
                geometry, _, _, cache_bytes = (
                    namespace["conditional_ridge_search"](
                        raw_state,
                        m_compare=np.array([m_true], dtype=int),
                        pitch_grid=np.arange(-45.0, -4.9, 1.0),
                    ))
                self.assertEqual(len(geometry), 1)
                recovered = geometry.iloc[0]
                phase_error_deg = np.rad2deg(abs(namespace["wrap_angle"](
                    recovered["Theta0"] - 0.7)))
                self.assertEqual(int(recovered["m_arms"]), m_true)
                self.assertLess(recovered["pitch_angle"], 0.0)
                self.assertLessEqual(
                    abs(recovered["pitch_angle"] - pitch_true), 2.0)
                self.assertLessEqual(phase_error_deg, 5.0)
                self.assertGreaterEqual(recovered["valid_held_out"], 10)
                self.assertGreater(recovered["validated_score"], 0.0)
                self.assertEqual(cache_bytes, 0)
        self.assertEqual(caught, [])

    def test_repeated_radii_are_rejected_before_endpoint_expansion(self):
        namespace = extracted_function_namespace([
            "preprocess_logpolar_raw",
            "build_log_polar_contrast",
        ])
        radius = np.array([1, 1, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9], dtype=float)
        repeated = pd.DataFrame({
            "radius_kpc": radius,
            "azimuth_rad": np.linspace(-2.5, 2.5, radius.size),
            "log_sfr": np.linspace(-2.0, -1.0, radius.size),
        })
        with self.assertRaisesRegex(ValueError, "strictly increasing"):
            namespace["build_log_polar_contrast"](
                repeated, {"lambda0_0": 0.1, "h_R": 2.0},
                n_u=4, n_phi=8)

    def test_preprocess_validates_raw_state_and_does_not_mutate_it(self):
        preprocess = extracted_function_namespace([
            "preprocess_logpolar_raw",
        ])["preprocess_logpolar_raw"]

        def valid_state():
            return {
                "raw_weighted_sum": np.arange(8, dtype=float).reshape(2, 4),
                "raw_coverage": np.full((2, 4), 2.0),
                "u": np.array([-0.5, 0.5]),
                "phi": np.linspace(-2.0, 2.0, 4),
                "u_edges": np.array([-1.0, 0.0, 1.0]),
                "phi_edges": np.linspace(-np.pi, np.pi, 5),
            }

        invalid_states = []
        state = valid_state()
        state["raw_weighted_sum"] = state["raw_weighted_sum"].ravel()
        state["raw_coverage"] = state["raw_coverage"].ravel()
        invalid_states.append(("one-dimensional maps", state, None))
        state = valid_state()
        state["raw_coverage"][0, 0] = -1.0
        invalid_states.append(("negative coverage", state, None))
        state = valid_state()
        state["raw_weighted_sum"][0, 0] = np.nan
        invalid_states.append(("NaN numerator", state, None))
        state = valid_state()
        state["u"] = state["u"][:1]
        invalid_states.append(("short u axis", state, None))
        state = valid_state()
        state["phi_edges"] = state["phi_edges"][:-1]
        invalid_states.append(("short phi edges", state, None))
        state = valid_state()
        state["phi"][0] = np.nan
        invalid_states.append(("NaN phi axis", state, None))
        invalid_states.append((
            "short allowed mask", valid_state(), np.array([True, False])))

        for label, state, allowed in invalid_states:
            with self.subTest(label=label):
                try:
                    preprocess(
                        state, allowed_phi=allowed, include_local=False,
                        smooth_sigma=1.0)
                except ValueError:
                    pass
                except Exception as exc:
                    self.fail(
                        f"Expected ValueError, got {type(exc).__name__}: {exc}")
                else:
                    self.fail("ValueError not raised")

        state = valid_state()
        original = {key: value.copy() for key, value in state.items()}
        preprocess(
            state, allowed_phi=np.array([True, False, True, False]),
            include_local=False)
        for key in state:
            np.testing.assert_array_equal(state[key], original[key])

    def test_build_logpolar_rejects_invalid_inputs_without_runtime_warning(self):
        namespace = extracted_function_namespace([
            "preprocess_logpolar_raw",
            "build_log_polar_contrast",
        ])
        build = namespace["build_log_polar_contrast"]

        def valid_table():
            return pd.DataFrame({
                "radius_kpc": np.array([0.5, 1.0, 1.5, 2.0]),
                "azimuth_rad": np.array([-2.0, -0.5, 0.5, 2.0]),
                "log_sfr": np.array([-2.0, -1.8, -1.6, -1.4]),
            })

        cases = []
        table = valid_table()
        table.loc[0, "radius_kpc"] = 0.0
        cases.append(("zero radius", table, {"lambda0_0": 0.1, "h_R": 2.0}))
        table = valid_table()
        table.loc[0, "azimuth_rad"] = np.nan
        cases.append(("NaN azimuth", table, {"lambda0_0": 0.1, "h_R": 2.0}))
        table = valid_table()
        table.loc[0, "log_sfr"] = np.nan
        cases.append(("NaN log SFR", table, {"lambda0_0": 0.1, "h_R": 2.0}))
        cases.append(("zero lambda0", valid_table(), {"lambda0_0": 0.0, "h_R": 2.0}))
        cases.append(("zero h_R", valid_table(), {"lambda0_0": 0.1, "h_R": 0.0}))
        cases.append(("infinite h_R", valid_table(), {"lambda0_0": 0.1, "h_R": np.inf}))

        for label, table, parameters in cases:
            with self.subTest(label=label):
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always", RuntimeWarning)
                    with self.assertRaises(ValueError):
                        build(table, parameters, n_u=1, n_phi=4)
                runtime_warnings = [
                    item for item in caught
                    if issubclass(item.category, RuntimeWarning)
                ]
                self.assertEqual(runtime_warnings, [])

    def test_conditional_search_uses_raw_fold_field_without_cache(self):
        self.assertIn("conditional_ridge_search", notebook_function_names())
        conditional_source = function_source("conditional_ridge_search")
        ridge_curve_source = function_source("ridge_curve_for_map")
        self.assertIn("ridge_curve_for_map", conditional_source)
        self.assertIn('field="radial_residual"', ridge_curve_source)
        self.assertIn("histogram_cache_payload_bytes = 0", conditional_source)
        self.assertNotIn('field="local_ridge"', ridge_curve_source)
        self.assertNotIn("histogram_cache = {}", conditional_source)

    def test_phase_histograms_and_variable_width_response_are_scale_invariant(self):
        namespace = extracted_function_namespace([
            "phase_row_histograms",
            "variable_width_ridge_response",
        ])
        n_u = 8
        n_phi = 36
        u = np.linspace(np.log(0.8), np.log(5.0), n_u)
        phi = np.linspace(-np.pi, np.pi, n_phi, endpoint=False)
        uu, pp = np.meshgrid(u, phi, indexing="ij")
        field = np.cos(
            (2.0 / np.tan(np.deg2rad(-20.0))) * uu - 2.0 * pp)
        logpolar_map = {
            "u": u,
            "phi": phi,
            "radial_residual": field,
            "coverage": np.full_like(field, 5.0),
            "valid": np.ones_like(field, dtype=bool),
        }
        weighted, support = namespace["phase_row_histograms"](
            logpolar_map, 2, -20.0, n_phase=72)
        self.assertEqual(weighted.shape, (n_u, 72))
        self.assertEqual(support.shape, (n_u, 72))
        response = namespace["variable_width_ridge_response"](
            weighted, support, u, 2, -20.0, n_radial_bands=4)
        scaled = namespace["variable_width_ridge_response"](
            6.0 * weighted, 6.0 * support, u, 2, -20.0,
            n_radial_bands=4)
        np.testing.assert_allclose(
            response["response"], scaled["response"],
            equal_nan=True, atol=1.0e-12)
        self.assertGreaterEqual(response["low_clamp_count"], 0)
        self.assertGreaterEqual(response["high_clamp_count"], 0)

    def test_conditional_search_returns_one_geometry_per_requested_mode(self):
        namespace = extracted_function_namespace([
            "preprocess_logpolar_raw",
            "sector_columns",
            "circular_dilate",
            "circular_erode",
            "logpolar_fold_masks",
            "build_sector_fold_maps",
            "phase_row_histograms",
            "variable_width_ridge_response",
            "aggregate_radial_response",
            "ridge_curve_for_map",
            "conditional_ridge_search",
        ])
        n_u = 8
        n_phi = 72
        u_edges = np.linspace(np.log(0.8), np.log(5.0), n_u + 1)
        phi_edges = np.linspace(-np.pi, np.pi, n_phi + 1)
        u = 0.5 * (u_edges[:-1] + u_edges[1:])
        phi = 0.5 * (phi_edges[:-1] + phi_edges[1:])
        uu, pp = np.meshgrid(u, phi, indexing="ij")
        signal = np.zeros((n_u, n_phi), dtype=float)
        for m_arms, pitch_angle, theta0 in [
                (2, -20.0, 0.5), (3, -30.0, 1.0)]:
            phase = ((m_arms / np.tan(np.deg2rad(pitch_angle))) * uu
                     - m_arms * pp + theta0)
            wrapped = np.angle(np.exp(1j * phase))
            signal += 0.7 * np.exp(-0.5 * (wrapped / 0.18)**2)
        coverage = np.full_like(signal, 40.0)
        raw_state = {
            "raw_weighted_sum": coverage * signal,
            "raw_coverage": coverage,
            "u": u,
            "phi": phi,
            "u_edges": u_edges,
            "phi_edges": phi_edges,
        }
        geometry, candidates, folds, cache_bytes = (
            namespace["conditional_ridge_search"](
                raw_state,
                m_compare=np.array([2, 3]),
                pitch_grid=np.array([-30.0, -20.0]),
                guard_bins=1,
            ))
        self.assertEqual(geometry["m_arms"].tolist(), [2, 3])
        self.assertTrue((geometry["pitch_angle"] < 0.0).all())
        self.assertTrue((geometry["geometry_source"]
                         == "leakage_controlled_negative_winding_ridge").all())
        self.assertTrue((geometry["valid_held_out"] >= 3).all())
        self.assertTrue((geometry["validated_score"] > 0.0).all())
        self.assertEqual(set(candidates["m_arms"]), {2, 3})
        self.assertEqual(set(folds["m_arms"]), {2, 3})
        self.assertEqual(cache_bytes, 0)

    def test_conditional_search_rejects_noninteger_modes_and_nonnegative_pitch(self):
        search = extracted_function_namespace([
            "conditional_ridge_search",
        ])["conditional_ridge_search"]
        with self.assertRaises(ValueError):
            search({}, m_compare=np.array([2.5]), pitch_grid=np.array([-20.0]))
        with self.assertRaises(ValueError):
            search({}, m_compare=np.array([2]), pitch_grid=np.array([20.0]))

    def test_blocked_nulls_are_descriptive_not_selective(self):
        self.assertIn("build_blocked_null_diagnostics", notebook_function_names())
        blocked_null_source = function_source("build_blocked_null_diagnostics")
        for required in [
            "null_block_summary",
            "pooled_null_summary",
            "descriptive_not_model_selection",
        ]:
            self.assertIn(required, blocked_null_source)
        self.assertNotIn("accepted.iloc[0].to_dict()", blocked_null_source)
        self.assertIn("RIDGE_NULL_BLOCKS_COMPLETE", notebook_code_source())

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

    def test_three_fixed_geometry_ktz_fits_are_declared(self):
        self.assertIn("fit_ktz_profile_fixed_geometry", notebook_function_names())
        fit_source = function_source("fit_ktz_profile_fixed_geometry")
        self.assertIn('"geometry_source": geometry["geometry_source"]', fit_source)
        code = notebook_code_source()
        for required in [
            "conditional_fits",
            "conditional_fit_tables",
            "profile_maxima_by_m",
            "KTZ_M234_PROFILE_FITS_COMPLETE",
            "HARMONIC_M234_MAXIMA_COMPLETE",
        ]:
            self.assertIn(required, code)

    def test_m234_comparative_figures_and_table_are_declared(self):
        code = notebook_code_source()
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
            self.assertIn(required, code)

    def test_adversarial_sector_signal_has_zero_training_field(self):
        namespace = extracted_function_namespace([
            "preprocess_logpolar_raw",
            "sector_columns",
            "circular_dilate",
            "circular_erode",
            "logpolar_fold_masks",
            "build_sector_fold_maps",
            "run_adversarial_leakage_check",
        ])
        result = namespace["run_adversarial_leakage_check"](
            n_u=12,
            n_phi=48,
            n_sectors=4,
            held_sector=1,
            guard_bins=2,
        )
        self.assertEqual(result["raw_training_signal_l1"], 0.0)
        self.assertLessEqual(result["max_abs_training_residual"], 1.0e-14)
        for required in [
            "train_valid_count",
            "max_abs_residual_difference",
            "max_abs_coverage_difference",
            "known_train_valid_count",
            "expected_training_residual",
            "mean_known_training_residual",
            "max_abs_known_training_error",
        ]:
            self.assertIn(required, result)
        self.assertGreater(result["train_valid_count"], 0)
        self.assertLessEqual(result["max_abs_residual_difference"], 1.0e-14)
        self.assertLessEqual(result["max_abs_coverage_difference"], 1.0e-14)
        self.assertGreater(result["known_train_valid_count"], 0)
        self.assertAlmostEqual(result["expected_training_residual"], 0.25)
        self.assertAlmostEqual(result["mean_known_training_residual"], 0.25)
        self.assertLessEqual(result["max_abs_known_training_error"], 1.0e-14)

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


if __name__ == "__main__":
    unittest.main()
