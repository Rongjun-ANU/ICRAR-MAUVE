import ast
import json
from pathlib import Path
import unittest

import numpy as np
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


def extracted_function_namespace(function_names):
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
        "RIDGE_GUARD_BINS": 4,
        "RIDGE_N_SECTORS": 12,
    }
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
        code = notebook_code_source()
        self.assertIn("ridge_geometry_table", code)
        self.assertIn("RIDGE_M234_REAL_COMPLETE", code)

    def test_leakage_controlled_fold_contract(self):
        function_names = notebook_function_names()
        for required_function in [
            "preprocess_logpolar_raw",
            "build_sector_fold_maps",
            "circular_dilate",
            "circular_erode",
            "conditional_ridge_search",
        ]:
            self.assertIn(required_function, function_names)
        preprocess_source = function_source("preprocess_logpolar_raw")
        for required in [
            "raw_weighted_sum",
            "raw_coverage",
        ]:
            self.assertIn(required, preprocess_source)
        conditional_source = function_source("conditional_ridge_search")
        self.assertIn('field="radial_residual"', conditional_source)
        self.assertIn("histogram_cache_payload_bytes = 0", conditional_source)
        self.assertNotIn('field="local_ridge"', conditional_source)
        self.assertNotIn("histogram_cache = {}", conditional_source)
        self.assertIn("RIDGE_LEAKAGE_GUARD_PASS", notebook_code_source())

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
