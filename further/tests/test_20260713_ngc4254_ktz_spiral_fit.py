import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "20260713_NGC4254_KTZ_spiral_fit.ipynb"


def load_notebook():
    assert NOTEBOOK.exists(), f"Notebook has not been created: {NOTEBOOK}"
    return json.loads(NOTEBOOK.read_text(encoding="utf-8"))


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

    def test_null_calibration_reuses_one_deterministic_shift_stream(self):
        text = notebook_source()
        self.assertIn("null_rng = np.random.default_rng(RNG_SEED)", text)
        self.assertIn("scramble_logpolar_rows(logpolar_map, null_rng)", text)
        self.assertNotIn("RNG_SEED + 1000 + null_index", text)

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


if __name__ == "__main__":
    unittest.main()
