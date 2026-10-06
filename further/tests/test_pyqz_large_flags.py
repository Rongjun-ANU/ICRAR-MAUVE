"""Preserve concatenated pyqz QC flags through the SFR FITS pipeline."""

import ast
from pathlib import Path

import numpy as np
from astropy.io import fits

from model_grid_diagnostics import BinSpectra, broadcast_bin_results, run_pyqz_spectra


def _sfr_definitions():
    source_path = Path(__file__).resolve().parents[1] / "SFR+Z.py"
    tree = ast.parse(source_path.read_text())
    names = {
        "selected_model_field",
        "append_ordered_image",
        "append_named_hii_sf_pair",
    }
    definitions = ast.Module(
        body=[node for node in tree.body if isinstance(node, ast.FunctionDef)
              and node.name in names],
        type_ignores=[],
    )
    scope = {"np": np, "fits": fits, "NOT_EVALUATED_FLAG": -99,
             "gas_header": fits.Header({"CTYPE1": "RA---TAN"})}
    exec(compile(definitions, str(source_path), "exec"), scope)
    return tree, scope


def test_pyqz_batch_retains_flag_91234_from_canfar_failure():
    class RawPyqzResult:
        def get_global_qz(self, *args, **kwargs):
            columns = ["<LogQ{KDE}>", "err(LogQ{KDE})", "<gas[O]+12{KDE}>",
                       "err(gas[O]+12{KDE})", "flag", "rs_offgrid"]
            return np.array([[7.2, 0.1, 8.4, 0.08, 91234, 2.5]]), columns

    spectra = BinSpectra(
        bin_ids=np.array([0]),
        fluxes=np.array([[10., 5., 28.6, 4., 3., 2.]]),
        errors=np.full((1, 6), 0.5),
        pixel_counts=np.array([1]),
    )
    result = run_pyqz_spectra(RawPyqzResult(), spectra)
    assert result.results["flag"].tolist() == [91234]
    assert result.results["flag"].dtype == np.int32
    assert result.results["valid"].dtype == np.int16


def test_broadcast_retains_large_flag_and_not_evaluated_sentinel():
    result = broadcast_bin_results(
        np.array([[0, 0, 1, -1]]), [0],
        {"flag": np.array([91234], dtype=np.int32)}, integer_fields={"flag"},
    )
    np.testing.assert_array_equal(result["flag"], [[91234, 91234, -99, -99]])
    assert result["flag"].dtype == np.int32


def test_sfr_region_mask_and_actual_flag_writer_retain_large_flags(tmp_path):
    tree, scope = _sfr_definitions()
    flag_map = np.array([[91234, 1234, -99]], dtype=np.int32)
    selection = np.array([[True, False, True]])
    selected = scope["selected_model_field"](flag_map, selection, integer=True)
    np.testing.assert_array_equal(selected, [[91234, -99, -99]])

    scope["ordered_hdul"] = fits.HDUList([fits.PrimaryHDU()])
    scope["MODEL_OUTPUT_MAPS"] = {"PYQZ_FLAG_HII": selected, "PYQZ_FLAG_SF": selected}
    scope["CARTA_DIMENSIONLESS_BUNIT"] = "1"
    scope["pyqz_references"] = ("pyqz raw flag",)
    writer = next(
        node for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
        and node.func.id == "append_named_hii_sf_pair"
        and len(node.args) > 1 and isinstance(node.args[1], ast.Constant)
        and node.args[1].value == "PYQZ_FLAG_HII"
    )
    eval(compile(ast.Expression(writer), "SFR+Z.py", "eval"), scope)
    output = tmp_path / "flags.fits"
    scope["ordered_hdul"].writeto(output)
    with fits.open(output) as hdul:
        for name in ("PYQZ_FLAG_HII", "PYQZ_FLAG_SF"):
            np.testing.assert_array_equal(hdul[name].data, selected)
            assert hdul[name].header["BITPIX"] == 32
            assert hdul[name].header["CTYPE1"] == "RA---TAN"


def test_sfr_schema_accepts_int32_pyqz_flags_and_int16_valid_masks():
    tree, scope = _sfr_definitions()
    assignment = next(
        node for node in ast.walk(tree) if isinstance(node, ast.Assign)
        and any(isinstance(target, ast.Name) and target.id == "expected_dtype"
                for target in node.targets)
    )
    scope["integer_model_names"] = {"PYQZ_FLAG_HII", "PYQZ_VALID_HII", "NB_FLAG_HII"}
    for name, expected in [("PYQZ_FLAG_HII", np.int32), ("PYQZ_VALID_HII", np.int16),
                           ("NB_FLAG_HII", np.int16), ("O_H_PYQZ_HII", np.float64)]:
        scope["name"] = name
        dtype = eval(compile(ast.Expression(assignment.value), "SFR+Z.py", "eval"), scope)
        assert dtype == np.dtype(expected)
