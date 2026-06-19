import ast
import glob
import os
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "make_ngist_masks_from_catalogs_VRI.py"
OPTIONAL_HELPERS = {
    "is_bare_galaxy_id",
    "input_patterns_for_argument",
}


def load_helpers(required_names):
    tree = ast.parse(SCRIPT.read_text(), filename=str(SCRIPT))
    wanted = set(required_names) | OPTIONAL_HELPERS
    nodes = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in wanted
    ]
    namespace = {"glob": glob, "os": os}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(SCRIPT), "exec"), namespace)
    missing = [name for name in required_names if name not in namespace]
    if missing:
        raise AssertionError(f"Missing helper(s): {', '.join(missing)}")
    return namespace


def test_fits_pattern_also_matches_gzip_counterpart():
    helpers = load_helpers(["fits_path_patterns", "expand_fits_input_patterns"])
    with tempfile.TemporaryDirectory() as tmpdir:
        root = Path(tmpdir)
        compressed = root / "NGC4064_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk_VRI.fits.gz"
        compressed.touch()

        matches = helpers["expand_fits_input_patterns"](
            [str(root / "*_DATACUBE*_VRI.fits")]
        )

        assert matches == [str(compressed)]


def test_bare_phangs_galid_matches_native_vri_fits():
    helpers = load_helpers(["fits_path_patterns", "expand_fits_input_patterns"])
    with tempfile.TemporaryDirectory() as tmpdir:
        root = Path(tmpdir)
        product = root / "NGC4254_PHANGS_DATACUBE_native_VRI.fits"
        product.touch()

        matches = helpers["expand_fits_input_patterns"]([str(root / "NGC4254")])

        assert matches == [str(product)]


def test_bare_phangs_galid_matches_native_vri_fits_gz():
    helpers = load_helpers(["fits_path_patterns", "expand_fits_input_patterns"])
    with tempfile.TemporaryDirectory() as tmpdir:
        root = Path(tmpdir)
        product = root / "NGC4321_PHANGS_DATACUBE_native_VRI.fits.gz"
        product.touch()

        matches = helpers["expand_fits_input_patterns"]([str(root / "NGC4321")])

        assert matches == [str(product)]


def test_bare_local_galid_matches_v3tk_vri_fits_gz_in_place():
    helpers = load_helpers(["fits_path_patterns", "expand_fits_input_patterns"])
    with tempfile.TemporaryDirectory() as tmpdir:
        root = Path(tmpdir)
        product = root / "NGC4380_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk_VRI.fits.gz"
        product.touch()

        matches = helpers["expand_fits_input_patterns"]([str(root / "NGC4380")])

        assert matches == [str(product)]
        assert product.exists()


def test_bare_multiple_phangs_galids_preserve_input_order():
    helpers = load_helpers(["fits_path_patterns", "expand_fits_input_patterns"])
    with tempfile.TemporaryDirectory() as tmpdir:
        root = Path(tmpdir)
        products = [
            root / "NGC4254_PHANGS_DATACUBE_native_VRI.fits",
            root / "NGC4321_PHANGS_DATACUBE_native_VRI.fits",
            root / "NGC4535_PHANGS_DATACUBE_native_VRI.fits",
        ]
        for product in products:
            product.touch()

        matches = helpers["expand_fits_input_patterns"](
            [str(root / "NGC4254"), str(root / "NGC4321"), str(root / "NGC4535")]
        )

        assert matches == [str(product) for product in products]


def test_safe_base_id_strips_phangs_native_suffix():
    helpers = load_helpers(["safe_base_id"])

    assert helpers["safe_base_id"]("NGC4254_PHANGS_DATACUBE_native_VRI.fits") == "NGC4254"
    assert (
        helpers["safe_base_id"]("NGC4064_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk_VRI.fits")
        == "NGC4064"
    )


if __name__ == "__main__":
    for name, func in sorted(globals().items()):
        if name.startswith("test_") and callable(func):
            func()
            print(f"PASS {name}")
