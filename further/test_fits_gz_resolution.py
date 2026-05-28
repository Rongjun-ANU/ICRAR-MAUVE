import ast
import tempfile
from pathlib import Path

import numpy as np
from astropy.io import fits


ROOT = Path(__file__).resolve().parent


OPTIONAL_HELPERS = {
    "_fits_path_alternates",
    "_glob_candidate_patterns",
    "_glob_matches_with_fits_gz",
}


def load_helpers(script_name, required_names):
    source_path = ROOT / script_name
    tree = ast.parse(source_path.read_text(), filename=str(source_path))
    wanted = set(required_names) | OPTIONAL_HELPERS
    nodes = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in wanted
    ]
    namespace = {"Path": Path}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(source_path), "exec"), namespace)
    return namespace


def test_astropy_reads_gzip_fits_without_manual_unzip():
    with tempfile.TemporaryDirectory() as tmpdir:
        path = Path(tmpdir) / "tiny.fits.gz"
        fits.PrimaryHDU(data=np.array([[1.5]])).writeto(path)

        with fits.open(path) as hdul:
            assert float(hdul[0].data[0, 0]) == 1.5


def test_mass_resolver_finds_compressed_exact_counterpart():
    helpers = load_helpers(
        "Mass.py",
        ["_unique_paths", "_has_glob_chars", "_candidate_paths", "resolve_existing_path"],
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        root = Path(tmpdir)
        compressed = root / "IC3392_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk.fits.gz"
        compressed.touch()

        resolved = helpers["resolve_existing_path"](
            "datacube FITS",
            root / "IC3392_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk.fits",
        )

        assert resolved == compressed.resolve()


def test_mass_resolver_finds_compressed_glob_counterpart():
    helpers = load_helpers(
        "Mass.py",
        ["_unique_paths", "_has_glob_chars", "_candidate_paths", "resolve_existing_path"],
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        root = Path(tmpdir)
        compressed = root / "IC3392_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk_v2.fits.gz"
        compressed.touch()

        resolved = helpers["resolve_existing_path"](
            "datacube FITS",
            root / "IC3392_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk*.fits",
        )

        assert resolved == compressed.resolve()


def test_sfr_keyworded_lookup_accepts_compressed_fits():
    helpers = load_helpers(
        "SFR+Z.py",
        [
            "_unique_paths",
            "_has_glob_chars",
            "_candidate_paths",
            "_input_search_dirs",
            "extract_shared_keyword",
            "collect_keyworded_input_paths",
            "find_keyworded_input_path",
        ],
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        root = Path(tmpdir)
        compressed = root / "NGC4064_gas_bin_maps.fits.gz"
        compressed.touch()

        resolved = helpers["find_keyworded_input_path"](
            root, None, Path("."), "NGC4064", ["gas_bin_maps.fits"], ""
        )

        assert resolved == compressed.resolve()


def test_output_path_from_compressed_input_stays_uncompressed_fits():
    for script_name in ("Mass.py", "SFR+Z.py"):
        helpers = load_helpers(script_name, ["build_further_output_path"])
        source = Path("/tmp") / "NGC4064_gas_bin_maps.fits.gz"

        output = helpers["build_further_output_path"](source)

        assert output == Path("/tmp") / "NGC4064_gas_bin_maps_further.fits"


def test_proxy_resolver_finds_compressed_counterpart():
    helpers = load_helpers("proxy_EWHa.py", ["_unique_paths", "resolve_existing_path"])
    with tempfile.TemporaryDirectory() as tmpdir:
        root = Path(tmpdir)
        compressed = root / "NGC4064_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk.fits.gz"
        compressed.touch()

        resolved = helpers["resolve_existing_path"](
            "continuum cube FITS",
            root / "NGC4064_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk.fits",
        )

        assert resolved == compressed.resolve()


def test_proxy_shell_delegates_input_resolution_to_python():
    shell_source = (ROOT / "proxy_EWHa.sh").read_text()

    assert "resolve_first_existing" not in shell_source
    assert "--bin-file" not in shell_source
    assert "--gas-file" not in shell_source
    assert "--cont-file" not in shell_source
    assert '--root "$ROOT_PRODUCT_BASE"' in shell_source
    assert '--fallback-root "$ROOT_LOCAL"' in shell_source
    assert '--cube-root "$CUBE_ROOT"' in shell_source


if __name__ == "__main__":
    for name, func in sorted(globals().items()):
        if name.startswith("test_") and callable(func):
            func()
            print(f"PASS {name}")
