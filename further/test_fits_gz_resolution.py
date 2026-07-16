import ast
import importlib.util
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
OPTIONAL_CONSTANTS = {
    "PHANGS_NATIVE_GALAXIES",
}


def load_helpers(script_name, required_names):
    source_path = ROOT / script_name
    tree = ast.parse(source_path.read_text(), filename=str(source_path))
    wanted = set(required_names) | OPTIONAL_HELPERS
    nodes = []
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name in wanted:
            nodes.append(node)
        elif isinstance(node, ast.Assign):
            target_names = {
                target.id for target in node.targets if isinstance(target, ast.Name)
            }
            if target_names & OPTIONAL_CONSTANTS:
                nodes.append(node)
    namespace = {"Path": Path, "np": np}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(source_path), "exec"), namespace)
    return namespace


def load_script_module(script_name):
    source_path = ROOT / script_name
    spec = importlib.util.spec_from_file_location(source_path.stem, source_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module spec for {source_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


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


def test_proxy_flux_density_to_cgs_handles_phangs_bunit_with_parenthesized_exponents():
    helpers = load_helpers("proxy_EWHa.py", ["flux_density_to_cgs"])
    bunit = "10**(-20)angstrom**(-1).cm**(-2).erg.s**(-1)"

    converted = helpers["flux_density_to_cgs"](np.array([2.5]), bunit)

    assert np.allclose(converted, np.array([2.5e-20]))


def test_mass_uses_phangs_native_cube_names_for_phangs_galaxies():
    helpers = load_helpers("Mass.py", ["datacube_filenames_for_galaxy"])

    for galid in ("NGC4254", "NGC4321", "NGC4535"):
        names = helpers["datacube_filenames_for_galaxy"](galid)

        assert names == [f"{galid}_PHANGS_DATACUBE_native.fits"]


def test_proxy_uses_phangs_native_cube_names_for_phangs_galaxies():
    helpers = load_helpers("proxy_EWHa.py", ["datacube_filenames_for_galaxy"])

    for galid in ("NGC4254", "NGC4321", "NGC4535"):
        names = helpers["datacube_filenames_for_galaxy"](galid)

        assert names == [f"{galid}_PHANGS_DATACUBE_native.fits"]


def test_shell_defaults_include_phangs_and_keep_cli_override():
    for shell_name in ("Mass.sh", "SFR.sh", "proxy_EWHa.sh"):
        shell_source = (ROOT / shell_name).read_text()

        for galid in ("NGC4254", "NGC4321", "NGC4535"):
            assert f'"{galid}"' in shell_source
        assert '[[ $# -gt 0 ]] && GALAXIES=("$@")' in shell_source


def test_cube_consuming_scripts_accept_phangs_cube_root_override():
    for script_name in ("Mass.py", "proxy_EWHa.py"):
        script_source = (ROOT / script_name).read_text()
        assert "--phangs-cube-root" in script_source

    for shell_name in ("Mass.sh", "proxy_EWHa.sh"):
        shell_source = (ROOT / shell_name).read_text()
        assert "PHANGS_CUBE_ROOT" in shell_source
        assert "--phangs-cube-root" in shell_source


def test_mass_filter_support_threshold_accepts_phangs_ao_gap():
    source_path = ROOT / "Mass.py"
    tree = ast.parse(source_path.read_text(), filename=str(source_path))

    threshold = None
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        if not any(
            isinstance(target, ast.Name) and target.id == "min_filter_support_fraction"
            for target in node.targets
        ):
            continue
        threshold = ast.literal_eval(node.value)
        break

    assert threshold is not None
    assert threshold <= 0.91


def test_proxy_combined_ngc4567_8_redshift_uses_member_mean():
    proxy_ewha = load_script_module("proxy_EWHa.py")
    with tempfile.TemporaryDirectory() as tmpdir:
        redshift_path = Path(tmpdir) / "new_redshifts"
        redshift_path.write_text(
            "NGC4567 0.007495\n"
            "NGC4568 0.007446\n",
            encoding="utf-8",
        )

        z = proxy_ewha.lookup_redshift("NGC4567_8", redshift_path)

        assert np.isclose(z, (0.007495 + 0.007446) / 2.0)


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
