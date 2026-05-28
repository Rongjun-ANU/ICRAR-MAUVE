import ast
import glob
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "make_ngist_masks_from_catalogs_VRI.py"


def load_helpers(required_names):
    tree = ast.parse(SCRIPT.read_text(), filename=str(SCRIPT))
    nodes = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in set(required_names)
    ]
    namespace = {"glob": glob}
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


if __name__ == "__main__":
    for name, func in sorted(globals().items()):
        if name.startswith("test_") and callable(func):
            func()
            print(f"PASS {name}")
