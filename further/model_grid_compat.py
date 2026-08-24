"""Runtime compatibility adapters for the installed model-grid packages.

The adapters in this module are deliberately version-gated and process-local.
They do not modify files in ``site-packages``.  This keeps the MAUVE workflow
portable to CANFAR while making unsupported package revisions fail loudly.
"""

from __future__ import annotations

import hashlib
import importlib
import importlib.util
from importlib import metadata
from pathlib import Path
import re
import shutil
import sys
import tempfile
from typing import Any, Sequence

import numpy as np


NEBULABAYES_VERSION = "0.9.9"
NEBULABAYES_GRID_PARAMETERS = ["log U", "log P/k", "12 + log O/H"]
PYQZ_VERSION = "0.8.4.0"
PYQZ_SOURCE_HASHES = {
    "__init__.py": "654bcee91a4abc68dbbcd2d36887d88a6852811a91e9da6c0ca44d47400b8207",
    "pyqz.py": "71fac51534eac1d32b40a068bfad9491f0c4754cc222b8b2b93d84938efe9cdb",
    "pyqz_metadata.py": "d7fe6c5ae1b54e8f9f1038f9b465621565e903ddeaa54a266283707a9119c9c8",
    "pyqz_tools.py": "2deb6a665831e51505d2538e18163e5a331402aa88a924e9f24f0dac9151755b",
}

_PYQZ_RUNTIME_DIRECTORY: tempfile.TemporaryDirectory[str] | None = None
_PYQZ_MODULE: Any = None


def _require_distribution_version(distribution: str, expected: str) -> None:
    """Fail before applying a compatibility patch to an unknown release."""

    try:
        installed = metadata.version(distribution)
    except metadata.PackageNotFoundError as exc:
        raise RuntimeError(
            f"Required package {distribution} is not installed in this interpreter."
        ) from exc
    if installed != expected:
        raise RuntimeError(
            f"The runtime adapter supports {distribution}=={expected}, but this "
            f"interpreter has {distribution}=={installed}."
        )


def _simps_compat(
    y: Any,
    x: Any = None,
    dx: float = 1.0,
    axis: int = -1,
    even: str | None = None,
) -> Any:
    """Expose SciPy's current Simpson implementation under its old API name."""

    del even
    from scipy import integrate

    return integrate.simpson(y, x=x, dx=dx, axis=axis)


def load_nebulabayes_compat() -> Any:
    """Import NebulaBayes 0.9.9 under current NumPy and SciPy releases."""

    _require_distribution_version("NebulaBayes", NEBULABAYES_VERSION)

    from scipy import integrate

    if "cumtrapz" not in integrate.__dict__:
        integrate.cumtrapz = integrate.cumulative_trapezoid
    if "simps" not in integrate.__dict__:
        integrate.simps = _simps_compat

    numpy_aliases = {
        "str": str,
        "float": float,
        "int": int,
        "product": np.prod,
    }
    for name, replacement in numpy_aliases.items():
        if name not in np.__dict__:
            setattr(np, name, replacement)
    # NebulaBayes 0.9.9 still calls the deprecated spelling repeatedly while
    # marginalising.  NumPy documents ``trapezoid`` as its direct replacement.
    np.trapz = np.trapezoid

    return importlib.import_module("NebulaBayes")


def load_nebulabayes_hii_grid() -> tuple[Any, list[str], Path]:
    """Load the bundled HII grid from its table HDU using float64 axes."""

    _require_distribution_version("NebulaBayes", NEBULABAYES_VERSION)

    from astropy.table import Table

    distribution = metadata.distribution("NebulaBayes")
    grid_path = Path(
        distribution.locate_file("NebulaBayes/grids/NB_HII_grid.fits.gz")
    )
    if not grid_path.is_file():
        raise RuntimeError(f"NebulaBayes HII grid is missing: {grid_path}")

    grid = Table.read(grid_path, hdu=1).to_pandas()
    for parameter in NEBULABAYES_GRID_PARAMETERS:
        if parameter not in grid.columns:
            raise RuntimeError(
                f"NebulaBayes HII grid is missing parameter column {parameter!r}."
            )
        grid[parameter] = grid[parameter].astype(np.float64)

    return grid, list(NEBULABAYES_GRID_PARAMETERS), grid_path


def make_nebulabayes_hii_model(
    line_list: Sequence[str],
    *,
    interpd_grid_shape: Sequence[int],
    interp_order: int = 1,
    grid_error: float = 0.10,
) -> Any:
    """Initialise the bundled NebulaBayes HII model through the safe loader."""

    nebula_bayes = load_nebulabayes_compat()
    grid, grid_parameters, _ = load_nebulabayes_hii_grid()
    return nebula_bayes.NB_Model(
        grid,
        grid_params=grid_parameters,
        line_list=list(line_list),
        interpd_grid_shape=list(interpd_grid_shape),
        interp_order=interp_order,
        grid_error=grid_error,
    )


def _replace_exact(
    source: str,
    old: str,
    new: str,
    label: str,
    *,
    expected: int = 1,
) -> str:
    count = source.count(old)
    if count != expected:
        raise RuntimeError(
            f"pyqz runtime patch {label!r} expected {expected} source matches, "
            f"but found {count}."
        )
    return source.replace(old, new)


def _replace_regex(
    source: str,
    pattern: str,
    replacement: str,
    label: str,
    *,
    expected: int,
) -> str:
    updated, count = re.subn(pattern, replacement, source, flags=re.MULTILINE)
    if count != expected:
        raise RuntimeError(
            f"pyqz runtime patch {label!r} expected {expected} source matches, "
            f"but found {count}."
        )
    return updated


def _validate_pyqz_installation() -> Path:
    _require_distribution_version("pyqz", PYQZ_VERSION)
    distribution = metadata.distribution("pyqz")
    source_directory = Path(distribution.locate_file("pyqz")).resolve()
    for filename, expected_hash in PYQZ_SOURCE_HASHES.items():
        source_path = source_directory / filename
        if not source_path.is_file():
            raise RuntimeError(f"Installed pyqz source file is missing: {source_path}")
        actual_hash = hashlib.sha256(source_path.read_bytes()).hexdigest()
        if actual_hash != expected_hash:
            raise RuntimeError(
                f"Installed pyqz source hash changed for {filename}: "
                f"expected {expected_hash}, found {actual_hash}. Refusing to apply "
                "the version-specific runtime patch."
            )
    return source_directory


def _patch_pyqz_init(package_directory: Path) -> None:
    path = package_directory / "__init__.py"
    source = path.read_text(encoding="utf-8")
    source = _replace_exact(
        source,
        "from pyqz import * # So that users only need to do import pyqz\n",
        "from .pyqz import * # Runtime Python-3 overlay\n",
        "package relative import",
    )
    source = _replace_exact(
        source,
        "from pyqz_metadata import __version__\n",
        "from .pyqz_metadata import __version__\n",
        "metadata relative import",
    )
    path.write_text(source, encoding="utf-8")


def _patch_pyqz_tools(package_directory: Path) -> None:
    path = package_directory / "pyqz_tools.py"
    source = path.read_text(encoding="utf-8")
    source = _replace_exact(
        source,
        "import pyqz_metadata as pyqzm\nfrom pyqz_metadata import __version__\n",
        "from . import pyqz_metadata as pyqzm\n"
        "from .pyqz_metadata import __version__\n",
        "tools relative imports",
    )
    source = _replace_regex(
        source,
        r"np\.int(?![A-Za-z0-9_])",
        "int",
        "tools NumPy int alias",
        expected=4,
    )
    source = _replace_regex(
        source,
        r"np\.float(?![A-Za-z0-9_])",
        "float",
        "tools NumPy float alias",
        expected=2,
    )
    source = _replace_regex(
        source,
        r"np\.str(?![A-Za-z0-9_])",
        "str",
        "tools NumPy str alias",
        expected=12,
    )
    source = _replace_exact(
        source,
        "    print ' '\n"
        "    print '  Success: %s resampled by a factor %ix%i and saved as %s' % \\\n"
        "           (fn.split('/')[-1],sampling,sampling,out_fn.split('/')[-1])\n",
        "    print(' ')\n"
        "    print('  Success: %s resampled by a factor %ix%i and saved as %s' %\n"
        "           (fn.split('/')[-1],sampling,sampling,out_fn.split('/')[-1]))\n",
        "tools success prints",
    )
    source = _replace_exact(
        source,
        "                print ' '\n"
        "                print str(counter)+'/'+str(loop_size)+': Running '+fn_out\n"
        "                print ' '\n",
        "                print(' ')\n"
        "                print(str(counter)+'/'+str(loop_size)+': Running '+fn_out)\n"
        "                print(' ')\n",
        "tools progress prints",
    )
    path.write_text(source, encoding="utf-8")


def _patch_pyqz_core(package_directory: Path) -> None:
    path = package_directory / "pyqz.py"
    source = path.read_text(encoding="utf-8")
    source = _replace_exact(
        source,
        "import matplotlib._cntr as cntr # To construct the KDE contours with plotting anything\n",
        "import contourpy\n",
        "contour backend import",
    )
    source = _replace_exact(
        source,
        "import pyqz_metadata as pyqzm\n"
        "from pyqz_metadata import __version__\n"
        "import pyqz_tools as pyqzt\n",
        "from . import pyqz_metadata as pyqzm\n"
        "from .pyqz_metadata import __version__\n"
        "from . import pyqz_tools as pyqzt\n",
        "core relative imports",
    )
    source = _replace_regex(
        source,
        r"np\.int(?![A-Za-z0-9_])",
        "int",
        "core NumPy int alias",
        expected=4,
    )
    source = _replace_regex(
        source,
        r"np\.float(?![A-Za-z0-9_])",
        "float",
        "core NumPy float alias",
        expected=1,
    )
    source = _replace_regex(
        source,
        r"np\.str(?![A-Za-z0-9_])",
        "str",
        "core NumPy str alias",
        expected=3,
    )
    source = _replace_regex(
        source,
        r"np\.infty(?![A-Za-z0-9_])",
        "np.inf",
        "core NumPy infinity alias",
        expected=2,
    )
    source = _replace_exact(
        source,
        "def get_global_qz_singlespec((j, final_cols, data, data_cols, which_grids,\n"
        "                              ids, qzs, Pk, kappa, struct, sampling, error_pdf, \n"
        "                              srs, flag_level, KDE_method, KDE_qz_sampling,\n"
        "                              KDE_do_singles, KDE_pickle_loc, verbose)):\n",
        "def get_global_qz_singlespec(args):\n"
        "    (j, final_cols, data, data_cols, which_grids,\n"
        "     ids, qzs, Pk, kappa, struct, sampling, error_pdf,\n"
        "     srs, flag_level, KDE_method, KDE_qz_sampling,\n"
        "     KDE_do_singles, KDE_pickle_loc, verbose) = args\n",
        "tuple parameter unpacking",
    )
    source = _replace_exact(
        source,
        "    nlines = len(data_cols)/2\n",
        "    nlines = len(data_cols)//2\n",
        "integer line count",
    )
    source = _replace_exact(
        source,
        "            results = map(get_global_qz_singlespec, jobs)    \n",
        "            results = list(map(get_global_qz_singlespec, jobs))\n",
        "eager map results",
    )
    source = _replace_exact(
        source,
        "            xv,yv = np.meshgrid(QZs_grid[qzs[0]],QZs_grid[qzs[1]])\n"
        "            # Now the magic function\n"
        "            c = cntr.Cntr(xv,yv, gridZ) \n"
        "            \n"
        "            paths = c.trace(pyqzm.PDF_cont_level) \n"
        "            npaths = len(paths)/2\n"
        "            paths = [Path(path,codes=paths[npaths+p]) for \n"
        "                       (p,path) in enumerate(paths[:npaths])] # Only keep the path nodes.\n",
        "            contour_generator = contourpy.contour_generator(\n"
        "                x=QZs_grid[qzs[0]], y=QZs_grid[qzs[1]], z=gridZ, name='serial'\n"
        "            )\n"
        "            paths = [Path(vertices) for vertices in\n"
        "                     contour_generator.lines(pyqzm.PDF_cont_level)]\n"
        "            npaths = len(paths)\n",
        "contour extraction",
    )
    source = _replace_exact(
        source,
        "        check1 = np.abs(mean_qz-mean_vert[qzs.index(qz)])/\\\n"
        "                                std_qz<=flag_level\n"
        "        check2 = np.abs(mean_qz-mean_vert[qzs.index(qz)])/\\\n"
        "                                err_vert[qzs.index(qz)]<=flag_level\n",
        "        with np.errstate(divide='ignore', invalid='ignore'):\n"
        "            check1 = np.abs(mean_qz-mean_vert[qzs.index(qz)])/\\\n"
        "                                    std_qz<=flag_level\n"
        "            check2 = np.abs(mean_qz-mean_vert[qzs.index(qz)])/\\\n"
        "                                    err_vert[qzs.index(qz)]<=flag_level\n",
        "zero-dispersion warning guard",
    )
    source = _replace_regex(
        source,
        r"^        print ' '$",
        "        print(' ')",
        "blank prints",
        expected=2,
    )
    print_replacements = {
        "            print '--> Received '+str(npoints)+' spectra ...'\n":
            "            print('--> Received '+str(npoints)+' spectra ...')\n",
        "            print '--> Received 1 spectrum ...'\n":
            "            print('--> Received 1 spectrum ...')\n",
        "                print '--> Dealing with them one at a time ... be patient now !'\n":
            "                print('--> Dealing with them one at a time ... be patient now !')\n",
        "                print '    (no status update until I am done ...)' \n":
            "                print('    (no status update until I am done ...)')\n",
        "                print '--> Launching the multiple processes ... ' \n":
            "                print('--> Launching the multiple processes ... ')\n",
        "        print 'All done in',dt.now()-starttime\n":
            "        print('All done in', dt.now()-starttime)\n",
    }
    for index, (old, new) in enumerate(print_replacements.items(), start=1):
        source = _replace_exact(source, old, new, f"core print {index}")
    path.write_text(source, encoding="utf-8")


def _build_pyqz_runtime_overlay(source_directory: Path) -> tempfile.TemporaryDirectory[str]:
    runtime_directory = tempfile.TemporaryDirectory(prefix="pyqz084-runtime-")
    package_directory = Path(runtime_directory.name) / "pyqz"
    try:
        shutil.copytree(
            source_directory,
            package_directory,
            ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
        )
        _patch_pyqz_init(package_directory)
        _patch_pyqz_tools(package_directory)
        _patch_pyqz_core(package_directory)
        for filename in ("__init__.py", "pyqz_metadata.py", "pyqz_tools.py", "pyqz.py"):
            source_path = package_directory / filename
            compile(source_path.read_text(encoding="utf-8"), str(source_path), "exec")
    except Exception:
        runtime_directory.cleanup()
        raise
    return runtime_directory


def load_pyqz_compat() -> Any:
    """Load pyqz 0.8.4 from a patched, disposable runtime copy.

    The installed package and its MAPPINGS grids are copied but never edited.
    Exact source hashes guard every transformation so a different build fails
    before scientific calculations begin.
    """

    global _PYQZ_MODULE, _PYQZ_RUNTIME_DIRECTORY

    if _PYQZ_MODULE is not None:
        return _PYQZ_MODULE
    if "pyqz" in sys.modules:
        existing = sys.modules["pyqz"]
        raise RuntimeError(
            "pyqz was imported before the runtime compatibility adapter "
            f"({getattr(existing, '__file__', 'unknown path')}). Start a fresh "
            "process and call load_pyqz_compat() before importing pyqz directly."
        )

    source_directory = _validate_pyqz_installation()
    runtime_directory = _build_pyqz_runtime_overlay(source_directory)
    package_directory = Path(runtime_directory.name) / "pyqz"
    spec = importlib.util.spec_from_file_location(
        "pyqz",
        package_directory / "__init__.py",
        submodule_search_locations=[str(package_directory)],
    )
    if spec is None or spec.loader is None:
        runtime_directory.cleanup()
        raise RuntimeError("Could not create an import specification for pyqz overlay.")

    module = importlib.util.module_from_spec(spec)
    sys.modules["pyqz"] = module
    try:
        spec.loader.exec_module(module)
    except Exception as exc:
        for module_name in list(sys.modules):
            if module_name == "pyqz" or module_name.startswith("pyqz."):
                del sys.modules[module_name]
        runtime_directory.cleanup()
        raise RuntimeError("Failed to import the pyqz runtime overlay.") from exc

    _PYQZ_RUNTIME_DIRECTORY = runtime_directory
    _PYQZ_MODULE = module
    return module
