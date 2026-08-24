from __future__ import annotations

import importlib
from importlib import metadata
import hashlib
from pathlib import Path
import warnings

import numpy as np
import pytest


def import_compat_module():
    try:
        return importlib.import_module("model_grid_compat")
    except ModuleNotFoundError as exc:
        pytest.fail(f"model-grid compatibility module is missing: {exc}")


def test_nebulabayes_runtime_compatibility_reaches_a_finite_hii_fit():
    compat = import_compat_module()

    nebula_bayes = compat.load_nebulabayes_compat()
    grid, grid_params, grid_path = compat.load_nebulabayes_hii_grid()

    assert metadata.version("NebulaBayes") == "0.9.9"
    assert nebula_bayes.__name__ == "NebulaBayes"
    assert np.trapz is np.trapezoid
    assert Path(grid_path).name == "NB_HII_grid.fits.gz"
    assert grid_params == ["log U", "log P/k", "12 + log O/H"]
    assert all(grid[param].dtype == np.dtype("float64") for param in grid_params)

    line_list = [
        "Hbeta",
        "OIII5007",
        "Halpha",
        "NII6583",
        "SII6716",
        "SII6731",
    ]
    model = compat.make_nebulabayes_hii_model(
        line_list,
        interpd_grid_shape=(40, 20, 160),
    )
    result = model(
        [1.0, 0.55, 2.86, 0.45, 0.35, 0.30],
        [0.05, 0.04, 0.10, 0.035, 0.03, 0.03],
        line_list,
        norm_line="Hbeta",
        deredden=False,
        prior="Uniform",
        likelihood_lines=None,
    )

    estimates = result.Posterior.DF_estimates
    assert model.Interpd_grids.shape == (40, 20, 160)
    assert estimates.loc["12 + log O/H", "Estimate"] == pytest.approx(
        8.274195,
        abs=1.0e-6,
    )
    assert estimates.loc["log U", "Estimate"] == pytest.approx(-3.384615, abs=1.0e-6)
    assert estimates.loc["log P/k", "Estimate"] == pytest.approx(6.747369, abs=1.0e-6)
    assert np.isfinite(result.Posterior.best_model["chi2"])


def test_pyqz_runtime_overlay_is_deterministic_and_keeps_site_packages_pristine():
    compat = import_compat_module()
    distribution = metadata.distribution("pyqz")
    installed_root = Path(distribution.locate_file("pyqz")).resolve()
    installed_hashes_before = {
        name: hashlib.sha256((installed_root / name).read_bytes()).hexdigest()
        for name in compat.PYQZ_SOURCE_HASHES
    }

    pyqz = compat.load_pyqz_compat()

    installed_hashes_after = {
        name: hashlib.sha256((installed_root / name).read_bytes()).hexdigest()
        for name in compat.PYQZ_SOURCE_HASHES
    }
    assert metadata.version("pyqz") == "0.8.4.0"
    assert pyqz.__version__ == "0.8.4"
    assert installed_hashes_before == compat.PYQZ_SOURCE_HASHES
    assert installed_hashes_after == installed_hashes_before
    assert not Path(pyqz.__file__).resolve().is_relative_to(installed_root)

    grid, columns, grid_metadata = pyqz.get_grid(
        "[NII]/[SII]+;[OIII]/[SII]+",
        Pk=5.0,
        kappa=np.inf,
        struct="pp",
        sampling=2,
    )
    assert grid.shape == (153, 7)
    assert columns == [
        "LogQ",
        "Tot[O]+12",
        "gas[O]+12",
        "[NII]/[SII]+",
        "[OIII]/[SII]+",
        "Mix_x",
        "Mix_y",
    ]
    assert "MAPPINGS V HII Region Grid" in grid_metadata["MV_id"]

    data = np.array(
        [[0.567, 0.0284, 0.511, 0.0255, 2.38, 0.119]],
        dtype=float,
    )
    data_cols = [
        "[NII]",
        "std[NII]",
        "[SII]+",
        "std[SII]+",
        "[OIII]",
        "std[OIII]",
    ]

    def run_once(srs=80):
        np.random.seed(20260823)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            return pyqz.get_global_qz(
                data,
                data_cols,
                ["[NII]/[SII]+;[OIII]/[SII]+"],
                qzs=["LogQ", "gas[O]+12"],
                Pk=5.0,
                kappa=np.inf,
                struct="pp",
                sampling=2,
                error_pdf="normal",
                srs=srs,
                KDE_method="multiv",
                KDE_qz_sampling=101j,
                KDE_do_singles=True,
                verbose=False,
                nproc=1,
            )

    values_1, output_cols_1 = run_once()
    values_2, output_cols_2 = run_once()
    np.testing.assert_array_equal(values_1, values_2)
    assert output_cols_1 == output_cols_2

    output = dict(zip(output_cols_1, values_1[0]))
    assert output["<LogQ{KDE}>"] == pytest.approx(7.370267749, abs=1.0e-8)
    assert output["err(LogQ{KDE})"] == pytest.approx(0.043950184, abs=1.0e-8)
    assert output["<gas[O]+12{KDE}>"] == pytest.approx(8.374295773, abs=1.0e-8)
    assert output["err(gas[O]+12{KDE})"] == pytest.approx(0.049729275, abs=1.0e-8)
    assert output["flag"] == 13
    assert output["rs_offgrid"] == 0

    production_1, production_columns_1 = run_once(srs=800)
    production_2, production_columns_2 = run_once(srs=800)
    np.testing.assert_array_equal(production_1, production_2)
    assert production_columns_1 == production_columns_2
    production = dict(zip(production_columns_1, production_1[0]))
    assert np.isfinite(production["<LogQ{KDE}>"])
    assert np.isfinite(production["err(LogQ{KDE})"])
    assert np.isfinite(production["<gas[O]+12{KDE}>"])
    assert np.isfinite(production["err(gas[O]+12{KDE})"])


def test_pyqz_upstream_grid_nodes_and_archived_example_remain_compatible():
    pyqz = import_compat_module().load_pyqz_compat()
    diagnostic = "[NII]/[SII]+;[OIII]/[SII]+"
    coefficients = pyqz.pyqzm.diagnostics[diagnostic]["coeffs"]
    grid, columns, _ = pyqz.get_grid(
        diagnostic,
        Pk=5.0,
        kappa=np.inf,
        struct="pp",
        coeffs=coefficients,
        sampling=1,
    )
    grid_nodes = grid[
        :,
        [
            columns.index("LogQ"),
            columns.index("Tot[O]+12"),
            columns.index("Mix_x"),
            columns.index("Mix_y"),
        ],
    ]
    interp_q = pyqz.interp_qz(
        "LogQ",
        [grid_nodes[:, -2], grid_nodes[:, -1]],
        diagnostic,
        coeffs=coefficients,
        Pk=5.0,
        kappa=np.inf,
        struct="pp",
        sampling=1,
    )
    interp_z = pyqz.interp_qz(
        "Tot[O]+12",
        [grid_nodes[:, -2], grid_nodes[:, -1]],
        diagnostic,
        coeffs=coefficients,
        Pk=5.0,
        kappa=np.inf,
        struct="pp",
        sampling=1,
    )
    assert np.all(np.round(interp_q, 2) == grid_nodes[:, 0])
    assert np.all(np.round(interp_z, 3) == grid_nodes[:, 1])

    off_grid = pyqz.interp_qz(
        "LogQ",
        [np.array(-1000.0), np.array(-1000.0)],
        diagnostic,
        coeffs=coefficients,
        Pk=5.0,
        kappa=np.inf,
        struct="pp",
        sampling=1,
    )
    assert np.isnan(off_grid)

    np.random.seed(20260823)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        bad_values, bad_columns = pyqz.get_global_qz(
            np.array([[0.567, 0.0284, 0.511, 0.0255, 2380.0, 119.0]]),
            ["[NII]", "std[NII]", "[SII]+", "std[SII]+", "[OIII]", "std[OIII]"],
            [diagnostic],
            qzs=["LogQ", "gas[O]+12"],
            Pk=5.0,
            kappa=np.inf,
            struct="pp",
            sampling=2,
            error_pdf="normal",
            srs=20,
            KDE_method="multiv",
            KDE_qz_sampling=51j,
            KDE_do_singles=True,
            verbose=False,
            nproc=1,
        )
    bad_output = dict(zip(bad_columns, bad_values[0]))
    assert bad_output["flag"] == 8
    assert np.isnan(bad_output["<LogQ>"])
    assert np.isnan(bad_output["<gas[O]+12>"])

    example_data = np.array(
        [[1.0, 0.05, 2.38, 0.119, 5.07, 0.2534, 0.567, 0.0284,
          0.511, 0.0255, 2.88, 0.144]],
        dtype=float,
    )
    example_columns = [
        "Hb", "stdHb", "[OIII]", "std[OIII]", "[OII]+", "std[OII]+",
        "[NII]", "std[NII]", "[SII]+", "std[SII]+", "Ha", "stdHa",
    ]
    example_diagnostics = [
        "[NII]/[SII]+;[OIII]/Hb",
        "[NII]/[OII]+;[OIII]/[SII]+",
    ]
    np.random.seed(20260823)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        values, output_columns = pyqz.get_global_qz(
            example_data,
            example_columns,
            example_diagnostics,
            qzs=["LogQ", "Tot[O]+12"],
            Pk=5.0,
            kappa=np.inf,
            struct="pp",
            sampling=1,
            error_pdf="normal",
            srs=400,
            KDE_method="multiv",
            KDE_qz_sampling=101j,
            KDE_do_singles=True,
            verbose=False,
            nproc=1,
        )
    output = dict(zip(output_columns, values[0]))

    archived_direct = {
        "[NII]/[SII]+;[OIII]/Hb|LogQ": 7.38832,
        "[NII]/[SII]+;[OIII]/Hb|Tot[O]+12": 8.48760,
        "[NII]/[OII]+;[OIII]/[SII]+|LogQ": 7.37681,
        "[NII]/[OII]+;[OIII]/[SII]+|Tot[O]+12": 8.48372,
        "<LogQ>": 7.38257,
        "<Tot[O]+12>": 8.48566,
    }
    for name, expected in archived_direct.items():
        assert output[name] == pytest.approx(expected, abs=5.0e-6)

    archived_kde = {
        "[NII]/[SII]+;[OIII]/Hb|LogQ{KDE}": 7.39892,
        "err([NII]/[SII]+;[OIII]/Hb|LogQ{KDE})": 0.04236,
        "[NII]/[SII]+;[OIII]/Hb|Tot[O]+12{KDE}": 8.48280,
        "err([NII]/[SII]+;[OIII]/Hb|Tot[O]+12{KDE})": 0.05038,
        "[NII]/[OII]+;[OIII]/[SII]+|LogQ{KDE}": 7.38182,
        "err([NII]/[OII]+;[OIII]/[SII]+|LogQ{KDE})": 0.02943,
        "[NII]/[OII]+;[OIII]/[SII]+|Tot[O]+12{KDE}": 8.48495,
        "err([NII]/[OII]+;[OIII]/[SII]+|Tot[O]+12{KDE})": 0.02561,
        "<LogQ{KDE}>": 7.38541,
        "err(LogQ{KDE})": 0.03434,
        "<Tot[O]+12{KDE}>": 8.48543,
        "err(Tot[O]+12{KDE})": 0.03251,
    }
    for name, expected in archived_kde.items():
        assert output[name] == pytest.approx(expected, abs=0.02)
    assert output["flag"] == 0
    assert output["rs_offgrid"] == 0
