from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest


JY22_GRID_PATH = (
    Path(__file__).resolve().parents[1]
    / "Peng2026"
    / "photoionization_models"
    / "photoionization_grid_interpolated.fits"
)


def _write_synthetic_jy22_grid(path, rows, *, omit_column=None):
    from astropy.io import fits

    names = (
        "grid_kind",
        "log_Z_Zsun",
        "log_U",
        "log_R3_OIII5007_Hbeta",
        "log_N2_NII6584_Halpha",
        "log_S2_SII6717_6731_Halpha",
    )
    columns = []
    for index, name in enumerate(names):
        if name == omit_column:
            continue
        values = np.asarray([row[index] for row in rows])
        if name == "grid_kind":
            columns.append(fits.Column(name=name, format="12A", array=values))
        else:
            columns.append(
                fits.Column(name=name, format="D", array=values.astype(float))
            )
    fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU.from_columns(columns)]).writeto(
        path
    )


def _make_linear_jy22_grid():
    from model_grid_diagnostics import JY22Grid

    o_h = np.array([8.0, 8.5, 9.0])
    log_u = np.array([-4.0, -3.0, -2.0])
    ratios = np.empty((o_h.size, log_u.size, 3), dtype=float)
    for z_index, abundance in enumerate(o_h):
        for u_index, ionisation in enumerate(log_u):
            ratios[z_index, u_index] = [
                -1.0 + 0.8 * (abundance - 8.5) + 0.05 * (ionisation + 3.0),
                -0.5 + 0.2 * (abundance - 8.5) + 0.25 * (ionisation + 3.0),
                0.2 - 0.4 * (abundance - 8.5) + 0.35 * (ionisation + 3.0),
            ]
    return JY22Grid(
        source_path=Path("synthetic.fits"),
        sha256="synthetic",
        log_z_zsun=o_h - 8.69,
        o_h=o_h,
        log_u=log_u,
        model_ratios=ratios,
        requested_log_u_bounds=(-4.0, -0.5),
        effective_log_u_bounds=(-4.0, -2.0),
        solar_o_h=8.69,
    )


def test_common_valid_mask_requires_positive_finite_six_line_fluxes_and_errors():
    from model_grid_diagnostics import MODEL_LINE_BASES, build_model_input_valid_mask

    shape = (2, 3)
    flux_maps = {name: np.full(shape, 1.0) for name in MODEL_LINE_BASES}
    error_maps = {name: np.full(shape, 10.0) for name in MODEL_LINE_BASES}
    region = np.array([[True, True, True], [True, True, False]])

    flux_maps["OIII5006"][0, 1] = 0.0
    error_maps["NII6583"][1, 0] = np.nan
    flux_maps["SII6730"][1, 1] = 0.01  # Low S/N remains valid by design.

    valid = build_model_input_valid_mask(flux_maps, error_maps, region)

    expected = np.array([[True, False, True], [False, True, False]])
    np.testing.assert_array_equal(valid, expected)


def test_sii_sum_uses_quadrature_error():
    from model_grid_diagnostics import sum_sii_doublet

    flux, error = sum_sii_doublet(3.0, 4.0, 0.3, 0.4)

    assert flux == pytest.approx(7.0)
    assert error == pytest.approx(0.5)


def test_nebulabayes_hbeta_normalisation_propagates_denominator_error():
    from model_grid_diagnostics import (
        NB_LINE_LIST,
        normalise_nebulabayes_to_hbeta,
    )

    fluxes = np.array([10.0, 5.0, 28.6, 4.0, 3.0, 2.0])
    errors = np.array([1.0, 0.5, 2.0, 0.4, 0.3, 0.2])

    ratios, ratio_errors = normalise_nebulabayes_to_hbeta(
        fluxes, errors, NB_LINE_LIST
    )

    expected_ratios = fluxes / fluxes[0]
    expected_errors = expected_ratios * np.sqrt(
        (errors / fluxes) ** 2 + (errors[0] / fluxes[0]) ** 2
    )
    np.testing.assert_allclose(ratios, expected_ratios)
    np.testing.assert_allclose(ratio_errors, expected_errors)
    assert ratios[0] == pytest.approx(1.0)


def test_logq_to_logu_uses_speed_of_light_in_cm_per_second():
    from model_grid_diagnostics import LOG10_C_CM_S, logq_to_logu

    assert LOG10_C_CM_S == pytest.approx(10.47712125472, abs=1.0e-11)
    assert logq_to_logu(7.5) == pytest.approx(-2.97712125472, abs=1.0e-11)


def test_jy22_empty_result_has_complete_typed_contract():
    from model_grid_diagnostics import (
        JY22_FLAG_COVARIANCE_INVALID,
        JY22_FLOAT_FIELDS,
        JY22_INTEGER_FIELDS,
        _records_to_arrays,
        empty_jy22_result,
    )

    result = empty_jy22_result(covariance_invalid=True)

    assert set(result) == set(JY22_FLOAT_FIELDS) | set(JY22_INTEGER_FIELDS)
    assert all(np.isnan(result[name]) for name in JY22_FLOAT_FIELDS)
    assert result["flag"] == JY22_FLAG_COVARIANCE_INVALID
    assert result["valid"] == 0

    arrays = _records_to_arrays(
        [result], JY22_FLOAT_FIELDS, JY22_INTEGER_FIELDS
    )
    assert all(arrays[name].dtype == np.float64 for name in JY22_FLOAT_FIELDS)
    assert all(arrays[name].dtype == np.int16 for name in JY22_INTEGER_FIELDS)


def test_jy22_loader_uses_released_grid_and_documented_logu_range():
    from model_grid_diagnostics import load_jy22_grid

    grid = load_jy22_grid(JY22_GRID_PATH)

    assert grid.sha256 == (
        "d7a219b60c9a1ea8b29339988c84b1832028b28b3c417d1c3c58b420831eb38a"
    )
    assert grid.log_z_zsun.shape == (40,)
    assert grid.o_h.shape == (40,)
    assert grid.log_u.shape == (28,)
    assert grid.model_ratios.shape == (40, 28, 3)
    assert grid.requested_log_u_bounds == (-4.0, -0.5)
    assert grid.effective_log_u_bounds == pytest.approx(
        (-4.0, -0.5384615384615388)
    )
    np.testing.assert_allclose(grid.o_h, grid.log_z_zsun + 8.69)
    assert grid.solar_o_h == pytest.approx(8.69)
    assert grid.source_path == JY22_GRID_PATH.resolve()


def test_jy22_resampling_matches_paper_grid_and_preserves_linear_surface():
    from model_grid_diagnostics import resample_jy22_grid

    source = _make_linear_jy22_grid()
    resampled = resample_jy22_grid(source, n_o_h=200, n_log_u=200)

    assert resampled.o_h.shape == (200,)
    assert resampled.log_z_zsun.shape == (200,)
    assert resampled.log_u.shape == (200,)
    assert resampled.model_ratios.shape == (200, 200, 3)
    assert resampled.o_h[[0, -1]] == pytest.approx(source.o_h[[0, -1]])
    assert resampled.log_u[[0, -1]] == pytest.approx(source.log_u[[0, -1]])
    assert resampled.requested_log_u_bounds == source.requested_log_u_bounds
    assert resampled.effective_log_u_bounds == source.effective_log_u_bounds
    assert resampled.source_path == source.source_path
    assert resampled.sha256 == source.sha256

    z_index, u_index = 73, 121
    abundance = resampled.o_h[z_index]
    ionisation = resampled.log_u[u_index]
    expected = np.array(
        [
            -1.0 + 0.8 * (abundance - 8.5) + 0.05 * (ionisation + 3.0),
            -0.5 + 0.2 * (abundance - 8.5) + 0.25 * (ionisation + 3.0),
            0.2 - 0.4 * (abundance - 8.5) + 0.35 * (ionisation + 3.0),
        ]
    )
    np.testing.assert_allclose(
        resampled.model_ratios[z_index, u_index], expected, atol=1.0e-12
    )


def test_jy22_resampling_uses_explicit_inference_logu_bounds():
    from model_grid_diagnostics import resample_jy22_grid

    source = _make_linear_jy22_grid()
    resampled = resample_jy22_grid(
        source,
        n_o_h=20,
        n_log_u=30,
        min_log_u=-3.8,
        max_log_u=-2.2,
    )

    assert resampled.log_u[[0, -1]] == pytest.approx([-3.8, -2.2])
    assert resampled.requested_log_u_bounds == (-3.8, -2.2)
    assert resampled.effective_log_u_bounds == (-3.8, -2.2)


@pytest.mark.parametrize(
    ("rows", "match"),
    [
        (
            [
                ("interpolated", -0.5, -4.0, 0.1, -1.0, -0.5),
                ("interpolated", -0.5, -3.0, 0.2, -0.9, -0.4),
                ("interpolated", 0.0, -4.0, 0.3, -0.8, -0.3),
            ],
            "rectangular",
        ),
        (
            [
                ("interpolated", -0.5, -4.0, 0.1, -1.0, -0.5),
                ("interpolated", -0.5, -4.0, 0.2, -0.9, -0.4),
                ("interpolated", 0.0, -4.0, 0.3, -0.8, -0.3),
                ("interpolated", 0.0, -3.0, 0.4, -0.7, -0.2),
            ],
            "duplicate",
        ),
        (
            [
                ("interpolated", -0.5, 0.0, 0.1, -1.0, -0.5),
                ("interpolated", 0.0, 0.5, 0.2, -0.9, -0.4),
            ],
            "requested log-U range",
        ),
    ],
)
def test_jy22_loader_rejects_malformed_or_out_of_range_grids(
    tmp_path, rows, match
):
    from model_grid_diagnostics import load_jy22_grid

    path = tmp_path / "bad_grid.fits"
    _write_synthetic_jy22_grid(path, rows)

    with pytest.raises(ValueError, match=match):
        load_jy22_grid(path)


def test_jy22_loader_rejects_missing_required_columns(tmp_path):
    from model_grid_diagnostics import load_jy22_grid

    rows = [
        ("interpolated", -0.5, -4.0, 0.1, -1.0, -0.5),
        ("interpolated", 0.0, -4.0, 0.2, -0.9, -0.4),
    ]
    path = tmp_path / "missing_column.fits"
    _write_synthetic_jy22_grid(path, rows, omit_column="log_U")

    with pytest.raises(KeyError, match="log_U"):
        load_jy22_grid(path)


def test_jy22_loader_rejects_mixed_grid_kind_in_requested_range(tmp_path):
    from model_grid_diagnostics import load_jy22_grid

    rows = [
        ("interpolated", -0.5, -4.0, 0.10, -1.00, -0.50),
        ("interpolated", -0.5, -3.0, 0.20, -0.90, -0.40),
        ("interpolated", 0.0, -4.0, 0.30, -0.80, -0.30),
        ("interpolated", 0.0, -3.0, 0.40, -0.70, -0.20),
        ("original", -0.5, -4.0, 9.99, 9.99, 9.99),
    ]
    path = tmp_path / "mixed_grid_kind.fits"
    _write_synthetic_jy22_grid(path, rows)

    with pytest.raises(ValueError, match="non-interpolated"):
        load_jy22_grid(path)


def test_jy22_loader_reorders_released_ratio_columns_to_n2_s2_r3(tmp_path):
    from model_grid_diagnostics import load_jy22_grid

    rows = [
        ("interpolated", -0.5, -4.0, 0.10, -1.00, -0.50),
        ("interpolated", -0.5, -3.0, 0.20, -0.90, -0.40),
        ("interpolated", 0.0, -4.0, 0.30, -0.80, -0.30),
        ("interpolated", 0.0, -3.0, 0.40, -0.70, -0.20),
    ]
    path = tmp_path / "ratio_order.fits"
    _write_synthetic_jy22_grid(path, rows)

    grid = load_jy22_grid(path)

    np.testing.assert_allclose(grid.model_ratios[0, 0], [-1.0, -0.5, 0.1])


def test_jy22_raw_ratio_jacobian_matches_finite_differences():
    from model_grid_diagnostics import jy22_log_ratios_and_covariance

    fluxes = np.array([10.0, 7.0, 35.0, 5.0, 3.0, 2.0])
    errors = np.array([0.8, 0.5, 1.2, 0.4, 0.3, 0.25])

    measurement = jy22_log_ratios_and_covariance(fluxes, errors)
    numerical = np.empty_like(measurement.jacobian)
    for index, value in enumerate(fluxes):
        step = value * 1.0e-6
        upper = fluxes.copy()
        lower = fluxes.copy()
        upper[index] += step
        lower[index] -= step
        upper_ratios = jy22_log_ratios_and_covariance(
            upper, errors
        ).log_ratios
        lower_ratios = jy22_log_ratios_and_covariance(
            lower, errors
        ).log_ratios
        numerical[:, index] = (upper_ratios - lower_ratios) / (2.0 * step)

    np.testing.assert_allclose(measurement.jacobian, numerical, rtol=2e-7, atol=2e-9)
    np.testing.assert_allclose(
        measurement.covariance,
        measurement.jacobian @ np.diag(errors**2) @ measurement.jacobian.T,
        rtol=1e-12,
        atol=1e-14,
    )
    np.testing.assert_allclose(
        measurement.covariance, measurement.covariance.T, atol=1e-15
    )
    assert np.all(np.linalg.eigvalsh(measurement.covariance) > 0.0)


def test_jy22_raw_ratios_keep_only_shared_halpha_covariance():
    from model_grid_diagnostics import jy22_log_ratios_and_covariance

    fluxes = np.array([10.0, 6.0, 25.0, 5.0, 3.0, 2.0])
    errors = np.array([1.0, 0.6, 2.0, 0.5, 0.3, 0.2])

    measurement = jy22_log_ratios_and_covariance(fluxes, errors)

    np.testing.assert_allclose(measurement.ratio_fluxes, fluxes)
    ln10 = np.log(10.0)
    expected_n2_s2_covariance = (errors[2] / fluxes[2] / ln10) ** 2
    assert measurement.covariance[0, 1] == pytest.approx(
        expected_n2_s2_covariance
    )
    assert measurement.covariance[0, 2] == pytest.approx(0.0, abs=1e-15)
    assert measurement.covariance[1, 2] == pytest.approx(0.0, abs=1e-15)


def test_jy22_raw_ratio_covariance_is_scale_invariant():
    from model_grid_diagnostics import jy22_log_ratios_and_covariance

    fluxes = np.array([10.0, 7.0, 35.0, 5.0, 3.0, 2.0])
    errors = 0.08 * fluxes

    original = jy22_log_ratios_and_covariance(fluxes, errors)
    scaled = jy22_log_ratios_and_covariance(
        fluxes * 1.0e5, errors * 1.0e5
    )

    np.testing.assert_allclose(original.log_ratios, scaled.log_ratios, atol=1e-14)
    np.testing.assert_allclose(original.covariance, scaled.covariance, atol=1e-14)
    expected_s2 = np.log10((fluxes[4] + fluxes[5]) / fluxes[2])
    assert original.log_ratios[1] == pytest.approx(expected_s2)
    np.testing.assert_allclose(original.ratio_fluxes, fluxes)


@pytest.mark.parametrize(
    ("fluxes", "errors"),
    [
        ([10, 5, 28.6, 4, 3, 0], [1, 1, 1, 1, 1, 1]),
        ([10, 5, 28.6, 4, 3, 2], [1, 1, -1, 1, 1, 1]),
    ],
)
def test_jy22_ratio_covariance_rejects_invalid_measurements(
    fluxes, errors
):
    from model_grid_diagnostics import jy22_log_ratios_and_covariance

    with pytest.raises(ValueError):
        jy22_log_ratios_and_covariance(fluxes, errors)


def test_jy22_posterior_recovers_grid_node_and_equal_tailed_marginals():
    from model_grid_diagnostics import JY22_DEPLETION_OFFSET_DEX, infer_jy22_posterior

    grid = _make_linear_jy22_grid()
    observed = grid.model_ratios[1, 1]

    result = infer_jy22_posterior(observed, np.eye(3) * 1.0e-6, grid)

    assert JY22_DEPLETION_OFFSET_DEX == pytest.approx(0.22)
    assert result["o_h_pre"] == pytest.approx(8.5, abs=1e-12)
    assert result["o_h_pre_16"] == pytest.approx(8.5, abs=1e-12)
    assert result["o_h_pre_84"] == pytest.approx(8.5, abs=1e-12)
    assert result["o_h"] == pytest.approx(8.28, abs=1e-12)
    assert result["o_h_16"] == pytest.approx(8.28, abs=1e-12)
    assert result["o_h_84"] == pytest.approx(8.28, abs=1e-12)
    assert result["log_u"] == pytest.approx(-3.0, abs=1e-12)
    assert result["log_u_16"] == pytest.approx(-3.0, abs=1e-12)
    assert result["log_u_84"] == pytest.approx(-3.0, abs=1e-12)
    assert result["chi2_min"] == pytest.approx(0.0, abs=1e-14)
    assert result["resid_n2"] == pytest.approx(0.0, abs=1e-14)
    assert result["resid_s2"] == pytest.approx(0.0, abs=1e-14)
    assert result["resid_r3"] == pytest.approx(0.0, abs=1e-14)
    assert result["resid_norm"] == pytest.approx(0.0, abs=1e-14)
    assert result["flag"] == 0
    assert result["valid"] == 1
    assert result["fit_ok"] == 1


def test_jy22_quantiles_interpolate_through_ordinary_cumulative_marginal():
    from model_grid_diagnostics import _weighted_axis_quantile

    axis = np.array([0.0, 1.0, 2.0])
    marginal = np.array([0.2, 0.3, 0.5])

    assert _weighted_axis_quantile(axis, marginal, 0.16) == pytest.approx(0.0)
    assert _weighted_axis_quantile(axis, marginal, 0.84) == pytest.approx(1.68)


def test_jy22_posterior_uses_full_covariance_and_broadens_with_error():
    from model_grid_diagnostics import JY22_DEPLETION_OFFSET_DEX, infer_jy22_posterior

    grid = _make_linear_jy22_grid()
    observed = grid.model_ratios[1, 1] + np.array([0.10, -0.08, 0.06])
    covariance = np.array(
        [[0.04, 0.025, -0.010], [0.025, 0.05, 0.012], [-0.010, 0.012, 0.03]]
    )

    result = infer_jy22_posterior(observed, covariance, grid)
    residual = grid.model_ratios - observed
    manual_chi2 = np.einsum(
        "...i,ij,...j->...", residual, np.linalg.inv(covariance), residual
    )
    weights = np.exp(-0.5 * (manual_chi2 - np.min(manual_chi2)))
    weights /= np.sum(weights)
    expected_o_h_pre = np.sum(weights * grid.o_h[:, None])
    expected_log_u = np.sum(weights * grid.log_u[None, :])
    best_index = np.unravel_index(np.argmin(manual_chi2), manual_chi2.shape)
    expected_residual = observed - grid.model_ratios[best_index]

    assert result["o_h_pre"] == pytest.approx(expected_o_h_pre, abs=1e-12)
    assert result["o_h"] == pytest.approx(
        expected_o_h_pre - JY22_DEPLETION_OFFSET_DEX, abs=1e-12
    )
    assert result["log_u"] == pytest.approx(expected_log_u, abs=1e-12)
    assert result["chi2_min"] == pytest.approx(np.min(manual_chi2), abs=1e-12)
    np.testing.assert_allclose(
        [result["resid_n2"], result["resid_s2"], result["resid_r3"]],
        expected_residual,
        atol=1e-12,
    )
    assert result["resid_norm"] == pytest.approx(
        np.linalg.norm(expected_residual), abs=1e-12
    )

    narrow = infer_jy22_posterior(
        grid.model_ratios[1, 1], np.eye(3) * 0.0025, grid
    )
    broad = infer_jy22_posterior(
        grid.model_ratios[1, 1], np.eye(3) * 0.25, grid
    )
    assert broad["o_h_84"] - broad["o_h_16"] > (
        narrow["o_h_84"] - narrow["o_h_16"]
    )
    assert broad["log_u_84"] - broad["log_u_16"] > (
        narrow["log_u_84"] - narrow["log_u_16"]
    )


def test_jy22_posterior_flags_edges_and_invalid_covariance():
    from model_grid_diagnostics import (
        JY22_FLAG_COVARIANCE_INVALID,
        JY22_FLAG_LOGU_EDGE,
        JY22_FLAG_OH_EDGE,
        infer_jy22_posterior,
    )

    grid = _make_linear_jy22_grid()
    edge = infer_jy22_posterior(
        grid.model_ratios[0, 0], np.eye(3) * 1.0e-6, grid
    )
    invalid = infer_jy22_posterior(
        grid.model_ratios[1, 1], np.ones((3, 3)), grid
    )

    assert edge["flag"] & JY22_FLAG_OH_EDGE
    assert edge["flag"] & JY22_FLAG_LOGU_EDGE
    assert edge["valid"] == 1
    assert edge["fit_ok"] == 1
    assert invalid["flag"] == JY22_FLAG_COVARIANCE_INVALID
    assert invalid["valid"] == 0
    assert invalid["fit_ok"] == 0


def test_jy22_posterior_marks_empirical_gross_mismatch_without_invalidating_result():
    from model_grid_diagnostics import (
        JY22_EMPIRICAL_CHI2_MAX,
        JY22_FLAG_POOR_FIT,
        JY22_FORMAL_CHI2_MAX,
        infer_jy22_posterior,
    )

    grid = _make_linear_jy22_grid()
    observed = grid.model_ratios[1, 1] + np.array([0.3, 0.3, 0.3])

    result = infer_jy22_posterior(observed, np.eye(3) * 1.0e-4, grid)

    assert JY22_EMPIRICAL_CHI2_MAX == pytest.approx(25.0)
    assert JY22_FORMAL_CHI2_MAX == JY22_EMPIRICAL_CHI2_MAX
    assert result["chi2_min"] > JY22_EMPIRICAL_CHI2_MAX
    assert result["flag"] & JY22_FLAG_POOR_FIT
    assert result["valid"] == 1
    assert result["fit_ok"] == 0


def test_jy22_batch_keeps_order_continues_after_bad_covariance_and_is_deterministic():
    from model_grid_diagnostics import (
        JY22RatioSpectra,
        JY22_FLAG_COVARIANCE_INVALID,
        run_jy22_spectra,
    )

    grid = _make_linear_jy22_grid()
    spectra = JY22RatioSpectra(
        bin_ids=np.array([3, 7, 11]),
        log_ratios=np.array(
            [grid.model_ratios[1, 1], grid.model_ratios[1, 1], grid.model_ratios[2, 1]]
        ),
        covariances=np.array(
            [np.eye(3) * 0.01, np.ones((3, 3)), np.eye(3) * 0.01]
        ),
        pixel_counts=np.array([2, 3, 1]),
    )

    first = run_jy22_spectra(grid, spectra)
    second = run_jy22_spectra(grid, spectra)

    assert first.results["valid"].tolist() == [1, 0, 1]
    assert first.results["flag"][1] == JY22_FLAG_COVARIANCE_INVALID
    assert first.failures == ()
    for name, values in first.results.items():
        np.testing.assert_array_equal(values, second.results[name])


def test_jy22_preprocessing_failure_is_classified_as_invalid_covariance():
    from model_grid_diagnostics import (
        BinSpectra,
        JY22_FLAG_COVARIANCE_INVALID,
        build_jy22_ratio_spectra,
        run_jy22_spectra,
    )

    raw = BinSpectra(
        bin_ids=np.array([17]),
        fluxes=np.array([[10.0, 7.0, 35.0, 5.0, 3.0, 0.0]]),
        errors=np.ones((1, 6)),
        pixel_counts=np.array([1]),
    )
    ratios = build_jy22_ratio_spectra(raw)

    result = run_jy22_spectra(_make_linear_jy22_grid(), ratios)

    assert result.results["flag"][0] == JY22_FLAG_COVARIANCE_INVALID
    assert result.results["valid"][0] == 0


def test_unique_bin_spectra_are_sorted_and_broadcast_back_to_pixels():
    from model_grid_diagnostics import (
        MODEL_LINE_BASES,
        broadcast_bin_results,
        extract_unique_bin_spectra,
    )

    binid = np.array([[2, 2, 0], [1, 1, -1]], dtype=np.int32)
    valid = binid >= 0
    base = np.array([[20.0, 20.0, 5.0], [10.0, 10.0, np.nan]])
    flux_maps = {
        name: base + index for index, name in enumerate(MODEL_LINE_BASES)
    }
    error_maps = {
        name: np.where(valid, 0.1 + index, np.nan)
        for index, name in enumerate(MODEL_LINE_BASES)
    }

    spectra = extract_unique_bin_spectra(binid, valid, flux_maps, error_maps)

    np.testing.assert_array_equal(spectra.bin_ids, [0, 1, 2])
    np.testing.assert_array_equal(spectra.pixel_counts, [1, 2, 2])
    assert spectra.fluxes.shape == (3, len(MODEL_LINE_BASES))
    assert spectra.fluxes[0, 0] == pytest.approx(5.0)
    assert spectra.fluxes[1, 0] == pytest.approx(10.0)
    assert spectra.fluxes[2, 0] == pytest.approx(20.0)

    maps = broadcast_bin_results(
        binid,
        spectra.bin_ids,
        {
            "value": np.array([100.0, 200.0, 300.0]),
            "flag": np.array([0, 8, 9], dtype=np.int16),
        },
        integer_fields={"flag"},
    )
    np.testing.assert_allclose(
        maps["value"],
        np.array([[300.0, 300.0, 100.0], [200.0, 200.0, np.nan]]),
        equal_nan=True,
    )
    np.testing.assert_array_equal(
        maps["flag"],
        np.array([[9, 9, 0], [8, 8, -99]], dtype=np.int16),
    )


def test_unique_bin_spectra_reject_inconsistent_values_within_one_bin():
    from model_grid_diagnostics import MODEL_LINE_BASES, extract_unique_bin_spectra

    binid = np.array([[0, 0]], dtype=np.int32)
    valid = np.ones_like(binid, dtype=bool)
    flux_maps = {name: np.ones((1, 2)) for name in MODEL_LINE_BASES}
    error_maps = {name: np.ones((1, 2)) for name in MODEL_LINE_BASES}
    flux_maps["NII6583"][0, 1] = 1.2

    with pytest.raises(ValueError, match="NII6583.*BINID 0"):
        extract_unique_bin_spectra(binid, valid, flux_maps, error_maps)


def test_pyqz_extractor_uses_combined_kde_columns_and_raw_qc():
    from model_grid_diagnostics import extract_pyqz_result

    columns = [
        "<LogQ{KDE}>",
        "err(LogQ{KDE})",
        "<gas[O]+12{KDE}>",
        "err(gas[O]+12{KDE})",
        "flag",
        "rs_offgrid",
    ]
    values = np.array([7.4, 0.08, 8.5, 0.06, 13.0, 4.2])

    result = extract_pyqz_result(values, columns)

    assert result["o_h"] == pytest.approx(8.5)
    assert result["o_h_err"] == pytest.approx(0.06)
    assert result["log_q"] == pytest.approx(7.4)
    assert result["log_u"] == pytest.approx(7.4 - 10.47712125472)
    assert result["log_u_err"] == pytest.approx(0.08)
    assert result["flag"] == 13
    assert result["rs_offgrid"] == pytest.approx(4.2)
    assert result["valid"] == 1


def test_nebulabayes_extractor_uses_continuous_means_and_preserves_modes_and_qc():
    from model_grid_diagnostics import extract_nebulabayes_result

    estimates = pd.DataFrame(
        {
            "Estimate": [8.55, -3.1, 6.4],
            "CI68_low": [8.45, -3.3, -np.inf],
            "CI68_high": [8.65, -2.9, 6.8],
            "Est_at_lower?": ["N", "Y", "N"],
            "Est_at_upper?": ["N", "N", "N"],
            "n_local_maxima": [1, 2, 1],
        },
        index=["12 + log O/H", "log U", "log P/k"],
    )
    posterior = SimpleNamespace(
        DF_estimates=estimates,
        best_model={"chi2": 1.75},
        Grid_spec=SimpleNamespace(
            paramName2paramValueArr={
                "12 + log O/H": np.array([8.4, 8.55, 8.7]),
                "log U": np.array([-3.4, -3.1, -2.8]),
                "log P/k": np.array([6.0, 6.4, 6.8]),
            }
        ),
        marginalised_1D={
            "12 + log O/H": np.array([1.0, 4.0, 2.0]),
            "log U": np.array([1.0, 3.0, 1.0]),
            "log P/k": np.array([1.0, 2.0, 1.0]),
        },
    )

    result = extract_nebulabayes_result(SimpleNamespace(Posterior=posterior))

    o_h_axis = posterior.Grid_spec.paramName2paramValueArr["12 + log O/H"]
    o_h_pdf = posterior.marginalised_1D["12 + log O/H"]
    expected_o_h_mean = np.trapezoid(o_h_axis * o_h_pdf, o_h_axis) / np.trapezoid(
        o_h_pdf, o_h_axis
    )
    assert result["o_h"] == pytest.approx(expected_o_h_mean)
    assert result["o_h"] != pytest.approx(result["o_h_mode"])
    assert result["o_h_mode"] == pytest.approx(8.55)
    assert result["log_u_mode"] == pytest.approx(-3.1)
    assert result["log_pk_mode"] == pytest.approx(6.4)
    assert result["o_h_ci68_low"] == pytest.approx(8.45)
    assert result["log_u_ci68_high"] == pytest.approx(-2.9)
    assert result["log_pk_ci68_low"] == -np.inf
    assert result["chi2_red"] == pytest.approx(1.75)
    assert result["n_localmax"] == 2
    assert result["flag"] == 2 | 8
    assert result["valid"] == 1
    assert result["o_h_reliable"] == 1
    assert result["log_u_reliable"] == 0
    assert result["log_pk_reliable"] == 0


def test_nebulabayes_empty_result_has_complete_typed_contract():
    from model_grid_diagnostics import (
        NB_FLOAT_FIELDS,
        NB_INTEGER_FIELDS,
        _records_to_arrays,
        empty_nebulabayes_result,
    )

    result = empty_nebulabayes_result(exception=True)

    assert set(result) == set(NB_FLOAT_FIELDS) | set(NB_INTEGER_FIELDS)
    assert all(np.isnan(result[name]) for name in NB_FLOAT_FIELDS)
    assert result["valid"] == 0
    assert result["o_h_reliable"] == 0
    assert result["log_u_reliable"] == 0
    assert result["log_pk_reliable"] == 0
    arrays = _records_to_arrays([result], NB_FLOAT_FIELDS, NB_INTEGER_FIELDS)
    assert all(arrays[name].dtype == np.float64 for name in NB_FLOAT_FIELDS)
    assert all(arrays[name].dtype == np.int16 for name in NB_INTEGER_FIELDS)


def test_integrated_lines_use_one_common_aperture_for_every_line():
    from model_grid_diagnostics import MODEL_LINE_BASES, integrate_line_maps

    shape = (2, 2)
    region = np.ones(shape, dtype=bool)
    flux_maps = {
        name: np.array([[1.0, 2.0], [3.0, 4.0]])
        for name in MODEL_LINE_BASES
    }
    error_maps = {name: np.ones(shape) for name in MODEL_LINE_BASES}
    flux_maps["SII6730"][0, 1] = np.nan
    binid = np.array([[0, 1], [2, 2]], dtype=np.int32)

    integrated = integrate_line_maps(
        flux_maps, error_maps, region, binid_map=binid
    )

    assert integrated.n_pixels == 3
    assert integrated.n_bins == 2
    for name in MODEL_LINE_BASES:
        if name == "SII6730":
            assert integrated.fluxes[name] == pytest.approx(8.0)
        else:
            assert integrated.fluxes[name] == pytest.approx(8.0)
        assert integrated.errors[name] == pytest.approx(np.sqrt(3.0))


def test_ordered_model_hdu_names_keep_immediate_hii_sf_pairs():
    from model_grid_diagnostics import ordered_model_hdu_names

    names = ordered_model_hdu_names()

    assert len(names) == len(set(names))
    for hii_name in [name for name in names if "_HII" in name]:
        sf_name = hii_name.replace("_HII", "_SF", 1)
        assert names.index(sf_name) == names.index(hii_name) + 1
    assert names[0:2] == ["O_H_PYQZ_HII", "O_H_PYQZ_SF"]
    assert "LOG_PK_NEBULABAYES_HII" in names
    assert "O_H_NEBULABAYES_HII_MODE" in names
    assert "NB_OH_RELIABLE_HII" in names
    assert "NB_LOGU_RELIABLE_HII" in names
    assert "NB_LOGPK_RELIABLE_HII" in names
    assert "PYQZ_RS_OFFGRID_SF" in names
    expected_jy22 = [
        "O_H_JY22_HII",
        "O_H_JY22_SF",
        "O_H_JY22_HII_16",
        "O_H_JY22_SF_16",
        "O_H_JY22_HII_84",
        "O_H_JY22_SF_84",
        "O_H_JY22_PREDEP_HII",
        "O_H_JY22_PREDEP_SF",
        "O_H_JY22_PREDEP_HII_16",
        "O_H_JY22_PREDEP_SF_16",
        "O_H_JY22_PREDEP_HII_84",
        "O_H_JY22_PREDEP_SF_84",
        "LOG_U_JY22_HII",
        "LOG_U_JY22_SF",
        "LOG_U_JY22_HII_16",
        "LOG_U_JY22_SF_16",
        "LOG_U_JY22_HII_84",
        "LOG_U_JY22_SF_84",
        "JY22_CHI2_MIN_HII",
        "JY22_CHI2_MIN_SF",
        "JY22_RESID_N2_HII",
        "JY22_RESID_N2_SF",
        "JY22_RESID_S2_HII",
        "JY22_RESID_S2_SF",
        "JY22_RESID_R3_HII",
        "JY22_RESID_R3_SF",
        "JY22_RESID_NORM_HII",
        "JY22_RESID_NORM_SF",
        "JY22_FLAG_HII",
        "JY22_FLAG_SF",
        "JY22_FIT_OK_HII",
        "JY22_FIT_OK_SF",
        "JY22_VALID_HII",
        "JY22_VALID_SF",
    ]
    assert names[-len(expected_jy22) :] == expected_jy22


class NarrowPdfFakePyqz:
    def __init__(self, *, central_on_grid=True, native_flag=1234):
        self.central_on_grid = central_on_grid
        self.native_flag = native_flag
        self.interp_calls = 0
        self.pyqzm = SimpleNamespace(
            M_version="MV",
            PDF_cont_level=0.61,
            QZs_lim={
                "MV": {
                    "LogQ": np.array([6.5, 8.5]),
                    "gas[O]+12": np.array([8.0, 8.875]),
                }
            },
            diagnostics={
                "[NII]/[SII]+;[OIII]/[SII]+": {
                    "coeffs": [[1, 0], [0, 1]]
                }
            },
        )

    def get_global_qz(self, data, data_cols, which_grids, **kwargs):
        del data, data_cols, which_grids, kwargs
        columns = [
            "<LogQ{KDE}>",
            "err(LogQ{KDE})",
            "<gas[O]+12{KDE}>",
            "err(gas[O]+12{KDE})",
            "flag",
            "rs_offgrid",
        ]
        values = np.array(
            [[np.nan, np.nan, np.nan, np.nan, float(self.native_flag), 0.0]]
        )
        return values, columns

    def interp_qz(self, qz, ratio_values, diagnostic, **kwargs):
        del diagnostic, kwargs
        self.interp_calls += 1
        first, second = (
            np.asarray(value, dtype=float) for value in ratio_values
        )
        if qz == "LogQ":
            output = 7.2 + 0.08 * first + 0.03 * second
        else:
            output = 8.5 + 0.02 * first - 0.04 * second
        if not self.central_on_grid:
            output = np.full_like(output, np.nan)
        return output


class RngReplayFakePyqz(NarrowPdfFakePyqz):
    def __init__(self, *, raise_after_sampling=False):
        super().__init__()
        self.raise_after_sampling = raise_after_sampling
        self.native_ratios = None
        self.fallback_ratio_calls = []

    def get_global_qz(self, data, data_cols, which_grids, **kwargs):
        from scipy import stats

        del which_grids
        srs = int(kwargs["srs"])
        line_names = [name for name in data_cols if not name.startswith("std")]
        sampled = np.full((srs + 1, len(line_names)), np.nan)
        for line_index, line_name in enumerate(line_names):
            flux = float(data[0, data_cols.index(line_name)])
            error = float(data[0, data_cols.index(f"std{line_name}")])
            sampled[0, line_index] = flux
            if error > 0.0:
                sampled[1:, line_index] = stats.truncnorm(
                    -flux / error, np.inf, loc=flux, scale=error
                ).rvs(srs)
            elif error == 0.0:
                sampled[1:, line_index] = flux
            else:
                raise ValueError("Unexpected synthetic line error.")
        self.native_ratios = [
            np.log10(sampled[:, 0] / sampled[:, 1]),
            np.log10(sampled[:, 2] / sampled[:, 1]),
        ]
        if self.raise_after_sampling:
            raise UnboundLocalError("synthetic native contour failure")
        return super().get_global_qz(data, data_cols, [], **kwargs)

    def interp_qz(self, qz, ratio_values, diagnostic, **kwargs):
        self.fallback_ratio_calls.append(
            [np.asarray(value).copy() for value in ratio_values]
        )
        return super().interp_qz(qz, ratio_values, diagnostic, **kwargs)


def _narrow_pdf_spectra(bin_id=4):
    from model_grid_diagnostics import BinSpectra

    return BinSpectra(
        bin_ids=np.array([bin_id]),
        fluxes=np.array([[10.0, 5.0, 28.6, 4.0, 3.0, 2.0]]),
        errors=np.full((1, 6), 0.01),
        pixel_counts=np.ones(1, dtype=int),
    )


def test_pyqz_adaptive_fallback_recovers_narrow_on_grid_pdf():
    from model_grid_diagnostics import run_pyqz_spectra

    run = run_pyqz_spectra(
        NarrowPdfFakePyqz(),
        _narrow_pdf_spectra(),
        random_seed=17,
        srs=800,
        adaptive_kde_fallback=True,
    )

    assert run.results["valid"].tolist() == [1]
    assert np.isfinite(run.results["o_h"][0])
    assert np.isfinite(run.results["o_h_err"][0])
    assert np.isfinite(run.results["log_q"][0])
    assert np.isfinite(run.results["log_q_err"][0])
    assert run.failures == ()
    assert run.recoveries and run.recoveries[0][0] == 4


def test_pyqz_adaptive_fallback_is_deterministic_and_restores_numpy_state():
    from model_grid_diagnostics import run_pyqz_spectra

    state_before = np.random.get_state()
    first = run_pyqz_spectra(
        NarrowPdfFakePyqz(),
        _narrow_pdf_spectra(),
        random_seed=17,
        srs=800,
        adaptive_kde_fallback=True,
    )
    state_after = np.random.get_state()
    second = run_pyqz_spectra(
        NarrowPdfFakePyqz(),
        _narrow_pdf_spectra(),
        random_seed=17,
        srs=800,
        adaptive_kde_fallback=True,
    )

    for name in first.results:
        np.testing.assert_array_equal(first.results[name], second.results[name])
    for before, after in zip(state_before, state_after):
        np.testing.assert_array_equal(before, after)


@pytest.mark.parametrize("raise_after_sampling", [False, True])
def test_pyqz_adaptive_fallback_replays_native_mc_draws_exactly(
    raise_after_sampling,
):
    from model_grid_diagnostics import run_pyqz_spectra

    fake = RngReplayFakePyqz(raise_after_sampling=raise_after_sampling)
    spectra = _narrow_pdf_spectra()
    spectra.errors[0, 1] = 0.0

    run = run_pyqz_spectra(
        fake,
        spectra,
        random_seed=29,
        srs=80,
        adaptive_kde_fallback=True,
    )

    assert run.results["valid"].tolist() == [1]
    assert run.results["flag"].tolist() == [13]
    assert run.results["rs_offgrid"].tolist() == [0.0]
    assert run.failures == ()
    assert run.recoveries and run.recoveries[0][0] == 4
    assert len(fake.fallback_ratio_calls) == 2
    for ratio_call in fake.fallback_ratio_calls:
        for expected, actual in zip(fake.native_ratios, ratio_call):
            np.testing.assert_array_equal(actual, expected)
    assert np.ptp(fake.native_ratios[1]) > 0.0


def test_pyqz_adaptive_fallback_keeps_native_flag8_without_replay():
    from model_grid_diagnostics import run_pyqz_spectra

    fake = NarrowPdfFakePyqz(central_on_grid=False, native_flag=8)
    run = run_pyqz_spectra(
        fake,
        _narrow_pdf_spectra(bin_id=8),
        adaptive_kde_fallback=True,
    )

    assert run.results["valid"].tolist() == [0]
    assert run.results["flag"].tolist() == [8]
    assert run.failures == ()
    assert run.recoveries == ()
    assert fake.interp_calls == 0


def test_pyqz_adaptive_fallback_fails_closed_when_direct_ratios_are_off_grid():
    from model_grid_diagnostics import run_pyqz_spectra

    run = run_pyqz_spectra(
        NarrowPdfFakePyqz(central_on_grid=False),
        _narrow_pdf_spectra(bin_id=9),
        adaptive_kde_fallback=True,
    )

    assert run.results["valid"].tolist() == [0]
    assert run.results["flag"].tolist() == [1234]
    assert run.recoveries == ()
    assert len(run.failures) == 1
    assert "Central diagnostic ratios are outside" in run.failures[0][1]


def test_pyqz_adaptive_fallback_is_opt_in():
    from model_grid_diagnostics import run_pyqz_spectra

    fake = NarrowPdfFakePyqz()
    run = run_pyqz_spectra(fake, _narrow_pdf_spectra())

    assert run.results["valid"].tolist() == [0]
    assert run.recoveries == ()
    assert fake.interp_calls == 0


def test_pyqz_batch_adapter_uses_locked_configuration_and_continues_after_failure():
    from model_grid_diagnostics import BinSpectra, run_pyqz_spectra

    class FakePyqz:
        def __init__(self):
            self.calls = []

        def get_global_qz(self, data, data_cols, which_grids, **kwargs):
            self.calls.append((data.copy(), list(data_cols), list(which_grids), kwargs))
            if data[0, 0] == 13.0:
                raise RuntimeError("synthetic off-grid failure")
            jitter = np.random.random()
            columns = [
                "<LogQ{KDE}>",
                "err(LogQ{KDE})",
                "<gas[O]+12{KDE}>",
                "err(gas[O]+12{KDE})",
                "flag",
                "rs_offgrid",
            ]
            values = np.array([[7.2 + jitter, 0.1, 8.4, 0.08, 3.0, 2.5]])
            return values, columns

    fluxes = np.array(
        [
            [10.0, 5.0, 28.6, 4.0, 3.0, 2.0],
            [10.0, 5.0, 28.6, 13.0, 3.0, 2.0],
            [20.0, 6.0, 57.2, 5.0, 4.0, 3.0],
        ]
    )
    errors = np.full_like(fluxes, 0.5)
    spectra = BinSpectra(
        bin_ids=np.array([2, 7, 9]),
        fluxes=fluxes,
        errors=errors,
        pixel_counts=np.ones(3, dtype=int),
    )
    fake = FakePyqz()

    state_before = np.random.get_state()
    run = run_pyqz_spectra(fake, spectra, random_seed=1234, srs=800)
    state_after = np.random.get_state()

    assert len(fake.calls) == 3
    first_data, first_cols, first_diags, first_kwargs = fake.calls[0]
    np.testing.assert_allclose(first_data[0], [4.0, 0.5, 5.0, np.sqrt(0.5), 5.0, 0.5])
    assert first_cols == [
        "[NII]", "std[NII]", "[SII]+", "std[SII]+", "[OIII]", "std[OIII]"
    ]
    assert first_diags == ["[NII]/[SII]+;[OIII]/[SII]+"]
    assert first_kwargs["qzs"] == ["LogQ", "gas[O]+12"]
    assert first_kwargs["Pk"] == 5.0
    assert first_kwargs["struct"] == "pp"
    assert np.isinf(first_kwargs["kappa"])
    assert first_kwargs["sampling"] == 2
    assert first_kwargs["error_pdf"] == "normal"
    assert first_kwargs["srs"] == 800
    assert first_kwargs["flag_level"] == 2.0
    assert first_kwargs["KDE_method"] == "multiv"
    assert first_kwargs["KDE_qz_sampling"] == 101j
    assert first_kwargs["KDE_do_singles"] is True
    assert first_kwargs["nproc"] == 1
    assert run.results["valid"].tolist() == [1, 0, 1]
    assert run.results["flag"].tolist() == [3, -99, 3]
    assert run.failures[0][0] == 7
    assert "synthetic off-grid failure" in run.failures[0][1]
    for before, after in zip(state_before, state_after):
        np.testing.assert_array_equal(before, after)


def test_nebulabayes_batch_adapter_prepropagates_hbeta_and_records_qc():
    from model_grid_diagnostics import BinSpectra, NB_LINE_LIST, run_nebulabayes_spectra

    class FakeNebulaBayesModel:
        def __init__(self):
            self.calls = []

        def __call__(self, fluxes, errors, names, **kwargs):
            self.calls.append((np.array(fluxes), np.array(errors), list(names), kwargs))
            if fluxes[1] > 0.75:
                raise RuntimeError("synthetic posterior failure")
            estimates = pd.DataFrame(
                {
                    "Estimate": [8.5, -3.0, 6.2],
                    "CI68_low": [8.4, -3.2, 5.9],
                    "CI68_high": [8.6, -2.8, 6.5],
                    "Est_at_lower?": ["N", "N", "N"],
                    "Est_at_upper?": ["N", "N", "N"],
                    "n_local_maxima": [1, 1, 2],
                },
                index=["12 + log O/H", "log U", "log P/k"],
            )
            posterior = SimpleNamespace(DF_estimates=estimates, best_model={"chi2": 1.2})
            posterior.Grid_spec = SimpleNamespace(
                paramName2paramValueArr={
                    "12 + log O/H": np.array([8.4, 8.5, 8.6]),
                    "log U": np.array([-3.2, -3.0, -2.8]),
                    "log P/k": np.array([5.9, 6.2, 6.5]),
                }
            )
            posterior.marginalised_1D = {
                "12 + log O/H": np.array([1.0, 4.0, 1.0]),
                "log U": np.array([1.0, 4.0, 1.0]),
                "log P/k": np.array([1.0, 4.0, 1.0]),
            }
            return SimpleNamespace(Posterior=posterior)

    spectra = BinSpectra(
        bin_ids=np.array([4, 5]),
        fluxes=np.array(
            [
                [10.0, 5.0, 28.6, 4.0, 3.0, 2.0],
                [10.0, 8.0, 28.6, 4.0, 3.0, 2.0],
            ]
        ),
        errors=np.array([[1.0, 0.5, 2.0, 0.4, 0.3, 0.2]] * 2),
        pixel_counts=np.ones(2, dtype=int),
    )
    model = FakeNebulaBayesModel()

    run = run_nebulabayes_spectra(model, spectra)

    expected_flux = spectra.fluxes[0] / spectra.fluxes[0, 0]
    expected_err = expected_flux * np.sqrt(
        (spectra.errors[0] / spectra.fluxes[0]) ** 2
        + (spectra.errors[0, 0] / spectra.fluxes[0, 0]) ** 2
    )
    np.testing.assert_allclose(model.calls[0][0], expected_flux)
    np.testing.assert_allclose(model.calls[0][1], expected_err)
    assert model.calls[0][2] == list(NB_LINE_LIST)
    assert model.calls[0][3] == {
        "norm_line": "Hbeta",
        "deredden": False,
        "prior": "Uniform",
        "likelihood_lines": None,
        "verbosity": "ERROR",
    }
    assert run.results["valid"].tolist() == [1, 0]
    assert run.results["n_localmax"].tolist() == [2, -99]
    assert run.results["flag"].tolist() == [0, 16]
    assert run.failures[0][0] == 5


def test_sfr_wires_fixed_model_grid_fields_and_provenance():
    source = (Path(__file__).resolve().parents[1] / "SFR+Z.py").read_text()

    required = (
        'f"O_H_NEBULABAYES_{region}_MODE": "o_h_mode"',
        'f"NB_OH_RELIABLE_{region}": "o_h_reliable"',
        'f"NB_LOGU_RELIABLE_{region}": "log_u_reliable"',
        'f"NB_LOGPK_RELIABLE_{region}": "log_pk_reliable"',
        'f"O_H_JY22_PREDEP_{region}": "o_h_pre"',
        'f"JY22_RESID_N2_{region}": "resid_n2"',
        'f"JY22_FIT_OK_{region}": "fit_ok"',
        "result['resid_norm']",
        "header['NBEST']",
        "header['JY22DEP']",
        "header['JY22OH']",
        "header['JY22FMAX']",
        "header['JY22SHP']",
        "JY22_INFERENCE_GRID_SHAPE = (200, 200)",
        "JY22_GRID = resample_jy22_grid(",
        'or "RELIABLE" in name',
        'or "FIT_OK" in name',
    )
    missing = [snippet for snippet in required if snippet not in source]
    assert not missing, f"SFR+Z.py is missing model-grid contract snippets: {missing}"
    assert "empirical chi2<=25 gross-mismatch cut" in source
    assert "chi2>25 empirical gross mismatch" in source
    assert source.count("adaptive_kde_fallback=True") == 1
    assert "adaptive-local-KDE recovery" in source
    spatial_pyqz_block = source[
        source.index("pyqz_bin_run = run_pyqz_spectra(") : source.index(
            'print_model_failures("pyqz bins"'
        )
    ]
    assert "adaptive_kde_fallback" not in spatial_pyqz_block
    assert "JY22 corrects SII6716 and SII6730 separately" not in source
    assert "shared Balmer correction" not in source
