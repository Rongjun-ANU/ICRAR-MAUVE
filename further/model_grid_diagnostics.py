"""Pure preparation and result helpers for MAUVE model-grid diagnostics.

The package-version and runtime compatibility work lives in
``model_grid_compat.py``.  This module contains the science-facing data
contract used by ``SFR+Z.py`` and is intentionally importable without loading
pyqz or NebulaBayes.
"""

from __future__ import annotations

import hashlib
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

import numpy as np


MODEL_LINE_BASES = (
    "HB4861",
    "OIII5006",
    "HA6562",
    "NII6583",
    "SII6716",
    "SII6730",
)
NB_LINE_LIST = (
    "Hbeta",
    "OIII5007",
    "Halpha",
    "NII6583",
    "SII6716",
    "SII6731",
)

NOT_EVALUATED_FLAG = -99
# Keep the conventional q = U c conversion used by the pyqz/MAPPINGS
# literature, where c is rounded to 3.0e10 cm s^-1.
LOG10_C_CM_S = 10.47712125472

NB_FLAG_OH_EDGE = 1
NB_FLAG_LOGU_EDGE = 2
NB_FLAG_LOGPK_EDGE = 4
NB_FLAG_OPEN_CI68 = 8
NB_FLAG_EXCEPTION = 16

JY22_FLAG_OH_EDGE = 1
JY22_FLAG_LOGU_EDGE = 2
JY22_FLAG_POSTERIOR_INVALID = 4
JY22_FLAG_COVARIANCE_INVALID = 8
JY22_FLAG_EXCEPTION = 16

PYQZ_FLOAT_FIELDS = (
    "o_h",
    "o_h_err",
    "log_q",
    "log_q_err",
    "log_u",
    "log_u_err",
    "rs_offgrid",
)
PYQZ_INTEGER_FIELDS = ("flag", "valid")
NB_FLOAT_FIELDS = (
    "o_h",
    "o_h_ci68_low",
    "o_h_ci68_high",
    "log_u",
    "log_u_ci68_low",
    "log_u_ci68_high",
    "log_pk",
    "log_pk_ci68_low",
    "log_pk_ci68_high",
    "chi2_red",
)
NB_INTEGER_FIELDS = ("n_localmax", "flag", "valid")
JY22_FLOAT_FIELDS = (
    "o_h",
    "o_h_16",
    "o_h_84",
    "log_u",
    "log_u_16",
    "log_u_84",
    "chi2_min",
)
JY22_INTEGER_FIELDS = ("flag", "valid")
JY22_RATIO_NAMES = ("N2", "S2", "R3")
JY22_REQUIRED_GRID_COLUMNS = (
    "grid_kind",
    "log_Z_Zsun",
    "log_U",
    "log_R3_OIII5007_Hbeta",
    "log_N2_NII6584_Halpha",
    "log_S2_SII6717_6731_Halpha",
)


@dataclass(frozen=True)
class BinSpectra:
    """One representative six-line spectrum for each eligible spatial bin."""

    bin_ids: np.ndarray
    fluxes: np.ndarray
    errors: np.ndarray
    pixel_counts: np.ndarray
    line_names: tuple[str, ...] = MODEL_LINE_BASES


@dataclass(frozen=True)
class IntegratedLines:
    """Raw line sums and quadrature errors over one common aperture."""

    fluxes: dict[str, float]
    errors: dict[str, float]
    valid_mask: np.ndarray
    n_pixels: int
    n_bins: int


@dataclass(frozen=True)
class ModelBatchRun:
    """One complete model pass, including isolated per-spectrum failures."""

    results: dict[str, np.ndarray]
    failures: tuple[tuple[int, str], ...]


@dataclass(frozen=True)
class JY22Grid:
    """Validated Peng2026/JY22 grid restricted to its documented log-U range."""

    source_path: Path
    sha256: str
    log_z_zsun: np.ndarray
    o_h: np.ndarray
    log_u: np.ndarray
    model_ratios: np.ndarray
    requested_log_u_bounds: tuple[float, float]
    effective_log_u_bounds: tuple[float, float]
    solar_o_h: float


@dataclass(frozen=True)
class JY22RatioMeasurement:
    """Dust-corrected JY22 ratios and their raw-flux Jacobian covariance."""

    log_ratios: np.ndarray
    covariance: np.ndarray
    ebv: float
    corrected_fluxes: np.ndarray
    jacobian: np.ndarray


@dataclass(frozen=True)
class JY22RatioSpectra:
    """JY22 N2/S2/R3 measurements for one row per spatial bin or aperture."""

    bin_ids: np.ndarray
    log_ratios: np.ndarray
    covariances: np.ndarray
    pixel_counts: np.ndarray
    ebv: np.ndarray


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_jy22_grid(
    path: str | Path,
    *,
    min_log_u: float = -4.0,
    max_log_u: float = -0.5,
    solar_o_h: float = 8.69,
) -> JY22Grid:
    """Load and validate the released interpolated JY22 photoionization grid."""

    from astropy.io import fits

    source_path = Path(path).expanduser().resolve()
    if not source_path.is_file():
        raise FileNotFoundError(f"JY22 grid does not exist: {source_path}")
    if not np.isfinite([min_log_u, max_log_u, solar_o_h]).all():
        raise ValueError("JY22 grid bounds and solar O/H must be finite.")
    if min_log_u > max_log_u:
        raise ValueError("JY22 minimum log U must not exceed its maximum.")

    with fits.open(source_path, memmap=True) as hdul:
        if len(hdul) < 2 or not hasattr(hdul[1], "columns"):
            raise ValueError("JY22 grid must contain a binary table in HDU 1.")
        columns = tuple(hdul[1].columns.names or ())
        _require_columns(columns, JY22_REQUIRED_GRID_COLUMNS, "JY22 grid")
        table = hdul[1].data
        grid_kind = np.char.strip(np.asarray(table["grid_kind"]).astype(str))
        log_z = np.asarray(table["log_Z_Zsun"], dtype=np.float64)
        log_u = np.asarray(table["log_U"], dtype=np.float64)
        interpolated = np.char.lower(grid_kind) == "interpolated"
        in_requested_range = (log_u >= min_log_u) & (log_u <= max_log_u)
        if np.any(in_requested_range & ~interpolated):
            unexpected = np.unique(grid_kind[in_requested_range & ~interpolated])
            raise ValueError(
                "JY22 grid contains non-interpolated rows in the requested "
                f"log-U range: {unexpected.tolist()}."
            )
        selected = interpolated & in_requested_range
        if not np.any(selected):
            raise ValueError(
                "JY22 grid has no interpolated rows in the requested log-U range "
                f"[{min_log_u}, {max_log_u}]."
            )
        selected_columns = {
            name: np.asarray(table[name][selected], dtype=np.float64)
            for name in JY22_REQUIRED_GRID_COLUMNS[1:]
        }

    finite_values = np.column_stack(tuple(selected_columns.values()))
    if not np.all(np.isfinite(finite_values)):
        raise ValueError("JY22 grid contains non-finite selected values.")

    selected_log_z = selected_columns["log_Z_Zsun"]
    selected_log_u = selected_columns["log_U"]
    parameter_pairs = np.column_stack((selected_log_z, selected_log_u))
    if np.unique(parameter_pairs, axis=0).shape[0] != parameter_pairs.shape[0]:
        raise ValueError("JY22 grid contains duplicate (log Z/Zsun, log U) rows.")

    log_z_axis = np.unique(selected_log_z)
    log_u_axis = np.unique(selected_log_u)
    expected_rows = log_z_axis.size * log_u_axis.size
    if parameter_pairs.shape[0] != expected_rows:
        raise ValueError(
            "JY22 grid is not a complete rectangular log Z/Zsun by log U grid."
        )

    model_ratios = np.full(
        (log_z_axis.size, log_u_axis.size, len(JY22_RATIO_NAMES)),
        np.nan,
        dtype=np.float64,
    )
    ratio_columns = (
        "log_N2_NII6584_Halpha",
        "log_S2_SII6717_6731_Halpha",
        "log_R3_OIII5007_Hbeta",
    )
    z_indices = np.searchsorted(log_z_axis, selected_log_z)
    u_indices = np.searchsorted(log_u_axis, selected_log_u)
    for ratio_index, column in enumerate(ratio_columns):
        model_ratios[z_indices, u_indices, ratio_index] = selected_columns[column]
    if not np.all(np.isfinite(model_ratios)):
        raise ValueError("JY22 grid is not a complete rectangular ratio grid.")

    return JY22Grid(
        source_path=source_path,
        sha256=_sha256_file(source_path),
        log_z_zsun=log_z_axis,
        o_h=log_z_axis + float(solar_o_h),
        log_u=log_u_axis,
        model_ratios=model_ratios,
        requested_log_u_bounds=(float(min_log_u), float(max_log_u)),
        effective_log_u_bounds=(float(log_u_axis[0]), float(log_u_axis[-1])),
        solar_o_h=float(solar_o_h),
    )


def jy22_log_ratios_and_covariance(
    fluxes: Sequence[float],
    errors: Sequence[float],
    extinction_coefficients: Sequence[float],
    *,
    intrinsic_balmer_ratio: float = 2.86,
) -> JY22RatioMeasurement:
    """Derive dust-corrected N2/S2/R3 and full first-order covariance.

    Inputs follow ``MODEL_LINE_BASES``.  The covariance is propagated directly
    from the independent raw six-line flux errors, including the shared Balmer
    decrement and denominator terms.  The Balmer-floor branch has zero dust
    derivative, matching the pipeline's ``max(Halpha/Hbeta, 2.86)`` rule.
    """

    flux = np.asarray(fluxes, dtype=np.float64)
    error = np.asarray(errors, dtype=np.float64)
    k_values = np.asarray(extinction_coefficients, dtype=np.float64)
    expected_shape = (len(MODEL_LINE_BASES),)
    if flux.shape != expected_shape:
        raise ValueError(f"JY22 fluxes must have shape {expected_shape}.")
    if error.shape != expected_shape:
        raise ValueError(f"JY22 errors must have shape {expected_shape}.")
    if k_values.shape != expected_shape:
        raise ValueError(
            f"JY22 extinction coefficients must have shape {expected_shape}."
        )
    if not np.all(np.isfinite(flux) & (flux > 0.0)):
        raise ValueError("JY22 raw fluxes must all be finite and positive.")
    if not np.all(np.isfinite(error) & (error > 0.0)):
        raise ValueError("JY22 raw flux errors must all be finite and positive.")
    if not np.all(np.isfinite(k_values)):
        raise ValueError("JY22 extinction coefficients must all be finite.")
    if not np.isfinite(intrinsic_balmer_ratio) or intrinsic_balmer_ratio <= 0.0:
        raise ValueError("JY22 intrinsic Balmer ratio must be finite and positive.")

    hb_index, ha_index = 0, 2
    k_difference = k_values[hb_index] - k_values[ha_index]
    if not np.isfinite(k_difference) or k_difference == 0.0:
        raise ValueError("JY22 requires distinct Hbeta and Halpha extinction values.")

    observed_balmer_ratio = flux[ha_index] / flux[hb_index]
    balmer_ratio = max(observed_balmer_ratio, float(intrinsic_balmer_ratio))
    ebv_coefficient = 2.5 / k_difference
    ebv = float(
        ebv_coefficient * np.log10(balmer_ratio / intrinsic_balmer_ratio)
    )
    d_ebv_d_flux = np.zeros(flux.size, dtype=np.float64)
    if observed_balmer_ratio > intrinsic_balmer_ratio:
        d_ebv_d_flux[hb_index] = -ebv_coefficient / (
            np.log(10.0) * flux[hb_index]
        )
        d_ebv_d_flux[ha_index] = ebv_coefficient / (
            np.log(10.0) * flux[ha_index]
        )

    correction = np.power(10.0, 0.4 * k_values * ebv)
    corrected = flux * correction
    corrected_jacobian = np.diag(correction) + (
        corrected[:, None]
        * (0.4 * np.log(10.0) * k_values[:, None])
        * d_ebv_d_flux[None, :]
    )

    hb, oiii, ha, nii, sii_6716, sii_6730 = corrected
    sii_sum = sii_6716 + sii_6730
    log_ratios = np.asarray(
        [
            np.log10(nii / ha),
            np.log10(sii_sum / ha),
            np.log10(oiii / hb),
        ],
        dtype=np.float64,
    )
    ratio_jacobian = np.empty((len(JY22_RATIO_NAMES), flux.size), dtype=np.float64)
    inverse_ln10 = 1.0 / np.log(10.0)
    ratio_jacobian[0] = inverse_ln10 * (
        corrected_jacobian[3] / nii - corrected_jacobian[2] / ha
    )
    ratio_jacobian[1] = inverse_ln10 * (
        (corrected_jacobian[4] + corrected_jacobian[5]) / sii_sum
        - corrected_jacobian[2] / ha
    )
    ratio_jacobian[2] = inverse_ln10 * (
        corrected_jacobian[1] / oiii - corrected_jacobian[0] / hb
    )
    covariance = ratio_jacobian @ np.diag(error**2) @ ratio_jacobian.T
    covariance = 0.5 * (covariance + covariance.T)
    if not np.all(np.isfinite(log_ratios)) or not np.all(np.isfinite(covariance)):
        raise ValueError("JY22 ratio or covariance calculation is non-finite.")

    return JY22RatioMeasurement(
        log_ratios=log_ratios,
        covariance=covariance,
        ebv=ebv,
        corrected_fluxes=corrected,
        jacobian=ratio_jacobian,
    )


def build_jy22_ratio_spectra(
    spectra: BinSpectra,
    extinction_coefficients: Sequence[float],
    *,
    intrinsic_balmer_ratio: float = 2.86,
) -> JY22RatioSpectra:
    """Convert raw six-line spectra into JY22 ratio/covariance measurements."""

    bin_ids = np.asarray(spectra.bin_ids, dtype=np.int64)
    fluxes = np.asarray(spectra.fluxes, dtype=np.float64)
    errors = np.asarray(spectra.errors, dtype=np.float64)
    pixel_counts = np.asarray(spectra.pixel_counts, dtype=np.int64)
    expected_flux_shape = (bin_ids.size, len(MODEL_LINE_BASES))
    if fluxes.shape != expected_flux_shape or errors.shape != expected_flux_shape:
        raise ValueError(
            "JY22 input spectra must have one six-line flux/error row per BINID."
        )
    if pixel_counts.shape != (bin_ids.size,):
        raise ValueError("JY22 pixel counts must have one value per BINID.")

    ratios = np.full((bin_ids.size, len(JY22_RATIO_NAMES)), np.nan, dtype=np.float64)
    covariances = np.full(
        (bin_ids.size, len(JY22_RATIO_NAMES), len(JY22_RATIO_NAMES)),
        np.nan,
        dtype=np.float64,
    )
    ebv = np.full(bin_ids.size, np.nan, dtype=np.float64)
    for index in range(bin_ids.size):
        try:
            measurement = jy22_log_ratios_and_covariance(
                fluxes[index],
                errors[index],
                extinction_coefficients,
                intrinsic_balmer_ratio=intrinsic_balmer_ratio,
            )
        except ValueError:
            continue
        ratios[index] = measurement.log_ratios
        covariances[index] = measurement.covariance
        ebv[index] = measurement.ebv
    return JY22RatioSpectra(
        bin_ids=bin_ids,
        log_ratios=ratios,
        covariances=covariances,
        pixel_counts=pixel_counts,
        ebv=ebv,
    )


def _normalise_map_dict(
    maps: Mapping[str, Any], label: str
) -> tuple[dict[str, np.ndarray], tuple[int, ...]]:
    missing = [name for name in MODEL_LINE_BASES if name not in maps]
    if missing:
        raise KeyError(f"Missing {label} maps: {missing}")

    arrays = {
        name: np.asarray(maps[name], dtype=np.float64) for name in MODEL_LINE_BASES
    }
    shape = arrays[MODEL_LINE_BASES[0]].shape
    for name, array in arrays.items():
        if array.shape != shape:
            raise ValueError(
                f"{label} map {name} has shape {array.shape}, expected {shape}."
            )
    return arrays, shape


def build_model_input_valid_mask(
    flux_maps: Mapping[str, Any],
    error_maps: Mapping[str, Any],
    region_mask: Any | None = None,
) -> np.ndarray:
    """Return the common finite-positive six-line model-input aperture."""

    fluxes, shape = _normalise_map_dict(flux_maps, "flux")
    errors, error_shape = _normalise_map_dict(error_maps, "error")
    if error_shape != shape:
        raise ValueError(f"Flux/error map shapes differ: {shape} versus {error_shape}.")

    valid = np.ones(shape, dtype=bool)
    for name in MODEL_LINE_BASES:
        valid &= np.isfinite(fluxes[name]) & (fluxes[name] > 0.0)
        valid &= np.isfinite(errors[name]) & (errors[name] > 0.0)

    if region_mask is not None:
        region = np.asarray(region_mask, dtype=bool)
        if region.shape != shape:
            raise ValueError(
                f"Region mask shape {region.shape} does not match line maps {shape}."
            )
        valid &= region
    return valid


def sum_sii_doublet(
    flux_6716: Any,
    flux_6730: Any,
    error_6716: Any,
    error_6730: Any,
) -> tuple[Any, Any]:
    """Return [S II] 6716+6730 flux and independent-error quadrature."""

    flux = np.asarray(flux_6716) + np.asarray(flux_6730)
    error = np.hypot(np.asarray(error_6716), np.asarray(error_6730))
    if flux.ndim == 0:
        return float(flux), float(error)
    return flux, error


def normalise_nebulabayes_to_hbeta(
    fluxes: Sequence[float],
    errors: Sequence[float],
    line_names: Sequence[str] = NB_LINE_LIST,
) -> tuple[np.ndarray, np.ndarray]:
    """Normalise to Hbeta and include denominator error in each ratio error."""

    flux = np.asarray(fluxes, dtype=np.float64)
    error = np.asarray(errors, dtype=np.float64)
    names = list(line_names)
    if flux.ndim != 1 or error.shape != flux.shape or len(names) != flux.size:
        raise ValueError("NebulaBayes flux, error, and line-name lengths must match.")
    if names.count("Hbeta") != 1:
        raise ValueError("NebulaBayes line list must contain Hbeta exactly once.")
    if not np.all(np.isfinite(flux) & (flux > 0.0)):
        raise ValueError("NebulaBayes fluxes must all be finite and positive.")
    if not np.all(np.isfinite(error) & (error > 0.0)):
        raise ValueError("NebulaBayes errors must all be finite and positive.")

    hb_index = names.index("Hbeta")
    hb_flux = flux[hb_index]
    hb_fractional_error = error[hb_index] / hb_flux
    ratios = flux / hb_flux
    ratio_errors = ratios * np.sqrt((error / flux) ** 2 + hb_fractional_error**2)
    return ratios, ratio_errors


def logq_to_logu(log_q: Any) -> Any:
    """Convert pyqz's dimensional LogQ [cm/s] to dimensionless log U."""

    converted = np.asarray(log_q, dtype=np.float64) - LOG10_C_CM_S
    return float(converted) if converted.ndim == 0 else converted


def extract_unique_bin_spectra(
    binid_map: Any,
    valid_mask: Any,
    flux_maps: Mapping[str, Any],
    error_maps: Mapping[str, Any],
    *,
    rtol: float = 1.0e-10,
    atol: float = 1.0e-12,
) -> BinSpectra:
    """Extract one verified representative spectrum per ascending BINID."""

    binid = np.asarray(binid_map)
    valid = np.asarray(valid_mask, dtype=bool)
    fluxes, shape = _normalise_map_dict(flux_maps, "flux")
    errors, error_shape = _normalise_map_dict(error_maps, "error")
    if binid.shape != shape or valid.shape != shape or error_shape != shape:
        raise ValueError("BINID, validity, flux, and error shapes must match.")

    eligible = valid & np.isfinite(binid) & (binid >= 0)
    eligible_flat_indices = np.flatnonzero(eligible)
    eligible_bin_ids = binid.ravel()[eligible_flat_indices].astype(np.int64)
    sort_order = np.argsort(eligible_bin_ids, kind="stable")
    sorted_flat_indices = eligible_flat_indices[sort_order]
    sorted_bin_ids = eligible_bin_ids[sort_order]
    bin_ids, group_starts, group_counts = np.unique(
        sorted_bin_ids, return_index=True, return_counts=True
    )
    flux_rows: list[list[float]] = []
    error_rows: list[list[float]] = []
    pixel_counts = group_counts.astype(np.int64)

    for bin_id, group_start, group_count in zip(
        bin_ids, group_starts, group_counts, strict=True
    ):
        pixel_indices = sorted_flat_indices[
            group_start : group_start + group_count
        ]
        flux_row: list[float] = []
        error_row: list[float] = []
        for name in MODEL_LINE_BASES:
            line_fluxes = fluxes[name].ravel()[pixel_indices]
            line_errors = errors[name].ravel()[pixel_indices]
            if line_fluxes.size == 0:
                raise ValueError(f"BINID {bin_id} has no eligible pixels.")
            if not np.allclose(
                line_fluxes, line_fluxes[0], rtol=rtol, atol=atol, equal_nan=False
            ):
                raise ValueError(f"{name} flux is inconsistent within BINID {bin_id}.")
            if not np.allclose(
                line_errors, line_errors[0], rtol=rtol, atol=atol, equal_nan=False
            ):
                raise ValueError(f"{name} error is inconsistent within BINID {bin_id}.")
            flux_row.append(float(line_fluxes[0]))
            error_row.append(float(line_errors[0]))
        flux_rows.append(flux_row)
        error_rows.append(error_row)

    n_lines = len(MODEL_LINE_BASES)
    return BinSpectra(
        bin_ids=np.asarray(bin_ids, dtype=np.int64),
        fluxes=np.asarray(flux_rows, dtype=np.float64).reshape((-1, n_lines)),
        errors=np.asarray(error_rows, dtype=np.float64).reshape((-1, n_lines)),
        pixel_counts=pixel_counts,
    )


def broadcast_bin_results(
    binid_map: Any,
    bin_ids: Sequence[int],
    results: Mapping[str, Sequence[Any]],
    *,
    integer_fields: set[str] | frozenset[str] | None = None,
) -> dict[str, np.ndarray]:
    """Broadcast one value per BINID into full-size float or int16 maps."""

    binid = np.asarray(binid_map)
    ids = np.asarray(bin_ids, dtype=np.int64)
    integers = set() if integer_fields is None else set(integer_fields)
    output: dict[str, np.ndarray] = {}

    if np.unique(ids).size != ids.size:
        raise ValueError("BINIDs for result broadcasting must be unique.")
    sort_order = np.argsort(ids)
    sorted_ids = ids[sort_order]
    eligible = np.isfinite(binid) & (binid >= 0)
    eligible_flat_indices = np.flatnonzero(eligible)
    eligible_bin_ids = binid.ravel()[eligible_flat_indices].astype(np.int64)
    positions = np.searchsorted(sorted_ids, eligible_bin_ids)
    matched = positions < sorted_ids.size
    matched_positions = np.flatnonzero(matched)
    matched[matched_positions] &= (
        sorted_ids[positions[matched_positions]]
        == eligible_bin_ids[matched_positions]
    )
    target_flat_indices = eligible_flat_indices[matched]
    result_positions = positions[matched]

    for name, field_values in results.items():
        values = np.asarray(field_values)
        if values.shape != (ids.size,):
            raise ValueError(
                f"Result field {name} has shape {values.shape}, expected {(ids.size,)}."
            )
        if name in integers:
            result_map = np.full(binid.shape, NOT_EVALUATED_FLAG, dtype=np.int16)
        else:
            result_map = np.full(binid.shape, np.nan, dtype=np.float64)
        sorted_values = values[sort_order]
        result_map.ravel()[target_flat_indices] = sorted_values[result_positions]
        output[name] = result_map
    return output


def _require_columns(columns: Sequence[str], required: Sequence[str], label: str) -> None:
    missing = [name for name in required if name not in columns]
    if missing:
        raise KeyError(f"{label} result is missing columns: {missing}")


def extract_pyqz_result(
    values: Sequence[float], columns: Sequence[str]
) -> dict[str, float | int]:
    """Extract the combined KDE pyqz values and their native QC fields."""

    required = (
        "<LogQ{KDE}>",
        "err(LogQ{KDE})",
        "<gas[O]+12{KDE}>",
        "err(gas[O]+12{KDE})",
        "flag",
        "rs_offgrid",
    )
    _require_columns(columns, required, "pyqz")
    row = np.asarray(values, dtype=np.float64)
    if row.shape != (len(columns),):
        raise ValueError(f"pyqz result row has shape {row.shape}; expected {(len(columns),)}.")
    lookup = {name: float(row[index]) for index, name in enumerate(columns)}

    log_q = lookup["<LogQ{KDE}>"]
    log_q_err = lookup["err(LogQ{KDE})"]
    o_h = lookup["<gas[O]+12{KDE}>"]
    o_h_err = lookup["err(gas[O]+12{KDE})"]
    valid = int(np.all(np.isfinite([log_q, log_q_err, o_h, o_h_err])))
    raw_flag = lookup["flag"]
    flag = int(raw_flag) if np.isfinite(raw_flag) else NOT_EVALUATED_FLAG
    return {
        "o_h": o_h,
        "o_h_err": o_h_err,
        "log_q": log_q,
        "log_q_err": log_q_err,
        "log_u": logq_to_logu(log_q),
        "log_u_err": log_q_err,
        "flag": flag,
        "rs_offgrid": lookup["rs_offgrid"],
        "valid": valid,
    }


def _truthy_table_value(value: Any) -> bool:
    if isinstance(value, str):
        return value.strip().upper() in {"Y", "YES", "TRUE", "T", "1"}
    return bool(value)


def extract_nebulabayes_result(result: Any) -> dict[str, float | int]:
    """Extract marginal peaks, 68% intervals, and objective NebulaBayes QC."""

    table = result.Posterior.DF_estimates
    parameters = {
        "o_h": "12 + log O/H",
        "log_u": "log U",
        "log_pk": "log P/k",
    }
    required_columns = (
        "Estimate",
        "CI68_low",
        "CI68_high",
        "Est_at_lower?",
        "Est_at_upper?",
        "n_local_maxima",
    )
    _require_columns(list(table.columns), required_columns, "NebulaBayes")
    missing_parameters = [name for name in parameters.values() if name not in table.index]
    if missing_parameters:
        raise KeyError(f"NebulaBayes result is missing parameters: {missing_parameters}")

    output: dict[str, float | int] = {}
    flag = 0
    edge_bits = {
        "o_h": NB_FLAG_OH_EDGE,
        "log_u": NB_FLAG_LOGU_EDGE,
        "log_pk": NB_FLAG_LOGPK_EDGE,
    }
    local_maxima: list[int] = []
    estimates: list[float] = []
    interval_bounds: list[float] = []
    for stem, parameter in parameters.items():
        row = table.loc[parameter]
        estimate = float(row["Estimate"])
        low = float(row["CI68_low"])
        high = float(row["CI68_high"])
        output[stem] = estimate
        output[f"{stem}_ci68_low"] = low
        output[f"{stem}_ci68_high"] = high
        estimates.append(estimate)
        interval_bounds.extend((low, high))
        local_maxima.append(int(row["n_local_maxima"]))
        if _truthy_table_value(row["Est_at_lower?"]) or _truthy_table_value(
            row["Est_at_upper?"]
        ):
            flag |= edge_bits[stem]

    if not np.all(np.isfinite(interval_bounds)):
        flag |= NB_FLAG_OPEN_CI68
    output["chi2_red"] = float(result.Posterior.best_model["chi2"])
    output["n_localmax"] = max(local_maxima)
    output["flag"] = flag
    output["valid"] = int(np.all(np.isfinite(estimates)))
    return output


def empty_pyqz_result(*, exception: bool = False) -> dict[str, float | int]:
    """Return a complete failed/not-evaluated pyqz result record."""

    del exception
    result: dict[str, float | int] = {name: np.nan for name in PYQZ_FLOAT_FIELDS}
    result.update({"flag": NOT_EVALUATED_FLAG, "valid": 0})
    return result


def empty_nebulabayes_result(*, exception: bool = False) -> dict[str, float | int]:
    """Return a complete failed/not-evaluated NebulaBayes result record."""

    result: dict[str, float | int] = {name: np.nan for name in NB_FLOAT_FIELDS}
    result.update(
        {
            "n_localmax": NOT_EVALUATED_FLAG,
            "flag": NB_FLAG_EXCEPTION if exception else NOT_EVALUATED_FLAG,
            "valid": 0,
        }
    )
    return result


def empty_jy22_result(
    *,
    exception: bool = False,
    covariance_invalid: bool = False,
    posterior_invalid: bool = False,
) -> dict[str, float | int]:
    """Return a complete failed/not-evaluated JY22 result record."""

    result: dict[str, float | int] = {name: np.nan for name in JY22_FLOAT_FIELDS}
    if exception:
        flag = JY22_FLAG_EXCEPTION
    elif covariance_invalid:
        flag = JY22_FLAG_COVARIANCE_INVALID
    elif posterior_invalid:
        flag = JY22_FLAG_POSTERIOR_INVALID
    else:
        flag = NOT_EVALUATED_FLAG
    result.update({"flag": flag, "valid": 0})
    return result


def _records_to_arrays(
    records: Sequence[Mapping[str, float | int]],
    float_fields: Sequence[str],
    integer_fields: Sequence[str],
) -> dict[str, np.ndarray]:
    return {
        **{
            name: np.asarray([record[name] for record in records], dtype=np.float64)
            for name in float_fields
        },
        **{
            name: np.asarray([record[name] for record in records], dtype=np.int16)
            for name in integer_fields
        },
    }


def _weighted_axis_quantile(
    axis: np.ndarray, marginal_probability: np.ndarray, quantile: float
) -> float:
    """Interpolate a quantile through the cumulative marginal posterior."""

    values = np.asarray(axis, dtype=np.float64)
    mass = np.asarray(marginal_probability, dtype=np.float64)
    positive = np.isfinite(values) & np.isfinite(mass) & (mass > 0.0)
    if not np.any(positive):
        return np.nan
    values = values[positive]
    mass = mass[positive]
    mass /= np.sum(mass)
    if values.size == 1:
        return float(values[0])
    cumulative_probability = np.cumsum(mass)
    return float(
        np.interp(
            quantile,
            cumulative_probability,
            values,
            left=values[0],
            right=values[-1],
        )
    )


def _marginal_touches_edge(
    axis: np.ndarray,
    marginal_probability: np.ndarray,
    low_quantile: float,
    high_quantile: float,
) -> bool:
    """Flag a marginal whose mode or equal-tailed interval touches a grid edge."""

    axis = np.asarray(axis, dtype=np.float64)
    marginal = np.asarray(marginal_probability, dtype=np.float64)
    if axis.size == 0 or not np.any(np.isfinite(marginal)):
        return True
    maximum = np.nanmax(marginal)
    modes = np.flatnonzero(
        np.isclose(marginal, maximum, rtol=1.0e-12, atol=1.0e-15)
    )
    mode_at_edge = bool(np.any((modes == 0) | (modes == axis.size - 1)))
    interval_at_edge = bool(
        np.isclose(low_quantile, axis[0], rtol=0.0, atol=1.0e-12)
        or np.isclose(high_quantile, axis[-1], rtol=0.0, atol=1.0e-12)
    )
    return mode_at_edge or interval_at_edge


def infer_jy22_posterior(
    log_ratios: Sequence[float],
    covariance: Any,
    grid: JY22Grid,
) -> dict[str, float | int]:
    """Infer posterior-mean O/H and log U using the full ratio covariance."""

    observed = np.asarray(log_ratios, dtype=np.float64)
    covariance_array = np.asarray(covariance, dtype=np.float64)
    if covariance_array.shape != (
        len(JY22_RATIO_NAMES),
        len(JY22_RATIO_NAMES),
    ):
        return empty_jy22_result(covariance_invalid=True)
    if not np.all(np.isfinite(covariance_array)) or not np.allclose(
        covariance_array,
        covariance_array.T,
        rtol=1.0e-10,
        atol=1.0e-14,
    ):
        return empty_jy22_result(covariance_invalid=True)
    try:
        cholesky = np.linalg.cholesky(covariance_array)
    except np.linalg.LinAlgError:
        return empty_jy22_result(covariance_invalid=True)
    if observed.shape != (len(JY22_RATIO_NAMES),) or not np.all(
        np.isfinite(observed)
    ):
        return empty_jy22_result(posterior_invalid=True)

    model_ratios = np.asarray(grid.model_ratios, dtype=np.float64)
    o_h_axis = np.asarray(grid.o_h, dtype=np.float64)
    log_u_axis = np.asarray(grid.log_u, dtype=np.float64)
    expected_shape = (o_h_axis.size, log_u_axis.size, len(JY22_RATIO_NAMES))
    if model_ratios.shape != expected_shape:
        raise ValueError(
            f"JY22 model ratios have shape {model_ratios.shape}; expected "
            f"{expected_shape}."
        )
    if not (
        np.all(np.isfinite(model_ratios))
        and np.all(np.isfinite(o_h_axis))
        and np.all(np.isfinite(log_u_axis))
    ):
        return empty_jy22_result(posterior_invalid=True)

    residuals = (model_ratios - observed).reshape((-1, len(JY22_RATIO_NAMES)))
    whitened = np.linalg.solve(cholesky, residuals.T).T
    chi2 = np.sum(whitened**2, axis=1).reshape(model_ratios.shape[:2])
    if not np.any(np.isfinite(chi2)):
        return empty_jy22_result(posterior_invalid=True)
    chi2_min = float(np.nanmin(chi2))
    unnormalised = np.exp(-0.5 * (chi2 - chi2_min))
    normalisation = float(np.sum(unnormalised))
    if not np.isfinite(normalisation) or normalisation <= 0.0:
        return empty_jy22_result(posterior_invalid=True)
    posterior = unnormalised / normalisation
    marginal_o_h = np.sum(posterior, axis=1)
    marginal_log_u = np.sum(posterior, axis=0)

    o_h_16 = _weighted_axis_quantile(o_h_axis, marginal_o_h, 0.16)
    o_h_84 = _weighted_axis_quantile(o_h_axis, marginal_o_h, 0.84)
    log_u_16 = _weighted_axis_quantile(log_u_axis, marginal_log_u, 0.16)
    log_u_84 = _weighted_axis_quantile(log_u_axis, marginal_log_u, 0.84)
    result: dict[str, float | int] = {
        "o_h": float(np.sum(marginal_o_h * o_h_axis)),
        "o_h_16": o_h_16,
        "o_h_84": o_h_84,
        "log_u": float(np.sum(marginal_log_u * log_u_axis)),
        "log_u_16": log_u_16,
        "log_u_84": log_u_84,
        "chi2_min": chi2_min,
        "flag": 0,
        "valid": 1,
    }
    if not np.all(
        np.isfinite([result[name] for name in JY22_FLOAT_FIELDS])
    ):
        return empty_jy22_result(posterior_invalid=True)
    if _marginal_touches_edge(o_h_axis, marginal_o_h, o_h_16, o_h_84):
        result["flag"] |= JY22_FLAG_OH_EDGE
    if _marginal_touches_edge(
        log_u_axis, marginal_log_u, log_u_16, log_u_84
    ):
        result["flag"] |= JY22_FLAG_LOGU_EDGE
    return result


def run_jy22_spectra(
    grid: JY22Grid,
    spectra: JY22RatioSpectra,
    *,
    progress: Callable[[int, int, int], None] | None = None,
) -> ModelBatchRun:
    """Run deterministic JY22 inference in the supplied stable spectrum order."""

    bin_ids = np.asarray(spectra.bin_ids, dtype=np.int64)
    log_ratios = np.asarray(spectra.log_ratios, dtype=np.float64)
    covariances = np.asarray(spectra.covariances, dtype=np.float64)
    if log_ratios.shape != (bin_ids.size, len(JY22_RATIO_NAMES)):
        raise ValueError("JY22 ratios must have one N2/S2/R3 row per BINID.")
    if covariances.shape != (
        bin_ids.size,
        len(JY22_RATIO_NAMES),
        len(JY22_RATIO_NAMES),
    ):
        raise ValueError("JY22 covariances must have one 3x3 matrix per BINID.")

    records: list[dict[str, float | int]] = []
    failures: list[tuple[int, str]] = []
    for index, bin_id in enumerate(bin_ids):
        try:
            records.append(
                infer_jy22_posterior(
                    log_ratios[index], covariances[index], grid
                )
            )
        except Exception as exc:
            records.append(empty_jy22_result(exception=True))
            failures.append((int(bin_id), f"{type(exc).__name__}: {exc}"))
        if progress is not None:
            progress(index + 1, bin_ids.size, int(bin_id))

    return ModelBatchRun(
        results=_records_to_arrays(
            records, JY22_FLOAT_FIELDS, JY22_INTEGER_FIELDS
        ),
        failures=tuple(failures),
    )


def run_pyqz_spectra(
    pyqz: Any,
    spectra: BinSpectra,
    *,
    diagnostic: str = "[NII]/[SII]+;[OIII]/[SII]+",
    qzs: Sequence[str] = ("LogQ", "gas[O]+12"),
    pk: float = 5.0,
    kappa: float = np.inf,
    struct: str = "pp",
    sampling: int = 2,
    error_pdf: str = "normal",
    random_seed: int = 20260823,
    srs: int = 800,
    flag_level: float = 2.0,
    kde_method: str = "multiv",
    kde_qz_sampling: complex = 101j,
    kde_do_singles: bool = True,
    nproc: int = 1,
    progress: Callable[[int, int, int], None] | None = None,
) -> ModelBatchRun:
    """Run the locked pyqz diagnostic once per ascending input BINID."""

    data_columns = [
        "[NII]",
        "std[NII]",
        "[SII]+",
        "std[SII]+",
        "[OIII]",
        "std[OIII]",
    ]
    records: list[dict[str, float | int]] = []
    failures: list[tuple[int, str]] = []
    random_state = np.random.get_state()
    np.random.seed(random_seed)
    try:
        for index, bin_id in enumerate(spectra.bin_ids):
            flux = spectra.fluxes[index]
            error = spectra.errors[index]
            sii_flux, sii_error = sum_sii_doublet(
                flux[4], flux[5], error[4], error[5]
            )
            data = np.asarray(
                [[flux[3], error[3], sii_flux, sii_error, flux[1], error[1]]],
                dtype=np.float64,
            )
            try:
                with warnings.catch_warnings():
                    # Off-grid spectra legitimately leave all-NaN intermediate
                    # samples; their information is retained in the raw pyqz
                    # flag and rs_offgrid fields rather than repeated warnings.
                    warnings.filterwarnings(
                        "ignore", message="Mean of empty slice", category=RuntimeWarning
                    )
                    warnings.filterwarnings(
                        "ignore",
                        message="Degrees of freedom <= 0 for slice.*",
                        category=RuntimeWarning,
                    )
                    values, columns = pyqz.get_global_qz(
                        data,
                        data_columns,
                        [diagnostic],
                        qzs=list(qzs),
                        Pk=pk,
                        kappa=kappa,
                        struct=struct,
                        sampling=sampling,
                        error_pdf=error_pdf,
                        srs=srs,
                        flag_level=flag_level,
                        KDE_method=kde_method,
                        KDE_qz_sampling=kde_qz_sampling,
                        KDE_do_singles=kde_do_singles,
                        verbose=False,
                        nproc=nproc,
                    )
                values = np.asarray(values)
                if values.shape != (1, len(columns)):
                    raise ValueError(
                        f"pyqz returned shape {values.shape}; expected "
                        f"{(1, len(columns))}."
                    )
                records.append(extract_pyqz_result(values[0], columns))
            except Exception as exc:
                records.append(empty_pyqz_result(exception=True))
                failures.append((int(bin_id), f"{type(exc).__name__}: {exc}"))
            if progress is not None:
                progress(index + 1, spectra.bin_ids.size, int(bin_id))
    finally:
        np.random.set_state(random_state)

    return ModelBatchRun(
        results=_records_to_arrays(records, PYQZ_FLOAT_FIELDS, PYQZ_INTEGER_FIELDS),
        failures=tuple(failures),
    )


def run_nebulabayes_spectra(
    model: Any,
    spectra: BinSpectra,
    *,
    norm_line: str = "Hbeta",
    deredden: bool = False,
    prior: str = "Uniform",
    likelihood_lines: Sequence[str] | None = None,
    progress: Callable[[int, int, int], None] | None = None,
) -> ModelBatchRun:
    """Run the locked NebulaBayes HII-grid inference per input BINID."""

    records: list[dict[str, float | int]] = []
    failures: list[tuple[int, str]] = []
    for index, bin_id in enumerate(spectra.bin_ids):
        try:
            fluxes, errors = normalise_nebulabayes_to_hbeta(
                spectra.fluxes[index], spectra.errors[index], NB_LINE_LIST
            )
            result = model(
                fluxes,
                errors,
                list(NB_LINE_LIST),
                norm_line=norm_line,
                deredden=deredden,
                prior=prior,
                likelihood_lines=(
                    None if likelihood_lines is None else list(likelihood_lines)
                ),
                verbosity="ERROR",
            )
            records.append(extract_nebulabayes_result(result))
        except Exception as exc:
            records.append(empty_nebulabayes_result(exception=True))
            failures.append((int(bin_id), f"{type(exc).__name__}: {exc}"))
        if progress is not None:
            progress(index + 1, spectra.bin_ids.size, int(bin_id))

    return ModelBatchRun(
        results=_records_to_arrays(records, NB_FLOAT_FIELDS, NB_INTEGER_FIELDS),
        failures=tuple(failures),
    )


def integrate_line_maps(
    flux_maps: Mapping[str, Any],
    error_maps: Mapping[str, Any],
    region_mask: Any,
    *,
    binid_map: Any | None = None,
) -> IntegratedLines:
    """Sum all six raw lines over one shared finite-positive region aperture."""

    fluxes, shape = _normalise_map_dict(flux_maps, "flux")
    errors, _ = _normalise_map_dict(error_maps, "error")
    valid = build_model_input_valid_mask(fluxes, errors, region_mask)
    n_pixels = int(np.count_nonzero(valid))
    summed_fluxes = {
        name: float(np.sum(fluxes[name][valid])) if n_pixels else np.nan
        for name in MODEL_LINE_BASES
    }
    summed_errors = {
        name: float(np.sqrt(np.sum(np.square(errors[name][valid]))))
        if n_pixels
        else np.nan
        for name in MODEL_LINE_BASES
    }

    n_bins = 0
    if binid_map is not None:
        binid = np.asarray(binid_map)
        if binid.shape != shape:
            raise ValueError(
                f"BINID shape {binid.shape} does not match line maps {shape}."
            )
        eligible_bins = binid[valid & np.isfinite(binid) & (binid >= 0)]
        n_bins = int(np.unique(eligible_bins.astype(np.int64)).size)

    return IntegratedLines(
        fluxes=summed_fluxes,
        errors=summed_errors,
        valid_mask=valid,
        n_pixels=n_pixels,
        n_bins=n_bins,
    )


def _append_region_pair(names: list[str], hii_name: str, sf_name: str) -> None:
    names.extend((hii_name, sf_name))


def ordered_model_hdu_names() -> list[str]:
    """Return the exact ordered FITS extension contract for model products."""

    names: list[str] = []
    for stem in ("O_H_PYQZ", "LOG_Q_PYQZ", "LOG_U_PYQZ"):
        _append_region_pair(names, f"{stem}_HII", f"{stem}_SF")
        _append_region_pair(names, f"{stem}_HII_ERR", f"{stem}_SF_ERR")
    for stem in ("PYQZ_FLAG", "PYQZ_RS_OFFGRID", "PYQZ_VALID"):
        _append_region_pair(names, f"{stem}_HII", f"{stem}_SF")

    for stem in (
        "O_H_NEBULABAYES",
        "LOG_U_NEBULABAYES",
        "LOG_PK_NEBULABAYES",
    ):
        _append_region_pair(names, f"{stem}_HII", f"{stem}_SF")
        _append_region_pair(
            names, f"{stem}_HII_CI68_LOW", f"{stem}_SF_CI68_LOW"
        )
        _append_region_pair(
            names, f"{stem}_HII_CI68_HIGH", f"{stem}_SF_CI68_HIGH"
        )
    for stem in ("NB_CHI2_RED", "NB_NLOCALMAX", "NB_FLAG", "NB_VALID"):
        _append_region_pair(names, f"{stem}_HII", f"{stem}_SF")

    for stem in ("O_H_JY22", "LOG_U_JY22"):
        _append_region_pair(names, f"{stem}_HII", f"{stem}_SF")
        _append_region_pair(names, f"{stem}_HII_16", f"{stem}_SF_16")
        _append_region_pair(names, f"{stem}_HII_84", f"{stem}_SF_84")
    for stem in ("JY22_CHI2_MIN", "JY22_FLAG", "JY22_VALID"):
        _append_region_pair(names, f"{stem}_HII", f"{stem}_SF")
    return names
