#!/usr/bin/env python
"""
proxy_EWHa.py - Build proxy and legacy pseudo EW(Halpha) maps for a MAUVE galaxy.

Changes (2026-03-31)
--------------------
* Replaced the old single-band proxy workflow with a continuum-cube proxy EW.
* The default output is now:
  - {gal}_proxy_EW_maps.fits
* The old pseudo-EW calculation is still preserved in separate comparison HDUs.
*
* Changes (2026-04-01)
* --------------------
* The proxy continuum window is now re-centered for each spaxel using the
* Halpha velocity map `HA6562_VEL` on top of the galaxy systemic redshift.
*
* Inputs:
  - Observed Halpha map from:
      {gal}_gas_BIN_maps_extended.fits, or
      {gal}_gas_BIN_maps.fits
    using extension `HA6562_FLUX` (or legacy fallbacks if needed).
  - Observed Halpha velocity map from the same gas FITS:
      `HA6562_VEL`
  - Continuum cube from:
      {gal}_CONTcube.fits
  - Legacy broad-band R proxy input from:
      {gal}_SPATIAL_BINNING_maps_extended.fits
    using extension `MAGNITUDE_R_UNCORRECTED`
  - Galaxy redshift from:
      new_redshifts
*
* Proxy EW(Halpha) calculation:
*
*   1. Start from the MaNGA DAP Halpha EW passband in vacuum:
*        lambda_vac,low  = 6557.6 A
*        lambda_vac,high = 6571.6 A
*        lambda_vac,ref  = 6564.608 A
*
*   2. Convert vacuum wavelengths to standard air using the NIST ASD / Peck &
*      Reeder (1972) refractive-index relation:
*
*        lambda_air = lambda_vac / n
*
*        (n - 1) * 1e8 =
*            8060.51
*          + 2480990 / (132.274  - sigma^2)
*          +   17455.7 / (39.32957 - sigma^2)
*
*      where sigma is the vacuum wavenumber in inverse microns.
*
*   3. Shift the observed continuum cube to the rest frame:
*
*        v_gal = Doppler velocity corresponding to the galaxy redshift z_gal
*        v_spaxel = v_gal + v_Halpha(x, y)
*        z_spaxel = Doppler redshift corresponding to v_spaxel
*
*        lambda_rest(x, y) = lambda_obs / (1 + z_spaxel(x, y))
*
*      For each spaxel, convert the observed flux density to the local
*      rest-frame flux density:
*
*        f_lambda,rest(x, y) = (1 + z_spaxel(x, y)) * f_lambda,obs
*
*   4. Measure the mean continuum flux density in the rest-frame Halpha window
*      for each spaxel:
*
*        <f_lambda,cont> =
*            (1 / Delta_lambda) * Integral[f_lambda,rest(lambda) d lambda]
*
*      In the sampled cube this is implemented as the arithmetic mean of all
*      continuum planes whose wavelength falls inside the selected Halpha EW
*      window after shifting each spaxel by its own z_spaxel(x, y).
*
*   5. Compute the proxy equivalent width:
*
*        EW_proxy(Halpha) = F_Halpha / <f_lambda,cont>
*
*      where F_Halpha is the observed Halpha line flux map converted to
*      erg s^-1 cm^-2.
*
* Legacy pseudo EW(Halpha) preserved for comparison:
*
*   1. Convert `MAGNITUDE_R_UNCORRECTED` to nanomaggies:
*
*        f_nmgy = 10^((22.5 - m_r) / 2.5)
*
*   2. Convert nanomaggies to f_nu:
*
*        f_nu = f_nmgy * 3.631e-29
*              [erg s^-1 cm^-2 Hz^-1]
*
*   3. Approximate the continuum flux density at Halpha using the old single
*      reference wavelength lambda_ref = 6562.8 A:
*
*        f_lambda,legacy = f_nu * c / lambda_ref^2
*
*      with c = 2.99792458e18 A s^-1.
*
*   4. Compute the legacy pseudo equivalent width:
*
*        EW_pseudo(Halpha) = F_Halpha / f_lambda,legacy
*
*      This legacy pseudo EW is kept in the output FITS only as a comparison
*      product against the new continuum-cube proxy EW.
"""

import argparse
import csv
from pathlib import Path
from typing import Any, cast

import numpy as np
from astropy.io import fits


HALPHA_RITZ_VAC_A = 6564.608
HALPHA_EW_VAC_RANGE_A = (6557.6, 6571.6)
LEGACY_HALPHA_REF_A = 6562.8
DEFAULT_REDSHIFT_FILE = "new_redshifts"
CONT_CGS_UNIT = "erg s-1 cm-2 Angstrom-1"
NMAGGY_TO_FNU = 3.631e-29
C_AA_PER_S = 2.99792458e18
C_KM_PER_S = 2.99792458e5


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create proxy EW(Halpha) maps from the observed Halpha flux map and "
            "the mean continuum flux density measured from the continuum cube "
            "over the rest-frame Halpha EW passband, while also preserving the "
            "legacy broad-band pseudo EW for comparison."
        )
    )
    parser.add_argument(
        "-g",
        "--galaxy",
        default="IC3392",
        help="Galaxy identifier, e.g. IC3392 (default: IC3392)",
    )
    parser.add_argument(
        "--root",
        default=".",
        help="Directory containing the input FITS files (default: current directory)",
    )
    parser.add_argument(
        "--fallback-root",
        default=None,
        help="Fallback directory searched when an input file is not found under --root",
    )
    parser.add_argument(
        "--bin-file",
        default=None,
        help="Optional explicit path to {gal}_SPATIAL_BINNING_maps_extended.fits",
    )
    parser.add_argument(
        "--gas-file",
        default=None,
        help=(
            "Optional explicit path to the gas maps FITS file; if omitted, the script "
            "tries {gal}_gas_BIN_maps_extended.fits and {gal}_gas_BIN_maps.fits"
        ),
    )
    parser.add_argument(
        "--cont-file",
        default=None,
        help="Optional explicit path to {gal}_CONTcube.fits",
    )
    parser.add_argument(
        "--redshift-file",
        default=DEFAULT_REDSHIFT_FILE,
        help=(
            "Path to the galaxy redshift table (default: new_redshifts). "
            "Whitespace two-column tables and CSV tables are both supported."
        ),
    )
    parser.add_argument(
        "--out",
        default=None,
        help="Optional output path (default: {gal}_proxy_EW_maps.fits in current directory)",
    )
    return parser.parse_args()


def _unique_paths(*paths: Path | None) -> list[Path]:
    unique: list[Path] = []
    seen: set[str] = set()

    for path in paths:
        if path is None:
            continue
        resolved = path.expanduser().resolve()
        if str(resolved) in seen:
            continue
        seen.add(str(resolved))
        unique.append(resolved)

    return unique


def resolve_existing_path(label: str, *paths: Path | None) -> Path:
    candidates = _unique_paths(*paths)
    for candidate in candidates:
        if candidate.exists():
            return candidate

    checked = "\n".join(f"  - {candidate}" for candidate in candidates)
    raise FileNotFoundError(f"Could not find {label}. Checked:\n{checked}")


def normalize_galaxy_id(galaxy: str) -> str:
    return galaxy.upper().replace(" ", "")


def vacuum_to_air_wavelength_angstrom(lam_vac_a: float) -> float:
    """
    Convert a vacuum wavelength in Angstrom to standard air using the
    Peck & Reeder (1972) formula adopted by NIST ASD.

    The refractive-index formula uses sigma in inverse microns.
    """
    lam_vac_um = lam_vac_a * 1e-4
    sigma_um_inv = 1.0 / lam_vac_um
    n_minus_1 = 1e-8 * (
        8060.51
        + 2480990.0 / (132.274 - sigma_um_inv**2)
        + 17455.7 / (39.32957 - sigma_um_inv**2)
    )
    return lam_vac_a / (1.0 + n_minus_1)


HALPHA_RITZ_AIR_A = vacuum_to_air_wavelength_angstrom(HALPHA_RITZ_VAC_A)
HALPHA_EW_AIR_RANGE_A = tuple(
    vacuum_to_air_wavelength_angstrom(lam) for lam in HALPHA_EW_VAC_RANGE_A
)


def halpha_flux_to_cgs(line_flux: np.ndarray, bunit: str) -> np.ndarray:
    """Convert Halpha map to erg s^-1 cm^-2 when BUNIT encodes a 1e-20 scale."""
    unit_norm = str(bunit).replace(" ", "").lower()
    scale = 1.0
    if ("1e-20" in unit_norm) or ("10^-20" in unit_norm) or ("10**-20" in unit_norm):
        scale = 1e-20
    return np.asarray(line_flux, dtype=np.float64) * scale


def flux_density_to_cgs(flux_density: np.ndarray, bunit: str) -> np.ndarray:
    """
    Convert continuum flux density to erg s^-1 cm^-2 Angstrom^-1 when BUNIT
    encodes a 1e-20 scale.

    This assumes the cube is already expressed per Angstrom, which matches the
    MAUVE CONTcube format seen locally.
    """
    unit_norm = str(bunit).replace(" ", "").lower()
    scale = 1.0
    if ("1e-20" in unit_norm) or ("10^-20" in unit_norm) or ("10**-20" in unit_norm):
        scale = 1e-20
    return np.asarray(flux_density, dtype=np.float64) * scale


def magnitude_to_nanomaggies(magnitude: np.ndarray) -> np.ndarray:
    return 10 ** ((22.5 - magnitude) / 2.5)


def nanomaggies_to_flambda(nmgy: np.ndarray, lam_a: float) -> np.ndarray:
    f_nu = np.asarray(nmgy, dtype=np.float64) * NMAGGY_TO_FNU
    return f_nu * C_AA_PER_S / (lam_a**2)


def read_observed_halpha(gas_path: Path) -> tuple[np.ndarray, fits.Header, str, str]:
    candidates = (
        "HA6562_FLUX",
        "HA6563_FLUX",
        "HALPHA6563_FLUX",
        "HALPHA_FLUX",
    )

    with fits.open(gas_path) as hdul:
        available = {
            str(getattr(hdu, "name", "")).upper()
            for hdu in hdul
            if getattr(hdu, "name", "")
        }

        for ext in candidates:
            if ext.upper() in available:
                hdu = cast(Any, hdul[ext])
                data = np.asarray(hdu.data, dtype=np.float64)
                header = hdu.header.copy()
                bunit = header.get("BUNIT", "1e-20 erg s-1 cm-2")
                return data, header, ext, bunit

    raise KeyError(
        f"Could not find observed Halpha extension in {gas_path.name}. "
        f"Tried: {candidates}"
    )


def velocity_to_kms(velocity: np.ndarray, bunit: str) -> np.ndarray:
    unit_norm = str(bunit).strip().lower().replace(" ", "")
    factor = 1.0
    if unit_norm in {"", "km/s", "kms-1", "km/s^-1", "kilometer/second"}:
        factor = 1.0
    elif unit_norm in {"m/s", "ms-1", "meter/second"}:
        factor = 1.0e-3
    return np.asarray(velocity, dtype=np.float64) * factor


def read_observed_halpha_velocity(
    gas_path: Path,
) -> tuple[np.ndarray, fits.Header, str, str]:
    candidates = (
        "HA6562_VEL",
        "HA6563_VEL",
        "HALPHA6563_VEL",
        "HALPHA_VEL",
    )

    with fits.open(gas_path) as hdul:
        available = {
            str(getattr(hdu, "name", "")).upper()
            for hdu in hdul
            if getattr(hdu, "name", "")
        }

        for ext in candidates:
            if ext.upper() in available:
                hdu = cast(Any, hdul[ext])
                data = np.asarray(hdu.data, dtype=np.float64)
                header = hdu.header.copy()
                bunit = header.get("BUNIT", "km/s")
                return data, header, ext, bunit

    raise KeyError(
        f"Could not find observed Halpha velocity extension in {gas_path.name}. "
        f"Tried: {candidates}"
    )


def read_observed_r_flux(bin_path: Path) -> tuple[np.ndarray, fits.Header, str, str]:
    with fits.open(bin_path) as hdul:
        available = {
            str(getattr(hdu, "name", "")).upper()
            for hdu in hdul
            if getattr(hdu, "name", "")
        }

        if "MAGNITUDE_R_UNCORRECTED" in available:
            hdu = cast(Any, hdul["magnitude_r_uncorrected"])
            magnitude = np.asarray(hdu.data, dtype=np.float64)
            data = magnitude_to_nanomaggies(magnitude)
            header = hdu.header.copy()
            return data, header, "magnitude_r_uncorrected", "nanomaggies"

    raise KeyError(
        f"Could not find R-band flux in {bin_path.name}. "
        "Expected MAGNITUDE_R_UNCORRECTED."
    )


def read_redshift_table(redshift_path: Path) -> dict[str, float]:
    with redshift_path.open("r", encoding="utf-8") as handle:
        first_line = handle.readline()
        handle.seek(0)

        if "," in first_line:
            reader = csv.DictReader(handle)
            galaxy_key = None
            redshift_key = None

            if reader.fieldnames is None:
                raise ValueError(f"CSV redshift file {redshift_path} has no header row.")

            for key in reader.fieldnames:
                key_norm = key.strip().lower()
                if galaxy_key is None and key_norm in {"galaxy", "galaxy_id", "id", "target"}:
                    galaxy_key = key
                if redshift_key is None and key_norm in {"redshift", "selected_redshift", "z"}:
                    redshift_key = key

            if galaxy_key is None or redshift_key is None:
                raise ValueError(
                    f"Could not identify galaxy/redshift columns in {redshift_path}."
                )

            table: dict[str, float] = {}
            for row in reader:
                galaxy = normalize_galaxy_id(row[galaxy_key].strip())
                if not galaxy:
                    continue
                table[galaxy] = float(row[redshift_key])
            return table

        table = {}
        for raw_line in handle:
            line = raw_line.strip()
            if (not line) or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            table[normalize_galaxy_id(parts[0])] = float(parts[1])
        return table


def lookup_redshift(galaxy: str, redshift_path: Path) -> float:
    table = read_redshift_table(redshift_path)
    galaxy_norm = normalize_galaxy_id(galaxy)
    if galaxy_norm not in table:
        available = ", ".join(sorted(table)[:10])
        raise KeyError(
            f"Galaxy {galaxy_norm} not found in redshift file {redshift_path}. "
            f"Example entries: {available}"
        )

    z = float(table[galaxy_norm])
    if 1.0 + z <= 0.0:
        raise ValueError(f"Invalid redshift for {galaxy_norm}: z={z}")
    return z


def redshift_to_velocity_kms(z: float) -> float:
    one_plus_z = 1.0 + float(z)
    beta = (one_plus_z**2 - 1.0) / (one_plus_z**2 + 1.0)
    return C_KM_PER_S * beta


def velocity_to_redshift(velocity_kms: np.ndarray | float) -> np.ndarray:
    beta = np.asarray(velocity_kms, dtype=np.float64) / C_KM_PER_S
    if np.any(np.abs(beta[np.isfinite(beta)]) >= 1.0):
        raise ValueError("Velocity implies |v| >= c; cannot convert to redshift.")
    return np.sqrt((1.0 + beta) / (1.0 - beta)) - 1.0


def build_spaxel_redshift_map(galaxy_z: float, ha_vel_kms: np.ndarray) -> np.ndarray:
    spaxel_z = np.full(ha_vel_kms.shape, galaxy_z, dtype=np.float64)
    valid_vel = np.isfinite(ha_vel_kms)
    if not np.any(valid_vel):
        return spaxel_z

    galaxy_vel_kms = redshift_to_velocity_kms(galaxy_z)
    # Follow the requested practical recipe: add the per-spaxel Halpha velocity
    # to the systemic galaxy velocity in km/s, then convert back to redshift.
    total_vel_kms = galaxy_vel_kms + ha_vel_kms[valid_vel]
    spaxel_z[valid_vel] = velocity_to_redshift(total_vel_kms)
    return spaxel_z


def find_continuum_hdu(hdul: fits.HDUList) -> tuple[np.ndarray, fits.Header, str, str]:
    for index, hdu in enumerate(hdul):
        data = getattr(hdu, "data", None)
        if data is None:
            continue
        data_array = np.asarray(data)
        if data_array.ndim != 3:
            continue
        header = hdu.header.copy()
        source = str(getattr(hdu, "name", "") or f"HDU{index}")
        bunit = header.get("BUNIT", "")
        return data_array, header, source, bunit

    raise KeyError("Could not find a 3D continuum cube HDU.")


def find_spectral_axis(header: fits.Header, ndim: int) -> tuple[int, int]:
    for fits_axis in range(1, ndim + 1):
        ctype = str(header.get(f"CTYPE{fits_axis}", "")).upper()
        if any(token in ctype for token in ("AWAV", "WAVE", "LAMB", "FREQ", "VELO")):
            numpy_axis = ndim - fits_axis
            return fits_axis, numpy_axis

    if (ndim == 3) and ("CRVAL3" in header):
        return 3, 0

    raise ValueError("Could not identify the spectral axis in the continuum cube header.")


def wavelength_unit_to_angstrom_factor(unit: str) -> float:
    unit_norm = unit.strip().lower().replace("angstroms", "angstrom")
    if unit_norm in {"", "a", "aa", "angstrom", "ang", "angs"}:
        return 1.0
    if unit_norm in {"nm", "nanometer", "nanometers"}:
        return 10.0
    if unit_norm in {"um", "micron", "microns", "micrometer", "micrometers"}:
        return 1.0e4
    if unit_norm in {"m", "meter", "meters"}:
        return 1.0e10
    raise ValueError(f"Unsupported wavelength unit for continuum cube: {unit!r}")


def build_wavelength_axis_angstrom(
    header: fits.Header, spectral_fits_axis: int, length: int
) -> np.ndarray:
    crval = float(header[f"CRVAL{spectral_fits_axis}"])
    crpix = float(header.get(f"CRPIX{spectral_fits_axis}", 1.0))

    cdelt_key = f"CDELT{spectral_fits_axis}"
    cd_key = f"CD{spectral_fits_axis}_{spectral_fits_axis}"
    if cdelt_key in header:
        delta = float(header[cdelt_key])
    elif cd_key in header:
        delta = float(header[cd_key])
    else:
        raise KeyError(
            f"Continuum cube header is missing {cdelt_key} and {cd_key}."
        )

    factor = wavelength_unit_to_angstrom_factor(
        str(header.get(f"CUNIT{spectral_fits_axis}", "Angstrom"))
    )
    pixel = np.arange(length, dtype=np.float64) + 1.0
    wavelength = crval + (pixel - crpix) * delta
    return wavelength * factor


def spectral_axis_uses_air_wavelength(header: fits.Header, spectral_fits_axis: int) -> bool:
    ctype = str(header.get(f"CTYPE{spectral_fits_axis}", "")).upper()
    return "AWAV" in ctype


def compute_mean_restframe_continuum(
    cont_path: Path,
    spaxel_z: np.ndarray,
    expected_spatial_shape: tuple[int, ...],
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    if spaxel_z.shape != expected_spatial_shape:
        raise ValueError(
            "Redshift-map shape mismatch: "
            f"expected {expected_spatial_shape} but got {spaxel_z.shape}."
        )

    with fits.open(cont_path, memmap=True) as hdul:
        cube, cont_header, cont_source, cont_bunit = find_continuum_hdu(hdul)

        spectral_fits_axis, spectral_numpy_axis = find_spectral_axis(
            cont_header, cube.ndim
        )
        cube_view = np.moveaxis(cube, spectral_numpy_axis, 0)

        if cube_view.shape[1:] != expected_spatial_shape:
            raise ValueError(
                "Spatial shape mismatch: "
                f"Halpha shape {expected_spatial_shape} vs continuum cube shape "
                f"{cube_view.shape[1:]}."
            )

        wave_obs_a = build_wavelength_axis_angstrom(
            cont_header, spectral_fits_axis, cube_view.shape[0]
        )

        use_air_window = spectral_axis_uses_air_wavelength(cont_header, spectral_fits_axis)
        if use_air_window:
            rest_low_a, rest_high_a = HALPHA_EW_AIR_RANGE_A
            rest_ref_a = HALPHA_RITZ_AIR_A
        else:
            rest_low_a, rest_high_a = HALPHA_EW_VAC_RANGE_A
            rest_ref_a = HALPHA_RITZ_VAC_A

        low = min(rest_low_a, rest_high_a)
        high = max(rest_low_a, rest_high_a)
        valid_z = np.isfinite(spaxel_z) & ((1.0 + spaxel_z) > 0.0)
        if not np.any(valid_z):
            raise ValueError(
                "Spaxel redshift map has no finite values with 1+z > 0."
            )

        one_plus_z = 1.0 + spaxel_z
        observed_low_a = low * one_plus_z
        observed_high_a = high * one_plus_z
        observed_ref_a = rest_ref_a * one_plus_z

        global_wave_mask = (wave_obs_a >= np.nanmin(observed_low_a[valid_z])) & (
            wave_obs_a <= np.nanmax(observed_high_a[valid_z])
        )
        nwave_window = int(np.count_nonzero(global_wave_mask))
        if nwave_window == 0:
            raise ValueError(
                "No continuum planes intersect the spaxel-dependent Halpha EW "
                f"window for {cont_path.name}."
            )

        valid_counts = np.zeros(expected_spatial_shape, dtype=np.int32)
        cont_sum = np.zeros(expected_spatial_shape, dtype=np.float64)

        for wave_idx in np.flatnonzero(global_wave_mask):
            plane_mask = (
                valid_z
                & (wave_obs_a[wave_idx] >= observed_low_a)
                & (wave_obs_a[wave_idx] <= observed_high_a)
            )
            if not np.any(plane_mask):
                continue

            plane_rest = flux_density_to_cgs(cube_view[wave_idx], cont_bunit) * one_plus_z
            finite = plane_mask & np.isfinite(plane_rest)
            if not np.any(finite):
                continue

            valid_counts[finite] += 1
            cont_sum[finite] += plane_rest[finite]

        cont_mean = np.full(expected_spatial_shape, np.nan, dtype=np.float64)
        valid = valid_counts > 0
        cont_mean[valid] = cont_sum[valid] / valid_counts[valid]

        meta: dict[str, Any] = {
            "source": cont_source,
            "bunit": cont_bunit,
            "spectral_fits_axis": spectral_fits_axis,
            "spectral_ctype": str(cont_header.get(f"CTYPE{spectral_fits_axis}", "")),
            "spectral_cunit": str(cont_header.get(f"CUNIT{spectral_fits_axis}", "")),
            "window_is_air": use_air_window,
            "rest_low_a": low,
            "rest_high_a": high,
            "rest_ref_a": rest_ref_a,
            "observed_low_min_a": float(np.nanmin(observed_low_a[valid_z])),
            "observed_low_max_a": float(np.nanmax(observed_low_a[valid_z])),
            "observed_high_min_a": float(np.nanmin(observed_high_a[valid_z])),
            "observed_high_max_a": float(np.nanmax(observed_high_a[valid_z])),
            "observed_ref_min_a": float(np.nanmin(observed_ref_a[valid_z])),
            "observed_ref_max_a": float(np.nanmax(observed_ref_a[valid_z])),
            "z_min": float(np.nanmin(spaxel_z[valid_z])),
            "z_max": float(np.nanmax(spaxel_z[valid_z])),
            "nspaxel_z": int(np.count_nonzero(valid_z)),
            "nwave_total": int(cube_view.shape[0]),
            "nwave_window": nwave_window,
        }
        return cont_mean, valid_counts.astype(np.int16), meta


def compute_legacy_pseudo_ew(
    obs_ha_cgs: np.ndarray,
    obs_r_nmgy: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    legacy_cont = nanomaggies_to_flambda(obs_r_nmgy, LEGACY_HALPHA_REF_A)
    legacy_pseudo_ew = np.full(obs_ha_cgs.shape, np.nan, dtype=np.float64)
    valid = np.isfinite(obs_ha_cgs) & np.isfinite(legacy_cont) & (legacy_cont > 0.0)
    legacy_pseudo_ew[valid] = obs_ha_cgs[valid] / legacy_cont[valid]
    return legacy_cont, legacy_pseudo_ew


def copy_basic_wcs(target_header: fits.Header, source_header: fits.Header) -> None:
    for key in (
        "CTYPE1",
        "CTYPE2",
        "CRVAL1",
        "CRVAL2",
        "CRPIX1",
        "CRPIX2",
        "CD1_1",
        "CD1_2",
        "CD2_1",
        "CD2_2",
        "CDELT1",
        "CDELT2",
        "CUNIT1",
        "CUNIT2",
    ):
        if key in source_header:
            target_header[key] = source_header[key]


def main() -> None:
    args = parse_args()
    gal = normalize_galaxy_id(args.galaxy)
    root = Path(args.root).expanduser().resolve()
    fallback_root = (
        Path(args.fallback_root).expanduser().resolve()
        if args.fallback_root is not None
        else None
    )

    bin_path = (
        Path(args.bin_file).expanduser().resolve()
        if args.bin_file
        else resolve_existing_path(
            "binning FITS",
            root / f"{gal}_SPATIAL_BINNING_maps_extended.fits",
            fallback_root / f"{gal}_SPATIAL_BINNING_maps_extended.fits"
            if fallback_root is not None
            else None,
        )
    )
    gas_path = (
        Path(args.gas_file).expanduser().resolve()
        if args.gas_file
        else resolve_existing_path(
            "gas FITS",
            root / f"{gal}_gas_BIN_maps_extended.fits",
            root / f"{gal}_gas_BIN_maps.fits",
            fallback_root / f"{gal}_gas_BIN_maps_extended.fits"
            if fallback_root is not None
            else None,
            fallback_root / f"{gal}_gas_BIN_maps.fits"
            if fallback_root is not None
            else None,
        )
    )
    cont_path = (
        Path(args.cont_file).expanduser().resolve()
        if args.cont_file
        else resolve_existing_path(
            "continuum cube FITS",
            root / f"{gal}_CONTcube.fits",
            fallback_root / f"{gal}_CONTcube.fits" if fallback_root is not None else None,
        )
    )

    redshift_arg = Path(args.redshift_file).expanduser()
    redshift_path = resolve_existing_path(
        "redshift table",
        redshift_arg if redshift_arg.is_absolute() else Path.cwd() / redshift_arg,
        root / redshift_arg if not redshift_arg.is_absolute() else None,
        fallback_root / redshift_arg
        if (fallback_root is not None) and (not redshift_arg.is_absolute())
        else None,
    )

    out_path = Path(args.out) if args.out else Path(f"{gal}_proxy_EW_maps.fits")
    z = lookup_redshift(gal, redshift_path)

    print("\n=== proxy_EWHa inputs/outputs ===")
    print("Galaxy          :", gal)
    print("Redshift (z)    :", z)
    print("Primary root    :", root)
    if fallback_root is not None:
        print("Fallback root   :", fallback_root)
    print("Binning FITS    :", bin_path)
    print("Gas FITS        :", gas_path)
    print("Continuum FITS  :", cont_path)
    print("Redshift file   :", redshift_path)
    print("Output FITS     :", out_path)
    print(
        "Halpha window (vacuum, rest) : "
        f"{HALPHA_EW_VAC_RANGE_A[0]:.3f} - {HALPHA_EW_VAC_RANGE_A[1]:.3f} A"
    )
    print(
        "Halpha window (air, rest)    : "
        f"{HALPHA_EW_AIR_RANGE_A[0]:.3f} - {HALPHA_EW_AIR_RANGE_A[1]:.3f} A"
    )
    print(
        "Halpha lambda (air, rest)    : "
        f"{HALPHA_RITZ_AIR_A:.3f} A"
    )
    print(
        "Halpha Delta lambda (air)    : "
        f"{HALPHA_EW_AIR_RANGE_A[1] - HALPHA_EW_AIR_RANGE_A[0]:.3f} A"
    )

    galaxy_vel_kms = redshift_to_velocity_kms(z)

    obs_ha, ha_header, ha_extname, ha_bunit = read_observed_halpha(gas_path)
    ha_vel, _vel_header, vel_extname, vel_bunit = read_observed_halpha_velocity(gas_path)
    obs_r_nmgy, r_header, r_source, r_bunit = read_observed_r_flux(bin_path)

    if obs_ha.shape != obs_r_nmgy.shape:
        raise ValueError(
            "Map shape mismatch: "
            f"Halpha shape {obs_ha.shape} vs R-band shape {obs_r_nmgy.shape}."
        )
    if obs_ha.shape != ha_vel.shape:
        raise ValueError(
            "Map shape mismatch: "
            f"Halpha flux shape {obs_ha.shape} vs Halpha velocity shape {ha_vel.shape}."
        )

    ha_vel_kms = velocity_to_kms(ha_vel, vel_bunit)
    spaxel_z = build_spaxel_redshift_map(z, ha_vel_kms)
    cont_mean, cont_counts, cont_meta = compute_mean_restframe_continuum(
        cont_path, spaxel_z, obs_ha.shape
    )

    obs_ha_cgs = halpha_flux_to_cgs(obs_ha, ha_bunit)

    proxy_ewha = np.full(obs_ha.shape, np.nan, dtype=np.float64)
    valid_proxy = np.isfinite(obs_ha_cgs) & np.isfinite(cont_mean) & (cont_mean > 0.0)
    proxy_ewha[valid_proxy] = obs_ha_cgs[valid_proxy] / cont_mean[valid_proxy]

    legacy_cont, pseudo_ewha = compute_legacy_pseudo_ew(obs_ha_cgs, obs_r_nmgy)

    primary = fits.PrimaryHDU(data=obs_ha.astype(np.float64))
    primary.header["EXTNAME"] = "OBS_HA6562_FLUX"
    primary.header["BUNIT"] = ha_bunit
    primary.header["GALAXY"] = gal
    primary.header["BINFILE"] = bin_path.name
    primary.header["GASFILE"] = gas_path.name
    primary.header["CNTFILE"] = cont_path.name
    primary.header["ZFILE"] = redshift_path.name
    primary.header["REDSHIFT"] = (z, "Galaxy systemic redshift")
    primary.header["GALVEL"] = (
        galaxy_vel_kms,
        "Systemic velocity from galaxy redshift [km/s]",
    )
    primary.header["HA_SRC"] = ha_extname
    primary.header["VELSRC"] = vel_extname
    primary.header["RSRC"] = r_source
    primary.header["CNTSRC"] = str(cont_meta["source"])
    primary.header["EWVAC"] = (HALPHA_RITZ_VAC_A, "MaNGA DAP Halpha Ritz wavelength [A]")
    primary.header["EWAIR"] = (
        HALPHA_RITZ_AIR_A,
        "Standard-air Halpha reference wavelength [A]",
    )
    primary.header["EWLOVAC"] = (HALPHA_EW_VAC_RANGE_A[0], "MaNGA DAP EW window low [A]")
    primary.header["EWHIVAC"] = (HALPHA_EW_VAC_RANGE_A[1], "MaNGA DAP EW window high [A]")
    primary.header["EWLOAIR"] = (HALPHA_EW_AIR_RANGE_A[0], "Air EW window low [A]")
    primary.header["EWHIAIR"] = (HALPHA_EW_AIR_RANGE_A[1], "Air EW window high [A]")
    primary.header["LEGHAREF"] = (
        LEGACY_HALPHA_REF_A,
        "Legacy pseudo-EW reference wavelength [A]",
    )
    copy_basic_wcs(primary.header, ha_header)

    hdu_cont = fits.ImageHDU(
        data=cont_mean.astype(np.float64), header=ha_header, name="CONT_HA_MEAN"
    )
    hdu_cont.header["BUNIT"] = CONT_CGS_UNIT
    hdu_cont.header["METHOD"] = "rest_mean"
    hdu_cont.header["REDSHIFT"] = (z, "Galaxy systemic redshift")
    hdu_cont.header["GALVEL"] = (
        galaxy_vel_kms,
        "Systemic velocity from galaxy redshift [km/s]",
    )
    hdu_cont.header["VELHDU"] = vel_extname
    hdu_cont.header["VELUNIT"] = "km/s"
    hdu_cont.header["CONTHDU"] = str(cont_meta["source"])
    hdu_cont.header["CNTBUNIT"] = str(cont_meta["bunit"])
    hdu_cont.header["WAVEAXIS"] = (
        int(cont_meta["spectral_fits_axis"]),
        "Continuum spectral FITS axis",
    )
    hdu_cont.header["WCTYPE"] = str(cont_meta["spectral_ctype"])
    hdu_cont.header["WCUNIT"] = str(cont_meta["spectral_cunit"])
    hdu_cont.header["WINAIR"] = (
        bool(cont_meta["window_is_air"]),
        "True if the selected EW window is in air wavelengths",
    )
    hdu_cont.header["EWLOWRF"] = (
        float(cont_meta["rest_low_a"]),
        "Selected rest-frame EW window low [A]",
    )
    hdu_cont.header["EWHIRF"] = (
        float(cont_meta["rest_high_a"]),
        "Selected rest-frame EW window high [A]",
    )
    hdu_cont.header["EWLOWOB"] = (
        float(cont_meta["rest_low_a"] * (1.0 + z)),
        "Observed-frame EW low [A] at galaxy z",
    )
    hdu_cont.header["EWHIOBS"] = (
        float(cont_meta["rest_high_a"] * (1.0 + z)),
        "Observed-frame EW high [A] at galaxy z",
    )
    hdu_cont.header["EWLOMIN"] = (
        float(cont_meta["observed_low_min_a"]),
        "Min spaxel observed-frame EW low [A]",
    )
    hdu_cont.header["EWLOMAX"] = (
        float(cont_meta["observed_low_max_a"]),
        "Max spaxel observed-frame EW low [A]",
    )
    hdu_cont.header["EWHIMIN"] = (
        float(cont_meta["observed_high_min_a"]),
        "Min spaxel observed-frame EW high [A]",
    )
    hdu_cont.header["EWHIMAX"] = (
        float(cont_meta["observed_high_max_a"]),
        "Max spaxel observed-frame EW high [A]",
    )
    hdu_cont.header["ZSPXMIN"] = (
        float(cont_meta["z_min"]),
        "Minimum spaxel redshift used for EW",
    )
    hdu_cont.header["ZSPXMAX"] = (
        float(cont_meta["z_max"]),
        "Maximum spaxel redshift used for EW",
    )
    hdu_cont.header["NWAVE"] = (int(cont_meta["nwave_total"]), "Total continuum planes")
    hdu_cont.header["NWIN"] = (int(cont_meta["nwave_window"]), "Planes in EW window")
    hdu_cont.header["RESTFSCL"] = (
        True,
        "Continuum density scaled by (1+z_spaxel) to rest-frame f_lambda",
    )
    hdu_cont.header.add_comment(
        "Mean continuum flux density measured from the continuum-only cube over the Halpha EW passband."
    )
    hdu_cont.header.add_comment(
        "The observed-frame EW window is re-centered per spaxel using the galaxy redshift plus HA6562_VEL."
    )

    hdu_cont_npix = fits.ImageHDU(
        data=cont_counts.astype(np.int16), header=ha_header, name="CONT_HA_NPIX"
    )
    hdu_cont_npix.header["BUNIT"] = "count"
    hdu_cont_npix.header.add_comment(
        "Number of finite continuum planes contributing to CONT_HA_MEAN per spaxel"
    )

    hdu_obs_r = fits.ImageHDU(
        data=obs_r_nmgy.astype(np.float64), header=r_header, name="OBS_R_FLUX"
    )
    hdu_obs_r.header["BUNIT"] = r_bunit
    hdu_obs_r.header.add_comment(
        "Observed broad-band R flux map from MAGNITUDE_R_UNCORRECTED converted to nanomaggies"
    )

    hdu_legacy_cont = fits.ImageHDU(
        data=legacy_cont.astype(np.float64), header=ha_header, name="CONT_HA_RBAND"
    )
    hdu_legacy_cont.header["BUNIT"] = CONT_CGS_UNIT
    hdu_legacy_cont.header["METHOD"] = "legacy_r"
    hdu_legacy_cont.header["RUNIT"] = "nanomaggies"
    hdu_legacy_cont.header["LREF"] = (
        LEGACY_HALPHA_REF_A,
        "Legacy single-wavelength continuum reference [A]",
    )
    hdu_legacy_cont.header.add_comment(
        "Legacy broad-band continuum approximation from the R map: f_lambda = f_nu * c / lambda^2."
    )

    hdu_proxy_ew = fits.ImageHDU(
        data=proxy_ewha.astype(np.float64), header=ha_header, name="proxy_EWHa"
    )
    hdu_proxy_ew.header["BUNIT"] = "Angstrom"
    hdu_proxy_ew.header["METHOD"] = "proxy"
    hdu_proxy_ew.header["REDSHIFT"] = (z, "Galaxy systemic redshift")
    hdu_proxy_ew.header["GALVEL"] = (
        galaxy_vel_kms,
        "Systemic velocity from galaxy redshift [km/s]",
    )
    hdu_proxy_ew.header["VELHDU"] = vel_extname
    hdu_proxy_ew.header["HAUNIT"] = "erg s-1 cm-2"
    hdu_proxy_ew.header["CNTUNIT"] = CONT_CGS_UNIT
    hdu_proxy_ew.header["CONTHDU"] = "CONT_HA_MEAN"
    hdu_proxy_ew.header.add_comment(
        "EW_proxy[A] = F_Ha / <f_lambda,cont> using the continuum cube over the Halpha EW passband."
    )
    hdu_proxy_ew.header.add_comment(
        "MaNGA DAP Halpha vacuum passband 6557.6-6571.6 A converted to air with Peck & Reeder (1972)."
    )
    hdu_proxy_ew.header.add_comment(
        "Continuum is shifted to the rest frame using the per-spaxel redshift from galaxy z plus HA6562_VEL and scaled as f_lambda,rest=(1+z_spaxel)*f_lambda,obs."
    )

    hdu_pseudo_ew = fits.ImageHDU(
        data=pseudo_ewha.astype(np.float64), header=ha_header, name="pseudo_EWHa"
    )
    hdu_pseudo_ew.header["BUNIT"] = "Angstrom"
    hdu_pseudo_ew.header["METHOD"] = "legacy_pseudo"
    hdu_pseudo_ew.header["HAUNIT"] = "erg s-1 cm-2"
    hdu_pseudo_ew.header["CNTUNIT"] = CONT_CGS_UNIT
    hdu_pseudo_ew.header["CONTHDU"] = "CONT_HA_RBAND"
    hdu_pseudo_ew.header["LREF"] = (
        LEGACY_HALPHA_REF_A,
        "Legacy single-wavelength continuum reference [A]",
    )
    hdu_pseudo_ew.header.add_comment(
        "Legacy pseudo EW retained for comparison only: EW_pseudo[A] = F_Ha / f_lambda,legacy."
    )
    hdu_pseudo_ew.header.add_comment(
        "f_lambda,legacy is derived from the broad-band R map after converting nanomaggies to f_nu and then to f_lambda at 6562.8 A."
    )

    hdul_out = fits.HDUList(
        [
            primary,
            hdu_cont,
            hdu_cont_npix,
            hdu_obs_r,
            hdu_legacy_cont,
            hdu_proxy_ew,
            hdu_pseudo_ew,
        ]
    )
    hdul_out.writeto(out_path, overwrite=True)

    print("\nSaved:", out_path.resolve())
    print("Halpha source extension       :", ha_extname)
    print("Halpha velocity extension     :", vel_extname)
    print("Legacy R-band source          :", r_source)
    print("Continuum cube source HDU     :", cont_meta["source"])
    print("Continuum spectral axis       :", cont_meta["spectral_ctype"])
    print("Galaxy systemic velocity      :", f"{galaxy_vel_kms:.3f} km/s")
    print("Finite Halpha velocity spaxels:", int(np.sum(np.isfinite(ha_vel_kms))))
    print(
        "Spaxel redshift range         : "
        f"{cont_meta['z_min']:.6f} - {cont_meta['z_max']:.6f}"
    )
    print(
        "Selected rest-frame window    : "
        f"{cont_meta['rest_low_a']:.3f} - {cont_meta['rest_high_a']:.3f} A"
    )
    print(
        "Observed-frame window range   : "
        f"{cont_meta['observed_low_min_a']:.3f} - {cont_meta['observed_high_max_a']:.3f} A"
    )
    print("Continuum planes in window    :", cont_meta["nwave_window"])
    print("Valid CONT_HA_MEAN pixels     :", int(np.sum(cont_counts > 0)))
    print("Valid proxy_EWHa pixels       :", int(np.sum(np.isfinite(proxy_ewha))))
    print("Valid pseudo_EWHa pixels      :", int(np.sum(np.isfinite(pseudo_ewha))))


if __name__ == "__main__":
    main()
