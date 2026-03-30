#!/usr/bin/env python
"""
pseudo_EWHa.py - Build proxy EW(Halpha) maps for a MAUVE galaxy.

Changes (2026-03-30)
--------------------
* Initial implementation.
* Inputs:
  - {gal}_SPATIAL_BINNING_maps_extended.fits
  - {gal}_gas_BIN_maps_extended.fits
* Output:
  - {gal}_pseudo_EW_maps.fits
        containing 3 HDUs total:
            PRIMARY(OBS_HA6562_FLUX), OBS_R_FLUX, pseudo_EWHa.

Notes:
* This script computes a proxy EW(Halpha) from broadband r-band continuum.
* The computation is:
        EW_proxy = F_Ha / f_lambda(6562.8 A)
    with:
        f_nu = R_nanomaggy * 3.631e-29  [erg s-1 cm-2 Hz-1]
        f_lambda = f_nu * c / lambda^2  [erg s-1 cm-2 A-1]
* No redshift correction is applied.
"""

import argparse
from pathlib import Path
from typing import Any, cast

import numpy as np
from astropy.io import fits


C_AA_PER_S = 2.99792458e18
NMAGGY_TO_FNU = 3.631e-29
HALPHA_REF_A = 6562.8


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create pseudo EW(Halpha) maps from observed Halpha and R-band flux maps."
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
        "--bin-file",
        default=None,
        help="Optional explicit path to {gal}_SPATIAL_BINNING_maps_extended.fits",
    )
    parser.add_argument(
        "--gas-file",
        default=None,
        help="Optional explicit path to {gal}_gas_BIN_maps_extended.fits",
    )
    parser.add_argument(
        "--out",
        default=None,
        help="Optional output path (default: {gal}_pseudo_EW_maps.fits in current directory)",
    )
    return parser.parse_args()


def magnitude_to_nanomaggies(magnitude: np.ndarray) -> np.ndarray:
    return 10 ** ((22.5 - magnitude) / 2.5)


def halpha_flux_to_cgs(line_flux: np.ndarray, bunit: str) -> np.ndarray:
    """Convert Halpha map to erg s^-1 cm^-2 when BUNIT encodes a 1e-20 scale."""
    unit_norm = str(bunit).replace(" ", "").lower()
    scale = 1.0
    if ("1e-20" in unit_norm) or ("10^-20" in unit_norm) or ("10**-20" in unit_norm):
        scale = 1e-20
    return line_flux * scale


def nanomaggies_to_flambda(nmgy: np.ndarray, lam_a: float) -> np.ndarray:
    """Convert nanomaggies (f_nu) to f_lambda at wavelength lam_a (Angstrom)."""
    f_nu = nmgy * NMAGGY_TO_FNU
    return f_nu * C_AA_PER_S / (lam_a ** 2)


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


def read_observed_r_flux(bin_path: Path) -> tuple[np.ndarray, fits.Header, str, str]:
    with fits.open(bin_path) as hdul:
        available = {
            str(getattr(hdu, "name", "")).upper()
            for hdu in hdul
            if getattr(hdu, "name", "")
        }

        if "MAGNITUDE_R_UNCORRECTED" in available:
            hdu = cast(Any, hdul["magnitude_r_uncorrected"])
            mag = np.asarray(hdu.data, dtype=np.float64)
            data = magnitude_to_nanomaggies(mag)
            header = hdu.header.copy()
            return data, header, "magnitude_r_uncorrected", "nanomaggies"

    raise KeyError(
        f"Could not find R-band flux in {bin_path.name}. "
        "Expected magnitude_r_uncorrected."
    )


def main() -> None:
    args = parse_args()
    gal = args.galaxy.upper()
    root = Path(args.root).expanduser().resolve()

    bin_path = (
        Path(args.bin_file)
        if args.bin_file
        else root / f"{gal}_SPATIAL_BINNING_maps_extended.fits"
    )
    gas_path = (
        Path(args.gas_file)
        if args.gas_file
        else root / f"{gal}_gas_BIN_maps_extended.fits"
    )
    out_path = Path(args.out) if args.out else Path(f"{gal}_pseudo_EW_maps.fits")

    print("\n=== pseudo_EWHa inputs/outputs ===")
    print("Galaxy      :", gal)
    print("Binning FITS:", bin_path)
    print("Gas FITS    :", gas_path)
    print("Output FITS :", out_path)

    if not bin_path.exists():
        raise FileNotFoundError(f"Input file not found: {bin_path}")
    if not gas_path.exists():
        raise FileNotFoundError(f"Input file not found: {gas_path}")

    obs_ha, ha_header, ha_extname, ha_bunit = read_observed_halpha(gas_path)
    obs_r, r_header, r_source, r_bunit = read_observed_r_flux(bin_path)

    if obs_ha.shape != obs_r.shape:
        raise ValueError(
            "Map shape mismatch: "
            f"Halpha shape {obs_ha.shape} vs R-band shape {obs_r.shape}."
        )

    obs_ha_cgs = halpha_flux_to_cgs(obs_ha, ha_bunit)
    cont_flam = nanomaggies_to_flambda(obs_r, HALPHA_REF_A)

    pseudo_ewha = np.full(obs_ha.shape, np.nan, dtype=np.float64)
    valid = np.isfinite(obs_ha_cgs) & np.isfinite(cont_flam) & (cont_flam > 0)
    pseudo_ewha[valid] = obs_ha_cgs[valid] / cont_flam[valid]

    primary = fits.PrimaryHDU(data=obs_ha.astype(np.float64))
    primary.header["EXTNAME"] = "OBS_HA6562_FLUX"
    primary.header["BUNIT"] = ha_bunit
    primary.header["GALAXY"] = gal
    primary.header["BINFILE"] = bin_path.name
    primary.header["GASFILE"] = gas_path.name
    primary.header["HA_SRC"] = ha_extname
    primary.header["R_SRC"] = r_source
    primary.header["EWLREF"] = (HALPHA_REF_A, "Lambda used for continuum conversion [A]")

    # Preserve basic WCS metadata from the Halpha extension when present.
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
        if key in ha_header:
            primary.header[key] = ha_header[key]

    hdu_obs_r = fits.ImageHDU(
        data=obs_r.astype(np.float64), header=r_header, name="OBS_R_FLUX"
    )
    hdu_obs_r.header["BUNIT"] = r_bunit
    hdu_obs_r.header["COMMENT"] = "Observed R-band flux map from spatial-binning extended FITS"

    hdu_pseudo_ew = fits.ImageHDU(
        data=pseudo_ewha.astype(np.float64), header=ha_header, name="pseudo_EWHa"
    )
    hdu_pseudo_ew.header["BUNIT"] = "Angstrom"
    hdu_pseudo_ew.header["METHOD"] = "proxy"
    hdu_pseudo_ew.header["HAUNIT"] = "erg s-1 cm-2"
    hdu_pseudo_ew.header["RUNIT"] = "nanomaggy"
    hdu_pseudo_ew.header.add_comment(
        "EW_proxy[A] = F_Ha / f_lambda(6562.8A); no redshift correction applied"
    )
    hdu_pseudo_ew.header.add_comment(
        "f_nu = R_nanomaggy * 3.631e-29 [erg s-1 cm-2 Hz-1]"
    )
    hdu_pseudo_ew.header.add_comment(
        "f_lambda = f_nu * c / lambda^2, c=2.99792458e18 [A s-1]"
    )

    hdul_out = fits.HDUList([primary, hdu_obs_r, hdu_pseudo_ew])
    hdul_out.writeto(out_path, overwrite=True)

    print("\nSaved:", out_path.resolve())
    print("Halpha source extension:", ha_extname)
    print("R-band source extension:", r_source)
    print("Valid pseudo_EWHa pixels:", int(np.sum(np.isfinite(pseudo_ewha))))


if __name__ == "__main__":
    main()

