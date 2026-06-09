#!/usr/bin/env python3
"""Build an nGIST-compatible fake spatial mask for a failed AO cube run.

This script exists for the NGC4254 PHANGS/MUSE AO failure mode where nGIST has
already written a mask-table template, but the run cannot proceed to Voronoi
binning. There were two separate problems:

1. A plain 2-D image mask such as ``GALID_mask.fits`` is not the format that
   nGIST's Voronoi module reads after spatial masking has been skipped. At that
   stage it expects extension 1 to be a binary table with a ``MASK`` column, so
   an image HDU fails with ``numpy.ndarray object has no attribute MASK``.
2. Preserving only the original spatial mask is not enough for AO cubes with
   bad variance/noise. Voronoi calls ``voronoi_2d_binning`` with
   ``cube["noise"][MASK == 0]`` and requires every selected noise value to be
   positive and finite. If a locally generated mask marks bad-noise pixels as
   unmasked, Voronoi fails with ``NOISE must be positive and finite``.

The script therefore uses ``GALID_mask_fail.fits`` as the structural template
for the output table and rewrites ``GALID_mask.fits`` with nGIST's mask-column
convention:

* ``0`` means unmasked/usable.
* ``1`` means masked/rejected.
* ``MASK_FILE`` stores the original spatial mask region.
* ``MASK_DEFUNCT`` stores pixels rejected because the cube signal or variance
  would produce invalid Voronoi noise.
* ``MASK`` is the union of ``MASK_FILE`` and ``MASK_DEFUNCT``.
* ``MASK_SNR`` is kept at zero because this script is not applying an S/N
  threshold mask.

When the input cube is available, the script reads ``CONFIG -> GENERAL.INPUT``,
shifts the wavelength grid to rest frame, detects global AO/LGS NaN wavelength
planes using a data-driven NaN-fraction threshold, excludes those planes from
the S/N/noise calculation, and marks pixels with non-finite signal or
non-positive/non-finite median variance as ``MASK_DEFUNCT``. This mirrors the
failure condition that matters for Voronoi without modifying nGIST itself.

Run this on the machine where the cube path exists, for example on Setonix. If
the cube is unavailable, the script refuses to write by default because a
spatial-only mask cannot fix a Voronoi noise assertion. Use ``--spatial-only``
only when intentionally recreating the old spatial-mask-only behavior.

This is a mask-file workaround only. It does not fill NaN planes in the cube,
patch ``MUSE_WFMAON.py``, or edit the pPXF spectral-mask files; those are
separate reader/input-level fixes.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from astropy.io import fits


def _parse_value(value: str):
    value = value.strip().strip("'\"")
    if value.lower() in {"true", "false"}:
        return value.lower() == "true"
    if value.lower() in {"none", "null"}:
        return None
    try:
        if any(char in value for char in ".eE"):
            return float(value)
        return int(value)
    except ValueError:
        return value


def _read_simple_yaml(path: Path) -> dict:
    config: dict[str, dict[str, object]] = {}
    section: str | None = None
    if not path.exists():
        return config

    for raw_line in path.read_text().splitlines():
        line = raw_line.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        if not line.startswith((" ", "\t")) and line.endswith(":"):
            section = line[:-1].strip()
            config[section] = {}
            continue
        if section is not None and ":" in line:
            key, value = line.split(":", 1)
            config[section][key.strip()] = _parse_value(value)

    return config


def _read_spatial_mask(path: Path) -> np.ndarray:
    with fits.open(path) as hdul:
        if len(hdul) > 1 and hdul[1].data is not None:
            data = hdul[1].data
        elif hdul[0].data is not None:
            data = hdul[0].data
        else:
            raise ValueError(f"No mask data found in {path}")

        names = getattr(data, "names", None)
        if names is not None:
            if "MASK_FILE" in names and np.any(np.asarray(data["MASK_FILE"]) != 0):
                data = data["MASK_FILE"]
            elif "MASK" in names:
                data = data["MASK"]

        return (np.asarray(data).reshape(-1) != 0).astype(np.int16)


def _config_float(config: dict, section: str, key: str, default: float) -> float:
    value = config.get(section, {}).get(key, default)
    if value is None:
        return default
    return float(value)


def _cube_noise_mask(
    config: dict,
    cube_path: Path,
    expected_size: int,
    *,
    global_nan_threshold: float,
    chunk_rows: int,
    wave_chunk: int,
) -> np.ndarray:
    redshift = _config_float(config, "GENERAL", "REDSHIFT", 0.0)
    lmin_tot = _config_float(config, "READ_DATA", "LMIN_TOT", 4800.0)
    lmax_tot = _config_float(config, "READ_DATA", "LMAX_TOT", 8900.0)
    lmin_snr = _config_float(config, "READ_DATA", "LMIN_SNR", lmin_tot)
    lmax_snr = _config_float(config, "READ_DATA", "LMAX_SNR", lmax_tot)

    with fits.open(cube_path, memmap=True) as hdul:
        data = hdul[1].data
        shape = data.shape
        if len(shape) != 3:
            raise ValueError(f"Expected 3D cube in HDU 1: {cube_path}")
        nwave, ny, nx = shape
        if ny * nx != expected_size:
            raise ValueError(
                f"Cube spatial size mismatch: {cube_path} has {ny * nx} pixels, "
                f"but the mask expects {expected_size}"
            )

        hdr = hdul[1].header
        cdelt = hdr.get("CD3_3", hdr.get("CDELT3"))
        if cdelt is None:
            raise ValueError(f"Cube header has no CD3_3/CDELT3 wavelength step: {cube_path}")
        wave = hdr["CRVAL3"] + np.arange(nwave) * cdelt
        wave = wave / (1.0 + redshift)

        idx_tot = np.where(np.logical_and(wave >= lmin_tot, wave <= lmax_tot))[0]
        global_nan = np.zeros(nwave, dtype=bool)
        for w0 in range(0, len(idx_tot), wave_chunk):
            idx_wave = idx_tot[w0 : w0 + wave_chunk]
            nan_fraction = np.mean(np.isnan(data[idx_wave, :, :]), axis=(1, 2))
            global_nan[idx_wave] = nan_fraction > global_nan_threshold

        idx_global_nan = np.where(global_nan)[0]
        if len(idx_global_nan) > 0:
            print(
                "Detected "
                f"{len(idx_global_nan)} global-NaN wavelength planes "
                f"(threshold>{global_nan_threshold:.2f})."
            )
            print(
                "Global-NaN rest-frame range: "
                f"{wave[idx_global_nan[0]]:.2f}A - {wave[idx_global_nan[-1]]:.2f}A"
            )
        else:
            print(
                "Detected 0 global-NaN wavelength planes; "
                "no AO/LGS planes were excluded from the noise calculation."
            )

        idx_snr = np.where(
            np.logical_and.reduce(
                (
                    wave >= lmin_tot,
                    wave <= lmax_tot,
                    wave >= lmin_snr,
                    wave <= lmax_snr,
                    ~global_nan,
                )
            )
        )[0]
        if len(idx_snr) == 0:
            raise ValueError("No wavelengths remain for the S/N noise check")

        bad = np.zeros(expected_size, dtype=bool)
        for y0 in range(0, ny, chunk_rows):
            y1 = min(y0 + chunk_rows, ny)
            signal = np.nanmedian(data[idx_snr, y0:y1, :], axis=0).reshape(-1)

            if len(hdul) >= 3 and hdul[2].data is not None:
                variance = np.nanmedian(hdul[2].data[idx_snr, y0:y1, :], axis=0).reshape(-1)
                noise_bad = ~np.isfinite(variance) | (variance <= 0.0)
            else:
                noise_bad = np.zeros_like(signal, dtype=bool)

            start = y0 * nx
            stop = y1 * nx
            bad[start:stop] = (~np.isfinite(signal)) | noise_bad

    return bad.astype(np.int16)


def _write_mask_in_template(
    mask_path: Path,
    template_path: Path,
    spatial_mask: np.ndarray,
    defunct_mask: np.ndarray,
) -> tuple[int, int, int, int]:
    with fits.open(template_path) as template:
        output = fits.HDUList([hdu.copy() for hdu in template])

    try:
        table = output[1].data
        names = table.names
    except (IndexError, AttributeError) as exc:
        raise ValueError(f"Template has no table extension with a MASK column: {template_path}") from exc

    if "MASK" not in names:
        raise ValueError(f"Template table has no MASK column: {template_path}")

    target = table["MASK"]
    if target.size != spatial_mask.size:
        raise ValueError(
            f"Mask size mismatch: {mask_path} has {spatial_mask.size} values, "
            f"but {template_path} expects {target.size}"
        )
    if target.size != defunct_mask.size:
        raise ValueError(
            f"Defunct mask size mismatch: got {defunct_mask.size} values, "
            f"but {template_path} expects {target.size}"
        )

    spatial_mask = spatial_mask.reshape(target.shape)
    defunct_mask = defunct_mask.reshape(target.shape)
    combined_mask = np.logical_or(spatial_mask != 0, defunct_mask != 0).astype(target.dtype)

    table["MASK"] = combined_mask
    if "MASK_DEFUNCT" in names:
        table["MASK_DEFUNCT"] = defunct_mask.astype(table["MASK_DEFUNCT"].dtype, copy=False)
    if "MASK_SNR" in names:
        table["MASK_SNR"] = np.zeros(target.shape, dtype=table["MASK_SNR"].dtype)
    if "MASK_FILE" in names:
        table["MASK_FILE"] = spatial_mask.astype(table["MASK_FILE"].dtype, copy=False)
    output.writeto(mask_path, overwrite=True)

    masked = int(np.sum(combined_mask != 0))
    unmasked = int(combined_mask.size - masked)
    spatial = int(np.sum(spatial_mask != 0))
    defunct = int(np.sum(defunct_mask != 0))
    return masked, unmasked, spatial, defunct


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Rewrite GALID_mask.fits using GALID_mask_fail.fits format, preserving "
            "the spatial mask and masking pixels with invalid Voronoi noise when "
            "the input cube is available."
        )
    )
    parser.add_argument("galid", help="Galaxy ID, for example NGC4254")
    parser.add_argument("--config", type=Path, help="Path to CONFIG. Defaults to GALID/CONFIG.")
    parser.add_argument("--cube", type=Path, help="Input datacube path. Defaults to GENERAL.INPUT in CONFIG.")
    parser.add_argument(
        "--spatial-only",
        action="store_true",
        help="Allow writing without checking the input cube. This cannot fix Voronoi noise failures.",
    )
    parser.add_argument("--chunk-rows", type=int, default=16, help="Spatial rows per FITS processing chunk.")
    parser.add_argument(
        "--wave-chunk",
        type=int,
        default=16,
        help="Wavelength planes per FITS processing chunk for global-NaN detection.",
    )
    parser.add_argument(
        "--global-nan-threshold",
        type=float,
        default=0.90,
        help="Mask wavelength planes with NaN fraction above this value from the noise calculation.",
    )
    args = parser.parse_args()

    base_dir = Path(__file__).resolve().parent
    galaxy_dir = base_dir / args.galid
    mask_path = galaxy_dir / f"{args.galid}_mask.fits"
    template_path = galaxy_dir / f"{args.galid}_mask_fail.fits"
    config_path = args.config or galaxy_dir / "CONFIG"

    config = _read_simple_yaml(config_path)
    spatial_mask = _read_spatial_mask(mask_path)
    defunct_mask = np.zeros_like(spatial_mask, dtype=np.int16)

    cube_value = args.cube or config.get("GENERAL", {}).get("INPUT")
    cube_path = Path(cube_value) if cube_value else None
    if cube_path is not None and cube_path.exists():
        defunct_mask = _cube_noise_mask(
            config,
            cube_path,
            spatial_mask.size,
            global_nan_threshold=args.global_nan_threshold,
            chunk_rows=args.chunk_rows,
            wave_chunk=args.wave_chunk,
        )
    else:
        if not args.spatial_only:
            raise FileNotFoundError(
                "Input cube is unavailable, so a Voronoi-safe fake mask cannot be generated. "
                f"Run this script where the CONFIG GENERAL.INPUT cube exists, or pass --cube explicitly. "
                f"Cube path: {cube_path}. Use --spatial-only only if you intentionally want the old "
                "spatial-mask-only behavior."
            )
        print(f"Warning: input cube unavailable; writing spatial mask only: {cube_path}")

    masked, unmasked, spatial, defunct = _write_mask_in_template(
        mask_path,
        template_path,
        spatial_mask,
        defunct_mask,
    )
    print(f"Wrote {mask_path}")
    print(f"MASK masked={masked} unmasked={unmasked}")
    print(f"MASK_FILE spatial_masked={spatial}")
    print(f"MASK_DEFUNCT invalid_noise_masked={defunct}")


if __name__ == "__main__":
    main()
