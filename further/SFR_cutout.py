#!/usr/bin/env python
"""Extract SFR surface-density HDUs from a MAUVE further gas-map FITS file."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from astropy.io import fits


SFR_HDU_PREFIX = "LOGSFR_SURFACE_DENSITY"
GALAXY_ID_PATTERN = re.compile(r"^(?:IC\d+|NGC\d+(?:_\d+)?)$", re.IGNORECASE)


def source_path_for(galaxy: str, root: Path, product_subdir: Path) -> Path:
    """Return the uncompressed or compressed further gas-map input path."""
    galaxy_dir = root / product_subdir / galaxy
    base = galaxy_dir / f"{galaxy}_gas_bin_maps_further.fits"
    for candidate in (base, Path(f"{base}.gz")):
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(
        f"Could not find {base} or {base}.gz for galaxy {galaxy}."
    )


def extract_sfr_hdus(source: Path, output: Path) -> list[str]:
    """Copy PRIMARY, BIN_ID, and every LOGSFR_SURFACE_DENSITY* HDU."""
    with fits.open(source, memmap=False) as source_hdul:
        binid_hdu = source_hdul["BIN_ID"].copy()
        selected = [
            hdu.copy()
            for hdu in source_hdul[1:]
            if hdu.name.upper().startswith(SFR_HDU_PREFIX)
        ]
        if not selected:
            raise ValueError(
                f"No {SFR_HDU_PREFIX}* HDUs were found in {source}."
            )

        output_hdul = fits.HDUList(
            [source_hdul[0].copy(), binid_hdu, *selected]
        )
        output_hdul.writeto(output, overwrite=True)

    return [hdu.name for hdu in selected]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract all LOGSFR_SURFACE_DENSITY* HDUs into a compact SFR FITS file."
        )
    )
    parser.add_argument("-g", "--galaxy", required=True, help="Galaxy ID")
    parser.add_argument(
        "--root",
        default=".",
        help="Root directory containing the selected product subdirectory",
    )
    parser.add_argument(
        "--product-subdir",
        default="v3tk_v7.6.8",
        help="Product subdirectory under --root (default: v3tk_v7.6.8)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    galaxy = args.galaxy.upper()
    if not GALAXY_ID_PATTERN.fullmatch(galaxy):
        raise ValueError(
            "Invalid galaxy ID; expected IC<digits>, NGC<digits>, "
            "or NGC<digits>_<digits>."
        )

    root = Path(args.root).expanduser().resolve()
    product_subdir = Path(args.product_subdir)
    source = source_path_for(galaxy, root, product_subdir)
    output = source.parent / f"{galaxy}_SFR_maps_further.fits"
    names = extract_sfr_hdus(source, output)

    print(f"Input : {source}")
    print(f"Output: {output}")
    print(f"Copied {len(names)} SFR HDUs: {', '.join(names)}")


if __name__ == "__main__":
    main()
