#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from astropy.io import fits


def _read_mask_values(path: Path) -> np.ndarray:
    with fits.open(path) as hdul:
        if len(hdul) > 1 and hdul[1].data is not None:
            data = hdul[1].data
        elif hdul[0].data is not None:
            data = hdul[0].data
        else:
            raise ValueError(f"No mask data found in {path}")

        names = getattr(data, "names", None)
        if names is not None and "MASK" in names:
            data = data["MASK"]

        return np.asarray(data).reshape(-1)


def _write_mask_in_template(mask_path: Path, template_path: Path, values: np.ndarray) -> None:
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
    if target.size != values.size:
        raise ValueError(
            f"Mask size mismatch: {mask_path} has {values.size} values, "
            f"but {template_path} expects {target.size}"
        )

    table["MASK"] = values.reshape(target.shape).astype(target.dtype, copy=False)
    output.writeto(mask_path, overwrite=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Rewrite GALID_mask.fits using GALID_mask_fail.fits format while preserving mask values."
    )
    parser.add_argument("galid", help="Galaxy ID, for example NGC4254")
    args = parser.parse_args()

    base_dir = Path(__file__).resolve().parent
    galaxy_dir = base_dir / args.galid
    mask_path = galaxy_dir / f"{args.galid}_mask.fits"
    template_path = galaxy_dir / f"{args.galid}_mask_fail.fits"

    values = _read_mask_values(mask_path)
    _write_mask_in_template(mask_path, template_path, values)


if __name__ == "__main__":
    main()
