#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
python="${PYTHON:-/opt/miniconda3/envs/ICRAR/bin/python}"
"$python" -B - "$script_dir" "$@" <<'PY'
import argparse
import csv
import math
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_area, proj_plane_pixel_scales
from PIL import Image, ImageDraw, ImageFont
from reproject import reproject_adaptive

script_dir = Path(sys.argv[1])
source = script_dir.parent / "v3tk_to_VRI"
sys.path.insert(0, str(source))
from v3tk_get_legacy import _compute_footprint_size_arcsec, _legacy_fits_cutout_url

SUFFIX = "_observed_VRI.png"
LEGACY_SCALE = 0.262
with (script_dir.parent / "further" / "simbad_object_crosswalk.csv").open(newline="") as handle:
    catalog = {row["component"]: row for row in csv.DictReader(handle)}


def legacy_background(wcs, shape, directory):
    height, width = shape
    center = wcs.pixel_to_world((width - 1) / 2, (height - 1) / 2).icrs
    width_arcsec, height_arcsec = _compute_footprint_size_arcsec(wcs, shape)
    # Cover the rotated footprint and provide room for the adaptive kernel.
    url = _legacy_fits_cutout_url(
        ra_deg=center.ra.deg, dec_deg=center.dec.deg, layer="ls-dr9",
        pixscale=LEGACY_SCALE, width_px=math.ceil(width_arcsec / LEGACY_SCALE) + 16,
        height_px=math.ceil(height_arcsec / LEGACY_SCALE) + 16, bands="grz",
    )
    path = directory / "legacy.fits"
    subprocess.run([
        "curl", "--fail", "--silent", "--show-error", "--location",
        "--retry", "4", "--connect-timeout", "15", "--max-time", "120",
        "--output", str(path), url,
    ], check=True)
    with fits.open(path) as hdus:
        hdu = hdus[0]
        if hdu.header.get("BANDS") != "grz" or hdu.data.shape[0] != 3:
            raise ValueError("Legacy cutout does not contain the requested g, r, z bands")
        input_wcs = WCS(hdu.header, naxis=2)
        flux, footprint = reproject_adaptive(
            (hdu.data, input_wcs), wcs, shape_out=(3, *shape), conserve_flux=True,
        )
    if not np.all((footprint > 0) & np.isfinite(flux)):
        raise ValueError("Legacy reprojection did not cover the full output canvas")

    # Convert flux per output pixel to surface brightness before stretching.
    surface_brightness = flux / (proj_plane_pixel_area(wcs) * 3600**2)
    # Use the Legacy Survey grz asinh stretch, in equivalent native-pixel units.
    # https://github.com/legacysurvey/legacypipe/blob/main/py/legacypipe/survey.py
    rgb = np.moveaxis(surface_brightness[::-1], 0, -1) * LEGACY_SCALE**2
    scaled = np.maximum(0, rgb * [2.2, 3.4, 6.0] + 0.03)
    intensity = scaled.mean(axis=-1, keepdims=True)
    stretch = np.divide(
        np.arcsinh(20 * intensity), np.sqrt(20) * intensity,
        out=np.zeros_like(intensity), where=intensity > 0,
    )
    return np.rint(np.clip(scaled * stretch, 0, 1) * 255).astype(np.uint8)


def centered_grid(wcs, a, b, galaxy):
    row = catalog[galaxy]
    sky = SkyCoord(float(row["ra_deg"]), float(row["dec_deg"]), unit="deg")
    x, y = map(float, wcs.world_to_pixel(sky))
    if not (0 <= x < a and 0 <= y < b):
        raise ValueError(f"{galaxy}: catalog center falls outside the observed grid")

    # Pixel centers span 0 to a-1 and 0 to b-1. The output midpoint may
    # differ from the projected nucleus by half a pixel with an integer shift.
    k = math.ceil(max(a / 4, b / 3, max(x, a - 1 - x) / 2, 2 * max(y, b - 1 - y) / 3))
    width, height = 4 * k, 3 * k
    left = min(max(round((width - 1) / 2 - x), 0), width - a)
    bottom = min(max(round((height - 1) / 2 - y), 0), height - b)
    target = wcs.deepcopy()
    target.wcs.crpix += [left, bottom]
    target.pixel_shape = (width, height)
    return target, (height, width), left, bottom


def render(png, output):
    cube = png.name.removesuffix(SUFFIX)
    [path] = source.glob(f"{cube}_*_VRI.fits")
    with Image.open(png) as image, fits.open(path) as hdus:
        a, b = image.size
        observed = np.asarray(image.convert("RGB"))[::-1]
        hdu = hdus["R_FLUX"]
        valid = np.isfinite(hdu.data) & (hdu.data > 0)
        wcs = WCS(hdu.header).celestial
    if valid.shape != (b, a):
        raise ValueError("Observed PNG and FITS dimensions differ")
    if not np.allclose(proj_plane_pixel_scales(wcs) * 3600, 0.2, atol=1e-5, rtol=0):
        raise ValueError("The observed grid does not have 0.2 arcsec pixels")

    labels = ("NGC4567", "NGC4568") if cube == "NGC4567_8" else (cube,)
    results = []
    with TemporaryDirectory(dir=output) as temporary:
        directory = Path(temporary)
        for galaxy in labels:
            target, (height, width), left, bottom = centered_grid(wcs, a, b, galaxy)
            rgb = legacy_background(target, (height, width), directory)
            rgb[bottom:bottom + b, left:left + a][valid] = observed[valid]
            composite = Image.fromarray(rgb[::-1])
            ImageDraw.Draw(composite).text(
                (12, 12), galaxy, font=ImageFont.load_default(size=24),
                fill="yellow", stroke_width=2, stroke_fill="black",
            )
            name = f"{galaxy}_MUSE_G4.jpg"
            composite.save(directory / name, quality=98, subsampling=0)
            results.append(f"{galaxy}: {a}x{b} -> {width}x{height}")
        for galaxy in labels:
            name = f"{galaxy}_MUSE_G4.jpg"
            (directory / name).replace(output / name)
    return "\n".join(results)


parser = argparse.ArgumentParser(description="Build native-scale 4:3 MUSE + Legacy images")
parser.add_argument("cubes", nargs="*", help="Observed cube IDs to process (default: all)")
parser.add_argument("--output-dir", type=Path, default=script_dir)
args = parser.parse_args(sys.argv[2:])
images = sorted(source.glob(f"*{SUFFIX}"))
if args.cubes:
    available = {p.name.removesuffix(SUFFIX): p for p in images}
    missing = set(args.cubes) - available.keys()
    if missing:
        parser.error(f"No observed PNG for: {', '.join(sorted(missing))}")
    images = [available[cube] for cube in dict.fromkeys(args.cubes)]
if not images:
    parser.error(f"No observed images in {source}")
args.output_dir.mkdir(parents=True, exist_ok=True)

failed = False
with ThreadPoolExecutor(max_workers=2) as pool:
    jobs = {pool.submit(render, png, args.output_dir): png for png in images}
    for job in as_completed(jobs):
        try:
            print(job.result(), flush=True)
        except Exception as exc:
            print(f"{jobs[job].name}: {exc}", file=sys.stderr, flush=True)
            failed = True
sys.exit(int(failed))
PY
