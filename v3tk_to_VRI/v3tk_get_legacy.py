#!/usr/bin/env python

from __future__ import annotations

import argparse
import concurrent.futures as cf
import math
import os
import pathlib
import time
from dataclasses import dataclass
import sys
import threading
import urllib.error


def _require_deps():
	try:
		import numpy as np  # noqa: F401
		from astropy.io import fits  # noqa: F401
		from astropy.wcs import WCS  # noqa: F401
		from astropy.coordinates import SkyCoord  # noqa: F401
		import astropy.units as u  # noqa: F401
		from reproject import reproject_adaptive  # noqa: F401
		from PIL import Image  # noqa: F401
		import urllib.request  # noqa: F401
		import urllib.parse  # noqa: F401
	except Exception as exc:  # pragma: no cover
		import sys as _sys
		raise RuntimeError(
			"Missing required dependency.\n\n"
			"Required: numpy astropy reproject pillow\n\n"
			f"Interpreter used: {_sys.executable}\n"
			f"Python version: {_sys.version.splitlines()[0]}\n\n"
			"If you already installed these in your conda env, you are likely running with a different Python.\n"
			"Try:\n"
			"  conda activate ICRAR\n"
			"  python v3tk_get_legacy.py --dry-run\n\n"
			"Install deps if needed:\n"
			"  conda install -c conda-forge numpy astropy reproject pillow\n"
			"  # or: python -m pip install numpy astropy reproject pillow\n"
		) from exc


@dataclass(frozen=True)
class Job:
	input_path: pathlib.Path
	galaxy_id: str
	observed_png: pathlib.Path
	legacy_cutout_jpg: pathlib.Path
	legacy_reprojected_jpg: pathlib.Path


def _parse_args(argv: list[str]) -> argparse.Namespace:
	p = argparse.ArgumentParser(
		prog="v3tk_get_legacy.py",
		description=(
			"For each '*_v3tk_VRI.fits' file, compute the MUSE WCS footprint, download a Legacy Survey cutout "
			"at native pixel scale (pixscale=0.262 arcsec/pix), reproject it to the MUSE grid using "
			"reproject.reproject_adaptive(conserve_flux=True), and save 'XXX_legacy_reprojected.jpg'."
		),
	)
	p.add_argument(
		"--input-dir",
		type=pathlib.Path,
		default=pathlib.Path("."),
		help="Directory containing *_v3tk_VRI.fits files (default: current directory)",
	)
	p.add_argument(
		"--pattern",
		default="*_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk_VRI.fits",
		help="Glob pattern to match input FITS files",
	)
	p.add_argument(
		"--layer",
		default="ls-dr9",
		help="Legacy viewer layer (default: ls-dr9)",
	)
	p.add_argument(
		"--pixscale",
		type=float,
		default=0.262,
		help="Legacy cutout pixel scale in arcsec/pixel (default: 0.262)",
	)
	p.add_argument(
		"--margin-factor",
		type=float,
		default=1.02,
		help="Multiply the MUSE footprint size by this factor when requesting the Legacy cutout (default: 1.02)",
	)
	p.add_argument(
		"--timeout",
		type=float,
		default=120.0,
		help="Download timeout in seconds (default: 120)",
	)
	p.add_argument(
		"--workers",
		type=int,
		default=None,
		help="Parallel workers across galaxies (default: min(2, os.cpu_count()))",
	)
	p.add_argument(
		"--max-concurrent-downloads",
		type=int,
		default=2,
		help="Max concurrent HTTP downloads (default: 2). Reduce if you hit HTTP 429.",
	)
	p.add_argument(
		"--overwrite",
		dest="overwrite",
		action="store_true",
		help="Overwrite existing outputs (default)",
	)
	p.add_argument(
		"--no-overwrite",
		dest="overwrite",
		action="store_false",
		help="Do not overwrite existing outputs",
	)
	p.set_defaults(overwrite=True)
	p.add_argument(
		"--keep-cutout",
		action="store_true",
		help="Keep the downloaded cutout JPEG (default: delete after successful reprojection)",
	)
	p.add_argument(
		"--dry-run",
		action="store_true",
		help="List planned work but do not download/reproject",
	)
	p.add_argument(
		"--quiet",
		action="store_true",
		help="Only print errors",
	)
	return p.parse_args(argv)


def _galaxy_id_from_filename(path: pathlib.Path) -> str:
	name = path.name
	marker = "_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk_VRI.fits"
	if name.endswith(marker):
		return name[: -len(marker)]
	if name.endswith("_v3tk_VRI.fits"):
		return name[: -len("_v3tk_VRI.fits")]
	return name.split("_", 1)[0]


def _select_flux_hdu(hdul):
	"""Return the HDU that contains the flux map in nanomaggy."""
	# Prefer explicit naming.
	for key in ("R_FLUX", "FLUX", "RFLUX"):
		try:
			return hdul[key]
		except Exception:
			pass

	# Heuristic: first 2D image HDU whose unit looks like nanomaggy.
	for hdu in hdul:
		data = getattr(hdu, "data", None)
		if data is None or getattr(data, "ndim", None) != 2:
			continue
		bunit = (hdu.header.get("BUNIT") or "").lower()
		extname = (hdu.header.get("EXTNAME") or "").upper()
		name = (getattr(hdu, "name", "") or "").upper()
		if "nanomag" in bunit or "nmgy" in bunit or extname == "R_FLUX" or name == "R_FLUX":
			return hdu

	# Fallback: first 2D image HDU.
	for hdu in hdul:
		data = getattr(hdu, "data", None)
		if data is None or getattr(data, "ndim", None) != 2:
			continue
		return hdu

	raise ValueError("Could not find a 2D flux image HDU")


def _muse_target_wcs_and_shape(fits_path: pathlib.Path):
	from astropy.io import fits
	from astropy.wcs import WCS

	with fits.open(fits_path, memmap=True) as hdul:
		h_flux = _select_flux_hdu(hdul)
		data = h_flux.data
		if data is None or getattr(data, "ndim", None) != 2:
			raise ValueError("Selected flux HDU does not contain a 2D image")
		w = WCS(h_flux.header).celestial
		ny, nx = data.shape
		return w, (ny, nx)


def _compute_footprint_size_arcsec(target_wcs, shape: tuple[int, int]) -> tuple[float, float]:
	"""Return (width_arcsec, height_arcsec) of the target footprint around CRVAL.

	Uses spherical offsets from the WCS reference position (CRVAL) to the four corners.
	This avoids RA wrap issues and works with rotated WCS.
	"""
	import astropy.units as u
	from astropy.coordinates import SkyCoord

	# Corners in pixel coordinates (0-based for WCS.pixel_to_world).
	ny, nx = shape
	pix = [(0, 0), (0, ny - 1), (nx - 1, 0), (nx - 1, ny - 1)]
	world = [target_wcs.pixel_to_world(x, y) for (x, y) in pix]
	coords = SkyCoord(world)

	# Center from the *image center pixel* (more robust than CRVAL/CRPIX for asymmetric WCS).
	ny, nx = shape
	xc = (nx - 1) / 2.0
	yc = (ny - 1) / 2.0
	center = target_wcs.pixel_to_world(xc, yc)
	center = SkyCoord(center).transform_to(coords.frame)

	# Convert corners into tangent-plane offsets from center.
	dlon, dlat = center.spherical_offsets_to(coords)  # SkyCoord -> Quantity
	max_lon = float(max(abs(dlon.to_value(u.arcsec))))
	max_lat = float(max(abs(dlat.to_value(u.arcsec))))
	width = 2.0 * max_lon
	height = 2.0 * max_lat
	return width, height


def _legacy_cutout_url(
	*,
	ra_deg: float,
	dec_deg: float,
	layer: str,
	pixscale: float,
	width_px: int,
	height_px: int,
) -> str:
	import urllib.parse

	params = {
		"ra": f"{ra_deg:.10f}",
		"dec": f"{dec_deg:.10f}",
		"layer": layer,
		"pixscale": f"{pixscale:.5f}",
		"width": str(int(width_px)),
		"height": str(int(height_px)),
	}
	return "https://www.legacysurvey.org/viewer/jpeg-cutout?" + urllib.parse.urlencode(params)


def _legacy_fits_cutout_url(
	*,
	ra_deg: float,
	dec_deg: float,
	layer: str,
	pixscale: float,
	width_px: int,
	height_px: int,
	bands: str = "r",
) -> str:
	import urllib.parse

	params = {
		"ra": f"{ra_deg:.10f}",
		"dec": f"{dec_deg:.10f}",
		"layer": layer,
		"pixscale": f"{pixscale:.5f}",
		"width": str(int(width_px)),
		"height": str(int(height_px)),
		"bands": bands,
	}
	return "https://www.legacysurvey.org/viewer/fits-cutout?" + urllib.parse.urlencode(params)


def _download_file(url: str, out_path: pathlib.Path, timeout_s: float) -> None:
	import urllib.request

	out_path.parent.mkdir(parents=True, exist_ok=True)
	req = urllib.request.Request(url, headers={"User-Agent": "ICRAR-MAUVE/legacy-cutout"})

	# Throttle concurrent downloads to avoid server rate-limits.
	sem = globals().get("_DOWNLOAD_SEM")
	if sem is None:
		# Fallback when run outside main.
		sem = threading.Semaphore(2)

	last_exc = None
	for attempt in range(20):
		with sem:
			try:
				with urllib.request.urlopen(req, timeout=timeout_s) as resp:
					data = resp.read()
				out_path.write_bytes(data)
				return
			except urllib.error.HTTPError as exc:
				# 429/503: exponential backoff with jitter
				last_exc = exc
				code = getattr(exc, "code", None)
				if code in (429, 503):
					# Random jitter to prevent thundering herd
					import random
					sleep_s = min(120.0, 5.0 + (1.5 ** attempt)) + random.uniform(0.0, 5.0)
					time.sleep(sleep_s)
					continue
				# Anything else: fail immediately with a picklable/printable error.
				raise RuntimeError(f"HTTP {code} for {url}") from None
			except Exception as exc:
				last_exc = exc
				import random
				sleep_s = min(30.0, 1.0 + attempt) + random.uniform(0.0, 2.0)
				time.sleep(sleep_s)

	raise RuntimeError(f"Failed to download after retries: {url} ({last_exc})") from None


def _target_corner_world_coords(target_wcs, shape: tuple[int, int]):
	from astropy.coordinates import SkyCoord

	ny, nx = shape
	pix = [(0, 0), (0, ny - 1), (nx - 1, 0), (nx - 1, ny - 1)]
	world = [target_wcs.pixel_to_world(x, y) for (x, y) in pix]
	return SkyCoord(world)


def _cutout_covers_target(
	*,
	legacy_wcs,
	legacy_shape: tuple[int, int],
	target_wcs,
	target_shape: tuple[int, int],
) -> bool:
	# legacy_shape is (h, w)
	coords = _target_corner_world_coords(target_wcs, target_shape)
	x, y = legacy_wcs.world_to_pixel(coords)
	w = legacy_shape[1]
	h = legacy_shape[0]
	# allow a tiny epsilon
	eps = 1e-3
	return (
		float(x.min()) >= -eps
		and float(y.min()) >= -eps
		and float(x.max()) <= (w - 1) + eps
		and float(y.max()) <= (h - 1) + eps
	)


def _legacy_wcs_for_cutout(
	*,
	ra_deg: float,
	dec_deg: float,
	pixscale_arcsec: float,
	width_px: int,
	height_px: int,
):
	"""Build a simple north-up TAN WCS matching the Legacy viewer cutout."""
	from astropy.wcs import WCS

	w = WCS(naxis=2)
	w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
	w.wcs.cunit = ["deg", "deg"]
	w.wcs.crval = [float(ra_deg), float(dec_deg)]
	# FITS CRPIX is 1-based; place reference at image center.
	w.wcs.crpix = [width_px / 2.0 + 0.5, height_px / 2.0 + 0.5]
	scale_deg = float(pixscale_arcsec) / 3600.0
	# East-left convention: RA decreases with +x.
	w.wcs.cdelt = [-scale_deg, scale_deg]
	return w


def _reproject_rgb_to_target(
	*,
	legacy_rgb_u8,
	legacy_wcs,
	target_wcs,
	shape_out: tuple[int, int],
):
	import numpy as np
	from reproject import reproject_adaptive

	if legacy_rgb_u8.ndim != 3 or legacy_rgb_u8.shape[2] != 3:
		raise ValueError(f"Expected RGB image (H,W,3), got {legacy_rgb_u8.shape}")

	out_channels: list[np.ndarray] = []
	for ch in range(3):
		arr = legacy_rgb_u8[:, :, ch].astype(np.float32)
		rep, _fp = reproject_adaptive(
			(arr, legacy_wcs),
			target_wcs,
			shape_out=shape_out,
			conserve_flux=True,
		)
		# Outside footprint will be NaN; use 0 so background is black.
		rep = np.nan_to_num(rep, nan=0.0, posinf=255.0, neginf=0.0)
		out_channels.append(rep)

	stack = np.stack(out_channels, axis=2)
	stack = np.clip(stack, 0.0, 255.0).astype(np.uint8)
	return stack


def _verify_expected_size(
	*,
	galaxy_id: str,
	shape: tuple[int, int],
	observed_png: pathlib.Path,
):
	# shape is (ny, nx)
	from PIL import Image

	ny, nx = shape
	if observed_png.exists():
		with Image.open(observed_png) as im:
			w, h = im.size
			if (w, h) != (nx, ny):
				raise ValueError(
					f"Size mismatch for {galaxy_id}: FITS is {nx}x{ny} but {observed_png.name} is {w}x{h}. "
					"Regenerate observed PNGs from the current FITS with: python v3tk_observed_VRI_image.py"
				)


def _process_one(job: Job, args: argparse.Namespace) -> tuple[str, str]:
	import numpy as np
	from PIL import Image
	from astropy.io import fits
	from astropy.wcs import WCS

	target_wcs, shape_out = _muse_target_wcs_and_shape(job.input_path)
	_verify_expected_size(galaxy_id=job.galaxy_id, shape=shape_out, observed_png=job.observed_png)

	# Determine center sky coordinate from the *image center pixel*.
	ny, nx = shape_out
	xc = (nx - 1) / 2.0
	yc = (ny - 1) / 2.0
	center = target_wcs.pixel_to_world(xc, yc)
	ra0 = float(center.ra.deg)
	dec0 = float(center.dec.deg)

	# Determine initial cutout size based on MUSE footprint.
	width_arcsec, height_arcsec = _compute_footprint_size_arcsec(target_wcs, shape_out)
	width_arcsec *= float(args.margin_factor)
	height_arcsec *= float(args.margin_factor)

	width_px0 = max(32, int(math.ceil(width_arcsec / float(args.pixscale))))
	height_px0 = max(32, int(math.ceil(height_arcsec / float(args.pixscale))))

	# Note: do not use CRVAL here; use true center computed above.

	if args.dry_run:
		return (
			job.galaxy_id,
			f"DRY-RUN: would download ~{width_px0}x{height_px0} @ {args.pixscale} arcsec/pix (auto-enlarge if needed)",
		)

	# Skip work if output exists and not overwriting.
	if job.legacy_reprojected_jpg.exists() and (not args.overwrite):
		return (job.galaxy_id, f"skip (exists): {job.legacy_reprojected_jpg.name}")

	# We download a FITS cutout to get an authoritative WCS, and a JPEG cutout for RGB values.
	# If the cutout does not fully cover the target footprint, enlarge and retry.
	last_w = width_px0
	last_h = height_px0
	legacy_wcs = None
	legacy_rgb_u8 = None
	last_debug = None
	for attempt in range(3):
		width_px = int(math.ceil(last_w))
		height_px = int(math.ceil(last_h))
		fits_tmp = job.legacy_cutout_jpg.with_suffix(".fits")
		jpg_url = _legacy_cutout_url(
			ra_deg=ra0,
			dec_deg=dec0,
			layer=args.layer,
			pixscale=float(args.pixscale),
			width_px=width_px,
			height_px=height_px,
		)
		fits_url = _legacy_fits_cutout_url(
			ra_deg=ra0,
			dec_deg=dec0,
			layer=args.layer,
			pixscale=float(args.pixscale),
			width_px=width_px,
			height_px=height_px,
			bands="r",
		)

		_download_file(url=fits_url, out_path=fits_tmp, timeout_s=float(args.timeout))
		with fits.open(fits_tmp, memmap=True) as hdul:
			# Legacy FITS cutout stores the image in primary.
			hdr = hdul[0].header
			legacy_wcs = WCS(hdr).celestial
			# Use NAXIS for shape.
			legacy_h = int(hdr.get("NAXIS2", height_px))
			legacy_w = int(hdr.get("NAXIS1", width_px))

		_download_file(url=jpg_url, out_path=job.legacy_cutout_jpg, timeout_s=float(args.timeout))
		with Image.open(job.legacy_cutout_jpg) as im:
			im = im.convert("RGB")
			legacy_rgb_u8 = np.asarray(im, dtype=np.uint8)

		# Sanity: JPEG and FITS should agree on dimensions; if not, trust JPEG.
		legacy_shape = (legacy_rgb_u8.shape[0], legacy_rgb_u8.shape[1])
		if (legacy_h, legacy_w) != legacy_shape:
			legacy_h, legacy_w = legacy_shape

		# Debug: check that the legacy WCS center pixel maps to the requested center.
		try:
			cx = (legacy_w - 1) / 2.0
			cy = (legacy_h - 1) / 2.0
			csky = legacy_wcs.pixel_to_world(cx, cy)
			dra = (float(csky.ra.deg) - ra0) * 3600.0 * math.cos(math.radians(dec0))
			ddec = (float(csky.dec.deg) - dec0) * 3600.0
			last_debug = (dra, ddec, legacy_w, legacy_h)
		except Exception:
			last_debug = None

		if _cutout_covers_target(
			legacy_wcs=legacy_wcs,
			legacy_shape=(legacy_h, legacy_w),
			target_wcs=target_wcs,
			target_shape=shape_out,
		):
			# Cleanup FITS temp immediately.
			try:
				fits_tmp.unlink(missing_ok=True)
			except Exception:
				pass
			break

		# Not covered: enlarge and retry.
		try:
			fits_tmp.unlink(missing_ok=True)
		except Exception:
			pass
		last_w = last_w * 1.25
		last_h = last_h * 1.25
		if attempt == 2:
			raise ValueError("Legacy cutout did not cover full MUSE footprint after retries; try increasing --margin-factor")

	if legacy_wcs is None or legacy_rgb_u8 is None:
		raise RuntimeError("Internal error: missing Legacy cutout/WCS")

	# Reproject to MUSE grid.
	stack_u8 = _reproject_rgb_to_target(
		legacy_rgb_u8=legacy_rgb_u8,
		legacy_wcs=legacy_wcs,
		target_wcs=target_wcs,
		shape_out=shape_out,
	)

	# Verify final image size matches target and observed PNG.
	ny, nx = shape_out
	if (stack_u8.shape[0], stack_u8.shape[1]) != (ny, nx):
		raise ValueError(
			f"Internal error for {job.galaxy_id}: reprojected image shape {stack_u8.shape} != expected {(ny, nx)}"
		)
	_verify_expected_size(galaxy_id=job.galaxy_id, shape=shape_out, observed_png=job.observed_png)

	out_im = Image.fromarray(stack_u8)
	out_im.save(job.legacy_reprojected_jpg, format="JPEG", quality=95)

	if not args.keep_cutout:
		try:
			job.legacy_cutout_jpg.unlink(missing_ok=True)
		except Exception:
			pass

	msg = f"wrote: {job.legacy_reprojected_jpg.name} ({nx}x{ny}) | muse-center RA={ra0:.6f} Dec={dec0:.6f}"
	if last_debug is not None:
		dra, ddec, lw, lh = last_debug
		msg += f" | legacy-center offset (arcsec): dRA={dra:+.3f} dDec={ddec:+.3f} | cutout={lw}x{lh}"
	return (job.galaxy_id, msg)


def _discover_jobs(input_dir: pathlib.Path, pattern: str) -> list[Job]:
	paths = sorted(input_dir.glob(pattern))
	j: list[Job] = []
	for p in paths:
		gid = _galaxy_id_from_filename(p)
		observed_png = input_dir / f"{gid}_observed_VRI.png"
		legacy_cutout_jpg = input_dir / f"{gid}_Legacy.jpg"
		legacy_reprojected_jpg = input_dir / f"{gid}_legacy_reprojected.jpg"
		j.append(
			Job(
				input_path=p,
				galaxy_id=gid,
				observed_png=observed_png,
				legacy_cutout_jpg=legacy_cutout_jpg,
				legacy_reprojected_jpg=legacy_reprojected_jpg,
			)
		)
	return j


def main(argv: list[str]) -> int:
	t0 = time.perf_counter()
	try:
		_require_deps()
		args = _parse_args(argv)
		input_dir = args.input_dir.resolve()
		jobs = _discover_jobs(input_dir=input_dir, pattern=args.pattern)
		if not jobs:
			raise FileNotFoundError(f"No files matched pattern '{args.pattern}' in {input_dir}")

		if not args.quiet:
			print(f"Found {len(jobs)} file(s)")
			print(f"Layer: {args.layer}")
			print(f"Legacy pixscale: {args.pixscale} arcsec/pix")
			print(f"Margin factor: {args.margin_factor}")

		workers = args.workers
		if workers is None:
			workers = min(2, os.cpu_count() or 1)
		workers = max(1, int(workers))

		max_dl = max(1, int(args.max_concurrent_downloads))
		globals()["_DOWNLOAD_SEM"] = threading.Semaphore(max_dl)

		n_ok = 0
		n_err = 0
		# Use threads: allows shared throttling and avoids pickling urllib exceptions.
		with cf.ThreadPoolExecutor(max_workers=workers) as ex:
			futs = [ex.submit(_process_one, job, args) for job in jobs]
			for fut in cf.as_completed(futs):
				try:
					gid, msg = fut.result()
					n_ok += 1
					if not args.quiet:
						print(f"[{gid}] {msg}")
				except Exception as exc:
					n_err += 1
					print(f"[ERROR] {exc}")

		if not args.quiet:
			dt = time.perf_counter() - t0
			print(f"Workers: {workers}")
			print(f"Max concurrent downloads: {max_dl}")
			print(f"Done: {n_ok}/{len(jobs)} processed in {dt:.1f}s")
		return 0 if n_err == 0 else 1
	except Exception as exc:
		print(f"[FATAL] {exc}")
		return 2


if __name__ == "__main__":
	raise SystemExit(main(sys.argv[1:]))

