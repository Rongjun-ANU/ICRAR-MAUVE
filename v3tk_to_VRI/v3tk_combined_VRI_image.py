#!/usr/bin/env python

from __future__ import annotations

import argparse
import concurrent.futures as cf
import os
import pathlib
import sys
import time
from dataclasses import dataclass


def _require_deps():
	try:
		import numpy as np  # noqa: F401
		from astropy.io import fits  # noqa: F401
		from PIL import Image  # noqa: F401
	except Exception as exc:  # pragma: no cover
		import sys as _sys
		raise RuntimeError(
			"Missing required dependency.\n\n"
			"Required: numpy astropy pillow\n\n"
			f"Interpreter used: {_sys.executable}\n"
			f"Python version: {_sys.version.splitlines()[0]}\n\n"
			"Try:\n"
			"  conda activate ICRAR\n"
			"  python v3tk_combined_VRI_image.py --dry-run\n\n"
			"Install deps if needed:\n"
			"  conda install -c conda-forge numpy astropy pillow\n"
			"  # or: python -m pip install numpy astropy pillow\n"
		) from exc


@dataclass(frozen=True)
class Job:
	input_fits: pathlib.Path
	galaxy_id: str
	observed_png: pathlib.Path
	legacy_reprojected_jpg: pathlib.Path
	output_png: pathlib.Path


def _parse_args(argv: list[str]) -> argparse.Namespace:
	p = argparse.ArgumentParser(
		prog="v3tk_combined_VRI_image.py",
		description=(
			"Combine MUSE observed R-band rendering with Legacy Survey background. "
			"Uses the FITS flux map to build a valid-pixel mask (finite & >0). "
			"Where pixels are invalid (NaN or <=0), it uses the Legacy reprojected image as background; "
			"otherwise it uses the colored observed_VRI.png."
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
		"--suffix",
		default="_combined_VRI.png",
		help="Output suffix (default: _combined_VRI.png)",
	)
	p.add_argument(
		"--workers",
		type=int,
		default=None,
		help="Parallel workers (default: os.cpu_count())",
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
		"--dry-run",
		action="store_true",
		help="List planned outputs but do not write files",
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
	# Prefer explicit naming.
	for key in ("R_FLUX", "FLUX", "RFLUX"):
		try:
			return hdul[key]
		except Exception:
			pass

	# Fallback: first 2D image HDU.
	for hdu in hdul:
		data = getattr(hdu, "data", None)
		if data is None or getattr(data, "ndim", None) != 2:
			continue
		return hdu

	raise ValueError("Could not find a 2D flux image HDU")


def _discover_jobs(input_dir: pathlib.Path, pattern: str, suffix: str) -> list[Job]:
	paths = sorted(input_dir.glob(pattern))
	j: list[Job] = []
	for p in paths:
		gid = _galaxy_id_from_filename(p)
		observed_png = input_dir / f"{gid}_observed_VRI.png"
		legacy_jpg = input_dir / f"{gid}_legacy_reprojected.jpg"
		out = input_dir / f"{gid}{suffix}"
		j.append(
			Job(
				input_fits=p,
				galaxy_id=gid,
				observed_png=observed_png,
				legacy_reprojected_jpg=legacy_jpg,
				output_png=out,
			)
		)
	return j


def _combine_one(job: Job, overwrite: bool) -> tuple[str, str]:
	import numpy as np
	from astropy.io import fits
	from PIL import Image

	if not job.observed_png.exists():
		raise FileNotFoundError(f"Missing observed PNG: {job.observed_png.name}")
	if not job.legacy_reprojected_jpg.exists():
		raise FileNotFoundError(f"Missing legacy reprojected JPG: {job.legacy_reprojected_jpg.name}")

	if job.output_png.exists() and (not overwrite):
		return (job.galaxy_id, f"skip (exists): {job.output_png.name}")

	with fits.open(job.input_fits, memmap=True) as hdul:
		h = _select_flux_hdu(hdul)
		data = h.data
		if data is None or getattr(data, "ndim", None) != 2:
			raise ValueError("Selected flux HDU does not contain a 2D image")
		arr = np.asarray(data, dtype=np.float32)

	finite = np.isfinite(arr)
	positive = finite & (arr > 0)
	# observed_VRI.png is written flipped vertically to match origin='lower'; match that here.
	mask = np.flipud(positive)

	with Image.open(job.observed_png) as im_obs:
		obs = im_obs.convert("RGB")
		obs_u8 = np.asarray(obs, dtype=np.uint8)
	with Image.open(job.legacy_reprojected_jpg) as im_leg:
		leg = im_leg.convert("RGB")
		leg_u8 = np.asarray(leg, dtype=np.uint8)

	ny, nx = arr.shape
	if obs_u8.shape[0] != ny or obs_u8.shape[1] != nx:
		raise ValueError(
			f"Size mismatch for {job.galaxy_id}: FITS is {nx}x{ny} but {job.observed_png.name} is {obs_u8.shape[1]}x{obs_u8.shape[0]}"
		)
	if leg_u8.shape[0] != ny or leg_u8.shape[1] != nx:
		raise ValueError(
			f"Size mismatch for {job.galaxy_id}: FITS is {nx}x{ny} but {job.legacy_reprojected_jpg.name} is {leg_u8.shape[1]}x{leg_u8.shape[0]}"
		)

	if mask.shape != (ny, nx):
		raise RuntimeError("Internal error: mask shape mismatch")

	out = leg_u8.copy()
	out[mask] = obs_u8[mask]
	Image.fromarray(out).save(job.output_png, format="PNG")
	return (job.galaxy_id, f"wrote: {job.output_png.name} ({nx}x{ny})")


def main(argv: list[str]) -> int:
	t0 = time.perf_counter()
	try:
		_require_deps()
		args = _parse_args(argv)
		input_dir = args.input_dir.resolve()
		jobs = _discover_jobs(input_dir=input_dir, pattern=args.pattern, suffix=args.suffix)
		if not jobs:
			raise FileNotFoundError(f"No files matched pattern '{args.pattern}' in {input_dir}")

		if args.dry_run:
			for j in jobs:
				print(
					f"{j.input_fits.name} + {j.observed_png.name} + {j.legacy_reprojected_jpg.name} -> {j.output_png.name}"
				)
			return 0

		workers = args.workers
		if workers is None:
			workers = os.cpu_count() or 1
		workers = max(1, int(workers))

		ok = 0
		fail = 0
		with cf.ProcessPoolExecutor(max_workers=workers) as ex:
			futs = [ex.submit(_combine_one, j, args.overwrite) for j in jobs]
			for fut in cf.as_completed(futs):
				try:
					gid, msg = fut.result()
					ok += 1
					if not args.quiet:
						print(f"[{gid}] {msg}")
				except Exception as exc:
					fail += 1
					print(f"[ERROR] {exc}", file=sys.stderr)

		dt = time.perf_counter() - t0
		if not args.quiet:
			print(f"Done. ok={ok} fail={fail} workers={workers} runtime={dt:.2f}s")
		return 0 if fail == 0 else 2
	except Exception as exc:
		dt = time.perf_counter() - t0
		print(f"ERROR: {exc}", file=sys.stderr)
		print(f"Runtime: {dt:.2f} s", file=sys.stderr)
		return 2


if __name__ == "__main__":
	raise SystemExit(main(sys.argv[1:]))

