#!/usr/bin/env python

from __future__ import annotations

import argparse
import concurrent.futures as cf
import os
import pathlib
import sys
import time
from dataclasses import dataclass

from fits_path_utils import expand_fits_glob, strip_optional_gzip_suffix


def _require_deps():
	try:
		import numpy as np  # noqa: F401
		from astropy.io import fits  # noqa: F401
		import matplotlib  # noqa: F401
		from PIL import Image  # noqa: F401
	except Exception as exc:  # pragma: no cover
		import sys as _sys
		raise RuntimeError(
			"Missing required dependency (numpy/astropy/matplotlib/pillow).\n\n"
			f"Interpreter used: {_sys.executable}\n"
			f"Python version: {_sys.version.splitlines()[0]}\n\n"
			"If you already installed astropy in your conda env, you are likely running the script with a different Python.\n"
			"Try:\n"
			"  conda activate ICRAR\n"
			"  python v3tk_observed_R_image.py --dry-run\n\n"
			"Install deps if needed:\n"
			"  conda install -c conda-forge numpy astropy matplotlib pillow\n"
			"  # or: python -m pip install numpy astropy matplotlib pillow"
		) from exc


@dataclass(frozen=True)
class Job:
	input_path: pathlib.Path
	galaxy_id: str
	output_png: pathlib.Path
	output_pdf: pathlib.Path


def _parse_args(argv: list[str]) -> argparse.Namespace:
	p = argparse.ArgumentParser(
		prog="v3tk_observed_R_image.py",
		description=(
			"Extract the observed R-band flux image (nanomaggy) from each '*_v3tk_R.fits' file "
			"and render it to 'XXX_observed_R.png' and 'XXX_observed_R.pdf'. Runs in parallel for efficiency. "
			"Rendering uses a log-scaled normalization per galaxy (LogNorm)."
		),
	)
	p.add_argument(
		"--input-dir",
		type=pathlib.Path,
		default=pathlib.Path("."),
		help="Directory containing *_v3tk_R.fits files (default: current directory)",
	)
	p.add_argument(
		"--pattern",
		default="*_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk_R.fits",
		help="Glob pattern to match R-band FITS files",
	)
	p.add_argument(
		"--workers",
		type=int,
		default=None,
		help="Number of parallel workers (default: os.cpu_count())",
	)
	p.add_argument(
		"--overwrite",
		dest="overwrite",
		action="store_true",
		help="Overwrite existing outputs (PNG/PDF) (default)",
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
	name = strip_optional_gzip_suffix(path)
	marker = "_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk_R.fits"
	if name.endswith(marker):
		return name[: -len(marker)]
	if name.endswith("_v3tk_R.fits"):
		return name[: -len("_v3tk_R.fits")]
	# Fallback: take prefix before first underscore.
	return name.split("_", 1)[0]


def _select_flux_hdu(hdul):
	"""Return the HDU that contains the flux map in nanomaggy.

	Handles common layouts:
	- [0]=flux, [1]=mag
	- [0]=primary(empty), [1]=flux (R_FLUX), [2]=mag (R_MAG)
	"""
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
		if "NANOMAG" in bunit.upper() or "NMGY" in bunit.upper() or extname == "R_FLUX" or name == "R_FLUX":
			return hdu

	# Fallback: first 2D image extension (skip empty primary).
	for idx, hdu in enumerate(hdul):
		data = getattr(hdu, "data", None)
		if data is None or getattr(data, "ndim", None) != 2:
			continue
		# Many FITS store images in extension 1+, but allow idx=0 if it really is an image.
		if idx == 0 and len(hdul) > 1:
			# If primary is image and there are other HDUs, still accept.
			return hdu
		return hdu

	raise ValueError("Could not find a 2D flux image HDU")


def _extract_one(job: Job, overwrite: bool) -> tuple[str, str]:
	# Import inside worker for ProcessPool compatibility.
	from astropy.io import fits
	import numpy as np

	import matplotlib
	matplotlib.use("Agg")
	from matplotlib import cm
	from matplotlib.colors import LogNorm
	from matplotlib.figure import Figure
	from PIL import Image

	with fits.open(job.input_path, memmap=True) as hdul:
		h_flux = _select_flux_hdu(hdul)
		data = h_flux.data
		if data is None or getattr(data, "ndim", None) != 2:
			raise ValueError("Selected flux HDU does not contain a 2D image")

	arr = np.asarray(data, dtype=np.float32)
	finite = np.isfinite(arr)
	if not np.any(finite):
		raise ValueError("Flux image contains no finite pixels")
	# Log scale is undefined for non-positive values: mask them to render as black.
	positive = finite & (arr > 0)
	if not np.any(positive):
		raise ValueError("Flux image contains no positive finite pixels (required for log colorbar)")
	img = np.ma.array(arr, mask=~positive)

	# Display scaling: per-galaxy full range of positive finite pixels (required for LogNorm).
	pos_vals = arr[positive]
	vmin = float(pos_vals.min())
	vmax = float(pos_vals.max())
	if not np.isfinite(vmin) or not np.isfinite(vmax):
		raise ValueError("Could not determine display range for image")
	if vmax <= vmin:
		# Constant-valued image: choose a tiny range so imshow doesn't warn.
		vmax = vmin + 1e-6
	norm = LogNorm(vmin=vmin, vmax=vmax)

	# "Native resolution": PNG dimensions match data shape exactly.
	ny, nx = arr.shape
	try:
		cmap = matplotlib.colormaps["plasma"].copy()
	except Exception:  # pragma: no cover
		cmap = cm.get_cmap("plasma")
		try:
			cmap = cmap.copy()
		except Exception:
			pass
	try:
		cmap.set_bad("black")
	except Exception:
		pass

	# Colormap to RGB bytes with exact pixel dimensions.
	rgba = cmap(norm(arr), bytes=True)  # (ny, nx, 4) uint8; NaNs map to 0
	rgb = rgba[:, :, :3].copy()
	rgb[~positive] = 0
	# Match the PDF/imshow convention (origin='lower'): flip vertically for image file output.
	rgb = np.flipud(rgb)
	if overwrite or (not job.output_png.exists()):
		Image.fromarray(rgb, mode="RGB").save(job.output_png, format="PNG")

	# Write PDF (not pixel-perfect, but preserves appearance).
	dpi = 100
	fig = Figure(figsize=(nx / dpi, ny / dpi), dpi=dpi)
	fig.patch.set_facecolor("black")
	fig.patch.set_alpha(1.0)
	ax = fig.add_axes([0, 0, 1, 1])
	ax.set_facecolor("black")
	ax.set_axis_off()
	ax.imshow(
		img,
		origin="lower",
		cmap=cmap,
		interpolation="nearest",
		norm=norm,
	)
	# Write PDF.
	if overwrite or (not job.output_pdf.exists()):
		fig.savefig(
			str(job.output_pdf),
			format="pdf",
			dpi=dpi,
			facecolor=fig.get_facecolor(),
			edgecolor="none",
		)

	return (job.galaxy_id, f"{job.output_png.name}, {job.output_pdf.name}")


def _discover_jobs(input_dir: pathlib.Path, pattern: str) -> list[Job]:
	paths = expand_fits_glob(input_dir, pattern)
	jobs: list[Job] = []
	for p in paths:
		gid = _galaxy_id_from_filename(p)
		out_png = input_dir / f"{gid}_observed_R.png"
		out_pdf = input_dir / f"{gid}_observed_R.pdf"
		jobs.append(Job(input_path=p, galaxy_id=gid, output_png=out_png, output_pdf=out_pdf))
	return jobs


def main(argv: list[str]) -> int:
	t0 = time.perf_counter()
	try:
		_require_deps()
		args = _parse_args(argv)
		input_dir = args.input_dir.resolve()
		jobs = _discover_jobs(input_dir=input_dir, pattern=args.pattern)

		if not jobs:
			raise FileNotFoundError(f"No files matched pattern '{args.pattern}' in {input_dir}")

		workers = args.workers
		if workers is None:
			workers = os.cpu_count() or 1
		workers = max(1, int(workers))

		if not args.quiet:
			print(f"Found {len(jobs)} file(s)")
			print(f"Workers: {workers}")
			print(f"Pattern: {args.pattern}")

		if args.dry_run:
			for j in jobs:
				print(f"{j.input_path.name} -> {j.output_png.name} + {j.output_pdf.name}")
			return 0

		# Use processes to avoid GIL overhead in FITS I/O / decompression.
		ok = 0
		fail = 0
		with cf.ProcessPoolExecutor(max_workers=workers) as ex:
			futs = [ex.submit(_extract_one, j, args.overwrite) for j in jobs]
			for fut in cf.as_completed(futs):
				try:
					gid, out = fut.result()
					ok += 1
					if not args.quiet:
						print(f"OK  {gid} -> {out}")
				except Exception as exc:
					fail += 1
					print(f"ERROR: {exc}", file=sys.stderr)

		dt = time.perf_counter() - t0
		if not args.quiet:
			print(f"Done. ok={ok} fail={fail} runtime={dt:.2f}s")
		return 0 if fail == 0 else 2
	except Exception as exc:
		dt = time.perf_counter() - t0
		print(f"ERROR: {exc}", file=sys.stderr)
		print(f"Runtime: {dt:.2f} s", file=sys.stderr)
		return 2


if __name__ == "__main__":
	raise SystemExit(main(sys.argv[1:]))
