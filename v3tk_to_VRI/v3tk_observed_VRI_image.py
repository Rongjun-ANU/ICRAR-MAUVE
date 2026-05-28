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
		from astropy.visualization import make_lupton_rgb  # noqa: F401
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
			"  python v3tk_observed_VRI_image.py --dry-run\n\n"
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


@dataclass(frozen=True)
class RenderOptions:
	percentile_low: float
	percentile_high: float
	stretch: float
	Q: float
	gamma: float
	post_boost: float


def _parse_args(argv: list[str]) -> argparse.Namespace:
	p = argparse.ArgumentParser(
		prog="v3tk_observed_VRI_image.py",
		description=(
			"Extract observed V/R/I-band flux images (nanomaggy) from each '*_v3tk_VRI.fits' file "
			"and render an RGB composite to 'XXX_observed_VRI.png' and 'XXX_observed_VRI.pdf'. Runs in parallel for efficiency. "
			"Rendering uses a Lupton RGB (asinh) stretch with percentile-based scaling per galaxy. "
			"Channel mapping: I→R, R→G, V→B."
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
		help="Glob pattern to match VRI FITS files",
	)
	p.add_argument(
		"--percentile-low",
		type=float,
		default=1.0,
		help="Lower percentile for background subtraction/scaling (default: 1)",
	)
	p.add_argument(
		"--percentile-high",
		type=float,
		default=99.9,
		help="Upper percentile for scaling (default: 99.9)",
	)
	p.add_argument(
		"--stretch",
		type=float,
		default=1.0,
		help="Lupton RGB stretch parameter (default: 1.0)",
	)
	p.add_argument(
		"--Q",
		type=float,
		default=8.0,
		help="Lupton RGB softening parameter Q (default: 8.0)",
	)
	p.add_argument(
		"--gamma",
		type=float,
		default=1.0,
		help=(
			"Gamma correction applied to scaled channels before Lupton RGB. "
			"<1 brightens mid-tones; >1 darkens. Default: 1.0"
		),
	)
	p.add_argument(
		"--post-boost",
		dest="post_boost",
		type=float,
		default=2.5,
		help="Global brightness multiplier applied to the final RGB image after rendering (default: 2.5)",
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
	marker = "_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk_VRI.fits"
	if name.endswith(marker):
		return name[: -len(marker)]
	if name.endswith("_v3tk_VRI.fits"):
		return name[: -len("_v3tk_VRI.fits")]
	# Fallback: take prefix before first underscore.
	return name.split("_", 1)[0]


def _select_band_flux_hdu(hdul, band: str):
	"""Return the HDU that contains the band flux map (nanomaggy) for a given band.

	Expected output layout from v3tk_to_VRI.py:
	- EXTNAME / HDU name: '{band}_FLUX' (e.g., 'V_FLUX', 'R_FLUX', 'I_FLUX')
	"""
	band = band.upper()
	# Prefer explicit naming.
	for key in (f"{band}_FLUX", f"{band}FLUX"):
		try:
			return hdul[key]
		except Exception:
			pass

	# Heuristic: any 2D image HDU with FILTER=band and nanomaggy units.
	for hdu in hdul:
		data = getattr(hdu, "data", None)
		if data is None or getattr(data, "ndim", None) != 2:
			continue
		bunit = (hdu.header.get("BUNIT") or "").lower()
		filt = (hdu.header.get("FILTER") or "").upper()
		if filt == band and ("nanomag" in bunit or "nmgy" in bunit):
			return hdu

	raise KeyError(f"Could not find {band}_FLUX HDU (or equivalent) in FITS")


def _log_scale_to_unit(arr, positive_mask, *, vmin=None, vmax=None):
	"""Log-scale an array to [0,1] using matplotlib LogNorm.

	Non-positive or non-finite pixels are returned as 0.
	"""
	from matplotlib.colors import LogNorm
	import numpy as np

	if vmin is None or vmax is None:
		pos_vals = arr[positive_mask]
		vmin = float(pos_vals.min())
		vmax = float(pos_vals.max())
		if not np.isfinite(vmin) or not np.isfinite(vmax):
			raise ValueError("Could not determine display range for image")
		if vmax <= vmin:
			vmax = vmin + 1e-6

	norm = LogNorm(vmin=vmin, vmax=vmax)
	# LogNorm returns masked array if input is masked; keep it simple and fill invalid with 0.
	out = np.zeros(arr.shape, dtype=np.float32)
	scaled = norm(arr)
	# norm(arr) yields floats; invalid/non-positive become <=0 or masked; use positive_mask.
	out[positive_mask] = np.asarray(scaled, dtype=np.float32)[positive_mask]
	out = np.clip(out, 0.0, 1.0)
	return out


def _prep_channel(arr) -> "object":
	"""Prepare an image channel for RGB composition.

	- Converts non-finite to 0
	- Clips negatives to 0 (flux should be >=0 for display)
	"""
	import numpy as np

	arr = np.asarray(arr, dtype=np.float32)
	arr = np.where(np.isfinite(arr), arr, 0.0).astype(np.float32, copy=False)
	arr[arr < 0] = 0.0
	return arr


def _scale_by_luminance_percentiles(r, g, b, p_low: float, p_high: float):
	"""Scale RGB channels using percentiles of a luminance image.

	This tends to preserve color differences better than normalizing each band independently.
	"""
	import numpy as np

	if not (0.0 <= p_low < p_high <= 100.0):
		raise ValueError("percentiles must satisfy 0 <= low < high <= 100")

	# Simple luminance proxy.
	lum = (r + g + b) / 3.0
	vals = lum[np.isfinite(lum) & (lum > 0)]
	if vals.size == 0:
		raise ValueError("RGB luminance has no positive finite pixels")
	lo = float(np.percentile(vals, p_low))
	hi = float(np.percentile(vals, p_high))
	if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
		hi = lo + 1e-6

	def _scale(x):
		y = x - lo
		y[y < 0] = 0
		y = y / (hi - lo)
		return np.clip(y, 0.0, 1.0).astype(np.float32, copy=False)

	return _scale(r), _scale(g), _scale(b)


def _gray_world_white_balance(r, g, b):
	"""Simple 'gray world' white balance.

	Scales channels so their robust medians match over moderately bright pixels.
	"""
	import numpy as np

	lum = (r + g + b) / 3.0
	vals = lum[np.isfinite(lum) & (lum > 0)]
	if vals.size == 0:
		return r, g, b

	lo = float(np.percentile(vals, 30.0))
	hi = float(np.percentile(vals, 99.5))
	mask = np.isfinite(lum) & (lum > lo) & (lum < hi)
	if int(mask.sum()) < 50:
		return r, g, b

	mr = float(np.median(r[mask]))
	mg = float(np.median(g[mask]))
	mb = float(np.median(b[mask]))

	# Avoid division by ~0
	eps = 1e-12
	mr = max(mr, eps)
	mg = max(mg, eps)
	mb = max(mb, eps)

	target = float(np.median([mr, mg, mb]))
	gr = target / mr
	gg = target / mg
	gb = target / mb

	# Keep memory / dtype stable.
	r = (r * gr).astype(np.float32, copy=False)
	g = (g * gg).astype(np.float32, copy=False)
	b = (b * gb).astype(np.float32, copy=False)
	return r, g, b


def _extract_one(job: Job, overwrite: bool, opts: RenderOptions) -> tuple[str, str]:
	# Import inside worker for ProcessPool compatibility.
	from astropy.io import fits
	from astropy.visualization import make_lupton_rgb
	import numpy as np
	import warnings

	import matplotlib
	matplotlib.use("Agg")
	from matplotlib.figure import Figure
	from PIL import Image

	with fits.open(job.input_path, memmap=True) as hdul:
		h_v = _select_band_flux_hdu(hdul, "V")
		h_r = _select_band_flux_hdu(hdul, "R")
		h_i = _select_band_flux_hdu(hdul, "I")
		v = h_v.data
		r = h_r.data
		i = h_i.data
		for band, data in (("V", v), ("R", r), ("I", i)):
			if data is None or getattr(data, "ndim", None) != 2:
				raise ValueError(f"Selected {band}_FLUX HDU does not contain a 2D image")

	# Prepare channels. Mapping: I->R, R->G, V->B.
	r_raw = _prep_channel(i)
	g_raw = _prep_channel(r)
	b_raw = _prep_channel(v)

	# Pre-balance channels (helps reduce a persistent warm/orange cast).
	r_raw, g_raw, b_raw = _gray_world_white_balance(r_raw, g_raw, b_raw)

	missing: list[str] = []
	# Track missing/empty inputs (after prep, empty means all zeros).
	for band, arr in (("I", r_raw), ("R", g_raw), ("V", b_raw)):
		if not np.any(arr > 0):
			missing.append(f"{band}(empty)")

	if not (np.any(r_raw > 0) or np.any(g_raw > 0) or np.any(b_raw > 0)):
		raise ValueError("All V/R/I bands are empty after masking; cannot render RGB")

	# Percentile scaling on luminance, optional gamma, then Lupton/asinh stretch.
	r_s, g_s, b_s = _scale_by_luminance_percentiles(r_raw, g_raw, b_raw, opts.percentile_low, opts.percentile_high)
	if not (opts.gamma > 0):
		raise ValueError("gamma must be > 0")
	if opts.gamma != 1.0:
		r_s = np.power(r_s, opts.gamma, dtype=np.float32)
		g_s = np.power(g_s, opts.gamma, dtype=np.float32)
		b_s = np.power(b_s, opts.gamma, dtype=np.float32)
	# Important: render in float space and only quantize at the end. This avoids
	# banding/ring artifacts when a large post-boost is applied.
	# make_lupton_rgb can emit benign RuntimeWarnings (e.g. 0/0 in chroma terms);
	# suppress them so --quiet stays quiet.
	with np.errstate(divide="ignore", invalid="ignore"), warnings.catch_warnings():
		warnings.filterwarnings(
			"ignore",
			category=RuntimeWarning,
			module=r"astropy\\.visualization\\.lupton_rgb",
		)
		rgb_float = make_lupton_rgb(
			r_s,
			g_s,
			b_s,
			Q=opts.Q,
			stretch=opts.stretch,
			output_dtype=float,
		)
	rgb_float = np.asarray(rgb_float, dtype=np.float32)
	# Lupton can emit NaNs in low-intensity regions (e.g. 0/0 in chroma terms).
	rgb_float = np.nan_to_num(rgb_float, nan=0.0, posinf=1.0, neginf=0.0)

	# Systematic brightness scale-up applied AFTER rendering (still in float).
	if not np.isfinite(opts.post_boost) or opts.post_boost <= 0:
		raise ValueError("post_boost must be a finite number > 0")
	if opts.post_boost != 1.0:
		rgb_float = np.clip(rgb_float * float(opts.post_boost), 0.0, 1.0)
	else:
		rgb_float = np.clip(rgb_float, 0.0, 1.0)

	# Quantize once at the very end.
	rgb_u8 = np.clip(rgb_float * 255.0 + 0.5, 0.0, 255.0).astype(np.uint8)

	# "Native resolution": PNG dimensions match data shape exactly.
	ny, nx = g_raw.shape
	# Match the PDF/imshow convention (origin='lower'): flip vertically for image file output.
	rgb_u8_png = np.flipud(rgb_u8)
	if overwrite or (not job.output_png.exists()):
		Image.fromarray(rgb_u8_png, mode="RGB").save(job.output_png, format="PNG")

	# Write PDF (not pixel-perfect, but preserves appearance).
	dpi = 100
	fig = Figure(figsize=(nx / dpi, ny / dpi), dpi=dpi)
	fig.patch.set_facecolor("black")
	fig.patch.set_alpha(1.0)
	ax = fig.add_axes([0, 0, 1, 1])
	ax.set_facecolor("black")
	ax.set_axis_off()
	ax.imshow(
		rgb_float,
		origin="lower",
		interpolation="nearest",
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

	out = f"{job.output_png.name}, {job.output_pdf.name}"
	if missing:
		out += f" (missing/empty: {', '.join(missing)})"
	return (job.galaxy_id, out)


def _discover_jobs(input_dir: pathlib.Path, pattern: str) -> list[Job]:
	paths = expand_fits_glob(input_dir, pattern)
	jobs: list[Job] = []
	for p in paths:
		gid = _galaxy_id_from_filename(p)
		out_png = input_dir / f"{gid}_observed_VRI.png"
		out_pdf = input_dir / f"{gid}_observed_VRI.pdf"
		jobs.append(Job(input_path=p, galaxy_id=gid, output_png=out_png, output_pdf=out_pdf))
	return jobs


def main(argv: list[str]) -> int:
	t0 = time.perf_counter()
	try:
		_require_deps()
		args = _parse_args(argv)
		opts = RenderOptions(
			percentile_low=float(args.percentile_low),
			percentile_high=float(args.percentile_high),
			stretch=float(args.stretch),
			Q=float(args.Q),
			gamma=float(args.gamma),
			post_boost=float(args.post_boost),
		)
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
			futs = [ex.submit(_extract_one, j, args.overwrite, opts) for j in jobs]
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
