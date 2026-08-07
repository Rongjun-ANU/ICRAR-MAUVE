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
	legacy_reprojected_jpg: pathlib.Path
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
	target_max_luminance: float
	center_radius_fraction: float
	legacy_transition_width: int


def _parse_args(argv: list[str]) -> argparse.Namespace:
	p = argparse.ArgumentParser(
		prog="v3tk_observed_VRI_image.py",
		description=(
			"Extract observed V/R/I-band flux images (nanomaggy) from each '*_DATACUBE*_VRI.fits' file "
			"and render an RGB composite to 'XXX_observed_VRI.png' and 'XXX_observed_VRI.pdf'. Runs in parallel for efficiency. "
			"Rendering uses a Lupton RGB (asinh) stretch with percentile-based scaling per galaxy, "
			"then matches its inner edge to the reprojected Legacy image and normalizes its display maximum. "
			"Channel mapping: I→R, R→G, V→B."
		),
	)
	p.add_argument(
		"--input-dir",
		type=pathlib.Path,
		default=pathlib.Path("."),
		help="Directory containing *_DATACUBE*_VRI.fits files (default: current directory)",
	)
	p.add_argument(
		"--pattern",
		default="*_DATACUBE*_VRI.fits",
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
		default=0.65,
		help=(
			"Color-preserving gamma applied to the final Lupton RGB intensity. "
			"<1 brightens faint and mid-tone detail while retaining RGB ratios; "
			">1 darkens. Default: 0.65"
		),
	)
	p.add_argument(
		"--post-boost",
		dest="post_boost",
		type=float,
		default=2.75,
		help=(
			"Strength of the color-preserving gamma lift in shadows and mid-tones. "
			"The lift fades to zero toward white to protect highlight detail. Default: 2.75"
		),
	)
	p.add_argument(
		"--target-max-luminance",
		dest="target_max_luminance",
		type=float,
		default=255.0,
		help=(
			"Common maximum Rec.709 display luminance for the galaxy-centre "
			"reference aperture, on the 0-255 pixel scale. Default: 255"
		),
	)
	p.add_argument(
		"--center-radius-fraction",
		dest="center_radius_fraction",
		type=float,
		default=0.22,
		help=(
			"Radius of the galaxy-centre brightness-reference aperture as a fraction "
			"of the shorter valid-footprint dimension. Foreground stars outside this "
			"aperture do not set the global VRI scale. Default: 0.22"
		),
	)
	p.add_argument(
		"--legacy-transition-width",
		"--legacy-feather-width",
		dest="legacy_transition_width",
		type=int,
		default=20,
		help=(
			"Width in pixels over which only the smooth Legacy color/brightness "
			"baseline transitions to pure VRI while VRI detail stays sharp. "
			"The old --legacy-feather-width spelling is retained as an alias. Default: 20"
		),
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


def _apply_color_preserving_brightness(rgb, *, gamma: float, boost: float):
	"""Lift shadows through one shared scale factor per RGB pixel.

	The gamma lift is weighted by ``1 - max(R, G, B)``, so its effect fades to
	zero in the highlights. All three channels receive the same per-pixel
	multiplier, retaining their ratios without flattening bright regions to white.
	"""
	import numpy as np

	if not np.isfinite(gamma) or gamma <= 0:
		raise ValueError("gamma must be a finite number > 0")
	if not np.isfinite(boost) or boost <= 0:
		raise ValueError("post_boost must be a finite number > 0")

	rgb = np.asarray(rgb, dtype=np.float32)
	if rgb.ndim != 3 or rgb.shape[-1] != 3:
		raise ValueError("rgb must have shape (height, width, 3)")
	rgb = np.clip(rgb, 0.0, None)
	maximum = np.max(rgb, axis=-1)
	gamma_lifted = np.power(maximum, gamma)
	target_maximum = np.clip(
		maximum + boost * (gamma_lifted - maximum) * (1.0 - maximum),
		0.0,
		1.0,
	)
	scale = np.divide(
		target_maximum,
		maximum,
		out=np.zeros_like(maximum, dtype=np.float32),
		where=maximum > 0,
	)
	return (rgb * scale[..., None]).astype(np.float32, copy=False)


def _inner_transition_alpha(valid_mask, width: int):
	"""Return VRI opacity that rises from 0 at the inner edge to 1 inward."""
	import numpy as np

	mask = np.asarray(valid_mask, dtype=bool)
	if mask.ndim != 2:
		raise ValueError("valid_mask must be a 2D array")
	if isinstance(width, bool) or not isinstance(width, (int, np.integer)) or width < 0:
		raise ValueError("transition width must be an integer >= 0")

	alpha = mask.astype(np.float32)
	if width == 0 or not np.any(mask):
		return alpha

	remaining = mask.copy()
	for layer in range(width):
		padded = np.pad(remaining, 1, mode="constant", constant_values=False)
		eroded = np.logical_and.reduce(
			[
				padded[dy : dy + mask.shape[0], dx : dx + mask.shape[1]]
				for dy in range(3)
				for dx in range(3)
			]
		)
		boundary = remaining & ~eroded
		alpha[boundary] = layer / float(width)
		remaining = eroded
		if not np.any(remaining):
			break

	return alpha


def _box_blur_2d(arr, radius: int, *, pad_mode: str):
	"""Fast square low-pass filter with float64 accumulation for stable edges."""
	import numpy as np

	values = np.asarray(arr, dtype=np.float32)
	if values.ndim != 2:
		raise ValueError("arr must be a 2D array")
	if isinstance(radius, bool) or not isinstance(radius, (int, np.integer)) or radius < 0:
		raise ValueError("blur radius must be an integer >= 0")
	if pad_mode not in {"constant", "edge"}:
		raise ValueError("pad_mode must be 'constant' or 'edge'")
	if radius == 0:
		return values.copy()

	kernel_size = 2 * int(radius) + 1
	padded = np.pad(values, int(radius), mode=pad_mode)
	# Accumulate in float64: float32 integral-image subtraction can create a
	# colored one-pixel contour in large images through cancellation error.
	integral = np.pad(padded, ((1, 0), (1, 0)), mode="constant")
	integral = np.cumsum(np.cumsum(integral, axis=0, dtype=np.float64), axis=1, dtype=np.float64)
	sums = (
		integral[kernel_size:, kernel_size:]
		- integral[:-kernel_size, kernel_size:]
		- integral[kernel_size:, :-kernel_size]
		+ integral[:-kernel_size, :-kernel_size]
	)
	return (sums / float(kernel_size * kernel_size)).astype(np.float32)


def _central_reference_mask(valid_mask, *, radius_fraction: float):
	"""Return a circular brightness-reference mask around the VRI footprint centre."""
	import numpy as np

	mask = np.asarray(valid_mask, dtype=bool)
	if mask.ndim != 2:
		raise ValueError("valid_mask must be a 2D array")
	if not np.any(mask):
		raise ValueError("valid_mask contains no VRI pixels")
	if (
		isinstance(radius_fraction, bool)
		or not np.isfinite(radius_fraction)
		or not (0.0 < radius_fraction <= 0.5)
	):
		raise ValueError("center_radius_fraction must be finite and in (0, 0.5]")

	y_valid, x_valid = np.nonzero(mask)
	x_min, x_max = int(x_valid.min()), int(x_valid.max())
	y_min, y_max = int(y_valid.min()), int(y_valid.max())
	width = x_max - x_min + 1
	height = y_max - y_min + 1
	x_center = 0.5 * (x_min + x_max)
	y_center = 0.5 * (y_min + y_max)
	radius = float(radius_fraction) * min(width, height)
	y_grid, x_grid = np.ogrid[: mask.shape[0], : mask.shape[1]]
	reference = mask & (
		(x_grid - x_center) ** 2 + (y_grid - y_center) ** 2 <= radius**2
	)
	if not np.any(reference):
		raise ValueError("centre reference aperture contains no valid VRI pixels")
	return reference


def _match_vri_edge_to_legacy(
	vri_rgb,
	legacy_rgb,
	valid_mask,
	*,
	transition_width: int,
	target_max_luminance: float,
	normalization_mask=None,
):
	"""Match the VRI edge illumination to Legacy without blurring VRI detail.

	Only Legacy pixels outside the valid VRI footprint contribute to the low-pass
	edge reference. That reference changes the low-frequency VRI baseline near the
	edge while the complete high-frequency VRI residual remains in the output.
	No raw Legacy pixel is copied or cross-faded inside the VRI footprint.
	"""
	import numpy as np

	vri = np.asarray(vri_rgb, dtype=np.float32)
	legacy = np.asarray(legacy_rgb, dtype=np.float32)
	mask = np.asarray(valid_mask, dtype=bool)
	if vri.ndim != 3 or vri.shape[-1] != 3:
		raise ValueError("vri_rgb must have shape (height, width, 3)")
	if legacy.shape != vri.shape:
		raise ValueError("legacy_rgb must have the same shape as vri_rgb")
	if mask.shape != vri.shape[:2]:
		raise ValueError("valid_mask must match the RGB image dimensions")
	if not np.isfinite(target_max_luminance) or not (0.0 < target_max_luminance <= 255.0):
		raise ValueError("target_max_luminance must be finite and in (0, 255]")
	if not np.any(mask):
		raise ValueError("valid_mask contains no VRI pixels")
	if normalization_mask is None:
		reference_mask = mask
	else:
		reference_mask = np.asarray(normalization_mask, dtype=bool)
		if reference_mask.shape != mask.shape:
			raise ValueError("normalization_mask must match valid_mask")
		reference_mask = mask & reference_mask
		if not np.any(reference_mask):
			raise ValueError("normalization_mask contains no valid VRI pixels")
	if (
		isinstance(transition_width, bool)
		or not isinstance(transition_width, (int, np.integer))
		or transition_width < 0
	):
		raise ValueError("legacy_transition_width must be an integer >= 0")

	vri = np.clip(np.nan_to_num(vri, nan=0.0, posinf=1.0, neginf=0.0), 0.0, 1.0)
	legacy = np.clip(np.nan_to_num(legacy, nan=0.0, posinf=1.0, neginf=0.0), 0.0, 1.0)
	target = float(target_max_luminance) / 255.0
	weights = np.array([0.2126, 0.7152, 0.0722], dtype=np.float32)

	# Scale the original VRI through one global multiplier set by the brightest
	# pixel in the galaxy-centre reference aperture. This preserves RGB ratios and
	# stops foreground stars outside the aperture from setting the galaxy scale.
	# Per-pixel shared gamut protection below prevents any highly colored pixel
	# from globally limiting the brightness.
	vri_luma = np.sum(vri * weights, axis=-1)
	current_max_luma = float(np.max(vri_luma[reference_mask]))
	if current_max_luma <= 0.0:
		raise ValueError("VRI image has no positive luminance inside the reference aperture")
	vri = vri * max(0.0, target / current_max_luma)

	matched = vri.copy()
	if transition_width > 0:
		transition_alpha = _inner_transition_alpha(mask, int(transition_width))
		edge_weight = 1.0 - transition_alpha
		mask_blur = _box_blur_2d(mask.astype(np.float32), int(transition_width), pad_mode="constant")
		outside_mask = (~mask).astype(np.float32)
		legacy_support = _box_blur_2d(
			outside_mask,
			int(transition_width),
			pad_mode="edge",
		)
		vri_baseline_channels = []
		legacy_baseline_channels = []
		for channel in range(3):
			vri_sum = _box_blur_2d(
				vri[..., channel] * mask,
				int(transition_width),
				pad_mode="constant",
			)
			vri_baseline_channels.append(
				np.divide(
					vri_sum,
					mask_blur,
					out=np.zeros_like(vri_sum),
					where=mask_blur > 1e-6,
				)
			)
			legacy_sum = _box_blur_2d(
				legacy[..., channel] * outside_mask,
				int(transition_width),
				pad_mode="edge",
			)
			legacy_baseline_channels.append(
				np.divide(
					legacy_sum,
					legacy_support,
					out=vri_baseline_channels[-1].copy(),
					where=legacy_support > 1e-6,
				)
			)
		vri_baseline = np.stack(vri_baseline_channels, axis=-1)
		legacy_baseline = np.stack(legacy_baseline_channels, axis=-1)
		matched = vri + edge_weight[..., None] * (legacy_baseline - vri_baseline)

	# Keep all color operations coupled while protecting the target luminance and
	# RGB gamut. High-frequency VRI detail remains present in ``matched``.
	matched = np.clip(matched, 0.0, None)
	maximum_channel = np.max(matched, axis=-1)
	matched = matched / np.maximum(maximum_channel, 1.0)[..., None]
	matched_luma = np.sum(matched * weights, axis=-1)
	luma_cap = np.minimum(
		1.0,
		np.divide(target, matched_luma, out=np.ones_like(matched_luma), where=matched_luma > 0.0),
	)
	matched = matched * luma_cap[..., None]

	out = np.zeros_like(vri, dtype=np.float32)
	out[mask] = matched[mask]
	return np.clip(out, 0.0, 1.0).astype(np.float32, copy=False)


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

	if not job.legacy_reprojected_jpg.exists():
		raise FileNotFoundError(
			f"Missing Legacy image required by observed renderer: {job.legacy_reprojected_jpg.name}. "
			"Run v3tk_get_legacy.py first."
		)

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
	valid_mask = g_raw > 0

	# Pre-balance channels (helps reduce a persistent warm/orange cast).
	r_raw, g_raw, b_raw = _gray_world_white_balance(r_raw, g_raw, b_raw)

	missing: list[str] = []
	# Track missing/empty inputs (after prep, empty means all zeros).
	for band, arr in (("I", r_raw), ("R", g_raw), ("V", b_raw)):
		if not np.any(arr > 0):
			missing.append(f"{band}(empty)")

	if not (np.any(r_raw > 0) or np.any(g_raw > 0) or np.any(b_raw > 0)):
		raise ValueError("All V/R/I bands are empty after masking; cannot render RGB")

	# Percentile scaling on luminance, then Lupton/asinh stretch.
	r_s, g_s, b_s = _scale_by_luminance_percentiles(r_raw, g_raw, b_raw, opts.percentile_low, opts.percentile_high)
	# Important: render in float space and only quantize at the end. This avoids
	# banding/ring artifacts when the post-render tone curve is applied.
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

	# Lift faint structure through one scale factor per pixel, with progressively
	# less gain toward white. This protects highlights and keeps RGB ratios coupled.
	rgb_float = _apply_color_preserving_brightness(
		rgb_float,
		gamma=opts.gamma,
		boost=opts.post_boost,
	)

	# "Native resolution": PNG dimensions match data shape exactly.
	ny, nx = g_raw.shape
	# Legacy JPG and PNG products use display orientation. Flip the VRI image and
	# mask first, then match its edge illumination to the Legacy background.
	with Image.open(job.legacy_reprojected_jpg) as im_legacy:
		legacy_rgb = np.asarray(im_legacy.convert("RGB"), dtype=np.float32) / 255.0
	if legacy_rgb.shape != (ny, nx, 3):
		raise ValueError(
			f"Size mismatch for {job.galaxy_id}: VRI FITS is {nx}x{ny} but "
			f"{job.legacy_reprojected_jpg.name} is {legacy_rgb.shape[1]}x{legacy_rgb.shape[0]}"
		)
	display_mask = np.flipud(valid_mask)
	normalization_mask = _central_reference_mask(
		display_mask,
		radius_fraction=opts.center_radius_fraction,
	)
	rgb_png_float = _match_vri_edge_to_legacy(
		np.flipud(rgb_float),
		legacy_rgb,
		display_mask,
		transition_width=opts.legacy_transition_width,
		target_max_luminance=opts.target_max_luminance,
		normalization_mask=normalization_mask,
	)

	# Quantize once at the very end.
	rgb_u8_png = np.clip(rgb_png_float * 255.0 + 0.5, 0.0, 255.0).astype(np.uint8)
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
		np.flipud(rgb_png_float),
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
		legacy_jpg = input_dir / f"{gid}_legacy_reprojected.jpg"
		out_png = input_dir / f"{gid}_observed_VRI.png"
		out_pdf = input_dir / f"{gid}_observed_VRI.pdf"
		jobs.append(
			Job(
				input_path=p,
				galaxy_id=gid,
				legacy_reprojected_jpg=legacy_jpg,
				output_png=out_png,
				output_pdf=out_pdf,
			)
		)
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
			target_max_luminance=float(args.target_max_luminance),
			center_radius_fraction=float(args.center_radius_fraction),
			legacy_transition_width=int(args.legacy_transition_width),
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
				print(
					f"{j.input_path.name} + {j.legacy_reprojected_jpg.name} "
					f"-> {j.output_png.name} + {j.output_pdf.name}"
				)
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
