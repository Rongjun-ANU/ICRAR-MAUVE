#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import pathlib
import sys
import time
from typing import Tuple


__version__ = "2026-01-15.1"


def _require_deps():
	try:
		import numpy as np  # noqa: F401
		from astropy.io import fits  # noqa: F401
		from astropy.wcs import WCS  # noqa: F401
		import astropy.units as u  # noqa: F401
	except Exception as exc:  # pragma: no cover
		raise RuntimeError(
			"Missing required dependency. Please install: numpy astropy\n"
			"Example (conda): conda install -c conda-forge numpy astropy\n"
			"Example (pip): pip install numpy astropy"
		) from exc

	try:
		from speclite import filters  # noqa: F401
	except Exception as exc:  # pragma: no cover
		raise RuntimeError(
			"Missing required dependency: speclite\n"
			"Example (conda): conda install -c conda-forge speclite\n"
			"Example (pip): pip install speclite"
		) from exc


def _parse_args(argv: list[str]) -> argparse.Namespace:
	p = argparse.ArgumentParser(
		prog="v3tk_to_VRI.py",
		description=(
			"Collapse a MUSE cube (DATA HDU) to observed V/R/I-band flux maps (nanomaggy) "
			"and AB magnitude maps, then write them to XXX_VRI.fits preserving spatial WCS."
		),
	)
	p.add_argument("fits_path", type=pathlib.Path, help="Input cube FITS path (e.g., XXX.fits)")
	p.add_argument(
		"--data-hdu",
		default="DATA",
		help="HDU name/index containing the cube data (default: DATA)",
	)
	p.add_argument(
		"--row-chunk",
		type=int,
		default=None,
		help=(
			"Number of spatial rows (Y) processed per chunk. "
			"Default: number of CPUs reported by os.cpu_count(). "
			"Lower this if you hit memory limits; raise it for speed if you have headroom."
		),
	)
	# Backward-compatible: --filter acts as an alias for --filter-r.
	p.add_argument(
		"--filter",
		dest="filter_r",
		default="bessell-R",
		help=argparse.SUPPRESS,
	)
	p.add_argument(
		"--filter-v",
		dest="filter_v",
		default="bessell-V",
		help="speclite filter name for V band (default: bessell-V)",
	)
	p.add_argument(
		"--filter-r",
		dest="filter_r",
		default="bessell-R",
		help="speclite filter name for R band (default: bessell-R)",
	)
	p.add_argument(
		"--filter-i",
		dest="filter_i",
		default="bessell-I",
		help="speclite filter name for I band (default: bessell-I)",
	)
	p.add_argument(
		"--flux-scale",
		type=float,
		default=1e-20,
		help=(
			"Scale factor applied to raw cube values to form F_lambda in erg/s/cm^2/Angstrom "
			"(default: 1e-20, typical for MUSE)")
		,
	)
	p.add_argument(
		"--output",
		type=pathlib.Path,
		default=None,
		help="Output FITS path. Default: input name with _VRI.fits",
	)
	p.add_argument(
		"--allow-partial-overlap",
		action="store_true",
		default=True,
		help=(
			"Allow integrating a filter even when the cube wavelength range does not fully cover the filter bandpass. "
			"This truncates/renormalizes the filter response to the cube range (so results are NOT standard Bessell V/R/I, "
			"but can be useful for visualization). Default: ON."
		),
	)
	p.add_argument(
		"--no-allow-partial-overlap",
		dest="allow_partial_overlap",
		action="store_false",
		help="Disable partial-overlap behavior (restore strict speclite behavior; may yield NaNs).",
	)
	p.add_argument(
		"--nan-policy",
		choices=("propagate", "zero"),
		default="zero",
		help=(
			"How to handle non-finite cube samples during filter integration. "
			"'propagate' keeps NaNs (can yield all-NaN band images if the cube has NaNs across that band). "
			"'zero' replaces non-finite cube samples with 0 before integrating (useful for visualization). "
			"Default: zero."
		),
	)
	return p.parse_args(argv)


def _spectral_wavelength_grid_aa(hdr, nz: int):
	import numpy as np
	import astropy.units as u
	from astropy.wcs import WCS

	w = WCS(hdr)
	spec_wcs = w.sub(["spectral"])  # 1D spectral WCS

	pix = np.arange(nz, dtype=float)

	# Prefer newer API when available.
	if hasattr(spec_wcs, "pixel_to_world_values"):
		wave_native = spec_wcs.pixel_to_world_values(pix)
	else:
		# all_pix2world wants shape (N, naxis) for naxis=1
		wave_native = spec_wcs.all_pix2world(pix[:, None], 0)[:, 0]

	cunit = None
	try:
		cunit = spec_wcs.wcs.cunit[0]
	except Exception:
		cunit = None

	if cunit is None or str(cunit).strip() == "":
		# Common for MUSE: spectral unit effectively Angstrom.
		wave = wave_native * u.AA
	else:
		wave = (wave_native * cunit).to(u.AA)

	return wave


def _spatial_wcs_header_2d(hdr):
	from astropy.wcs import WCS

	w2 = WCS(hdr).celestial
	# relax=True helps preserve non-standard WCS cards when possible.
	return w2.to_header(relax=True)


def _ensure_wave_increasing(
	wave_aa, data_3d
) -> Tuple["object", "object"]:
	# wave_aa is an astropy Quantity of shape (nz,)
	# data_3d is numpy array shape (nz, ny, nx)
	import numpy as np

	if wave_aa.shape[0] < 2:
		return wave_aa, data_3d

	if np.isfinite(wave_aa.value[0]) and np.isfinite(wave_aa.value[-1]) and (wave_aa.value[-1] < wave_aa.value[0]):
		return wave_aa[::-1], data_3d[::-1, :, :]
	return wave_aa, data_3d


def _wave_is_decreasing(wave_aa) -> bool:
	import numpy as np

	if wave_aa.shape[0] < 2:
		return False
	if np.isfinite(wave_aa.value[0]) and np.isfinite(wave_aa.value[-1]) and (wave_aa.value[-1] < wave_aa.value[0]):
		return True
	return False


def _filter_bandpass_range_aa(filt) -> tuple[float, float]:
	"""Return an approximate (min,max) wavelength range in Angstrom for a speclite filter.

	Uses the region where response > 0 to estimate bandpass coverage.
	"""
	import numpy as np
	import astropy.units as u

	w = getattr(filt, "wavelength", None)
	if w is None:
		raise ValueError("Filter has no wavelength attribute")
	# speclite may store wavelength as an astropy Quantity or as a plain ndarray.
	try:
		w_aa = np.asarray(w.to_value(u.AA), dtype=float)  # type: ignore[attr-defined]
	except Exception:
		w_aa = np.asarray(w, dtype=float)
	r = np.asarray(getattr(filt, "response", None), dtype=float)
	if r is None:
		raise ValueError("Filter has no response attribute")
	mask = np.isfinite(w) & np.isfinite(r) & (r > 0)
	if not np.any(mask):
		raise ValueError("Filter response has no positive samples")
	return float(w_aa[mask].min()), float(w_aa[mask].max())


def _truncate_filter_to_cube_range(filt, cube_min_aa: float, cube_max_aa: float, filters_mod):
	"""Return a new filter truncated to [cube_min_aa, cube_max_aa], or None if no overlap."""
	import numpy as np
	import astropy.units as u

	w = getattr(filt, "wavelength", None)
	if w is None:
		return None
	try:
		w_aa = np.asarray(w.to_value(u.AA), dtype=float)  # type: ignore[attr-defined]
	except Exception:
		w_aa = np.asarray(w, dtype=float)
	r = np.asarray(getattr(filt, "response", None), dtype=float)
	if r is None:
		return None
	mask = (
		np.isfinite(w_aa)
		& np.isfinite(r)
		& (r > 0)
		& (w_aa >= cube_min_aa)
		& (w_aa <= cube_max_aa)
	)
	if int(mask.sum()) < 2:
		return None

	FilterResponse = getattr(filters_mod, "FilterResponse", None)
	if FilterResponse is None:
		from speclite.filters import FilterResponse  # type: ignore

	name = getattr(filt, "name", None) or "FILTER"
	try:
		meta = dict(getattr(filt, "meta", {}) or {})
	except Exception:
		meta = {}
	# Record provenance without relying on FilterResponse(name=...), which is not supported
	# by some speclite versions.
	meta = meta or {}
	meta.setdefault("source_name", name)
	meta["truncated_to_cube"] = True
	meta["cube_min_aa"] = float(cube_min_aa)
	meta["cube_max_aa"] = float(cube_max_aa)

	# speclite versions differ in accepted kwargs; try the most common signatures.
	try:
		return FilterResponse(wavelength=w_aa[mask] * u.AA, response=r[mask], meta=meta)
	except TypeError:
		return FilterResponse(w_aa[mask] * u.AA, r[mask], meta)


def _print_coverage_report(wave_aa, filters_by_band: list[tuple[str, str, "object"]]) -> tuple[float, float]:
	"""Print cube spectral range and filter bandpass ranges, plus overlap fraction."""
	import numpy as np
	import astropy.units as u

	wave_vals = np.asarray(wave_aa.to_value(u.AA), dtype=float)
	finite = np.isfinite(wave_vals)
	if not np.any(finite):
		print("WARNING: cube wavelength grid has no finite values")
		return (float("nan"), float("nan"))

	cube_min = float(wave_vals[finite].min())
	cube_max = float(wave_vals[finite].max())
	if cube_max < cube_min:
		cube_min, cube_max = cube_max, cube_min

	print(f"v3tk_to_VRI.py version: {__version__}")
	print(f"Cube spectral coverage: {cube_min:.1f}–{cube_max:.1f} Å")
	for band, name, filt in filters_by_band:
		try:
			fmin, fmax = _filter_bandpass_range_aa(filt)
		except Exception as exc:
			print(f"  {band} ({name}): could not determine bandpass range: {exc}")
			continue
		overlap_min = max(cube_min, fmin)
		overlap_max = min(cube_max, fmax)
		span = max(0.0, fmax - fmin)
		overlap = max(0.0, overlap_max - overlap_min)
		frac = (overlap / span) if span > 0 else 0.0
		print(
			f"  {band} ({name}): bandpass {fmin:.1f}–{fmax:.1f} Å; "
			f"overlap {overlap_min:.1f}–{overlap_max:.1f} Å ({100.0 * frac:.1f}%)"
		)
		if overlap <= 0:
			print(f"  WARNING: {band}-band has no wavelength overlap with cube; output will be NaN")
	return (cube_min, cube_max)


def compute_vri_maps(
	fits_path: pathlib.Path,
	data_hdu: str,
	filter_v: str,
	filter_r: str,
	filter_i: str,
	flux_scale: float,
	row_chunk: int,
	allow_partial_overlap: bool,
	nan_policy: str,
) -> Tuple[
	"object",
	"object",
	"object",
	"object",
	"object",
	"object",
	"object",
	"object",
]:
	import numpy as np
	import astropy.units as u
	from astropy.io import fits
	from speclite import filters

	if row_chunk <= 0:
		raise ValueError("--row-chunk must be a positive integer")

	with fits.open(fits_path, memmap=True) as hdul:
		hdu = hdul[data_hdu]
		hdr = hdu.header
		data = hdu.data

		if data is None:
			raise ValueError(f"No data found in HDU '{data_hdu}'")
		if getattr(data, "ndim", None) != 3:
			raise ValueError(f"Expected 3D cube in HDU '{data_hdu}', got shape {getattr(data, 'shape', None)}")

		nz, ny, nx = data.shape
		wave = _spectral_wavelength_grid_aa(hdr, nz)
		reverse_spec = _wave_is_decreasing(wave)
		if reverse_spec:
			wave = wave[::-1]
		spec_slice = slice(None, None, -1) if reverse_spec else slice(None)

		primary_hdr = hdul[0].header.copy()
		spatial_hdr = _spatial_wcs_header_2d(hdr)

		v_flux_nmgy = np.full((ny, nx), np.nan, dtype=np.float32)
		v_mag_ab = np.full((ny, nx), np.nan, dtype=np.float32)
		r_flux_nmgy = np.full((ny, nx), np.nan, dtype=np.float32)
		r_mag_ab = np.full((ny, nx), np.nan, dtype=np.float32)
		i_flux_nmgy = np.full((ny, nx), np.nan, dtype=np.float32)
		i_mag_ab = np.full((ny, nx), np.nan, dtype=np.float32)

		f_v = filters.load_filter(filter_v)
		f_r = filters.load_filter(filter_r)
		f_i = filters.load_filter(filter_i)
		cube_min, cube_max = _print_coverage_report(wave, [("V", filter_v, f_v), ("R", filter_r, f_r), ("I", filter_i, f_i)])
		if allow_partial_overlap and (cube_min == cube_min) and (cube_max == cube_max):
			# Only truncate when the cube does NOT fully cover the filter bandpass.
			changed = False
			for band, filt_name, filt_obj in (("V", filter_v, f_v), ("R", filter_r, f_r), ("I", filter_i, f_i)):
				try:
					fmin, fmax = _filter_bandpass_range_aa(filt_obj)
				except Exception:
					# If we can't determine bandpass, skip truncation.
					continue
				fully_covered = (cube_min <= fmin) and (cube_max >= fmax)
				if fully_covered:
					continue
				new_filt = _truncate_filter_to_cube_range(filt_obj, cube_min, cube_max, filters_mod=filters)
				if new_filt is None:
					continue
				if band == "V":
					f_v = new_filt
				elif band == "R":
					f_r = new_filt
				else:
					f_i = new_filt
				changed = True
			if changed:
				print("NOTE: partial-overlap enabled: truncating/renormalizing filters to cube range where needed")
		scale_unit = flux_scale * u.erg / (u.s * u.cm**2 * u.AA)
		if nan_policy not in ("propagate", "zero"):
			raise ValueError("nan_policy must be 'propagate' or 'zero'")

		stats = {
			"V": {"finite": 0, "pos": 0},
			"R": {"finite": 0, "pos": 0},
			"I": {"finite": 0, "pos": 0},
		}

		for y0 in range(0, ny, row_chunk):
			y1 = min(ny, y0 + row_chunk)
			slab = np.asarray(data[spec_slice, y0:y1, :], dtype=np.float32)
			if nan_policy == "zero":
				slab = np.where(np.isfinite(slab), slab, 0.0).astype(np.float32, copy=False)
			f_lambda = slab * scale_unit
			f_lambda = np.moveaxis(f_lambda, 0, -1)  # (ychunk, nx, nz)

			for band, filt, flux_out, mag_out in (
				("V", f_v, v_flux_nmgy, v_mag_ab),
				("R", f_r, r_flux_nmgy, r_mag_ab),
				("I", f_i, i_flux_nmgy, i_mag_ab),
			):
				maggies = filt.get_ab_maggies(f_lambda, wavelength=wave)
				maggies = np.asarray(maggies, dtype=np.float64)  # (ychunk, nx)
				finite = np.isfinite(maggies)
				stats[band]["finite"] += int(finite.sum())
				pos = finite & (maggies > 0)
				stats[band]["pos"] += int(pos.sum())
				flux_out[y0:y1, :] = (maggies * 1e9).astype(np.float32)
				good = pos
				if np.any(good):
					mag_chunk = np.full(maggies.shape, np.nan, dtype=np.float32)
					mag_chunk[good] = (-2.5 * np.log10(maggies[good])).astype(np.float32)
					mag_out[y0:y1, :] = mag_chunk

		for band in ("V", "R", "I"):
			if stats[band]["finite"] == 0:
				print(f"WARNING: {band}-band produced 0 finite pixels (all NaN)")
			elif stats[band]["pos"] == 0:
				print(f"WARNING: {band}-band produced no positive pixels")

	return (
		primary_hdr,
		spatial_hdr,
		v_flux_nmgy,
		v_mag_ab,
		r_flux_nmgy,
		r_mag_ab,
		i_flux_nmgy,
		i_mag_ab,
	)


def write_output_fits(
	output_path: pathlib.Path,
	primary_hdr,
	spatial_hdr,
	v_flux_nmgy,
	v_mag_ab,
	r_flux_nmgy,
	r_mag_ab,
	i_flux_nmgy,
	i_mag_ab,
	filter_v: str,
	filter_r: str,
	filter_i: str,
):
	from astropy.io import fits

	hdu0 = fits.PrimaryHDU(header=primary_hdr)

	def _band_hdus(band: str, filt_name: str, flux, mag):
		h_flux = fits.ImageHDU(data=flux, header=spatial_hdr.copy(), name=f"{band}_FLUX")
		h_flux.header["BUNIT"] = ("nanomaggy", f"{band}-band flux density integrated via filter")
		h_flux.header["FILTER"] = (band, "Photometric band")
		h_flux.header["FILTNAM"] = (filt_name, "speclite filter used")
		h_flux.header["MAGZP"] = (22.5, "AB zeropoint for nanomaggy convention")

		h_mag = fits.ImageHDU(data=mag, header=spatial_hdr.copy(), name=f"{band}_MAG")
		h_mag.header["BUNIT"] = ("mag", "AB magnitude")
		h_mag.header["FILTER"] = (band, "Photometric band")
		h_mag.header["FILTNAM"] = (filt_name, "speclite filter used")
		return h_flux, h_mag

	v_hflux, v_hmag = _band_hdus("V", filter_v, v_flux_nmgy, v_mag_ab)
	r_hflux, r_hmag = _band_hdus("R", filter_r, r_flux_nmgy, r_mag_ab)
	i_hflux, i_hmag = _band_hdus("I", filter_i, i_flux_nmgy, i_mag_ab)

	hdul = fits.HDUList([hdu0, v_hflux, v_hmag, r_hflux, r_hmag, i_hflux, i_hmag])
	hdul.writeto(output_path, overwrite=True)


def main(argv: list[str]) -> int:
	t0 = time.perf_counter()
	try:
		_require_deps()
		args = _parse_args(argv)

		row_chunk = args.row_chunk
		if row_chunk is None:
			row_chunk = os.cpu_count() or 1
			row_chunk = max(1, int(row_chunk))
		print(f"row_chunk: {row_chunk}")

		fits_path = args.fits_path
		if not fits_path.exists():
			raise FileNotFoundError(f"Input FITS not found: {fits_path}")

		out = args.output
		if out is None:
			out = fits_path.with_name(f"{fits_path.stem}_VRI.fits")

		(
			primary_hdr,
			spatial_hdr,
			v_flux_nmgy,
			v_mag_ab,
			r_flux_nmgy,
			r_mag_ab,
			i_flux_nmgy,
			i_mag_ab,
		) = compute_vri_maps(
			fits_path=fits_path,
			data_hdu=args.data_hdu,
			filter_v=args.filter_v,
			filter_r=args.filter_r,
			filter_i=args.filter_i,
			flux_scale=args.flux_scale,
			row_chunk=row_chunk,
			allow_partial_overlap=args.allow_partial_overlap,
			nan_policy=args.nan_policy,
		)
		write_output_fits(
			out,
			primary_hdr,
			spatial_hdr,
			v_flux_nmgy,
			v_mag_ab,
			r_flux_nmgy,
			r_mag_ab,
			i_flux_nmgy,
			i_mag_ab,
			filter_v=args.filter_v,
			filter_r=args.filter_r,
			filter_i=args.filter_i,
		)
		print(f"Wrote: {out}")
		print("HDUs: V_FLUX/V_MAG, R_FLUX/R_MAG, I_FLUX/I_MAG")
		dt = time.perf_counter() - t0
		print(f"Runtime: {dt:.2f} s")
		return 0
	except Exception as exc:
		dt = time.perf_counter() - t0
		print(f"ERROR: {exc}", file=sys.stderr)
		print(f"Runtime: {dt:.2f} s", file=sys.stderr)
		return 2


if __name__ == "__main__":
	raise SystemExit(main(sys.argv[1:]))
