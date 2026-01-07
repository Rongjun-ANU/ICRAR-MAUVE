#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import pathlib
import sys
import time
from typing import Tuple


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
		prog="v3tk_to_R.py",
		description=(
			"Collapse a MUSE cube (DATA HDU) to an observed R-band flux map (nanomaggy) "
			"and AB magnitude map, then write them to XXX_R.fits preserving spatial WCS."
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
	p.add_argument(
		"--filter",
		dest="filter_name",
		default="bessell-R",
		help=(
			"speclite filter name (default: bessell-R). Examples: sdss2010-r, sdss2010noatm-r, "
			"decamDR1noatm-r"
		),
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
		help="Output FITS path. Default: input name with _R.fits",
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


def compute_rband_maps(
	fits_path: pathlib.Path,
	data_hdu: str,
	filter_name: str,
	flux_scale: float,
	row_chunk: int,
) -> Tuple["object", "object", "object", "object"]:
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

		flux_nmgy = np.full((ny, nx), np.nan, dtype=np.float32)
		mag_ab = np.full((ny, nx), np.nan, dtype=np.float32)

		f_r = filters.load_filter(filter_name)
		scale_unit = flux_scale * u.erg / (u.s * u.cm**2 * u.AA)

		for y0 in range(0, ny, row_chunk):
			y1 = min(ny, y0 + row_chunk)
			slab = np.asarray(data[spec_slice, y0:y1, :], dtype=np.float32)
			f_lambda = slab * scale_unit
			f_lambda = np.moveaxis(f_lambda, 0, -1)  # (ychunk, nx, nz)
			maggies = f_r.get_ab_maggies(f_lambda, wavelength=wave)
			maggies = np.asarray(maggies, dtype=np.float64)  # (ychunk, nx)

			flux_nmgy[y0:y1, :] = (maggies * 1e9).astype(np.float32)
			good = np.isfinite(maggies) & (maggies > 0)
			if np.any(good):
				mag_chunk = np.full(maggies.shape, np.nan, dtype=np.float32)
				mag_chunk[good] = (-2.5 * np.log10(maggies[good])).astype(np.float32)
				mag_ab[y0:y1, :] = mag_chunk

	return primary_hdr, spatial_hdr, flux_nmgy, mag_ab


def write_output_fits(
	output_path: pathlib.Path,
	primary_hdr,
	spatial_hdr,
	flux_nmgy,
	mag_ab,
):
	from astropy.io import fits

	hdu0 = fits.PrimaryHDU(header=primary_hdr)

	h_flux = fits.ImageHDU(data=flux_nmgy, header=spatial_hdr.copy(), name="R_FLUX")
	h_flux.header["BUNIT"] = ("nanomaggy", "R-band flux density integrated via filter")
	h_flux.header["FILTER"] = ("R", "Photometric band")
	h_flux.header["MAGZP"] = (22.5, "AB zeropoint for nanomaggy convention")
	# Note: AB mags still computed as -2.5 log10(maggies)

	h_mag = fits.ImageHDU(data=mag_ab, header=spatial_hdr.copy(), name="R_MAG")
	h_mag.header["BUNIT"] = ("mag", "AB magnitude")
	h_mag.header["FILTER"] = ("R", "Photometric band")

	hdul = fits.HDUList([hdu0, h_flux, h_mag])
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
			out = fits_path.with_name(f"{fits_path.stem}_R.fits")

		primary_hdr, spatial_hdr, flux_nmgy, mag_ab = compute_rband_maps(
			fits_path=fits_path,
			data_hdu=args.data_hdu,
			filter_name=args.filter_name,
			flux_scale=args.flux_scale,
			row_chunk=row_chunk,
		)
		write_output_fits(out, primary_hdr, spatial_hdr, flux_nmgy, mag_ab)
		print(f"Wrote: {out}")
		print("HDU1: R_FLUX (nanomaggy), HDU2: R_MAG (AB mag)")
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
