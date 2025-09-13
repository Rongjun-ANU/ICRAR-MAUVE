#!/usr/bin/env python
"""
Foreground_star.py – flag foreground stars in MAUVE kinematic maps.

Logic
-----
1. Read *_KIN_maps.fits for one galaxy.
2. Identify candidate foreground-star spaxels as 5-σ outliers in either
   velocity (V) or velocity dispersion (SIGMA).
3. Refine that mask with an iterative minimum-enclosing-circle algorithm.
4. Write *_KIN_maps_extended.fits with a new HDU called FOREGROUND_STAR.

Parameters
----------
max_radius     = 50   # pixels
max_iterations = 10
"""

import argparse, logging, sys, time
from pathlib import Path

import numpy as np
from astropy.io import fits

# ------------------------------------------------------------------
# 0.  Command-line interface
# ------------------------------------------------------------------
def cli():
    p = argparse.ArgumentParser(description="Flag foreground stars in KIN maps")
    p.add_argument("-g", "--galaxy", default="IC3392",
                   help="Galaxy identifier (default IC3392)")
    p.add_argument("--root", default="/arc/projects/mauve",
                   help="Root of MAUVE directory tree")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Verbose logging")
    return p.parse_args()

args   = cli()
loglvl = logging.INFO if args.verbose else logging.WARNING
logging.basicConfig(level=loglvl,
                    format="%(asctime)s  %(levelname)-8s  %(message)s",
                    datefmt="%H:%M:%S")
t0     = time.perf_counter()

galaxy = args.galaxy.upper()
root   = Path(args.root).expanduser().resolve()
kin_path = root / "products/v0.6" / galaxy / f"{galaxy}_KIN_maps.fits"
out_path = Path(f"{galaxy}_KIN_maps_extended.fits")

logging.info("Loading kinematic map from %s", kin_path)
if not kin_path.is_file():
    logging.error("File not found: %s", kin_path)
    sys.exit(1)

# ------------------------------------------------------------------
# 1.  Load every extension (except PRIMARY) into globals()
# ------------------------------------------------------------------
with fits.open(kin_path) as hdul:
    ext_names = [hdul[i].name for i in range(1, len(hdul))]
    for name in ext_names:
        globals()[name] = hdul[name].data
        logging.debug("Loaded %-12s %s", name, globals()[name].shape)
    kin_header = hdul[0].header.copy()       # keep WCS metadata

# ------------------------------------------------------------------
# 2.  Initial foreground-star mask: 5-σ outliers in V or SIGMA
# ------------------------------------------------------------------
V, SIGMA = globals().get("V"), globals().get("SIGMA")
if V is None or SIGMA is None:
    logging.error("V and/or SIGMA layers missing in %s", kin_path)
    sys.exit(1)

foreground_idx = np.where(
    (np.abs(V - np.nanmedian(V)) >= 5 * np.nanstd(V)) |
    (np.abs(SIGMA - np.nanmedian(SIGMA)) >= 5 * np.nanstd(SIGMA))
)
foreground_mask = np.zeros_like(V, dtype=bool)
foreground_mask[foreground_idx] = True
logging.info("Initial 5-σ selection: %d candidate spaxels",
             foreground_mask.sum())

# ------------------------------------------------------------------
# 3.  Iterative circle-refinement helper
# ------------------------------------------------------------------
def refine_mask(mask, max_radius=50, max_iterations=10):
    current = mask.copy()
    for it in range(max_iterations):
        ys, xs = np.where(current)
        if xs.size == 0:
            logging.info("No spaxels remain after %d iterations", it)
            return current

        cx, cy = xs.mean(), ys.mean()             # provisional centre
        d = np.hypot(xs - cx, ys - cy)
        drop = d > 2 * max_radius                 # beyond 2×radius

        if not np.any(drop):
            logging.debug("Converged in %d iterations", it)
            break

        current[ys[drop], xs[drop]] = False
        logging.debug("Iter %d: removed %d spaxels", it + 1, drop.sum())

    ys, xs = np.where(current)
    if xs.size == 0:
        return current

    cx, cy = xs.mean(), ys.mean()
    final_r = np.max(np.hypot(xs - cx, ys - cy))
    logging.info("Final circle centre (x,y) = (%.1f, %.1f); radius = %.2f",
                 cx, cy, final_r)

    ny, nx = mask.shape
    y_idx, x_idx = np.indices((ny, nx))          # orientation-safe grids
    return np.hypot(x_idx - cx, y_idx - cy) <= final_r

foreground_mask = refine_mask(foreground_mask)

# ---- YOUR exact NaN-preserving line ------------------------------
foreground_mask = np.where(np.isfinite(V), foreground_mask, np.nan)

logging.info("Refined mask: %d spaxels flagged",
             np.nansum(foreground_mask).astype(int))

# ------------------------------------------------------------------
# 4.  Copy original HDUs, append mask, and write to disk
# ------------------------------------------------------------------
with fits.open(kin_path) as hdul_in:
    hdul_out = fits.HDUList([hdu.copy() for hdu in hdul_in])

# Preserve NaNs → store as float so they survive in FITS
hdul_out.append(
    fits.ImageHDU(foreground_mask.astype(np.float32),
                  header=kin_header, name="FOREGROUND_STAR")
)

hdul_out.writeto(out_path, overwrite=True)
logging.info("Extended file written ➜ %s", out_path.resolve())
logging.info("Total runtime: %.1f s", time.perf_counter() - t0)
