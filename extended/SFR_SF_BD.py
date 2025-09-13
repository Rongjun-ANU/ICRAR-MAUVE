#!/usr/bin/env python
"""
SFR.py – produce Balmer-decrement, dust-corrected Hα maps and SFR
surface-density maps for a MAUVE galaxy, following exactly the
command-line and file-path style of Mass.py.

Changes (2025-06-13)
--------------------
* No initial S/N cut: `cut = 0`
* Star-forming mask now further restricted to spaxels with BD > 2.5
* Negative colour excess values are forced to zero:
      E_BV_BD[E_BV_BD < 0] = 0.0
"""

import argparse, logging, time, sys
from pathlib import Path

# ------------------------------------------------------------------
# 0.  Command-line interface
# ------------------------------------------------------------------
def cli():
    p = argparse.ArgumentParser(description="Generate SFR maps for a MAUVE galaxy")
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

# ------------------------------------------------------------------
# 1.  Resolve file paths & sanity checks
# ------------------------------------------------------------------
gas_path = root / "products/v0.6" / galaxy / f"{galaxy}_gas_BIN_maps.fits"
SFH_path = root / "products/v0.6" / galaxy / f"{galaxy}_SFH_maps.fits"
out_path = Path(f"{galaxy}_gas_BIN_maps_extended.fits")

print("\n=== Using the following files ===")
print("Gas lines    :", gas_path)
print("SFH map      :", SFH_path)
print("Output       :", out_path, "\n")

missing = [p for p in (gas_path, SFH_path) if not p.is_file()]
if missing:
    for p in missing:
        logging.error("Missing file: %s", p)
    sys.exit(1)

# ------------------------------------------------------------------
# 2.  Notebook-derived code starts
# ------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from urllib import request
from scipy.interpolate import interp1d
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import speclite.filters as sp
from speclite import filters
from scipy.ndimage import sum_labels, mean

from astropy.io import fits
from astropy import units as u
from astropy import constants as c
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

from ppxf.ppxf import ppxf, rebin
import ppxf.ppxf_util as util
from ppxf import sps_util as lib

# ------------------------------------------------------------------
# 2.1  Load emission-line maps
# ------------------------------------------------------------------
print(f"Loading gas line map from {gas_path}")
with fits.open(gas_path) as hdul:
    HB4861_FLUX      = hdul['HB4861_FLUX'].data
    HB4861_FLUX_ERR  = hdul['HB4861_FLUX_ERR'].data
    HA6562_FLUX      = hdul['HA6562_FLUX'].data
    HA6562_FLUX_ERR  = hdul['HA6562_FLUX_ERR'].data
    OIII5006_FLUX    = hdul['OIII5006_FLUX'].data
    OIII5006_FLUX_ERR= hdul['OIII5006_FLUX_ERR'].data
    NII6583_FLUX     = hdul['NII6583_FLUX'].data
    NII6583_FLUX_ERR = hdul['NII6583_FLUX_ERR'].data
    SII6716_FLUX     = hdul['SII6716_FLUX'].data
    SII6716_FLUX_ERR = hdul['SII6716_FLUX_ERR'].data
    SII6730_FLUX     = hdul['SII6730_FLUX'].data
    SII6730_FLUX_ERR = hdul['SII6730_FLUX_ERR'].data
    gas_header       = hdul[0].header.copy()          # keep WCS for new layers

# ------------------------------------------------------------------
# 3.  (No) S/N cut – set threshold to zero
# ------------------------------------------------------------------
cut = 0  # ← new default: accept every spaxel, even low-S/N ones

mask_HB = HB4861_FLUX / HB4861_FLUX_ERR >= cut
mask_HA = HA6562_FLUX / HA6562_FLUX_ERR >= cut
mask_combined = mask_HB & mask_HA           # practically all True now

HB4861_FLUX_cut = np.where(mask_combined, HB4861_FLUX, np.nan)
HA6562_FLUX_cut = np.where(mask_combined, HA6562_FLUX, np.nan)

# ------------------------------------------------------------------
# 4.  Balmer Decrement
# ------------------------------------------------------------------
BD = HA6562_FLUX_cut / HB4861_FLUX_cut

# ------------------------------------------------------------------
# 5.  BPT masks (unchanged logic)
# ------------------------------------------------------------------
logN2 = np.log10(NII6583_FLUX / HA6562_FLUX)
logS2 = np.log10((SII6716_FLUX + SII6730_FLUX) / HA6562_FLUX)
logO3 = np.log10(OIII5006_FLUX / HB4861_FLUX)

mask_N2 = NII6583_FLUX / NII6583_FLUX_ERR >= cut
mask_S2 = (SII6716_FLUX + SII6730_FLUX) / (SII6716_FLUX_ERR + SII6730_FLUX_ERR) >= cut
mask_O3 = OIII5006_FLUX / OIII5006_FLUX_ERR >= cut
mask_combinedd = mask_combined & mask_N2 & mask_S2 & mask_O3

#  BPT demarcations
def kewley01_N2(x):   return 0.61 / (x - 0.47) + 1.19
def kauff03_N2(x):    return 0.61 / (x - 0.05) + 1.30
def kewley01_S2(x):   return 0.72 / (x - 0.32) + 1.30
def kewley06_Sy_LIN(x): return 1.89 * x + 0.76

#  Classification
N2_HII   = logO3 <= kauff03_N2(logN2)
N2_Comp  = (logO3 > kauff03_N2(logN2)) & (logO3 <= kewley01_N2(logN2))
S2_HII   = logO3 <= kewley01_S2(logS2)

mask_SF  = (N2_HII | N2_Comp) & S2_HII    # star-forming in BOTH diagrams

# ------------------------------------------------------------------
# 6.  Additional “good” criterion: BD > 2.5
# ------------------------------------------------------------------
mask_good = mask_SF & (BD > 2.5)

# ------------------------------------------------------------------
# 7.  Reddening from BD, with negative values floored at 0
# ------------------------------------------------------------------
# Calzetti (2000) coefficients
k_Hb, k_Ha, R_int = 4.598, 3.325, 2.86
E_BV_BD = 2.5 / (k_Hb - k_Ha) * np.log10( (HA6562_FLUX_cut / HB4861_FLUX_cut) / R_int )
E_BV_BD[E_BV_BD < 0] = 0.0    # ← new behaviour

# ------------------------------------------------------------------
# 8.  Dust-corrected Hα (good SF spaxels only)
# ------------------------------------------------------------------
HA6562_FLUX_SF         = np.where(mask_good, HA6562_FLUX_cut, np.nan)
HA6562_FLUX_SF_corr    = HA6562_FLUX_SF * 10**(0.4 * k_Ha * E_BV_BD)

# ------------------------------------------------------------------
# 9.  Flux → luminosity → SFR
# ------------------------------------------------------------------
def flux_to_luminosity(flux, distance=16.5):          # distance in Mpc
    return (flux * 1e-20 * u.erg / u.s / u.cm**2 * 4 * np.pi * (distance * u.Mpc)**2).cgs

def calzetti_sfr(lum):                                # Calzetti ’07 + IMF tweak
    return 5.3e-42 * lum.value / 0.67 * 0.63          # M⊙ yr⁻¹

HA6562_L_SF_corr       = flux_to_luminosity(HA6562_FLUX_SF_corr)
total_HA6562_L_SF_corr = np.nansum(HA6562_L_SF_corr)

HA6562_SFR_SF_corr     = calzetti_sfr(HA6562_L_SF_corr)
total_HA6562_SFR_SF_corr = calzetti_sfr(total_HA6562_L_SF_corr)

print(f"Total SFR (good, dust-corrected): "
      f"{total_HA6562_SFR_SF_corr:.2f} M_sun/yr "
      f"(log = {np.log10(total_HA6562_SFR_SF_corr):.2f})")

# ------------------------------------------------------------------
# 10.  Σ_SFR (surface density)
# ------------------------------------------------------------------
legacy_wcs2    = WCS(gas_header).celestial
pix_scale      = (proj_plane_pixel_scales(legacy_wcs2) * u.deg).to(u.arcsec)
pix_area_kpc2  = ((pix_scale[0].to(u.rad) * 16.5 * u.Mpc) *
                  (pix_scale[1].to(u.rad) * 16.5 * u.Mpc)).to(u.kpc**2)

SFR_surface_density     = HA6562_SFR_SF_corr / pix_area_kpc2   # M⊙ yr⁻¹ kpc⁻²
log_SFR_surface_density = np.log10(SFR_surface_density.value)

# ------------------------------------------------------------------
# 11.  Write new layers
# ------------------------------------------------------------------
with fits.open(gas_path) as hdul:
    new_hdul = fits.HDUList([hdu.copy() for hdu in hdul])

# Balmer Decrement
new_hdul.append(fits.ImageHDU(BD.astype(np.float64),
                              header=gas_header, name="Balmer_Decrement"))
# E(B-V)_BD
new_hdul.append(fits.ImageHDU(E_BV_BD.astype(np.float64),
                              header=gas_header, name="E(B-V)_BD"))
# Corrected Hα luminosity
new_hdul.append(fits.ImageHDU(HA6562_L_SF_corr.value.astype(np.float64),
                              header=gas_header, name="Halpha_Luminosity_SF_corr",
                              bunit="erg/s"))
# Corrected SFR
new_hdul.append(fits.ImageHDU(HA6562_SFR_SF_corr.astype(np.float64),
                              header=gas_header, name="Halpha_SFR_SF_corr",
                              bunit="M_sun/yr"))
# log Σ_SFR
new_hdul.append(fits.ImageHDU(log_SFR_surface_density.astype(np.float64),
                              header=gas_header, name="LOGSFR_SURFACE_DENSITY",
                              bunit="log(M_sun/yr/kpc2)"))

new_hdul.writeto(out_path, overwrite=True)
print("Extended file written ➜", out_path.resolve())

print(f"Total runtime: {time.perf_counter() - t0:.1f} s")
