
#!/usr/bin/env python
"""
SFR.py – produce Balmer‐decrement, dust‐corrected Hα maps and SFR surface–density maps
for a MAUVE galaxy, following exactly the command‑line and file‑path style of Mass.py.

Usage
-----
    python SFR.py -g IC3392 [--root /arc/projects/mauve]

Outputs
-------
{GAL}_gas_BIN_maps_extended.fits  (existing layers retained, new SFR‑related layers appended)
"""

import argparse, logging, time, sys
from pathlib import Path

def cli():
    p = argparse.ArgumentParser(description="Generate SFR maps for a MAUVE galaxy")
    p.add_argument("-g", "--galaxy", default="IC3392",
                   help="Galaxy identifier (default IC3392)")
    p.add_argument("--root", default="/arc/projects/mauve",
                   help="Root of MAUVE directory tree")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Verbose logging")
    return p.parse_args()

args = cli()
loglvl = logging.INFO if args.verbose else logging.WARNING
logging.basicConfig(level=loglvl, format="%(asctime)s  %(levelname)-8s  %(message)s",
                    datefmt="%H:%M:%S")

t0 = time.perf_counter()
galaxy = args.galaxy.upper()
root = Path(args.root).expanduser().resolve()

# Resolve file paths
gas_path  = root / "products/v0.6" / galaxy / f"{galaxy}_gas_BIN_maps.fits"
SFH_path  = root / "products/v0.6" / galaxy / f"{galaxy}_SFH_maps.fits"
out_path  = Path(f"{galaxy}_gas_BIN_maps_extended.fits")

print("\n=== Using the following files ===")
print("Gas lines    :", gas_path)
print("SFH map      :", SFH_path)
print("Output       :", out_path, "\n")

missing = [p for p in (gas_path, SFH_path) if not p.is_file()]
if missing:
    for p in missing:
        logging.error("Missing file: %s", p)
    sys.exit(1)

# ---- notebook‑derived code starts ----

from pathlib import Path
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

import os
import sys
import glob


# Load gas line map IC3392_gas_BIN_maps.fits 
print(f"Loading gas line map from {gas_path}")
with fits.open(gas_path) as hdul:
    HB4861_FLUX = hdul['HB4861_FLUX'].data
    HB4861_FLUX_ERR = hdul['HB4861_FLUX_ERR'].data
    HA6562_FLUX = hdul['HA6562_FLUX'].data
    HA6562_FLUX_ERR = hdul['HA6562_FLUX_ERR'].data
    OIII5006_FLUX = hdul['OIII5006_FLUX'].data
    OIII5006_FLUX_ERR = hdul['OIII5006_FLUX_ERR'].data
    NII6583_FLUX = hdul['NII6583_FLUX'].data
    NII6583_FLUX_ERR = hdul['NII6583_FLUX_ERR'].data
    SII6716_FLUX = hdul['SII6716_FLUX'].data
    SII6716_FLUX_ERR = hdul['SII6716_FLUX_ERR'].data
    SII6730_FLUX = hdul['SII6730_FLUX'].data
    SII6730_FLUX_ERR = hdul['SII6730_FLUX_ERR'].data
    gas_header = hdul['HB4861_FLUX'].header  # Use header from one of the lines
    hdul.close()

# Apply a first cut of FLUX/ERR ≥ 5 on every line, 
# at least 22 to get min BD>2.86
# at least 15 to have all SF in BPT diagram
cut = 15
mask_HB = HB4861_FLUX / HB4861_FLUX_ERR >= cut
mask_HA = HA6562_FLUX / HA6562_FLUX_ERR >= cut
# Combine masks for both lines
mask_combined = mask_HB & mask_HA
# Apply the mask to the flux maps
HB4861_FLUX_cut = np.where(mask_combined, HB4861_FLUX, np.nan)
HA6562_FLUX_cut = np.where(mask_combined, HA6562_FLUX, np.nan)

# Balmer-decrement map – H α/H β (start with all spaxels).
BD = HA6562_FLUX_cut / HB4861_FLUX_cut

# Print the lowest and highest values of the Balmer decrement map
print(f"Lowest Balmer Decrement: {np.nanmin(BD)}")
print(f"Highest Balmer Decrement: {np.nanmax(BD)}")
# Print lowest 5 non-NaN values (unique values only)
bd_valid = BD[~np.isnan(BD)]
bd_unique = np.unique(bd_valid)
lowest_5_unique = bd_unique[:5]
print(f"Lowest 5 unique non-NaN Balmer Decrement values: {lowest_5_unique}")



# 1) read the whole HDUList from disk
with fits.open(gas_path) as hdul:
    # 2) clone all existing HDUs into a new list
    new_hdul = fits.HDUList([hdu.copy() for hdu in hdul])

# 3) build the new Balmer Decrement image HDU
BD_hdu = fits.ImageHDU(
    data=BD.astype(np.float64), name="Balmer_Decrement",)

# 4) keep WCS and pixel-scale info by copying the original BINID header
BD_hdu.header.update(gas_header)         # you created gas_hdr earlier
BD_hdu.header["EXTNAME"] = "Balmer_Decrement"
BD_hdu.header["BUNIT"]   = ""  # physical units

# 5) append and write to disk
new_hdul.append(BD_hdu)
new_hdul.writeto(out_path, overwrite=True)

print(f"Saved extended file to {out_path.resolve()}")

# O'Donnell (1994) update
# k_Hb = 3.609   # 4861 Å
# k_Ha = 2.535   # 6563 Å
# Calzetti (2000)
k_Hb = 4.598 
k_Ha = 3.325
R_int = 2.86

E_BV_BD = 2.5/(k_Hb - k_Ha) * np.log10( (HA6562_FLUX_cut/HB4861_FLUX_cut) / R_int )

with fits.open(out_path, mode="append") as hdul:             # open existing file
    E_BV_BD_hdu = fits.ImageHDU(
        data=E_BV_BD.astype(np.float64),  # like others
        header=gas_header, name="E(B-V)_BD")  # new HDU name
    E_BV_BD_hdu.header["BUNIT"] = ""  # units keyword
    hdul.append(E_BV_BD_hdu)                                    # add as new HDU 
    hdul.flush()                                             # write in-place
print("E(B-V)_BD layer saved ➜", out_path.resolve())

# ---- line ratios --------------------------------------------------
logN2  = np.log10(NII6583_FLUX / HA6562_FLUX)        # [N II]/Hα
logS2  = np.log10((SII6716_FLUX+SII6730_FLUX) / HA6562_FLUX)   # Σ[S II]/Hα
logO3  = np.log10(OIII5006_FLUX / HB4861_FLUX)       # [O III]/Hβ         

mask_N2 = NII6583_FLUX / NII6583_FLUX_ERR >= cut
mask_S2 = (SII6716_FLUX + SII6730_FLUX) / (SII6716_FLUX_ERR + SII6730_FLUX_ERR) >= cut
mask_O3 = OIII5006_FLUX / OIII5006_FLUX_ERR >= cut
mask_combinedd = mask_combined & mask_N2 & mask_S2 & mask_O3

logN2_cut = np.where(mask_combinedd, logN2, np.nan)
logS2_cut = np.where(mask_combinedd, logS2, np.nan)
logO3_cut = np.where(mask_combinedd, logO3, np.nan)

#  N II BPT -----------------------------------------
def kewley01_N2(x):   # max-starburst
    return 0.61/(x-0.47) + 1.19
def kauff03_N2(x):    # empirical SF upper envelope
    return 0.61/(x-0.05) + 1.30                            

#  S II BPT -----------------------------------------
def kewley01_S2(x):
    return 0.72/(x-0.32) + 1.30                           
def kewley06_Sy_LIN(x):   # Seyfert/LINER division
    return 1.89*x + 0.76        

# Count the number of spaxels in each region
N2_HII = logO3 <= kauff03_N2(logN2)
N2_Comp = (logO3 > kauff03_N2(logN2)) & (logO3 <= kewley01_N2(logN2))
N2_AGN = logO3 > kewley01_N2(logN2)
S2_HII = logO3 <= kewley01_S2(logS2)
S2_Seyfert = (logO3 > kewley01_S2(logS2)) & (logO3 > kewley06_Sy_LIN(logS2))
S2_LINER = (logO3 > kewley01_S2(logS2)) & (logO3 <= kewley06_Sy_LIN(logS2))
# Count the number of spaxels in each region
N2_HII_count = np.sum(N2_HII)
N2_Comp_count = np.sum(N2_Comp)
N2_AGN_count = np.sum(N2_AGN)
S2_HII_count = np.sum(S2_HII)
S2_Seyfert_count = np.sum(S2_Seyfert)
S2_LINER_count = np.sum(S2_LINER)
print(f"Number of spaxels in [N II] BPT regions:")
print(f"HII: {N2_HII_count}, Comp: {N2_Comp_count}, AGN: {N2_AGN_count}")
print(f"Number of spaxels in [S II] BPT regions:")
print(f"HII: {S2_HII_count}, Seyfert: {S2_Seyfert_count}, LINER: {S2_LINER_count}")

# Purely SF spaxels in both BPT diagram
mask_SF = (N2_HII+N2_Comp) & (S2_HII)
# Apply the mask to the Halpha flux map
HA6562_FLUX_SF = np.where(mask_SF, HA6562_FLUX_cut, np.nan)

# Corrected Halpha map with E(B-V) from gas lines
HA6562_FLUX_corr = HA6562_FLUX_cut * 10**(0.4 * k_Ha * E_BV_BD)

HA6562_FLUX_SF_corr = HA6562_FLUX_SF * 10**(0.4 * k_Ha * E_BV_BD)

# Convert the corrected Halpha map ($10^{-20}erg/(s cm^2)$) to luminosity (erg/s)
def flux_to_luminosity(flux, distance=16.5):
    """
    Convert flux to luminosity.
    
    Parameters:
    flux : array-like
        Flux in erg/(s * Angstrom * cm^2).
    distance : float
        Distance in parsecs.
        
    Returns:
    luminosity : array-like
        Luminosity in erg/s.
    """
    return (flux*1e-20*u.erg/u.s/u.cm**2 * 4*np.pi*(distance*u.Mpc)**2).cgs

HA6562_Luminosity_corr = flux_to_luminosity(HA6562_FLUX_corr)
HA6562_Luminosity = flux_to_luminosity(HA6562_FLUX)
HA6562_Luminosity_SF_corr = flux_to_luminosity(HA6562_FLUX_SF_corr)
total_HA6562_Luminosity_corr = np.nansum(HA6562_Luminosity_corr)
total_HA6562_Luminosity = np.nansum(HA6562_Luminosity)
total_HA6562_Luminosity_SF_corr = np.nansum(HA6562_Luminosity_SF_corr)
log_halpha_corr = np.log10(HA6562_Luminosity_corr/c.L_sun.cgs)
total_log_halpha_corr = np.log10(total_HA6562_Luminosity_corr/c.L_sun.cgs)
total_log_halpha_SF_corr = np.log10(total_HA6562_Luminosity_SF_corr/c.L_sun.cgs)
print(f"Total corrected Hα Luminosity: {total_HA6562_Luminosity_corr:.2e} erg/s")
print(f"Total Hα Luminosity: {total_HA6562_Luminosity:.2e} ")
print(f"Total corrected Hα Luminosity for purely SF spaxels: {total_HA6562_Luminosity_SF_corr:.2e} ")
print(f"Total corrected log Hα Luminosity: {total_log_halpha_corr:.2f} L_sun")
print(f"Total corrected log Hα Luminosity for purely SF spaxels: {total_log_halpha_SF_corr:.2f} L_sun")

with fits.open(out_path, mode="append") as hdul:             # open existing file
    HA6562_Luminosity_SF_corr_hdu = fits.ImageHDU(
        data=HA6562_Luminosity_SF_corr.value.astype(np.float64),  # like others
        header=gas_header, name="Halpha_Luminosity_SF_corr")  # new HDU name
    HA6562_Luminosity_SF_corr_hdu.header["BUNIT"] = "erg/s"  # units keyword
    hdul.append(HA6562_Luminosity_SF_corr_hdu)                                    
    hdul.flush()                                             # write in-place
print("Corrected Halpha luminosity from SF layer saved ➜", out_path.resolve())


# SFR map from Halpha luminosity, using Calzetti 2007
def calzetti_sfr(luminosity):
    """
    Convert Halpha luminosity to SFR using Calzetti 2007.
    But it is assuming the Kroupa IMF, 
    so we need to times a coefficient to go to Chabrier IMF.
    
    Parameters:
    luminosity : array-like
        Halpha luminosity in erg/s.
        
    Returns:
    sfr : array-like
        Star formation rate in solar masses per year.
    """
    return 5.3e-42 * luminosity.cgs.value / 0.67 * 0.63 # SFR in M_sun/yr
HA6562_SFR_corr = calzetti_sfr(HA6562_Luminosity_corr)
total_HA6562_SFR_corr = calzetti_sfr(total_HA6562_Luminosity_corr)
HA6562_SFR_SF_corr = calzetti_sfr(HA6562_Luminosity_SF_corr)
total_HA6562_SFR_SF_corr = calzetti_sfr(total_HA6562_Luminosity_SF_corr)

# Print the total SFR in the map
print(f"Total SFR from corrected Hα Luminosity: {total_HA6562_SFR_corr:.2f} M_sun/yr or log(SFR) = {np.log10(total_HA6562_SFR_corr):.2f} M_sun/yr")
print(f"Total SFR from corrected Hα Luminosity for purely SF spaxels: {total_HA6562_SFR_SF_corr:.2f} M_sun/yr or log(SFR) = {np.log10(total_HA6562_SFR_SF_corr):.2f} M_sun/yr")


with fits.open(out_path, mode="append") as hdul:             # open existing file
    HA6562_SFR_SF_corr_hdu = fits.ImageHDU(
        data=HA6562_SFR_SF_corr.astype(np.float64),  # like others
        header=gas_header, name="Halpha_SFR_SF_corr")  # new HDU name
    HA6562_SFR_SF_corr_hdu.header["BUNIT"] = "M_sun/yr"  # units keyword
    hdul.append(HA6562_SFR_SF_corr_hdu)                                    
    hdul.flush()                                             # write in-place
print("Corrected Halpha SFR from SF layer saved ➜", out_path.resolve())


# Getting the SFR surface density 
# Convert to surface density in M☉/pc²
# 1. Convert pixel area to physical area in pc²
legacy_wcs2 = WCS(gas_header).celestial  # strip spectral axis
pixel_scale = (proj_plane_pixel_scales(legacy_wcs2) * u.deg).to(u.arcsec)
pixel_area_Mpc = ((pixel_scale[0]).to(u.rad).value*16.5*u.Mpc)*(((pixel_scale[1]).to(u.rad).value*16.5*u.Mpc))
pixel_area_kpc = pixel_area_Mpc.to(u.kpc**2)
# 2. Convert SFR to surface density
SFR_surface_density = HA6562_SFR_SF_corr / pixel_area_kpc  # M☉/yr/kpc²
# 3. Convert to log10 scale
log_SFR_surface_density = np.log10(SFR_surface_density.value)


with fits.open(out_path, mode="append") as hdul:             # open existing file
    log_SFR_surface_density_hdu = fits.ImageHDU(
        data=log_SFR_surface_density.astype(np.float64),  # like others
        header=gas_header, name="LOGSFR_SURFACE_DENSITY")  # new HDU name
    log_SFR_surface_density_hdu.header["BUNIT"] = "log(M_sun/yr/kpc2)"  # units keyword
    hdul.append(log_SFR_surface_density_hdu)                                    
    hdul.flush()                                             # write in-place
print("Corrected Halpha SFR surface density from SF layer saved ➜", out_path.resolve())


print(f'Total runtime: {time.perf_counter() - t0:.1f} s')
