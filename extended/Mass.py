#!/usr/bin/env python
"""
Generate extended Voronoi‑binning maps for a MAUVE galaxy.

This script is adapted from `map_MILES_brief.py` and now takes the galaxy
identifier (e.g. IC3392) as a command‑line argument:

    python map_MILES_brief_param.py -g IC3392

The MAUVE directory tree is assumed to look like::

    /arc/projects/mauve/
      ├─ cubes/v3.0/{galaxy}_DATACUBE_FINAL_WCS_Pall_mad_red_v3.fits
      └─ products/v0.6/{galaxy}/
             ├─ {galaxy}_SPATIAL_BINNING_maps.fits
             ├─ {galaxy}_SFH_maps.fits
             └─ {galaxy}_sfh-weights.fits

Pass a different root path with ``--root`` if your tree lives elsewhere.

Changes (2025-09-14)
-----------------------
* Added inclination correction for stellar mass surface density calculation.
* Implemented read_galaxy_inclination() function to read inclination angles from MAUVE_Inclination.dat.
* Applied cos(θ) correction factor to stellar_mass_surface_density where θ is the galaxy inclination angle.
* Enhanced logging to show inclination values and correction factors applied.

Changes (2025-09-15)
-----------------------
* Added user-configurable inclination correction parameter.
* Users can now enable/disable inclination correction by setting apply_inclination_correction = True/False.

Changes (2025-09-17)
-----------------------
* Enhanced inclination correction methodology from simple cos(θ) to account for b/a factor.
* More physically accurate correction accounting for finite disc thickness rather than infinitely thin discs.
* Implemented b/a = sqrt((1-q₀²)*cos²(i) + q₀²) correction where q₀ = 0.2 (intrinsic disc thickness).
* Updated logging to report inclination angle, cos(θ), b/a factor, and adopted q₀ parameter.

"""

# ------------------------------------------------------------------
# User Configuration Parameters
# ------------------------------------------------------------------

# Inclination correction toggle
# Set to True to apply cos(θ) inclination correction, False to disable
apply_inclination_correction = True

# ------------------------------------------------------------------
# 0.  Command‑line interface
# ------------------------------------------------------------------
import argparse
from pathlib import Path

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate extended Voronoi‑binning maps for a MAUVE galaxy"
    )
    p.add_argument(
        "-g", "--galaxy", default="IC3392",
        help="Galaxy identifier, e.g. IC3392 (default: IC3392)"
    )
    p.add_argument(
        "--root", default="/arc/projects/mauve",
        help="Root directory of the MAUVE data tree (default: /arc/projects/mauve)"
    )
    return p.parse_args()

args    = parse_args()
galaxy  = args.galaxy.upper()   # ensure consistent capitalisation
rootdir = Path(args.root).expanduser().resolve()

# ------------------------------------------------------------------
# 1.  File paths derived from CLI args
# ------------------------------------------------------------------
cube_path   = rootdir / "cubes/v3.0"            / f"{galaxy}_DATACUBE_FINAL_WCS_Pall_mad_red_v3.fits"
bin_path    = rootdir / "products/v0.6" / galaxy / f"{galaxy}_SPATIAL_BINNING_maps.fits"
sfh_path    = rootdir / "products/v0.6" / galaxy / f"{galaxy}_SFH_maps.fits"
weight_path = rootdir / "products/v0.6" / galaxy / f"{galaxy}_sfh-weights.fits"
phot_path   = Path("BaSTI+Chabrier.dat")                         # local file
out_path    = Path(f"{galaxy}_SPATIAL_BINNING_maps_extended.fits")   # output in CWD

# For backwards compatibility with variable names in the original notebook
vor_path     = bin_path
binning_path = bin_path

print("\n=== Using the following files ===")
print("Cube         :", cube_path)
print("Binning map  :", bin_path)
print("SFH map      :", sfh_path)
print("Weights      :", weight_path)
print("Photometry   :", phot_path)
print("Output       :", out_path, "\n")

# ------------------------------------------------------------------
# 2.  Imports (same as original)
# ------------------------------------------------------------------
import os, sys, glob, warnings
import numpy as np
import matplotlib.pyplot as plt
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
# 3.  Helper function for inclination correction
# ------------------------------------------------------------------

def read_galaxy_inclination(galaxy_name, inclination_file="MAUVE_Inclination.dat"):
    """
    Read galaxy inclination from MAUVE_Inclination.dat file.
    
    Parameters:
    -----------
    galaxy_name : str
        Name of the galaxy (e.g., 'IC3392')
    inclination_file : str
        Path to the inclination data file
        
    Returns:
    --------
    float
        Inclination angle in degrees, or None if galaxy not found
    """
    try:
        with open(inclination_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) == 2 and parts[0].upper() == galaxy_name.upper():
                    return float(parts[1])
        print(f"Warning: Galaxy {galaxy_name} not found in {inclination_file}")
        return None
    except FileNotFoundError:
        print(f"Warning: Inclination file {inclination_file} not found")
        return None
    except Exception as e:
        print(f"Warning: Error reading inclination file: {e}")
        return None

# ------------------------------------------------------------------
# 4.  Begin main workflow 
# ------------------------------------------------------------------

# Load spatial binning map IC3392_individual.fits 
# --------- file location (edit if needed) ----------
print("Loading:", binning_path.resolve())
with fits.open(binning_path) as hdul:
    # check data structure and header
    print(hdul.info())
    binning_primary = hdul[0]
    binning_BINID   = hdul[1].data
    binning_FLUX    = hdul[2].data
    binning_hdr     = hdul[1].header
    hdul.close()


# Load SFH and weights data IC3392_sfh-weights.fits
# --------- file location (edit if needed) ----------
print("Loading:", weight_path.resolve())
with fits.open(weight_path) as hdul:
    # check data structure and header
    print(hdul.info())
    weights_data = hdul[1].data
    grid_data = hdul[2].data
    weights_hdr  = hdul[1].header
    grid_hdr  = hdul[2].header

    hdul.close()

# ---------- 2.1  Column names (from HEADER_out_phot) ----------
names = [
    'IMF','slope','MH','Age','U','B','V','R','I','J','H','K',
    'UminusV','BminusV','VminusR','VminusI','VminusJ','VminusH','VminusK',
    'ML_U','ML_B','ML_V','ML_R','ML_I','ML_J','ML_H','ML_K',
    'F439W','F555W','F675W','F814W','C439_555','C555_675','C555_814'
]

fname = Path("BaSTI+Chabrier.dat")

# ---------- 2.2  Load data, skip the two header lines ----------
tbl = np.genfromtxt(
    fname, dtype=None, encoding=None, names=names,
    comments='#', skip_header=2, autostrip=True)

# ---------- 2.3  Keep only Chabrier rows ----------
mask = (tbl['IMF'] == 'Ch')
phot = tbl[mask]

print(f"Loaded {len(phot)} out of {len(tbl)} rows of data from {fname.name}")


# --- 3.1  Build a lookup dict keyed by (logAge, MH) rounded to 2 dec —--
key_ml = {}

for row in phot:
    age_gyr = round(row['Age'], 2)                     # e.g. 0.03 → 0.03
    mh_dex  = round(row['MH'],  2)                     # e.g. –2.27 → –2.27
    mlr     = row['ML_R']                   # keep 2 dp as requested
    key_ml[(age_gyr, mh_dex)] = mlr

grid = grid_data  # the FITS_rec you already loaded

# --- 4.1  Prepare new array with an extra ML_R column -------------
mlr_values = np.full(len(grid), np.nan, dtype=np.float32)

for i, (logage, mh, _) in enumerate(grid):
    age_gyr = round(10**logage, 2)   # yrs → Gyr, 2 dp
    mh_dex  = round(mh, 2)                 # already dex
    mlr_values[i] = key_ml.get((age_gyr, mh_dex), np.nan)

# --- 4.2  Build a new structured array including ML_R -------------
ml_dtype = grid.dtype.descr + [('ML_R', 'f4')]
grid_mlr  = np.empty(len(grid), dtype=ml_dtype)

for name in grid.dtype.names:
    grid_mlr[name] = grid[name]
grid_mlr['ML_R'] = mlr_values


# 1) convert the opaque FITS_rec into a plain ndarray
w      = weights_data['WEIGHTS'].astype(np.float32)        # (4077, 477)
ml_ssp = grid_mlr['ML_R'].astype(np.float32)               # (477,)

# 2) light-weighted M/L_R per Voronoi bin (shape 4077)
ml_bin = (w * ml_ssp).sum(axis=1)     # or: np.dot(w, ml_ssp)

# 3) optional sanity check: every bin should return a finite, positive value
assert np.all(np.isfinite(ml_bin)) and (ml_bin > 0).all()


# --- 0)  inputs already in memory ------------------------------------------
# binning_BINID  -> (ny, nx) float array   (NaN   = originally masked pixel
#                                            <0    = masked, but *belongs to* |id|)
# ml_bin         -> (N_bin,) float array   (your 4 077 zone M/L_R values)

# --- 1)  create blank map, same shape & dtype ------------------------------
binning_MLR = np.full_like(binning_BINID, np.nan, dtype=np.float32)

# --- 2)  fill *valid* Voronoi zones ----------------------------------------
valid = binning_BINID >= 0                     # True where BINID is a real zone
binning_MLR[valid] = ml_bin[binning_BINID[valid].astype(int)]


# ---------------------------------------------------------------------
# 1.  R-band magnitude map per spaxel using speclite
# ---------------------------------------------------------------------
print("Loading:", cube_path.resolve())
cube = fits.open(cube_path, memmap=False)
data = cube["DATA"].data.astype(np.float32)          # (nz, ny, nx)
hdr  = cube["DATA"].header
nz, ny, nx = data.shape

# wavelength grid (native header unit → Å)
spec_wcs   = WCS(hdr).sub(["spectral"])              # 1-axis WCS
wave_native = spec_wcs.all_pix2world(
                 np.arange(nz)[:, None], 0)[:, 0]    # numeric values
wave = (wave_native * spec_wcs.wcs.cunit[0]).to(u.AA)   # use cunit[0]

# Load filter using speclite (more robust than synphot for this purpose)
f_r = filters.load_filter('bessell-R')  
# decamDR1noatm-r, bessell-R or 'sdss2010noatm-r'/'sdss2010-r' for SDSS r-band

# Convert flux units: MUSE cube is in 1e-20 erg s⁻¹ cm⁻² Å⁻¹
F_lambda = data * (1e-20 * u.erg / (u.s * u.cm**2 * u.AA))  # Add units 

# Calculate AB magnitude for each spaxel
ny, nx = data.shape[1], data.shape[2]
m_r_map = np.full((ny, nx), np.nan, dtype=np.float32)
maggies = np.full((ny, nx), np.nan, dtype=np.float32)

print("Computing AB magnitudes for each spaxel...")
for j in range(ny):
    if j % 50 == 0:  # Progress indicator
        print(f"Processing row {j}/{ny}")
    for i in range(nx):
        # Extract spectrum for this spaxel
        spectrum = F_lambda[:, j, i]
        m_r_map[j, i] = f_r.get_ab_magnitude(spectrum, wavelength=wave)
        maggies[j, i] = f_r.get_ab_maggies(spectrum, wavelength=wave)
cube.close()

# ---------------------------------------------------------------------
# 2.  Collapse to Voronoi bins
# ---------------------------------------------------------------------
with fits.open(vor_path) as hd_v:
    # ---------- PATCH BEGIN ----------
    # keep BINID as float so NaNs survive the read
    BINID_f32 = hd_v["BINID"].data.astype(np.float32)
    muse_hdr2 = hd_v["BINID"].header
    nan_mask  = ~np.isfinite(hd_v["FLUX"].data)

    # build an *integer* copy with sentinel -1 for “no bin”
    bad_pix = (~np.isfinite(BINID_f32)) | (BINID_f32 < 0)
    BINID   = BINID_f32.astype(np.int32)   # safe: all finite →
    BINID[bad_pix] = -1                   # mark blanks with −1
    # ---------- PATCH END ----------

uniq = np.unique(BINID[BINID >= 0])                 # keep this line

# Convert magnitude map to flux for averaging
F0_ref = 3631e-23  # Reference flux in erg s⁻¹ cm⁻² Hz⁻¹ for AB magnitude zero point
flux_map = F0_ref * 10**(-0.4 * m_r_map)  # Convert mag to flux
flux_map[nan_mask] = np.nan  # Keep NaNs where needed

# Average flux in each bin (excluding NaN pixels)
flux_map_clean = np.where(nan_mask | ~np.isfinite(m_r_map), np.nan, flux_map)
valid_pixels = sum_labels(~np.isnan(flux_map_clean), BINID, uniq)
sum_flux = sum_labels(np.nan_to_num(flux_map_clean), BINID, uniq)
mean_flux = np.divide(sum_flux, valid_pixels, 
                     out=np.full_like(sum_flux, np.nan), 
                     where=valid_pixels > 0)

# Convert back to magnitudes
mean_mag = -2.5 * np.log10(mean_flux / F0_ref)

# Create lookup table for bin-averaged magnitudes
lut = np.full(int(BINID.max()) + 1, np.nan, dtype=np.float32)
lut[uniq] = mean_mag
m_r_binned = lut[BINID]                                    # (ny, nx)
m_r_binned[nan_mask] = np.nan  # Keep NaNs where needed
flux_map_binned = F0_ref * 10**(-0.4 * m_r_binned)  # Convert mag to flux

# ---------------------------------------------------------------------
# 3.  Galactic-extinction correction
# ---------------------------------------------------------------------
EBV = fits.getdata(sfh_path, "EBV").astype(np.float32)

# Choose extinction coefficient based on your filter choice:
if 'bessell' in f_r.name.lower():
    A_r = 2.32 * EBV  # Bessell R-band coefficient (Fitzpatrick 1999)
    M_r_sun = 4.65    # Solar absolute magnitude in Bessell R (AB magnitude), https://mips.as.arizona.edu/~cnaw/sun_2006.html
elif 'sdss' in f_r.name.lower():
    A_r = 2.285 * EBV  # SDSS r-band coefficient (Schlafly & Finkbeiner 2011)
    M_r_sun = 4.65     # Solar absolute magnitude in SDSS r (AB magnitude), https://mips.as.arizona.edu/~cnaw/sun_2006.html
else:
    print(f"Warning: Unknown filter {f_r.name}, using SDSS coefficients")
    A_r = 2.285 * EBV
    M_r_sun = 4.65

# Apply extinction correction to magnitudes (subtract extinction)
m_r_corr = m_r_binned - A_r  # Magnitude correction
m_r_corr[nan_mask] = np.nan

# magnitude back to nanomaggies in Legacy survey format
def magnitude_to_nanomaggies(magnitude):
    return 10**((22.5 - magnitude) / 2.5)
FLUX_R_corr = magnitude_to_nanomaggies(m_r_corr)  # Convert to flux

# ---------------------------------------------------------------------
# 4.  Luminosity → stellar-mass map
# ---------------------------------------------------------------------
# Distance modulus for 16.5 Mpc
distmod = 5 * np.log10((16.5 * u.Mpc).to(u.pc).value / 10)

# Absolute magnitude
M_r = m_r_corr - distmod

# Luminosity in solar units
L_Lsun = 10**(-0.4 * (M_r - M_r_sun))

# Stellar mass (uncomment when binning_MLR is available)
M_star = L_Lsun * binning_MLR
logM_star = np.where(M_star > 0, np.log10(M_star), np.nan)
logM_star[nan_mask] = np.nan

print("R-band magnitude calculation completed!")
print(f"Filter used: {f_r.name}")
print(f"Magnitude range: {np.nanmin(m_r_corr):.2f} to {np.nanmax(m_r_corr):.2f}")
print(f"Distance modulus: {distmod:.2f} mag")

# ──────────────────────────────────────────────────────────────────
#  NEW summary lines
# ──────────────────────────────────────────────────────────────────
tot_mag  = -2.5 * np.log10(np.nansum(flux_map) / F0_ref)
tot_L_R  = np.log10(np.nansum(L_Lsun))          # log₁₀(L/L☉)
tot_M_R  = np.log10(np.nansum(M_star))          # log₁₀(M/M☉)

print(f"Total R-band magnitude   : {tot_mag:.3f} mag (AB)")
print(f"Total R-band luminosity  : {tot_L_R:.3f} log10(L☉)")
print(f"Total stellar mass (R)   : {tot_M_R:.3f} log10(M☉)")


# 1) read the whole HDUList from disk
with fits.open(binning_path) as hdul:
    # 2) clone all existing HDUs into a new list
    new_hdul = fits.HDUList([hdu.copy() for hdu in hdul])

# 3) build the new FLUX_R_corr image HDU
FLUX_R_corr_hdu = fits.ImageHDU(
    data=FLUX_R_corr.astype(np.float64), name="FLUX_R_corr")

# 4) keep WCS and pixel-scale info by copying the original BINID header
FLUX_R_corr_hdu.header.update(binning_hdr)         # you created binning_hdr earlier
FLUX_R_corr_hdu.header["EXTNAME"] = "FLUX_R_corr"
FLUX_R_corr_hdu.header["BUNIT"]   = "nanomaggies"  # physical units

# 5) append and write to disk
new_hdul.append(FLUX_R_corr_hdu)
new_hdul.writeto(out_path, overwrite=True)

print(f"Saved extended file to {out_path.resolve()}")

with fits.open(out_path, mode="append") as hdul:             # open existing file
    new_hdu = fits.ImageHDU(data=binning_MLR.astype(np.float64),  # like others
                             header=binning_hdr, name="ML_R")
    new_hdu.header["EXTNAME"] = "ML_R"                       # name keyword
    new_hdu.header["BUNIT"] = "Msol/Lsol_R"                   # units keyword
    hdul.append(new_hdu)                                    # add as 9-th HDU
    hdul.flush()                                             # write in-place

print("M/L_R layer saved ➜", out_path.resolve())


with fits.open(out_path, mode="append") as hdul:             # open existing file
    mass_hdu = fits.ImageHDU(data=logM_star.astype(np.float64),  # like others
                             header=binning_hdr, name="LOGMSTAR")
    mass_hdu.header["BUNIT"] = "log(Msol)"                   # units keyword
    hdul.append(mass_hdu)                                    # add as 9-th HDU
    hdul.flush()                                             # write in-place

print("Stellar-mass layer saved ➜", out_path.resolve())


# Getting the stellar mass surface density 
# Convert to surface density in M☉/pc²
# 1. Convert pixel area to physical area in pc²
legacy_wcs2 = WCS(binning_hdr).celestial  # strip spectral axis
pixel_scale = (proj_plane_pixel_scales(legacy_wcs2) * u.deg).to(u.arcsec)
pixel_area_Mpc = ((pixel_scale[0]).to(u.rad).value*16.5*u.Mpc)*(((pixel_scale[1]).to(u.rad).value*16.5*u.Mpc))
pixel_area_kpc = pixel_area_Mpc.to(u.kpc**2)

# 2. Read galaxy inclination and calculate correction factor
if apply_inclination_correction:
    galaxy_inclination = read_galaxy_inclination(galaxy)
    if galaxy_inclination is not None:
        inclination_rad = np.deg2rad(galaxy_inclination)
        cos_inclination = np.cos(inclination_rad)
        # Calculate b/a factor: sqrt((1-q0^2)*cos^2(i) + q0^2) where q0=0.2
        ba_factor = np.abs(np.sqrt((1-0.2**2)*cos_inclination**2 + 0.2**2))
        print(f"Galaxy {galaxy} inclination: {galaxy_inclination}° (cos θ = {cos_inclination:.3f})")
        print(f"Inclination correction ENABLED: applying b/a = {ba_factor:.3f} (adopting intrinsic thickness q₀ = 0.2 for disc galaxy)")
    else:
        ba_factor = 1.0
        print(f"No inclination data found for {galaxy}, using ba_factor = 1.0")
else:
    ba_factor = 1.0
    print(f"Inclination correction DISABLED: using ba_factor = 1.0")

# 3. Convert stellar mass to surface density with inclination correction
stellar_mass_surface_density = M_star / pixel_area_kpc  # M☉/kpc²
stellar_mass_surface_density_corrected = stellar_mass_surface_density * ba_factor  # Apply inclination correction
# stellar_mass_surface_density_corrected = stellar_mass_surface_density 

# 4. Convert to log10 scale
log_stellar_mass_surface_density = np.log10(stellar_mass_surface_density_corrected.value)


with fits.open(out_path, mode="append") as hdul:             # open existing file
    mass_density_hdu = fits.ImageHDU(
        data=log_stellar_mass_surface_density.astype(np.float64),  # like others
        header=binning_hdr, name="LOGMASS_SURFACE_DENSITY")
    mass_density_hdu.header["BUNIT"] = "log(Msol/kpc2)"  # units keyword
    hdul.append(mass_density_hdu)                                    # add as 10-th HDU 
    hdul.flush()                                             # write in-place

print("Stellar mass surface density layer saved ➜", out_path.resolve())


with fits.open(out_path, mode="append") as hdul:             # open existing file
    m_r_hdu = fits.ImageHDU(data=m_r_corr.astype(np.float64),  # like others
                             header=binning_hdr, name="magnitude_r")
    m_r_hdu.header["BUNIT"] = "mag_AB"                   # units keyword
    hdul.append(m_r_hdu)                                    # add as 9-th HDU
    hdul.flush()                                             # write in-place

print("r-band magnitude layer saved ➜", out_path.resolve())


with fits.open(out_path, mode="append") as hdul:             # open existing file
    m_r_uncorrected_hdu = fits.ImageHDU(data=m_r_binned.astype(np.float64),  # like others
                             header=binning_hdr, name="magnitude_r_uncorrected")
    m_r_uncorrected_hdu.header["BUNIT"] = "mag_AB"                   # units keyword
    hdul.append(m_r_uncorrected_hdu)                                    # add as 9-th HDU
    hdul.flush()                                             # write in-place

print("r-band uncorrected magnitude layer saved ➜", out_path.resolve())


