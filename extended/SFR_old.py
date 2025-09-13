#!/usr/bin/env python
"""
SFR.py – produce Balmer-decrement, dust-corrected Hα maps and SFR
surface-density maps for a MAUVE galaxy.

Changes (2025-06-17)
-----------------------
* Foreground-star mask from {gal}_KIN_maps_extended.fits (if present).
* Calzetti (2000) extinction curve coded explicitly.
* Correct definition of “upper-limit’’ E(B–V) and dereddened fluxes:
      – BD_upper = 2.86 for spaxels that fail the S/N≥15 cut but still
        have a finite stellar velocity (V_STARS2).
* Two Σ_SFR layers written:
      LOGSFR_SURFACE_DENSITY          – pure SF spaxels
      LOGSFR_SURFACE_DENSITY_UPPER    – all spaxels with S/N≥15
"""
import argparse, logging, time, sys
from pathlib import Path
import numpy as np

from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales


# ------------------------------------------------------------------
# 0.  Command-line interface  (exactly as requested)
# ------------------------------------------------------------------
import argparse, logging, sys, time
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

p = argparse.ArgumentParser(description="Generate SFR maps for a MAUVE galaxy")
p.add_argument("-g", "--galaxy", default="IC3392", help="Galaxy ID (default IC3392)")
p.add_argument("--root", default="/arc/projects/mauve", help="MAUVE root path")
p.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
args = p.parse_args()

loglvl = logging.INFO if args.verbose else logging.WARNING
logging.basicConfig(level=loglvl, format="%(levelname)s %(message)s", stream=sys.stdout)

t0   = time.perf_counter()
gal  = args.galaxy.upper()
root = Path(args.root)

gas_path = root / "products" / "v0.6" / gal / f"{gal}_gas_BIN_maps.fits"
out_path = Path(f"{gal}_gas_BIN_maps_extended.fits")
kin_path = Path(f"{gal}_KIN_maps_extended.fits")          # foreground + stellar V

# ------------------------------------------------------------------
# 1.  Load gas-line maps (extended version file if it already exists)
# ------------------------------------------------------------------
if gas_path.exists():
    src = gas_path
elif Path(out_path).exists():
    src = Path(out_path)
else:
    sys.exit(f"Cannot find {gas_path} or {out_path}")

print(f"Reading gas-line FITS ➜ {src}")
with fits.open(src) as hdul:
    V_STARS2 = hdul['V_STARS2'].data
    SIGMA_STARS2 = hdul['SIGMA_STARS2'].data
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
    gas_header = hdul[5].header
    hdul.close()

gas_header

# ------------------------------------------------------------------
# 2.  Foreground-star removal  (verbatim from notebook snippet)
# ------------------------------------------------------------------
if kin_path.exists():
    print(f"Loading kinematic map from {kin_path}")
    with fits.open(kin_path) as hdul:
        kin_info = hdul.info()
        
        # Read data from all extensions except PRIMARY
        extension_names = [hdul[i].name for i in range(1, len(hdul))]
        print(f"Available extensions: {extension_names}")
        
        # Read each extension's data before closing the file
        for ext_name in extension_names:
            if ext_name and ext_name != "PRIMARY":
                globals()[ext_name] = hdul[ext_name].data
                print(f"Loaded {ext_name}: shape {globals()[ext_name].shape}")
    
    print("All data loaded successfully!")
    hdul.close()
    # Invert the true/false in FOREGROUND_STAR, but except the nan values
    non_FOREGROUND_STAR = np.where(np.isnan(FOREGROUND_STAR), np.nan, ~FOREGROUND_STAR.astype(bool))

    V_STARS2 = np.where(non_FOREGROUND_STAR, V, np.nan)
    SIGMA_STARS2 = np.where(non_FOREGROUND_STAR, SIGMA, np.nan)
    HB4861_FLUX = np.where(non_FOREGROUND_STAR, HB4861_FLUX, np.nan)
    HB4861_FLUX_ERR = np.where(non_FOREGROUND_STAR, HB4861_FLUX_ERR, np.nan)
    HA6562_FLUX = np.where(non_FOREGROUND_STAR, HA6562_FLUX, np.nan)
    HA6562_FLUX_ERR = np.where(non_FOREGROUND_STAR, HA6562_FLUX_ERR, np.nan)
    OIII5006_FLUX = np.where(non_FOREGROUND_STAR, OIII5006_FLUX, np.nan)
    OIII5006_FLUX_ERR = np.where(non_FOREGROUND_STAR, OIII5006_FLUX_ERR, np.nan)
    NII6583_FLUX = np.where(non_FOREGROUND_STAR, NII6583_FLUX, np.nan)
    NII6583_FLUX_ERR = np.where(non_FOREGROUND_STAR, NII6583_FLUX_ERR, np.nan)
    SII6716_FLUX = np.where(non_FOREGROUND_STAR, SII6716_FLUX, np.nan)
    SII6716_FLUX_ERR = np.where(non_FOREGROUND_STAR, SII6716_FLUX_ERR, np.nan)
    SII6730_FLUX = np.where(non_FOREGROUND_STAR, SII6730_FLUX, np.nan)
    SII6730_FLUX_ERR = np.where(non_FOREGROUND_STAR, SII6730_FLUX_ERR, np.nan)
    print("Foreground stars are removed successfully!")

else:
    print(f"File not found: {kin_path}")

# ------------------------------------------------------------------
# 3.  S/N≥15 cut on both Balmer lines
# ------------------------------------------------------------------
cut = 15
mask_HB = HB4861_FLUX / HB4861_FLUX_ERR >= cut
mask_HA = HA6562_FLUX / HA6562_FLUX_ERR >= cut
# Combine masks for both lines
mask_combined = mask_HB & mask_HA
# Apply the mask to the flux maps
HB4861_FLUX_cut = np.where(mask_combined, HB4861_FLUX, np.nan)
HA6562_FLUX_cut = np.where(mask_combined, HA6562_FLUX, np.nan)
HB4861_FLUX_out = np.where(~mask_combined, HB4861_FLUX, np.nan)
HA6562_FLUX_out = np.where(~mask_combined, HA6562_FLUX, np.nan)
OIII5006_FLUX_cut = np.where(mask_combined, OIII5006_FLUX, np.nan)
OIII5006_FLUX_out = np.where(~mask_combined, OIII5006_FLUX, np.nan)
NII6583_FLUX_cut = np.where(mask_combined, NII6583_FLUX, np.nan)
NII6583_FLUX_out = np.where(~mask_combined, NII6583_FLUX, np.nan)
SII6716_FLUX_cut = np.where(mask_combined, SII6716_FLUX, np.nan)
SII6716_FLUX_out = np.where(~mask_combined, SII6716_FLUX, np.nan)
SII6730_FLUX_cut = np.where(mask_combined, SII6730_FLUX, np.nan)
SII6730_FLUX_out = np.where(~mask_combined, SII6730_FLUX, np.nan)

# Balmer-decrement map – H α/H β (start with all spaxels).
BD =  HA6562_FLUX / HB4861_FLUX
BD_cut = HA6562_FLUX_cut / HB4861_FLUX_cut
BD_out = HA6562_FLUX_out / HB4861_FLUX_out
BD_upper = np.where(~mask_combined & np.isfinite(V_STARS2), 2.86, BD)

# Print the lowest and highest values of the Balmer decrement map
print(f"Lowest Balmer Decrement: {np.nanmin(BD_cut)}")
print(f"Highest Balmer Decrement: {np.nanmax(BD_cut)}")
# Print lowest 5 non-NaN values (unique values only)
bd_valid = BD_cut[np.isfinite(BD_cut)]
bd_unique = np.unique(bd_valid)
lowest_5_unique = bd_unique[:5]
print(f"Lowest 5 unique non-NaN Balmer Decrement values: {lowest_5_unique}")

# ------------------------------------------------------------------
# 4.  Calzetti (2000) k(λ) and E(B–V)
# ------------------------------------------------------------------
# Calzetti (2000) curve
def calzetti_curve(wavelengths):
    """Calzetti (2000) extinction curve."""
    # Convert single values to array
    if np.isscalar(wavelengths):
        wavelengths = np.array([wavelengths])
        
    # Calzetti (2000) parameters
    Rv = 4.05  # R_V for Calzetti law
    A_lambda = np.zeros_like(wavelengths, dtype=float)
    
    # Calculate the extinction for each wavelength
    mask_short = wavelengths < 6300
    
    # Short wavelengths
    A_lambda[mask_short] = 2.659 * (-2.156 + 1.509/wavelengths[mask_short] 
                                   - 0.198/wavelengths[mask_short]**2 
                                   + 0.011/wavelengths[mask_short]**3) + Rv
    # Long wavelengths
    A_lambda[~mask_short] = 2.659 * (-1.857 + 1.040/wavelengths[~mask_short]) + Rv
    
    return A_lambda[0] if np.isscalar(wavelengths) else A_lambda

# Calculate k values for Hβ and Hα
k_HB4861 = calzetti_curve(0.4861)  # Hβ
k_HA6562 = calzetti_curve(0.6562)  # Hα
k_OIII5006 = calzetti_curve(0.5006)  # [OIII] 5006
k_NII6583 = calzetti_curve(0.6583)  # [NII] 6583
k_SII6716 = calzetti_curve(0.6716)  # [SII] 6716
k_SII6730 = calzetti_curve(0.6730)  # [SII] 6730

R_int = 2.86

E_BV_BD_cut = 2.5/(k_HB4861 - k_HA6562) * np.log10( (HA6562_FLUX_cut/HB4861_FLUX_cut) / R_int )
E_BV_BD_upper = 2.5/(k_HB4861 - k_HA6562) * np.log10( BD_upper / R_int )

HB4861_FLUX_upper = HB4861_FLUX * 10**(0.4 * k_HB4861 * E_BV_BD_upper)
HA6562_FLUX_upper = HA6562_FLUX * 10**(0.4 * k_HA6562 * E_BV_BD_upper)
OIII5006_FLUX_upper = OIII5006_FLUX * 10**(0.4 * k_OIII5006 * E_BV_BD_upper)
NII6583_FLUX_upper = NII6583_FLUX * 10**(0.4 * k_NII6583 * E_BV_BD_upper)
SII6716_FLUX_upper = SII6716_FLUX * 10**(0.4 * k_SII6716 * E_BV_BD_upper)
SII6730_FLUX_upper = SII6730_FLUX * 10**(0.4 * k_SII6730 * E_BV_BD_upper)

# ------------------------------------------------------------------
# 5.  BPT masks
# ------------------------------------------------------------------
logN2  = np.log10(NII6583_FLUX_upper / HA6562_FLUX_upper)        # [N II]/Hα
logS2  = np.log10((SII6716_FLUX_upper+SII6730_FLUX_upper) / HA6562_FLUX_upper)   # Σ[S II]/Hα
logO3  = np.log10(OIII5006_FLUX_upper / HB4861_FLUX_upper)       # [O III]/Hβ         

mask_N2 = NII6583_FLUX_upper / NII6583_FLUX_ERR >= cut
mask_S2 = (SII6716_FLUX_upper + SII6730_FLUX_upper) / (SII6716_FLUX_ERR + SII6730_FLUX_ERR) >= cut
mask_O3 = OIII5006_FLUX_upper / OIII5006_FLUX_ERR >= cut
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

# Apply the mask to the flux maps
HB4861_FLUX_SF = np.where(mask_SF*mask_combinedd, HB4861_FLUX_upper, np.nan)
HA6562_FLUX_SF = np.where(mask_SF*mask_combinedd, HA6562_FLUX_upper, np.nan)
OIII5006_FLUX_SF = np.where(mask_SF*mask_combinedd, OIII5006_FLUX_upper, np.nan)
NII6583_FLUX_SF = np.where(mask_SF*mask_combinedd, NII6583_FLUX_upper, np.nan)
SII6716_FLUX_SF = np.where(mask_SF*mask_combinedd, SII6716_FLUX_upper, np.nan)
SII6730_FLUX_SF = np.where(mask_SF*mask_combinedd, SII6730_FLUX_upper, np.nan)

HB4861_FLUX_out_SF = np.where(mask_SF*(~mask_combinedd), HB4861_FLUX_upper, np.nan)
HA6562_FLUX_out_SF = np.where(mask_SF*(~mask_combinedd), HA6562_FLUX_upper, np.nan)
OIII5006_FLUX_out_SF = np.where(mask_SF*(~mask_combinedd), OIII5006_FLUX_upper, np.nan)
NII6583_FLUX_out_SF = np.where(mask_SF*(~mask_combinedd), NII6583_FLUX_upper, np.nan)
SII6716_FLUX_out_SF = np.where(mask_SF*(~mask_combinedd), SII6716_FLUX_upper, np.nan)
SII6730_FLUX_out_SF = np.where(mask_SF*(~mask_combinedd), SII6730_FLUX_upper, np.nan)

# ------------------------------------------------------------------
# 6.  Luminosity, SFR, Σ_SFR
# ------------------------------------------------------------------
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

# Calculate luminosity for all the flux maps
lum_HB4861_upper = flux_to_luminosity(HB4861_FLUX_upper)
lum_HA6562_upper = flux_to_luminosity(HA6562_FLUX_upper)
lum_OIII5006_upper = flux_to_luminosity(OIII5006_FLUX_upper)
lum_NII6583_upper = flux_to_luminosity(NII6583_FLUX_upper)
lum_SII6716_upper = flux_to_luminosity(SII6716_FLUX_upper)
lum_SII6730_upper = flux_to_luminosity(SII6730_FLUX_upper)

lum_HB4861_SF = flux_to_luminosity(HB4861_FLUX_SF)
lum_HA6562_SF = flux_to_luminosity(HA6562_FLUX_SF)
lum_OIII5006_SF = flux_to_luminosity(OIII5006_FLUX_SF)
lum_NII6583_SF = flux_to_luminosity(NII6583_FLUX_SF)
lum_SII6716_SF = flux_to_luminosity(SII6716_FLUX_SF)
lum_SII6730_SF = flux_to_luminosity(SII6730_FLUX_SF)

lum_HB4861_out_SF = flux_to_luminosity(HB4861_FLUX_out_SF)
lum_HA6562_out_SF = flux_to_luminosity(HA6562_FLUX_out_SF)
lum_OIII5006_out_SF = flux_to_luminosity(OIII5006_FLUX_out_SF)
lum_NII6583_out_SF = flux_to_luminosity(NII6583_FLUX_out_SF)    
lum_SII6716_out_SF = flux_to_luminosity(SII6716_FLUX_out_SF)
lum_SII6730_out_SF = flux_to_luminosity(SII6730_FLUX_out_SF)

total_lum_HA6562_upper = np.nansum(lum_HA6562_upper)
total_lum_HA6562_SF = np.nansum(lum_HA6562_SF)

print(f"Total luminosity of Hα (corrected): {total_lum_HA6562_upper:.3e} erg/s")
print(f"Total luminosity of Hα (SF spaxels): {total_lum_HA6562_SF:.3e} erg/s")

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
    return 5.3e-42 * luminosity.cgs.value / 0.67 *0.63 # SFR in M_sun/yr

# Calculate SFR for all the Halpha maps
sfr_HA6562_upper = calzetti_sfr(lum_HA6562_upper)
sfr_HA6562_SF = calzetti_sfr(lum_HA6562_SF)
sfr_HA6562_out_SF = calzetti_sfr(lum_HA6562_out_SF)

# Getting the SFR surface density 
# Convert to surface density in M☉/pc²
# 1. Convert pixel area to physical area in pc²
legacy_wcs2 = WCS(gas_header).celestial  # strip spectral axis
pixel_scale = (proj_plane_pixel_scales(legacy_wcs2) * u.deg).to(u.arcsec)
pixel_area_Mpc = ((pixel_scale[0]).to(u.rad).value*16.5*u.Mpc)*(((pixel_scale[1]).to(u.rad).value*16.5*u.Mpc))
pixel_area_kpc = pixel_area_Mpc.to(u.kpc**2)
# 2. Convert SFR to surface density
sfr_surface_density_HA6562_upper = sfr_HA6562_upper / pixel_area_kpc.value
sfr_surface_density_HA6562_SF = sfr_HA6562_SF / pixel_area_kpc.value
sfr_surface_density_HA6562_out_SF = sfr_HA6562_out_SF / pixel_area_kpc.value
# 3. Convert to log10 scale
log_sfr_surface_density_HA6562_upper = np.log10(sfr_surface_density_HA6562_upper)
log_sfr_surface_density_HA6562_SF = np.log10(sfr_surface_density_HA6562_SF)
log_sfr_surface_density_HA6562_out_SF = np.log10(sfr_surface_density_HA6562_out_SF)

# ------------------------------------------------------------------
# 7.  Append Σ_SFR layers
# ------------------------------------------------------------------
with fits.open(gas_path) as hdul:
    new_hdul = fits.HDUList([hdu.copy() for hdu in hdul])

# Balmer Decrement
new_hdul.append(fits.ImageHDU(BD_upper.astype(np.float64),
                              header=gas_header, name="Balmer_Decrement_upper"))
# E(B-V)_BD
new_hdul.append(fits.ImageHDU(E_BV_BD_upper.astype(np.float64),
                              header=gas_header, name="E(B-V)_BD_upper"))
# Corrected Hα luminosity
hdu_halpha_lum = fits.ImageHDU(lum_HA6562_SF.value.astype(np.float64),
                               header=gas_header, name="Halpha_Luminosity_SF_corr")
hdu_halpha_lum.header['BUNIT'] = 'erg/s'
new_hdul.append(hdu_halpha_lum)
# Corrected SFR
hdu_sfr = fits.ImageHDU(sfr_HA6562_SF.astype(np.float64),
                        header=gas_header, name="Halpha_SFR_SF_corr")
hdu_sfr.header['BUNIT'] = 'M_sun/yr'
new_hdul.append(hdu_sfr)
# log Σ_SFR
hdu_logsfr_sf = fits.ImageHDU(log_sfr_surface_density_HA6562_SF.astype(np.float64),
                              header=gas_header, name="LOGSFR_SURFACE_DENSITY_SF")
hdu_logsfr_sf.header['BUNIT'] = 'log(M_sun/yr/kpc2)'
new_hdul.append(hdu_logsfr_sf)
hdu_logsfr_upper = fits.ImageHDU(log_sfr_surface_density_HA6562_upper.astype(np.float64),
                                 header=gas_header, name="LOGSFR_SURFACE_DENSITY_UPPER")
hdu_logsfr_upper.header['BUNIT'] = 'log(M_sun/yr/kpc2)'
new_hdul.append(hdu_logsfr_upper)
hdu_logsfr_out_sf = fits.ImageHDU(log_sfr_surface_density_HA6562_out_SF.astype(np.float64),
                                 header=gas_header, name="LOGSFR_SURFACE_DENSITY_OUT_SF")
hdu_logsfr_out_sf.header['BUNIT'] = 'log(M_sun/yr/kpc2)'
new_hdul.append(hdu_logsfr_out_sf)

new_hdul.writeto(out_path, overwrite=True)
print("Extended file written ➜", out_path.resolve())

print(f"Total runtime: {time.perf_counter() - t0:.1f} s")

