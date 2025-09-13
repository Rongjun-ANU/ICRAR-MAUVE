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

Changes (2025-06-30)
-----------------------
* Major refactoring of S/N cut methodology: changed from fixed S/N≥15 to configurable parameters (cut=3, noise=20).
* Complete restructuring of BPT diagram analysis with comprehensive mask system:
  - Added error propagation for BPT ratios and constraint validation
  - Implemented detailed classification: SF, upper-limit, unconstrained, and unknown spaxels
  - Added both [N II] and [S II] BPT diagram analysis with "both" and "either" logic
* Enhanced flux correction methodology for undetected Balmer lines
* Refactored code structure with modular functions and detailed roadmap documentation
* Expanded output with four new SFR surface density maps: SF, upper-limit, unconstrained, and unknown
* Added comprehensive quality control masks for all emission lines

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

# two key cut parameters
cut = 3 # FLUX/FLUX_ERR
noise = 20 # detection limit of FLUX, in the unit of 10^-20 erg/s

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

    V_STARS2 = np.where(non_FOREGROUND_STAR, V_STARS2, np.nan)
    SIGMA_STARS2 = np.where(non_FOREGROUND_STAR, SIGMA_STARS2, np.nan)
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
# 3.  Roadmap
# ------------------------------------------------------------------

# Road map:
# 1. Calculate the Balmer Decrement (BD) from Hβ and Hα
# 2. Convert BD to gas E(B-V) using the Calzetti (2000) extinction curve
# 3. Use E(B-V) to correct the fluxes of the gas lines
# 4. Convert the corrected Hα flux to luminosity
# 5. Calculate the star formation rate (SFR) from the Hα luminosity using the Calzetti (2007) relation
# 6. Calculate the SFR surface density from the SFR map

# Define a function to calculate the Balmer Decrement, 
# with an argument to decide calculate the raw BD, or the corrected BD (i.e., if raw BD < 2.86, then corrected BD = 2.86)
def calculate_balmer_decrement(HB4861_FLUX, HA6562_FLUX, corrected=True):
    BD = HA6562_FLUX / HB4861_FLUX
    # check if an element in BD is NaN or infinite, but it is finite in V_STARS2, then set this element to 2.86
    BD[(~np.isfinite(BD)*np.isfinite(V_STARS2))] = 2.86
    if corrected:
        BD = np.where(BD < 2.86, 2.86, BD)
    return BD

# Calculate the Balmer Decrement
BD = calculate_balmer_decrement(HB4861_FLUX, HA6562_FLUX, corrected=True)
# 

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

# Define a function to convert the BD to gas E(B-V) 
# using the formula E(B-V)_BD = 2.5/(k_HB4861 - k_HA6562) * np.log10(BD/R_int)
def convert_bd_to_ebv(BD, k_HB4861, k_HA6562, R_int=2.86):
    E_BV_BD = 2.5 / (k_HB4861 - k_HA6562) * np.log10(BD / R_int)
    return E_BV_BD

# Calculate the gas E(B-V) from BD
E_BV_BD = convert_bd_to_ebv(BD, k_HB4861, k_HA6562, R_int)

# Use E(B-V)_BD to correct the fluxes
def correct_flux_with_ebv(flux, ebv, k):
    """Correct flux with gas E(B-V) and extinction coefficient k."""
    return flux * 10**(0.4 * k * ebv)

# Correct the fluxes with E(B-V)_BD
HB4861_FLUX_corr = correct_flux_with_ebv(HB4861_FLUX, E_BV_BD, k_HB4861)
HA6562_FLUX_corr = correct_flux_with_ebv(HA6562_FLUX, E_BV_BD, k_HA6562)
OIII5006_FLUX_corr = correct_flux_with_ebv(OIII5006_FLUX, E_BV_BD, k_OIII5006)
NII6583_FLUX_corr = correct_flux_with_ebv(NII6583_FLUX, E_BV_BD, k_NII6583)
SII6716_FLUX_corr = correct_flux_with_ebv(SII6716_FLUX, E_BV_BD, k_SII6716)
SII6730_FLUX_corr = correct_flux_with_ebv(SII6730_FLUX, E_BV_BD, k_SII6730)

###################
# Modify the the corrected Flux map to deal with the case that Halpha and/or Hbeta are not detected. 

# Balmer detection masks: (HB4861_FLUX/HB4861_FLUX_ERR>=cut) & (HA6562_FLUX/HA6562_FLUX_ERR>=cut)
Balmer_detected = ((((HB4861_FLUX / HB4861_FLUX_ERR) >= cut) & (HB4861_FLUX >= noise)) & ((HA6562_FLUX / HA6562_FLUX_ERR) >= cut) & (HA6562_FLUX >= noise))
Balmer_not_detected = ~Balmer_detected

# If there is a spaxel that Halpha and/or Hbeta are not detected (Balmer_not_detected), all lines' fluxes in that spaxel are set to max(noise, FLUX_Corr) in the unit of 10^-20 erg/s
def modify_Balmer_not_detected_map(flux_map, flux_raw_map, mask=Balmer_not_detected, noise=noise): 
    """
    Apply a mask to the flux map based on Balmer detection.

    Parameters:
    flux_map : array-like
        The flux map to be modified.
    mask : array-like, optional
        The mask indicating where to apply the correction (default is Balmer_not_detected).
    noise : float, optional
        The noise level to set for undetected regions (default is 20).
        
    Returns:
    modified_flux_map : array-like
        The modified flux map with undetected regions set to max(noise, FLUX_Corr).
    """
    modified_flux_map = flux_map.copy()
    # For spaxels where Balmer lines are not detected, set flux to max(noise, original corrected flux)
    modified_flux_map[mask] = np.maximum(noise, flux_raw_map[mask])
    
    return modified_flux_map

# Apply the Balmer detection mask to the corrected flux maps for Further calculation of SFR
HA6562_FLUX_Corr = modify_Balmer_not_detected_map(flux_map=HA6562_FLUX_corr, flux_raw_map=HA6562_FLUX, mask=Balmer_not_detected, noise=noise)

###################

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
    return (flux*1e-20*u.erg/u.s/u.cm**2 * 4*np.pi*(distance*u.Mpc)**2).cgs.value

# Calculate the luminosity of Halpha
HA6562_LUM = flux_to_luminosity(HA6562_FLUX_Corr)

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
    return 5.3e-42 * luminosity / 0.67 * 0.63  # SFR in M_sun/yr

# Calculate the SFR map from Halpha luminosity
SFR_map = calzetti_sfr(HA6562_LUM)

# Getting the SFR surface density 
# Convert to surface density in M☉/pc²
# 1. Convert pixel area to physical area in pc²
legacy_wcs2 = WCS(gas_header).celestial  # strip spectral axis
pixel_scale = (proj_plane_pixel_scales(legacy_wcs2) * u.deg).to(u.arcsec)
pixel_area_Mpc = ((pixel_scale[0]).to(u.rad).value*16.5*u.Mpc)*(((pixel_scale[1]).to(u.rad).value*16.5*u.Mpc))
pixel_area_kpc = pixel_area_Mpc.to(u.kpc**2)
# 2. Convert SFR to surface density
SFR_surface_density_map = SFR_map / pixel_area_kpc.value

# 3. Convert to log10 scale
LOG_SFR_surface_density_map = np.log10(SFR_surface_density_map)

# ------------------------------------------------------------------
# 4.  Masks: basic QC cut
# ------------------------------------------------------------------

# Define a function to apply signal to noise cut at lines, then return the masks
def apply_QC(cut=cut, noise=noise): 
    QC_good = {
        'HB4861': ((HB4861_FLUX / HB4861_FLUX_ERR) >= cut) & (HB4861_FLUX >= noise),
        'HA6562': ((HA6562_FLUX / HA6562_FLUX_ERR) >= cut) & (HA6562_FLUX >= noise),
        'OIII5006': ((OIII5006_FLUX / OIII5006_FLUX_ERR) >= cut) & (OIII5006_FLUX >= noise),
        'NII6583': ((NII6583_FLUX / NII6583_FLUX_ERR) >= cut) & (NII6583_FLUX >= noise),
        'SII6716': ((SII6716_FLUX / SII6716_FLUX_ERR) >= cut) & (SII6716_FLUX >= noise),
        'SII6730': ((SII6730_FLUX / SII6730_FLUX_ERR) >= cut) & (SII6730_FLUX >= noise)
    }
    QC_bad = {
        'HB4861': ((HB4861_FLUX / HB4861_FLUX_ERR) < cut) | (HB4861_FLUX < noise),
        'HA6562': ((HA6562_FLUX / HA6562_FLUX_ERR) < cut) | (HA6562_FLUX < noise),
        'OIII5006': ((OIII5006_FLUX / OIII5006_FLUX_ERR) < cut) | (OIII5006_FLUX < noise),
        'NII6583': ((NII6583_FLUX / NII6583_FLUX_ERR) < cut) | (NII6583_FLUX < noise),
        'SII6716': ((SII6716_FLUX / SII6716_FLUX_ERR) < cut) | (SII6716_FLUX < noise),
        'SII6730': ((SII6730_FLUX / SII6730_FLUX_ERR) < cut) | (SII6730_FLUX < noise)
    }
    return QC_good, QC_bad

# Apply the SNR cut to each line
QC_good, QC_bad = apply_QC(cut=cut, noise=noise)

# Extract individual masks
HB4861_QC_good = QC_good['HB4861']
HB4861_QC_bad = QC_bad['HB4861']
HA6562_QC_good = QC_good['HA6562']
HA6562_QC_bad = QC_bad['HA6562']
OIII5006_QC_good = QC_good['OIII5006']
OIII5006_QC_bad = QC_bad['OIII5006']
NII6583_QC_good = QC_good['NII6583']
NII6583_QC_bad = QC_bad['NII6583']
SII6716_QC_good = QC_good['SII6716']
SII6716_QC_bad = QC_bad['SII6716']
SII6730_QC_good = QC_good['SII6730']
SII6730_QC_bad = QC_bad['SII6730']

# ------------------------------------------------------------------
# 5.  Masks: BPT selection: HII, Comp, AGN in [NII] BPT; HII, LINER, Seyfert in [SII] BPT.  
# ------------------------------------------------------------------

# ---- line ratios --------------------------------------------------
logN2  = np.log10(NII6583_FLUX_corr / HA6562_FLUX_corr)        # [N II]/Hα
logS2  = np.log10((SII6716_FLUX_corr+SII6730_FLUX_corr) / HA6562_FLUX_corr)   # Σ[S II]/Hα
logO3  = np.log10(OIII5006_FLUX_corr / HB4861_FLUX_corr)       # [O III]/Hβ

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

# Create x arrays for the theoretical lines
x_kewley_N2 = np.linspace(-2.0, 0.3, 200)
x_kauff_N2 = np.linspace((286-np.sqrt(2871561))/1100, 0.0, 200)
x_kewley_S2 = np.linspace(-2.0, 0.3, 200)
x_kewley06_S2 = np.linspace((159-np.sqrt(105081))/525, 0.5, 200)

# Define a function to apply the BPT masks, 
# the BPT masks are to find the HII, Comp, and AGN regions in NII BPT, 
# while the HII, LINER, and Seyfert regions in SII BPT, respectively.
def apply_bpt_masks(logN2, logS2, logO3):
    # NII BPT masks
    mask_N2_HII = (logO3 < kauff03_N2(logN2)) & (logO3 < kewley01_N2(logN2)) & (logN2 < 0.05)
    mask_N2_Comp = (logO3 >= kauff03_N2(logN2)) & (logO3 < kewley01_N2(logN2)) & (logN2 < 0.47)
    mask_N2_AGN = (logO3 >= kewley01_N2(logN2)) | (logN2 >= 0.47)

    # SII BPT masks
    mask_S2_HII = (logO3 < kewley01_S2(logS2)) & (logS2 < 0.32)
    mask_S2_LINER = (((logO3 >= kewley01_S2(logS2)) & (logS2 < 0.32)) | (logS2 >= 0.32)) & (logO3 < kewley06_Sy_LIN(logS2))
    mask_S2_Seyfert = (((logO3 >= kewley01_S2(logS2)) & (logS2 < 0.32)) | (logS2 >= 0.32)) & (logO3 >= kewley06_Sy_LIN(logS2))

    return (mask_N2_HII, mask_N2_Comp, mask_N2_AGN), (mask_S2_HII, mask_S2_LINER, mask_S2_Seyfert)

# Apply the BPT masks
masks_N2, masks_S2 = apply_bpt_masks(logN2, logS2, logO3)
mask_N2_HII, mask_N2_Comp, mask_N2_AGN = masks_N2
mask_S2_HII, mask_S2_LINER, mask_S2_Seyfert = masks_S2

# NII SF and non-SF masks
mask_N2_SF = mask_N2_HII | mask_N2_Comp
mask_N2_nonSF = mask_N2_AGN
# SII SF and non-SF masks
mask_S2_SF = mask_S2_HII
mask_S2_nonSF = mask_S2_LINER | mask_S2_Seyfert

# ------------------------------------------------------------------
# 6.  Masks: for [NII] BPT
# ------------------------------------------------------------------

# Flag 4 final cases that we want to track

# definite SF spaxels: or HA_detected_HB_detected & mask_N2_constrainted & mask_N2_SF
mask_SF_N2 = ((HA_detected_HB_detected_NII_detected_OIII_detected & mask_N2_constrainted & mask_N2_SF) | 
                    (HA_detected_HB_detected_NII_not_detected_OIII_not_detected & mask_N2_constrainted & mask_N2_SF) | 
                    (HA_detected_HB_detected_NII_not_detected_OIII_detected & mask_N2_constrainted & mask_N2_SF) | 
                    (HA_detected_HB_detected_NII_detected_OIII_not_detected & mask_N2_constrainted & mask_N2_SF))
# get SFR as nonSF: : or HA_detected_HB_detected & mask_N2_constrainted & mask_N2_nonSF
mask_nonSF_N2 = ((HA_detected_HB_detected_NII_detected_OIII_detected & mask_N2_constrainted & mask_N2_nonSF) |
              (HA_detected_HB_detected_NII_not_detected_OIII_not_detected & mask_N2_constrainted & mask_N2_nonSF) |
              (HA_detected_HB_detected_NII_not_detected_OIII_detected & mask_N2_constrainted & mask_N2_nonSF) |
              (HA_detected_HB_detected_NII_detected_OIII_not_detected & mask_N2_constrainted & mask_N2_nonSF))
# all the unclassified spaxels: : or HA_detected_HB_detected & mask_N2_unconstrainted
mask_unclassified_N2 = ((HA_detected_HB_detected_NII_not_detected_OIII_not_detected & mask_N2_unconstrainted) | 
                       (HA_detected_HB_detected_NII_not_detected_OIII_detected & mask_N2_unconstrainted) | 
                       (HA_detected_HB_detected_NII_detected_OIII_not_detected & mask_N2_unconstrainted) | 
                       (HA_detected_HB_detected_NII_detected_OIII_detected & mask_N2_unconstrainted))
# the rest are upper spaxels: 
mask_upper = (HA_not_detected | HA_detected_HB_not_detected)

# Something else might be useful

# all the constrainted spaxels
mask_constrainted_N2 = mask_SF_N2 | mask_nonSF_N2

# ------------------------------------------------------------------
# 7.  Masks: for [SII] BPT
# ------------------------------------------------------------------

# Flag 4 final cases that we want to track for [SII] BPT

# definite SF spaxels: or HA_detected_HB_detected & mask_S2_constrainted & mask_S2_SF
mask_SF_S2 = ((HA_detected_HB_detected_SII_detected_OIII_detected & mask_S2_constrainted & mask_S2_SF) | 
              (HA_detected_HB_detected_SII_not_detected_OIII_not_detected & mask_S2_constrainted & mask_S2_SF) | 
              (HA_detected_HB_detected_SII_not_detected_OIII_detected & mask_S2_constrainted & mask_S2_SF) | 
              (HA_detected_HB_detected_SII_detected_OIII_not_detected & mask_S2_constrainted & mask_S2_SF))
# get SFR as nonSF: or HA_detected_HB_detected & mask_S2_constrainted & mask_S2_nonSF
mask_nonSF_S2 = ((HA_detected_HB_detected_SII_detected_OIII_detected & mask_S2_constrainted & mask_S2_nonSF) |
                  (HA_detected_HB_detected_SII_not_detected_OIII_not_detected & mask_S2_constrainted & mask_S2_nonSF) |
                  (HA_detected_HB_detected_SII_not_detected_OIII_detected & mask_S2_constrainted & mask_S2_nonSF) |
                  (HA_detected_HB_detected_SII_detected_OIII_not_detected & mask_S2_constrainted & mask_S2_nonSF))
# all the unclassified spaxels: or HA_detected_HB_detected & mask_S2_unconstrainted
mask_unclassified_S2 = ((HA_detected_HB_detected_SII_not_detected_OIII_not_detected & mask_S2_unconstrainted) | 
                           (HA_detected_HB_detected_SII_not_detected_OIII_detected & mask_S2_unconstrainted) | 
                           (HA_detected_HB_detected_SII_detected_OIII_not_detected & mask_S2_unconstrainted) | 
                           (HA_detected_HB_detected_SII_detected_OIII_detected & mask_S2_unconstrainted))
# # the rest are upper spaxels: 
# mask_upper = (HA_not_detected | HA_detected_HB_not_detected)

# Something else might be useful

# all the constrained spaxels
mask_constrainted_S2 = mask_SF_S2 | mask_nonSF_S2

# ------------------------------------------------------------------
# 8.  Masks: Combine two BPT
# ------------------------------------------------------------------

# both:

# SF: SF in both N2 and S2 BPT diagrams:
mask_SF_both = mask_SF_N2 & mask_SF_S2
# nonSF: constrained in both N2 and S2 BPT diagrams, but not SF in either or both:
mask_nonSF_both = ((mask_constrainted_N2 & mask_constrainted_S2) & ~mask_SF_both)
# unclassified: unconstrained in either N2 or S2 BPT diagrams: or 
# mask_unclassified_both = (mask_unclassified_N2 | mask_unclassified_S2)
mask_unclassified_both = ((~(mask_constrainted_N2 & mask_constrainted_S2)) & HA_detected_HB_detected)
# # upper: the rest are upper spaxels:
# mask_upper = (HA_not_detected | HA_detected_HB_not_detected)

# Something else might be useful
# all the constrained spaxels: or mask_constrainted_both = ~mask_unclassified_both
mask_constrainted_both = (mask_constrainted_N2 & mask_constrainted_S2)

# either:

# SF: SF in either N2 or S2 BPT diagrams:
mask_SF_either = mask_SF_N2 | mask_SF_S2
# nonSF: constrained in either N2 or S2 BPT diagrams, but not SF in either or both:
mask_nonSF_either = ((mask_constrainted_N2 | mask_constrainted_S2) & ~mask_SF_either)
# unclassified: unconstrained in either N2 or S2 BPT diagrams: or 
# mask_unclassified_either = (mask_unclassified_N2 & mask_unclassified_S2)
mask_unclassified_either = ((~(mask_constrainted_N2 | mask_constrainted_S2)) & HA_detected_HB_detected)
# # upper: the rest are upper spaxels:
# mask_upper = (HA_not_detected | HA_detected_HB_not_detected)

# Something else might be useful
# all the constrained spaxels: or mask_constrainted_either = ~mask_unclassified_either
mask_constrainted_either = (mask_constrainted_N2 | mask_constrainted_S2)

# ------------------------------------------------------------------
# 9.  Append Σ_SFR layers (choose 'both' for now, but can be changed to 'either' or just fall back to 'N2' or 'S2')
# ------------------------------------------------------------------

def choose_BPT(choice='both'):
    """
    Choose BPT classification method and return SFR surface density maps and masks.
    
    Parameters:
    choice : str, optional
        BPT classification choice. Options: 'both', 'either', 'N2', 'S2' (default: 'both')
        
    Returns:
    tuple : (SFR_maps, masks) where:
        - SFR_maps: tuple of four arrays (SF, nonSF, unclassified, upper) for the chosen BPT method
        - masks: tuple of four boolean arrays (mask_SF, mask_nonSF, mask_unclassified, mask_upper)
    """
    # Get the appropriate masks based on choice
    if choice == 'both':
        mask_SF = mask_SF_both
        mask_nonSF = mask_nonSF_both
        mask_unclassified = mask_unclassified_both
    elif choice == 'either':
        mask_SF = mask_SF_either
        mask_nonSF = mask_nonSF_either
        mask_unclassified = mask_unclassified_either
    elif choice == 'N2':
        mask_SF = mask_SF_N2
        mask_nonSF = mask_nonSF_N2
        mask_unclassified = mask_unclassified_N2
    elif choice == 'S2':
        mask_SF = mask_SF_S2
        mask_nonSF = mask_nonSF_S2
        mask_unclassified = mask_unclassified_S2
    else:
        raise ValueError(f"Invalid choice '{choice}'. Options: 'both', 'either', 'N2', 'S2'")
    
    # Apply masks to create SFR surface density maps
    LOG_SFR_surface_density_map_SF = np.where(mask_SF, LOG_SFR_surface_density_map, np.nan)
    LOG_SFR_surface_density_map_nonSF = np.where(mask_nonSF, LOG_SFR_surface_density_map, np.nan)
    LOG_SFR_surface_density_map_unclassified = np.where(mask_unclassified, LOG_SFR_surface_density_map, np.nan)
    LOG_SFR_surface_density_map_upper = np.where(mask_upper, LOG_SFR_surface_density_map, np.nan)
    
    # Return both SFR maps and masks
    sfr_maps = (LOG_SFR_surface_density_map_SF, LOG_SFR_surface_density_map_nonSF, 
                LOG_SFR_surface_density_map_unclassified, LOG_SFR_surface_density_map_upper)
    masks = (mask_SF, mask_nonSF, mask_unclassified, mask_upper)
    
    return sfr_maps, masks

# Get the SFR surface density maps and masks using the default 'both' choice
(LOG_SFR_surface_density_map_SF, LOG_SFR_surface_density_map_nonSF, 
 LOG_SFR_surface_density_map_unclassified, LOG_SFR_surface_density_map_upper), (mask_SF, mask_nonSF, mask_unclassified, mask_upper) = choose_BPT()

with fits.open(gas_path) as hdul:
    new_hdul = fits.HDUList([hdu.copy() for hdu in hdul])

# Corrected Hα luminosity
hdu_halpha_lum = fits.ImageHDU(HA6562_LUM.astype(np.float64),
                               header=gas_header, name="Halpha_Luminosity_corr")
hdu_halpha_lum.header['BUNIT'] = 'erg/s'
new_hdul.append(hdu_halpha_lum)
# Corrected SFR
hdu_sfr = fits.ImageHDU(SFR_map.astype(np.float64),
                        header=gas_header, name="Halpha_SFR_corr")
hdu_sfr.header['BUNIT'] = 'M_sun/yr'
new_hdul.append(hdu_sfr)
# log Σ_SFR
hdu_logsfr = fits.ImageHDU(LOG_SFR_surface_density_map.astype(np.float64),
                           header=gas_header, name="LOGSFR_SURFACE_DENSITY")
hdu_logsfr.header['BUNIT'] = 'log(M_sun/yr/kpc2)'
new_hdul.append(hdu_logsfr)
hdu_logsfr_sf = fits.ImageHDU(LOG_SFR_surface_density_map_SF.astype(np.float64),
                              header=gas_header, name="LOGSFR_SURFACE_DENSITY_SF")
hdu_logsfr_sf.header['BUNIT'] = 'log(M_sun/yr/kpc2)'
new_hdul.append(hdu_logsfr_sf)
hdu_logsfr_nonSF = fits.ImageHDU(LOG_SFR_surface_density_map_nonSF.astype(np.float64),
                                 header=gas_header, name="LOGSFR_SURFACE_DENSITY_NONSF")
hdu_logsfr_nonSF.header['BUNIT'] = 'log(M_sun/yr/kpc2)'
new_hdul.append(hdu_logsfr_nonSF)
hdu_logsfr_unclassified = fits.ImageHDU(LOG_SFR_surface_density_map_unclassified.astype(np.float64),
                                           header=gas_header, name="LOGSFR_SURFACE_DENSITY_UNCLASSIFIED")
hdu_logsfr_unclassified.header['BUNIT'] = 'log(M_sun/yr/kpc2)'
new_hdul.append(hdu_logsfr_unclassified)
hdu_logsfr_upper = fits.ImageHDU(LOG_SFR_surface_density_map_upper.astype(np.float64),
                                   header=gas_header, name="LOGSFR_SURFACE_DENSITY_UPPER")
hdu_logsfr_upper.header['BUNIT'] = 'log(M_sun/yr/kpc2)'
new_hdul.append(hdu_logsfr_upper)

new_hdul.writeto(out_path, overwrite=True)
print("Extended file written ➜", out_path.resolve())

# ------------------------------------------------------------------
# 10.  Print some useful information
# ------------------------------------------------------------------

# Print the total non-nan spaxels
print("--------------------------------------------------------------")
total_spaxels = np.sum(np.isfinite(V_STARS2))
print("Total non-nan spaxels:", total_spaxels)
# Print the number of 6 cases that need number, 2 upper cases, and 4 unclassified cases
print("Number of pixels with Halpha not detected:", np.sum(HA_not_detected))
print("Number of pixels with Halpha detected, Hbeta not detected:", np.sum(HA_detected_HB_not_detected))
print("Number of pixels with Halpha detected, Hbeta detected, NII not detected, OIII not detected and unclassified:", 
      np.sum(HA_detected_HB_detected_NII_not_detected_OIII_not_detected & mask_N2_unconstrainted))
print("Number of pixels with Halpha detected, Hbeta detected, NII not detected, OIII detected and unclassified:",
      np.sum(HA_detected_HB_detected_NII_not_detected_OIII_detected & mask_N2_unconstrainted))
print("Number of pixels with Halpha detected, Hbeta detected, NII detected, OIII not detected and unclassified:",
      np.sum(HA_detected_HB_detected_NII_detected_OIII_not_detected & mask_N2_unconstrainted))
print("Number of pixels with Halpha detected, Hbeta detected, NII detected, OIII detected and unclassified:",
      np.sum(HA_detected_HB_detected_NII_detected_OIII_detected & mask_N2_unconstrainted))
print("Number of pixels with Halpha detected, Hbeta detected, SII not detected, OIII not detected and unclassified:", 
      np.sum(HA_detected_HB_detected_SII_not_detected_OIII_not_detected & mask_N2_unconstrainted))
print("Number of pixels with Halpha detected, Hbeta detected, SII not detected, OIII detected and unclassified:",
      np.sum(HA_detected_HB_detected_SII_not_detected_OIII_detected & mask_N2_unconstrainted))
print("Number of pixels with Halpha detected, Hbeta detected, SII detected, OIII not detected and unclassified:",
      np.sum(HA_detected_HB_detected_SII_detected_OIII_not_detected & mask_N2_unconstrainted))
print("Number of pixels with Halpha detected, Hbeta detected, SII detected, OIII detected and unclassified:",
      np.sum(HA_detected_HB_detected_SII_detected_OIII_detected & mask_N2_unconstrainted))
print("--------------------------------------------------------------")
print(f"Total corrected Halpha luminosity: {np.nansum(HA6562_LUM):.2e} erg/s")
print(f"Total corrected Halpha luminosity from SF region: {np.nansum(SFR_map[mask_SF]):.2e} erg/s")
print(f"Total Halpha SFR: {np.nansum(SFR_map):.2f} M☉/yr or in log10 scale: {np.log10(np.nansum(SFR_map)):.2f} log(M☉/yr)")
print(f"Total Halpha SFR from SF region: {np.nansum(SFR_map[mask_SF]):.2f} M☉/yr or in log10 scale: {np.log10(np.nansum(SFR_map[mask_SF])):.2f} log(M☉/yr)")
print("--------------------------------------------------------------")

print(f"Total runtime: {time.perf_counter() - t0:.1f} s")
