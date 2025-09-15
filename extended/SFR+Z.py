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
  - Added error propagation for BPT ratios and classified validation
  - Implemented detailed classification: SF, upper-limit, unconstrained, and upper spaxels
  - Added both [N II] and [S II] BPT diagram analysis with "both" and "either" logic
* Enhanced flux correction methodology for undetected Balmer lines
* Refactored code structure with modular functions and detailed roadmap documentation
* Expanded output with four new SFR surface density maps: SF, upper-limit, unconstrained, and upper
* Added comprehensive quality control masks for all emission lines

Changes (2025-07-28)
-----------------------
* Metallicity [O/H] calculation (12+log(O/H)) added using different methods: Dopita et al. (2016), Pilyugin & Grebel (2016). 

Changes (2025-09-04)
-----------------------
* Fix the Error Propagation in BPT Analysis
* Default to CANFAR path for gas maps; keep local-testing path commented out.
* Replaced `calzetti_curve` with vectorized helper `calzetti_k(w_um)` for Calzetti (2000) k(λ).
* Extended `choose_BPT()` to also return SF-masked, dust-corrected line-flux maps
  (HB4861, HA6562, OIII5006, NII6583, SII6716, SII6730).
* Added SF-integrated metallicities computed from SF-summed lines:
  - `O_H_D16_SF_total` (Dopita+2016),
  - `O_H_PG16_SF_total` (Pilyugin & Grebel 2016; branch via total log N2).
* Expanded terminal summary to include SF-region totals.
* No FITS schema changes; integrated information are kept in log file (not written to the FITS).

Changes (2025-09-08)
-----------------------
* Added O3N2-M13 (Marino et al. 2013) metallicity calibration: [O/H] = 8.533 - 0.214 * O3N2
* Added N2-M13 (Marino et al. 2013) metallicity calibration: [O/H] = 8.743 + 0.462 * N2
* Added O3N2-PP04 (Pettini & Pagel 2004) metallicity calibration: [O/H] = 8.73 - 0.32 * O3N2
* Added N2-PP04 (Pettini & Pagel 2004) metallicity calibration: [O/H] = 9.37 + 2.03*N2 + 1.26*N2^2 + 0.32*N2^3
* Added comprehensive Curti et al. (2020) C20 metallicity calibration suite:
  - O3N2-C20: Quadratic equation solver for O3N2 index
  - O3S2-C20: Quartic polynomial root-finding for O3S2 index  
  - RS32-C20: Quartic polynomial for RS32 = log([OIII]/Hβ + [SII]/Hα)
  - R3-C20: Cubic polynomial for R3 = log([OIII]/Hβ)
  - N2-C20: Quartic polynomial with strict range selection for N2 index
  - S2-C20: Quartic polynomial for S2 = log([SII]/Hα)
* Added Combined-C20 metallicity using priority-based method selection per spaxel
  (Note: Generally N/A to MAUVE data due to limited line coverage, but calculated for completeness)
* Optimized C20 calibrations to be independent without cross-dependencies on D16 metallicity
* Integrated all metallicity maps and total calculations for SF regions with comprehensive FITS output
* Updated terminal summary to include all metallicity method totals with method usage statistics
* Enhanced polynomial root-finding with sophisticated tolerance systems and range enforcement

Changes (2025-09-11)
-----------------------
* Added comprehensive metallicity calculations for total available regions (Section 11):
  - Extended D16, PG16, M13, PP04, and all C20 calibrations to total flux calculations
  - Implemented proper flux summing and extinction correction for integrated measurements
  - Added all C20 method calculations (O3N2, O3S2, RS32, R3, N2, S2, Combined) for total regions
* Enhanced terminal output with expanded metallicity reporting:
  - Added detailed metallicity totals for both SF regions and total available regions
  - Improved summary statistics with comprehensive method comparisons
* Restructured code organization with clearer section numbering (12 → 13 for final output)
* Added robust error handling and validation for all total metallicity calculations
* Maintained backwards compatibility with existing FITS output structure

Changes (2025-09-14)
-----------------------
* Added inclination correction for SFR surface density calculation.
* Implemented read_galaxy_inclination() function to read inclination angles from MAUVE_Inclination.dat.
* Applied cos(θ) correction factor to SFR_surface_density_map where θ is the galaxy inclination angle.
* Enhanced logging to show inclination values and correction factors applied.

Changes (2025-09-15)
-----------------------
* Added user-configurable inclination correction parameter.
* Users can now enable/disable inclination correction by setting apply_inclination_correction = True/False.

"""

# ------------------------------------------------------------------
# User Configuration Parameters
# ------------------------------------------------------------------

# Inclination correction toggle
# Set to True to apply cos(θ) inclination correction, False to disable
apply_inclination_correction = False

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

# ------------------------------------------------------------------
# Helper function for inclination correction
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

gas_path = root / "products" / "v0.6" / gal / f"{gal}_gas_BIN_maps.fits" # For CANFAR
gas_path = Path(f"{gal}_gas_BIN_maps.fits") # For local testing
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
# 3. Use E(B-V) to correct the fluxes of the gas lines, then use different methods to calculate the metallicity [O/H] (12+log(O/H))
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

# Calzetti (2000) curve
def calzetti_k(w_um):
    """Return k(λ) = A(λ)/E(B−V) for Calzetti (2000); wavelengths in microns."""
    import numpy as np
    w = np.asarray(w_um, dtype=float)
    Rv = 4.05
    k = np.empty_like(w, dtype=float)

    short = (w >= 0.12) & (w < 0.63)
    long  = (w >= 0.63) & (w <= 2.2)

    k[short] = 2.659 * (-2.156 + 1.509/w[short] - 0.198/w[short]**2 + 0.011/w[short]**3) + Rv
    k[long]  = 2.659 * (-1.857 + 1.040/w[long]) + Rv
    return k.item() if k.ndim == 1 and k.size == 1 else k

k_HB4861  = calzetti_k(0.4861)  # ≈ 4.598
k_HA6562  = calzetti_k(0.6562)  # ≈ 3.326
k_OIII5006= calzetti_k(0.5006)  # ≈ 4.465
k_NII6583 = calzetti_k(0.6583)  # ≈ 3.313
k_SII6716 = calzetti_k(0.6716)  # ≈ 3.230
k_SII6730 = calzetti_k(0.6730)  # ≈ 3.221

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

# ------------------------------------------------------------------
# Metallicity [O/H] calculation (12+log(O/H)) using different methods
# ------------------------------------------------------------------

# Error propogation for BPT diagrams (sigma of log_10(numerator/denominator))
def bpt_error_propagation(numerator, denominator, numerator_err, denominator_err):
    """
    Calculate the propagated error for the BPT ratio log10(numerator/denominator).
    
    Parameters:
    numerator (np.ndarray): The numerator values.
    denominator (np.ndarray): The denominator values.
    numerator_err (np.ndarray): The error in the numerator.
    denominator_err (np.ndarray): The error in the denominator.
    
    Returns:
    np.ndarray: The propagated error for the BPT ratio.
    """
    # Avoid division by zero
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = numerator / denominator
        log_ratio = np.log10(ratio)
        log_ratio_err = 1/(np.log(10)) * np.sqrt((numerator_err / numerator)**2 + (denominator_err / denominator)**2)
        return log_ratio_err

# Dopita et al. (2016) metallicity calculation
y = np.log10(NII6583_FLUX_corr / (SII6716_FLUX_corr + SII6730_FLUX_corr)) + 0.264*np.log10(NII6583_FLUX_corr / HA6562_FLUX_corr)
O_H_D16 = 8.77 + y + 0.45*(y + 0.3)**5
# Set O_H_D16 to be nan if outside the range of 7.63 and 9.23
O_H_D16 = np.where((O_H_D16 < 7.63) | (O_H_D16 > 9.23), np.nan, O_H_D16)

# Pilyugin & Grebel (2016) metallicity calculation (the S calibration)
# note that here we assume [O III] = 1.33 [O III] 5007, [N II] = 1.34 [N II] 6583, see watts et al. (2024) for details
# PG16 set different coefficients for different branches (logN_2>=-0.6 and logN_2<-0.6)
OIII_scaled = 1.33 * OIII5006_FLUX_corr  # [O III] = 1.33 * [O III] 5006
NII_scaled = 1.34 * NII6583_FLUX_corr     # [N II] = 1.34 * [N II] 6583
# Calculate the line ratios needed for PG16
N2 = NII_scaled / HB4861_FLUX_corr   # N2 = I([N II]λ6548 + λ6584)/I(Hβ)
S2 = (SII6716_FLUX_corr + SII6730_FLUX_corr) / HB4861_FLUX_corr  # S2 = I([S II]λ6717 + λ6731)/I(Hβ)
R3 = OIII_scaled / HB4861_FLUX_corr  # R3 = I([O III]λ4959 + λ5007)/I(Hβ) (same value as R2 in this case)
# Calculate log values
log_R3_S2 = np.log10(R3/S2)
log_N2 = np.log10(N2)
log_S2 = np.log10(S2)
# Determine which branch to use based on log(N2)
# Upper branch: log(N2) >= -0.6
# Lower branch: log(N2) < -0.6
# Initialize arrays for the results - preserve original shape and fill with NaN
O_H_PG16 = np.full_like(log_N2, np.nan)
# Upper branch coefficients (log N2 >= -0.6)
upper_mask = log_N2 >= -0.6
a1_upper = 8.424
a2_upper = 0.030
a3_upper = 0.751
a4_upper = -0.349
a5_upper = 0.182
a6_upper = 0.508
# Lower branch coefficients (log N2 < -0.6)
lower_mask = log_N2 < -0.6
a1_lower = 8.072
a2_lower = 0.789
a3_lower = 0.726
a4_lower = 1.069
a5_lower = -0.170
a6_lower = 0.022
# Calculate (O/H)S,U for upper branch
O_H_PG16[upper_mask] = (a1_upper + a2_upper * log_R3_S2[upper_mask] + a3_upper * log_N2[upper_mask] + 
                        (a4_upper + a5_upper * log_R3_S2[upper_mask] + a6_upper * log_N2[upper_mask]) * log_S2[upper_mask])
# Calculate (O/H)S,L for lower branch  
O_H_PG16[lower_mask] = (a1_lower + a2_lower * log_R3_S2[lower_mask] + a3_lower * log_N2[lower_mask] + 
                        (a4_lower + a5_lower * log_R3_S2[lower_mask] + a6_lower * log_N2[lower_mask]) * log_S2[lower_mask])
# Set O_H_PG16 to be nan if outside the range of 7.63 and 9.23
O_H_PG16 = np.where((O_H_PG16 < 7.63) | (O_H_PG16 > 9.23), np.nan, O_H_PG16)

# O3N2-M13 (Marino et al. 2013) metallicity calculation function
def calculate_o3n2_m13_metallicity(hb4861_flux, oiii5006_flux, nii6583_flux, ha6562_flux, oh_d16_sf):
    """Calculate [O/H] using O3N2-M13 (Marino et al. 2013) calibration: [O/H] = 8.533 - 0.214 * O3N2"""
    # Use basic finite checks on emission lines
    good_mask = (np.isfinite(hb4861_flux) & np.isfinite(oiii5006_flux) &
                np.isfinite(nii6583_flux) & np.isfinite(ha6562_flux) &
                (hb4861_flux > 0) & (oiii5006_flux > 0) &
                (nii6583_flux > 0) & (ha6562_flux > 0))
    
    # Calculate O3N2 ratio and then [O/H] metallicity using M13 calibration
    oh_o3n2_m13 = np.full_like(hb4861_flux, np.nan)
    oiii_hb = oiii5006_flux / hb4861_flux
    nii_ha = nii6583_flux / ha6562_flux
    o3n2_ratio = np.log10(oiii_hb / nii_ha)
    # Apply O3N2-M13 (Marino et al. 2013) calibration: [O/H] = 8.533 - 0.214 * O3N2
    oh_o3n2_m13[good_mask] = 8.533 - 0.214 * o3n2_ratio[good_mask]
    
    return oh_o3n2_m13, good_mask

# Calculate O3N2-M13 metallicity
O_H_O3N2_M13, o3n2_m13_good_mask = calculate_o3n2_m13_metallicity(HB4861_FLUX_corr, OIII5006_FLUX_corr, 
                                                                  NII6583_FLUX_corr, HA6562_FLUX_corr, O_H_D16)
# Set O_H_O3N2_M13 to be nan if outside the range of 7.63 and 9.23
O_H_O3N2_M13 = np.where((O_H_O3N2_M13 < 7.63) | (O_H_O3N2_M13 > 9.23), np.nan, O_H_O3N2_M13)

# N2-M13 (Marino et al. 2013) metallicity calculation function
def calculate_n2_m13_metallicity(nii6583_flux, ha6562_flux, oh_d16_sf):
    """Calculate [O/H] using N2-M13 (Marino et al. 2013) calibration: [O/H] = 8.743 + 0.462*N2"""
    # Use basic finite checks on emission lines
    good_mask = (np.isfinite(nii6583_flux) & np.isfinite(ha6562_flux) &
                (nii6583_flux > 0) & (ha6562_flux > 0))
    
    # Calculate N2 ratio and then [O/H] metallicity using M13 calibration
    oh_n2_m13 = np.full_like(nii6583_flux, np.nan)
    n2_ratio = np.log10(nii6583_flux / ha6562_flux)
    # Apply N2-M13 (Marino et al. 2013) calibration: [O/H] = 8.743 + 0.462*N2
    oh_n2_m13[good_mask] = 8.743 + 0.462 * n2_ratio[good_mask]
    
    return oh_n2_m13, good_mask

# Calculate N2-M13 metallicity
O_H_N2_M13, n2_m13_good_mask = calculate_n2_m13_metallicity(NII6583_FLUX_corr, HA6562_FLUX_corr, O_H_D16)
# Set O_H_N2_M13 to be nan if outside the range of 7.63 and 9.23
O_H_N2_M13 = np.where((O_H_N2_M13 < 7.63) | (O_H_N2_M13 > 9.23), np.nan, O_H_N2_M13)

# O3N2-PP04 (Pettini & Pagel 2004) metallicity calculation function
def calculate_o3n2_pp04_metallicity(hb4861_flux, oiii5006_flux, nii6583_flux, ha6562_flux, oh_d16_sf):
    """Calculate [O/H] using O3N2-PP04 (Pettini & Pagel 2004) calibration: [O/H] = 8.73 - 0.32 * O3N2"""
    # Use basic finite checks on emission lines
    good_mask = (np.isfinite(hb4861_flux) & np.isfinite(oiii5006_flux) &
                np.isfinite(nii6583_flux) & np.isfinite(ha6562_flux) &
                (hb4861_flux > 0) & (oiii5006_flux > 0) &
                (nii6583_flux > 0) & (ha6562_flux > 0))
    
    # Calculate O3N2 ratio and then [O/H] metallicity using PP04 calibration
    oh_o3n2_pp04 = np.full_like(hb4861_flux, np.nan)
    oiii_hb = oiii5006_flux / hb4861_flux
    nii_ha = nii6583_flux / ha6562_flux
    o3n2_ratio = np.log10(oiii_hb / nii_ha)
    # Apply O3N2-PP04 (Pettini & Pagel 2004) calibration: [O/H] = 8.73 - 0.32 * O3N2
    oh_o3n2_pp04[good_mask] = 8.73 - 0.32 * o3n2_ratio[good_mask]
    
    return oh_o3n2_pp04, good_mask

# Calculate O3N2-PP04 metallicity
O_H_O3N2_PP04, o3n2_pp04_good_mask = calculate_o3n2_pp04_metallicity(HB4861_FLUX_corr, OIII5006_FLUX_corr, 
                                                                     NII6583_FLUX_corr, HA6562_FLUX_corr, O_H_D16)
# Set O_H_O3N2_PP04 to be nan if outside the range of 7.63 and 9.23
O_H_O3N2_PP04 = np.where((O_H_O3N2_PP04 < 7.63) | (O_H_O3N2_PP04 > 9.23), np.nan, O_H_O3N2_PP04)

# N2-PP04 (Pettini & Pagel 2004) metallicity calculation function
def calculate_n2_pp04_metallicity(nii6583_flux, ha6562_flux, oh_d16_sf):
    """Calculate [O/H] using N2-PP04 (Pettini & Pagel 2004) calibration: [O/H] = 9.37 + 2.03*N2 + 1.26*N2^2 + 0.32*N2^3"""
    # Use basic finite checks on emission lines
    good_mask = (np.isfinite(nii6583_flux) & np.isfinite(ha6562_flux) &
                (nii6583_flux > 0) & (ha6562_flux > 0))
    
    # Calculate N2 ratio and then [O/H] metallicity using PP04 calibration
    oh_n2_pp04 = np.full_like(nii6583_flux, np.nan)
    n2_ratio = np.log10(nii6583_flux / ha6562_flux)
    # Apply N2-PP04 (Pettini & Pagel 2004) calibration: [O/H] = 9.37 + 2.03*N2 + 1.26*N2^2 + 0.32*N2^3
    oh_n2_pp04[good_mask] = (9.37 + 2.03 * n2_ratio[good_mask] + 
                            1.26 * n2_ratio[good_mask]**2 + 
                            0.32 * n2_ratio[good_mask]**3)
    
    return oh_n2_pp04, good_mask

# Calculate N2-PP04 metallicity
O_H_N2_PP04, n2_pp04_good_mask = calculate_n2_pp04_metallicity(NII6583_FLUX_corr, HA6562_FLUX_corr, O_H_D16)
# Set O_H_N2_PP04 to be nan if outside the range of 7.63 and 9.23
O_H_N2_PP04 = np.where((O_H_N2_PP04 < 7.63) | (O_H_N2_PP04 > 9.23), np.nan, O_H_N2_PP04)

# O3N2-C20 (Curti et al. 2020) metallicity calculation function
def calculate_o3n2_c20_metallicity(hb4861_flux, oiii5006_flux, nii6583_flux, ha6562_flux, 
                                   hb4861_flux_err, oiii5006_flux_err, nii6583_flux_err, ha6562_flux_err, 
                                   oh_d16_sf):
    """Calculate [O/H] using O3N2-C20 calibration from Curti+2020 with error propagation"""
    # Use basic finite checks on emission lines
    good_mask = (np.isfinite(hb4861_flux) & np.isfinite(oiii5006_flux) &
                np.isfinite(nii6583_flux) & np.isfinite(ha6562_flux) &
                (hb4861_flux > 0) & (oiii5006_flux > 0) &
                (nii6583_flux > 0) & (ha6562_flux > 0))
    
    # Additional checks for positive fluxes and finite errors where O/H is valid
    good_mask = (good_mask & (hb4861_flux > 0) & (oiii5006_flux > 0) & (nii6583_flux > 0) & (ha6562_flux > 0) &
                 np.isfinite(hb4861_flux_err) & np.isfinite(oiii5006_flux_err) & 
                 np.isfinite(nii6583_flux_err) & np.isfinite(ha6562_flux_err))
    
    # Initialize output arrays
    oh_o3n2_c20 = np.full_like(hb4861_flux, np.nan)
    oh_o3n2_c20_err = np.full_like(hb4861_flux, np.nan)
    
    # Calculate O3N2 ratio and errors
    oiii_hb = oiii5006_flux / hb4861_flux
    nii_ha = nii6583_flux / ha6562_flux
    o3n2_ratio = np.log10(oiii_hb / nii_ha)
    
    # Calculate errors for the line ratios using error propagation
    oiii_hb_err = bpt_error_propagation(oiii5006_flux, hb4861_flux, oiii5006_flux_err, hb4861_flux_err)
    nii_ha_err = bpt_error_propagation(nii6583_flux, ha6562_flux, nii6583_flux_err, ha6562_flux_err)
    
    # Error for O3N2 = log10(OIII/Hb / NII/Ha) = log10(OIII/Hb) - log10(NII/Ha)
    # Error propagation: sqrt(err1^2 + err2^2) for difference of independent variables
    o3n2_ratio_err = np.sqrt(oiii_hb_err**2 + nii_ha_err**2)
    
    # Apply O3N2-C20 calibration (Curti+2020)
    # Step 1: Compute R = O3N2 and y = log10(R)
    R = o3n2_ratio  # This is already log10(O3N2)
    y = R
    y_err = o3n2_ratio_err
    
    # Step 2: Solve quadratic equation y - (c0 + c1*x + c2*x^2) = 0 for x
    # Coefficients from Curti+2020
    c0 = 0.281
    c1 = -4.765
    c2 = -2.268
    
    # Rearrange to standard form: c2*x^2 + c1*x + (c0 - y) = 0
    # Using quadratic formula: x = (-b ± sqrt(b^2 - 4ac)) / (2a)
    a = c2
    b = c1
    c = c0 - y
    
    # Calculate discriminant
    discriminant = b**2 - 4*a*c
    
    # Only calculate where discriminant is positive
    valid_discriminant = discriminant >= 0
    combined_mask = good_mask & valid_discriminant
    
    if np.any(combined_mask):
        # Use the appropriate root (typically the positive one for metallicity)
        x_solution1 = (-b + np.sqrt(discriminant[combined_mask])) / (2*a)
        x_solution2 = (-b - np.sqrt(discriminant[combined_mask])) / (2*a)
        
        # Choose the solution that gives reasonable metallicity values
        # Typically x should be in the range corresponding to 12+log(O/H) ~ 7.6 to 8.85
        x_final = np.where((x_solution1 >= -1.1) & (x_solution1 <= 1.25), x_solution1, x_solution2)
        
        # Calculate error in x using derivative approach
        # For equation f(x,y) = y - (c0 + c1*x + c2*x^2) = 0
        # df/dx = -(c1 + 2*c2*x), df/dy = 1
        # x_err = |df/dy| * y_err / |df/dx| = y_err / (|c1| + |2*c2*x|)
        derivative_x = np.abs(c1 + 2*c2*x_final)
        x_err = y_err[combined_mask] / derivative_x
        
        # Step 3: Return 12 + log(O/H) = x + 8.69
        oh_o3n2_c20[combined_mask] = x_final + 8.69
        # Add intrinsic fitting error from O3S2-C20 calibration (0.09 dex)
        fitting_err = 0.09  # dex
        oh_o3n2_c20_err[combined_mask] = np.sqrt(x_err**2)

    return oh_o3n2_c20, oh_o3n2_c20_err, combined_mask

# Calculate O3N2-C20 metallicity with error propagation
O_H_O3N2_C20, O_H_O3N2_C20_ERR, o3n2_c20_good_mask = calculate_o3n2_c20_metallicity(
    HB4861_FLUX_corr, OIII5006_FLUX_corr, NII6583_FLUX_corr, HA6562_FLUX_corr,
    HB4861_FLUX_ERR, OIII5006_FLUX_ERR, NII6583_FLUX_ERR, HA6562_FLUX_ERR, 
    O_H_D16)
# Set O_H_O3N2_C20 to be nan if outside the range of 7.63 and 9.23
O_H_O3N2_C20 = np.where((O_H_O3N2_C20 < 7.63) | (O_H_O3N2_C20 > 9.23), np.nan, O_H_O3N2_C20)
O_H_O3N2_C20_ERR = np.where((O_H_O3N2_C20 < 7.63) | (O_H_O3N2_C20 > 9.23), np.nan, O_H_O3N2_C20_ERR)

# O3S2-C20 (Curti et al. 2020) metallicity calculation function
def calculate_o3s2_c20_metallicity(hb4861_flux, oiii5006_flux, sii6716_flux, sii6730_flux, 
                                   hb4861_flux_err, oiii5006_flux_err, sii6716_flux_err, sii6730_flux_err, 
                                   oh_d16_sf):
    """Calculate [O/H] using O3S2-C20 calibration from Curti+2020 with error propagation"""
    # Use basic finite checks on emission lines
    good_mask = (np.isfinite(hb4861_flux) & np.isfinite(oiii5006_flux) &
                np.isfinite(sii6716_flux) & np.isfinite(sii6730_flux) &
                (hb4861_flux > 0) & (oiii5006_flux > 0) &
                (sii6716_flux > 0) & (sii6730_flux > 0))
    
    # Additional checks for positive fluxes and finite errors where O/H is valid
    good_mask = (good_mask & (hb4861_flux > 0) & (oiii5006_flux > 0) & (sii6716_flux > 0) & (sii6730_flux > 0) &
                 np.isfinite(hb4861_flux_err) & np.isfinite(oiii5006_flux_err) & 
                 np.isfinite(sii6716_flux_err) & np.isfinite(sii6730_flux_err))
    
    # Initialize output arrays
    oh_o3s2_c20 = np.full_like(hb4861_flux, np.nan)
    oh_o3s2_c20_err = np.full_like(hb4861_flux, np.nan)
    
    # Calculate line ratios and errors
    oiii_hb = oiii5006_flux / hb4861_flux
    sii_total = sii6716_flux + sii6730_flux  # Total [SII] flux
    sii_total_err = np.sqrt(sii6716_flux_err**2 + sii6730_flux_err**2)  # Error for sum
    sii_hb = sii_total / hb4861_flux
    
    # Calculate O3S2 ratio: ([OIII]/Hβ) / ([SII]/Hβ) = [OIII]/[SII]
    o3s2_ratio = np.log10(oiii_hb / sii_hb)
    
    # Calculate errors for the line ratios using error propagation
    oiii_hb_err = bpt_error_propagation(oiii5006_flux, hb4861_flux, oiii5006_flux_err, hb4861_flux_err)
    sii_hb_err = bpt_error_propagation(sii_total, hb4861_flux, sii_total_err, hb4861_flux_err)
    
    # Error for O3S2 = log10(OIII/Hb / SII/Hb) = log10(OIII/Hb) - log10(SII/Hb)
    # Error propagation: sqrt(err1^2 + err2^2) for difference of independent variables
    o3s2_ratio_err = np.sqrt(oiii_hb_err**2 + sii_hb_err**2)
    
    # Apply O3S2-C20 calibration (Curti+2020)
    # Step 1: Compute R = O3S2 and y = log10(R)
    R = o3s2_ratio  # This is already log10(O3S2)
    y = R
    y_err = o3s2_ratio_err
    
    # Step 2: Solve polynomial equation y - (c0 + c1*x + c2*x^2 + c3*x^3 + c4*x^4) = 0 for x
    # Coefficients from Curti+2020 for O3S2
    c0 = 0.191
    c1 = -4.292
    c2 = -2.538
    c3 = 0.053
    c4 = 0.332
    
    # This is now a 4th order polynomial: c4*x^4 + c3*x^3 + c2*x^2 + c1*x + (c0 - y) = 0
    # We need to solve this numerically for each valid spaxel
    combined_mask = np.copy(good_mask)
    
    if np.any(good_mask):
        valid_indices = np.where(good_mask)
        for i in range(len(valid_indices[0])):
            idx_y, idx_x = valid_indices[0][i], valid_indices[1][i]
            y_val = y[idx_y, idx_x]
            y_err_val = y_err[idx_y, idx_x]
            
            # Polynomial coefficients for numpy.roots (highest degree first)
            poly_coeffs = [c4, c3, c2, c1, (c0 - y_val)]
            roots = np.roots(poly_coeffs)
            
            # Select the real root that gives reasonable metallicity values
            real_roots = roots[np.isreal(roots)].real
            if len(real_roots) > 0:
                # Choose root that gives reasonable x values (around -1 to 1 for typical metallicities)
                reasonable_roots = real_roots[(real_roots >= -2) & (real_roots <= 2)]
                if len(reasonable_roots) > 0:
                    x_final = reasonable_roots[0]  # Take first reasonable root
                    
                    # Calculate error in x using derivative approach
                    # For equation f(x,y) = y - (c0 + c1*x + c2*x^2 + c3*x^3 + c4*x^4) = 0
                    # df/dx = -(c1 + 2*c2*x + 3*c3*x^2 + 4*c4*x^3), df/dy = 1
                    # x_err = |df/dy| * y_err / |df/dx| = y_err / (|c1| + |2*c2*x| + |3*c3*x^2| + |4*c4*x^3|)
                    derivative_x = (np.abs(c1 + 2*c2*x_final + 3*c3*x_final**2 + 4*c4*x_final**3))
                    x_err = y_err_val / derivative_x
                    
                    oh_o3s2_c20[idx_y, idx_x] = x_final + 8.69
                    # Add intrinsic fitting error from O3S2-C20 calibration (0.11 dex)
                    fitting_err = 0.11  # dex
                    oh_o3s2_c20_err[idx_y, idx_x] = np.sqrt(x_err**2)
                else:
                    combined_mask[idx_y, idx_x] = False
            else:
                combined_mask[idx_y, idx_x] = False
    
    # Update the combined mask to only include spaxels where we found valid solutions
    combined_mask = combined_mask & np.isfinite(oh_o3s2_c20)
    
    return oh_o3s2_c20, oh_o3s2_c20_err, combined_mask

# Calculate O3S2-C20 metallicity with error propagation
O_H_O3S2_C20, O_H_O3S2_C20_ERR, o3s2_c20_good_mask = calculate_o3s2_c20_metallicity(
    HB4861_FLUX_corr, OIII5006_FLUX_corr, SII6716_FLUX_corr, SII6730_FLUX_corr,
    HB4861_FLUX_ERR, OIII5006_FLUX_ERR, SII6716_FLUX_ERR, SII6730_FLUX_ERR, 
    O_H_D16)
# Set O_H_O3S2_C20 to be nan if outside the range of 7.63 and 9.23
O_H_O3S2_C20 = np.where((O_H_O3S2_C20 < 7.63) | (O_H_O3S2_C20 > 9.23), np.nan, O_H_O3S2_C20)
O_H_O3S2_C20_ERR = np.where((O_H_O3S2_C20 < 7.63) | (O_H_O3S2_C20 > 9.23), np.nan, O_H_O3S2_C20_ERR)

# RS32-C20 (Curti et al. 2020) metallicity calculation function
def calculate_rs32_c20_metallicity(hb4861_flux, ha6563_flux,
                                   oiii5006_flux, sii6716_flux, sii6730_flux,
                                   hb4861_flux_err, ha6563_flux_err,
                                   oiii5006_flux_err, sii6716_flux_err, sii6730_flux_err,
                                   oh_d16_sf,
                                   coeffs=(-0.054, -2.546, -1.970, 0.082, 0.222)):
    """
    RS32–C20 calibration (Curti+2020; user-provided coefficients) with error propagation:
      RS32 = log10( [OIII]/Hβ + ([SII]6716+6730)/Hα )
      Let y = RS32 and x = (12+log(O/H)) - 8.69
      Then: y = c0 + c1 x + c2 x^2 + c3 x^3 + c4 x^4
      Solve per spaxel for x, return 12+log(O/H) = x + 8.69
    """
    c0, c1, c2, c3, c4 = coeffs

    # Good spaxels: use basic finite checks on emission lines
    good_mask = np.ones_like(hb4861_flux, dtype=bool)

    pos = (
        np.isfinite(hb4861_flux) & np.isfinite(ha6563_flux) &
        np.isfinite(oiii5006_flux) & np.isfinite(sii6716_flux) & np.isfinite(sii6730_flux) &
        np.isfinite(hb4861_flux_err) & np.isfinite(ha6563_flux_err) &
        np.isfinite(oiii5006_flux_err) & np.isfinite(sii6716_flux_err) & np.isfinite(sii6730_flux_err) &
        (hb4861_flux > 0) & (ha6563_flux > 0) &
        (oiii5006_flux > 0) & (sii6716_flux > 0) & (sii6730_flux > 0)
    )
    good_mask &= pos

    oh_rs32_c20 = np.full_like(hb4861_flux, np.nan, dtype=float)
    oh_rs32_c20_err = np.full_like(hb4861_flux, np.nan, dtype=float)

    if np.any(good_mask):
        # RS32 (linear inside the log): [OIII]/Hβ + [SII]/Hα
        oiii_hb = oiii5006_flux[good_mask] / hb4861_flux[good_mask]
        sii_total = sii6716_flux[good_mask] + sii6730_flux[good_mask]
        sii_total_err = np.sqrt(sii6716_flux_err[good_mask]**2 + sii6730_flux_err[good_mask]**2)
        sii_ha = sii_total / ha6563_flux[good_mask]
        
        r_lin = oiii_hb + sii_ha
        r_lin = np.where(r_lin > 0, r_lin, np.nan)
        y = np.log10(r_lin)
        
        # Calculate errors for the line ratios using error propagation
        oiii_hb_err = bpt_error_propagation(oiii5006_flux[good_mask], hb4861_flux[good_mask], 
                                            oiii5006_flux_err[good_mask], hb4861_flux_err[good_mask])
        sii_ha_err = bpt_error_propagation(sii_total, ha6563_flux[good_mask], 
                                           sii_total_err, ha6563_flux_err[good_mask])
        
        # Error for RS32 = log10(OIII/Hb + SII/Ha)
        # For f = A + B, df = sqrt(dA^2 + dB^2)
        # For g = log10(f), dg = (1/ln(10)) * df/f
        r_lin_err = np.sqrt(oiii_hb_err**2 + sii_ha_err**2)
        y_err = (1/np.log(10)) * (r_lin_err / r_lin)

        # Solve quartic per valid pixel:
        # c4*x^4 + c3*x^3 + c2*x^2 + c1*x + (c0 - y) = 0
        idxs = np.argwhere(good_mask)
        for idx_in_good, (iy, ix) in enumerate(idxs):
            y_val = y[idx_in_good]
            y_err_val = y_err[idx_in_good]
            if not np.isfinite(y_val):
                continue
            roots = np.roots([c4, c3, c2, c1, (c0 - y_val)])
            real = roots[np.isreal(roots)].real
            if real.size:
                # Sensible metallicity range for Te-anchored scales:
                # 12+log(O/H) ~ 8.0–8.9 ⇒ x = (12+log(O/H))-8.69 ∈ [-0.7, +0.3]
                phys = real[(real >= -0.7) & (real <= 0.3)]
                cand = phys if phys.size else real  # fallback if nothing in phys range
                # choose the root that best reproduces y
                y_pred = c0 + c1*cand + c2*cand**2 + c3*cand**3 + c4*cand**4
                x_final = cand[np.argmin(np.abs(y_pred - y_val))]
                
                # Calculate error in x using derivative approach
                # For equation f(x,y) = y - (c0 + c1*x + c2*x^2 + c3*x^3 + c4*x^4) = 0
                # df/dx = -(c1 + 2*c2*x + 3*c3*x^2 + 4*c4*x^3), df/dy = 1
                # x_err = |df/dy| * y_err / |df/dx| = y_err / (|c1| + |2*c2*x| + |3*c3*x^2| + |4*c4*x^3|)
                derivative_x = (np.abs(c1 + 2*c2*x_final + 3*c3*x_final**2 + 4*c4*x_final**3))
                x_err = y_err_val / derivative_x
                
                oh_rs32_c20[iy, ix] = x_final + 8.69
                # Add intrinsic fitting error from RS32-C20 calibration (0.08 dex)
                fitting_err = 0.08  # dex
                oh_rs32_c20_err[iy, ix] = np.sqrt(x_err**2)

    combined_mask = good_mask & np.isfinite(oh_rs32_c20)
    return oh_rs32_c20, oh_rs32_c20_err, combined_mask

# Calculate RS32-C20 metallicity
O_H_RS32_C20, O_H_RS32_C20_ERR, rs32_c20_good_mask = calculate_rs32_c20_metallicity(HB4861_FLUX_corr, HA6562_FLUX_corr, 
                                                                  OIII5006_FLUX_corr, SII6716_FLUX_corr, SII6730_FLUX_corr,
                                                                  HB4861_FLUX_ERR, HA6562_FLUX_ERR,
                                                                  OIII5006_FLUX_ERR, SII6716_FLUX_ERR, SII6730_FLUX_ERR,
                                                                  None)
# Set O_H_RS32_C20 to be nan if outside the range of 7.63 and 9.23
O_H_RS32_C20 = np.where((O_H_RS32_C20 < 7.63) | (O_H_RS32_C20 > 9.23), np.nan, O_H_RS32_C20)
O_H_RS32_C20_ERR = np.where((O_H_RS32_C20 < 7.63) | (O_H_RS32_C20 > 9.23), np.nan, O_H_RS32_C20_ERR)

# R3-C20 (Curti et al. 2020) metallicity calculation function
def calculate_r3_c20_metallicity(hb4861_flux, hb4861_flux_err,
                                 oiii5006_flux, oiii5006_flux_err,
                                 oh_d16_sf,
                                 coeffs=(-0.277, -3.549, -3.593, -0.981),
                                 fitting_error=0.07):
    """
    R3–C20 calibration (Curti+2020; user-provided coefficients):
      R3 = log10( [OIII]5007 / Hβ )
      Let y = R3 and x = (12+log(O/H)) - 8.69
      Then: y = c0 + c1 x + c2 x^2 + c3 x^3
      Solve per spaxel for x, return 12+log(O/H) = x + 8.69
    """
    c0, c1, c2, c3 = coeffs

    # Good spaxels: use basic finite checks on emission lines
    good_mask = np.ones_like(hb4861_flux, dtype=bool)

    pos = (
        np.isfinite(hb4861_flux) & np.isfinite(oiii5006_flux) &
        (hb4861_flux > 0) & (oiii5006_flux > 0) &
        np.isfinite(hb4861_flux_err) & np.isfinite(oiii5006_flux_err) &
        (hb4861_flux_err > 0) & (oiii5006_flux_err > 0)
    )
    good_mask &= pos

    oh_r3_c20 = np.full_like(hb4861_flux, np.nan, dtype=float)
    oh_r3_c20_err = np.full_like(hb4861_flux, np.nan, dtype=float)

    if np.any(good_mask):
        # R3 = log10([OIII]/Hβ) and its error
        r_lin = (oiii5006_flux[good_mask] / hb4861_flux[good_mask])
        r_lin = np.where(r_lin > 0, r_lin, np.nan)
        y = np.log10(r_lin)
        
        # Calculate error in R3 using BPT error propagation
        r3_error = bpt_error_propagation(
            oiii5006_flux[good_mask], hb4861_flux[good_mask],
            oiii5006_flux_err[good_mask], hb4861_flux_err[good_mask]
        )

        # Solve cubic per valid pixel and calculate error:
        # c3*x^3 + c2*x^2 + c1*x + (c0 - y) = 0
        idxs = np.argwhere(good_mask)
        for idx, ((iy, ix), y_val) in enumerate(zip(idxs, y)):
            if not np.isfinite(y_val):
                continue
            roots = np.roots([c3, c2, c1, (c0 - y_val)])
            real = roots[np.isreal(roots)].real
            if real.size:
                # Sensible metallicity range for Te-anchored scales:
                # 12+log(O/H) ~ 8.0–8.9 ⇒ x = (12+log(O/H))-8.69 ∈ [-0.7, +0.3]
                phys = real[(real >= -0.7) & (real <= 0.3)]
                cand = phys if phys.size else real  # fallback if nothing in phys range
                # choose the root that best reproduces y
                y_pred = c0 + c1*cand + c2*cand**2 + c3*cand**3
                x_final = cand[np.argmin(np.abs(y_pred - y_val))]
                oh_r3_c20[iy, ix] = x_final + 8.69
                
                # Error propagation: derivative of polynomial with respect to y
                # dy/dx = c1 + 2*c2*x + 3*c3*x^2
                derivative_y = np.abs(c1 + 2*c2*x_final + 3*c3*x_final**2)
                
                if derivative_y > 0:
                    # dx/dy = 1/(dy/dx)
                    derivative_x = 1.0 / derivative_y
                    
                    # Error in metallicity from observational error in R3
                    obs_error = derivative_x * r3_error[idx]
                    
                    # Combine observational error with fitting error
                    total_error = np.sqrt(obs_error**2)
                    oh_r3_c20_err[iy, ix] = total_error

    combined_mask = good_mask & np.isfinite(oh_r3_c20)
    return oh_r3_c20, oh_r3_c20_err, combined_mask

# Calculate R3-C20 metallicity
O_H_R3_C20, O_H_R3_C20_ERR, r3_c20_good_mask = calculate_r3_c20_metallicity(HB4861_FLUX_corr, HB4861_FLUX_ERR,
                                                                             OIII5006_FLUX_corr, OIII5006_FLUX_ERR,
                                                                             None)
# Set O_H_R3_C20 to be nan if outside the range of 7.63 and 9.23
O_H_R3_C20 = np.where((O_H_R3_C20 < 7.63) | (O_H_R3_C20 > 9.23), np.nan, O_H_R3_C20)
O_H_R3_C20_ERR = np.where((O_H_R3_C20 < 7.63) | (O_H_R3_C20 > 9.23), np.nan, O_H_R3_C20_ERR)

# N2-C20 (Curti et al. 2020) metallicity calculation function
def calculate_n2_c20_metallicity(ha6563_flux, ha6563_flux_err,
                                 nii6584_flux, nii6584_flux_err,
                                 oh_d16_sf,
                                 coeffs=(-0.489, 1.513, -2.554, -5.293, -2.867),
                                 fitting_error=0.10):
    """
    N2–C20 calibration (Curti+2020; user-provided coefficients):
      N2 = log10( [NII]6584 / Hα )
      Let y = N2 and x = (12+log(O/H)) - 8.69
      Then: y = c0 + c1 x + c2 x^2 + c3 x^3 + c4 x^4
      Solve per spaxel for x, return 12+log(O/H) = x + 8.69

    Selection rule (as requested):
      • Get ALL (near-)real roots of the quartic.
      • Keep only roots with x ∈ [-0.7, 0.3].
      • If multiple such roots exist, pick the SMALLEST one.
      • If none exist, discard the spaxel (leave NaN).
      • No post-hoc clipping.
    """
    c0, c1, c2, c3, c4 = coeffs

    # Good spaxels: use basic finite checks on emission lines
    good_mask = np.ones_like(ha6563_flux, dtype=bool)

    pos = (
        np.isfinite(ha6563_flux) & np.isfinite(nii6584_flux) &
        (ha6563_flux > 0) & (nii6584_flux > 0) &
        np.isfinite(ha6563_flux_err) & np.isfinite(nii6584_flux_err) &
        (ha6563_flux_err > 0) & (nii6584_flux_err > 0)
    )
    good_mask &= pos

    oh_n2_c20 = np.full_like(ha6563_flux, np.nan, dtype=float)
    oh_n2_c20_err = np.full_like(ha6563_flux, np.nan, dtype=float)

    if np.any(good_mask):
        # N2 (linear inside the log): [NII]6584 / Hα
        n2_lin = nii6584_flux[good_mask] / ha6563_flux[good_mask]
        n2_lin = np.where(n2_lin > 0, n2_lin, np.nan)
        y = np.log10(n2_lin)
        
        # Calculate error in N2 using BPT error propagation
        n2_error = bpt_error_propagation(
            nii6584_flux[good_mask], ha6563_flux[good_mask],
            nii6584_flux_err[good_mask], ha6563_flux_err[good_mask]
        )

        idxs = np.argwhere(good_mask)

        # tolerances for "almost-real" roots and range comparison
        REAL_ATOL = 1e-8
        RANGE_EPS = 0.0  # exact bounds as requested

        for idx, ((iy, ix), y_val) in enumerate(zip(idxs, y)):
            if not np.isfinite(y_val):
                continue

            roots = np.roots([c4, c3, c2, c1, (c0 - y_val)])

            # treat tiny-imag roots as real
            realish = roots[np.abs(roots.imag) <= REAL_ATOL].real
            if realish.size == 0:
                # fall back to the least-imag root if imag part is tiny-ish
                k = np.argmin(np.abs(roots.imag))
                if np.abs(roots[k].imag) <= 1e-6:
                    realish = np.array([roots[k].real])

            if realish.size == 0:
                continue  # no usable real roots

            # STRICT selection inside [-0.7, 0.3]
            in_rng = realish[(realish >= -0.7 - RANGE_EPS) & (realish <= 0.3 + RANGE_EPS)]
            if in_rng.size == 0:
                # No in-range root → discard this spaxel
                continue

            # If multiple, pick the smallest one
            x_final = np.min(in_rng)
            oh_n2_c20[iy, ix] = x_final + 8.69
            
            # Error propagation: derivative of 4th-order polynomial with respect to y
            # dy/dx = c1 + 2*c2*x + 3*c3*x^2 + 4*c4*x^3
            derivative_y = np.abs(c1 + 2*c2*x_final + 3*c3*x_final**2 + 4*c4*x_final**3)
            
            if derivative_y > 0:
                # dx/dy = 1/(dy/dx)
                derivative_x = 1.0 / derivative_y
                
                # Error in metallicity from observational error in N2
                obs_error = derivative_x * n2_error[idx]
                
                # Combine observational error with fitting error
                total_error = np.sqrt(obs_error**2)
                oh_n2_c20_err[iy, ix] = total_error

    combined_mask = good_mask & np.isfinite(oh_n2_c20)
    return oh_n2_c20, oh_n2_c20_err, combined_mask

# Calculate N2-C20 metallicity
O_H_N2_C20, O_H_N2_C20_ERR, n2_c20_good_mask = calculate_n2_c20_metallicity(HA6562_FLUX_corr, HA6562_FLUX_ERR,
                                                                             NII6583_FLUX_corr, NII6583_FLUX_ERR,
                                                                             None)
# Set O_H_N2_C20 to be nan if outside the range of 7.63 and 9.23
O_H_N2_C20 = np.where((O_H_N2_C20 < 7.63) | (O_H_N2_C20 > 9.23), np.nan, O_H_N2_C20)
O_H_N2_C20_ERR = np.where((O_H_N2_C20 < 7.63) | (O_H_N2_C20 > 9.23), np.nan, O_H_N2_C20_ERR)

def s2_error_propagation(sii6716_flux, sii6716_flux_err, sii6730_flux, sii6730_flux_err, ha6563_flux, ha6563_flux_err):
    """Calculate propagated error for log10(([SII]6716 + [SII]6730) / Hα)"""
    # Error in numerator (sum of [SII] lines)
    numerator = sii6716_flux + sii6730_flux
    numerator_err = np.sqrt(sii6716_flux_err**2 + sii6730_flux_err**2)
    
    # Error propagation for log10(numerator/denominator)
    ratio_rel_err = np.sqrt((numerator_err / numerator)**2 + (ha6563_flux_err / ha6563_flux)**2)
    log_ratio_err = ratio_rel_err / np.log(10)
    return log_ratio_err

# S2-C20 (Curti et al. 2020) metallicity calculation function
def calculate_s2_c20_metallicity(ha6563_flux, ha6563_flux_err,
                                 sii6716_flux, sii6716_flux_err,
                                 sii6730_flux, sii6730_flux_err,
                                 oh_d16_sf,
                                 coeffs=(-0.442, -0.360, -6.271, -8.339, -3.559),
                                 fitting_error=0.06):
    """
    S2–C20 calibration (Curti+2020; user-provided coefficients):
      S2 = log10( ([SII]6716 + [SII]6730) / Hα )
      Let y = S2 and x = (12+log(O/H)) - 8.69
      Then: y = c0 + c1 x + c2 x^2 + c3 x^3 + c4 x^4
      Solve per spaxel for x, return 12+log(O/H) = x + 8.69

    Root selection (strict):
      • Collect all (near-)real roots.
      • Keep only roots with x ∈ [-0.7, 0.3].
      • If multiple, choose the smallest.
      • If none in range, discard spaxel (NaN).
      • No post-hoc clipping.
    """
    c0, c1, c2, c3, c4 = coeffs

    # Good spaxels: use basic finite checks on emission lines
    good_mask = np.ones_like(ha6563_flux, dtype=bool)

    pos = (
        np.isfinite(ha6563_flux) & np.isfinite(sii6716_flux) & np.isfinite(sii6730_flux) &
        (ha6563_flux > 0) & (sii6716_flux > 0) & (sii6730_flux > 0) &
        np.isfinite(ha6563_flux_err) & np.isfinite(sii6716_flux_err) & np.isfinite(sii6730_flux_err) &
        (ha6563_flux_err > 0) & (sii6716_flux_err > 0) & (sii6730_flux_err > 0)
    )
    good_mask &= pos

    oh_s2_c20 = np.full_like(ha6563_flux, np.nan, dtype=float)
    oh_s2_c20_err = np.full_like(ha6563_flux, np.nan, dtype=float)

    if np.any(good_mask):
        # S2 (linear inside the log): ([SII]6716+6730)/Hα
        s2_lin = (sii6716_flux[good_mask] + sii6730_flux[good_mask]) / ha6563_flux[good_mask]
        s2_lin = np.where(s2_lin > 0, s2_lin, np.nan)
        y = np.log10(s2_lin)
        
        # Calculate error in S2 using specialized error propagation
        s2_error = s2_error_propagation(
            sii6716_flux[good_mask], sii6716_flux_err[good_mask],
            sii6730_flux[good_mask], sii6730_flux_err[good_mask],
            ha6563_flux[good_mask], ha6563_flux_err[good_mask]
        )

        idxs = np.argwhere(good_mask)

        REAL_ATOL = 1e-8  # accept roots with tiny imaginary part as real

        for idx, ((iy, ix), y_val) in enumerate(zip(idxs, y)):
            if not np.isfinite(y_val):
                continue

            # Solve: c4*x^4 + c3*x^3 + c2*x^2 + c1*x + (c0 - y) = 0
            roots = np.roots([c4, c3, c2, c1, (c0 - y_val)])

            # Treat tiny-imag roots as real
            realish = roots[np.abs(roots.imag) <= REAL_ATOL].real
            if realish.size == 0:
                # fallback: least-imag root if imag part is still tiny-ish
                k = np.argmin(np.abs(roots.imag))
                if np.abs(roots[k].imag) <= 1e-6:
                    realish = np.array([roots[k].real])
                else:
                    continue  # no usable real roots

            # STRICT in-range selection: x ∈ [-0.7, 0.3]
            in_range = realish[(realish >= -0.7) & (realish <= 0.3)]
            if in_range.size == 0:
                continue  # discard spaxel if no valid root in range

            # If multiple, pick the smallest
            x_final = np.min(in_range)
            oh_s2_c20[iy, ix] = x_final + 8.69
            
            # Error propagation: derivative of 4th-order polynomial with respect to y
            # dy/dx = c1 + 2*c2*x + 3*c3*x^2 + 4*c4*x^3
            derivative_y = np.abs(c1 + 2*c2*x_final + 3*c3*x_final**2 + 4*c4*x_final**3)
            
            if derivative_y > 0:
                # dx/dy = 1/(dy/dx)
                derivative_x = 1.0 / derivative_y
                
                # Error in metallicity from observational error in S2
                obs_error = derivative_x * s2_error[idx]
                
                # Combine observational error with fitting error
                total_error = np.sqrt(obs_error**2)
                oh_s2_c20_err[iy, ix] = total_error

    combined_mask = good_mask & np.isfinite(oh_s2_c20)
    return oh_s2_c20, oh_s2_c20_err, combined_mask

# Calculate S2-C20 metallicity
O_H_S2_C20, O_H_S2_C20_ERR, s2_c20_good_mask = calculate_s2_c20_metallicity(HA6562_FLUX_corr, HA6562_FLUX_ERR,
                                                                             SII6716_FLUX_corr, SII6716_FLUX_ERR,
                                                                             SII6730_FLUX_corr, SII6730_FLUX_ERR,
                                                                             None)
# Set O_H_S2_C20 to be nan if outside the range of 7.63 and 9.23
O_H_S2_C20 = np.where((O_H_S2_C20 < 7.63) | (O_H_S2_C20 > 9.23), np.nan, O_H_S2_C20)
O_H_S2_C20_ERR = np.where((O_H_S2_C20 < 7.63) | (O_H_S2_C20 > 9.23), np.nan, O_H_S2_C20_ERR)

# Combined C20 metallicity calculation function

def calculate_combined_c20_metallicity(gal):
    """
    Calculate combined C20 metallicity by selecting the method with smallest error for each spaxel.
    
    For each spaxel, we:
    1. Calculate metallicity and error for all 6 C20 methods
    2. Select the method with the smallest error
    3. Return the corresponding metallicity value
    
    Returns:
        oh_combined_c20: Combined metallicity map
        oh_combined_c20_err: Combined error map
        method_map: Map showing which method was used for each spaxel (0-5)
        combined_mask: Combined valid spaxel mask
    """
    # Load all required fluxes
    with fits.open(f'{gal}_SPATIAL_BINNING_maps_extended.fits') as h:
        sigM = h['LOGMASS_SURFACE_DENSITY'].data
    with fits.open(f'{gal}_gas_BIN_maps_extended.fits') as h:
        sigSFR = h['LOGSFR_SURFACE_DENSITY_SF'].data
        
        # Load emission line fluxes
        hb4861_flux = h['HB4861_FLUX_corr'].data
        oiii5006_flux = h['OIII5006_FLUX_corr'].data
        sii6716_flux = h['SII6716_FLUX_corr'].data
        sii6730_flux = h['SII6730_FLUX_corr'].data
        
        # Load emission line flux errors
        hb4861_flux_err = h['HB4861_FLUX_ERR'].data
        oiii5006_flux_err = h['OIII5006_FLUX_ERR'].data
        sii6716_flux_err = h['SII6716_FLUX_ERR'].data
        sii6730_flux_err = h['SII6730_FLUX_ERR'].data
        
        # Load Halpha and NII with multiple naming conventions
        ha_key_candidates = ('HA6563_FLUX_corr', 'HA6562_FLUX_corr', 'HALPHA6563_FLUX_corr', 'HALPHA_FLUX_corr')
        ha_err_candidates = ('HA6563_FLUX_ERR', 'HA6562_FLUX_ERR', 'HALPHA6563_FLUX_ERR', 'HALPHA_FLUX_ERR')
        nii_key_candidates = ('NII6584_FLUX_corr', 'NII6583_FLUX_corr', 'NII6584_FLUX', 'NII6583_FLUX')
        nii_err_candidates = ('NII6584_FLUX_ERR', 'NII6583_FLUX_ERR')
        
        ha6563_flux = None
        for k in ha_key_candidates:
            if k in h:
                ha6563_flux = h[k].data
                break
        
        ha6563_flux_err = None
        for k in ha_err_candidates:
            if k in h:
                ha6563_flux_err = h[k].data
                break
                
        nii6584_flux = None
        for k in nii_key_candidates:
            if k in h:
                nii6584_flux = h[k].data
                break
                
        nii6584_flux_err = None
        for k in nii_err_candidates:
            if k in h:
                nii6584_flux_err = h[k].data
                break
        
        # Load reference O/H data
        try:
            oh_d16_sf = h['O_H_D16_SF'].data
        except KeyError:
            oh_d16_sf = None
    
    # Calculate metallicity for all 6 methods
    print(f"Calculating all 6 C20 metallicities for {gal}...")
    
    # Method 0: O3N2-C20
    oh_o3n2_c20, oh_o3n2_c20_err, mask_o3n2 = calculate_o3n2_c20_metallicity(
        hb4861_flux, oiii5006_flux, nii6584_flux, ha6563_flux,
        hb4861_flux_err, oiii5006_flux_err, nii6584_flux_err, ha6563_flux_err, oh_d16_sf
    )
    
    # Method 1: O3S2-C20
    oh_o3s2_c20, oh_o3s2_c20_err, mask_o3s2 = calculate_o3s2_c20_metallicity(
        hb4861_flux, oiii5006_flux, sii6716_flux, sii6730_flux,
        hb4861_flux_err, oiii5006_flux_err, sii6716_flux_err, sii6730_flux_err, oh_d16_sf
    )
    
    # Method 2: RS32-C20
    oh_rs32_c20, oh_rs32_c20_err, mask_rs32 = calculate_rs32_c20_metallicity(
        hb4861_flux, ha6563_flux, oiii5006_flux, sii6716_flux, sii6730_flux,
        hb4861_flux_err, ha6563_flux_err, oiii5006_flux_err, sii6716_flux_err, sii6730_flux_err, oh_d16_sf
    )
    
    # Method 3: R3-C20
    oh_r3_c20, oh_r3_c20_err, mask_r3 = calculate_r3_c20_metallicity(
        hb4861_flux, hb4861_flux_err, oiii5006_flux, oiii5006_flux_err, oh_d16_sf
    )
    
    # Method 4: N2-C20
    oh_n2_c20, oh_n2_c20_err, mask_n2 = calculate_n2_c20_metallicity(
        ha6563_flux, ha6563_flux_err, nii6584_flux, nii6584_flux_err, oh_d16_sf
    )
    
    # Method 5: S2-C20
    oh_s2_c20, oh_s2_c20_err, mask_s2 = calculate_s2_c20_metallicity(
        ha6563_flux, ha6563_flux_err, sii6716_flux, sii6716_flux_err, sii6730_flux, sii6730_flux_err, oh_d16_sf
    )
    
    # Stack all metallicities and errors
    all_metallicities = np.stack([oh_o3n2_c20, oh_o3s2_c20, oh_rs32_c20, oh_r3_c20, oh_n2_c20, oh_s2_c20], axis=0)
    all_errors = np.stack([oh_o3n2_c20_err, oh_o3s2_c20_err, oh_rs32_c20_err, oh_r3_c20_err, oh_n2_c20_err, oh_s2_c20_err], axis=0)
    all_masks = np.stack([mask_o3n2, mask_o3s2, mask_rs32, mask_r3, mask_n2, mask_s2], axis=0)
    
    # Initialize output arrays
    oh_combined_c20 = np.full_like(hb4861_flux, np.nan)
    oh_combined_c20_err = np.full_like(hb4861_flux, np.nan)
    method_map = np.full_like(hb4861_flux, -1, dtype=int)  # -1 indicates no valid method
    
    # For each spaxel, find the method with the smallest error
    for i in range(hb4861_flux.shape[0]):
        for j in range(hb4861_flux.shape[1]):
            # Get valid methods for this spaxel
            valid_methods = all_masks[:, i, j]
            
            if np.any(valid_methods):
                # Get errors for valid methods only
                valid_errors = all_errors[valid_methods, i, j]
                valid_metallicities = all_metallicities[valid_methods, i, j]
                
                # Find method with minimum error
                min_error_idx = np.nanargmin(valid_errors)
                
                # Map back to original method index
                method_indices = np.where(valid_methods)[0]
                best_method = method_indices[min_error_idx]
                
                # Store results
                oh_combined_c20[i, j] = valid_metallicities[min_error_idx]
                oh_combined_c20_err[i, j] = valid_errors[min_error_idx]
                method_map[i, j] = best_method
    
    # Create combined mask
    combined_mask = np.isfinite(oh_combined_c20)
    
    return oh_combined_c20, oh_combined_c20_err, method_map, combined_mask

# Calculate combined C20 metallicity
O_H_COMBINED_C20, O_H_COMBINED_C20_ERR, combined_c20_method_map, combined_c20_good_mask = calculate_combined_c20_metallicity(gal)
# Set combined C20 to be nan if outside the range of 7.63 and 9.23
O_H_COMBINED_C20 = np.where((O_H_COMBINED_C20 < 7.63) | (O_H_COMBINED_C20 > 9.23), np.nan, O_H_COMBINED_C20)
O_H_COMBINED_C20_ERR = np.where((O_H_COMBINED_C20 < 7.63) | (O_H_COMBINED_C20 > 9.23), np.nan, O_H_COMBINED_C20_ERR)

print(f"Combined C20 metallicity: median = {np.nanmedian(O_H_COMBINED_C20):.3f}, range = ({np.nanmin(O_H_COMBINED_C20):.3f}, {np.nanmax(O_H_COMBINED_C20):.3f})")
print(f"Combined C20 method usage: O3N2={np.sum(combined_c20_method_map==0)}, O3S2={np.sum(combined_c20_method_map==1)}, RS32={np.sum(combined_c20_method_map==2)}, R3={np.sum(combined_c20_method_map==3)}, N2={np.sum(combined_c20_method_map==4)}, S2={np.sum(combined_c20_method_map==5)}")

# # For D16 and PG16, select the finite values in both maps (O3N2-M13, N2-M13, O3N2-PP04, N2-PP04, O3N2-C20, O3S2-C20, RS32-C20, R3-C20, N2-C20 and S2-C20 will be calculated where D16/PG16 are valid)
# valid_mask = np.isfinite(O_H_D16) & np.isfinite(O_H_PG16) & np.isfinite(O_H_O3N2_M13) & np.isfinite(O_H_N2_M13) & np.isfinite(O_H_O3N2_PP04) & np.isfinite(O_H_N2_PP04) & np.isfinite(O_H_O3N2_C20) & np.isfinite(O_H_O3S2_C20) & np.isfinite(O_H_RS32_C20) & np.isfinite(O_H_R3_C20) & np.isfinite(O_H_N2_C20) & np.isfinite(O_H_S2_C20)
# O_H_D16 = np.where(valid_mask, O_H_D16, np.nan)
# O_H_PG16 = np.where(valid_mask, O_H_PG16, np.nan)
# # Apply the same mask to O3N2-M13, N2-M13, O3N2-PP04, N2-PP04, O3N2-C20, O3S2-C20, RS32-C20, R3-C20, N2-C20 and S2-C20 to ensure consistency
# O_H_O3N2_M13 = np.where(valid_mask, O_H_O3N2_M13, np.nan)
# O_H_N2_M13 = np.where(valid_mask, O_H_N2_M13, np.nan)
# O_H_O3N2_PP04 = np.where(valid_mask, O_H_O3N2_PP04, np.nan)
# O_H_N2_PP04 = np.where(valid_mask, O_H_N2_PP04, np.nan)
# O_H_O3N2_C20 = np.where(valid_mask, O_H_O3N2_C20, np.nan)
# O_H_O3S2_C20 = np.where(valid_mask, O_H_O3S2_C20, np.nan)
# O_H_RS32_C20 = np.where(valid_mask, O_H_RS32_C20, np.nan)
# O_H_R3_C20 = np.where(valid_mask, O_H_R3_C20, np.nan)
# O_H_N2_C20 = np.where(valid_mask, O_H_N2_C20, np.nan)
# O_H_S2_C20 = np.where(valid_mask, O_H_S2_C20, np.nan)
# O_H_COMBINED_C20 = np.where(valid_mask, O_H_COMBINED_C20, np.nan)

# ------------------------------------------------------------------
# End of Metallicity [O/H] calculation (12+log(O/H)) using different methods
# ------------------------------------------------------------------

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

# 2. Read galaxy inclination and calculate correction factor
if apply_inclination_correction:
    galaxy_inclination = read_galaxy_inclination(gal)
    if galaxy_inclination is not None:
        inclination_rad = np.deg2rad(galaxy_inclination)
        cos_inclination = np.cos(inclination_rad)
        print(f"Galaxy {gal} inclination: {galaxy_inclination}° (cos θ = {cos_inclination:.3f})")
        print(f"Inclination correction ENABLED: applying cos(θ) = {cos_inclination:.3f}")
    else:
        cos_inclination = 1.0
        print(f"No inclination data found for {gal}, using cos_inclination = 1.0")
else:
    cos_inclination = 1.0
    print(f"Inclination correction DISABLED: using cos_inclination = 1.0")

# 3. Convert SFR to surface density with inclination correction
SFR_surface_density_map = SFR_map / pixel_area_kpc.value
SFR_surface_density_map_corrected = SFR_surface_density_map * cos_inclination  # Apply inclination correction
# SFR_surface_density_map_corrected = SFR_surface_density_map 

# 4. Convert to log10 scale
LOG_SFR_surface_density_map = np.log10(SFR_surface_density_map_corrected)

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
# 6.  Masks: classified or not?  
# ------------------------------------------------------------------

# now I want to use the error bars to determine the mask called mask_classified, 
# which is for each point, its value +/- errorbars, are all still inside the same region of on the BPT. 
# These regions are HII+Comp and AGN for NII BPT; and HII, LINER and Seyfert for SII BPT. 

# Error propogation for BPT diagrams (sigma of log_10(numerator/denominator))
# def bpt_error_propagation(numerator, denominator, numerator_err, denominator_err):
#     """
#     Calculate the propagated error for the BPT ratio log10(numerator/denominator).
    
#     Parameters:
#     numerator (np.ndarray): The numerator values.
#     denominator (np.ndarray): The denominator values.
#     numerator_err (np.ndarray): The error in the numerator.
#     denominator_err (np.ndarray): The error in the denominator.
    
#     Returns:
#     np.ndarray: The propagated error for the BPT ratio.
#     """
#     # Avoid division by zero
#     with np.errstate(divide='ignore', invalid='ignore'):
#         ratio = numerator / denominator
#         log_ratio = np.log10(ratio)
#         log_ratio_err = 1/(np.log(10)) * np.sqrt((numerator_err / numerator)**2 + (denominator_err / denominator)**2)
#         return log_ratio_err
    
# Calculate the errors for the BPT ratios
logN2_err = bpt_error_propagation(NII6583_FLUX, HA6562_FLUX, 
                                   NII6583_FLUX_ERR, HA6562_FLUX_ERR)
logS2_err = bpt_error_propagation(SII6716_FLUX + SII6730_FLUX, HA6562_FLUX, 
                                    np.sqrt(SII6716_FLUX_ERR**2 + SII6730_FLUX_ERR**2), HA6562_FLUX_ERR)
logO3_err = bpt_error_propagation(OIII5006_FLUX, HB4861_FLUX, 
                                   OIII5006_FLUX_ERR, HB4861_FLUX_ERR)

mask_N2_left, mask_S2_left = apply_bpt_masks(logN2=logN2-logN2_err, logS2=logS2-logS2_err, logO3=logO3)
mask_N2_right, mask_S2_right = apply_bpt_masks(logN2=logN2+logN2_err, logS2=logS2+logS2_err, logO3=logO3)
mask_N2_down, mask_S2_down = apply_bpt_masks(logN2=logN2, logS2=logS2, logO3=logO3-logO3_err)
mask_N2_up, mask_S2_up = apply_bpt_masks(logN2=logN2, logS2=logS2, logO3=logO3+logO3_err)

mask_N2_left_HII, mask_N2_left_Comp, mask_N2_left_AGN = mask_N2_left
mask_N2_right_HII, mask_N2_right_Comp, mask_N2_right_AGN = mask_N2_right
mask_N2_down_HII, mask_N2_down_Comp, mask_N2_down_AGN = mask_N2_down
mask_N2_up_HII, mask_N2_up_Comp, mask_N2_up_AGN = mask_N2_up
mask_S2_left_HII, mask_S2_left_LINER, mask_S2_left_Seyfert = mask_S2_left
mask_S2_right_HII, mask_S2_right_LINER, mask_S2_right_Seyfert = mask_S2_right
mask_S2_down_HII, mask_S2_down_LINER, mask_S2_down_Seyfert = mask_S2_down
mask_S2_up_HII, mask_S2_up_LINER, mask_S2_up_Seyfert = mask_S2_up

mask_classified_N2_HII_Comp = ((mask_N2_left_HII | mask_N2_left_Comp) & 
                               (mask_N2_right_HII | mask_N2_right_Comp) & 
                                 (mask_N2_down_HII | mask_N2_down_Comp) &
                                 (mask_N2_up_HII | mask_N2_up_Comp))
mask_classified_N2_AGN = (mask_N2_left_AGN & mask_N2_right_AGN &
                          mask_N2_down_AGN & mask_N2_up_AGN)
mask_classified_S2_HII = (mask_S2_left_HII & mask_S2_right_HII &
                          mask_S2_down_HII & mask_S2_up_HII)
mask_classified_S2_LINER = (mask_S2_left_LINER & mask_S2_right_LINER &
                            mask_S2_down_LINER & mask_S2_up_LINER)
mask_classified_S2_Seyfert = (mask_S2_left_Seyfert & mask_S2_right_Seyfert &
                             mask_S2_down_Seyfert & mask_S2_up_Seyfert)

# NII classified and unclassified masks
mask_N2_classified = (mask_classified_N2_HII_Comp | mask_classified_N2_AGN)
mask_N2_unclassified = ~mask_N2_classified
# SII classified and unclassified masks
mask_S2_classified = (mask_classified_S2_HII | mask_classified_S2_LINER | mask_classified_S2_Seyfert)
mask_S2_unclassified = ~mask_S2_classified

# ------------------------------------------------------------------
# 6.  Masks: for [NII] BPT
# ------------------------------------------------------------------

# Halpha detected
HA_detected = HA6562_QC_good
# Halpha not detected
HA_not_detected = HA6562_QC_bad

# Halpha detected, Hbeta detected
HA_detected_HB_detected = HA6562_QC_good & HB4861_QC_good
# Halpha detected, Hbeta not detected
HA_detected_HB_not_detected = HA6562_QC_good & HB4861_QC_bad

# Halpha detected, Hbeta detected, NII detected
HA_detected_HB_detected_NII_detected = HA6562_QC_good & HB4861_QC_good & NII6583_QC_good
# Halpha detected, Hbeta detected, NII not detected
HA_detected_HB_detected_NII_not_detected = HA6562_QC_good & HB4861_QC_good & NII6583_QC_bad

# Halpha detected, Hbeta detected, NII detected, OIII detected
HA_detected_HB_detected_NII_detected_OIII_detected = (HA6562_QC_good & 
                                                      HB4861_QC_good & 
                                                      NII6583_QC_good &
                                                      OIII5006_QC_good)
# Halpha detected, Hbeta detected, NII detected, OIII not detected
HA_detected_HB_detected_NII_detected_OIII_not_detected = (HA6562_QC_good & 
                                                          HB4861_QC_good & 
                                                          NII6583_QC_good &
                                                          OIII5006_QC_bad)

# Halpha detected, Hbeta detected, NII not detected, OIII detected
HA_detected_HB_detected_NII_not_detected_OIII_detected = (HA6562_QC_good & 
                                                          HB4861_QC_good & 
                                                          NII6583_QC_bad &
                                                          OIII5006_QC_good)
# Halpha detected, Hbeta detected, NII not detected, OIII not detected
HA_detected_HB_detected_NII_not_detected_OIII_not_detected = (HA6562_QC_good & 
                                                              HB4861_QC_good & 
                                                              NII6583_QC_bad &
                                                              OIII5006_QC_bad)


# Flag 4 final cases that we want to track

# definite SF spaxels: or HA_detected_HB_detected & mask_N2_classified & mask_N2_SF
mask_SF_N2 = ((HA_detected_HB_detected_NII_detected_OIII_detected & mask_N2_classified & mask_N2_SF) | 
                    (HA_detected_HB_detected_NII_not_detected_OIII_not_detected & mask_N2_classified & mask_N2_SF) | 
                    (HA_detected_HB_detected_NII_not_detected_OIII_detected & mask_N2_classified & mask_N2_SF) | 
                    (HA_detected_HB_detected_NII_detected_OIII_not_detected & mask_N2_classified & mask_N2_SF))
# get SFR as non-SF: : or HA_detected_HB_detected & mask_N2_classified & mask_N2_nonSF
mask_nonSF_N2 = ((HA_detected_HB_detected_NII_detected_OIII_detected & mask_N2_classified & mask_N2_nonSF) |
              (HA_detected_HB_detected_NII_not_detected_OIII_not_detected & mask_N2_classified & mask_N2_nonSF) |
              (HA_detected_HB_detected_NII_not_detected_OIII_detected & mask_N2_classified & mask_N2_nonSF) |
              (HA_detected_HB_detected_NII_detected_OIII_not_detected & mask_N2_classified & mask_N2_nonSF))
# all the unclassified spaxels: : or HA_detected_HB_detected & mask_N2_unclassified
mask_unclassified_N2 = ((HA_detected_HB_detected_NII_not_detected_OIII_not_detected & mask_N2_unclassified) | 
                       (HA_detected_HB_detected_NII_not_detected_OIII_detected & mask_N2_unclassified) | 
                       (HA_detected_HB_detected_NII_detected_OIII_not_detected & mask_N2_unclassified) | 
                       (HA_detected_HB_detected_NII_detected_OIII_detected & mask_N2_unclassified))
# the rest are upper spaxels: 
mask_upper = (HA_not_detected | HA_detected_HB_not_detected)

# Something else might be useful

# all the classified spaxels
mask_classified_N2 = mask_SF_N2 | mask_nonSF_N2

# ------------------------------------------------------------------
# 7.  Masks: for [SII] BPT
# ------------------------------------------------------------------

# # Halpha detected
# HA_detected = HA6562_QC_good
# # Halpha not detected
# HA_not_detected = HA6562_QC_bad

# # Halpha detected, Hbeta detected
# HA_detected_HB_detected = HA6562_QC_good & HB4861_QC_good
# # Halpha detected, Hbeta not detected
# HA_detected_HB_not_detected = HA6562_QC_good & HB4861_QC_bad

# Halpha detected, Hbeta detected, SII detected
HA_detected_HB_detected_SII_detected = HA6562_QC_good & HB4861_QC_good & (SII6716_QC_good & SII6730_QC_good)
# Halpha detected, Hbeta detected, SII not detected
HA_detected_HB_detected_SII_not_detected = HA6562_QC_good & HB4861_QC_good & ~(SII6716_QC_good & SII6730_QC_good)

# Halpha detected, Hbeta detected, SII detected, OIII detected
HA_detected_HB_detected_SII_detected_OIII_detected = (HA6562_QC_good & 
                                                      HB4861_QC_good & 
                                                      (SII6716_QC_good & SII6730_QC_good) &
                                                      OIII5006_QC_good)
# Halpha detected, Hbeta detected, SII detected, OIII not detected
HA_detected_HB_detected_SII_detected_OIII_not_detected = (HA6562_QC_good & 
                                                          HB4861_QC_good & 
                                                          (SII6716_QC_good & SII6730_QC_good) &
                                                          OIII5006_QC_bad)

# Halpha detected, Hbeta detected, SII not detected, OIII detected
HA_detected_HB_detected_SII_not_detected_OIII_detected = (HA6562_QC_good & 
                                                          HB4861_QC_good & 
                                                          ~(SII6716_QC_good & SII6730_QC_good) &
                                                          OIII5006_QC_good)
# Halpha detected, Hbeta detected, SII not detected, OIII not detected
HA_detected_HB_detected_SII_not_detected_OIII_not_detected = (HA6562_QC_good & 
                                                              HB4861_QC_good & 
                                                              ~(SII6716_QC_good & SII6730_QC_good) &
                                                              OIII5006_QC_bad)


# Flag 4 final cases that we want to track for [SII] BPT

# definite SF spaxels: or HA_detected_HB_detected & mask_S2_classified & mask_S2_SF
mask_SF_S2 = ((HA_detected_HB_detected_SII_detected_OIII_detected & mask_S2_classified & mask_S2_SF) | 
              (HA_detected_HB_detected_SII_not_detected_OIII_not_detected & mask_S2_classified & mask_S2_SF) | 
              (HA_detected_HB_detected_SII_not_detected_OIII_detected & mask_S2_classified & mask_S2_SF) | 
              (HA_detected_HB_detected_SII_detected_OIII_not_detected & mask_S2_classified & mask_S2_SF))
# get SFR as non-SF: or HA_detected_HB_detected & mask_S2_classified & mask_S2_nonSF
mask_nonSF_S2 = ((HA_detected_HB_detected_SII_detected_OIII_detected & mask_S2_classified & mask_S2_nonSF) |
                  (HA_detected_HB_detected_SII_not_detected_OIII_not_detected & mask_S2_classified & mask_S2_nonSF) |
                  (HA_detected_HB_detected_SII_not_detected_OIII_detected & mask_S2_classified & mask_S2_nonSF) |
                  (HA_detected_HB_detected_SII_detected_OIII_not_detected & mask_S2_classified & mask_S2_nonSF))
# all the unconstrained spaxels: or HA_detected_HB_detected & mask_S2_unclassified
mask_unclassified_S2 = ((HA_detected_HB_detected_SII_not_detected_OIII_not_detected & mask_S2_unclassified) | 
                           (HA_detected_HB_detected_SII_not_detected_OIII_detected & mask_S2_unclassified) | 
                           (HA_detected_HB_detected_SII_detected_OIII_not_detected & mask_S2_unclassified) | 
                           (HA_detected_HB_detected_SII_detected_OIII_detected & mask_S2_unclassified))
# # the rest are upper spaxels: 
# mask_upper = (HA_not_detected | HA_detected_HB_not_detected)

# Something else might be useful

# all the constrained spaxels
mask_classified_S2 = mask_SF_S2 | mask_nonSF_S2

# ------------------------------------------------------------------
# 8.  Masks: Combine two BPT
# ------------------------------------------------------------------

# both:

# SF: SF in both N2 and S2 BPT diagrams:
mask_SF_both = mask_SF_N2 & mask_SF_S2
# non-SF: constrained in both N2 and S2 BPT diagrams, but not SF in either or both:
mask_nonSF_both = ((mask_classified_N2 & mask_classified_S2) & ~mask_SF_both)
# Unclassified: unconstrained in either N2 or S2 BPT diagrams: or 
# mask_unclassified_both = (mask_unclassified_N2 | mask_unclassified_S2)
mask_unclassified_both = ((~(mask_classified_N2 & mask_classified_S2)) & HA_detected_HB_detected)
# # Upper: the rest are upper spaxels:
# mask_upper = (HA_not_detected | HA_detected_HB_not_detected)

# Something else might be useful
# all the constrained spaxels: or mask_classified_both = ~mask_unclassified_both
mask_classified_both = (mask_classified_N2 & mask_classified_S2)

# either:

# SF: SF in either N2 or S2 BPT diagrams:
mask_SF_either = mask_SF_N2 | mask_SF_S2
# non-SF: constrained in either N2 or S2 BPT diagrams, but not SF in either or both:
mask_nonSF_either = ((mask_classified_N2 | mask_classified_S2) & ~mask_SF_either)
# Unclassified: unconstrained in either N2 or S2 BPT diagrams: or 
# mask_unclassified_either = (mask_unclassified_N2 & mask_unclassified_S2)
mask_unclassified_either = ((~(mask_classified_N2 | mask_classified_S2)) & HA_detected_HB_detected)
# # Upper: the rest are upper spaxels:
# mask_upper = (HA_not_detected | HA_detected_HB_not_detected)

# Something else might be useful
# all the constrained spaxels: or mask_classified_either = ~mask_unclassified_either
mask_classified_either = (mask_classified_N2 | mask_classified_S2)

# ------------------------------------------------------------------
# 9.  Append Σ_SFR layers (choose 'both' for now, but can be changed to 'either' or just fall back to 'N2' or 'S2')
# ------------------------------------------------------------------

def choose_BPT(choice='both'):
    """
    Choose BPT classification method and return SFR surface density maps, metallicity maps, line maps, and masks.
    
    Parameters:
    choice : str, optional
        BPT classification choice. Options: 'both', 'either', 'N2', 'S2' (default: 'both')
        
    Returns:
    tuple : (SFR_maps, metallicity_maps, line_maps, masks) where:
        - SFR_maps: tuple of four arrays (SF, nonSF, unconstrained, upper) for the chosen BPT method
        - metallicity_maps: tuple of thirteen arrays (O_H_D16_SF, O_H_PG16_SF, O_H_O3N2_M13_SF, O_H_N2_M13_SF, O_H_O3N2_PP04_SF, O_H_N2_PP04_SF, O_H_O3N2_C20_SF, O_H_O3S2_C20_SF, O_H_RS32_C20_SF, O_H_R3_C20_SF, O_H_N2_C20_SF, O_H_S2_C20_SF, O_H_COMBINED_C20_SF) for SF regions only
        - line_maps: tuple of six arrays (HB4861, HA6562, OIII5006, NII6583, SII6716, SII6730) for SF regions only
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
    
    # Apply SF mask to create metallicity maps (only for SF regions)
    O_H_D16_SF = np.where(mask_SF, O_H_D16, np.nan)
    O_H_PG16_SF = np.where(mask_SF, O_H_PG16, np.nan)
    O_H_O3N2_M13_SF = np.where(mask_SF, O_H_O3N2_M13, np.nan)
    O_H_N2_M13_SF = np.where(mask_SF, O_H_N2_M13, np.nan)
    O_H_O3N2_PP04_SF = np.where(mask_SF, O_H_O3N2_PP04, np.nan)
    O_H_N2_PP04_SF = np.where(mask_SF, O_H_N2_PP04, np.nan)
    O_H_O3N2_C20_SF = np.where(mask_SF, O_H_O3N2_C20, np.nan)
    O_H_O3S2_C20_SF = np.where(mask_SF, O_H_O3S2_C20, np.nan)
    O_H_RS32_C20_SF = np.where(mask_SF, O_H_RS32_C20, np.nan)
    O_H_R3_C20_SF = np.where(mask_SF, O_H_R3_C20, np.nan)
    O_H_N2_C20_SF = np.where(mask_SF, O_H_N2_C20, np.nan)
    O_H_S2_C20_SF = np.where(mask_SF, O_H_S2_C20, np.nan)
    O_H_COMBINED_C20_SF = np.where(mask_SF, O_H_COMBINED_C20, np.nan)
    O_H_COMBINED_C20_SF_ERR = np.where(mask_SF, O_H_COMBINED_C20_ERR, np.nan)
    
    # Create SF error maps for all individual C20 methods
    O_H_O3N2_C20_SF_ERR = np.where(mask_SF, O_H_O3N2_C20_ERR, np.nan)
    O_H_O3S2_C20_SF_ERR = np.where(mask_SF, O_H_O3S2_C20_ERR, np.nan)
    O_H_RS32_C20_SF_ERR = np.where(mask_SF, O_H_RS32_C20_ERR, np.nan)
    O_H_R3_C20_SF_ERR = np.where(mask_SF, O_H_R3_C20_ERR, np.nan)
    O_H_N2_C20_SF_ERR = np.where(mask_SF, O_H_N2_C20_ERR, np.nan)
    O_H_S2_C20_SF_ERR = np.where(mask_SF, O_H_S2_C20_ERR, np.nan)

    # Apply SF mask to create line maps in SF regions
    HB4861_FLUX_corr_SF = np.where(mask_SF, HB4861_FLUX_corr, np.nan)
    HA6562_FLUX_corr_SF = np.where(mask_SF, HA6562_FLUX_corr, np.nan)
    OIII5006_FLUX_corr_SF = np.where(mask_SF, OIII5006_FLUX_corr, np.nan)
    NII6583_FLUX_corr_SF = np.where(mask_SF, NII6583_FLUX_corr, np.nan)
    SII6716_FLUX_corr_SF = np.where(mask_SF, SII6716_FLUX_corr, np.nan)
    SII6730_FLUX_corr_SF = np.where(mask_SF, SII6730_FLUX_corr, np.nan)
    
    # Apply SF mask to create line error maps in SF regions  
    HB4861_FLUX_ERR_SF = np.where(mask_SF, HB4861_FLUX_ERR, np.nan)
    HA6562_FLUX_ERR_SF = np.where(mask_SF, HA6562_FLUX_ERR, np.nan)
    OIII5006_FLUX_ERR_SF = np.where(mask_SF, OIII5006_FLUX_ERR, np.nan)
    NII6583_FLUX_ERR_SF = np.where(mask_SF, NII6583_FLUX_ERR, np.nan)
    SII6716_FLUX_ERR_SF = np.where(mask_SF, SII6716_FLUX_ERR, np.nan)
    SII6730_FLUX_ERR_SF = np.where(mask_SF, SII6730_FLUX_ERR, np.nan)

    # Return SFR maps, metallicity maps, metallicity error maps, line maps, and masks
    sfr_maps = (LOG_SFR_surface_density_map_SF, LOG_SFR_surface_density_map_nonSF, 
                LOG_SFR_surface_density_map_unclassified, LOG_SFR_surface_density_map_upper)
    metallicity_maps = (O_H_D16_SF, O_H_PG16_SF, O_H_O3N2_M13_SF, O_H_N2_M13_SF, O_H_O3N2_PP04_SF, O_H_N2_PP04_SF, O_H_O3N2_C20_SF, O_H_O3S2_C20_SF, O_H_RS32_C20_SF, O_H_R3_C20_SF, O_H_N2_C20_SF, O_H_S2_C20_SF, O_H_COMBINED_C20_SF)
    metallicity_error_maps = (O_H_O3N2_C20_SF_ERR, O_H_O3S2_C20_SF_ERR, O_H_RS32_C20_SF_ERR, O_H_R3_C20_SF_ERR, O_H_N2_C20_SF_ERR, O_H_S2_C20_SF_ERR, O_H_COMBINED_C20_SF_ERR)
    line_maps = (HB4861_FLUX_corr_SF, HA6562_FLUX_corr_SF, OIII5006_FLUX_corr_SF, 
                 NII6583_FLUX_corr_SF, SII6716_FLUX_corr_SF, SII6730_FLUX_corr_SF)
    masks = (mask_SF, mask_nonSF, mask_unclassified, mask_upper)
    
    return sfr_maps, metallicity_maps, metallicity_error_maps, line_maps, masks

# Get the SFR surface density maps, metallicity maps, metallicity error maps, line maps, and masks using the default 'both' choice
(LOG_SFR_surface_density_map_SF, LOG_SFR_surface_density_map_nonSF, 
 LOG_SFR_surface_density_map_unclassified, LOG_SFR_surface_density_map_upper), (O_H_D16_SF, O_H_PG16_SF, O_H_O3N2_M13_SF, O_H_N2_M13_SF, O_H_O3N2_PP04_SF, O_H_N2_PP04_SF, O_H_O3N2_C20_SF, O_H_O3S2_C20_SF, O_H_RS32_C20_SF, O_H_R3_C20_SF, O_H_N2_C20_SF, O_H_S2_C20_SF, O_H_COMBINED_C20_SF), (O_H_O3N2_C20_SF_ERR, O_H_O3S2_C20_SF_ERR, O_H_RS32_C20_SF_ERR, O_H_R3_C20_SF_ERR, O_H_N2_C20_SF_ERR, O_H_S2_C20_SF_ERR, O_H_COMBINED_C20_SF_ERR), (HB4861_FLUX_corr_SF, HA6562_FLUX_corr_SF, OIII5006_FLUX_corr_SF, NII6583_FLUX_corr_SF, SII6716_FLUX_corr_SF, SII6730_FLUX_corr_SF), (mask_SF, mask_nonSF, mask_unclassified, mask_upper) = choose_BPT()

# ------------------------------------------------------------------
# 10.  Calculate the total Metallicity in SF regions
# ------------------------------------------------------------------

# nansum the line maps in SF regions
HB4861_FLUX_corr_SF_total = np.nansum(HB4861_FLUX_corr_SF)
HA6562_FLUX_corr_SF_total = np.nansum(HA6562_FLUX_corr_SF)
OIII5006_FLUX_corr_SF_total = np.nansum(OIII5006_FLUX_corr_SF)
NII6583_FLUX_corr_SF_total = np.nansum(NII6583_FLUX_corr_SF)
SII6716_FLUX_corr_SF_total = np.nansum(SII6716_FLUX_corr_SF)
SII6730_FLUX_corr_SF_total = np.nansum(SII6730_FLUX_corr_SF)

# Dopita et al. (2016) metallicity calculation (total)
y_SF_total = np.log10(NII6583_FLUX_corr_SF_total / (SII6716_FLUX_corr_SF_total + SII6730_FLUX_corr_SF_total)) + 0.264*np.log10(NII6583_FLUX_corr_SF_total / HA6562_FLUX_corr_SF_total)
O_H_D16_SF_total = 8.77 + y_SF_total + 0.45*(y_SF_total + 0.3)**5

# Pilyugin & Grebel (2016) metallicity calculation (the S calibration)
# note that here we assume [O III] = 1.33 [O III] 5007, [N II] = 1.34 [N II] 6583, see watts et al. (2024) for details
# PG16 set different coefficients for different branches (logN_2>=-0.6 and logN_2<-0.6)
OIII_scaled_SF_total = 1.33 * OIII5006_FLUX_corr_SF_total  # [O III] = 1.33 * [O III] 5006
NII_scaled_SF_total = 1.34 * NII6583_FLUX_corr_SF_total  # [N II] = 1.34 * [N II] 6583
# Calculate the line ratios needed for PG16
N2_SF_total = NII_scaled_SF_total / HB4861_FLUX_corr_SF_total   # N2 = I([N II]λ6548 + λ6584)/I(Hβ)
S2_SF_total = (SII6716_FLUX_corr_SF_total + SII6730_FLUX_corr_SF_total) / HB4861_FLUX_corr_SF_total  # S2 = I([S II]λ6717 + λ6731)/I(Hβ)
R3_SF_total = OIII_scaled_SF_total / HB4861_FLUX_corr_SF_total  # R3 = I([O III]λ4959 + λ5007)/I(Hβ) (same value as R2 in this case)
# Calculate log values
log_R3_S2_SF_total = np.log10(R3_SF_total/S2_SF_total)
log_N2_SF_total = np.log10(N2_SF_total)
log_S2_SF_total = np.log10(S2_SF_total)
# Determine which branch to use based on log(N2)
# Upper branch: log(N2) >= -0.6
# Lower branch: log(N2) < -0.6
O_H_PG16_SF_total = []
if log_N2_SF_total >= -0.6:
    O_H_PG16_SF_total = (a1_upper + a2_upper * log_R3_S2_SF_total + a3_upper * log_N2_SF_total + 
                      (a4_upper + a5_upper * log_R3_S2_SF_total + a6_upper * log_N2_SF_total) * log_S2_SF_total)
# O3N2-M13 (Marino et al. 2013) metallicity calculation (total)
# Calculate O3N2 ratio for total SF region using M13 calibration
oiii_hb_SF_total = OIII5006_FLUX_corr_SF_total / HB4861_FLUX_corr_SF_total
nii_ha_SF_total = NII6583_FLUX_corr_SF_total / HA6562_FLUX_corr_SF_total
o3n2_ratio_SF_total = np.log10(oiii_hb_SF_total / nii_ha_SF_total)
# Apply O3N2-M13 (Marino et al. 2013) calibration: [O/H] = 8.533 - 0.214 * O3N2
O_H_O3N2_M13_SF_total = 8.533 - 0.214 * o3n2_ratio_SF_total

# N2-M13 (Marino et al. 2013) metallicity calculation (total)
# Calculate N2 ratio for total SF region using M13 calibration
n2_ratio_SF_total = np.log10(NII6583_FLUX_corr_SF_total / HA6562_FLUX_corr_SF_total)
# Apply N2-M13 (Marino et al. 2013) calibration: [O/H] = 8.743 + 0.462 * N2
O_H_N2_M13_SF_total = 8.743 + 0.462 * n2_ratio_SF_total

# O3N2-PP04 (Pettini & Pagel 2004) metallicity calculation (total)
# Calculate O3N2 ratio for total SF region using PP04 calibration
# (reuse previously calculated ratios from O3N2-M13)
# Apply O3N2-PP04 (Pettini & Pagel 2004) calibration: [O/H] = 8.73 - 0.32 * O3N2
O_H_O3N2_PP04_SF_total = 8.73 - 0.32 * o3n2_ratio_SF_total

# N2-PP04 (Pettini & Pagel 2004) metallicity calculation (total)
# Calculate N2 ratio for total SF region using PP04 calibration
# (reuse previously calculated ratio from N2-M13)
# Apply N2-PP04 (Pettini & Pagel 2004) calibration: [O/H] = 9.37 + 2.03*N2 + 1.26*N2^2 + 0.32*N2^3
O_H_N2_PP04_SF_total = (9.37 + 2.03 * n2_ratio_SF_total + 
                       1.26 * n2_ratio_SF_total**2 + 
                       0.32 * n2_ratio_SF_total**3)

# O3N2-C20 (Curti et al. 2020) metallicity calculation (SF total) with error propagation
if (np.isfinite(HB4861_FLUX_corr_SF_total) and np.isfinite(OIII5006_FLUX_corr_SF_total) and
    np.isfinite(NII6583_FLUX_corr_SF_total) and np.isfinite(HA6562_FLUX_corr_SF_total) and
    HB4861_FLUX_corr_SF_total > 0 and OIII5006_FLUX_corr_SF_total > 0 and
    NII6583_FLUX_corr_SF_total > 0 and HA6562_FLUX_corr_SF_total > 0):
    
    # For SF total, use simplified error calculation (similar to spaxel approach)
    # Calculate error propagation for SF total fluxes 
    # For now, use estimated relative errors (can be improved with proper error propagation)
    rel_err_estimate = 0.05  # 5% relative error estimate for total fluxes
    
    oiii_hb_err_sf_total = rel_err_estimate * np.sqrt(2)  # Combined relative error for ratio
    nii_ha_err_sf_total = rel_err_estimate * np.sqrt(2)   # Combined relative error for ratio
    
    # Error for O3N2 = log10(OIII/Hb) - log10(NII/Ha)
    o3n2_ratio_err_sf_total = np.sqrt(oiii_hb_err_sf_total**2 + nii_ha_err_sf_total**2)
    
    # Apply O3N2-C20 calibration
    y_total = o3n2_ratio_SF_total
    y_err_total = o3n2_ratio_err_sf_total
    c0 = 0.281
    c1 = -4.765
    c2 = -2.268
    a = c2
    b = c1
    c = c0 - y_total
    
    # Calculate discriminant for total
    discriminant_total = b**2 - 4*a*c
    
    # Calculate O3N2-C20 total metallicity if discriminant is positive
    if discriminant_total >= 0:
        x_solution1_total = (-b + np.sqrt(discriminant_total)) / (2*a)
        x_solution2_total = (-b - np.sqrt(discriminant_total)) / (2*a)
        
        # Choose the solution that gives reasonable metallicity values
        if (x_solution1_total >= -1.1) and (x_solution1_total <= 1.25):
            x_final_total = x_solution1_total
        else:
            x_final_total = x_solution2_total
        
        # Calculate error in x using derivative approach
        derivative_x_sf_total = np.abs(c1 + 2*c2*x_final_total)
        x_err_sf_total = y_err_total / derivative_x_sf_total
        
        O_H_O3N2_C20_SF_total = x_final_total + 8.69
        O_H_O3N2_C20_SF_total_ERR = np.sqrt(x_err_sf_total**2)  # Can add fitting error here if needed
    else:
        O_H_O3N2_C20_SF_total = np.nan
        O_H_O3N2_C20_SF_total_ERR = np.nan
else:
    O_H_O3N2_C20_SF_total = np.nan
    O_H_O3N2_C20_SF_total_ERR = np.nan

# O3S2-C20 (Curti et al. 2020) SF total metallicity calculation with error propagation
if (np.isfinite(OIII5006_FLUX_corr_SF_total) and np.isfinite(HB4861_FLUX_corr_SF_total) and
    np.isfinite(SII6716_FLUX_corr_SF_total) and np.isfinite(SII6730_FLUX_corr_SF_total) and
    OIII5006_FLUX_corr_SF_total > 0 and HB4861_FLUX_corr_SF_total > 0 and
    SII6716_FLUX_corr_SF_total > 0 and SII6730_FLUX_corr_SF_total > 0):
    
    # Calculate line ratios
    oiii_hb_total = OIII5006_FLUX_corr_SF_total / HB4861_FLUX_corr_SF_total
    sii_total_total = SII6716_FLUX_corr_SF_total + SII6730_FLUX_corr_SF_total
    sii_hb_total = sii_total_total / HB4861_FLUX_corr_SF_total
    o3s2_ratio_total = np.log10(oiii_hb_total / sii_hb_total)
    
    # For SF total, use simplified error calculation (similar to O3N2 approach)
    rel_err_estimate = 0.05  # 5% relative error estimate for SF total fluxes
    
    oiii_hb_err_sf = rel_err_estimate * np.sqrt(2)  # Combined relative error for ratio
    sii_hb_err_sf = rel_err_estimate * np.sqrt(2)   # Combined relative error for ratio
    
    # Error for O3S2 = log10(OIII/Hb) - log10(SII/Hb)
    o3s2_ratio_err_total = np.sqrt(oiii_hb_err_sf**2 + sii_hb_err_sf**2)
    
    # Coefficients from Curti+2020 for O3S2
    c0 = 0.191
    c1 = -4.292
    c2 = -2.538
    c3 = 0.053
    c4 = 0.332
    
    y_total = o3s2_ratio_total
    y_err_total = o3s2_ratio_err_total
    poly_coeffs = [c4, c3, c2, c1, (c0 - y_total)]
    roots = np.roots(poly_coeffs)
    
    # Select the real root that gives reasonable metallicity values
    real_roots = roots[np.isreal(roots)].real
    if len(real_roots) > 0:
        # Choose root that gives reasonable x values
        reasonable_roots = real_roots[(real_roots >= -2) & (real_roots <= 2)]
        if len(reasonable_roots) > 0:
            x_final_total = reasonable_roots[0]
            
            # Calculate error in x using derivative approach 
            # For 4th order polynomial: df/dx = -(c1 + 2*c2*x + 3*c3*x^2 + 4*c4*x^3)
            derivative_x_sf = np.abs(c1 + 2*c2*x_final_total + 3*c3*x_final_total**2 + 4*c4*x_final_total**3)
            x_err_sf = y_err_total / derivative_x_sf
            
            O_H_O3S2_C20_SF_total = x_final_total + 8.69
            O_H_O3S2_C20_SF_total_ERR = np.sqrt(x_err_sf**2)  # Can add fitting error here if needed
        else:
            O_H_O3S2_C20_SF_total = np.nan
            O_H_O3S2_C20_SF_total_ERR = np.nan
    else:
        O_H_O3S2_C20_SF_total = np.nan
        O_H_O3S2_C20_SF_total_ERR = np.nan
else:
    O_H_O3S2_C20_SF_total = np.nan
    O_H_O3S2_C20_SF_total_ERR = np.nan

# RS32-C20 (Curti et al. 2020) total metallicity calculation
if np.isfinite(OIII5006_FLUX_corr_SF_total) and np.isfinite(HB4861_FLUX_corr_SF_total) and \
   np.isfinite(HA6562_FLUX_corr_SF_total) and np.isfinite(SII6716_FLUX_corr_SF_total) and \
   np.isfinite(SII6730_FLUX_corr_SF_total) and \
   OIII5006_FLUX_corr_SF_total > 0 and HB4861_FLUX_corr_SF_total > 0 and \
   HA6562_FLUX_corr_SF_total > 0 and SII6716_FLUX_corr_SF_total > 0 and SII6730_FLUX_corr_SF_total > 0:
    
    # RS32 = log10([OIII]/Hβ + [SII]/Hα)
    oiii_hb_total = OIII5006_FLUX_corr_SF_total / HB4861_FLUX_corr_SF_total
    sii_total = SII6716_FLUX_corr_SF_total + SII6730_FLUX_corr_SF_total
    sii_ha_total = sii_total / HA6562_FLUX_corr_SF_total
    r_lin_total = oiii_hb_total + sii_ha_total
    
    if r_lin_total > 0:
        y_total = np.log10(r_lin_total)
        
        # For SF region, use simplified error estimates where detailed flux errors are not readily available
        # Approximate error as 10% of the flux values (typical for integrated measurements)
        oiii_hb_err_sf = 0.1 * oiii_hb_total
        sii_ha_err_sf = 0.1 * sii_ha_total
        r_lin_err_sf = np.sqrt(oiii_hb_err_sf**2 + sii_ha_err_sf**2)
        y_err_sf = (1/np.log(10)) * (r_lin_err_sf / r_lin_total)
        
        # Coefficients from Curti+2020 for RS32
        c0 = -0.054
        c1 = -2.546
        c2 = -1.970
        c3 = 0.082
        c4 = 0.222
        
        poly_coeffs = [c4, c3, c2, c1, (c0 - y_total)]
        roots = np.roots(poly_coeffs)
        
        # Select the real root that gives reasonable metallicity values
        real_roots = roots[np.isreal(roots)].real
        if len(real_roots) > 0:
            # Sensible metallicity range for Te-anchored scales:
            # 12+log(O/H) ~ 8.0–8.9 ⇒ x = (12+log(O/H))-8.69 ∈ [-0.7, +0.3]
            reasonable_roots = real_roots[(real_roots >= -0.7) & (real_roots <= 0.3)]
            if len(reasonable_roots) > 0:
                # Choose the root that best reproduces y
                y_pred = c0 + c1*reasonable_roots + c2*reasonable_roots**2 + c3*reasonable_roots**3 + c4*reasonable_roots**4
                x_final_total = reasonable_roots[np.argmin(np.abs(y_pred - y_total))]
                O_H_RS32_C20_SF_total = x_final_total + 8.69
                
                # Calculate error in x using derivative approach
                derivative_x = (np.abs(c1 + 2*c2*x_final_total + 3*c3*x_final_total**2 + 4*c4*x_final_total**3))
                x_err_sf = y_err_sf / derivative_x
                O_H_RS32_C20_SF_total_ERR = x_err_sf
            else:
                # Fallback to any real root if none in reasonable range
                y_pred = c0 + c1*real_roots + c2*real_roots**2 + c3*real_roots**3 + c4*real_roots**4
                x_final_total = real_roots[np.argmin(np.abs(y_pred - y_total))]
                O_H_RS32_C20_SF_total = x_final_total + 8.69
                
                # Calculate error in x using derivative approach
                derivative_x = (np.abs(c1 + 2*c2*x_final_total + 3*c3*x_final_total**2 + 4*c4*x_final_total**3))
                x_err_sf = y_err_sf / derivative_x
                O_H_RS32_C20_SF_total_ERR = x_err_sf
        else:
            O_H_RS32_C20_SF_total = np.nan
            O_H_RS32_C20_SF_total_ERR = np.nan
    else:
        O_H_RS32_C20_SF_total = np.nan
        O_H_RS32_C20_SF_total_ERR = np.nan
else:
    O_H_RS32_C20_SF_total = np.nan
    O_H_RS32_C20_SF_total_ERR = np.nan

# R3-C20 (Curti et al. 2020) total metallicity calculation
if np.isfinite(OIII5006_FLUX_corr_SF_total) and np.isfinite(HB4861_FLUX_corr_SF_total) and \
   OIII5006_FLUX_corr_SF_total > 0 and HB4861_FLUX_corr_SF_total > 0:
    
    # R3 = log10([OIII]/Hβ)
    r_lin_total = OIII5006_FLUX_corr_SF_total / HB4861_FLUX_corr_SF_total
    
    if r_lin_total > 0:
        y_total = np.log10(r_lin_total)
        
        # For SF region, use simplified error estimates where detailed flux errors are not readily available
        # Approximate error as 10% of the flux values (typical for integrated measurements)
        r3_error_sf = 0.1  # simplified error estimate for R3 ratio
        
        # Coefficients from Curti+2020 for R3
        c0 = -0.277
        c1 = -3.549
        c2 = -3.593
        c3 = -0.981
        
        poly_coeffs = [c3, c2, c1, (c0 - y_total)]
        roots = np.roots(poly_coeffs)
        
        # Select the real root that gives reasonable metallicity values
        real_roots = roots[np.isreal(roots)].real
        if len(real_roots) > 0:
            # Sensible metallicity range for Te-anchored scales:
            # 12+log(O/H) ~ 8.0–8.9 ⇒ x = (12+log(O/H))-8.69 ∈ [-0.7, +0.3]
            reasonable_roots = real_roots[(real_roots >= -0.7) & (real_roots <= 0.3)]
            if len(reasonable_roots) > 0:
                # Choose the root that best reproduces y
                y_pred = c0 + c1*reasonable_roots + c2*reasonable_roots**2 + c3*reasonable_roots**3
                x_final_total = reasonable_roots[np.argmin(np.abs(y_pred - y_total))]
                O_H_R3_C20_SF_total = x_final_total + 8.69
                
                # Error propagation: derivative of polynomial with respect to y
                derivative_y = np.abs(c1 + 2*c2*x_final_total + 3*c3*x_final_total**2)
                if derivative_y > 0:
                    derivative_x = 1.0 / derivative_y
                    obs_error = derivative_x * r3_error_sf
                    O_H_R3_C20_SF_total_ERR = obs_error
                else:
                    O_H_R3_C20_SF_total_ERR = np.nan
            else:
                # Fallback to any real root if none in reasonable range
                y_pred = c0 + c1*real_roots + c2*real_roots**2 + c3*real_roots**3
                x_final_total = real_roots[np.argmin(np.abs(y_pred - y_total))]
                O_H_R3_C20_SF_total = x_final_total + 8.69
                
                # Error propagation: derivative of polynomial with respect to y
                derivative_y = np.abs(c1 + 2*c2*x_final_total + 3*c3*x_final_total**2)
                if derivative_y > 0:
                    derivative_x = 1.0 / derivative_y
                    obs_error = derivative_x * r3_error_sf
                    O_H_R3_C20_SF_total_ERR = obs_error
                else:
                    O_H_R3_C20_SF_total_ERR = np.nan
        else:
            O_H_R3_C20_SF_total = np.nan
            O_H_R3_C20_SF_total_ERR = np.nan
    else:
        O_H_R3_C20_SF_total = np.nan
        O_H_R3_C20_SF_total_ERR = np.nan
else:
    O_H_R3_C20_SF_total = np.nan
    O_H_R3_C20_SF_total_ERR = np.nan

# N2-C20 (Curti et al. 2020) total metallicity calculation
if np.isfinite(NII6583_FLUX_corr_SF_total) and np.isfinite(HA6562_FLUX_corr_SF_total) and \
   NII6583_FLUX_corr_SF_total > 0 and HA6562_FLUX_corr_SF_total > 0:
    
    # N2 = log10([NII]/Hα)
    n2_lin_total = NII6583_FLUX_corr_SF_total / HA6562_FLUX_corr_SF_total
    
    if n2_lin_total > 0:
        y_total = np.log10(n2_lin_total)
        
        # For SF region, use simplified error estimates where detailed flux errors are not readily available
        # Approximate error as 10% of the flux values (typical for integrated measurements)
        n2_error_sf = 0.1  # simplified error estimate for N2 ratio
        
        # Coefficients from Curti+2020 for N2
        c0 = -0.489
        c1 = 1.513
        c2 = -2.554
        c3 = -5.293
        c4 = -2.867
        
        poly_coeffs = [c4, c3, c2, c1, (c0 - y_total)]
        roots = np.roots(poly_coeffs)
        
        # tolerances for "almost-real" roots and range comparison
        REAL_ATOL = 1e-8
        RANGE_EPS = 0.0  # exact bounds as requested
        
        # treat tiny-imag roots as real
        realish = roots[np.abs(roots.imag) <= REAL_ATOL].real
        if realish.size == 0:
            # fall back to the least-imag root if imag part is tiny-ish
            k = np.argmin(np.abs(roots.imag))
            if np.abs(roots[k].imag) <= 1e-6:
                realish = np.array([roots[k].real])
        
        if realish.size > 0:
            # STRICT selection inside [-0.7, 0.3]
            in_rng = realish[(realish >= -0.7 - RANGE_EPS) & (realish <= 0.3 + RANGE_EPS)]
            if in_rng.size > 0:
                # If multiple, pick the smallest one
                x_final_total = np.min(in_rng)
                O_H_N2_C20_SF_total = x_final_total + 8.69
                
                # Error propagation: derivative of 4th-order polynomial with respect to y
                derivative_y = np.abs(c1 + 2*c2*x_final_total + 3*c3*x_final_total**2 + 4*c4*x_final_total**3)
                if derivative_y > 0:
                    derivative_x = 1.0 / derivative_y
                    obs_error = derivative_x * n2_error_sf
                    O_H_N2_C20_SF_total_ERR = obs_error
                else:
                    O_H_N2_C20_SF_total_ERR = np.nan
            else:
                # No in-range root → discard this spaxel
                O_H_N2_C20_SF_total = np.nan
                O_H_N2_C20_SF_total_ERR = np.nan
        else:
            O_H_N2_C20_SF_total = np.nan
            O_H_N2_C20_SF_total_ERR = np.nan
    else:
        O_H_N2_C20_SF_total = np.nan
        O_H_N2_C20_SF_total_ERR = np.nan
else:
    O_H_N2_C20_SF_total = np.nan
    O_H_N2_C20_SF_total_ERR = np.nan

# S2-C20 (Curti et al. 2020) total metallicity calculation
if np.isfinite(SII6716_FLUX_corr_SF_total) and np.isfinite(SII6730_FLUX_corr_SF_total) and \
   np.isfinite(HA6562_FLUX_corr_SF_total) and \
   SII6716_FLUX_corr_SF_total > 0 and SII6730_FLUX_corr_SF_total > 0 and HA6562_FLUX_corr_SF_total > 0:
    
    # S2 = log10(([SII]6716+6730)/Hα)
    s2_lin_total = (SII6716_FLUX_corr_SF_total + SII6730_FLUX_corr_SF_total) / HA6562_FLUX_corr_SF_total
    
    if s2_lin_total > 0:
        y_total = np.log10(s2_lin_total)
        
        # For SF region, use simplified error estimates where detailed flux errors are not readily available
        # Approximate error as 10% of the flux values (typical for integrated measurements)
        s2_error_sf = 0.1  # simplified error estimate for S2 ratio
        
        # Coefficients from Curti+2020 for S2
        c0 = -0.442
        c1 = -0.360
        c2 = -6.271
        c3 = -8.339
        c4 = -3.559
        
        poly_coeffs = [c4, c3, c2, c1, (c0 - y_total)]
        roots = np.roots(poly_coeffs)
        
        REAL_ATOL = 1e-8  # accept roots with tiny imaginary part as real
        
        # Treat tiny-imag roots as real
        realish = roots[np.abs(roots.imag) <= REAL_ATOL].real
        if realish.size == 0:
            # fallback: least-imag root if imag part is still tiny-ish
            k = np.argmin(np.abs(roots.imag))
            if np.abs(roots[k].imag) <= 1e-6:
                realish = np.array([roots[k].real])
        
        if realish.size > 0:
            # STRICT in-range selection: x ∈ [-0.7, 0.3]
            in_range = realish[(realish >= -0.7) & (realish <= 0.3)]
            if in_range.size > 0:
                # If multiple, pick the smallest
                x_final_total = np.min(in_range)
                O_H_S2_C20_SF_total = x_final_total + 8.69
                
                # Error propagation: derivative of 4th-order polynomial with respect to y
                derivative_y = np.abs(c1 + 2*c2*x_final_total + 3*c3*x_final_total**2 + 4*c4*x_final_total**3)
                if derivative_y > 0:
                    derivative_x = 1.0 / derivative_y
                    obs_error = derivative_x * s2_error_sf
                    O_H_S2_C20_SF_total_ERR = obs_error
                else:
                    O_H_S2_C20_SF_total_ERR = np.nan
            else:
                # discard if no valid root in range
                O_H_S2_C20_SF_total = np.nan
                O_H_S2_C20_SF_total_ERR = np.nan
        else:
            O_H_S2_C20_SF_total = np.nan
            O_H_S2_C20_SF_total_ERR = np.nan
    else:
        O_H_S2_C20_SF_total = np.nan
        O_H_S2_C20_SF_total_ERR = np.nan
else:
    O_H_S2_C20_SF_total = np.nan
    O_H_S2_C20_SF_total_ERR = np.nan

# Total Combined C20 metallicity in SF region
# Calculate metallicity using integrated line fluxes for each C20 method and select minimum error
if np.any(np.isfinite(O_H_COMBINED_C20_SF)):
    # Get integrated line fluxes from SF regions (these are calculated above)
    
    # Check if we have valid integrated fluxes for C20 calculations
    has_valid_fluxes = (
        np.isfinite(HB4861_FLUX_corr_SF_total) and np.isfinite(HA6562_FLUX_corr_SF_total) and
        np.isfinite(OIII5006_FLUX_corr_SF_total) and np.isfinite(NII6583_FLUX_corr_SF_total) and
        np.isfinite(SII6716_FLUX_corr_SF_total) and np.isfinite(SII6730_FLUX_corr_SF_total) and
        HB4861_FLUX_corr_SF_total > 0 and HA6562_FLUX_corr_SF_total > 0 and
        OIII5006_FLUX_corr_SF_total > 0 and NII6583_FLUX_corr_SF_total > 0 and
        SII6716_FLUX_corr_SF_total > 0 and SII6730_FLUX_corr_SF_total > 0
    )
    
    if has_valid_fluxes:
        # Calculate metallicity and error for each C20 method using integrated fluxes
        c20_methods = []
        c20_errors = []
        c20_names = ['O3N2-C20', 'O3S2-C20', 'RS32-C20', 'R3-C20', 'N2-C20', 'S2-C20']
        
        # Method 0: O3N2-C20 (use calculated total values)
        if np.isfinite(O_H_O3N2_C20_SF_total) and np.isfinite(O_H_O3N2_C20_SF_total_ERR):
            c20_methods.append(O_H_O3N2_C20_SF_total)
            c20_errors.append(O_H_O3N2_C20_SF_total_ERR)
        else:
            c20_methods.append(np.nan)
            c20_errors.append(np.inf)
        
        # Method 1: O3S2-C20 (use calculated total values)
        if np.isfinite(O_H_O3S2_C20_SF_total) and np.isfinite(O_H_O3S2_C20_SF_total_ERR):
            c20_methods.append(O_H_O3S2_C20_SF_total)
            c20_errors.append(O_H_O3S2_C20_SF_total_ERR)
        else:
            c20_methods.append(np.nan)
            c20_errors.append(np.inf)
        
        # Method 2: RS32-C20 (use calculated total values)
        if np.isfinite(O_H_RS32_C20_SF_total) and np.isfinite(O_H_RS32_C20_SF_total_ERR):
            c20_methods.append(O_H_RS32_C20_SF_total)
            c20_errors.append(O_H_RS32_C20_SF_total_ERR)
        else:
            c20_methods.append(np.nan)
            c20_errors.append(np.inf)
        
        # Method 3: R3-C20 (use calculated total values)
        if np.isfinite(O_H_R3_C20_SF_total) and np.isfinite(O_H_R3_C20_SF_total_ERR):
            c20_methods.append(O_H_R3_C20_SF_total)
            c20_errors.append(O_H_R3_C20_SF_total_ERR)
        else:
            c20_methods.append(np.nan)
            c20_errors.append(np.inf)
        
        # Method 4: N2-C20 (use calculated total values)
        if np.isfinite(O_H_N2_C20_SF_total) and np.isfinite(O_H_N2_C20_SF_total_ERR):
            c20_methods.append(O_H_N2_C20_SF_total)
            c20_errors.append(O_H_N2_C20_SF_total_ERR)
        else:
            c20_methods.append(np.nan)
            c20_errors.append(np.inf)
        
        # Method 5: S2-C20 (use calculated total values)
        if np.isfinite(O_H_S2_C20_SF_total) and np.isfinite(O_H_S2_C20_SF_total_ERR):
            c20_methods.append(O_H_S2_C20_SF_total)
            c20_errors.append(O_H_S2_C20_SF_total_ERR)
        else:
            c20_methods.append(np.nan)
            c20_errors.append(np.inf)
        
        # Find method with minimum error among valid results
        valid_indices = [i for i, (oh, err) in enumerate(zip(c20_methods, c20_errors)) if np.isfinite(oh) and np.isfinite(err)]
        
        if valid_indices:
            min_error_idx = min(valid_indices, key=lambda i: c20_errors[i])
            O_H_COMBINED_C20_SF_total = c20_methods[min_error_idx]
            best_method_name = c20_names[min_error_idx]
            print(f"Combined C20 total: Selected {best_method_name} with error {c20_errors[min_error_idx]:.3f}")
            print(f"Available methods: {[(c20_names[i], c20_methods[i], c20_errors[i]) for i in valid_indices]}")
        else:
            O_H_COMBINED_C20_SF_total = np.nan
    else:
        O_H_COMBINED_C20_SF_total = np.nan
else:
    O_H_COMBINED_C20_SF_total = np.nan

# ------------------------------------------------------------------
# 11.  Calculate the total Metallicity in total available regions
# ------------------------------------------------------------------

# nansum the line maps in total available regions
HB4861_FLUX_total = np.nansum(HB4861_FLUX)
HA6562_FLUX_total = np.nansum(HA6562_FLUX)
OIII5006_FLUX_total = np.nansum(OIII5006_FLUX)
NII6583_FLUX_total = np.nansum(NII6583_FLUX)
SII6716_FLUX_total = np.nansum(SII6716_FLUX)
SII6730_FLUX_total = np.nansum(SII6730_FLUX)

# Calculate total flux errors by error propagation (sqrt of sum of squares)
HB4861_FLUX_ERR_total = np.sqrt(np.nansum(HB4861_FLUX_ERR**2))
HA6562_FLUX_ERR_total = np.sqrt(np.nansum(HA6562_FLUX_ERR**2))
OIII5006_FLUX_ERR_total = np.sqrt(np.nansum(OIII5006_FLUX_ERR**2))
NII6583_FLUX_ERR_total = np.sqrt(np.nansum(NII6583_FLUX_ERR**2))
SII6716_FLUX_ERR_total = np.sqrt(np.nansum(SII6716_FLUX_ERR**2))
SII6730_FLUX_ERR_total = np.sqrt(np.nansum(SII6730_FLUX_ERR**2))

# Calculate Balmer Decrement for total flux
BD_total = HA6562_FLUX_total / HB4861_FLUX_total
if BD_total < 2.86:
    BD_total = 2.86

# Calculate E(B-V) for total flux
E_BV_BD_total = convert_bd_to_ebv(BD_total, k_HB4861, k_HA6562, R_int)

# Correct total fluxes with the uniform E(B-V)
HB4861_FLUX_corr_total = correct_flux_with_ebv(HB4861_FLUX_total, E_BV_BD_total, k_HB4861)
HA6562_FLUX_corr_total = correct_flux_with_ebv(HA6562_FLUX_total, E_BV_BD_total, k_HA6562)
OIII5006_FLUX_corr_total = correct_flux_with_ebv(OIII5006_FLUX_total, E_BV_BD_total, k_OIII5006)
NII6583_FLUX_corr_total = correct_flux_with_ebv(NII6583_FLUX_total, E_BV_BD_total, k_NII6583)
SII6716_FLUX_corr_total = correct_flux_with_ebv(SII6716_FLUX_total, E_BV_BD_total, k_SII6716)
SII6730_FLUX_corr_total = correct_flux_with_ebv(SII6730_FLUX_total, E_BV_BD_total, k_SII6730)

# Dopita et al. (2016) metallicity calculation (total)
y_total = np.log10(NII6583_FLUX_corr_total / (SII6716_FLUX_corr_total + SII6730_FLUX_corr_total)) + 0.264*np.log10(NII6583_FLUX_corr_total / HA6562_FLUX_corr_total)
O_H_D16_total = 8.77 + y_total + 0.45*(y_total + 0.3)**5

# Pilyugin & Grebel (2016) metallicity calculation (the S calibration)
OIII_scaled_total = 1.33 * OIII5006_FLUX_corr_total
NII_scaled_total = 1.34 * NII6583_FLUX_corr_total
N2_total = NII_scaled_total / HB4861_FLUX_corr_total
S2_total = (SII6716_FLUX_corr_total + SII6730_FLUX_corr_total) / HB4861_FLUX_corr_total
R3_total = OIII_scaled_total / HB4861_FLUX_corr_total
log_R3_S2_total = np.log10(R3_total/S2_total)
log_N2_total = np.log10(N2_total)
log_S2_total = np.log10(S2_total)
O_H_PG16_total = []
if log_N2_total >= -0.6:
    O_H_PG16_total = (a1_upper + a2_upper * log_R3_S2_total + a3_upper * log_N2_total + 
                      (a4_upper + a5_upper * log_R3_S2_total + a6_upper * log_N2_total) * log_S2_total)
else:
    O_H_PG16_total = (a1_lower + a2_lower * log_R3_S2_total + a3_lower * log_N2_total +
                        (a4_lower + a5_lower * log_R3_S2_total + a6_lower * log_N2_total) * log_S2_total)

# O3N2-M13 (Marino et al. 2013) metallicity calculation (total)
# Calculate O3N2 ratio for total SF region using M13 calibration
oiii_hb_total = OIII5006_FLUX_corr_total / HB4861_FLUX_corr_total
nii_ha_total = NII6583_FLUX_corr_total / HA6562_FLUX_corr_total
o3n2_ratio_total = np.log10(oiii_hb_total / nii_ha_total)
# Apply O3N2-M13 (Marino et al. 2013) calibration: [O/H] = 8.533 - 0.214 * O3N2
O_H_O3N2_M13_total = 8.533 - 0.214 * o3n2_ratio_total

# N2-M13 (Marino et al. 2013) metallicity calculation (total)
# Calculate N2 ratio for total SF region using M13 calibration
n2_ratio_total = np.log10(NII6583_FLUX_corr_total / HA6562_FLUX_corr_total)
# Apply N2-M13 (Marino et al. 2013) calibration: [O/H] = 8.743 + 0.462 * N2
O_H_N2_M13_total = 8.743 + 0.462 * n2_ratio_total

# O3N2-PP04 (Pettini & Pagel 2004) metallicity calculation (total)
O_H_O3N2_PP04_total = 8.73 - 0.32 * o3n2_ratio_total

# N2-PP04 (Pettini & Pagel 2004) metallicity calculation (total)
O_H_N2_PP04_total = (9.37 + 2.03 * n2_ratio_total + 
                       1.26 * n2_ratio_total**2 + 
                       0.32 * n2_ratio_total**3)

# O3N2-C20 (Curti et al. 2020) metallicity calculation (total) with error propagation
if (np.isfinite(HB4861_FLUX_corr_total) and np.isfinite(OIII5006_FLUX_corr_total) and
    np.isfinite(NII6583_FLUX_corr_total) and np.isfinite(HA6562_FLUX_corr_total) and
    HB4861_FLUX_corr_total > 0 and OIII5006_FLUX_corr_total > 0 and
    NII6583_FLUX_corr_total > 0 and HA6562_FLUX_corr_total > 0 and
    np.isfinite(HB4861_FLUX_ERR_total) and np.isfinite(OIII5006_FLUX_ERR_total) and
    np.isfinite(NII6583_FLUX_ERR_total) and np.isfinite(HA6562_FLUX_ERR_total)):
    
    # Calculate error propagation for total fluxes
    oiii_hb_err_total = bpt_error_propagation(OIII5006_FLUX_corr_total, HB4861_FLUX_corr_total, 
                                             OIII5006_FLUX_ERR_total, HB4861_FLUX_ERR_total)
    nii_ha_err_total = bpt_error_propagation(NII6583_FLUX_corr_total, HA6562_FLUX_corr_total, 
                                            NII6583_FLUX_ERR_total, HA6562_FLUX_ERR_total)
    
    # Error for O3N2 = log10(OIII/Hb) - log10(NII/Ha)
    o3n2_ratio_err_total = np.sqrt(oiii_hb_err_total**2 + nii_ha_err_total**2)
    
    # Apply O3N2-C20 calibration
    y_total_c20 = o3n2_ratio_total
    y_err_total_c20 = o3n2_ratio_err_total
    c0 = 0.281
    c1 = -4.765
    c2 = -2.268
    a = c2
    b = c1
    c = c0 - y_total_c20
    discriminant_total_c20 = b**2 - 4*a*c
    
    if discriminant_total_c20 >= 0:
        x_solution1_total = (-b + np.sqrt(discriminant_total_c20)) / (2*a)
        x_solution2_total = (-b - np.sqrt(discriminant_total_c20)) / (2*a)
        if (x_solution1_total >= -1.1) and (x_solution1_total <= 1.25):
            x_final_total = x_solution1_total
        else:
            x_final_total = x_solution2_total
        
        # Calculate error in x using derivative approach
        derivative_x_total = np.abs(c1 + 2*c2*x_final_total)
        x_err_total = y_err_total_c20 / derivative_x_total
        
        O_H_O3N2_C20_total = x_final_total + 8.69
        O_H_O3N2_C20_total_ERR = np.sqrt(x_err_total**2)  # Can add fitting error here if needed
    else:
        O_H_O3N2_C20_total = np.nan
        O_H_O3N2_C20_total_ERR = np.nan
else:
    O_H_O3N2_C20_total = np.nan
    O_H_O3N2_C20_total_ERR = np.nan

# O3S2-C20 (Curti et al. 2020) total metallicity calculation with error propagation
if (np.isfinite(OIII5006_FLUX_corr_total) and np.isfinite(HB4861_FLUX_corr_total) and
    np.isfinite(SII6716_FLUX_corr_total) and np.isfinite(SII6730_FLUX_corr_total) and
    OIII5006_FLUX_corr_total > 0 and HB4861_FLUX_corr_total > 0 and
    SII6716_FLUX_corr_total > 0 and SII6730_FLUX_corr_total > 0 and
    np.isfinite(OIII5006_FLUX_ERR_total) and np.isfinite(HB4861_FLUX_ERR_total) and
    np.isfinite(SII6716_FLUX_ERR_total) and np.isfinite(SII6730_FLUX_ERR_total)):
    
    # Calculate line ratios
    oiii_hb_total_c20 = OIII5006_FLUX_corr_total / HB4861_FLUX_corr_total
    sii_total_total_c20 = SII6716_FLUX_corr_total + SII6730_FLUX_corr_total
    sii_hb_total_c20 = sii_total_total_c20 / HB4861_FLUX_corr_total
    o3s2_ratio_total_c20 = np.log10(oiii_hb_total_c20 / sii_hb_total_c20)
    
    # Calculate error propagation for total fluxes
    sii_total_err_total = np.sqrt(SII6716_FLUX_ERR_total**2 + SII6730_FLUX_ERR_total**2)
    oiii_hb_err_total = bpt_error_propagation(OIII5006_FLUX_corr_total, HB4861_FLUX_corr_total,
                                             OIII5006_FLUX_ERR_total, HB4861_FLUX_ERR_total)
    sii_hb_err_total = bpt_error_propagation(sii_total_total_c20, HB4861_FLUX_corr_total,
                                            sii_total_err_total, HB4861_FLUX_ERR_total)
    
    # Error for O3S2 = log10(OIII/Hb) - log10(SII/Hb)
    o3s2_ratio_err_total = np.sqrt(oiii_hb_err_total**2 + sii_hb_err_total**2)
    
    # Apply O3S2-C20 calibration
    c0 = 0.191; c1 = -4.292; c2 = -2.538; c3 = 0.053; c4 = 0.332
    y_total_c20 = o3s2_ratio_total_c20
    y_err_total_c20 = o3s2_ratio_err_total
    poly_coeffs = [c4, c3, c2, c1, (c0 - y_total_c20)]
    roots = np.roots(poly_coeffs)
    real_roots = roots[np.isreal(roots)].real
    
    if len(real_roots) > 0:
        reasonable_roots = real_roots[(real_roots >= -2) & (real_roots <= 2)]
        if len(reasonable_roots) > 0:
            x_final_total = reasonable_roots[0]
            
            # Calculate error in x using derivative approach 
            # For 4th order polynomial: df/dx = -(c1 + 2*c2*x + 3*c3*x^2 + 4*c4*x^3)
            derivative_x_total = np.abs(c1 + 2*c2*x_final_total + 3*c3*x_final_total**2 + 4*c4*x_final_total**3)
            x_err_total = y_err_total_c20 / derivative_x_total
            
            O_H_O3S2_C20_total = x_final_total + 8.69
            O_H_O3S2_C20_total_ERR = np.sqrt(x_err_total**2)  # Can add fitting error here if needed
        else:
            O_H_O3S2_C20_total = np.nan
            O_H_O3S2_C20_total_ERR = np.nan
    else:
        O_H_O3S2_C20_total = np.nan
        O_H_O3S2_C20_total_ERR = np.nan
else:
    O_H_O3S2_C20_total = np.nan
    O_H_O3S2_C20_total_ERR = np.nan

# RS32-C20 (Curti et al. 2020) total metallicity calculation
if np.isfinite(OIII5006_FLUX_corr_total) and np.isfinite(HB4861_FLUX_corr_total) and \
   np.isfinite(HA6562_FLUX_corr_total) and np.isfinite(SII6716_FLUX_corr_total) and \
   np.isfinite(SII6730_FLUX_corr_total) and \
   OIII5006_FLUX_corr_total > 0 and HB4861_FLUX_corr_total > 0 and \
   HA6562_FLUX_corr_total > 0 and SII6716_FLUX_corr_total > 0 and SII6730_FLUX_corr_total > 0:
    
    # Calculate total flux errors
    HB4861_FLUX_ERR_total = np.sqrt(np.nansum(HB4861_FLUX_ERR[~np.isnan(HB4861_FLUX_corr)]**2))
    HA6562_FLUX_ERR_total = np.sqrt(np.nansum(HA6562_FLUX_ERR[~np.isnan(HA6562_FLUX_corr)]**2))
    OIII5006_FLUX_ERR_total = np.sqrt(np.nansum(OIII5006_FLUX_ERR[~np.isnan(OIII5006_FLUX_corr)]**2))
    SII6716_FLUX_ERR_total = np.sqrt(np.nansum(SII6716_FLUX_ERR[~np.isnan(SII6716_FLUX_corr)]**2))
    SII6730_FLUX_ERR_total = np.sqrt(np.nansum(SII6730_FLUX_ERR[~np.isnan(SII6730_FLUX_corr)]**2))

    oiii_hb_total_c20 = OIII5006_FLUX_corr_total / HB4861_FLUX_corr_total
    sii_total = SII6716_FLUX_corr_total + SII6730_FLUX_corr_total
    sii_total_err = np.sqrt(SII6716_FLUX_ERR_total**2 + SII6730_FLUX_ERR_total**2)
    sii_ha_total_c20 = sii_total / HA6562_FLUX_corr_total
    r_lin_total_c20 = oiii_hb_total_c20 + sii_ha_total_c20
    
    if r_lin_total_c20 > 0:
        y_total_c20 = np.log10(r_lin_total_c20)
        
        # Calculate errors for the line ratios using error propagation
        oiii_hb_err_total = bpt_error_propagation(OIII5006_FLUX_corr_total, HB4861_FLUX_corr_total, 
                                                  OIII5006_FLUX_ERR_total, HB4861_FLUX_ERR_total)
        sii_ha_err_total = bpt_error_propagation(sii_total, HA6562_FLUX_corr_total, 
                                                 sii_total_err, HA6562_FLUX_ERR_total)
        
        # Error for RS32 = log10(OIII/Hb + SII/Ha)
        r_lin_err_total = np.sqrt(oiii_hb_err_total**2 + sii_ha_err_total**2)
        y_err_total = (1/np.log(10)) * (r_lin_err_total / r_lin_total_c20)
        
        c0 = -0.054; c1 = -2.546; c2 = -1.970; c3 = 0.082; c4 = 0.222
        poly_coeffs = [c4, c3, c2, c1, (c0 - y_total_c20)]
        roots = np.roots(poly_coeffs)
        real_roots = roots[np.isreal(roots)].real
        if len(real_roots) > 0:
            reasonable_roots = real_roots[(real_roots >= -0.7) & (real_roots <= 0.3)]
            if len(reasonable_roots) > 0:
                y_pred = c0 + c1*reasonable_roots + c2*reasonable_roots**2 + c3*reasonable_roots**3 + c4*reasonable_roots**4
                x_final_total = reasonable_roots[np.argmin(np.abs(y_pred - y_total_c20))]
                O_H_RS32_C20_total = x_final_total + 8.69
                
                # Calculate error in x using derivative approach
                derivative_x = (np.abs(c1 + 2*c2*x_final_total + 3*c3*x_final_total**2 + 4*c4*x_final_total**3))
                x_err_total = y_err_total / derivative_x
                O_H_RS32_C20_total_ERR = x_err_total
            else:
                y_pred = c0 + c1*real_roots + c2*real_roots**2 + c3*real_roots**3 + c4*real_roots**4
                x_final_total = real_roots[np.argmin(np.abs(y_pred - y_total_c20))]
                O_H_RS32_C20_total = x_final_total + 8.69
                
                # Calculate error in x using derivative approach
                derivative_x = (np.abs(c1 + 2*c2*x_final_total + 3*c3*x_final_total**2 + 4*c4*x_final_total**3))
                x_err_total = y_err_total / derivative_x
                O_H_RS32_C20_total_ERR = x_err_total
        else:
            O_H_RS32_C20_total = np.nan
            O_H_RS32_C20_total_ERR = np.nan
    else:
        O_H_RS32_C20_total = np.nan
        O_H_RS32_C20_total_ERR = np.nan
else:
    O_H_RS32_C20_total = np.nan
    O_H_RS32_C20_total_ERR = np.nan

# R3-C20 (Curti et al. 2020) total metallicity calculation
if np.isfinite(OIII5006_FLUX_corr_total) and np.isfinite(HB4861_FLUX_corr_total) and \
   OIII5006_FLUX_corr_total > 0 and HB4861_FLUX_corr_total > 0:
    
    # Calculate total flux errors
    HB4861_FLUX_ERR_total = np.sqrt(np.nansum(HB4861_FLUX_ERR[~np.isnan(HB4861_FLUX_corr)]**2))
    OIII5006_FLUX_ERR_total = np.sqrt(np.nansum(OIII5006_FLUX_ERR[~np.isnan(OIII5006_FLUX_corr)]**2))
    
    r_lin_total_c20 = OIII5006_FLUX_corr_total / HB4861_FLUX_corr_total
    if r_lin_total_c20 > 0:
        y_total_c20 = np.log10(r_lin_total_c20)
        
        # Calculate error in R3 using BPT error propagation
        r3_error_total = bpt_error_propagation(OIII5006_FLUX_corr_total, HB4861_FLUX_corr_total,
                                               OIII5006_FLUX_ERR_total, HB4861_FLUX_ERR_total)
        
        c0 = -0.277; c1 = -3.549; c2 = -3.593; c3 = -0.981
        poly_coeffs = [c3, c2, c1, (c0 - y_total_c20)]
        roots = np.roots(poly_coeffs)
        real_roots = roots[np.isreal(roots)].real
        if len(real_roots) > 0:
            reasonable_roots = real_roots[(real_roots >= -0.7) & (real_roots <= 0.3)]
            if len(reasonable_roots) > 0:
                y_pred = c0 + c1*reasonable_roots + c2*reasonable_roots**2 + c3*reasonable_roots**3
                x_final_total = reasonable_roots[np.argmin(np.abs(y_pred - y_total_c20))]
                O_H_R3_C20_total = x_final_total + 8.69
                
                # Error propagation: derivative of polynomial with respect to y
                derivative_y = np.abs(c1 + 2*c2*x_final_total + 3*c3*x_final_total**2)
                if derivative_y > 0:
                    derivative_x = 1.0 / derivative_y
                    obs_error = derivative_x * r3_error_total
                    O_H_R3_C20_total_ERR = obs_error
                else:
                    O_H_R3_C20_total_ERR = np.nan
            else:
                y_pred = c0 + c1*real_roots + c2*real_roots**2 + c3*real_roots**3
                x_final_total = real_roots[np.argmin(np.abs(y_pred - y_total_c20))]
                O_H_R3_C20_total = x_final_total + 8.69
                
                # Error propagation: derivative of polynomial with respect to y
                derivative_y = np.abs(c1 + 2*c2*x_final_total + 3*c3*x_final_total**2)
                if derivative_y > 0:
                    derivative_x = 1.0 / derivative_y
                    obs_error = derivative_x * r3_error_total
                    O_H_R3_C20_total_ERR = obs_error
                else:
                    O_H_R3_C20_total_ERR = np.nan
        else:
            O_H_R3_C20_total = np.nan
            O_H_R3_C20_total_ERR = np.nan
    else:
        O_H_R3_C20_total = np.nan
        O_H_R3_C20_total_ERR = np.nan
else:
    O_H_R3_C20_total = np.nan
    O_H_R3_C20_total_ERR = np.nan

# N2-C20 (Curti et al. 2020) total metallicity calculation
if np.isfinite(NII6583_FLUX_corr_total) and np.isfinite(HA6562_FLUX_corr_total) and \
   NII6583_FLUX_corr_total > 0 and HA6562_FLUX_corr_total > 0:
    
    # Calculate total flux errors
    HA6562_FLUX_ERR_total = np.sqrt(np.nansum(HA6562_FLUX_ERR[~np.isnan(HA6562_FLUX_corr)]**2))
    NII6583_FLUX_ERR_total = np.sqrt(np.nansum(NII6583_FLUX_ERR[~np.isnan(NII6583_FLUX_corr)]**2))
    
    n2_lin_total_c20 = NII6583_FLUX_corr_total / HA6562_FLUX_corr_total
    if n2_lin_total_c20 > 0:
        y_total_c20 = np.log10(n2_lin_total_c20)
        
        # Calculate error in N2 using BPT error propagation
        n2_error_total = bpt_error_propagation(NII6583_FLUX_corr_total, HA6562_FLUX_corr_total,
                                               NII6583_FLUX_ERR_total, HA6562_FLUX_ERR_total)
        
        c0 = -0.489; c1 = 1.513; c2 = -2.554; c3 = -5.293; c4 = -2.867
        poly_coeffs = [c4, c3, c2, c1, (c0 - y_total_c20)]
        roots = np.roots(poly_coeffs)
        realish = roots[np.abs(roots.imag) <= 1e-8].real
        if realish.size > 0:
            in_rng = realish[(realish >= -0.7) & (realish <= 0.3)]
            if in_rng.size > 0:
                x_final_total = np.min(in_rng)
                O_H_N2_C20_total = x_final_total + 8.69
                
                # Error propagation: derivative of 4th-order polynomial with respect to y
                derivative_y = np.abs(c1 + 2*c2*x_final_total + 3*c3*x_final_total**2 + 4*c4*x_final_total**3)
                if derivative_y > 0:
                    derivative_x = 1.0 / derivative_y
                    obs_error = derivative_x * n2_error_total
                    O_H_N2_C20_total_ERR = obs_error
                else:
                    O_H_N2_C20_total_ERR = np.nan
            else:
                O_H_N2_C20_total = np.nan
                O_H_N2_C20_total_ERR = np.nan
        else:
            O_H_N2_C20_total = np.nan
            O_H_N2_C20_total_ERR = np.nan
    else:
        O_H_N2_C20_total = np.nan
        O_H_N2_C20_total_ERR = np.nan
else:
    O_H_N2_C20_total = np.nan
    O_H_N2_C20_total_ERR = np.nan

# S2-C20 (Curti et al. 2020) total metallicity calculation
if np.isfinite(SII6716_FLUX_corr_total) and np.isfinite(SII6730_FLUX_corr_total) and \
   np.isfinite(HA6562_FLUX_corr_total) and \
   SII6716_FLUX_corr_total > 0 and SII6730_FLUX_corr_total > 0 and HA6562_FLUX_corr_total > 0:
    
    # Calculate total flux errors
    HA6562_FLUX_ERR_total = np.sqrt(np.nansum(HA6562_FLUX_ERR[~np.isnan(HA6562_FLUX_corr)]**2))
    SII6716_FLUX_ERR_total = np.sqrt(np.nansum(SII6716_FLUX_ERR[~np.isnan(SII6716_FLUX_corr)]**2))
    SII6730_FLUX_ERR_total = np.sqrt(np.nansum(SII6730_FLUX_ERR[~np.isnan(SII6730_FLUX_corr)]**2))
    
    s2_lin_total_c20 = (SII6716_FLUX_corr_total + SII6730_FLUX_corr_total) / HA6562_FLUX_corr_total
    if s2_lin_total_c20 > 0:
        y_total_c20 = np.log10(s2_lin_total_c20)
        
        # Calculate error in S2 using specialized error propagation
        s2_error_total = s2_error_propagation(SII6716_FLUX_corr_total, SII6716_FLUX_ERR_total,
                                              SII6730_FLUX_corr_total, SII6730_FLUX_ERR_total,
                                              HA6562_FLUX_corr_total, HA6562_FLUX_ERR_total)
        
        c0 = -0.442; c1 = -0.360; c2 = -6.271; c3 = -8.339; c4 = -3.559
        poly_coeffs = [c4, c3, c2, c1, (c0 - y_total_c20)]
        roots = np.roots(poly_coeffs)
        realish = roots[np.abs(roots.imag) <= 1e-8].real
        if realish.size > 0:
            in_range = realish[(realish >= -0.7) & (realish <= 0.3)]
            if in_range.size > 0:
                x_final_total = np.min(in_range)
                O_H_S2_C20_total = x_final_total + 8.69
                
                # Error propagation: derivative of 4th-order polynomial with respect to y
                derivative_y = np.abs(c1 + 2*c2*x_final_total + 3*c3*x_final_total**2 + 4*c4*x_final_total**3)
                if derivative_y > 0:
                    derivative_x = 1.0 / derivative_y
                    obs_error = derivative_x * s2_error_total
                    O_H_S2_C20_total_ERR = obs_error
                else:
                    O_H_S2_C20_total_ERR = np.nan
            else:
                O_H_S2_C20_total = np.nan
                O_H_S2_C20_total_ERR = np.nan
        else:
            O_H_S2_C20_total = np.nan
            O_H_S2_C20_total_ERR = np.nan
    else:
        O_H_S2_C20_total = np.nan
        O_H_S2_C20_total_ERR = np.nan
else:
    O_H_S2_C20_total = np.nan
    O_H_S2_C20_total_ERR = np.nan

# Total Combined C20 metallicity
# Calculate metallicity using integrated line fluxes for each C20 method and select minimum error
if np.any(np.isfinite(O_H_COMBINED_C20)):
    # Get integrated line fluxes from total regions (these are calculated above)
    
    # Check if we have valid integrated fluxes for C20 calculations
    has_valid_fluxes = (
        np.isfinite(HB4861_FLUX_corr_total) and np.isfinite(HA6562_FLUX_corr_total) and
        np.isfinite(OIII5006_FLUX_corr_total) and np.isfinite(NII6583_FLUX_corr_total) and
        np.isfinite(SII6716_FLUX_corr_total) and np.isfinite(SII6730_FLUX_corr_total) and
        HB4861_FLUX_corr_total > 0 and HA6562_FLUX_corr_total > 0 and
        OIII5006_FLUX_corr_total > 0 and NII6583_FLUX_corr_total > 0 and
        SII6716_FLUX_corr_total > 0 and SII6730_FLUX_corr_total > 0
    )
    
    if has_valid_fluxes:
        # Calculate metallicity and error for each C20 method using integrated fluxes
        c20_methods = []
        c20_errors = []
        c20_names = ['O3N2-C20', 'O3S2-C20', 'RS32-C20', 'R3-C20', 'N2-C20', 'S2-C20']
        
        # Method 0: O3N2-C20 (use calculated total values)
        if np.isfinite(O_H_O3N2_C20_total) and np.isfinite(O_H_O3N2_C20_total_ERR):
            c20_methods.append(O_H_O3N2_C20_total)
            c20_errors.append(O_H_O3N2_C20_total_ERR)
        else:
            c20_methods.append(np.nan)
            c20_errors.append(np.inf)
        
        # Method 1: O3S2-C20 (use calculated total values)
        if np.isfinite(O_H_O3S2_C20_total) and np.isfinite(O_H_O3S2_C20_total_ERR):
            c20_methods.append(O_H_O3S2_C20_total)
            c20_errors.append(O_H_O3S2_C20_total_ERR)
        else:
            c20_methods.append(np.nan)
            c20_errors.append(np.inf)
        
        # Method 2: RS32-C20 (use calculated total values)
        if np.isfinite(O_H_RS32_C20_total) and np.isfinite(O_H_RS32_C20_total_ERR):
            c20_methods.append(O_H_RS32_C20_total)
            c20_errors.append(O_H_RS32_C20_total_ERR)
        else:
            c20_methods.append(np.nan)
            c20_errors.append(np.inf)
        
        # Method 3: R3-C20 (use calculated total values)
        if np.isfinite(O_H_R3_C20_total) and np.isfinite(O_H_R3_C20_total_ERR):
            c20_methods.append(O_H_R3_C20_total)
            c20_errors.append(O_H_R3_C20_total_ERR)
        else:
            c20_methods.append(np.nan)
            c20_errors.append(np.inf)
        
        # Method 4: N2-C20 (use calculated total values)
        if np.isfinite(O_H_N2_C20_total) and np.isfinite(O_H_N2_C20_total_ERR):
            c20_methods.append(O_H_N2_C20_total)
            c20_errors.append(O_H_N2_C20_total_ERR)
        else:
            c20_methods.append(np.nan)
            c20_errors.append(np.inf)
        
        # Method 5: S2-C20 (use calculated total values)
        if np.isfinite(O_H_S2_C20_total) and np.isfinite(O_H_S2_C20_total_ERR):
            c20_methods.append(O_H_S2_C20_total)
            c20_errors.append(O_H_S2_C20_total_ERR)
        else:
            c20_methods.append(np.nan)
            c20_errors.append(np.inf)
        
        # Find method with minimum error among valid results
        valid_indices = [i for i, (oh, err) in enumerate(zip(c20_methods, c20_errors)) if np.isfinite(oh) and np.isfinite(err)]
        
        if valid_indices:
            min_error_idx = min(valid_indices, key=lambda i: c20_errors[i])
            O_H_COMBINED_C20_total = c20_methods[min_error_idx]
            best_method_name = c20_names[min_error_idx]
            print(f"Combined C20 total: Selected {best_method_name} with error {c20_errors[min_error_idx]:.3f}")
            print(f"Available methods: {[(c20_names[i], c20_methods[i], c20_errors[i]) for i in valid_indices]}")
        else:
            O_H_COMBINED_C20_total = np.nan
    else:
        O_H_COMBINED_C20_total = np.nan
else:
    O_H_COMBINED_C20_total = np.nan

# ------------------------------------------------------------------
# 12.  Output the results
# ------------------------------------------------------------------

with fits.open(gas_path) as hdul:
    new_hdul = fits.HDUList([hdu.copy() for hdu in hdul])

# Add provenance information to primary header
new_hdul[0].header['BPTMODE'] = 'both'
new_hdul[0].header['CUT_SN'] = cut
new_hdul[0].header['NOISE'] = noise  # in 1e-20 erg s-1 cm-2
new_hdul[0].header['DIST_MPC'] = 16.5

# Gas E(B-V)
hdu_E_BV_BD = fits.ImageHDU(E_BV_BD.astype(np.float64),
                           header=gas_header, name="Gas_E_BV_BD")
hdu_E_BV_BD.header['BUNIT'] = 'mag'
new_hdul.append(hdu_E_BV_BD)
# Corrected line fluxes
hdu_HB4861_FLUX_corr = fits.ImageHDU(HB4861_FLUX_corr.astype(np.float64),
                                     header=gas_header, name="HB4861_FLUX_corr")
hdu_HB4861_FLUX_corr.header['BUNIT'] = '1e-20 erg s-1 cm-2'
new_hdul.append(hdu_HB4861_FLUX_corr)
hdu_HA6562_FLUX_corr = fits.ImageHDU(HA6562_FLUX_corr.astype(np.float64),
                                     header=gas_header, name="HA6562_FLUX_corr")
hdu_HA6562_FLUX_corr.header['BUNIT'] = '1e-20 erg s-1 cm-2'
new_hdul.append(hdu_HA6562_FLUX_corr)
hdu_OIII5006_FLUX_corr = fits.ImageHDU(OIII5006_FLUX_corr.astype(np.float64),
                                       header=gas_header, name="OIII5006_FLUX_corr")
hdu_OIII5006_FLUX_corr.header['BUNIT'] = '1e-20 erg s-1 cm-2'
new_hdul.append(hdu_OIII5006_FLUX_corr)
hdu_NII6583_FLUX_corr = fits.ImageHDU(NII6583_FLUX_corr.astype(np.float64),
                                      header=gas_header, name="NII6583_FLUX_corr")
hdu_NII6583_FLUX_corr.header['BUNIT'] = '1e-20 erg s-1 cm-2'
new_hdul.append(hdu_NII6583_FLUX_corr)
hdu_SII6716_FLUX_corr = fits.ImageHDU(SII6716_FLUX_corr.astype(np.float64),
                                      header=gas_header, name="SII6716_FLUX_corr")
hdu_SII6716_FLUX_corr.header['BUNIT'] = '1e-20 erg s-1 cm-2'
new_hdul.append(hdu_SII6716_FLUX_corr)
hdu_SII6730_FLUX_corr = fits.ImageHDU(SII6730_FLUX_corr.astype(np.float64),
                                      header=gas_header, name="SII6730_FLUX_corr")
hdu_SII6730_FLUX_corr.header['BUNIT'] = '1e-20 erg s-1 cm-2'
new_hdul.append(hdu_SII6730_FLUX_corr)
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
# [O/H]
hdu_O_H_D16_SF = fits.ImageHDU(O_H_D16_SF.astype(np.float64),
                             header=gas_header, name="O_H_D16_SF")
hdu_O_H_D16_SF.header['BUNIT'] = '12+log(O/H)'
new_hdul.append(hdu_O_H_D16_SF)
hdu_O_H_PG16_SF = fits.ImageHDU(O_H_PG16_SF.astype(np.float64),
                             header=gas_header, name="O_H_PG16_SF")
hdu_O_H_PG16_SF.header['BUNIT'] = '12+log(O/H)'
new_hdul.append(hdu_O_H_PG16_SF)
hdu_O_H_O3N2_M13_SF = fits.ImageHDU(O_H_O3N2_M13_SF.astype(np.float64),
                             header=gas_header, name="O_H_O3N2_M13_SF")
hdu_O_H_O3N2_M13_SF.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_O3N2_M13_SF.header['COMMENT'] = 'O3N2-M13 (Marino et al. 2013) metallicity in SF regions'
new_hdul.append(hdu_O_H_O3N2_M13_SF)
hdu_O_H_N2_M13_SF = fits.ImageHDU(O_H_N2_M13_SF.astype(np.float64),
                             header=gas_header, name="O_H_N2_M13_SF")
hdu_O_H_N2_M13_SF.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_N2_M13_SF.header['COMMENT'] = 'N2-M13 (Marino et al. 2013) metallicity in SF regions'
new_hdul.append(hdu_O_H_N2_M13_SF)
hdu_O_H_O3N2_PP04_SF = fits.ImageHDU(O_H_O3N2_PP04_SF.astype(np.float64),
                             header=gas_header, name="O_H_O3N2_PP04_SF")
hdu_O_H_O3N2_PP04_SF.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_O3N2_PP04_SF.header['COMMENT'] = 'O3N2-PP04 (Pettini & Pagel 2004) metallicity in SF regions'
new_hdul.append(hdu_O_H_O3N2_PP04_SF)
hdu_O_H_N2_PP04_SF = fits.ImageHDU(O_H_N2_PP04_SF.astype(np.float64),
                             header=gas_header, name="O_H_N2_PP04_SF")
hdu_O_H_N2_PP04_SF.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_N2_PP04_SF.header['COMMENT'] = 'N2-PP04 (Pettini & Pagel 2004) metallicity in SF regions'
new_hdul.append(hdu_O_H_N2_PP04_SF)
hdu_O_H_O3N2_C20_SF = fits.ImageHDU(O_H_O3N2_C20_SF.astype(np.float64),
                             header=gas_header, name="O_H_O3N2_C20_SF")
hdu_O_H_O3N2_C20_SF.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_O3N2_C20_SF.header['COMMENT'] = 'O3N2-C20 (Curti et al. 2020) metallicity in SF regions'
new_hdul.append(hdu_O_H_O3N2_C20_SF)
hdu_O_H_O3S2_C20_SF = fits.ImageHDU(O_H_O3S2_C20_SF.astype(np.float64),
                             header=gas_header, name="O_H_O3S2_C20_SF")
hdu_O_H_O3S2_C20_SF.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_O3S2_C20_SF.header['COMMENT'] = 'O3S2-C20 (Curti et al. 2020) metallicity in SF regions'
new_hdul.append(hdu_O_H_O3S2_C20_SF)
hdu_O_H_RS32_C20_SF = fits.ImageHDU(O_H_RS32_C20_SF.astype(np.float64),
                             header=gas_header, name="O_H_RS32_C20_SF")
hdu_O_H_RS32_C20_SF.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_RS32_C20_SF.header['COMMENT'] = 'RS32-C20 (Curti et al. 2020) metallicity in SF regions'
new_hdul.append(hdu_O_H_RS32_C20_SF)
hdu_O_H_R3_C20_SF = fits.ImageHDU(O_H_R3_C20_SF.astype(np.float64),
                             header=gas_header, name="O_H_R3_C20_SF")
hdu_O_H_R3_C20_SF.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_R3_C20_SF.header['COMMENT'] = 'R3-C20 (Curti et al. 2020) metallicity in SF regions'
new_hdul.append(hdu_O_H_R3_C20_SF)
hdu_O_H_N2_C20_SF = fits.ImageHDU(O_H_N2_C20_SF.astype(np.float64),
                             header=gas_header, name="O_H_N2_C20_SF")
hdu_O_H_N2_C20_SF.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_N2_C20_SF.header['COMMENT'] = 'N2-C20 (Curti et al. 2020) metallicity in SF regions'
new_hdul.append(hdu_O_H_N2_C20_SF)
hdu_O_H_S2_C20_SF = fits.ImageHDU(O_H_S2_C20_SF.astype(np.float64),
                             header=gas_header, name="O_H_S2_C20_SF")
hdu_O_H_S2_C20_SF.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_S2_C20_SF.header['COMMENT'] = 'S2-C20 (Curti et al. 2020) metallicity in SF regions'
new_hdul.append(hdu_O_H_S2_C20_SF)

hdu_O_H_COMBINED_C20_SF = fits.ImageHDU(O_H_COMBINED_C20_SF.astype(np.float64),
                             header=gas_header, name="O_H_COMBINED_C20_SF")
hdu_O_H_COMBINED_C20_SF.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_COMBINED_C20_SF.header['COMMENT'] = 'Combined C20 metallicity (best method per spaxel) in SF regions'
new_hdul.append(hdu_O_H_COMBINED_C20_SF)

hdu_COMBINED_C20_METHOD = fits.ImageHDU(combined_c20_method_map.astype(np.float64),
                             header=gas_header, name="COMBINED_C20_METHOD")
hdu_COMBINED_C20_METHOD.header['BUNIT'] = 'method_index'
hdu_COMBINED_C20_METHOD.header['COMMENT'] = 'Method used for Combined C20: 0=O3N2, 1=O3S2, 2=RS32, 3=R3, 4=N2, 5=S2'
new_hdul.append(hdu_COMBINED_C20_METHOD)

new_hdul.writeto(out_path, overwrite=True)
print("Extended file written ➜", out_path.resolve())

# ------------------------------------------------------------------
# 13.  Print some useful information
# ------------------------------------------------------------------

# Print the total non-nan spaxels
print("--------------------------------------------------------------")
total_spaxels = np.sum(np.isfinite(V_STARS2))
print("Total non-nan spaxels:", total_spaxels)
# Print the number of 6 cases that need number, 2 upper cases, and 4 unclassified cases
print("Number of pixels with Halpha not detected:", np.sum(HA_not_detected))
print("Number of pixels with Halpha detected, Hbeta not detected:", np.sum(HA_detected_HB_not_detected))
print("Number of pixels with Halpha detected, Hbeta detected, NII not detected, OIII not detected and unclassified:", 
      np.sum(HA_detected_HB_detected_NII_not_detected_OIII_not_detected & mask_N2_unclassified))
print("Number of pixels with Halpha detected, Hbeta detected, NII not detected, OIII detected and unclassified:",
      np.sum(HA_detected_HB_detected_NII_not_detected_OIII_detected & mask_N2_unclassified))
print("Number of pixels with Halpha detected, Hbeta detected, NII detected, OIII not detected and unclassified:",
      np.sum(HA_detected_HB_detected_NII_detected_OIII_not_detected & mask_N2_unclassified))
print("Number of pixels with Halpha detected, Hbeta detected, NII detected, OIII detected and unclassified:",
      np.sum(HA_detected_HB_detected_NII_detected_OIII_detected & mask_N2_unclassified))
print("Number of pixels with Halpha detected, Hbeta detected, SII not detected, OIII not detected and unclassified:", 
      np.sum(HA_detected_HB_detected_SII_not_detected_OIII_not_detected & mask_N2_unclassified))
print("Number of pixels with Halpha detected, Hbeta detected, SII not detected, OIII detected and unclassified:",
      np.sum(HA_detected_HB_detected_SII_not_detected_OIII_detected & mask_N2_unclassified))
print("Number of pixels with Halpha detected, Hbeta detected, SII detected, OIII not detected and unclassified:",
      np.sum(HA_detected_HB_detected_SII_detected_OIII_not_detected & mask_N2_unclassified))
print("Number of pixels with Halpha detected, Hbeta detected, SII detected, OIII detected and unclassified:",
      np.sum(HA_detected_HB_detected_SII_detected_OIII_detected & mask_N2_unclassified))
print("--------------------------------------------------------------")
print(f"Total corrected Halpha luminosity: {np.nansum(HA6562_LUM):.2e} erg/s")
print(f"Total corrected Halpha luminosity from SF region: {np.nansum(HA6562_LUM[mask_SF]):.2e} erg/s")
print(f"Total Halpha SFR: {np.nansum(SFR_map):.2f} M☉/yr or in log10 scale: {np.log10(np.nansum(SFR_map)):.2f} log(M☉/yr)")
print(f"Total Halpha SFR from SF region: {np.nansum(SFR_map[mask_SF]):.2f} M☉/yr or in log10 scale: {np.log10(np.nansum(SFR_map[mask_SF])):.2f} log(M☉/yr)")
print("--------------------------------------------------------------")
print("[O/H] D16 SF: Total metallicity in SF region: ", O_H_D16_SF_total)
print("[O/H] PG16 SF: Total metallicity in SF region: ", O_H_PG16_SF_total)
print("[O/H] O3N2-M13 SF: Total metallicity in SF region: ", O_H_O3N2_M13_SF_total)
print("[O/H] N2-M13 SF: Total metallicity in SF region: ", O_H_N2_M13_SF_total)
print("[O/H] O3N2-PP04 SF: Total metallicity in SF region: ", O_H_O3N2_PP04_SF_total)
print("[O/H] N2-PP04 SF: Total metallicity in SF region: ", O_H_N2_PP04_SF_total)
print("[O/H] O3N2-C20 SF: Total metallicity in SF region: ", O_H_O3N2_C20_SF_total)
print("[O/H] O3S2-C20 SF: Total metallicity in SF region: ", O_H_O3S2_C20_SF_total)
print("[O/H] RS32-C20 SF: Total metallicity in SF region: ", O_H_RS32_C20_SF_total)
print("[O/H] R3-C20 SF: Total metallicity in SF region: ", O_H_R3_C20_SF_total)
print("[O/H] N2-C20 SF: Total metallicity in SF region: ", O_H_N2_C20_SF_total)
print("[O/H] S2-C20 SF: Total metallicity in SF region: ", O_H_S2_C20_SF_total)
print("[O/H] Combined-C20 SF: Total metallicity in SF region: ", O_H_COMBINED_C20_SF_total)
print("--------------------------------------------------------------")
print("[O/H] D16: Total metallicity in total region: ", O_H_D16_total)
print("[O/H] PG16: Total metallicity in total region: ", O_H_PG16_total)
print("[O/H] O3N2-M13: Total metallicity in total region: ", O_H_O3N2_M13_total)
print("[O/H] N2-M13: Total metallicity in total region: ", O_H_N2_M13_total)
print("[O/H] O3N2-PP04: Total metallicity in total region: ", O_H_O3N2_PP04_total)
print("[O/H] N2-PP04: Total metallicity in total region: ", O_H_N2_PP04_total)
print("[O/H] O3N2-C20: Total metallicity in total region: ", O_H_O3N2_C20_total) 
print("[O/H] O3S2-C20: Total metallicity in total region: ", O_H_O3S2_C20_total)
print("[O/H] RS32-C20: Total metallicity in total region: ", O_H_RS32_C20_total)
print("[O/H] R3-C20: Total metallicity in total region: ", O_H_R3_C20_total)
print("[O/H] N2-C20: Total metallicity in total region: ", O_H_N2_C20_total)
print("[O/H] S2-C20: Total metallicity in total region: ", O_H_S2_C20_total)
print("[O/H] Combined-C20: Total metallicity in total region: ", O_H_COMBINED_C20_total)
print("--------------------------------------------------------------")


print(f"Total runtime: {time.perf_counter() - t0:.1f} s")

# End of file