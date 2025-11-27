#!/usr/bin/env python
"""
PHANGS_Z.py - Adapted for PHANGS data (1D catalogs) to calculate metallicity in HII regions.
"""

import argparse
import logging
import sys
import time
import glob
import os
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy import units as u
import warnings

# Suppress warnings
warnings.filterwarnings('ignore')

# ------------------------------------------------------------------
# Helper Functions
# ------------------------------------------------------------------

# Calzetti (2000) curve
def calzetti_k(w_um):
    """Return k(λ) = A(λ)/E(B−V) for Calzetti (2000); wavelengths in microns."""
    w = np.asarray(w_um, dtype=float)
    Rv = 4.05
    k = np.empty_like(w, dtype=float)

    short = (w >= 0.12) & (w < 0.63)
    long  = (w >= 0.63) & (w <= 2.2)

    if w.ndim == 0:
        if short:
            return 2.659 * (-2.156 + 1.509/w - 0.198/w**2 + 0.011/w**3) + Rv
        elif long:
            return 2.659 * (-1.857 + 1.040/w) + Rv
        else:
            return 0.0
    
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

def convert_bd_to_ebv(BD, k_HB4861, k_HA6562, R_int=2.86):
    E_BV_BD = 2.5 / (k_HB4861 - k_HA6562) * np.log10(BD / R_int)
    return E_BV_BD

def correct_flux_with_ebv(flux, ebv, k):
    """Correct flux with gas E(B-V) and extinction coefficient k."""
    return flux * 10**(0.4 * k * ebv)

def bpt_error_propagation(numerator, denominator, numerator_err, denominator_err):
    """
    Calculate the propagated error for the BPT ratio log10(numerator/denominator).
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = numerator / denominator
        # log_ratio = np.log10(ratio) # Not used directly
        log_ratio_err = 1/(np.log(10)) * np.sqrt((numerator_err / numerator)**2 + (denominator_err / denominator)**2)
        return log_ratio_err

def s2_error_propagation(sii6716_flux, sii6716_flux_err, sii6730_flux, sii6730_flux_err, ha6563_flux, ha6563_flux_err):
    """Calculate propagated error for log10(([SII]6716 + [SII]6730) / Hα)"""
    numerator = sii6716_flux + sii6730_flux
    numerator_err = np.sqrt(sii6716_flux_err**2 + sii6730_flux_err**2)
    ratio_rel_err = np.sqrt((numerator_err / numerator)**2 + (ha6563_flux_err / ha6563_flux)**2)
    log_ratio_err = ratio_rel_err / np.log(10)
    return log_ratio_err

# ------------------------------------------------------------------
# Metallicity Calculation Functions (Dimension Agnostic)
# ------------------------------------------------------------------

def calculate_n2s2_n06_metallicity(nii6583_flux, ha6562_flux, sii6716_flux, sii6730_flux):
    good_mask = (np.isfinite(nii6583_flux) & np.isfinite(ha6562_flux) &
                np.isfinite(sii6716_flux) & np.isfinite(sii6730_flux) &
                (nii6583_flux > 0) & (ha6562_flux > 0) &
                (sii6716_flux > 0) & (sii6730_flux > 0))
    
    oh_n2s2_n06 = np.full_like(nii6583_flux, np.nan)
    sii_total = sii6716_flux + sii6730_flux
    n2s2_ratio = np.log10(nii6583_flux / sii_total)
    
    c3 = 0.17963; c2 = 0.58181; c1 = 0.74100; c0 = -0.25214
    
    if np.any(good_mask):
        indices = np.where(good_mask)
        # Handle 1D vs ND
        if len(indices) == 1:
             iterator = indices[0]
             for i in iterator:
                n2s2_val = n2s2_ratio[i]
                poly_coeffs = [c3, c2, c1, (c0 - n2s2_val)]
                roots = np.roots(poly_coeffs)
                real_roots = roots[np.isreal(roots)].real
                if len(real_roots) > 0:
                    oh_n2s2_n06[i] = real_roots[0] + 8.69
        else:
             for idx in zip(*indices):
                n2s2_val = n2s2_ratio[idx]
                poly_coeffs = [c3, c2, c1, (c0 - n2s2_val)]
                roots = np.roots(poly_coeffs)
                real_roots = roots[np.isreal(roots)].real
                if len(real_roots) > 0:
                    oh_n2s2_n06[idx] = real_roots[0] + 8.69
    
    return oh_n2s2_n06, good_mask

def calculate_o3n2_m13_metallicity(hb4861_flux, oiii5006_flux, nii6583_flux, ha6562_flux):
    good_mask = (np.isfinite(hb4861_flux) & np.isfinite(oiii5006_flux) &
                np.isfinite(nii6583_flux) & np.isfinite(ha6562_flux) &
                (hb4861_flux > 0) & (oiii5006_flux > 0) &
                (nii6583_flux > 0) & (ha6562_flux > 0))
    
    oh_o3n2_m13 = np.full_like(hb4861_flux, np.nan)
    oiii_hb = oiii5006_flux / hb4861_flux
    nii_ha = nii6583_flux / ha6562_flux
    o3n2_ratio = np.log10(oiii_hb / nii_ha)
    oh_o3n2_m13[good_mask] = 8.533 - 0.214 * o3n2_ratio[good_mask]
    return oh_o3n2_m13, good_mask

def calculate_n2_m13_metallicity(nii6583_flux, ha6562_flux):
    good_mask = (np.isfinite(nii6583_flux) & np.isfinite(ha6562_flux) &
                (nii6583_flux > 0) & (ha6562_flux > 0))
    oh_n2_m13 = np.full_like(nii6583_flux, np.nan)
    n2_ratio = np.log10(nii6583_flux / ha6562_flux)
    oh_n2_m13[good_mask] = 8.743 + 0.462 * n2_ratio[good_mask]
    return oh_n2_m13, good_mask

def calculate_o3n2_pp04_metallicity(hb4861_flux, oiii5006_flux, nii6583_flux, ha6562_flux):
    good_mask = (np.isfinite(hb4861_flux) & np.isfinite(oiii5006_flux) &
                np.isfinite(nii6583_flux) & np.isfinite(ha6562_flux) &
                (hb4861_flux > 0) & (oiii5006_flux > 0) &
                (nii6583_flux > 0) & (ha6562_flux > 0))
    oh_o3n2_pp04 = np.full_like(hb4861_flux, np.nan)
    oiii_hb = oiii5006_flux / hb4861_flux
    nii_ha = nii6583_flux / ha6562_flux
    o3n2_ratio = np.log10(oiii_hb / nii_ha)
    oh_o3n2_pp04[good_mask] = 8.73 - 0.32 * o3n2_ratio[good_mask]
    return oh_o3n2_pp04, good_mask

def calculate_n2_pp04_metallicity(nii6583_flux, ha6562_flux):
    good_mask = (np.isfinite(nii6583_flux) & np.isfinite(ha6562_flux) &
                (nii6583_flux > 0) & (ha6562_flux > 0))
    oh_n2_pp04 = np.full_like(nii6583_flux, np.nan)
    n2_ratio = np.log10(nii6583_flux / ha6562_flux)
    oh_n2_pp04[good_mask] = (9.37 + 2.03 * n2_ratio[good_mask] + 
                            1.26 * n2_ratio[good_mask]**2 + 
                            0.32 * n2_ratio[good_mask]**3)
    return oh_n2_pp04, good_mask

def calculate_o3n2_c20_metallicity(hb4861_flux, oiii5006_flux, nii6583_flux, ha6562_flux, 
                                   hb4861_flux_err, oiii5006_flux_err, nii6583_flux_err, ha6562_flux_err):
    good_mask = (np.isfinite(hb4861_flux) & np.isfinite(oiii5006_flux) &
                np.isfinite(nii6583_flux) & np.isfinite(ha6562_flux) &
                (hb4861_flux > 0) & (oiii5006_flux > 0) &
                (nii6583_flux > 0) & (ha6562_flux > 0) &
                 np.isfinite(hb4861_flux_err) & np.isfinite(oiii5006_flux_err) & 
                 np.isfinite(nii6583_flux_err) & np.isfinite(ha6562_flux_err))
    
    oh_o3n2_c20 = np.full_like(hb4861_flux, np.nan)
    oh_o3n2_c20_err = np.full_like(hb4861_flux, np.nan)
    
    oiii_hb = oiii5006_flux / hb4861_flux
    nii_ha = nii6583_flux / ha6562_flux
    o3n2_ratio = np.log10(oiii_hb / nii_ha)
    
    oiii_hb_err = bpt_error_propagation(oiii5006_flux, hb4861_flux, oiii5006_flux_err, hb4861_flux_err)
    nii_ha_err = bpt_error_propagation(nii6583_flux, ha6562_flux, nii6583_flux_err, ha6562_flux_err)
    o3n2_ratio_err = np.sqrt(oiii_hb_err**2 + nii_ha_err**2)
    
    R = o3n2_ratio
    y = R
    y_err = o3n2_ratio_err
    c0 = 0.281; c1 = -4.765; c2 = -2.268
    a = c2; b = c1; c = c0 - y
    discriminant = b**2 - 4*a*c
    valid_discriminant = discriminant >= 0
    combined_mask = good_mask & valid_discriminant
    
    if np.any(combined_mask):
        x_solution1 = (-b + np.sqrt(discriminant[combined_mask])) / (2*a)
        x_solution2 = (-b - np.sqrt(discriminant[combined_mask])) / (2*a)
        x_final = np.where((x_solution1 >= -1.1) & (x_solution1 <= 1.25), x_solution1, x_solution2)
        derivative_x = np.abs(c1 + 2*c2*x_final)
        x_err = y_err[combined_mask] / derivative_x
        oh_o3n2_c20[combined_mask] = x_final + 8.69
        oh_o3n2_c20_err[combined_mask] = np.sqrt(x_err**2)

    return oh_o3n2_c20, oh_o3n2_c20_err, combined_mask

def calculate_o3s2_c20_metallicity(hb4861_flux, oiii5006_flux, sii6716_flux, sii6730_flux, 
                                   hb4861_flux_err, oiii5006_flux_err, sii6716_flux_err, sii6730_flux_err):
    good_mask = (np.isfinite(hb4861_flux) & np.isfinite(oiii5006_flux) &
                np.isfinite(sii6716_flux) & np.isfinite(sii6730_flux) &
                (hb4861_flux > 0) & (oiii5006_flux > 0) &
                (sii6716_flux > 0) & (sii6730_flux > 0) &
                 np.isfinite(hb4861_flux_err) & np.isfinite(oiii5006_flux_err) & 
                 np.isfinite(sii6716_flux_err) & np.isfinite(sii6730_flux_err))
    
    oh_o3s2_c20 = np.full_like(hb4861_flux, np.nan)
    oh_o3s2_c20_err = np.full_like(hb4861_flux, np.nan)
    
    oiii_hb = oiii5006_flux / hb4861_flux
    sii_total = sii6716_flux + sii6730_flux
    sii_total_err = np.sqrt(sii6716_flux_err**2 + sii6730_flux_err**2)
    sii_hb = sii_total / hb4861_flux
    o3s2_ratio = np.log10(oiii_hb / sii_hb)
    
    oiii_hb_err = bpt_error_propagation(oiii5006_flux, hb4861_flux, oiii5006_flux_err, hb4861_flux_err)
    sii_hb_err = bpt_error_propagation(sii_total, hb4861_flux, sii_total_err, hb4861_flux_err)
    o3s2_ratio_err = np.sqrt(oiii_hb_err**2 + sii_hb_err**2)
    
    y = o3s2_ratio
    y_err = o3s2_ratio_err
    c0 = 0.191; c1 = -4.292; c2 = -2.538; c3 = 0.053; c4 = 0.332
    
    combined_mask = np.copy(good_mask)
    
    if np.any(good_mask):
        indices = np.where(good_mask)
        if len(indices) == 1:
             iterator = indices[0]
             for i in iterator:
                y_val = y[i]
                y_err_val = y_err[i]
                poly_coeffs = [c4, c3, c2, c1, (c0 - y_val)]
                roots = np.roots(poly_coeffs)
                real_roots = roots[np.isreal(roots)].real
                if len(real_roots) > 0:
                    reasonable_roots = real_roots[(real_roots >= -2) & (real_roots <= 2)]
                    if len(reasonable_roots) > 0:
                        x_final = reasonable_roots[0]
                        derivative_x = (np.abs(c1 + 2*c2*x_final + 3*c3*x_final**2 + 4*c4*x_final**3))
                        x_err = y_err_val / derivative_x
                        oh_o3s2_c20[i] = x_final + 8.69
                        oh_o3s2_c20_err[i] = np.sqrt(x_err**2)
                    else:
                        combined_mask[i] = False
                else:
                    combined_mask[i] = False
        else:
             for idx in zip(*indices):
                y_val = y[idx]
                y_err_val = y_err[idx]
                poly_coeffs = [c4, c3, c2, c1, (c0 - y_val)]
                roots = np.roots(poly_coeffs)
                real_roots = roots[np.isreal(roots)].real
                if len(real_roots) > 0:
                    reasonable_roots = real_roots[(real_roots >= -2) & (real_roots <= 2)]
                    if len(reasonable_roots) > 0:
                        x_final = reasonable_roots[0]
                        derivative_x = (np.abs(c1 + 2*c2*x_final + 3*c3*x_final**2 + 4*c4*x_final**3))
                        x_err = y_err_val / derivative_x
                        oh_o3s2_c20[idx] = x_final + 8.69
                        oh_o3s2_c20_err[idx] = np.sqrt(x_err**2)
                    else:
                        combined_mask[idx] = False
                else:
                    combined_mask[idx] = False
    
    combined_mask = combined_mask & np.isfinite(oh_o3s2_c20)
    return oh_o3s2_c20, oh_o3s2_c20_err, combined_mask

def calculate_rs32_c20_metallicity(hb4861_flux, ha6563_flux, oiii5006_flux, sii6716_flux, sii6730_flux,
                                   hb4861_flux_err, ha6563_flux_err, oiii5006_flux_err, sii6716_flux_err, sii6730_flux_err,
                                   coeffs=(-0.054, -2.546, -1.970, 0.082, 0.222)):
    c0, c1, c2, c3, c4 = coeffs
    good_mask = (np.isfinite(hb4861_flux) & np.isfinite(ha6563_flux) &
        np.isfinite(oiii5006_flux) & np.isfinite(sii6716_flux) & np.isfinite(sii6730_flux) &
        np.isfinite(hb4861_flux_err) & np.isfinite(ha6563_flux_err) &
        np.isfinite(oiii5006_flux_err) & np.isfinite(sii6716_flux_err) & np.isfinite(sii6730_flux_err) &
        (hb4861_flux > 0) & (ha6563_flux > 0) &
        (oiii5006_flux > 0) & (sii6716_flux > 0) & (sii6730_flux > 0))

    oh_rs32_c20 = np.full_like(hb4861_flux, np.nan)
    oh_rs32_c20_err = np.full_like(hb4861_flux, np.nan)

    if np.any(good_mask):
        oiii_hb = oiii5006_flux[good_mask] / hb4861_flux[good_mask]
        sii_total = sii6716_flux[good_mask] + sii6730_flux[good_mask]
        sii_total_err = np.sqrt(sii6716_flux_err[good_mask]**2 + sii6730_flux_err[good_mask]**2)
        sii_ha = sii_total / ha6563_flux[good_mask]
        r_lin = oiii_hb + sii_ha
        r_lin = np.where(r_lin > 0, r_lin, np.nan)
        y = np.log10(r_lin)
        
        oiii_hb_err = bpt_error_propagation(oiii5006_flux[good_mask], hb4861_flux[good_mask], 
                                            oiii5006_flux_err[good_mask], hb4861_flux_err[good_mask])
        sii_ha_err = bpt_error_propagation(sii_total, ha6563_flux[good_mask], 
                                           sii_total_err, ha6563_flux_err[good_mask])
        r_lin_err = np.sqrt(oiii_hb_err**2 + sii_ha_err**2)
        y_err = (1/np.log(10)) * (r_lin_err / r_lin)

        indices = np.where(good_mask)
        if len(indices) == 1:
             iterator = indices[0]
             for idx_in_good, i in enumerate(iterator):
                y_val = y[idx_in_good]
                y_err_val = y_err[idx_in_good]
                if not np.isfinite(y_val): continue
                roots = np.roots([c4, c3, c2, c1, (c0 - y_val)])
                real = roots[np.isreal(roots)].real
                if real.size:
                    phys = real[(real >= -0.7) & (real <= 0.3)]
                    cand = phys if phys.size else real
                    y_pred = c0 + c1*cand + c2*cand**2 + c3*cand**3 + c4*cand**4
                    x_final = cand[np.argmin(np.abs(y_pred - y_val))]
                    derivative_x = (np.abs(c1 + 2*c2*x_final + 3*c3*x_final**2 + 4*c4*x_final**3))
                    x_err = y_err_val / derivative_x
                    oh_rs32_c20[i] = x_final + 8.69
                    oh_rs32_c20_err[i] = np.sqrt(x_err**2)
        else:
             for idx_in_good, idx in enumerate(zip(*indices)):
                y_val = y[idx_in_good]
                y_err_val = y_err[idx_in_good]
                if not np.isfinite(y_val): continue
                roots = np.roots([c4, c3, c2, c1, (c0 - y_val)])
                real = roots[np.isreal(roots)].real
                if real.size:
                    phys = real[(real >= -0.7) & (real <= 0.3)]
                    cand = phys if phys.size else real
                    y_pred = c0 + c1*cand + c2*cand**2 + c3*cand**3 + c4*cand**4
                    x_final = cand[np.argmin(np.abs(y_pred - y_val))]
                    derivative_x = (np.abs(c1 + 2*c2*x_final + 3*c3*x_final**2 + 4*c4*x_final**3))
                    x_err = y_err_val / derivative_x
                    oh_rs32_c20[idx] = x_final + 8.69
                    oh_rs32_c20_err[idx] = np.sqrt(x_err**2)

    combined_mask = good_mask & np.isfinite(oh_rs32_c20)
    return oh_rs32_c20, oh_rs32_c20_err, combined_mask

def calculate_r3_c20_metallicity(hb4861_flux, hb4861_flux_err, oiii5006_flux, oiii5006_flux_err,
                                 coeffs=(-0.277, -3.549, -3.593, -0.981)):
    c0, c1, c2, c3 = coeffs
    good_mask = (np.isfinite(hb4861_flux) & np.isfinite(oiii5006_flux) &
        (hb4861_flux > 0) & (oiii5006_flux > 0) &
        np.isfinite(hb4861_flux_err) & np.isfinite(oiii5006_flux_err) &
        (hb4861_flux_err > 0) & (oiii5006_flux_err > 0))

    oh_r3_c20 = np.full_like(hb4861_flux, np.nan)
    oh_r3_c20_err = np.full_like(hb4861_flux, np.nan)

    if np.any(good_mask):
        r_lin = (oiii5006_flux[good_mask] / hb4861_flux[good_mask])
        r_lin = np.where(r_lin > 0, r_lin, np.nan)
        y = np.log10(r_lin)
        r3_error = bpt_error_propagation(oiii5006_flux[good_mask], hb4861_flux[good_mask],
            oiii5006_flux_err[good_mask], hb4861_flux_err[good_mask])

        indices = np.where(good_mask)
        if len(indices) == 1:
             iterator = indices[0]
             for idx_in_good, i in enumerate(iterator):
                y_val = y[idx_in_good]
                if not np.isfinite(y_val): continue
                roots = np.roots([c3, c2, c1, (c0 - y_val)])
                real = roots[np.isreal(roots)].real
                if real.size:
                    phys = real[(real >= -0.7) & (real <= 0.3)]
                    cand = phys if phys.size else real
                    y_pred = c0 + c1*cand + c2*cand**2 + c3*cand**3
                    x_final = cand[np.argmin(np.abs(y_pred - y_val))]
                    oh_r3_c20[i] = x_final + 8.69
                    derivative_y = np.abs(c1 + 2*c2*x_final + 3*c3*x_final**2)
                    if derivative_y > 0:
                        derivative_x = 1.0 / derivative_y
                        obs_error = derivative_x * r3_error[idx_in_good]
                        oh_r3_c20_err[i] = np.sqrt(obs_error**2)
        else:
             for idx_in_good, idx in enumerate(zip(*indices)):
                y_val = y[idx_in_good]
                if not np.isfinite(y_val): continue
                roots = np.roots([c3, c2, c1, (c0 - y_val)])
                real = roots[np.isreal(roots)].real
                if real.size:
                    phys = real[(real >= -0.7) & (real <= 0.3)]
                    cand = phys if phys.size else real
                    y_pred = c0 + c1*cand + c2*cand**2 + c3*cand**3
                    x_final = cand[np.argmin(np.abs(y_pred - y_val))]
                    oh_r3_c20[idx] = x_final + 8.69
                    derivative_y = np.abs(c1 + 2*c2*x_final + 3*c3*x_final**2)
                    if derivative_y > 0:
                        derivative_x = 1.0 / derivative_y
                        obs_error = derivative_x * r3_error[idx_in_good]
                        oh_r3_c20_err[idx] = np.sqrt(obs_error**2)

    combined_mask = good_mask & np.isfinite(oh_r3_c20)
    return oh_r3_c20, oh_r3_c20_err, combined_mask

def calculate_n2_c20_metallicity(ha6563_flux, ha6563_flux_err, nii6584_flux, nii6584_flux_err,
                                 coeffs=(-0.489, 1.513, -2.554, -5.293, -2.867)):
    c0, c1, c2, c3, c4 = coeffs
    good_mask = (np.isfinite(ha6563_flux) & np.isfinite(nii6584_flux) &
        (ha6563_flux > 0) & (nii6584_flux > 0) &
        np.isfinite(ha6563_flux_err) & np.isfinite(nii6584_flux_err) &
        (ha6563_flux_err > 0) & (nii6584_flux_err > 0))

    oh_n2_c20 = np.full_like(ha6563_flux, np.nan)
    oh_n2_c20_err = np.full_like(ha6563_flux, np.nan)

    if np.any(good_mask):
        n2_lin = nii6584_flux[good_mask] / ha6563_flux[good_mask]
        n2_lin = np.where(n2_lin > 0, n2_lin, np.nan)
        y = np.log10(n2_lin)
        n2_error = bpt_error_propagation(nii6584_flux[good_mask], ha6563_flux[good_mask],
            nii6584_flux_err[good_mask], ha6563_flux_err[good_mask])

        indices = np.where(good_mask)
        REAL_ATOL = 1e-8
        RANGE_EPS = 0.0

        if len(indices) == 1:
             iterator = indices[0]
             for idx_in_good, i in enumerate(iterator):
                y_val = y[idx_in_good]
                if not np.isfinite(y_val): continue
                roots = np.roots([c4, c3, c2, c1, (c0 - y_val)])
                realish = roots[np.abs(roots.imag) <= REAL_ATOL].real
                if realish.size == 0:
                    k = np.argmin(np.abs(roots.imag))
                    if np.abs(roots[k].imag) <= 1e-6: realish = np.array([roots[k].real])
                if realish.size == 0: continue
                in_rng = realish[(realish >= -0.7 - RANGE_EPS) & (realish <= 0.3 + RANGE_EPS)]
                if in_rng.size == 0: continue
                x_final = np.min(in_rng)
                oh_n2_c20[i] = x_final + 8.69
                derivative_y = np.abs(c1 + 2*c2*x_final + 3*c3*x_final**2 + 4*c4*x_final**3)
                if derivative_y > 0:
                    derivative_x = 1.0 / derivative_y
                    obs_error = derivative_x * n2_error[idx_in_good]
                    oh_n2_c20_err[i] = np.sqrt(obs_error**2)
        else:
             for idx_in_good, idx in enumerate(zip(*indices)):
                y_val = y[idx_in_good]
                if not np.isfinite(y_val): continue
                roots = np.roots([c4, c3, c2, c1, (c0 - y_val)])
                realish = roots[np.abs(roots.imag) <= REAL_ATOL].real
                if realish.size == 0:
                    k = np.argmin(np.abs(roots.imag))
                    if np.abs(roots[k].imag) <= 1e-6: realish = np.array([roots[k].real])
                if realish.size == 0: continue
                in_rng = realish[(realish >= -0.7 - RANGE_EPS) & (realish <= 0.3 + RANGE_EPS)]
                if in_rng.size == 0: continue
                x_final = np.min(in_rng)
                oh_n2_c20[idx] = x_final + 8.69
                derivative_y = np.abs(c1 + 2*c2*x_final + 3*c3*x_final**2 + 4*c4*x_final**3)
                if derivative_y > 0:
                    derivative_x = 1.0 / derivative_y
                    obs_error = derivative_x * n2_error[idx_in_good]
                    oh_n2_c20_err[idx] = np.sqrt(obs_error**2)

    combined_mask = good_mask & np.isfinite(oh_n2_c20)
    return oh_n2_c20, oh_n2_c20_err, combined_mask

def calculate_s2_c20_metallicity(ha6563_flux, ha6563_flux_err, sii6716_flux, sii6716_flux_err,
                                 sii6730_flux, sii6730_flux_err,
                                 coeffs=(-0.442, -0.360, -6.271, -8.339, -3.559)):
    c0, c1, c2, c3, c4 = coeffs
    good_mask = (np.isfinite(ha6563_flux) & np.isfinite(sii6716_flux) & np.isfinite(sii6730_flux) &
        (ha6563_flux > 0) & (sii6716_flux > 0) & (sii6730_flux > 0) &
        np.isfinite(ha6563_flux_err) & np.isfinite(sii6716_flux_err) & np.isfinite(sii6730_flux_err) &
        (ha6563_flux_err > 0) & (sii6716_flux_err > 0) & (sii6730_flux_err > 0))

    oh_s2_c20 = np.full_like(ha6563_flux, np.nan)
    oh_s2_c20_err = np.full_like(ha6563_flux, np.nan)

    if np.any(good_mask):
        s2_lin = (sii6716_flux[good_mask] + sii6730_flux[good_mask]) / ha6563_flux[good_mask]
        s2_lin = np.where(s2_lin > 0, s2_lin, np.nan)
        y = np.log10(s2_lin)
        s2_error = s2_error_propagation(sii6716_flux[good_mask], sii6716_flux_err[good_mask],
            sii6730_flux[good_mask], sii6730_flux_err[good_mask],
            ha6563_flux[good_mask], ha6563_flux_err[good_mask])

        indices = np.where(good_mask)
        REAL_ATOL = 1e-8
        
        if len(indices) == 1:
             iterator = indices[0]
             for idx_in_good, i in enumerate(iterator):
                y_val = y[idx_in_good]
                if not np.isfinite(y_val): continue
                roots = np.roots([c4, c3, c2, c1, (c0 - y_val)])
                realish = roots[np.abs(roots.imag) <= REAL_ATOL].real
                if realish.size == 0:
                    k = np.argmin(np.abs(roots.imag))
                    if np.abs(roots[k].imag) <= 1e-6: realish = np.array([roots[k].real])
                if realish.size == 0: continue
                in_range = realish[(realish >= -0.7) & (realish <= 0.3)]
                if in_range.size == 0: continue
                x_final = np.min(in_range)
                oh_s2_c20[i] = x_final + 8.69
                derivative_y = np.abs(c1 + 2*c2*x_final + 3*c3*x_final**2 + 4*c4*x_final**3)
                if derivative_y > 0:
                    derivative_x = 1.0 / derivative_y
                    obs_error = derivative_x * s2_error[idx_in_good]
                    oh_s2_c20_err[i] = np.sqrt(obs_error**2)
        else:
             for idx_in_good, idx in enumerate(zip(*indices)):
                y_val = y[idx_in_good]
                if not np.isfinite(y_val): continue
                roots = np.roots([c4, c3, c2, c1, (c0 - y_val)])
                realish = roots[np.abs(roots.imag) <= REAL_ATOL].real
                if realish.size == 0:
                    k = np.argmin(np.abs(roots.imag))
                    if np.abs(roots[k].imag) <= 1e-6: realish = np.array([roots[k].real])
                if realish.size == 0: continue
                in_range = realish[(realish >= -0.7) & (realish <= 0.3)]
                if in_range.size == 0: continue
                x_final = np.min(in_range)
                oh_s2_c20[idx] = x_final + 8.69
                derivative_y = np.abs(c1 + 2*c2*x_final + 3*c3*x_final**2 + 4*c4*x_final**3)
                if derivative_y > 0:
                    derivative_x = 1.0 / derivative_y
                    obs_error = derivative_x * s2_error[idx_in_good]
                    oh_s2_c20_err[idx] = np.sqrt(obs_error**2)

    combined_mask = good_mask & np.isfinite(oh_s2_c20)
    return oh_s2_c20, oh_s2_c20_err, combined_mask

# ------------------------------------------------------------------
# Main Execution
# ------------------------------------------------------------------

if __name__ == "__main__":
    # Find all FITS files in current directory
    fits_files = glob.glob("*.fits")
    # Exclude extended files
    fits_files = [f for f in fits_files if "_extended.fits" not in f]
    
    print(f"Found {len(fits_files)} FITS files to process.")
    
    for fits_file in fits_files:
        print(f"Processing {fits_file}...")
        
        try:
            with fits.open(fits_file) as hdul:
                # Check if PHASE3CATALOG extension exists
                if 'PHASE3CATALOG' not in hdul:
                    print(f"Skipping {fits_file}: PHASE3CATALOG extension not found.")
                    continue
                
                data = hdul['PHASE3CATALOG'].data
                header = hdul['PHASE3CATALOG'].header
                
                # Extract columns (using _CORR as requested)
                # Use .copy() to ensure we have numpy arrays and not FITS_rec objects which can be tricky
                try:
                    HB4861_FLUX = data['HB4861_FLUX_CORR'].copy()
                    HB4861_FLUX_ERR = data['HB4861_FLUX_CORR_ERR'].copy()
                    
                    HA6562_FLUX = data['HA6562_FLUX_CORR'].copy()
                    HA6562_FLUX_ERR = data['HA6562_FLUX_CORR_ERR'].copy()
                    
                    OIII5006_FLUX = data['OIII5006_FLUX_CORR'].copy()
                    OIII5006_FLUX_ERR = data['OIII5006_FLUX_CORR_ERR'].copy()
                    
                    NII6583_FLUX = data['NII6583_FLUX_CORR'].copy()
                    NII6583_FLUX_ERR = data['NII6583_FLUX_CORR_ERR'].copy()
                    
                    SII6716_FLUX = data['SII6716_FLUX_CORR'].copy()
                    SII6716_FLUX_ERR = data['SII6716_FLUX_CORR_ERR'].copy()
                    
                    SII6730_FLUX = data['SII6730_FLUX_CORR'].copy()
                    SII6730_FLUX_ERR = data['SII6730_FLUX_CORR_ERR'].copy()
                except KeyError as e:
                    print(f"Skipping {fits_file}: Missing required column {e}")
                    continue

                # Calculate Metallicities
                print("  Calculating metallicities...")
                
                # Dopita et al. (2016)
                y = np.log10(NII6583_FLUX / (SII6716_FLUX + SII6730_FLUX)) + 0.264*np.log10(NII6583_FLUX / HA6562_FLUX)
                O_H_D16 = 8.77 + y + 0.45*(y + 0.3)**5
                O_H_D16 = np.where((O_H_D16 < 7.63) | (O_H_D16 > 9.23), np.nan, O_H_D16)
                
                # Pilyugin & Grebel (2016)
                OIII_scaled = 1.33 * OIII5006_FLUX
                NII_scaled = 1.34 * NII6583_FLUX
                N2 = NII_scaled / HB4861_FLUX
                S2 = (SII6716_FLUX + SII6730_FLUX) / HB4861_FLUX
                R3 = OIII_scaled / HB4861_FLUX
                
                log_R3_S2 = np.log10(R3/S2)
                log_N2 = np.log10(N2)
                log_S2 = np.log10(S2)
                
                O_H_PG16 = np.full_like(log_N2, np.nan)
                upper_mask = log_N2 >= -0.6
                lower_mask = log_N2 < -0.6
                
                # Coefficients
                a1_u=8.424; a2_u=0.030; a3_u=0.751; a4_u=-0.349; a5_u=0.182; a6_u=0.508
                a1_l=8.072; a2_l=0.789; a3_l=0.726; a4_l=1.069; a5_l=-0.170; a6_l=0.022
                
                if np.any(upper_mask):
                    O_H_PG16[upper_mask] = (a1_u + a2_u * log_R3_S2[upper_mask] + a3_u * log_N2[upper_mask] + 
                                            (a4_u + a5_u * log_R3_S2[upper_mask] + a6_u * log_N2[upper_mask]) * log_S2[upper_mask])
                if np.any(lower_mask):
                    O_H_PG16[lower_mask] = (a1_l + a2_l * log_R3_S2[lower_mask] + a3_l * log_N2[lower_mask] + 
                                            (a4_l + a5_l * log_R3_S2[lower_mask] + a6_l * log_N2[lower_mask]) * log_S2[lower_mask])
                
                O_H_PG16 = np.where((O_H_PG16 < 7.63) | (O_H_PG16 > 9.23), np.nan, O_H_PG16)
                
                # N2S2-N06
                O_H_N2S2_N06, _ = calculate_n2s2_n06_metallicity(NII6583_FLUX, HA6562_FLUX, SII6716_FLUX, SII6730_FLUX)
                
                # O3N2-M13
                O_H_O3N2_M13, _ = calculate_o3n2_m13_metallicity(HB4861_FLUX, OIII5006_FLUX, NII6583_FLUX, HA6562_FLUX)
                O_H_O3N2_M13 = np.where((O_H_O3N2_M13 < 7.63) | (O_H_O3N2_M13 > 9.23), np.nan, O_H_O3N2_M13)
                
                # N2-M13
                O_H_N2_M13, _ = calculate_n2_m13_metallicity(NII6583_FLUX, HA6562_FLUX)
                O_H_N2_M13 = np.where((O_H_N2_M13 < 7.63) | (O_H_N2_M13 > 9.23), np.nan, O_H_N2_M13)
                
                # O3N2-PP04
                O_H_O3N2_PP04, _ = calculate_o3n2_pp04_metallicity(HB4861_FLUX, OIII5006_FLUX, NII6583_FLUX, HA6562_FLUX)
                O_H_O3N2_PP04 = np.where((O_H_O3N2_PP04 < 7.63) | (O_H_O3N2_PP04 > 9.23), np.nan, O_H_O3N2_PP04)
                
                # N2-PP04
                O_H_N2_PP04, _ = calculate_n2_pp04_metallicity(NII6583_FLUX, HA6562_FLUX)
                O_H_N2_PP04 = np.where((O_H_N2_PP04 < 7.63) | (O_H_N2_PP04 > 9.23), np.nan, O_H_N2_PP04)
                
                # C20 Suite
                O_H_O3N2_C20, O_H_O3N2_C20_ERR, _ = calculate_o3n2_c20_metallicity(
                    HB4861_FLUX, OIII5006_FLUX, NII6583_FLUX, HA6562_FLUX,
                    HB4861_FLUX_ERR, OIII5006_FLUX_ERR, NII6583_FLUX_ERR, HA6562_FLUX_ERR)
                O_H_O3N2_C20 = np.where((O_H_O3N2_C20 < 7.63) | (O_H_O3N2_C20 > 9.23), np.nan, O_H_O3N2_C20)
                
                O_H_O3S2_C20, O_H_O3S2_C20_ERR, _ = calculate_o3s2_c20_metallicity(
                    HB4861_FLUX, OIII5006_FLUX, SII6716_FLUX, SII6730_FLUX,
                    HB4861_FLUX_ERR, OIII5006_FLUX_ERR, SII6716_FLUX_ERR, SII6730_FLUX_ERR)
                O_H_O3S2_C20 = np.where((O_H_O3S2_C20 < 7.63) | (O_H_O3S2_C20 > 9.23), np.nan, O_H_O3S2_C20)
                
                O_H_RS32_C20, O_H_RS32_C20_ERR, _ = calculate_rs32_c20_metallicity(
                    HB4861_FLUX, HA6562_FLUX, OIII5006_FLUX, SII6716_FLUX, SII6730_FLUX,
                    HB4861_FLUX_ERR, HA6562_FLUX_ERR, OIII5006_FLUX_ERR, SII6716_FLUX_ERR, SII6730_FLUX_ERR)
                O_H_RS32_C20 = np.where((O_H_RS32_C20 < 7.63) | (O_H_RS32_C20 > 9.23), np.nan, O_H_RS32_C20)
                
                O_H_R3_C20, O_H_R3_C20_ERR, _ = calculate_r3_c20_metallicity(
                    HB4861_FLUX, HB4861_FLUX_ERR, OIII5006_FLUX, OIII5006_FLUX_ERR)
                O_H_R3_C20 = np.where((O_H_R3_C20 < 7.63) | (O_H_R3_C20 > 9.23), np.nan, O_H_R3_C20)
                
                O_H_N2_C20, O_H_N2_C20_ERR, _ = calculate_n2_c20_metallicity(
                    HA6562_FLUX, HA6562_FLUX_ERR, NII6583_FLUX, NII6583_FLUX_ERR)
                O_H_N2_C20 = np.where((O_H_N2_C20 < 7.63) | (O_H_N2_C20 > 9.23), np.nan, O_H_N2_C20)
                
                O_H_S2_C20, O_H_S2_C20_ERR, _ = calculate_s2_c20_metallicity(
                    HA6562_FLUX, HA6562_FLUX_ERR, SII6716_FLUX, SII6716_FLUX_ERR, SII6730_FLUX, SII6730_FLUX_ERR)
                O_H_S2_C20 = np.where((O_H_S2_C20 < 7.63) | (O_H_S2_C20 > 9.23), np.nan, O_H_S2_C20)
                
                # Combined C20
                # Stack arrays
                all_metallicities = np.stack([O_H_O3N2_C20, O_H_O3S2_C20, O_H_RS32_C20, O_H_R3_C20, O_H_N2_C20, O_H_S2_C20], axis=0)
                all_errors = np.stack([O_H_O3N2_C20_ERR, O_H_O3S2_C20_ERR, O_H_RS32_C20_ERR, O_H_R3_C20_ERR, O_H_N2_C20_ERR, O_H_S2_C20_ERR], axis=0)
                
                O_H_COMBINED_C20 = np.full_like(HB4861_FLUX, np.nan)
                O_H_COMBINED_C20_ERR = np.full_like(HB4861_FLUX, np.nan)
                
                # Find min error index along axis 0
                # Handle all-NaN slices by filling with inf temporarily
                temp_errors = all_errors.copy()
                temp_errors[~np.isfinite(temp_errors)] = np.inf
                min_error_idx = np.argmin(temp_errors, axis=0)
                    
                # Create a mask for valid entries (where at least one method worked)
                valid_mask = np.any(np.isfinite(all_errors), axis=0)
                
                if np.any(valid_mask):
                    # Use advanced indexing to select the best values
                    # We need indices for the other dimensions
                    if all_metallicities.ndim == 2: # (6, N)
                        idx_0 = min_error_idx[valid_mask]
                        idx_1 = np.where(valid_mask)[0]
                        O_H_COMBINED_C20[valid_mask] = all_metallicities[idx_0, idx_1]
                        O_H_COMBINED_C20_ERR[valid_mask] = all_errors[idx_0, idx_1]
                    else:
                        # Fallback for generic shape, though we expect 1D
                        # Flatten and reshape if needed, or just loop (slow but safe)
                        # Since we are 1D, the above should work.
                        pass

                # Create new columns
                new_cols = []
                new_cols.append(fits.Column(name='O_H_D16', format='D', array=O_H_D16))
                new_cols.append(fits.Column(name='O_H_PG16', format='D', array=O_H_PG16))
                new_cols.append(fits.Column(name='O_H_N2S2_N06', format='D', array=O_H_N2S2_N06))
                new_cols.append(fits.Column(name='O_H_O3N2_M13', format='D', array=O_H_O3N2_M13))
                new_cols.append(fits.Column(name='O_H_N2_M13', format='D', array=O_H_N2_M13))
                new_cols.append(fits.Column(name='O_H_O3N2_PP04', format='D', array=O_H_O3N2_PP04))
                new_cols.append(fits.Column(name='O_H_N2_PP04', format='D', array=O_H_N2_PP04))
                
                new_cols.append(fits.Column(name='O_H_O3N2_C20', format='D', array=O_H_O3N2_C20))
                new_cols.append(fits.Column(name='O_H_O3N2_C20_ERR', format='D', array=O_H_O3N2_C20_ERR))
                
                new_cols.append(fits.Column(name='O_H_O3S2_C20', format='D', array=O_H_O3S2_C20))
                new_cols.append(fits.Column(name='O_H_O3S2_C20_ERR', format='D', array=O_H_O3S2_C20_ERR))
                
                new_cols.append(fits.Column(name='O_H_RS32_C20', format='D', array=O_H_RS32_C20))
                new_cols.append(fits.Column(name='O_H_RS32_C20_ERR', format='D', array=O_H_RS32_C20_ERR))
                
                new_cols.append(fits.Column(name='O_H_R3_C20', format='D', array=O_H_R3_C20))
                new_cols.append(fits.Column(name='O_H_R3_C20_ERR', format='D', array=O_H_R3_C20_ERR))
                
                new_cols.append(fits.Column(name='O_H_N2_C20', format='D', array=O_H_N2_C20))
                new_cols.append(fits.Column(name='O_H_N2_C20_ERR', format='D', array=O_H_N2_C20_ERR))
                
                new_cols.append(fits.Column(name='O_H_S2_C20', format='D', array=O_H_S2_C20))
                new_cols.append(fits.Column(name='O_H_S2_C20_ERR', format='D', array=O_H_S2_C20_ERR))
                
                new_cols.append(fits.Column(name='O_H_COMBINED_C20', format='D', array=O_H_COMBINED_C20))
                new_cols.append(fits.Column(name='O_H_COMBINED_C20_ERR', format='D', array=O_H_COMBINED_C20_ERR))
                
                # Create new HDU
                orig_cols = hdul['PHASE3CATALOG'].columns
                hdu_new = fits.BinTableHDU.from_columns(orig_cols + fits.ColDefs(new_cols))
                hdu_new.name = 'PHASE3CATALOG'
                
                # Save to new file
                out_name = fits_file.replace(".fits", "_extended.fits")
                print(f"  Saving to {out_name}...")
                hdu_new.writeto(out_name, overwrite=True)
                
        except Exception as e:
            print(f"Error processing {fits_file}: {e}")
            import traceback
            traceback.print_exc()
