#!/usr/bin/env python
"""
nGIST-Compatible Spatial Masking Tool (v3tk)

This script generates binary spatial masks (FITS) and diagnostic overlays (PNG) for
MUSE galaxy data cubes. It is designed as a pre-processing step for nGIST (or other
pipelines), automatically masking foreground stars and (conservatively) background
galaxies while preserving the target galaxy emission.

Usage:
    python create_masks.py [pattern ...]

    If one or more [pattern] arguments are provided (e.g., "NGC*.fits"), it processes
    all matching files (deduplicated, order-preserving).
    Otherwise, it defaults to processing all "*_DATACUBE*_VRI.fits" files in the current
    directory.

--------------------------------------------------------------------------------
I. FILE INPUTS & OUTPUTS
--------------------------------------------------------------------------------
Inputs (per galaxy XXX):
  1. XXX_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk_VRI.fits  (Required)
     - 2D image or 3D cube (collapsed to 2D via nan-median) used for WCS and pixel geometry.

  2. XXX_combined_VRI.png                               (Optional)
     - High-res visual reference used ONLY if it matches the FITS pixel grid exactly
       (same nx, ny). If the PNG dimensions differ, the overlay is rendered on the
       FITS background to guarantee alignment.

Outputs (per galaxy XXX):
  1. XXX_mask.fits
     - Binary mask with the same spatial WCS/dimensions as the input image.
     - 0 = Unmasked (target/sky), 1 = Masked (foreground star / background object).

  2. XXX_combined_VRI_mask.png
     - Diagnostic overlay (pixel-locked to the chosen background).
     - Green outlines = masked stars. Brown outlines = masked background galaxies.
     - Optional: a dashed contour of the target-galaxy footprint (blue) if enabled.
     - Optional: a dashed/colored contour of the MUSE FoV footprint (yellow) if enabled.

  3. Optional DS9 region files (Legacy Surveys PA validation; OFF by default):
     - XXX_legacy_PA_EofN.reg
     - XXX_legacy_phi_from_E.reg

Logging:
  - The script prints detailed progress to stdout. If you want a persistent log file,
    run e.g.:
        python create_masks.py ... 2>&1 | tee v3tk_masking.log

--------------------------------------------------------------------------------
II. MASKING ALGORITHM
--------------------------------------------------------------------------------
The masking is hierarchical and conservative by default:

0) Load data + WCS
   - Reads the first HDU with >=2D data; collapses 3D to 2D (nan-median along axis=0).
   - Builds a 2D celestial WCS (prefers data-HDU WCS; falls back to primary if needed).
   - Computes pixel scale (arcsec/pix) from the WCS and the FoV radius from the corners
     (plus a small padding of ~10").

0b) MUSE FoV footprint (recommended; enabled by default)
   - If fov_use_mask=True, reads the R_FLUX extension to define the true MUSE footprint:
       FoV = finite & non-zero pixels (with optional morphological cleanup).
   - All star/galaxy rasterization can be gated to this FoV mask to prevent masking
     outside the valid MUSE area.
   - A small edge-buffer option can keep partial edge objects from being over-masked.
   - The FoV contour can be drawn on the diagnostic overlay (fov_draw_contour=True).

1) FOREGROUND STARS (Gaia DR3)
   - Query: Gaia DR3 cone search over the full FoV via ADQL (astroquery.gaia).
   - Modes (Config.gaia_star_mode):
     a) "foreground" (default): mask ONLY sources that look like Milky Way stars via
        kinematics (significant parallax and/or proper motion). This reduces the risk
        of masking Virgo-distance compact sources that happen to appear in Gaia.
     b) "strict": mask all Gaia sources in the FoV, but require reasonable astrometric
        quality (RUWE/IPD/excess-noise thresholds).
     c) "loose": mask any Gaia detection in the FoV (most complete; highest risk of
        masking extragalactic compact sources).
   - Size model: power-law radius vs Gaia G magnitude, with a seeing floor (>= 1×FWHM),
     plus a bright-star boost for very bright stars, then padded by gaia_margin_arcsec:
        r = max(r_min, r_ref * 10^(-0.2*(G - G_ref))) capped at r_max
   - Implementation note: circles are rasterized directly into the output mask; objects
     fully outside the FITS footprint are skipped robustly; FoV gating is applied if enabled.

2) TARGET-GALAXY FOOTPRINT (optional but enabled by default as a veto)
     - Goal: avoid masking “background candidates” that fall inside a pragmatic target
         footprint (helps suppress midplane artifacts / HII-region contamination).
     - Built directly from the R_MAG extension (surface-brightness map):
         * footprint = finite R_MAG spaxels with R_MAG < target_mu_lim
     - No foreground-star exclusion, no component filtering, no dilation.
     - If draw_target_iso_contour=True, the footprint is drawn as a dashed contour
         on the overlay (all connected regions are outlined, including single spaxels).

3) BACKGROUND GALAXIES (layered strategy; evidence-based by default)

   A) Legacy Surveys DR9 (default first-pass; best morphology + photo-z)
      - Query mechanism: NOIRLab Data Lab TAP via pyvo (if available).
      - Tables:
          ls_dr9.tractor  (morphology/type + optional shape_r, shape_e1/e2)
          ls_dr9.photo_z  (photometric redshift summary)
        joined by ls_id.
      - Selection (conservative):
        * Reject PSF-type detections (type == "PSF").
        * Require z_phot_l95 > legacy_z_l95_min.
        * Require photo-z significance using (z_mean, z_l95, z_u95):
            sigma_z ~ (z_u95 - z_l95)/3.92
            z_snr  ~ (z_mean - z_cut)/sigma_z   (z_cut defaults to legacy_z_l95_min)
          and enforce z_snr >= legacy_z_snr_min (unless disabled).
        * Optional maximum allowed z-width (legacy_z_width_max > 0).
        * Optional rejection near “good” Gaia point sources (legacy_reject_if_near_gaia_arcsec).
        * Optional veto inside the target footprint (reject_bg_inside_target_footprint=True).
      - Shape / sizing:
        * Prefer Tractor intrinsic size shape_r (treated as angular), scaled by
          legacy_shape_r_scale, then floored at legacy_r_min_arcsec and capped
          at legacy_r_max_arcsec (or global cap).
        * If enabled (legacy_use_ellipses=True) and e1/e2 are present, derive axis ratio
          and a position angle; otherwise fall back to a circle.
        * If legacy_wcs_sample_ellipses=True, ellipses are rasterized by sampling in the
          local (east,north) tangent plane and mapping through the FITS WCS to avoid
          common PA/parity mistakes on rotated/flipped WCS grids.
        * If you prefer maximum robustness, force circles for all Legacy objects with
          legacy_force_circles=True.
      - If Legacy masking succeeds (i.e., at least one object is masked), the script
        skips all lower-fidelity fallback catalogs (PS1/SkyMapper/SDSS/NED).

   B) Evidence-based confirmation (SDSS spec-z / NED spec-z; used when Legacy did NOT succeed)
      - If require_nonvirgo_confirmation_for_galaxy_mask=True (default), PS1/SkyMapper/SDSS
        photometric candidates are masked ONLY when there is spectroscopic evidence that
        they are background:
          * cz <= virgo_keep_cz_max_kms       -> treated as Virgo/nearby (NOT masked)
          * cz >= background_mask_cz_min_kms  -> treated as background (masked)
        SDSS spectroscopy is preferred; NED “SPEC” redshifts are used as fallback.

      - Optional “no-spec” fallback (still conservative; OFF/ON via flags):
        * If distance is unknown, the script may still mask only objects that are VERY
          extended and bright (ps1_allow_photometric_fallback, ps1_fallback_ext_min,
          ps1_fallback_rmag_max). SDSS has a parallel (off-by-default) fallback.

   C) Pan-STARRS (PS1) / SkyMapper (photometric fallback; only if Legacy did not succeed)
      - PS1 queried via MAST if possible; falls back to VizieR (II/349/ps1) if needed.
      - Galaxy-like selection primarily via extendedness:
          (PSF mag - Kron mag) > ps1_ext_thresh, with a sanity cap ps1_ext_max.
      - Additional guards (when VizieR PS1 provides them):
        * Minimum detections (Nr >= ps1_min_Nr).
        * Qual bitmask filtering to require extendedness/quality and reject suspect stack objects
          (ps1_qual_* flags).
        * Optional strict color cuts (ps1_enable_color_cuts) to reduce contamination from compact
          blue sources in the target galaxy, at the cost of completeness.
      - Sizes:
        * PS1 VizieR lacks robust size/shape columns; when extendedness is satisfied, a small
          fallback radius gal_fallback_arcsec is used.
        * SkyMapper may provide a/b/PA; if present, ellipses are rasterized.

   D) SDSS photometry (supplemental; footprint-limited)
      - Used only if Legacy did not succeed and enable_sdss=True.
      - SDSS query_region is limited to <= 3 arcmin radius by the service.
      - Uses morphology-aware sizing:
        * point-like if (psfMag_r - modelMag_r) < sdss_pointlike_dmag_max -> star-like radius model
        * extended: prefer petroR90, then petroR50, then model radii (deV/exp), with configurable scales
      - High-z PSF override: if a confirmed spec-z exceeds highz_psf_override_zmin, use a seeing-based
        radius (highz_psf_k_fwhm × FWHM), capped by highz_psf_rmax_arcsec if enabled.
      - Optional SDSS photometric fallback exists but is OFF by default.

4) STAR FALLBACK (SDSS; only if Gaia yields zero stars)
   - If Gaia is unavailable/empty and sdss_star_fallback=True, SDSS “STAR” detections are
     masked using the same star-radius model (using an r-band magnitude proxy).

--------------------------------------------------------------------------------
III. KEY CONFIGURATION PARAMETERS (Config dataclass)
--------------------------------------------------------------------------------
[Global]
  fwhm_arcsec ............................: seeing FWHM used for floors/overrides
  gaia_margin_arcsec .....................: padding added to star/galaxy radii
  exclude_center_arcsec ..................: inner radius to force UNMASKED (default 0)

[FoV gating (R_FLUX footprint)]
  fov_use_mask ...........................: enable FoV gating from R_FLUX extension
  fov_extname ............................: extension name (default "R_FLUX")
  fov_close_size_pix .....................: cleanup kernel size for closing/filling
  fov_min_abs ............................: minimum absolute flux for valid FoV pixels
  fov_edge_star_buffer_arcsec ............: edge buffer to prevent over-masking outside FoV
  fov_draw_contour .......................: draw FoV contour on overlay
  fov_flip_y .............................: flip FoV contour for PNG background orientation

[Gaia Stars]
  gaia_star_mode .........................: "foreground" | "strict" | "loose"
  gaia_gmag_max ..........................: faint limit for Gaia query
  (Foreground-by-kinematics thresholds)
    gaia_parallax_snr_min, gaia_parallax_min_mas
    gaia_pm_snr_min, gaia_pm_min_masyr
  (Strict-mode quality thresholds)
    gaia_ruwe_max, gaia_ipd_frac_multi_peak_max, gaia_astrometric_excess_noise_sig_max
  (Star radius model)
    star_r_min_arcsec, star_r_ref_arcsec, star_g_ref, star_r_max_arcsec

[Target Footprint Veto / Contour]
  reject_bg_inside_target_footprint ......: veto masking inside target footprint
  target_mu_lim ..........................: isophote (mag/arcsec^2) for footprint
  target_mu_min_area_pix, target_mu_min_area_frac
  target_mu_close_radius_pix .............: morphological closing radius (pixels)
  target_iso_dilate_fwhm .................: dilation by ~FWHM (in pixel units internally)
  draw_target_iso_contour ................: draw dashed contour on overlay

[Legacy Surveys DR9]
  enable_legacy ..........................: enable DR9 TAP queries (pyvo required)
  legacy_tap_url .........................: TAP endpoint (default NOIRLab Data Lab)
  legacy_tractor_table, legacy_photoz_table
  legacy_z_l95_min .......................: photo-z lower 95% bound cut
  legacy_z_snr_min .......................: photo-z significance cut
  legacy_z_snr_use_threshold .............: use (z_mean - z_cut)/sigma_z
  legacy_z_width_max .....................: optional maximum (z_u95 - z_l95)
  legacy_use_ellipses ....................: use Tractor e1/e2 to form ellipses
  legacy_force_circles ...................: mask all Legacy objects as circles
  legacy_wcs_sample_ellipses .............: robust WCS-sampled ellipse rasterization
  legacy_pa_east_of_north ................: convention toggle for sampled ellipses
  legacy_pa_offset_deg ...................: fixed extra rotation (degrees)
  legacy_ellipse_npts ....................: polygon sampling resolution
  legacy_r_min_arcsec, legacy_r_max_arcsec
  legacy_shape_r_scale ...................: scale for Tractor shape_r
  legacy_seeing_scale ....................: fallback when shape_r missing/invalid
  legacy_reject_if_near_gaia_arcsec ......: reject Legacy detections near Gaia sources
  legacy_write_ds9_regions ...............: write DS9 region files for PA debugging

[Evidence-based masking (Spec-z)]
  require_nonvirgo_confirmation_for_galaxy_mask ..: if True, only mask with cz evidence
  virgo_keep_cz_max_kms ..................: keep (NOT mask) if cz <= this
  background_mask_cz_min_kms .............: mask if cz >= this
  sdss_spec_match_arcsec, ned_match_arcsec: crossmatch radii

[PS1 / SkyMapper Photometric Galaxies]
  ps1_ext_thresh, ps1_ext_max ............: PSF-Kron extendedness thresholds
  ps1_rmag_max ...........................: optional magnitude cut
  ps1_min_Nr .............................: minimum detections (VizieR PS1)
  ps1_qual_extended_required_bits ........: required Qual bits (VizieR PS1)
  ps1_require_qual_good / _primary_best ..: quality/primary constraints
  ps1_reject_qual_suspect ................: reject suspect/poor stack bits
  ps1_enable_color_cuts ..................: optional g-r / r-i cuts
  ps1_e_mag_max ..........................: max allowed mag error for color cuts
  ps1_allow_photometric_fallback .........: allow no-cz fallback for very extended/bright
  ps1_fallback_ext_min, ps1_fallback_rmag_max
  ps1_reject_if_near_gaia_arcsec .........: avoid masking Gaia-like sources as galaxies
  gal_fallback_arcsec ....................: radius when no size is available
  gal_r_min_arcsec, gal_r_max_arcsec .....: global radius floor/cap

[SDSS Supplemental]
  enable_sdss ............................: enable SDSS photometry/spectroscopy
  sdss_data_release ......................: SDSS DR to query
  sdss_pointlike_dmag_max ................: psf-model threshold for point-like
  sdss_petro_r90_scale, sdss_petro_r50_scale, sdss_model_radius_scale
  highz_psf_override_zmin ................: if spec-z > this, use seeing-based radius
  highz_psf_k_fwhm, highz_psf_rmax_arcsec
  sdss_reject_if_near_gaia_arcsec ........: Gaia-near rejection
  sdss_star_fallback .....................: use SDSS stars only if Gaia returns none
  sdss_allow_photometric_fallback ........: SDSS no-cz fallback (OFF by default)

[Output]
  use_png_background .....................: use *_combined_VRI.png ONLY if pixel-matched
  output_dpi .............................: DPI for overlays (pixel-locked sizing)
  log_each_star / log_each_galaxy ........: per-object logging toggles
  log_max_galaxies .......................: cap per-object galaxy logs (<=0 unlimited)

--------------------------------------------------------------------------------
IV. REQUIREMENTS
--------------------------------------------------------------------------------
Core:
  pip install numpy astropy matplotlib

Optional (for PNG background):
  pip install pillow

Catalog services (Gaia/PS1/VizieR/NED/SDSS):
  pip install astroquery

Legacy DR9 TAP (recommended for best background masking):
  pip install pyvo

FoV/target-footprint morphology (recommended; otherwise these steps are skipped):
  pip install scipy
"""

from __future__ import annotations

import glob
import os
import sys
import traceback
from dataclasses import dataclass
import numpy as np
from datetime import datetime
from time import perf_counter

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs.utils import proj_plane_pixel_scales

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, Polygon
from matplotlib.path import Path

try:
    from PIL import Image
except Exception:
    Image = None

# Astroquery imports (kept inside try so the script can still run without some services)
try:
    from astroquery.gaia import Gaia
except Exception:
    Gaia = None

try:
    from astroquery.mast import Catalogs
except Exception:
    Catalogs = None

try:
    from astroquery.vizier import Vizier
except Exception:
    Vizier = None

try:
    from astroquery.ipac.ned import Ned
except Exception:
    try:
        from astroquery.ned import Ned
    except Exception:
        Ned = None

try:
    from astroquery.sdss import SDSS
except Exception:
    SDSS = None

# Optional: NOIRLab Astro Data Lab TAP queries (Legacy Surveys DR9)
try:
    import pyvo
except Exception:
    pyvo = None


@dataclass
class Config:
    fwhm_arcsec: float = 1.0                # typical MUSE seeing FWHM, tune if needed
    gaia_gmag_max: float = 21.0             # ignore very faint Gaia sources
    gaia_margin_arcsec: float = 1.0         # extra padding on radii (registration / wings)

    # MUSE FoV detection from R_FLUX (non-NaN region)
    fov_use_mask: bool = True
    fov_extname: str = "R_FLUX"
    fov_close_size_pix: int = 5
    fov_min_abs: float = 0.0
    fov_edge_star_buffer_arcsec: float = 5.0
    fov_draw_contour: bool = True
    fov_contour_color: str = "yellow"
    fov_contour_linestyle: str = ":"
    fov_contour_linewidth: float = 1.0
    fov_flip_y: bool = True

    # Foreground-star selection in Gaia:
    # - "loose": mask any Gaia detection within the FOV (most complete; can mask some non-foreground knots)
    # - "strict": like loose but also requires good astrometric-quality metrics (cleaner; can miss some stars)
    # - "foreground": mask only sources that look like Milky Way stars via *kinematics*
    #                 (significant parallax and/or proper motion). This is the closest practical
    #                 proxy to "not at Virgo distance" because Gaia cannot measure 16.5 Mpc parallaxes.
    gaia_star_mode: str = "foreground"

    # Gaia quality thresholds used in "strict" mode
    gaia_ruwe_max: float = 1.4
    gaia_ipd_frac_multi_peak_max: float = 2.0
    gaia_astrometric_excess_noise_sig_max: float = 2.0

    # Gaia kinematics thresholds used in "foreground" mode
    # Require either significant positive parallax OR significant proper motion.
    gaia_parallax_snr_min: float = 5.0
    gaia_parallax_min_mas: float = 0.002
    gaia_pm_snr_min: float = 5.0
    gaia_pm_min_masyr: float = 0.02

    # NOTE: We intentionally do NOT exclude the inner galaxy by default, because
    # you may still want to mask true foreground stars/background galaxies there.
    # Keep this at 0 unless you explicitly want an inner "no-mask" zone.
    exclude_center_arcsec: float = 0.0

    # PS1(VizieR) quality cuts (to avoid spurious detections / bad photometry).
    # Defaults are "looser but still sane" (tune as needed).
    ps1_min_Nr: int = 2
    ps1_e_mag_max: float = 0.20

    # PS1(VizieR) qualityFlag (Qual) bitmask filtering (see II/349 ReadMe, Note 3)
    #  1: extended in our data
    #  2: extended in external data
    #  4: good-quality measurement in our data
    # 16: good-quality object in the stack (>1 good stack measurement)
    # 64: suspect object in the stack
    #128: poor-quality stack object
    # For very strict selection, require bit 1 (extended in PS1) and do not accept
    # bit-2-only (extended only in external catalogs).
    ps1_qual_extended_required_bits: int = 1
    ps1_require_qual_good: bool = True
    ps1_reject_qual_suspect: bool = True
    ps1_require_qual_primary_best: bool = True

    # Optional very-strict color cuts to reduce contamination from compact blue
    # sources in the target galaxy (HII regions / some PNe). This will also drop
    # some real blue background galaxies.
    ps1_enable_color_cuts: bool = True
    ps1_g_r_min: float = 0.2
    ps1_r_i_min: float = 0.0

    # VizieR PS1 table (II/349/ps1) does not provide robust size/shape columns.
    # When an object passes the *extendedness* test (PSF - Kron), we apply this
    # small fallback radius.
    gal_fallback_arcsec: float = 3.0

    # star mask radius model (in arcsec) as function of Gaia G magnitude:
    # r = max(r_min, r_ref * 10^(-0.2*(G-G_ref))) capped at r_max
    star_r_min_arcsec: float = 1.5
    star_r_ref_arcsec: float = 5.0
    star_g_ref: float = 15.0
    star_r_max_arcsec: float = 25.0

    # galaxy-like selection (Pan-STARRS): extended if (PSF - Kron) > threshold
    ps1_ext_thresh: float = 0.25
    # Upper bound to guard against pathological photometry (blends/saturation can
    # produce huge PSF-Kron differences that are not real galaxies).
    ps1_ext_max: float = 1.5
    ps1_rmag_max: float = 22.0              # ignore very faint objects (optional)
    ps1_require_ri_extended: bool = False    # if i-band mags exist, require extendedness in both r and i
    gal_r_min_arcsec: float = 2.0
    gal_r_max_arcsec: float = 30.0

    # Photometric fallback when no cz/z evidence exists (PS1/SkyMapper)
    ps1_allow_photometric_fallback: bool = True
    ps1_fallback_ext_min: float = 0.6
    ps1_fallback_rmag_max: float = 21.0

    # Reject PS1 objects that coincide with Gaia sources (usually stars / blends)
    # so we don't double-count stars as "background galaxies".
    ps1_reject_if_near_gaia_arcsec: float = 0.8

    # Virgo-distance veto for *galaxy-like* objects (when a redshift/velocity is known).
    # If enabled and a candidate matches a catalog object consistent with Virgo distance,
    # we do NOT mask it.
    enable_virgo_distance_veto: bool = True
    virgo_distance_mpc: float = 16.5
    virgo_distance_tolerance_mpc: float = 5.0
    hubble_km_s_mpc: float = 70.0
    virgo_match_arcsec: float = 1.0

    # Logging controls (galaxy fields can be huge; per-object prints can look "stuck")
    log_each_star: bool = True
    log_each_galaxy: bool = True
    # <=0 means unlimited
    log_max_galaxies: int = 0

    # Target-galaxy footprint veto (pragmatic guard against midplane artifacts)
    reject_bg_inside_target_footprint: bool = True
    target_mu_lim: float = 26.0             # R-band surface brightness isophote (mag/arcsec^2)
    target_mu_min_area_pix: int = 500       # minimum area to keep a component
    target_mu_min_area_frac: float = 0.002  # minimum area as a fraction of the FoV
    target_mu_close_radius_pix: int = 3     # morphological closing radius (pixels)
    target_iso_dilate_fwhm: float = 1.0     # dilate footprint by ~FWHM
    draw_target_iso_contour: bool = True

    use_png_background: bool = True         # if False, uses FITS as background for overlay
    output_dpi: int = 200

    # SDSS integration (supplemental; footprint-limited)
    enable_sdss: bool = True
    sdss_data_release: int = 18
    log_sdss_colnames: bool = True
    # Legacy (kept for compatibility): older versions used petroRad_* × scale.
    # The current SDSS masking uses morphology-aware sizing (psfMag-modelMag) and
    # prefers petroR90/petroR50/model radii, so this is typically unused.
    sdss_petro_scale: float = 2.5
    # SDSS morphology proxy: point-like if (psfMag_r - modelMag_r) < threshold
    sdss_pointlike_dmag_max: float = 0.145
    # SDSS extended-object sizing factors (arcsec)
    sdss_petro_r90_scale: float = 1.2
    sdss_petro_r50_scale: float = 1.8
    sdss_model_radius_scale: float = 3.0

    # Photometric fallback when no cz/z evidence exists (SDSS)
    sdss_allow_photometric_fallback: bool = False
    sdss_fallback_dmag_min: float = 0.30
    sdss_fallback_petroR90_min_arcsec: float = 2.0
    sdss_fallback_petroR90_max_arcsec: float = 8.0
    sdss_fallback_rmag_max: float = 21.0

    # If a candidate has a confirmed high redshift, it is typically PSF-limited at MUSE resolution.
    # In that case, ignore SDSS size proxies and use a seeing-based radius.
    highz_psf_override_zmin: float = 0.1
    highz_psf_k_fwhm: float = 1.5
    # <=0 disables the cap
    highz_psf_rmax_arcsec: float = 3.0
    # Reject SDSS "galaxies" near good Gaia sources to avoid stellar contamination
    sdss_reject_if_near_gaia_arcsec: float = 0.8
    # Only use SDSS stars if Gaia returns no stars (fallback classifier)
    sdss_star_fallback: bool = True

    # === Legacy Surveys DR9 background-galaxy masking (default first-pass) ===
    enable_legacy: bool = True

    # Photo-z gating: require lower 95% bound above this redshift
    legacy_z_l95_min: float = 0.01

    # Photo-z significance control (Legacy DR9)
    legacy_z_snr_min: float = 5.0
    legacy_z_snr_use_threshold: bool = True
    legacy_z_width_max: float = 0.0

    # Masking size model
    # Legacy sizing policy:
    # - Prefer Tractor intrinsic size (shape_r) when available.
    # - Enforce a minimum (legacy_r_min_arcsec).
    # - Only fall back to seeing when shape_r is missing/invalid.
    legacy_use_ellipses: bool = True        # use shape_e1/e2 when available
    # If True, ignore PA/axis ratio and mask Legacy detections as circles.
    # This is the most robust option if you want to avoid all PA/parity conventions.
    legacy_force_circles: bool = False
    # If True, rasterize Tractor ellipses by sampling in the local sky tangent plane
    # (east/north offsets) and transforming those points through the cube WCS.
    # This is robust to WCS rotation/parity and prevents PA mismatches when the
    # cube pixel axes are rotated relative to north/east.
    legacy_wcs_sample_ellipses: bool = True
    # Optional convention tweak:
    # - Default assumes the (e1,e2)->angle half-angle is measured from +east toward +north.
    # - If your reference overlay expects PA measured East of North, set this True.
    legacy_pa_east_of_north: bool = True
    # Additional fixed rotation to apply (degrees) after any convention conversion.
    legacy_pa_offset_deg: float = 0.0
    # Number of points used to approximate an ellipse polygon when sampling.
    legacy_ellipse_npts: int = 96

    # DS9 region export (for validating PA convention).
    # Writes FK5 regions with ellipse(ra,dec,a",b",PA) where PA is degrees East of North.
    legacy_write_ds9_regions: bool = False
    # If True, include circular Legacy objects as DS9 circles as well.
    legacy_ds9_include_circles: bool = False
    legacy_r_min_arcsec: float = 1.0        # floor on Legacy semi-major/minor axes
    legacy_shape_r_scale: float = 1.0       # "just use it" by default (shape_r is already angular)
    legacy_reject_if_near_gaia_arcsec: float = 0.8
    # Optional separate Legacy cap (independent from global gal_r_max_arcsec)
    legacy_r_max_arcsec: float = 15.0
    # Only used when shape_r is missing/invalid
    legacy_seeing_scale: float = 1.0        # multiply fwhm_arcsec by this

    # Data Lab TAP settings (table names may vary; keep configurable)
    legacy_tap_url: str = "https://datalab.noirlab.edu/tap"
    legacy_tractor_table: str = "ls_dr9.tractor"
    legacy_photoz_table: str = "ls_dr9.photo_z"

    # === Evidence-based galaxy masking ===
    # Only mask galaxies when spectroscopic redshift confirms they are NOT at Virgo distance
    require_nonvirgo_confirmation_for_galaxy_mask: bool = True
    # Treat as "possible Virgo/nearby" (do NOT mask) if cz <= this
    virgo_keep_cz_max_kms: float = 3500.0
    # Treat as "definitely background" (mask) if cz >= this
    background_mask_cz_min_kms: float = 5000.0
    # Matching radii for redshift catalogs
    sdss_spec_match_arcsec: float = 1.0
    ned_match_arcsec: float = 1.0


def safe_base_id(rfits_path: str) -> str:
    # From XXX_DATACUBE..._VRI.fits => XXX
    bn = os.path.basename(rfits_path)
    return bn.split("_DATACUBE")[0]


def format_radec_hmsdms(sc: SkyCoord, precision: int = 2) -> tuple[str, str]:
    """Return RA/Dec strings in sexagesimal format (NED-friendly).

    RA:  hh:mm:ss.ss
    Dec: ±dd:mm:ss.ss
    """
    ra_hms = sc.ra.to_string(unit=u.hour, sep=":", precision=precision, pad=True)
    dec_dms = sc.dec.to_string(unit=u.deg, sep=":", precision=precision, pad=True, alwayssign=True)
    return ra_hms, dec_dms


def iau_coord_name(prefix: str, sc: SkyCoord, ra_precision: int = 2, dec_precision: int = 1) -> str:
    """Return an IAU-style coordinate name like: SDSS J122544.87+123947.6"""
    ra = sc.ra.to_string(unit=u.hour, sep="", precision=ra_precision, pad=True)
    dec = sc.dec.to_string(unit=u.deg, sep="", precision=dec_precision, pad=True, alwayssign=True)
    return f"{prefix} J{ra}{dec}"


def load_r_image_and_wcs(rfits_path: str):
    with fits.open(rfits_path) as hdul:
        # Some of these products store WCS in the primary header (HDU0) but the
        # actual image/cube data in an extension, so hdul[0].data can be None.
        primary_hdr = hdul[0].header

        data = None
        data_hdr = None
        for hdu in hdul:
            if hdu.data is None:
                continue
            if getattr(hdu.data, "ndim", 0) >= 2:
                data = hdu.data
                data_hdr = hdu.header
                break

        if data is None:
            raise ValueError(f"No 2D/3D image data found in {rfits_path}")

    # allow for accidental 3D arrays (collapse if needed)
    if data.ndim == 3:
        data2d = np.nanmedian(data, axis=0)
    elif data.ndim == 2:
        data2d = data
    else:
        raise ValueError(f"Unexpected data ndim={data.ndim} in {rfits_path}")

    # Prefer WCS from the primary header if it is valid; otherwise fall back to
    # the data HDU header. Force a 2D WCS (spatial/celestial) regardless of any
    # higher-dimensional keywords.
    hdr_for_wcs = data_hdr
    try:
        w0 = WCS(data_hdr, naxis=2)
        if not w0.has_celestial:
            raise ValueError("Data header WCS has no celestial component")
        w = w0
    except Exception:
        hdr_for_wcs = primary_hdr
        w = WCS(hdr_for_wcs, naxis=2)

    ny, nx = data2d.shape
    return data2d, w, hdr_for_wcs, nx, ny


def build_muse_fov_mask(
    rfits_path: str,
    fov_extname: str,
    closing_size: int = 5,
    min_abs: float = 0.0,
) -> np.ndarray | None:
    """Return boolean FoV mask from the R_FLUX HDU (finite & non-zero pixels), with optional cleanup."""
    try:
        with fits.open(rfits_path) as hdul:
            fov_hdu = None
            for hdu in hdul:
                if getattr(hdu, "name", "").strip().upper() == str(fov_extname).strip().upper():
                    fov_hdu = hdu
                    break
            if fov_hdu is None:
                print(f"WARNING: FoV extension '{fov_extname}' not found in {rfits_path}; FoV gating disabled.")
                return None
            data = fov_hdu.data
            if data is None:
                print(f"WARNING: FoV extension '{fov_extname}' has no data; FoV gating disabled.")
                return None
    except Exception as e:
        print(f"WARNING: failed to read FoV extension '{fov_extname}' from {rfits_path}: {e}")
        return None

    if data.ndim == 3:
        img = np.nanmedian(data, axis=0)
    elif data.ndim == 2:
        img = data
    else:
        print(f"WARNING: FoV extension '{fov_extname}' has ndim={data.ndim}; FoV gating disabled.")
        return None

    try:
        min_abs_val = float(min_abs)
    except Exception:
        min_abs_val = 0.0
    if min_abs_val < 0:
        min_abs_val = 0.0
    fov = np.isfinite(img) & (np.abs(img) > min_abs_val)

    if closing_size and int(closing_size) > 1:
        try:
            from scipy.ndimage import binary_closing, binary_fill_holes

            st = np.ones((int(closing_size), int(closing_size)), dtype=bool)
            fov = binary_closing(fov, structure=st)
            fov = binary_fill_holes(fov)
        except Exception:
            print("WARNING: scipy.ndimage unavailable; FoV mask cleanup skipped.")

    return np.asarray(fov, dtype=bool)


def pixel_scale_arcsec(w: WCS) -> float:
    # robust pixel scale estimate (arcsec/pix)
    # for 2D WCS, scales[0], scales[1] in deg/pix
    scales = proj_plane_pixel_scales(w) * u.deg
    return float(np.mean(scales.to_value(u.arcsec)))


def fov_center_and_radius(w: WCS, nx: int, ny: int):
    # Center: use the true pixel center for stability (avoid RA wrap / corner-mean issues)
    cx = (nx - 1) / 2.0
    cy = (ny - 1) / 2.0
    center = w.pixel_to_world(cx, cy)

    # Radius: max separation to the four corners
    pix = np.array([[0, 0], [nx - 1, 0], [nx - 1, ny - 1], [0, ny - 1]], dtype=float)
    sky = w.pixel_to_world(pix[:, 0], pix[:, 1])
    radius = center.separation(sky).max()
    return center, radius


def star_radius_arcsec_from_g(cfg: Config, gmag: float) -> float:
    # sanitize
    try:
        g = float(gmag)
    except Exception:
        g = float("nan")
    if (not np.isfinite(g)) or (g < -5.0) or (g > 40.0):
        g = 18.0  # safe fallback

    exp = -0.2 * (g - float(cfg.star_g_ref))
    exp = float(np.clip(exp, -10.0, 10.0))  # hard overflow guard

    r = float(cfg.star_r_ref_arcsec) * (10.0 ** exp)
    r = max(float(cfg.star_r_min_arcsec), min(float(cfg.star_r_max_arcsec), float(r)))
    r = max(float(r), 1.0 * float(cfg.fwhm_arcsec))

    # --- NEW: bright-star boost ---
    if g < 10.0:
        r *= 1.5
    elif g < 14.0:
        r *= 1.25

    r += float(cfg.gaia_margin_arcsec)
    return float(r)


def query_gaia_sources(center: SkyCoord, radius: u.Quantity, cfg: Config):
        if Gaia is None:
                print("WARNING: astroquery.gaia not available; skipping Gaia query.")
                return None

        # Use ADQL to keep it stable across astroquery versions
        # radius is in deg for CIRCLE
        rad_deg = radius.to(u.deg).value

        mode = str(getattr(cfg, "gaia_star_mode", "strict")).lower().strip()

        where_quality = ""
        if mode == "strict":
            where_quality = f"""
                AND (ruwe IS NULL OR ruwe < {float(cfg.gaia_ruwe_max)})
                AND (ipd_frac_multi_peak IS NULL OR ipd_frac_multi_peak <= {float(cfg.gaia_ipd_frac_multi_peak_max)})
                AND (astrometric_excess_noise_sig IS NULL OR astrometric_excess_noise_sig <= {float(cfg.gaia_astrometric_excess_noise_sig_max)})
                """

        query = f"""
            SELECT
              source_id, ra, dec,
              phot_g_mean_mag,
              ruwe,
              ipd_frac_multi_peak,
              astrometric_excess_noise_sig,
              parallax, parallax_error,
              pmra, pmra_error,
              pmdec, pmdec_error
            FROM gaiadr3.gaia_source
            WHERE 1=CONTAINS(
              POINT('ICRS', ra, dec),
              CIRCLE('ICRS', {center.ra.deg}, {center.dec.deg}, {rad_deg})
            )
            AND phot_g_mean_mag IS NOT NULL
            AND phot_g_mean_mag < {cfg.gaia_gmag_max}
            {where_quality}
        """

        import time
        max_retries = 3
        for attempt in range(max_retries):
            try:
                job = Gaia.launch_job_async(query, dump_to_file=False)
                return job.get_results()
            except Exception as e:
                # Catch general exceptions including requests.exceptions.HTTPError
                error_msg = str(e)
                if "500" in error_msg or "503" in error_msg or "504" in error_msg or "timeout" in error_msg.lower():
                    if attempt < max_retries - 1:
                        print(f"WARNING: Gaia query failed with server error ({error_msg}). Retrying {attempt+1}/{max_retries}...")
                        time.sleep(2 * (attempt + 1))
                        continue
                print(f"WARNING: Gaia query failed: {e}")
                return None


def _gaia_row_is_foreground_by_kinematics(row, cfg: Config) -> bool:
    """Best-effort foreground-star flag using Gaia parallax and proper motion.

    Gaia cannot measure Virgo distances (parallax at 16.5 Mpc is ~0.00006 mas),
    so we instead tag Milky Way stars via significant parallax / proper motion.
    """

    def _get_float(name: str) -> float:
        try:
            v = row[name]
            # Astroquery/astropy tables may use masked scalars.
            if v is None or getattr(v, "mask", False):
                return float("nan")
            return float(v)
        except Exception:
            return float("nan")

    plx = _get_float("parallax")
    e_plx = _get_float("parallax_error")
    pmra = _get_float("pmra")
    e_pmra = _get_float("pmra_error")
    pmdec = _get_float("pmdec")
    e_pmdec = _get_float("pmdec_error")

    parallax_snr = plx / e_plx if (np.isfinite(plx) and np.isfinite(e_plx) and e_plx > 0) else float("nan")
    pm = np.hypot(pmra, pmdec) if (np.isfinite(pmra) and np.isfinite(pmdec)) else float("nan")
    pm_err = np.hypot(e_pmra, e_pmdec) if (np.isfinite(e_pmra) and np.isfinite(e_pmdec)) else float("nan")
    pm_snr = pm / pm_err if (np.isfinite(pm) and np.isfinite(pm_err) and pm_err > 0) else float("nan")

    is_fg_parallax = (
        np.isfinite(parallax_snr)
        and parallax_snr >= float(cfg.gaia_parallax_snr_min)
        and np.isfinite(plx)
        and plx >= float(cfg.gaia_parallax_min_mas)
    )
    is_fg_pm = (
        np.isfinite(pm_snr)
        and pm_snr >= float(cfg.gaia_pm_snr_min)
        and np.isfinite(pm)
        and pm >= float(cfg.gaia_pm_min_masyr)
    )

    return bool(is_fg_parallax or is_fg_pm)


def gaia_foreground_reason(row, cfg: Config) -> str:
    """Return a human-readable explanation of the Gaia foreground selection."""

    def gf(name: str) -> float:
        try:
            v = row[name]
            if v is None or getattr(v, "mask", False):
                return float("nan")
            return float(v)
        except Exception:
            return float("nan")

    plx = gf("parallax")
    eplx = gf("parallax_error")
    pmra = gf("pmra")
    epmra = gf("pmra_error")
    pmdec = gf("pmdec")
    epmdec = gf("pmdec_error")

    plx_snr = plx / eplx if (np.isfinite(plx) and np.isfinite(eplx) and eplx > 0) else float("nan")
    pm = np.hypot(pmra, pmdec) if (np.isfinite(pmra) and np.isfinite(pmdec)) else float("nan")
    pm_err = np.hypot(epmra, epmdec) if (np.isfinite(epmra) and np.isfinite(epmdec)) else float("nan")
    pm_snr = pm / pm_err if (np.isfinite(pm) and np.isfinite(pm_err) and pm_err > 0) else float("nan")

    is_fg_parallax = (
        np.isfinite(plx_snr)
        and plx_snr >= float(cfg.gaia_parallax_snr_min)
        and np.isfinite(plx)
        and plx >= float(cfg.gaia_parallax_min_mas)
    )
    is_fg_pm = (
        np.isfinite(pm_snr)
        and pm_snr >= float(cfg.gaia_pm_snr_min)
        and np.isfinite(pm)
        and pm >= float(cfg.gaia_pm_min_masyr)
    )

    parts = []
    parts.append(f"plx={plx:.3f}±{eplx:.3f} mas (SNR={plx_snr:.1f})")
    parts.append(f"pm={pm:.2f}±{pm_err:.2f} mas/yr (SNR={pm_snr:.1f})")

    triggers: list[str] = []
    if is_fg_parallax:
        triggers.append(
            f"PARALLAX: plx≥{cfg.gaia_parallax_min_mas} mas and SNR≥{cfg.gaia_parallax_snr_min}"
        )
    if is_fg_pm:
        triggers.append(f"PM: pm≥{cfg.gaia_pm_min_masyr} mas/yr and SNR≥{cfg.gaia_pm_snr_min}")

    if len(triggers) == 0:
        triggers_str = "NOT foreground by kinematics thresholds"
    else:
        triggers_str = "foreground because " + " OR ".join(triggers)

    return " | ".join(parts) + " | " + triggers_str


def _virgo_velocity_range_kms(cfg: Config) -> tuple[float, float]:
    d0 = float(cfg.virgo_distance_mpc)
    dd = float(cfg.virgo_distance_tolerance_mpc)
    h0 = float(cfg.hubble_km_s_mpc)
    vmin = h0 * max(0.0, d0 - dd)
    vmax = h0 * (d0 + dd)
    return vmin, vmax


def _is_virgo_distance_from_z_or_v(cfg: Config, z: float | None, v_kms: float | None) -> bool:
    vmin, vmax = _virgo_velocity_range_kms(cfg)
    if v_kms is not None and np.isfinite(v_kms):
        return bool(vmin <= float(v_kms) <= vmax)
    if z is not None and np.isfinite(z):
        # Non-relativistic cz is fine at Virgo.
        cz = 299792.458 * float(z)
        return bool(vmin <= cz <= vmax)
    return False


def query_ned_redshifts(center: SkyCoord, radius: u.Quantity):
    """Query NED around the field, returning a table with positions and redshift/velocity when available."""
    if Ned is None:
        return None
    try:
        # Extragalactic objects only; returns columns including RA, DEC and often Redshift.
        tab = Ned.query_region(center, radius=radius, equinox="J2000.0")
        return tab
    except Exception:
        return None


def query_ps1_galaxy_like(center: SkyCoord, radius: u.Quantity, cfg: Config):
    if Catalogs is None:
        print("WARNING: astroquery.mast not available; skipping Pan-STARRS query.")
        return query_ps1_vizier(center, radius, cfg)

    # Pan-STARRS coverage is Dec >= -30 deg. We'll still try; if empty, fallback elsewhere.
    try:
        tab = Catalogs.query_region(
            center,
            radius=radius,
            catalog="Panstarrs",
            data_release="dr2",
        )
        return tab
    except Exception as e:
        print(f"WARNING: Pan-STARRS (MAST) query failed; falling back to VizieR. ({e})")
        return query_ps1_vizier(center, radius, cfg)


def query_ps1_vizier(center: SkyCoord, radius: u.Quantity, cfg: Config):
    if Vizier is None:
        print("WARNING: astroquery.vizier not available; cannot query Pan-STARRS via VizieR.")
        return None

    # VizieR Pan-STARRS1 catalog (provides PSF mags and Kron-like mags as *Kmag*)
    cols = [
        "RAJ2000", "DEJ2000", "Nr", "Qual",
        "gmag", "rmag", "imag", "gKmag", "rKmag", "iKmag",
        "e_gmag", "e_rmag", "e_imag", "e_gKmag", "e_rKmag", "e_iKmag",
    ]
    v = Vizier(
        columns=cols,
        row_limit=300000,
        column_filters={
            "rmag": f"<{cfg.ps1_rmag_max}",
            "Nr": f">={cfg.ps1_min_Nr}",
        },
    )
    try:
        res = v.query_region(center, radius=radius, catalog="II/349/ps1")
        if len(res) == 0:
            return None
        tab = res[0]
    except Exception as e:
        print(f"WARNING: Pan-STARRS (VizieR) query failed; skipping. ({e})")
        return None
    return tab


def query_sdss_photoobj(center: SkyCoord, radius: u.Quantity, cfg: Config):
    """Query SDSS photometric objects in the field (supplemental; footprint-limited)."""
    if SDSS is None:
        return None
    # SDSS query_region enforces radius <= 3 arcmin
    r = radius.to(u.arcmin)
    if r > 3.0 * u.arcmin:
        r = 3.0 * u.arcmin
    try:
        # Try a richer set of fields (morphology proxy + robust radii).
        # If SDSS rejects any field name (schema/DR differences), fall back.
        photo_fields_full = [
            "ra", "dec", "objid", "type",
            "psfMag_r", "modelMag_r",
            "petroR50_r", "petroR90_r",
            "expRad_r", "deVRad_r",
            "petroRad_r",
        ]
        try:
            tab = SDSS.query_region(
                center,
                radius=r,
                spectro=False,
                photoobj_fields=photo_fields_full,
                data_release=int(getattr(cfg, "sdss_data_release", 17)),
            )
            return tab
        except Exception:
            photo_fields_fallback = [
                "ra", "dec", "objid", "type",
                "petroRad_r", "petroRad_g", "petroRad_i",
                "psfMag_r", "modelMag_r",
            ]
            tab = SDSS.query_region(
                center,
                radius=r,
                spectro=False,
                photoobj_fields=photo_fields_fallback,
                data_release=int(getattr(cfg, "sdss_data_release", 17)),
            )
            return tab
    except Exception:
        return None


def query_sdss_spectro(center: SkyCoord, radius: u.Quantity, cfg: Config):
    """Get SDSS spectroscopic objects (with z) within the field, if SDSS is available."""
    if SDSS is None:
        return None
    r = radius.to(u.arcmin)
    if r > 3.0 * u.arcmin:
        r = 3.0 * u.arcmin
    try:
        tab = SDSS.query_region(
            center,
            radius=r,
            spectro=True,
            data_release=int(getattr(cfg, "sdss_data_release", 17)),
        )
        return tab
    except Exception:
        return None


def _cz_from_z(z: float) -> float:
    return 299792.458 * float(z)


def _get_col_float(row, name: str):
    try:
        v = row[name]
        if v is None or getattr(v, "mask", False):
            return None

        # Some astroquery/astropy tables can yield 0-d numpy scalars, masked
        # scalars, or (rarely) 1-element arrays; normalize those to a Python float.
        if isinstance(v, np.ma.MaskedArray):
            if v is np.ma.masked or np.any(getattr(v, "mask", False)):
                return None
            v = v.data
        if isinstance(v, np.ndarray) and v.ndim > 0:
            if v.size != 1:
                return None
            v = v.reshape(-1)[0]

        vf = float(v)
        return vf if np.isfinite(vf) else None
    except Exception:
        return None


def _as_str(x) -> str:
    if x is None:
        return ""
    if isinstance(x, (bytes, bytearray)):
        try:
            return x.decode(errors="ignore")
        except Exception:
            return ""
    return str(x)


def query_legacy_dr9_tractor_and_photoz(center: SkyCoord, radius: u.Quantity, cfg: Config):
    """Query Legacy Surveys DR9 via NOIRLab Data Lab TAP.

        Returns an astropy Table joined via `ls_id` with columns:
            ls_id, ra, dec, type, release, brickid, objid, z_phot_l95, z_phot_u95, z_phot_mean
        and (optionally) shape_r, shape_e1, shape_e2 if present.

    If TAP/pyvo/tables are unavailable, returns None.
    """

    if not bool(getattr(cfg, "enable_legacy", True)):
        return None
    if pyvo is None:
        print("WARNING: pyvo not available; skipping Legacy DR9 TAP queries.")
        return None

    rad_deg = float(radius.to_value(u.deg))
    ra0 = float(center.ra.deg)
    dec0 = float(center.dec.deg)

    def _tap_service():
        # Reuse service across calls to reduce overhead.
        if not hasattr(_tap_service, "svc"):
            _tap_service.svc = pyvo.dal.TAPService(str(getattr(cfg, "legacy_tap_url", "https://datalab.noirlab.edu/tap")))
        return _tap_service.svc

    def _run_sync(query: str, maxrec: int | None = None):
        import time
        max_retries = 3
        for attempt in range(max_retries):
            try:
                if maxrec is None:
                    return _tap_service().run_sync(query).to_table()
                return _tap_service().run_sync(query, maxrec=int(maxrec)).to_table()
            except Exception as e:
                # Catch general exceptions including requests.exceptions.HTTPError
                error_msg = str(e)
                if "500" in error_msg or "503" in error_msg or "504" in error_msg or "timeout" in error_msg.lower() or "connection" in error_msg.lower() or "Error 500" in error_msg:
                    if attempt < max_retries - 1:
                        print(f"WARNING: Legacy TAP query failed ({error_msg}). Retrying {attempt+1}/{max_retries}...")
                        time.sleep(2 * (attempt + 1))
                        continue
                raise

    # 1) tractor-like table (morphology/type + potential sizes)
    tractor_table = str(getattr(cfg, "legacy_tractor_table", "ls_dr9.tractor"))

    # NOTE: This particular service supports `q3c_radial_query` but (in practice)
    # does not accept ADQL POINT/CIRCLE geometry functions.
    # Also, its SQL parser does not handle boolean literals well, so we compare
    # the boolean-returning function to the Postgres-friendly string 't'.
    tractor_q_with_shape = f"""
    SELECT ls_id, ra, dec, type, release, brickid, objid,
           shape_r, shape_e1, shape_e2
    FROM {tractor_table}
    WHERE q3c_radial_query(ra, dec, {ra0}, {dec0}, {rad_deg}) = 't'
    """.strip()
    tractor_q_min = f"""
    SELECT ls_id, ra, dec, type, release, brickid, objid
    FROM {tractor_table}
    WHERE q3c_radial_query(ra, dec, {ra0}, {dec0}, {rad_deg}) = 't'
    """.strip()

    ttab = None
    tractor_err = None
    try:
        try:
            ttab = _run_sync(tractor_q_with_shape, maxrec=20000)
        except Exception:
            ttab = _run_sync(tractor_q_min, maxrec=20000)
    except Exception as e:
        tractor_err = e
        ttab = None

    if ttab is None:
        print(f"WARNING: Legacy tractor TAP query failed: {tractor_err}")
        return None

    if ttab is None or len(ttab) == 0:
        return None

    # 2) photo-z table (z columns). NOTE: `ls_dr9.photo_z` does not expose ra/dec,
    # so we join via `ls_id`.
    photoz_table = str(getattr(cfg, "legacy_photoz_table", "ls_dr9.photo_z"))

    ls_ids = []
    for r in ttab:
        try:
            ls_ids.append(int(r["ls_id"]))
        except Exception:
            continue
    if len(ls_ids) == 0:
        return None

    pkey = {}
    last_err = None
    # Chunk the IN-list to avoid overlong queries.
    for i0 in range(0, len(ls_ids), 500):
        chunk = ls_ids[i0 : i0 + 500]
        in_list = ",".join(str(int(x)) for x in chunk)
        q = f"""
        SELECT ls_id,
               z_phot_l95 AS z_phot_l95,
               z_phot_u95 AS z_phot_u95,
               z_phot_mean AS z_phot_mean
        FROM {photoz_table}
        WHERE ls_id IN ({in_list})
        """.strip()
        try:
            ptab = _run_sync(q, maxrec=20000)
            last_err = None
        except Exception as e:
            last_err = e
            continue
        for pr in ptab:
            try:
                pkey[int(pr["ls_id"])] = pr
            except Exception:
                continue

    if len(pkey) == 0:
        print(f"WARNING: Legacy photo-z TAP query failed: {last_err}")
        return None

    rows = []
    for r in ttab:
        try:
            lid = int(r["ls_id"])
        except Exception:
            continue
        pr = pkey.get(lid)
        if pr is None:
            continue
        rows.append((r, pr))

    if len(rows) == 0:
        return None

    from astropy.table import Table

    out = Table()
    out["ls_id"] = [int(r["ls_id"]) for r, _ in rows]
    out["ra"] = [float(r["ra"]) for r, _ in rows]
    out["dec"] = [float(r["dec"]) for r, _ in rows]
    out["type"] = [_as_str(r["type"]).strip() for r, _ in rows]
    out["release"] = [int(r["release"]) for r, _ in rows]
    out["brickid"] = [int(r["brickid"]) for r, _ in rows]
    out["objid"] = [int(r["objid"]) for r, _ in rows]

    for c in ["shape_r", "shape_e1", "shape_e2"]:
        if c in ttab.colnames:
            out[c] = [r[c] for r, _ in rows]

    out["z_phot_l95"] = [float(p["z_phot_l95"]) for _, p in rows]
    out["z_phot_u95"] = [float(p["z_phot_u95"]) for _, p in rows]
    out["z_phot_mean"] = [float(p["z_phot_mean"]) for _, p in rows]

    return out


def sdss_is_pointlike(row, cfg: Config) -> bool:
    psf = _get_col_float(row, "psfMag_r")
    mod = _get_col_float(row, "modelMag_r")
    if psf is None or mod is None:
        return False
    if (not np.isfinite(psf)) or (not np.isfinite(mod)):
        return False
    # SDSS can use sentinel/unphysical magnitudes; don't let those drive morphology.
    if (psf < -5.0) or (psf > 40.0) or (mod < -5.0) or (mod > 40.0):
        return False
    return bool((float(psf) - float(mod)) < float(getattr(cfg, "sdss_pointlike_dmag_max", 0.145)))


def sdss_mask_radius_arcsec(row, cfg: Config) -> float | None:
    # If point-like (star/QSO), mask like a star using r-band magnitude as a proxy.
    if sdss_is_pointlike(row, cfg):
        mag = _get_col_float(row, "psfMag_r")
        if mag is None:
            mag = _get_col_float(row, "modelMag_r")
        if mag is None:
            mag = 18.0
        return star_radius_arcsec_from_g(cfg, float(mag))

    # Extended: prefer petroR90, then petroR50, then model radii.
    r90 = _get_col_float(row, "petroR90_r")
    if r90 is not None:
        return float(getattr(cfg, "sdss_petro_r90_scale", 1.2)) * float(r90) + float(cfg.gaia_margin_arcsec)

    r50 = _get_col_float(row, "petroR50_r")
    if r50 is not None:
        return float(getattr(cfg, "sdss_petro_r50_scale", 1.8)) * float(r50) + float(cfg.gaia_margin_arcsec)

    dev = _get_col_float(row, "deVRad_r")
    exp = _get_col_float(row, "expRad_r")
    model_rad = None
    if dev is not None and exp is not None:
        model_rad = max(float(dev), float(exp))
    elif dev is not None:
        model_rad = float(dev)
    elif exp is not None:
        model_rad = float(exp)
    if model_rad is not None:
        return float(getattr(cfg, "sdss_model_radius_scale", 3.0)) * float(model_rad) + float(cfg.gaia_margin_arcsec)

    return None


def _extract_cz_kms_from_sdss_specrow(row):
    z = _get_col_float(row, "z") or _get_col_float(row, "Z")
    if z is None:
        return None
    return _cz_from_z(z)


def _extract_cz_kms_from_ned_row(row):
    def _ned_is_spec(r) -> bool:
        # NED includes a "Redshift Flag" column for many entries.
        # We only trust spectroscopic redshifts for definitive masking.
        try:
            if "Redshift Flag" not in r.colnames:
                return False
            flag = r["Redshift Flag"]
            flag = flag.decode() if isinstance(flag, (bytes, bytearray)) else str(flag)
            flag = flag.upper()
            return "SPEC" in flag
        except Exception:
            return False

    if not _ned_is_spec(row):
        return None

    v = None
    for vname in ["Velocity", "cz", "Vel", "V"]:
        v = _get_col_float(row, vname)
        if v is not None:
            return float(v)
    z = None
    for zname in ["Redshift", "z", "Z"]:
        z = _get_col_float(row, zname)
        if z is not None:
            return _cz_from_z(z)
    return None


def is_definitely_background(
    sc: SkyCoord,
    cfg: Config,
    sdss_spec_sky: SkyCoord | None,
    sdss_spec_tab,
    ned_sky: SkyCoord | None,
    ned_tab,
) -> bool | None:
    """
    Returns:
      True  -> definitely background (mask)
      False -> definitely Virgo/nearby or uncertain (do not mask)
      None  -> no distance info available (do not mask, for safety)
    """
    # 1) SDSS spectroscopy (preferred)
    if sdss_spec_sky is not None and sdss_spec_tab is not None and len(sdss_spec_tab) > 0:
        try:
            idx, sep2d, _ = sc.match_to_catalog_sky(sdss_spec_sky)
            if sep2d < (float(cfg.sdss_spec_match_arcsec) * u.arcsec):
                cz = _extract_cz_kms_from_sdss_specrow(sdss_spec_tab[idx])
                if cz is None:
                    return None
                if cz <= float(cfg.virgo_keep_cz_max_kms):
                    return False
                if cz >= float(cfg.background_mask_cz_min_kms):
                    return True
                return False
        except Exception:
            pass
    # 2) NED (fallback)
    if ned_sky is not None and ned_tab is not None and len(ned_tab) > 0:
        try:
            idx, sep2d, _ = sc.match_to_catalog_sky(ned_sky)
            if sep2d < (float(cfg.ned_match_arcsec) * u.arcsec):
                cz = _extract_cz_kms_from_ned_row(ned_tab[idx])
                if cz is None:
                    return None
                if cz <= float(cfg.virgo_keep_cz_max_kms):
                    return False
                if cz >= float(cfg.background_mask_cz_min_kms):
                    return True
                return False
        except Exception:
            pass
    return None


def get_best_cz_info(
    sc: SkyCoord,
    cfg: Config,
    sdss_spec_sky: SkyCoord | None,
    sdss_spec_tab,
    ned_sky: SkyCoord | None,
    ned_tab,
):
    """
    Returns (cz_kms, z, source, sep_arcsec) or (None, None, None, None)
    """
    # SDSS spectroscopy first
    if sdss_spec_sky is not None and sdss_spec_tab is not None and len(sdss_spec_tab) > 0:
        try:
            idx, sep2d, _ = sc.match_to_catalog_sky(sdss_spec_sky)
            if sep2d < (float(cfg.sdss_spec_match_arcsec) * u.arcsec):
                cz = _extract_cz_kms_from_sdss_specrow(sdss_spec_tab[idx])
                if cz is not None:
                    z = float(cz) / 299792.458
                    sep_arcsec = sep2d.to_value(u.arcsec)
                    if isinstance(sep_arcsec, np.ndarray) and sep_arcsec.ndim > 0:
                        if sep_arcsec.size == 1:
                            sep_arcsec = sep_arcsec.reshape(-1)[0]
                    return float(cz), z, "SDSS(spec)", float(sep_arcsec)
        except Exception:
            pass

    # NED fallback
    if ned_sky is not None and ned_tab is not None and len(ned_tab) > 0:
        try:
            idx, sep2d, _ = sc.match_to_catalog_sky(ned_sky)
            if sep2d < (float(cfg.ned_match_arcsec) * u.arcsec):
                cz = _extract_cz_kms_from_ned_row(ned_tab[idx])
                if cz is not None:
                    z = float(cz) / 299792.458
                    sep_arcsec = sep2d.to_value(u.arcsec)
                    if isinstance(sep_arcsec, np.ndarray) and sep_arcsec.ndim > 0:
                        if sep_arcsec.size == 1:
                            sep_arcsec = sep_arcsec.reshape(-1)[0]
                    return float(cz), z, "NED(spec)", float(sep_arcsec)
        except Exception:
            pass

    return None, None, None, None


def query_skymapper(center: SkyCoord, radius: u.Quantity):
    if Vizier is None:
        print("WARNING: astroquery.vizier not available; skipping SkyMapper query.")
        return None

    v = Vizier(columns=["**"], row_limit=200000)
    try:
        # SkyMapper DR4 Vizier table
        # II/379/smssdr4 (very large table; cone searches are OK but can be slower)
        res = v.query_region(center, radius=radius, catalog="II/379/smssdr4")
        if len(res) == 0:
            return None
        return res[0]
    except Exception as e:
        print(f"WARNING: SkyMapper (VizieR) query failed; skipping. ({e})")
        return None


def pick_first_existing_col(table, candidates):
    for c in candidates:
        if c in table.colnames:
            return c
    return None


def _to_float_or_nan(x):
    try:
        if np.ma.is_masked(x):
            return float("nan")
    except Exception:
        pass
    try:
        return float(x)
    except Exception:
        return float("nan")


def rasterize_circle(mask: np.ndarray, xi: float, yi: float, r_pix: float, fov_mask: np.ndarray | None = None) -> None:
    ny, nx = mask.shape
    if not np.isfinite(xi) or not np.isfinite(yi) or not np.isfinite(r_pix) or r_pix <= 0:
        return
    x0 = max(0, int(np.floor(xi - r_pix)))
    x1 = min(nx - 1, int(np.ceil(xi + r_pix)))
    y0 = max(0, int(np.floor(yi - r_pix)))
    y1 = min(ny - 1, int(np.ceil(yi + r_pix)))
    if x1 < x0 or y1 < y0:
        return
    yy, xx = np.ogrid[y0 : y1 + 1, x0 : x1 + 1]
    rr2 = (xx - xi) ** 2 + (yy - yi) ** 2
    inside = rr2 <= (r_pix**2)
    if fov_mask is not None:
        fsub = fov_mask[y0 : y1 + 1, x0 : x1 + 1]
        if fsub.shape == inside.shape:
            mask[y0 : y1 + 1, x0 : x1 + 1][inside & fsub] = 1
        return
    mask[y0 : y1 + 1, x0 : x1 + 1][inside] = 1


def circle_intersects_fov(xi: float, yi: float, r_pix: float, nx: int, ny: int) -> bool:
    """True if the circle overlaps the image rectangle, even partially."""
    if not (np.isfinite(xi) and np.isfinite(yi) and np.isfinite(r_pix)) or r_pix <= 0:
        return False
    # Use pixel-edge coordinates [-0.5, nx-0.5] etc. for correct partial overlap logic
    return not (
        (xi + r_pix) < -0.5 or (xi - r_pix) > (nx - 0.5) or
        (yi + r_pix) < -0.5 or (yi - r_pix) > (ny - 0.5)
    )


def ellipse_intersects_fov(xi: float, yi: float, a_pix: float, b_pix: float, nx: int, ny: int) -> bool:
    """Conservative overlap test using the ellipse bounding circle."""
    r = float(max(a_pix, b_pix))
    return circle_intersects_fov(xi, yi, r, nx, ny)


def rasterize_ellipse(mask: np.ndarray, xi: float, yi: float, a_pix: float, b_pix: float, angle_deg: float, fov_mask: np.ndarray | None = None) -> None:
    ny, nx = mask.shape
    if (
        not np.isfinite(xi)
        or not np.isfinite(yi)
        or not np.isfinite(a_pix)
        or not np.isfinite(b_pix)
        or a_pix <= 0
        or b_pix <= 0
    ):
        return
    r = float(max(a_pix, b_pix))
    x0 = max(0, int(np.floor(xi - r)))
    x1 = min(nx - 1, int(np.ceil(xi + r)))
    y0 = max(0, int(np.floor(yi - r)))
    y1 = min(ny - 1, int(np.ceil(yi + r)))
    if x1 < x0 or y1 < y0:
        return
    yy, xx = np.ogrid[y0 : y1 + 1, x0 : x1 + 1]
    th = np.deg2rad(angle_deg)
    xp = (xx - xi) * np.cos(th) + (yy - yi) * np.sin(th)
    yp = -(xx - xi) * np.sin(th) + (yy - yi) * np.cos(th)
    inside = (xp / a_pix) ** 2 + (yp / b_pix) ** 2 <= 1.0
    if fov_mask is not None:
        fsub = fov_mask[y0 : y1 + 1, x0 : x1 + 1]
        if fsub.shape == inside.shape:
            mask[y0 : y1 + 1, x0 : x1 + 1][inside & fsub] = 1
        return
    mask[y0 : y1 + 1, x0 : x1 + 1][inside] = 1


def rasterize_polygon(mask: np.ndarray, xverts: np.ndarray, yverts: np.ndarray, fov_mask: np.ndarray | None = None) -> None:
    ny, nx = mask.shape
    x = np.asarray(xverts, dtype=float)
    y = np.asarray(yverts, dtype=float)
    if x.size < 3 or y.size < 3:
        return
    if x.size != y.size:
        return
    if not (np.all(np.isfinite(x)) and np.all(np.isfinite(y))):
        return

    verts = np.column_stack([x, y])
    if not np.allclose(verts[0], verts[-1]):
        verts = np.vstack([verts, verts[0]])

    x0 = max(0, int(np.floor(np.min(verts[:, 0]))))
    x1 = min(nx - 1, int(np.ceil(np.max(verts[:, 0]))))
    y0 = max(0, int(np.floor(np.min(verts[:, 1]))))
    y1 = min(ny - 1, int(np.ceil(np.max(verts[:, 1]))))
    if x1 < x0 or y1 < y0:
        return

    xx, yy = np.meshgrid(np.arange(x0, x1 + 1), np.arange(y0, y1 + 1))
    pts = np.column_stack([xx.ravel(), yy.ravel()])
    inside = np.asarray(Path(verts).contains_points(pts), dtype=bool)
    inside = inside.reshape(yy.shape)
    if fov_mask is not None:
        fsub = fov_mask[y0 : y1 + 1, x0 : x1 + 1]
        if fsub.shape == inside.shape:
            mask[y0 : y1 + 1, x0 : x1 + 1][inside & fsub] = 1
        return
    mask[y0 : y1 + 1, x0 : x1 + 1][inside] = 1


def circle_overlaps_mask(fov_mask: np.ndarray | None, xi: float, yi: float, r_pix: float) -> bool:
    if fov_mask is None:
        return True
    if not np.isfinite(xi) or not np.isfinite(yi) or not np.isfinite(r_pix) or r_pix <= 0:
        return False
    ny, nx = fov_mask.shape
    x0 = max(0, int(np.floor(xi - r_pix)))
    x1 = min(nx - 1, int(np.ceil(xi + r_pix)))
    y0 = max(0, int(np.floor(yi - r_pix)))
    y1 = min(ny - 1, int(np.ceil(yi + r_pix)))
    if x1 < x0 or y1 < y0:
        return False
    sub = fov_mask[y0 : y1 + 1, x0 : x1 + 1]
    if not np.any(sub):
        return False
    yy, xx = np.ogrid[y0 : y1 + 1, x0 : x1 + 1]
    rr2 = (xx - xi) ** 2 + (yy - yi) ** 2
    inside = rr2 <= (r_pix**2)
    return bool(np.any(sub & inside))


def circle_within_mask(fov_mask: np.ndarray | None, xi: float, yi: float, r_pix: float) -> bool:
    if fov_mask is None:
        return True
    if not np.isfinite(xi) or not np.isfinite(yi) or not np.isfinite(r_pix) or r_pix <= 0:
        return False
    ny, nx = fov_mask.shape
    x0 = max(0, int(np.floor(xi - r_pix)))
    x1 = min(nx - 1, int(np.ceil(xi + r_pix)))
    y0 = max(0, int(np.floor(yi - r_pix)))
    y1 = min(ny - 1, int(np.ceil(yi + r_pix)))
    if x1 < x0 or y1 < y0:
        return False
    yy, xx = np.ogrid[y0 : y1 + 1, x0 : x1 + 1]
    rr2 = (xx - xi) ** 2 + (yy - yi) ** 2
    inside = rr2 <= (r_pix**2)
    fsub = fov_mask[y0 : y1 + 1, x0 : x1 + 1]
    if fsub.shape != inside.shape:
        return False
    return bool(np.all(fsub[inside]))


def ellipse_overlaps_mask(fov_mask: np.ndarray | None, xi: float, yi: float, a_pix: float, b_pix: float, angle_deg: float) -> bool:
    if fov_mask is None:
        return True
    if (
        not np.isfinite(xi)
        or not np.isfinite(yi)
        or not np.isfinite(a_pix)
        or not np.isfinite(b_pix)
        or a_pix <= 0
        or b_pix <= 0
    ):
        return False
    ny, nx = fov_mask.shape
    r = float(max(a_pix, b_pix))
    x0 = max(0, int(np.floor(xi - r)))
    x1 = min(nx - 1, int(np.ceil(xi + r)))
    y0 = max(0, int(np.floor(yi - r)))
    y1 = min(ny - 1, int(np.ceil(yi + r)))
    if x1 < x0 or y1 < y0:
        return False
    sub = fov_mask[y0 : y1 + 1, x0 : x1 + 1]
    if not np.any(sub):
        return False
    yy, xx = np.ogrid[y0 : y1 + 1, x0 : x1 + 1]
    th = np.deg2rad(angle_deg)
    xp = (xx - xi) * np.cos(th) + (yy - yi) * np.sin(th)
    yp = -(xx - xi) * np.sin(th) + (yy - yi) * np.cos(th)
    inside = (xp / a_pix) ** 2 + (yp / b_pix) ** 2 <= 1.0
    return bool(np.any(sub & inside))


def polygon_overlaps_mask(fov_mask: np.ndarray | None, xverts: np.ndarray, yverts: np.ndarray) -> bool:
    if fov_mask is None:
        return True
    x = np.asarray(xverts, dtype=float)
    y = np.asarray(yverts, dtype=float)
    if x.size < 3 or y.size < 3:
        return False
    if x.size != y.size:
        return False
    if not (np.all(np.isfinite(x)) and np.all(np.isfinite(y))):
        return False

    ny, nx = fov_mask.shape
    verts = np.column_stack([x, y])
    if not np.allclose(verts[0], verts[-1]):
        verts = np.vstack([verts, verts[0]])

    x0 = max(0, int(np.floor(np.min(verts[:, 0]))))
    x1 = min(nx - 1, int(np.ceil(np.max(verts[:, 0]))))
    y0 = max(0, int(np.floor(np.min(verts[:, 1]))))
    y1 = min(ny - 1, int(np.ceil(np.max(verts[:, 1]))))
    if x1 < x0 or y1 < y0:
        return False
    sub = fov_mask[y0 : y1 + 1, x0 : x1 + 1]
    if not np.any(sub):
        return False

    xx, yy = np.meshgrid(np.arange(x0, x1 + 1), np.arange(y0, y1 + 1))
    pts = np.column_stack([xx.ravel(), yy.ravel()])
    inside = np.asarray(Path(verts).contains_points(pts), dtype=bool)
    inside = inside.reshape(yy.shape)
    return bool(np.any(sub & inside))


def _ds9_region_header() -> list[str]:
    return [
        "# Region file format: DS9 version 4.1",
        'global color=cyan dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1',
        "fk5",
    ]


def _ds9_fmt_angle_deg(x: float) -> float:
    # DS9 ellipses treat PA modulo 180 as equivalent.
    a = float(x) % 180.0
    if a < 0:
        a += 180.0
    return a


def write_ds9_legacy_regions(
    base: str,
    entries: list[dict],
    suffix: str,
    *,
    color: str = "cyan",
) -> str:
    """Write a DS9 FK5 region file for the given entries.

    Each entry dict must contain: ra_deg, dec_deg, a_arcsec, b_arcsec, pa_deg.
    """
    out = f"{base}_legacy_{suffix}.reg"
    lines = _ds9_region_header()
    # Override color per file for quick visual A/B.
    lines[1] = lines[1].replace("color=cyan", f"color={color}")
    for e in entries:
        ra = float(e["ra_deg"])
        dec = float(e["dec_deg"])
        a = float(e["a_arcsec"])
        b = float(e["b_arcsec"])
        pa = _ds9_fmt_angle_deg(float(e["pa_deg"]))
        tag = str(e.get("tag", ""))
        # DS9 expects a and b as semi-axes in arcsec.
        s = f"ellipse({ra:.8f},{dec:.8f},{a:.4f}\",{b:.4f}\",{pa:.4f})"
        if tag:
            s += f" # text={{{tag}}}"
        lines.append(s)
    with open(out, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    return out


def sample_ellipse_via_wcs(
    w: WCS,
    center: SkyCoord,
    a_arcsec: float,
    b_arcsec: float,
    angle_deg: float,
    npts: int = 96,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample an ellipse in the local tangent plane and map it through the WCS.

    The tangent-plane basis is (east, north) about `center`. `angle_deg` is
    interpreted as the major-axis angle measured CCW from +east toward +north.
    """
    n = int(max(16, npts))
    t = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    x_maj = float(a_arcsec) * np.cos(t)
    y_min = float(b_arcsec) * np.sin(t)
    phi = np.deg2rad(float(angle_deg))

    east = x_maj * np.cos(phi) - y_min * np.sin(phi)
    north = x_maj * np.sin(phi) + y_min * np.cos(phi)

    sky = center.spherical_offsets_by(east * u.arcsec, north * u.arcsec)
    xp, yp = w.world_to_pixel(sky)
    return np.asarray(xp, dtype=float), np.asarray(yp, dtype=float)


def make_target_footprint_from_rmag(
    rfits_path: str,
    cfg: Config,
    shape: tuple[int, int] | None = None,
) -> tuple[np.ndarray | None, float | None]:
    """Build a target-galaxy footprint mask directly from the R_MAG HDU.

    Returns (footprint_mask, mu_lim). If R_MAG is missing/unreadable, returns (None, None).
    """
    try:
        with fits.open(rfits_path) as hdul:
            rmag_hdu = None
            for hdu in hdul:
                if getattr(hdu, "name", "").strip().upper() == "R_MAG":
                    rmag_hdu = hdu
                    break
            if rmag_hdu is None:
                print(f"WARNING: R_MAG extension not found in {rfits_path}; skipping target footprint.")
                return None, None
            data = rmag_hdu.data
            if data is None:
                print(f"WARNING: R_MAG extension has no data; skipping target footprint.")
                return None, None
    except Exception as e:
        print(f"WARNING: failed to read R_MAG from {rfits_path}: {e}")
        return None, None

    if data.ndim == 3:
        img = np.nanmedian(data, axis=0)
    elif data.ndim == 2:
        img = data
    else:
        print(f"WARNING: R_MAG has ndim={data.ndim}; skipping target footprint.")
        return None, None

    if shape is not None and img.shape != shape:
        if img.shape == (shape[1], shape[0]):
            img = img.T
        else:
            print(f"WARNING: R_MAG shape {img.shape} != expected {shape}; skipping target footprint.")
            return None, None

    mu_lim = float(getattr(cfg, "target_mu_lim", 26.0))
    footprint = np.isfinite(img) & (img < mu_lim)
    return footprint, mu_lim


def build_masks_for_one(rfits_path: str, cfg: Config):
    base = safe_base_id(rfits_path)
    png_path = f"{base}_combined_VRI.png"
    out_mask_fits = f"{base}_mask.fits"
    out_overlay_png = f"{base}_combined_VRI_mask.png"

    use_png_bg = bool(cfg.use_png_background and os.path.exists(png_path) and Image is not None)

    data2d, w, hdr, nx, ny = load_r_image_and_wcs(rfits_path)
    pixscale = pixel_scale_arcsec(w)
    center, rad = fov_center_and_radius(w, nx, ny)
    # small padding so we don't miss edge objects
    rad = rad + (10.0 * u.arcsec)

    try:
        print(f"[FOV] {base}: nx={nx} ny={ny} pixscale={pixscale:.4f}\"/pix rad={rad.to(u.arcmin).value:.3f} arcmin")
    except Exception:
        pass

    mask = np.zeros((ny, nx), dtype=np.uint8)
    exclude_center = cfg.exclude_center_arcsec * u.arcsec

    fov_mask = None
    fov_mask_plot = None
    if bool(getattr(cfg, "fov_use_mask", True)):
        fov_mask = build_muse_fov_mask(
            rfits_path,
            getattr(cfg, "fov_extname", "R_FLUX"),
            getattr(cfg, "fov_close_size_pix", 5),
            getattr(cfg, "fov_min_abs", 0.0),
        )
        if fov_mask is not None and fov_mask.shape != (ny, nx):
            if fov_mask.shape == (nx, ny):
                print("WARNING: FoV mask shape appears transposed; applying transpose to match data.")
                fov_mask = fov_mask.T
            else:
                print(f"WARNING: FoV mask shape {fov_mask.shape} != data shape {(ny, nx)}; FoV gating disabled.")
                fov_mask = None
        if fov_mask is not None and not np.any(fov_mask):
            print("WARNING: FoV mask is empty (all False). Gating will suppress all masks.")
        if fov_mask is not None:
            fov_mask_plot = np.flipud(fov_mask) if bool(getattr(cfg, "fov_flip_y", False)) else fov_mask
            try:
                frac = float(np.sum(fov_mask)) / float(fov_mask.size)
                print(f"[FOV] FoV mask coverage: {100.0 * frac:.2f}%")
            except Exception:
                pass
    if fov_mask is None:
        fov_mask = np.isfinite(data2d)
        print("WARNING: FoV mask unavailable; using finite data2d as FoV to prevent masking outside image.")
        fov_mask_plot = fov_mask

    star_edge_buffer_pix = None
    if fov_mask is not None:
        star_edge_buffer_pix = float(getattr(cfg, "fov_edge_star_buffer_arcsec", 5.0)) / float(pixscale)

    # ---------- Gaia stars ----------
    gaia = query_gaia_sources(center, rad, cfg)
    star_patches = []
    n_star_masked = 0
    star_exclude = np.zeros((ny, nx), dtype=np.uint8)
    gaia_sky = None
    gaia_sky_for_ps1_reject = None
    if gaia is not None and len(gaia) > 0:
        ra = np.array(gaia["ra"])
        dec = np.array(gaia["dec"])
        gmag = np.array(gaia["phot_g_mean_mag"])
        gaia_sky = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")

        # For PS1 galaxy rejection we only want high-quality Gaia point sources.
        # In loose mode, Gaia contains more dubious detections (and even some
        # galaxy cores), which would otherwise suppress real background galaxies.
        try:
            ruwe = np.array(gaia["ruwe"], dtype=float)
            ipd = np.array(gaia["ipd_frac_multi_peak"], dtype=float)
            exsig = np.array(gaia["astrometric_excess_noise_sig"], dtype=float)

            ok_ruwe = np.isnan(ruwe) | (ruwe < float(cfg.gaia_ruwe_max))
            ok_ipd = np.isnan(ipd) | (ipd <= float(cfg.gaia_ipd_frac_multi_peak_max))
            ok_exsig = np.isnan(exsig) | (exsig <= float(cfg.gaia_astrometric_excess_noise_sig_max))
            good_for_reject = ok_ruwe & ok_ipd & ok_exsig
            gaia_sky_for_ps1_reject = gaia_sky[good_for_reject]
        except Exception:
            gaia_sky_for_ps1_reject = gaia_sky

        x, y = w.world_to_pixel(gaia_sky)  # float pixel coords

        mode = str(getattr(cfg, "gaia_star_mode", "strict")).lower().strip()

        if cfg.log_each_star:
            print(
                f"[Gaia] mode={mode} (foreground thresholds: "
                f"plx≥{cfg.gaia_parallax_min_mas} mas & SNR≥{cfg.gaia_parallax_snr_min} "
                f"OR pm≥{cfg.gaia_pm_min_masyr} mas/yr & SNR≥{cfg.gaia_pm_snr_min})"
            )

        for row, xi, yi, gi, sc in zip(gaia, x, y, gmag, gaia_sky):
            if not np.isfinite(xi) or not np.isfinite(yi):
                continue

            # Optional: only mask sources that look like Milky Way stars.
            # This avoids masking Virgo-distance / extragalactic Gaia detections.
            if mode == "foreground":
                try:
                    if not _gaia_row_is_foreground_by_kinematics(row, cfg):
                        continue
                except Exception:
                    # If kinematics are unavailable/unparsable, don't mask it in foreground mode.
                    continue

            if sc.separation(center) < exclude_center:
                continue
            r_arcsec = star_radius_arcsec_from_g(cfg, float(gi))
            r_pix = r_arcsec / pixscale

            # Skip if it does not touch the FITS image at all
            if not circle_intersects_fov(float(xi), float(yi), float(r_pix), nx, ny):
                continue

            # FoV gating: ignore stars fully outside FoV; if edge-buffer only partially overlaps,
            # mask only the buffer-sized circle.
            if fov_mask is not None:
                if not circle_overlaps_mask(fov_mask, float(xi), float(yi), float(r_pix)):
                    continue
                r_mask_pix = float(r_pix)
                if star_edge_buffer_pix is not None and float(star_edge_buffer_pix) > 0:
                    if not circle_overlaps_mask(fov_mask, float(xi), float(yi), float(star_edge_buffer_pix)):
                        continue
                    if not circle_within_mask(fov_mask, float(xi), float(yi), float(star_edge_buffer_pix)):
                        r_mask_pix = float(min(r_pix, float(star_edge_buffer_pix)))
            else:
                r_mask_pix = float(r_pix)

            # rasterize into mask
            rasterize_circle(mask, xi, yi, r_mask_pix, fov_mask=fov_mask)

            # Rasterize an expanded star-exclusion mask for footprint estimation
            r_arcsec_fp = float(r_arcsec) + (2.0 * float(cfg.fwhm_arcsec))
            r_pix_fp = r_arcsec_fp / float(pixscale)
            if fov_mask is not None and r_mask_pix < float(r_pix):
                r_pix_fp = min(float(r_pix_fp), float(r_mask_pix))
            rasterize_circle(star_exclude, xi, yi, r_pix_fp, fov_mask=fov_mask)

            sid = row["source_id"] if "source_id" in row.colnames else "?"
            if cfg.log_each_star:
                ra_hms, dec_dms = format_radec_hmsdms(sc, precision=2)
                gaia_name = f"Gaia DR3 {sid}"
                gaia_iau = iau_coord_name("GAIA", sc, ra_precision=2, dec_precision=1)
                reason = ""
                if mode == "foreground":
                    try:
                        reason = " | " + gaia_foreground_reason(row, cfg)
                    except Exception:
                        reason = ""
                print(
                    f"[STAR] {gaia_name} ({gaia_iau}) ra={sc.ra.deg:.6f} dec={sc.dec.deg:.6f} "
                    f"RA={ra_hms} DEC={dec_dms} "
                    f"GaiaG={float(gi):.2f} r_arcsec={r_arcsec:.2f}{reason}"
                )
            n_star_masked += 1

            # The provided PNG background images are typically rendered with a
            # different vertical origin than raw FITS array indices.
            y_plot = (ny - 1 - yi) if use_png_bg else yi
            star_patches.append(Circle((xi, y_plot), r_mask_pix, fill=False))

    target_fp = None
    target_thresh = None
    if bool(getattr(cfg, "reject_bg_inside_target_footprint", False)) or bool(getattr(cfg, "draw_target_iso_contour", False)):
        target_fp, target_thresh = make_target_footprint_from_rmag(
            rfits_path,
            cfg,
            shape=(ny, nx),
        )

    # ---------- Background galaxies (Pan-STARRS preferred) ----------
    gal_patches = []
    gal_tab = None
    gal_catalog = None
    n_gal_masked = 0
    legacy_success = False
    legacy_attempted = False
    legacy_ds9_entries_pa_eofn: list[dict] = []
    legacy_ds9_entries_phi_from_east: list[dict] = []

    # ---------- Legacy Surveys DR9 photo-z-gated background objects (DEFAULT first-pass) ----------
    if bool(getattr(cfg, "enable_legacy", True)):
        legacy_attempted = True
        legacy = query_legacy_dr9_tractor_and_photoz(center, rad, cfg)
        if legacy is not None and len(legacy) > 0:
            n_gal_before_legacy = int(n_gal_masked)
            try:
                print(f"[LEGACY] N(joined)={len(legacy)} cols={list(legacy.colnames)}")
            except Exception:
                pass

            legacy_sky = SkyCoord(
                ra=np.array(legacy["ra"], dtype=float) * u.deg,
                dec=np.array(legacy["dec"], dtype=float) * u.deg,
                frame="icrs",
            )
            xL, yL = w.world_to_pixel(legacy_sky)

            for row, xi, yi, sc in zip(legacy, xL, yL, legacy_sky):
                if not (np.isfinite(xi) and np.isfinite(yi)):
                    continue
                if exclude_center.value > 0 and sc.separation(center) < exclude_center:
                    continue

                if bool(getattr(cfg, "reject_bg_inside_target_footprint", False)) and (target_fp is not None):
                    xi_i = int(round(float(xi)))
                    yi_i = int(round(float(yi)))
                    if (0 <= xi_i < nx) and (0 <= yi_i < ny) and bool(target_fp[yi_i, xi_i]):
                        continue

                typ = _as_str(row["type"]).strip().upper() if "type" in row.colnames else ""
                z_l95 = _get_col_float(row, "z_phot_l95")
                if z_l95 is None or (not np.isfinite(z_l95)):
                    continue

                # Rule: non-PSF AND lower-95% bound > threshold
                if typ == "PSF":
                    continue
                if float(z_l95) <= float(getattr(cfg, "legacy_z_l95_min", 0.01)):
                    continue

                z_mean = _get_col_float(row, "z_phot_mean")
                z_u95 = _get_col_float(row, "z_phot_u95")
                if z_mean is None or z_u95 is None:
                    continue

                z_width = float(z_u95) - float(z_l95)
                if (not np.isfinite(z_width)) or (z_width <= 0):
                    continue
                if float(getattr(cfg, "legacy_z_width_max", 0.0)) > 0 and z_width > float(cfg.legacy_z_width_max):
                    continue

                sigma_z = z_width / 3.92
                if (not np.isfinite(sigma_z)) or (sigma_z <= 0):
                    continue

                z_cut = float(getattr(cfg, "legacy_z_l95_min", 0.2))
                if bool(getattr(cfg, "legacy_z_snr_use_threshold", True)):
                    z_snr = (float(z_mean) - z_cut) / sigma_z
                else:
                    z_snr = float(z_mean) / sigma_z

                if (not np.isfinite(z_snr)) or (z_snr < float(getattr(cfg, "legacy_z_snr_min", 3.0))):
                    continue

                # Optional: reject Legacy detections sitting on Gaia point sources
                # (often star halos / bad fits misclassified as extended)
                if gaia_sky_for_ps1_reject is not None and float(getattr(cfg, "legacy_reject_if_near_gaia_arcsec", 0.0)) > 0:
                    try:
                        _, sepg, _ = sc.match_to_catalog_sky(gaia_sky_for_ps1_reject)
                        if sepg < (float(cfg.legacy_reject_if_near_gaia_arcsec) * u.arcsec):
                            continue
                    except Exception:
                        pass

                sr = _get_col_float(row, "shape_r")
                e1 = _get_col_float(row, "shape_e1")
                e2 = _get_col_float(row, "shape_e2")

                # --- Choose intrinsic size ---
                if sr is not None and np.isfinite(sr) and float(sr) > 0:
                    a_arcsec = float(cfg.legacy_shape_r_scale) * float(sr)
                    a_arcsec = max(float(cfg.legacy_r_min_arcsec), a_arcsec)

                    if bool(getattr(cfg, "legacy_use_ellipses", False)) and (e1 is not None) and (e2 is not None) and np.isfinite(e1) and np.isfinite(e2):
                        e = float(np.hypot(float(e1), float(e2)))
                        e = float(np.clip(e, 0.0, 0.85))
                        q = (1.0 - e) / (1.0 + e)  # axis ratio b/a
                        q = float(np.clip(q, 0.2, 1.0))
                        b_arcsec = max(float(cfg.legacy_r_min_arcsec), a_arcsec * q)
                        angle_deg = float(np.degrees(0.5 * np.arctan2(float(e2), float(e1))))
                    else:
                        b_arcsec = a_arcsec
                        angle_deg = 0.0
                else:
                    # No Tractor size -> fall back to seeing (but still respect the 1" floor)
                    a_arcsec = float(getattr(cfg, "legacy_seeing_scale", 1.0)) * float(cfg.fwhm_arcsec)
                    a_arcsec = max(float(cfg.legacy_r_min_arcsec), a_arcsec)
                    b_arcsec = a_arcsec
                    angle_deg = 0.0

                # --- Apply caps (optional separate Legacy cap; otherwise use global) ---
                rmax = float(getattr(cfg, "legacy_r_max_arcsec", cfg.gal_r_max_arcsec))
                a_arcsec = float(np.clip(a_arcsec, float(cfg.legacy_r_min_arcsec), rmax))
                b_arcsec = float(np.clip(b_arcsec, float(cfg.legacy_r_min_arcsec), rmax))

                # --- Add padding margin ---
                a_arcsec += float(cfg.gaia_margin_arcsec)
                b_arcsec += float(cfg.gaia_margin_arcsec)

                # --- DS9 region export (FK5; for PA convention validation) ---
                # Tractor phi is commonly interpreted as angle CCW from +east toward +north.
                # DS9 expects PA degrees East of North.
                if bool(getattr(cfg, "legacy_write_ds9_regions", False)):
                    try:
                        tag = f"{typ} z_l95={float(z_l95):.3f}"
                    except Exception:
                        tag = f"{typ}"
                    try:
                        legacy_ds9_entries_phi_from_east.append(
                            {
                                "ra_deg": float(sc.ra.deg),
                                "dec_deg": float(sc.dec.deg),
                                "a_arcsec": float(a_arcsec),
                                "b_arcsec": float(b_arcsec),
                                "pa_deg": float(angle_deg),
                                "tag": tag,
                            }
                        )
                        pa_eofn = 90.0 - float(angle_deg)
                        legacy_ds9_entries_pa_eofn.append(
                            {
                                "ra_deg": float(sc.ra.deg),
                                "dec_deg": float(sc.dec.deg),
                                "a_arcsec": float(a_arcsec),
                                "b_arcsec": float(b_arcsec),
                                "pa_deg": float(pa_eofn),
                                "tag": tag,
                            }
                        )
                    except Exception:
                        pass

                a_pix = a_arcsec / float(pixscale)
                b_pix = b_arcsec / float(pixscale)

                # --- Force Legacy to circles (robust to PA convention issues) ---
                if bool(getattr(cfg, "legacy_force_circles", False)):
                    r_arcsec = float(max(a_arcsec, b_arcsec))
                    r_pix = r_arcsec / float(pixscale)
                    if not circle_intersects_fov(float(xi), float(yi), float(r_pix), nx, ny):
                        continue
                    if fov_mask is not None and not circle_overlaps_mask(fov_mask, float(xi), float(yi), float(r_pix)):
                        continue
                    rasterize_circle(mask, xi, yi, r_pix, fov_mask=fov_mask)
                    y_plot = (ny - 1 - yi) if use_png_bg else yi
                    gal_patches.append(Circle((xi, y_plot), r_pix, fill=False))
                else:
                    # Original behavior: ellipses when available, otherwise circles.
                    if bool(getattr(cfg, "legacy_use_ellipses", False)) and (abs(float(a_arcsec) - float(b_arcsec)) > 1e-6):
                        if not ellipse_intersects_fov(float(xi), float(yi), float(a_pix), float(b_pix), nx, ny):
                            continue
                        if fov_mask is not None and not ellipse_overlaps_mask(fov_mask, float(xi), float(yi), float(a_pix), float(b_pix), float(angle_deg)):
                            continue
                        use_wcs_sampling = bool(getattr(cfg, "legacy_wcs_sample_ellipses", True))
                        if use_wcs_sampling:
                            try:
                                angle_for_sampling = float(angle_deg)
                                if bool(getattr(cfg, "legacy_pa_east_of_north", False)):
                                    angle_for_sampling = 90.0 - angle_for_sampling
                                angle_for_sampling = angle_for_sampling + float(getattr(cfg, "legacy_pa_offset_deg", 0.0))
                                xv, yv = sample_ellipse_via_wcs(
                                    w,
                                    sc,
                                    float(a_arcsec),
                                    float(b_arcsec),
                                    float(angle_for_sampling),
                                    npts=int(getattr(cfg, "legacy_ellipse_npts", 96)),
                                )
                                if fov_mask is not None and not polygon_overlaps_mask(fov_mask, xv, yv):
                                    continue
                                rasterize_polygon(mask, xv, yv, fov_mask=fov_mask)
                                if use_png_bg:
                                    yv_plot = (ny - 1) - yv
                                else:
                                    yv_plot = yv
                                gal_patches.append(Polygon(np.column_stack([xv, yv_plot]), closed=True, fill=False))
                            except Exception:
                                rasterize_ellipse(mask, xi, yi, a_pix, b_pix, angle_deg, fov_mask=fov_mask)
                                y_plot = (ny - 1 - yi) if use_png_bg else yi
                                angle_plot = (-angle_deg) if use_png_bg else angle_deg
                                gal_patches.append(Ellipse((xi, y_plot), 2 * a_pix, 2 * b_pix, angle=angle_plot, fill=False))
                        else:
                            rasterize_ellipse(mask, xi, yi, a_pix, b_pix, angle_deg, fov_mask=fov_mask)
                            y_plot = (ny - 1 - yi) if use_png_bg else yi
                            angle_plot = (-angle_deg) if use_png_bg else angle_deg
                            gal_patches.append(Ellipse((xi, y_plot), 2 * a_pix, 2 * b_pix, angle=angle_plot, fill=False))
                    else:
                        r_pix = float(max(a_pix, b_pix))
                        if not circle_intersects_fov(float(xi), float(yi), float(r_pix), nx, ny):
                            continue
                        if fov_mask is not None and not circle_overlaps_mask(fov_mask, float(xi), float(yi), float(r_pix)):
                            continue
                        rasterize_circle(mask, xi, yi, r_pix, fov_mask=fov_mask)
                        y_plot = (ny - 1 - yi) if use_png_bg else yi
                        gal_patches.append(Circle((xi, y_plot), r_pix, fill=False))
                n_gal_masked += 1

                if cfg.log_each_galaxy:
                    ra_hms, dec_dms = format_radec_hmsdms(sc, precision=2)
                    legacy_iau = iau_coord_name("LS", sc, ra_precision=2, dec_precision=1)
                    if bool(getattr(cfg, "legacy_force_circles", False)):
                        r_arcsec = float(max(a_arcsec, b_arcsec))
                        print(
                            f"[GAL][LEGACY] ({legacy_iau}) ra={sc.ra.deg:.6f} dec={sc.dec.deg:.6f} "
                            f"RA={ra_hms} DEC={dec_dms} type={typ} z_l95={float(z_l95):.3f} "
                            f"r={r_arcsec:.2f}\" (circle from max axis) z_mean={float(z_mean):.3f} "
                            f"z_u95={float(z_u95):.3f} sigma_z~{float(sigma_z):.3f} SNR={float(z_snr):.1f}"
                        )
                    else:
                        print(
                            f"[GAL][LEGACY] ({legacy_iau}) ra={sc.ra.deg:.6f} dec={sc.dec.deg:.6f} "
                            f"RA={ra_hms} DEC={dec_dms} type={typ} z_l95={float(z_l95):.3f} "
                            f"a={a_arcsec:.2f}\" b={b_arcsec:.2f}\" pa={float(angle_deg):.1f} "
                            f"z_mean={float(z_mean):.3f} z_u95={float(z_u95):.3f} "
                            f"sigma_z~{float(sigma_z):.3f} SNR={float(z_snr):.1f}"
                        )

            legacy_success = bool(n_gal_masked > n_gal_before_legacy)
            if legacy_success:
                print("[LEGACY] Legacy masking succeeded; skipping PS1/SDSS/NED fallback catalogs.")
                if bool(getattr(cfg, "legacy_write_ds9_regions", False)) and len(legacy_ds9_entries_pa_eofn) > 0:
                    try:
                        p1 = write_ds9_legacy_regions(base, legacy_ds9_entries_pa_eofn, "PA_EofN", color="cyan")
                        p2 = write_ds9_legacy_regions(base, legacy_ds9_entries_phi_from_east, "phi_from_E", color="magenta")
                        print(f"[DS9] Wrote {p1}")
                        print(f"[DS9] Wrote {p2}")
                        print("[DS9] Open in DS9 over the FITS; PA_EofN should match if the convention is correct.")
                    except Exception as e:
                        print(f"WARNING: failed to write DS9 region files: {e}")
        else:
            print("[LEGACY] No DR9 photo-z-gated objects found (or query unavailable).")

    # Optional: NED crossmatch for redshift-based background determination
    ned_tab = None
    ned_sky = None
    if (not legacy_attempted) and (
        bool(getattr(cfg, "enable_virgo_distance_veto", False))
        or bool(getattr(cfg, "require_nonvirgo_confirmation_for_galaxy_mask", False))
    ):
        ned_tab = query_ned_redshifts(center, rad)
        try:
            if ned_tab is not None and len(ned_tab) > 0 and "RA" in ned_tab.colnames and "DEC" in ned_tab.colnames:
                ned_sky = SkyCoord(ra=np.array(ned_tab["RA"]) * u.deg, dec=np.array(ned_tab["DEC"]) * u.deg, frame="icrs")
        except Exception:
            ned_sky = None

    try:
        if not legacy_attempted:
            print(f"[NED] N={0 if ned_tab is None else len(ned_tab)} cols={None if ned_tab is None else ned_tab.colnames}")
    except Exception:
        pass

    # SDSS spectroscopy table (has redshifts; preferred for evidence-based masking)
    sdss_spec_tab = query_sdss_spectro(center, rad, cfg) if ((not legacy_attempted) and bool(getattr(cfg, "enable_sdss", True))) else None
    sdss_spec_sky = None
    if sdss_spec_tab is not None and len(sdss_spec_tab) > 0:
        ra_c = pick_first_existing_col(sdss_spec_tab, ["ra", "RA"])
        dec_c = pick_first_existing_col(sdss_spec_tab, ["dec", "DEC"])
        if ra_c and dec_c:
            sdss_spec_sky = SkyCoord(
                ra=np.array(sdss_spec_tab[ra_c]) * u.deg,
                dec=np.array(sdss_spec_tab[dec_c]) * u.deg,
                frame="icrs",
            )

    try:
        if not legacy_attempted:
            print(
                f"[SDSS-spec] N={0 if sdss_spec_tab is None else len(sdss_spec_tab)} cols={None if sdss_spec_tab is None else sdss_spec_tab.colnames}"
            )
    except Exception:
        pass

    if (not legacy_attempted) and center.dec.deg >= -30:
        gal_tab = query_ps1_galaxy_like(center, rad, cfg)
        if gal_tab is not None and len(gal_tab) > 0:
            gal_catalog = "PS1"

    # If PS1 failed/empty and SkyMapper is available, try SkyMapper (useful mainly in the south)
    if (not legacy_attempted) and (gal_tab is None or len(gal_tab) == 0) and center.dec.deg <= +28:
        gal_tab = query_skymapper(center, rad)
        if gal_tab is not None and len(gal_tab) > 0:
            gal_catalog = "SkyMapper"

    if gal_catalog is None:
        gal_catalog = "CAT"

    if gal_tab is not None and len(gal_tab) > 0:
        # Dynamic column picking (PS1 and SkyMapper differ)
        ra_col = pick_first_existing_col(gal_tab, ["raMean", "RAJ2000", "ra", "_VRIAJ2000"])
        dec_col = pick_first_existing_col(gal_tab, ["decMean", "DEJ2000", "dec", "_DEJ2000"])

        if ra_col and dec_col:
            sky = SkyCoord(ra=np.array(gal_tab[ra_col]) * u.deg,
                           dec=np.array(gal_tab[dec_col]) * u.deg,
                           frame="icrs")
            x, y = w.world_to_pixel(sky)

            # Try PS1-style extendedness. Depending on source, columns differ:
            # - MAST PS1: rMeanPSFMag / rMeanKronMag
            # - VizieR PS1 (II/349/ps1): rmag (PSF-like) and rKmag (Kron-like)
            psf_col = pick_first_existing_col(
                gal_tab,
                [
                    "rMeanPSFMag", "iMeanPSFMag", "gMeanPSFMag",
                    "rmag", "imag", "gmag",
                ],
            )
            kron_col = pick_first_existing_col(
                gal_tab,
                [
                    "rMeanKronMag", "iMeanKronMag", "gMeanKronMag",
                    "rKmag", "iKmag", "gKmag",
                ],
            )

            g_psf_col = pick_first_existing_col(gal_tab, ["gMeanPSFMag", "gmag"])
            r_psf_col = pick_first_existing_col(gal_tab, ["rMeanPSFMag", "rmag"])
            i_psf_col = pick_first_existing_col(gal_tab, ["iMeanPSFMag", "imag"])

            g_err_col = pick_first_existing_col(gal_tab, ["gMeanPSFMagErr", "e_gmag"])
            r_err_col = pick_first_existing_col(gal_tab, ["rMeanPSFMagErr", "e_rmag"])
            i_err_col = pick_first_existing_col(gal_tab, ["iMeanPSFMagErr", "e_imag"])
            psf_err_col = pick_first_existing_col(
                gal_tab,
                [
                    "rMeanPSFMagErr", "iMeanPSFMagErr", "gMeanPSFMagErr",
                    "e_rmag", "e_imag", "e_gmag",
                ],
            )
            kron_err_col = pick_first_existing_col(
                gal_tab,
                [
                    "rMeanKronMagErr", "iMeanKronMagErr", "gMeanKronMagErr",
                    "e_rKmag", "e_iKmag", "e_gKmag",
                ],
            )
            rmag_col = pick_first_existing_col(
                gal_tab,
                [
                    "rMeanKronMag", "rMeanPSFMag",
                    "rKmag", "rmag",
                    "rPSF", "r",
                ],
            )

            # Try size columns that might exist:
            kronrad_col = pick_first_existing_col(gal_tab, ["rKronRad", "rMeanKronRad", "iKronRad", "KronRad"])
            hlr_col = pick_first_existing_col(gal_tab, ["rHalfLightRad", "iHalfLightRad", "halfLightRadius"])

            # If SkyMapper: try semi-major/minor & PA if present
            a_col = pick_first_existing_col(gal_tab, ["a", "A", "semimajor", "aWorld"])
            b_col = pick_first_existing_col(gal_tab, ["b", "B", "semiminor", "bWorld"])
            pa_col = pick_first_existing_col(gal_tab, ["pa", "PA", "theta", "posang"])

            i_kron_col = pick_first_existing_col(gal_tab, ["iMeanKronMag", "iKmag"])

            gal_logged = 0
            gal_log_limit = int(getattr(cfg, "log_max_galaxies", 0) or 0)

            for i in range(len(gal_tab)):
                xi, yi = float(x[i]), float(y[i])
                if not np.isfinite(xi) or not np.isfinite(yi):
                    continue

                sc = sky[i]
                cat_iau = iau_coord_name(gal_catalog, sc, ra_precision=2, dec_precision=1)
                if exclude_center.value > 0 and sc.separation(center) < exclude_center:
                    continue

                if bool(getattr(cfg, "reject_bg_inside_target_footprint", False)) and (target_fp is not None):
                    xi_i = int(round(float(xi)))
                    yi_i = int(round(float(yi)))
                    if (0 <= xi_i < nx) and (0 <= yi_i < ny) and bool(target_fp[yi_i, xi_i]):
                        continue

                # If it matches a Gaia source, treat it as stellar (already handled by Gaia masks)
                # and do NOT count it as a background galaxy.
                if gaia_sky_for_ps1_reject is not None and cfg.ps1_reject_if_near_gaia_arcsec > 0:
                    try:
                        idx, sep2d, _ = sc.match_to_catalog_sky(gaia_sky_for_ps1_reject)
                        if sep2d < (cfg.ps1_reject_if_near_gaia_arcsec * u.arcsec):
                            continue
                    except Exception:
                        pass

                # Evidence-based masking: only mask if we have definitive redshift proof it's background
                if cfg.require_nonvirgo_confirmation_for_galaxy_mask:
                    bg = is_definitely_background(
                        sc=sc,
                        cfg=cfg,
                        sdss_spec_sky=sdss_spec_sky,
                        sdss_spec_tab=sdss_spec_tab,
                        ned_sky=ned_sky,
                        ned_tab=ned_tab,
                    )
                    if bg is True:
                        pass  # mask
                    elif bg is None:
                        # --- fallback for unknown distance: mask only VERY extended and bright ---
                        # (tune these to taste)
                        if psf_col and kron_col:
                            psf = _to_float_or_nan(gal_tab[psf_col][i])
                            kron = _to_float_or_nan(gal_tab[kron_col][i])
                            ext = (psf - kron) if (np.isfinite(psf) and np.isfinite(kron)) else -np.inf
                        else:
                            ext = -np.inf

                        rmag = _to_float_or_nan(gal_tab[r_psf_col][i]) if r_psf_col else np.inf

                        if not bool(getattr(cfg, "ps1_allow_photometric_fallback", True)):
                            continue

                        if not (
                            float(ext) > float(getattr(cfg, "ps1_fallback_ext_min", 0.8))
                            and float(rmag) < float(getattr(cfg, "ps1_fallback_rmag_max", 20.0))
                        ):
                            continue
                    else:
                        # cz says Virgo/nearby (or otherwise not definitely background) -> do not mask
                        continue
                # Fallback: old Virgo-distance veto (kept for compatibility if new flag is disabled)
                elif ned_sky is not None and cfg.virgo_match_arcsec > 0:
                    try:
                        nidx, nsep, _ = sc.match_to_catalog_sky(ned_sky)
                        if nsep < (float(cfg.virgo_match_arcsec) * u.arcsec):
                            z = None
                            v_kms = None
                            for zname in ["Redshift", "z", "Z"]:
                                if zname in ned_tab.colnames:
                                    try:
                                        zv = ned_tab[zname][nidx]
                                        if zv is None or getattr(zv, "mask", False):
                                            z = None
                                        else:
                                            z = float(zv)
                                    except Exception:
                                        z = None
                                    break
                            for vname in ["Velocity", "cz", "Vel", "V"]:
                                if vname in ned_tab.colnames:
                                    try:
                                        vv = ned_tab[vname][nidx]
                                        if vv is None or getattr(vv, "mask", False):
                                            v_kms = None
                                        else:
                                            v_kms = float(vv)
                                    except Exception:
                                        v_kms = None
                                    break
                            if _is_virgo_distance_from_z_or_v(cfg, z=z, v_kms=v_kms):
                                continue
                    except Exception:
                        pass

                # Catalog-based quality filters (VizieR PS1 provides Nd/Nr and per-band errors).
                if "Nr" in gal_tab.colnames:
                    try:
                        if int(gal_tab["Nr"][i]) < int(cfg.ps1_min_Nr):
                            continue
                    except Exception:
                        pass

                # If available, prefer the PS1 qualityFlag (Qual) bitmask to ensure
                # we're selecting real extended objects (i.e., background galaxies).
                if "Qual" in gal_tab.colnames:
                    try:
                        q = int(gal_tab["Qual"][i])
                        if int(getattr(cfg, "ps1_qual_extended_required_bits", 0)):
                            if (q & int(cfg.ps1_qual_extended_required_bits)) != int(cfg.ps1_qual_extended_required_bits):
                                continue
                        if cfg.ps1_require_qual_good and (q & (4 | 16)) != (4 | 16):
                            continue
                        if cfg.ps1_require_qual_primary_best and (q & 32) == 0:
                            continue
                        if cfg.ps1_reject_qual_suspect and (q & (64 | 128)) != 0:
                            continue
                    except Exception:
                        pass

                # Very strict color cuts (if the needed bands are present)
                if cfg.ps1_enable_color_cuts and g_psf_col and r_psf_col and i_psf_col:
                    g = _to_float_or_nan(gal_tab[g_psf_col][i])
                    r = _to_float_or_nan(gal_tab[r_psf_col][i])
                    ii = _to_float_or_nan(gal_tab[i_psf_col][i])
                    if not (np.isfinite(g) and np.isfinite(r) and np.isfinite(ii)):
                        continue
                    # Require reasonably good errors if available
                    if g_err_col and r_err_col and i_err_col:
                        eg = _to_float_or_nan(gal_tab[g_err_col][i])
                        er = _to_float_or_nan(gal_tab[r_err_col][i])
                        ei = _to_float_or_nan(gal_tab[i_err_col][i])
                        if np.isfinite(eg) and eg > float(cfg.ps1_e_mag_max):
                            continue
                        if np.isfinite(er) and er > float(cfg.ps1_e_mag_max):
                            continue
                        if np.isfinite(ei) and ei > float(cfg.ps1_e_mag_max):
                            continue
                    if (g - r) < float(cfg.ps1_g_r_min):
                        continue
                    if (r - ii) < float(cfg.ps1_r_i_min):
                        continue

                if psf_err_col and kron_err_col:
                    try:
                        e1 = _to_float_or_nan(gal_tab[psf_err_col][i])
                        e2 = _to_float_or_nan(gal_tab[kron_err_col][i])
                        if np.isfinite(e1) and e1 > float(cfg.ps1_e_mag_max):
                            continue
                        if np.isfinite(e2) and e2 > float(cfg.ps1_e_mag_max):
                            continue
                    except Exception:
                        pass

                # Optional magnitude cut
                if rmag_col and np.isfinite(gal_tab[rmag_col][i]):
                    if float(gal_tab[rmag_col][i]) > cfg.ps1_rmag_max:
                        continue

                # Decide if "galaxy-like"
                is_gal_like = False
                if psf_col and kron_col:
                    psf = _to_float_or_nan(gal_tab[psf_col][i])
                    kron = _to_float_or_nan(gal_tab[kron_col][i])
                    if np.isfinite(psf) and np.isfinite(kron):
                        ext_val_r = float(psf - kron)
                        ext_r = (ext_val_r > float(cfg.ps1_ext_thresh)) and (ext_val_r < float(cfg.ps1_ext_max))
                    else:
                        ext_r = False

                    # If i-band equivalents exist, optionally require extendedness in i as well.
                    if cfg.ps1_require_ri_extended:
                        if i_psf_col and i_kron_col:
                            ipsf = _to_float_or_nan(gal_tab[i_psf_col][i])
                            ikron = _to_float_or_nan(gal_tab[i_kron_col][i])
                            if np.isfinite(ipsf) and np.isfinite(ikron):
                                ext_val_i = float(ipsf - ikron)
                                ext_i = (ext_val_i > float(cfg.ps1_ext_thresh)) and (ext_val_i < float(cfg.ps1_ext_max))
                            else:
                                ext_i = False
                            is_gal_like = bool(ext_r and ext_i)
                        else:
                            is_gal_like = bool(ext_r)
                    else:
                        is_gal_like = bool(ext_r)
                else:
                    # If we can't classify, be conservative but not crazy: treat as gal-like only if size exists
                    is_gal_like = (kronrad_col is not None) or (hlr_col is not None) or (a_col is not None)

                # If Qual is present, it already encodes "extended" and real-vs-false-positive.
                # Keep the PSF-Kron extendedness requirement as an additional guard when available,
                # but do not allow objects that fail Qual cuts above.

                if not is_gal_like:
                    continue

                # Determine mask size.
                if kronrad_col and np.isfinite(gal_tab[kronrad_col][i]):
                    r_arcsec = float(gal_tab[kronrad_col][i])
                elif hlr_col and np.isfinite(gal_tab[hlr_col][i]):
                    r_arcsec = 2.0 * float(gal_tab[hlr_col][i])
                elif a_col and b_col and np.isfinite(gal_tab[a_col][i]) and np.isfinite(gal_tab[b_col][i]):
                    # We'll use ellipse below; just set a representative radius for logging.
                    r_arcsec = float(max(gal_tab[a_col][i], gal_tab[b_col][i]))
                elif psf_col and kron_col:
                    # VizieR PS1 lacks size; use a small fallback ONLY when extendedness is available.
                    r_arcsec = float(cfg.gal_fallback_arcsec)
                else:
                    continue

                r_arcsec = float(np.clip(r_arcsec, cfg.gal_r_min_arcsec, cfg.gal_r_max_arcsec))
                r_pix = r_arcsec / pixscale

                objid_col = pick_first_existing_col(gal_tab, ["objID", "ObjID", "source_id", "ID", "Name"])
                objid = "?"
                if objid_col is not None:
                    try:
                        objid = str(gal_tab[objid_col][i])
                    except Exception:
                        objid = "?"

                # Ellipse if axis ratio exists; otherwise circle
                if a_col and b_col and np.isfinite(gal_tab[a_col][i]) and np.isfinite(gal_tab[b_col][i]):
                    a_arcsec = float(gal_tab[a_col][i])
                    b_arcsec = float(gal_tab[b_col][i])
                    a_arcsec = float(np.clip(a_arcsec, cfg.gal_r_min_arcsec, cfg.gal_r_max_arcsec))
                    b_arcsec = float(np.clip(b_arcsec, cfg.gal_r_min_arcsec, cfg.gal_r_max_arcsec))
                    a_pix = a_arcsec / pixscale
                    b_pix = b_arcsec / pixscale
                    angle = 0.0
                    if pa_col and np.isfinite(gal_tab[pa_col][i]):
                        angle = float(gal_tab[pa_col][i])

                    # Skip if ellipse does not touch the FITS image at all
                    if not ellipse_intersects_fov(float(xi), float(yi), float(a_pix), float(b_pix), nx, ny):
                        continue
                    if fov_mask is not None and not ellipse_overlaps_mask(fov_mask, float(xi), float(yi), float(a_pix), float(b_pix), float(angle)):
                        continue

                    # Plot coords (only affects overlay, not the FITS mask)
                    y_plot = (ny - 1 - yi) if use_png_bg else yi
                    angle_plot = (-angle) if use_png_bg else angle

                    # rasterize ellipse (local cutout)
                    rasterize_ellipse(mask, xi, yi, a_pix, b_pix, angle, fov_mask=fov_mask)

                    # Get cz info for logging
                    cz, zval, zsrc, zsep = get_best_cz_info(sc, cfg, sdss_spec_sky, sdss_spec_tab, ned_sky, ned_tab)
                    d_mpc = (cz / float(cfg.hubble_km_s_mpc)) if (cz is not None) else None

                    if cz is None and cfg.log_each_galaxy and (gal_log_limit <= 0 or gal_logged < gal_log_limit):
                        # optional: show nearest-match separations for debugging
                        if sdss_spec_sky is not None and len(sdss_spec_sky) > 0:
                            try:
                                idx, sep2d, _ = sc.match_to_catalog_sky(sdss_spec_sky)
                                print(f"[DBG] nearest SDSS(spec) sep={sep2d.to_value(u.arcsec):.2f}\"")
                            except Exception:
                                pass
                        if ned_sky is not None and len(ned_sky) > 0:
                            try:
                                idx, sep2d, _ = sc.match_to_catalog_sky(ned_sky)
                                print(f"[DBG] nearest NED sep={sep2d.to_value(u.arcsec):.2f}\"")
                            except Exception:
                                pass

                    if cfg.log_each_galaxy and (gal_log_limit <= 0 or gal_logged < gal_log_limit):
                        ra_hms, dec_dms = format_radec_hmsdms(sc, precision=2)
                        qv = gal_tab["Qual"][i] if "Qual" in gal_tab.colnames else "?"
                        if cz is not None:
                            print(
                                f"[GAL] {objid} ({cat_iau}) ra={sc.ra.deg:.6f} dec={sc.dec.deg:.6f} "
                                f"RA={ra_hms} DEC={dec_dms} "
                                f"cz={cz:.0f}km/s z={zval:.5f} D~{d_mpc:.1f}Mpc ({zsrc}) "
                                f"a={a_arcsec:.2f}\" b={b_arcsec:.2f}\" pa={angle:.1f}"
                            )
                        else:
                            print(
                                f"[GAL] {objid} ({cat_iau}) ra={sc.ra.deg:.6f} dec={sc.dec.deg:.6f} "
                                f"RA={ra_hms} DEC={dec_dms} "
                                f"a={a_arcsec:.2f}\" b={b_arcsec:.2f}\" pa={angle:.1f} Qual={qv}"
                            )
                        gal_logged += 1
                    n_gal_masked += 1

                    gal_patches.append(Ellipse((xi, y_plot), 2 * a_pix, 2 * b_pix, angle=angle_plot, fill=False))
                else:
                    # Skip if circle does not touch the FITS image at all
                    if not circle_intersects_fov(float(xi), float(yi), float(r_pix), nx, ny):
                        continue

                    if fov_mask is not None and not circle_overlaps_mask(fov_mask, float(xi), float(yi), float(r_pix)):
                        continue

                    if fov_mask is not None and not circle_overlaps_mask(fov_mask, float(xi), float(yi), float(r_pix)):
                        continue

                    # rasterize circle
                    rasterize_circle(mask, xi, yi, r_pix, fov_mask=fov_mask)

                    # Get cz info for logging
                    cz, zval, zsrc, zsep = get_best_cz_info(sc, cfg, sdss_spec_sky, sdss_spec_tab, ned_sky, ned_tab)
                    d_mpc = (cz / float(cfg.hubble_km_s_mpc)) if (cz is not None) else None

                    if cz is None and cfg.log_each_galaxy and (gal_log_limit <= 0 or gal_logged < gal_log_limit):
                        # optional: show nearest-match separations for debugging
                        if sdss_spec_sky is not None and len(sdss_spec_sky) > 0:
                            try:
                                idx, sep2d, _ = sc.match_to_catalog_sky(sdss_spec_sky)
                                print(f"[DBG] nearest SDSS(spec) sep={sep2d.to_value(u.arcsec):.2f}\"")
                            except Exception:
                                pass
                        if ned_sky is not None and len(ned_sky) > 0:
                            try:
                                idx, sep2d, _ = sc.match_to_catalog_sky(ned_sky)
                                print(f"[DBG] nearest NED sep={sep2d.to_value(u.arcsec):.2f}\"")
                            except Exception:
                                pass

                    if cfg.log_each_galaxy and (gal_log_limit <= 0 or gal_logged < gal_log_limit):
                        ra_hms, dec_dms = format_radec_hmsdms(sc, precision=2)
                        qv = gal_tab["Qual"][i] if "Qual" in gal_tab.colnames else "?"
                        if cz is not None:
                            print(
                                f"[GAL] {objid} ({cat_iau}) ra={sc.ra.deg:.6f} dec={sc.dec.deg:.6f} "
                                f"RA={ra_hms} DEC={dec_dms} "
                                f"cz={cz:.0f}km/s z={zval:.5f} D~{d_mpc:.1f}Mpc ({zsrc}) "
                                f"r_arcsec={r_arcsec:.2f}"
                            )
                        else:
                            # if PSF/Kron mags available, print extendedness for debugging
                            ext = "?"
                            if psf_col and kron_col:
                                psf = _to_float_or_nan(gal_tab[psf_col][i])
                                kron = _to_float_or_nan(gal_tab[kron_col][i])
                                if np.isfinite(psf) and np.isfinite(kron):
                                    ext = f"{(psf-kron):.3f}"
                            nr = gal_tab["Nr"][i] if "Nr" in gal_tab.colnames else "?"
                            print(
                                f"[GAL] {objid} ({cat_iau}) ra={sc.ra.deg:.6f} dec={sc.dec.deg:.6f} "
                                f"RA={ra_hms} DEC={dec_dms} "
                                f"r_arcsec={r_arcsec:.2f} ext={ext} Nr={nr} Qual={qv}"
                            )
                        gal_logged += 1
                    n_gal_masked += 1
                    y_plot = (ny - 1 - yi) if use_png_bg else yi
                    gal_patches.append(Circle((xi, y_plot), r_pix, fill=False))

            if cfg.log_each_galaxy and gal_log_limit > 0 and n_gal_masked > gal_logged:
                print(f"[GAL] ... suppressed {n_gal_masked - gal_logged} more galaxy logs")

    # ---------- Supplemental galaxies from SDSS (if available) ----------
    if (not legacy_attempted) and bool(getattr(cfg, "enable_sdss", True)):
        sdss_tab = query_sdss_photoobj(center, rad, cfg)
        if sdss_tab is not None and len(sdss_tab) > 0:
            if bool(getattr(cfg, "log_sdss_colnames", False)):
                try:
                    print(f"[SDSS] photo columns={list(sdss_tab.colnames)}")
                except Exception:
                    pass
            # SDSS columns commonly available: 'ra', 'dec', 'type', 'petroRad_r', 'petroRadErr_r'
            ra_col = pick_first_existing_col(sdss_tab, ["ra", "RA", "_VRIAJ2000"])
            dec_col = pick_first_existing_col(sdss_tab, ["dec", "DEC", "_DEJ2000"])
            type_col = pick_first_existing_col(sdss_tab, ["type", "Type"])

            if ra_col and dec_col:
                sdss_sky = SkyCoord(ra=np.array(sdss_tab[ra_col]) * u.deg,
                                     dec=np.array(sdss_tab[dec_col]) * u.deg,
                                     frame="icrs")
                x_sdss, y_sdss = w.world_to_pixel(sdss_sky)
                objid_col = pick_first_existing_col(sdss_tab, ["objid", "objID", "ObjID", "ID"])  # SDSS uses objid

                gal_logged = 0
                gal_log_limit = int(getattr(cfg, "log_max_galaxies", 0) or 0)

                for i in range(len(sdss_tab)):
                    xi, yi = float(x_sdss[i]), float(y_sdss[i])
                    if not np.isfinite(xi) or not np.isfinite(yi):
                        continue

                    sc = sdss_sky[i]
                    sdss_iau = iau_coord_name("SDSS", sc, ra_precision=2, dec_precision=1)
                    if exclude_center.value > 0 and sc.separation(center) < exclude_center:
                        continue

                    if bool(getattr(cfg, "reject_bg_inside_target_footprint", False)) and (target_fp is not None):
                        xi_i = int(round(float(xi)))
                        yi_i = int(round(float(yi)))
                        if (0 <= xi_i < nx) and (0 <= yi_i < ny) and bool(target_fp[yi_i, xi_i]):
                            continue

                    # Only take SDSS photometric galaxies; skip stars
                    if type_col and sdss_tab[type_col][i] is not None:
                        t = sdss_tab[type_col][i]
                        is_gal = False
                        try:
                            is_gal = (int(t) == 3)
                        except Exception:
                            ts = t.decode() if isinstance(t, (bytes, bytearray)) else str(t)
                            is_gal = (ts.strip().upper() == "GALAXY")
                        if not is_gal:
                            continue

                    # Reject near good Gaia sources (stellar contaminants)
                    if gaia_sky_for_ps1_reject is not None and cfg.sdss_reject_if_near_gaia_arcsec > 0:
                        try:
                            idx, sep2d, _ = sc.match_to_catalog_sky(gaia_sky_for_ps1_reject)
                            if sep2d < (cfg.sdss_reject_if_near_gaia_arcsec * u.arcsec):
                                continue
                        except Exception:
                            pass

                    # Get cz/z once (used both for sizing override and for logging).
                    cz, zval, zsrc, zsep = get_best_cz_info(sc, cfg, sdss_spec_sky, sdss_spec_tab, ned_sky, ned_tab)
                    d_mpc = (cz / float(cfg.hubble_km_s_mpc)) if (cz is not None) else None

                    # High-z objects are expected to be PSF-limited in MUSE; apply this *before* any
                    # SDSS morphology-based sizing to avoid sentinel-mag overflows.
                    highz_psf_override = False
                    r_arcsec = None
                    try:
                        zsrc_is_spec = isinstance(zsrc, str) and ("SDSS(spec)" in zsrc or "NED(spec)" in zsrc)
                        if zsrc_is_spec and zval is not None and np.isfinite(zval) and float(zval) > float(cfg.highz_psf_override_zmin):
                            r_arcsec = (float(cfg.highz_psf_k_fwhm) * float(cfg.fwhm_arcsec)) + float(cfg.gaia_margin_arcsec)
                            if float(cfg.highz_psf_rmax_arcsec) > 0:
                                r_arcsec = min(float(r_arcsec), float(cfg.highz_psf_rmax_arcsec))
                            highz_psf_override = True
                    except Exception:
                        r_arcsec = None

                    if r_arcsec is None:
                        try:
                            r_arcsec = sdss_mask_radius_arcsec(sdss_tab[i], cfg)
                        except OverflowError:
                            psf = _get_col_float(sdss_tab[i], "psfMag_r")
                            mod = _get_col_float(sdss_tab[i], "modelMag_r")
                            objid = str(sdss_tab[objid_col][i]) if objid_col else "?"
                            print(f"[SDSS-OVERFLOW] objid={objid} psfMag_r={psf} modelMag_r={mod} -> using r_max")
                            r_arcsec = float(cfg.star_r_max_arcsec + cfg.gaia_margin_arcsec)

                    if r_arcsec is None or not np.isfinite(r_arcsec):
                        continue

                    r_arcsec = float(np.clip(float(r_arcsec), cfg.gal_r_min_arcsec, cfg.gal_r_max_arcsec))
                    r_pix = r_arcsec / pixscale

                    # Skip if circle does not touch the FITS image at all
                    if not circle_intersects_fov(float(xi), float(yi), float(r_pix), nx, ny):
                        continue

                    # Evidence-based masking for SDSS galaxies too
                    if cfg.require_nonvirgo_confirmation_for_galaxy_mask:
                        bg = is_definitely_background(
                            sc=sc,
                            cfg=cfg,
                            sdss_spec_sky=sdss_spec_sky,
                            sdss_spec_tab=sdss_spec_tab,
                            ned_sky=ned_sky,
                            ned_tab=ned_tab,
                        )
                        if bg is True:
                            pass
                        elif bg is None:
                            # Optional SDSS photometric fallback (parallel to PS1 fallback)
                            if not bool(getattr(cfg, "sdss_allow_photometric_fallback", False)):
                                continue

                            psf = _get_col_float(sdss_tab[i], "psfMag_r")
                            mod = _get_col_float(sdss_tab[i], "modelMag_r")
                            dmag = (float(psf) - float(mod)) if (psf is not None and mod is not None) else None
                            r90 = _get_col_float(sdss_tab[i], "petroR90_r")

                            mag_ok = (mod is not None and float(mod) < float(getattr(cfg, "sdss_fallback_rmag_max", 21.0)))
                            morph_ok = (dmag is not None and float(dmag) >= float(getattr(cfg, "sdss_fallback_dmag_min", 0.30)))
                            size_ok = (
                                r90 is not None
                                and float(r90) >= float(getattr(cfg, "sdss_fallback_petroR90_min_arcsec", 2.0))
                                and float(r90) <= float(getattr(cfg, "sdss_fallback_petroR90_max_arcsec", float("inf")))
                            )

                            if not (mag_ok and morph_ok and size_ok):
                                continue
                        else:
                            # cz says Virgo/nearby (or otherwise not definitely background) -> do not mask
                            continue
                    # Fallback: old Virgo-distance veto
                    elif ned_sky is not None and cfg.virgo_match_arcsec > 0:
                        try:
                            nidx, nsep, _ = sc.match_to_catalog_sky(ned_sky)
                            if nsep < (float(cfg.virgo_match_arcsec) * u.arcsec):
                                z = None
                                v_kms = None
                                for zname in ["Redshift", "z", "Z"]:
                                    if zname in ned_tab.colnames:
                                        try:
                                            zv = ned_tab[zname][nidx]
                                            z = None if (zv is None or getattr(zv, "mask", False)) else float(zv)
                                        except Exception:
                                            z = None
                                        break
                                for vname in ["Velocity", "cz", "Vel", "V"]:
                                    if vname in ned_tab.colnames:
                                        try:
                                            vv = ned_tab[vname][nidx]
                                            v_kms = None if (vv is None or getattr(vv, "mask", False)) else float(vv)
                                        except Exception:
                                            v_kms = None
                                        break
                                if _is_virgo_distance_from_z_or_v(cfg, z=z, v_kms=v_kms):
                                    continue
                        except Exception:
                            pass

                    # Rasterize circle
                    rasterize_circle(mask, xi, yi, r_pix, fov_mask=fov_mask)

                    objid = "?"
                    if objid_col:
                        try:
                            objid = str(sdss_tab[objid_col][i])
                        except Exception:
                            objid = "?"

                    if cfg.log_each_galaxy and (gal_log_limit <= 0 or gal_logged < gal_log_limit):
                        ra_hms, dec_dms = format_radec_hmsdms(sc, precision=2)
                        dmag = None
                        try:
                            psf = _get_col_float(sdss_tab[i], "psfMag_r")
                            mod = _get_col_float(sdss_tab[i], "modelMag_r")
                            if psf is not None and mod is not None:
                                dmag = float(psf) - float(mod)
                        except Exception:
                            dmag = None
                        morph = "pt" if sdss_is_pointlike(sdss_tab[i], cfg) else "ext"
                        morph_str = f" morph={morph}" + (f" dmag={dmag:.3f}" if dmag is not None else "")
                        if highz_psf_override:
                            morph_str += " highz_psf"
                        if cz is not None:
                            print(
                                f"[GAL] SDSS {objid} ({sdss_iau}) ra={sc.ra.deg:.6f} dec={sc.dec.deg:.6f} "
                                f"RA={ra_hms} DEC={dec_dms} "
                                f"cz={cz:.0f}km/s z={zval:.5f} D~{d_mpc:.1f}Mpc ({zsrc}) "
                                f"r_arcsec={r_arcsec:.2f}{morph_str}"
                            )
                        else:
                            print(
                                f"[GAL] SDSS {objid} ({sdss_iau}) ra={sc.ra.deg:.6f} dec={sc.dec.deg:.6f} "
                                f"RA={ra_hms} DEC={dec_dms} r_arcsec={r_arcsec:.2f}{morph_str}"
                            )
                        gal_logged += 1
                    n_gal_masked += 1
                    y_plot = (ny - 1 - yi) if use_png_bg else yi
                    gal_patches.append(Circle((xi, y_plot), r_pix, fill=False))

            if cfg.log_each_galaxy and gal_log_limit > 0 and n_gal_masked > gal_logged:
                print(f"[GAL] ... suppressed {n_gal_masked - gal_logged} more galaxy logs")

    # ---------- SDSS star fallback (only if Gaia yielded none) ----------
    if (gaia is None or len(gaia) == 0) and bool(getattr(cfg, "sdss_star_fallback", True)) and bool(getattr(cfg, "enable_sdss", True)):
        sdss_tab = query_sdss_photoobj(center, rad, cfg)
        if sdss_tab is not None and len(sdss_tab) > 0:
            ra_col = pick_first_existing_col(sdss_tab, ["ra", "RA", "_VRIAJ2000"])
            dec_col = pick_first_existing_col(sdss_tab, ["dec", "DEC", "_DEJ2000"])
            type_col = pick_first_existing_col(sdss_tab, ["type", "Type"])
            mag_col = pick_first_existing_col(sdss_tab, ["psfMag_r", "modelMag_r", "r"])
            if ra_col and dec_col and type_col:
                sdss_sky = SkyCoord(ra=np.array(sdss_tab[ra_col]) * u.deg,
                                     dec=np.array(sdss_tab[dec_col]) * u.deg,
                                     frame="icrs")
                x_sdss, y_sdss = w.world_to_pixel(sdss_sky)
                for i in range(len(sdss_tab)):
                    t = sdss_tab[type_col][i]
                    is_star = False
                    try:
                        is_star = (int(t) == 6)
                    except Exception:
                        ts = t.decode() if isinstance(t, (bytes, bytearray)) else str(t)
                        is_star = (ts.strip().upper() == "STAR")
                    if not is_star:
                        continue
                    xi, yi = float(x_sdss[i]), float(y_sdss[i])
                    if not np.isfinite(xi) or not np.isfinite(yi):
                        continue
                    sc = sdss_sky[i]
                    sdss_iau = iau_coord_name("SDSS", sc, ra_precision=2, dec_precision=1)
                    if exclude_center.value > 0 and sc.separation(center) < exclude_center:
                        continue
                    g_like_mag = _to_float_or_nan(sdss_tab[mag_col][i]) if mag_col else np.nan
                    r_arcsec = star_radius_arcsec_from_g(cfg, g_like_mag if np.isfinite(g_like_mag) else 18.0)
                    r_pix = r_arcsec / pixscale

                    # Skip if circle does not touch the FITS image at all
                    if not circle_intersects_fov(float(xi), float(yi), float(r_pix), nx, ny):
                        continue

                    if fov_mask is not None:
                        if not circle_overlaps_mask(fov_mask, float(xi), float(yi), float(r_pix)):
                            continue
                        r_mask_pix = float(r_pix)
                        if star_edge_buffer_pix is not None and float(star_edge_buffer_pix) > 0:
                            if not circle_overlaps_mask(fov_mask, float(xi), float(yi), float(star_edge_buffer_pix)):
                                continue
                            if not circle_within_mask(fov_mask, float(xi), float(yi), float(star_edge_buffer_pix)):
                                r_mask_pix = float(min(r_pix, float(star_edge_buffer_pix)))
                    else:
                        r_mask_pix = float(r_pix)

                    rasterize_circle(mask, xi, yi, r_mask_pix, fov_mask=fov_mask)
                    if cfg.log_each_star:
                        ra_hms, dec_dms = format_radec_hmsdms(sc, precision=2)
                        print(
                            f"[STAR] SDSS ({sdss_iau}) ra={sc.ra.deg:.6f} dec={sc.dec.deg:.6f} "
                            f"RA={ra_hms} DEC={dec_dms} r_arcsec={r_arcsec:.2f}"
                        )
                    n_star_masked += 1
                    y_plot = (ny - 1 - yi) if use_png_bg else yi
                    star_patches.append(Circle((xi, y_plot), r_mask_pix, fill=False))

    # ---------- Write mask FITS (nGIST expects 0=unmasked, 1=masked) ----------
    if fov_mask is not None and fov_mask.shape == mask.shape:
        mask = np.where(fov_mask, mask, 0)
    # Same spatial dims as input image/cube. Ensure extension name is MASK.
    hdr_mask = hdr.copy()
    hdr_mask["EXTNAME"] = "MASK"
    hdul_out = fits.HDUList([
        fits.PrimaryHDU(),
        fits.ImageHDU(mask.astype(np.uint8), header=hdr_mask, name="MASK"),
    ])
    hdul_out.writeto(out_mask_fits, overwrite=True)
    print(
        f"[OK] Wrote {out_mask_fits} (shape={mask.shape}, masked={(mask>0).sum()} px; "
        f"stars={n_star_masked}, galaxies={n_gal_masked})"
    )

    # ---------- Overlay PNG (pixel-locked output size) ----------
    if use_png_bg:
        assert Image is not None
        im_obj = Image.open(png_path)
        im = np.asarray(im_obj)
        H, W = im.shape[0], im.shape[1]

        # If PNG size differs from FITS, resize PNG to match FITS so all outputs align.
        if (H, W) != (ny, nx):
            print(
                f"WARNING: {png_path} shape=({H}, {W}) != FITS shape=({ny}, {nx}); resizing PNG to match FITS for overlay."
            )
            try:
                im_obj = im_obj.resize((nx, ny), resample=Image.BILINEAR)
                im = np.asarray(im_obj)
                H, W = im.shape[0], im.shape[1]
                try:
                    im_obj.save(png_path)
                    print(f"[PNG] Resized and overwrote {png_path} to match FITS size.")
                except Exception as e:
                    print(f"WARNING: failed to write resized PNG back to disk ({e}).")
            except Exception as e:
                print(f"WARNING: failed to resize PNG ({e}); using FITS background instead.")
                use_png_bg = False

    dpi = int(cfg.output_dpi)

    if use_png_bg:
        # Output exactly matches the input PNG pixel dimensions
        fig = plt.figure(figsize=(nx / dpi, ny / dpi), dpi=dpi)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.imshow(im, origin="upper", interpolation="nearest")
    else:
        # Output exactly matches the FITS pixel dimensions
        fig = plt.figure(figsize=(nx / dpi, ny / dpi), dpi=dpi)
        ax = fig.add_axes([0, 0, 1, 1])
        v = np.nanpercentile(data2d, [2, 98])
        ax.imshow(data2d, origin="upper", vmin=v[0], vmax=v[1], interpolation="nearest")

    # Fix axes limits to the image pixel grid
    ax.set_xlim(-0.5, nx - 0.5)
    ax.set_ylim(ny - 0.5, -0.5)
    ax.set_axis_off()

    if bool(getattr(cfg, "fov_draw_contour", False)) and (fov_mask_plot is not None):
        ax.contour(
            fov_mask_plot.astype(float),
            levels=[0.5],
            colors=str(getattr(cfg, "fov_contour_color", "yellow")),
            linestyles=str(getattr(cfg, "fov_contour_linestyle", ":")),
            linewidths=float(getattr(cfg, "fov_contour_linewidth", 1.0)),
        )

    if bool(getattr(cfg, "draw_target_iso_contour", False)) and (target_fp is not None):
        ax.contour(
            target_fp.astype(float),
            levels=[0.5],
            colors=["blue"],
            linestyles=["--"],
            linewidths=[1.2],
            origin="upper",
        )

    # Draw outlines (clip to axes so nothing expands the saved canvas)
    for p in star_patches:
        p.set_edgecolor("green")
        p.set_linewidth(1.2)
        p.set_clip_on(True)
        ax.add_patch(p)

    for p in gal_patches:
        p.set_edgecolor("brown")
        p.set_linewidth(1.2)
        p.set_clip_on(True)
        ax.add_patch(p)

    # IMPORTANT: no bbox_inches="tight" here; it changes output size
    fig.savefig(out_overlay_png, dpi=dpi)
    plt.close(fig)
    if Image is not None:
        try:
            out_im = Image.open(out_overlay_png)
            if out_im.size != (nx, ny):
                print(
                    f"WARNING: {out_overlay_png} size={out_im.size} != FITS size=({nx}, {ny}); resizing output PNG."
                )
                out_im = out_im.resize((nx, ny), resample=Image.BILINEAR)
                out_im.save(out_overlay_png)
        except Exception as e:
            print(f"WARNING: failed to verify/resize output PNG size ({e}).")
    print(f"[OK] Wrote {out_overlay_png}")

    return out_mask_fits, out_overlay_png


def main():
    cfg = Config()

    # Force "all on-screen log" on (and unlimited per-object logging)
    cfg.log_each_star = True
    cfg.log_each_galaxy = True
    cfg.log_sdss_colnames = True
    cfg.log_max_galaxies = 0

    # Default: do not write DS9 region files (avoid clutter).
    cfg.legacy_write_ds9_regions = False

    # Tee ALL stdout/stderr to a log file as well as the console
    log_path = "v3tk_masking.log"

    class Tee:
        def __init__(self, *streams):
            self.streams = streams

        def write(self, data):
            for s in self.streams:
                s.write(data)

        def flush(self):
            for s in self.streams:
                s.flush()

        def isatty(self):
            return any(getattr(s, "isatty", lambda: False)() for s in self.streams)

    logf = open(log_path, "w", buffering=1)
    orig_stdout, orig_stderr = sys.stdout, sys.stderr
    sys.stdout = Tee(orig_stdout, logf)
    sys.stderr = Tee(orig_stderr, logf)
    try:
        print(f"=== v3tk masking run start: {datetime.now().isoformat(timespec='seconds')} ===")

        t_total0 = perf_counter()

        if len(sys.argv) > 1:
            rfits_list: list[str] = []
            for pat in sys.argv[1:]:
                rfits_list.extend(sorted(glob.glob(pat)))
            # de-dup while preserving order
            seen = set()
            rfits_list = [p for p in rfits_list if not (p in seen or seen.add(p))]
        else:
            rfits_list = sorted(glob.glob("*_DATACUBE*_VRI.fits"))
        if len(rfits_list) == 0:
            raise SystemExit("No *_DATACUBE*_VRI.fits files found in the current directory.")

        for rfits in rfits_list:
            t0 = perf_counter()
            try:
                build_masks_for_one(rfits, cfg)
            except Exception as e:
                print(f"[FAIL] {rfits}: {e}")
                traceback.print_exc()
            finally:
                dt = perf_counter() - t0
                print(f"[TIME] {rfits} runtime_s={dt:.2f}")

        t_total = perf_counter() - t_total0
        print(f"[TIME] total runtime_s={t_total:.2f}")

        print(f"=== v3tk masking run end:   {datetime.now().isoformat(timespec='seconds')} ===")
    finally:
        sys.stdout = orig_stdout
        sys.stderr = orig_stderr
        logf.close()


if __name__ == "__main__":
    main()
