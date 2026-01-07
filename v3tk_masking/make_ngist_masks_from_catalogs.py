#!/usr/bin/env python
"""
Create nGIST-compatible spatial masks from catalogs.

Inputs expected in the working directory (per galaxy XXX):
  - XXX_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk_R.fits   (2D image, WCS in HDU0)
  - XXX_combined_R.png                               (optional, same pixel grid)

Outputs:
  - XXX_mask.fits                 (0=unmasked, 1=masked; same spatial dims as FITS)
  - XXX_combined_R_mask.png       (overlay: green=Gaia stars, brown=galaxy-like)

Notes:
  - Gaia query: foreground stars (uses Gaia DR3 via astroquery.gaia)
  - Galaxies: Pan-STARRS (MAST) if Dec >= -30 deg; otherwise tries SkyMapper DR4 (Vizier).
  - You will probably still need to manually add a few artefact regions.

Requirements:
  pip install astropy astroquery matplotlib numpy pillow
"""

from __future__ import annotations

import glob
import os
import sys
from dataclasses import dataclass
import numpy as np

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs.utils import proj_plane_pixel_scales

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse

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


@dataclass
class Config:
    fwhm_arcsec: float = 1.0                # typical MUSE seeing FWHM, tune if needed
    gaia_gmag_max: float = 21.0             # ignore very faint Gaia sources
    gaia_margin_arcsec: float = 1.0         # extra padding on radii (registration / wings)

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
    gaia_parallax_snr_min: float = 3.0
    gaia_parallax_min_mas: float = 0.2
    gaia_pm_snr_min: float = 5.0
    gaia_pm_min_masyr: float = 2.0

    # NOTE: We intentionally do NOT exclude the inner galaxy by default, because
    # you may still want to mask true foreground stars/background galaxies there.
    # Keep this at 0 unless you explicitly want an inner "no-mask" zone.
    exclude_center_arcsec: float = 0.0

    # PS1(VizieR) quality cuts (to avoid spurious detections / bad photometry).
    # These are intentionally conservative; relax if you miss too many background galaxies.
    # VERY STRICT preset (few galaxies, minimal false positives)
    ps1_min_Nr: int = 3
    ps1_e_mag_max: float = 0.10

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
    ps1_ext_thresh: float = 0.40
    # Upper bound to guard against pathological photometry (blends/saturation can
    # produce huge PSF-Kron differences that are not real galaxies).
    ps1_ext_max: float = 1.5
    ps1_rmag_max: float = 22.0              # ignore very faint objects (optional)
    ps1_require_ri_extended: bool = True    # if i-band mags exist, require extendedness in both r and i
    gal_r_min_arcsec: float = 2.0
    gal_r_max_arcsec: float = 30.0

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

    use_png_background: bool = True         # if False, uses FITS as background for overlay
    output_dpi: int = 200


def safe_base_id(rfits_path: str) -> str:
    # From XXX_DATACUBE..._R.fits => XXX
    bn = os.path.basename(rfits_path)
    return bn.split("_DATACUBE")[0]


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
    r = cfg.star_r_ref_arcsec * (10.0 ** (-0.2 * (gmag - cfg.star_g_ref)))
    r = max(cfg.star_r_min_arcsec, r)
    r = min(cfg.star_r_max_arcsec, r)
    # also ensure at least ~1 FWHM
    r = max(r, 1.0 * cfg.fwhm_arcsec)
    # margin
    r += cfg.gaia_margin_arcsec
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
        job = Gaia.launch_job_async(query, dump_to_file=False)
        return job.get_results()


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
    v = Vizier(columns=["**"], row_limit=-1)
    try:
        res = v.query_region(center, radius=radius, catalog="II/349/ps1")
        if len(res) == 0:
            return None
        tab = res[0]
    except Exception as e:
        print(f"WARNING: Pan-STARRS (VizieR) query failed; skipping. ({e})")
        return None

    # Optional magnitude cut to keep table sizes reasonable
    if "rmag" in tab.colnames:
        try:
            tab = tab[tab["rmag"] < cfg.ps1_rmag_max]
        except Exception:
            pass
    return tab


def query_skymapper(center: SkyCoord, radius: u.Quantity):
    if Vizier is None:
        print("WARNING: astroquery.vizier not available; skipping SkyMapper query.")
        return None

    v = Vizier(columns=["**"], row_limit=-1)
    # SkyMapper DR4 Vizier table
    # II/379/smssdr4 (very large table; cone searches are OK but can be slower)
    res = v.query_region(center, radius=radius, catalog="II/379/smssdr4")
    if len(res) == 0:
        return None
    return res[0]


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


def rasterize_circle(mask: np.ndarray, xi: float, yi: float, r_pix: float) -> None:
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
    mask[y0 : y1 + 1, x0 : x1 + 1][rr2 <= (r_pix**2)] = 1


def rasterize_ellipse(mask: np.ndarray, xi: float, yi: float, a_pix: float, b_pix: float, angle_deg: float) -> None:
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
    mask[y0 : y1 + 1, x0 : x1 + 1][inside] = 1


def build_masks_for_one(rfits_path: str, cfg: Config):
    base = safe_base_id(rfits_path)
    png_path = f"{base}_combined_R.png"
    out_mask_fits = f"{base}_mask.fits"
    out_overlay_png = f"{base}_combined_R_mask.png"

    use_png_bg = bool(cfg.use_png_background and os.path.exists(png_path) and Image is not None)

    data2d, w, hdr, nx, ny = load_r_image_and_wcs(rfits_path)
    pixscale = pixel_scale_arcsec(w)
    center, rad = fov_center_and_radius(w, nx, ny)
    # small padding so we don't miss edge objects
    rad = rad + (10.0 * u.arcsec)

    mask = np.zeros((ny, nx), dtype=np.uint8)
    exclude_center = cfg.exclude_center_arcsec * u.arcsec

    # ---------- Gaia stars ----------
    gaia = query_gaia_sources(center, rad, cfg)
    star_patches = []
    n_star_masked = 0
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

            # rasterize into mask
            rasterize_circle(mask, xi, yi, r_pix)

            sid = row["source_id"] if "source_id" in row.colnames else "?"
            if cfg.log_each_star:
                print(
                    f"[STAR] Gaia {sid} ra={sc.ra.deg:.6f} dec={sc.dec.deg:.6f} "
                    f"G={float(gi):.2f} r_arcsec={r_arcsec:.2f}"
                )
            n_star_masked += 1

            # The provided PNG background images are typically rendered with a
            # different vertical origin than raw FITS array indices.
            y_plot = (ny - 1 - yi) if use_png_bg else yi
            star_patches.append(Circle((xi, y_plot), r_pix, fill=False))

    # ---------- Background galaxies (Pan-STARRS preferred) ----------
    gal_patches = []
    gal_tab = None
    n_gal_masked = 0

    # Optional: NED crossmatch for Virgo-distance veto (only used when enabled and available)
    ned_tab = None
    ned_sky = None
    if bool(getattr(cfg, "enable_virgo_distance_veto", False)):
        ned_tab = query_ned_redshifts(center, rad)
        try:
            if ned_tab is not None and len(ned_tab) > 0 and "RA" in ned_tab.colnames and "DEC" in ned_tab.colnames:
                ned_sky = SkyCoord(ra=np.array(ned_tab["RA"]) * u.deg, dec=np.array(ned_tab["DEC"]) * u.deg, frame="icrs")
        except Exception:
            ned_sky = None

    if center.dec.deg >= -30:
        gal_tab = query_ps1_galaxy_like(center, rad, cfg)

    # If PS1 failed/empty and SkyMapper is available, try SkyMapper (useful mainly in the south)
    if (gal_tab is None or len(gal_tab) == 0) and center.dec.deg <= +28:
        gal_tab = query_skymapper(center, rad)

    if gal_tab is not None and len(gal_tab) > 0:
        # Dynamic column picking (PS1 and SkyMapper differ)
        ra_col = pick_first_existing_col(gal_tab, ["raMean", "RAJ2000", "ra", "_RAJ2000"])
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
                if exclude_center.value > 0 and sc.separation(center) < exclude_center:
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

                # If this candidate is likely at Virgo distance (using NED z/velocity when available),
                # do NOT mask it.
                if ned_sky is not None and cfg.virgo_match_arcsec > 0:
                    try:
                        nidx, nsep, _ = sc.match_to_catalog_sky(ned_sky)
                        if nsep < (float(cfg.virgo_match_arcsec) * u.arcsec):
                            z = None
                            v_kms = None
                            # NED column names vary; try common ones.
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

                    # Plot coords (only affects overlay, not the FITS mask)
                    y_plot = (ny - 1 - yi) if use_png_bg else yi
                    angle_plot = (-angle) if use_png_bg else angle

                    # rasterize ellipse (local cutout)
                    rasterize_ellipse(mask, xi, yi, a_pix, b_pix, angle)

                    if cfg.log_each_galaxy and (gal_log_limit <= 0 or gal_logged < gal_log_limit):
                        qv = gal_tab["Qual"][i] if "Qual" in gal_tab.colnames else "?"
                        print(
                            f"[GAL] {objid} ra={sc.ra.deg:.6f} dec={sc.dec.deg:.6f} "
                            f"a_arcsec={a_arcsec:.2f} b_arcsec={b_arcsec:.2f} pa={angle:.1f} Qual={qv}"
                        )
                        gal_logged += 1
                    n_gal_masked += 1

                    gal_patches.append(Ellipse((xi, y_plot), 2 * a_pix, 2 * b_pix, angle=angle_plot, fill=False))
                else:
                    # rasterize circle
                    rasterize_circle(mask, xi, yi, r_pix)
                    if cfg.log_each_galaxy and (gal_log_limit <= 0 or gal_logged < gal_log_limit):
                        qv = gal_tab["Qual"][i] if "Qual" in gal_tab.colnames else "?"
                        # if PSF/Kron mags available, print extendedness for debugging
                        ext = "?"
                        if psf_col and kron_col:
                            psf = _to_float_or_nan(gal_tab[psf_col][i])
                            kron = _to_float_or_nan(gal_tab[kron_col][i])
                            if np.isfinite(psf) and np.isfinite(kron):
                                ext = f"{(psf-kron):.3f}"
                        nr = gal_tab["Nr"][i] if "Nr" in gal_tab.colnames else "?"
                        print(
                            f"[GAL] {objid} ra={sc.ra.deg:.6f} dec={sc.dec.deg:.6f} "
                            f"r_arcsec={r_arcsec:.2f} ext={ext} Nr={nr} Qual={qv}"
                        )
                        gal_logged += 1
                    n_gal_masked += 1
                    y_plot = (ny - 1 - yi) if use_png_bg else yi
                    gal_patches.append(Circle((xi, y_plot), r_pix, fill=False))

            if cfg.log_each_galaxy and gal_log_limit > 0 and n_gal_masked > gal_logged:
                print(f"[GAL] ... suppressed {n_gal_masked - gal_logged} more galaxy logs")

    # ---------- Write mask FITS (nGIST expects 0=unmasked, 1=masked) ----------
    # Same spatial dims as input image/cube. :contentReference[oaicite:1]{index=1}
    fits.writeto(out_mask_fits, mask.astype(np.uint8), header=hdr, overwrite=True)
    print(
        f"[OK] Wrote {out_mask_fits} (shape={mask.shape}, masked={(mask>0).sum()} px; "
        f"stars={n_star_masked}, galaxies={n_gal_masked})"
    )

    # ---------- Overlay PNG ----------
    fig, ax = plt.subplots(figsize=(6, 6))

    if use_png_bg:
        assert Image is not None
        im = np.asarray(Image.open(png_path))
        ax.imshow(im, origin="upper")
    else:
        # Use FITS as background (more robust; ensures same orientation as WCS pixel coords)
        v = np.nanpercentile(data2d, [2, 98])
        ax.imshow(data2d, origin="upper", vmin=v[0], vmax=v[1])

    # Draw outlines
    for p in star_patches:
        p.set_edgecolor("green")
        p.set_linewidth(1.2)
        ax.add_patch(p)

    for p in gal_patches:
        p.set_edgecolor("brown")
        p.set_linewidth(1.2)
        ax.add_patch(p)

    ax.set_axis_off()
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    fig.savefig(out_overlay_png, dpi=cfg.output_dpi, bbox_inches="tight", pad_inches=0)
    plt.close(fig)
    print(f"[OK] Wrote {out_overlay_png}")

    return out_mask_fits, out_overlay_png


def main():
    cfg = Config()

    if len(sys.argv) > 1:
        rfits_list: list[str] = []
        for pat in sys.argv[1:]:
            rfits_list.extend(sorted(glob.glob(pat)))
        # de-dup while preserving order
        seen = set()
        rfits_list = [p for p in rfits_list if not (p in seen or seen.add(p))]
    else:
        rfits_list = sorted(glob.glob("*_DATACUBE*_R.fits"))
    if len(rfits_list) == 0:
        raise SystemExit("No *_DATACUBE*_R.fits files found in the current directory.")

    for rfits in rfits_list:
        try:
            build_masks_for_one(rfits, cfg)
        except Exception as e:
            print(f"[FAIL] {rfits}: {e}")


if __name__ == "__main__":
    main()
