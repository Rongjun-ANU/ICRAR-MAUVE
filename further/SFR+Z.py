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

Changes (2025-09-17)
-----------------------
* Enhanced inclination correction methodology from simple cos(θ) to account for b/a factor.
* More physically accurate correction accounting for finite disc thickness rather than infinitely thin discs.
* Implemented b/a = sqrt((1-q₀²)*cos²(i) + q₀²) correction where q₀ = 0.2 (intrinsic disc thickness).
* Updated logging to report inclination angle, cos(θ), b/a factor, and adopted q₀ parameter.

Changes (2025-09-23)
-----------------------
* Implemented dual BPT classification system for comparative analysis:
  - Classification 1 (Liberal): SF = HII + Composite regions (existing approach)
  - Classification 2 (Conservative): SF = HII regions only (new conservative approach)
* Enhanced choose_BPT() function with classification parameter (classification=1/2):
  - Returns both log and regular SFR surface density maps for each classification
  - Maintains complete parallel outputs for direct comparison between approaches
* Added comprehensive Classification 2 outputs:
  - Complete mask hierarchy: mask_classified2_N2_HII, mask_classified2_N2_Comp_AGN, etc.
  - HII-specific FITS outputs: LOGSFR_SURFACE_DENSITY_HII, Halpha_SFR_corr_HII, etc.
  - Full metallicity analysis for HII regions: O_H_D16_HII, O_H_PG16_HII, all M13/PP04/C20 methods
  - Metallicity error maps for HII regions: O_H_O3N2_C20_HII_ERR, O_H_COMBINED_C20_HII_ERR, etc.
  - Line flux maps for HII regions: HB4861_FLUX_corr_HII, HA6562_FLUX_corr_HII, etc.
* Added missing LOGSFR_SURFACE_DENSITY_UNCLASSIFIED2 for complete Classification 2 coverage
* Enhanced terminal output with parallel reporting for both classifications:
  - Total Halpha SFR from SF region vs HII region comparison
  - Complete metallicity totals for both SF and HII regions using all calibration methods
* Maintained backwards compatibility while providing new conservative analysis option
* Updated variable naming convention: '_classified' → '_classified1', '_unclassified' → '_unclassified1'

Changes (2025-09-24 & 2025-09-25)
-----------------------
* Added N2S2-N06 metallicity calibration (Nakajima & Ouchi 2014) positioned between PG16 and O3N2-M13:
  - Implemented calculate_n2s2_n06_metallicity() function with cubic polynomial equation
  - Equation: log([N II]λ6584/([S II]λ6716+λ6731)) = -0.25214 + 0.74100·x + 0.58181·x² + 0.17963·x³
  - Uses numpy.roots() for accurate 3rd-order polynomial root solving
* Complete integration across dual BPT classification systems:
  - Added O_H_N2S2_N06_SF and O_H_N2S2_N06_HII maps for both classifications
  - Integrated total region calculations: O_H_N2S2_N06_SF_total, O_H_N2S2_N06_HII_total, O_H_N2S2_N06_total
  - Added FITS HDU extensions with descriptive headers and metadata
  - Enhanced terminal reporting with N2S2-N06 metallicity summaries for all regions
* Robust implementation features:
  - Comprehensive error handling for invalid flux data (NaN, negative, zero values)
  - Maintains full backward compatibility with existing analysis pipeline

Changes (2026-01-15)
-----------------------
* Added Milky Way (CCM89; Cardelli, Clayton & Mathis 1989) extinction curve as
  an optional k(λ)=A(λ)/E(B−V) alternative.
* Kept Calzetti (2000) available via `extinction_k(..., law="calzetti")`.

Changes (2026-03-31)
-----------------------
* Added `--fallback-root` for secondary input lookup.
* Gas, kinematic, and extended input FITS files are now resolved per file,
  checking the primary root first and then the fallback root.
* Removed the hardwired local gas-path override so mixed-root runs work
  consistently when some inputs are on CANFAR and others are local.

Changes (2026-04-14)
-----------------------
* For SF/HII/total integrated properties, totals now follow the same integrated flow
    as the total-region block: sum masked raw line maps first, then apply one
    integrated Balmer-decrement correction, then derive metallicities from
    integrated corrected line fluxes.

Changes (2026-04-26)
-----------------------
* Added flexible input filename resolution for gas, kinematic, and extended
    binning products, including upper/lower-case variants and `_bin_maps`
    aliases, while preserving normalized `_extended` output naming.
* Propagated emission-line flux uncertainties through the Balmer-decrement dust
    correction before BPT and C20 metallicity error calculations.
* Fixed C20 metallicity-error masking, included the configured calibration
    fitting scatter in C20 errors, and corrected RS32 error propagation in
    linear-ratio space.
* Fixed non-log classified SFR maps to store inclination-corrected SFR surface
    density, matching their FITS units.

Changes (2026-04-27)
-----------------------
* Non-cube product inputs are now resolved as a matched optional-keyword group,
    e.g. `{galaxy}_{keyword}_gas_BIN_maps.fits`,
    `{galaxy}_{keyword}_SPATIAL_BINNING_maps_extended.fits`, and
    `{galaxy}_{keyword}_KIN_maps_extended.fits`.

Changes (2026-05-24)
-----------------------
* Fixed Combined-C20 so it is an inverse-variance weighted combination of all
    finite C20 calibrations, with method-to-method scatter added to the formal
    uncertainty. The diagnostic method map now records the dominant-weight
    contributor instead of defining the combined metallicity as one selected
    calibration.
* Added independent `NII_BPT` and `SII_BPT` HDUs using dust-corrected line
    fluxes. Both maps use -1=unknown/non-detection and 0=unclassified.
    `NII_BPT` uses 1=HII, 2=Comp, 3=AGN; `SII_BPT` uses
    1=HII, 2=LINER, 3=Seyfert.

Changes (2026-05-25)
-----------------------
* Updated the default project layout for the second MAUVE project:
  - raw v3tk gas products are read from `v3tk_v7.6.8/{galaxy}` under the working directory.
  - the Mass.py input is now `{galaxy}_spatial_binning_maps_further.fits` in the same galaxy folder.
* New SFR/metallicity products now use the `_further` suffix and are written
  back into each galaxy product folder, e.g. `v3tk_v7.6.8/{galaxy}/{galaxy}_gas_bin_maps_further.fits`.
* Replaced the old `FOREGROUND_STAR` HDU masking dependency with the nGIST
  spatial mask file `{galaxy}_mask.fits`. For `MASKFILE[MASK]`, row-major
  `MASK == 0` is treated as usable/unmasked and `MASK == 1` is excluded,
  matching the finite footprint in the current v3tk products.

Changes (2026-05-26)
-----------------------
* Added PyNeb-based electron-density and auroral-line metallicity products:
  - `NE_SII` from the [S II] 6716/6730 density ratio, with finite values below
    20 cm^-3 clamped to the Brazzini+2024 low-density value.
  - HII-gated `TE_NII_HII`, `TE_SIII_BR24_HII`, `TE_OIII_BR24_HII`, and
    `TE_OIII_NII_CHAIN_HII` temperature maps with first-order uncertainties.
  - HII-gated `O_H_BR24_DIRECT_HII`, `O_H_NII_OII7325_HII`, and
    `O_H_NII_K25_HII` oxygen-abundance maps with propagated uncertainties.
* Extra auroral/red-line inputs are optional per FITS extension, but each Te
  method only evaluates where all of its required lines pass the configured
  detection cuts; multi-line complexes such as [O II] 7325 and [O III] 4959+5007
  use summed-flux detection gates. Missing required lines produce empty/NaN
  method maps.
* The Kreckel+2025/Mendez-Delgado+2023 NII-only proxy is masked outside the
  published 8000-13000 K Te([N II]) calibration range to avoid extrapolated
  low-metallicity artifacts from very hot/noisy auroral detections.
* Te products now fail fast when enabled and the active Python environment
  cannot import PyNeb, preventing silent science-product degradation.

Changes (2026-05-26, density/log refinement)
-----------------------
* Added `NE_SII_ALL`, which keeps measured `NE_SII` where available but assigns
  20 cm^-3 to spatially usable pixels that fail the [S II] density detection
  gate. Spatially masked pixels remain NaN.
* The fast [S II] density lookup now falls back to exact PyNeb `getTemDen` for
  suspicious high-density ratios, while keeping the lookup path for ordinary
  low-density pixels.
* Te logs now report the fixed-20 density count as both a count and fraction.

Changes (2026-05-26, SF Te products and integrated Te logs)
-----------------------
* Added SF-region versions of the PyNeb Te/ionic-abundance maps, matching the
  existing HII-region products but gated by the Classification-1 SF mask.
* Added SF and HII valid-mask HDUs for the new Te metallicity methods.
* Added integrated PyNeb Te-method log summaries for SF and HII regions. These
  summaries sum raw line fluxes over the method-valid pixels, apply one
  integrated Balmer-decrement correction, then run the same PyNeb equations as
  the spaxel maps.

Changes (2026-05-26, density lookup cache and HII valid-mask names)
-----------------------
* Renamed the HII Te-method validity-mask HDUs to `*_VALID_HII`, so they are
  explicit counterparts to the existing `*_VALID_SF` maps.
* Added a reusable PyNeb [S II] density lookup cache in the current working
  directory. The script loads the `.npz` table when present, otherwise creates
  it and writes a companion plasma-colormap `.png` heatmap.

Changes (2026-05-26, multi-temperature density lookup plot)
-----------------------
* Expanded the [S II] density lookup cache to store temperature rows from
  8000 K to 13000 K in 1000 K steps. The density maps still use the 10000 K
  row, matching the previous calculation.
* Replaced the lookup visualization with a line plot of electron density on the
  x-axis and [S II] 6716/6731 ratio on the y-axis, with one plasma-colored line
  per temperature.

Changes (2026-05-26, corrected-line schema and CASA-friendly units)
-----------------------
* Emission lines are now discovered from every available `*_FLUX` /
  `*_FLUX_ERR` pair in the input gas-line FITS, so dust-corrected flux and
  error products are written for all fitted lines, not only the legacy strong
  lines and Te lines.
* Added HII and SF versions of every corrected line-flux and line-error map,
  ordered as immediate HII/SF pairs in the written FITS schema.
* Reordered the appended products into EBV, corrected lines, BPT/SFR,
  strong-line metallicities, and then PyNeb density/Te products.
* Abundance, ratio, class, method, mask, and logarithmic diagnostic maps now
  use dimensionless `BUNIT = 1`, with the physical/logarithmic meaning kept in
  FITS comments to avoid CASA/CARTA unknown-unit warnings.
* The output FITS now preserves every original input gas-map HDU first; all new
  post-processing products are appended after the original input schema.

Changes (2026-05-26, Calzetti default restored)
-----------------------
* Restored Calzetti (2000) as the default Balmer-decrement extinction curve for
  consistency with the original MAUVE SFR+Z products. CCM89 remains available as
  an explicit optional law through `extinction_k(..., law="mw")`.

Changes (2026-05-28)
-----------------------
* FITS input lookup now accepts gzip-compressed `.fits.gz` counterparts wherever
  a `.fits` path or glob is requested, while keeping new `_further` outputs as
  normal `.fits` files.

Changes (2026-06-29)
-----------------------
* Added `--product-subdir` so wrappers can run either `v3tk_v7.6.8` or
  `v3tk_v7.6.8_7000` under the selected products root and write outputs back to
  the matching run folder.

Changes (2026-08-23)
-----------------------
* Added pyqz 0.8.4.0 MAPPINGS-V model-grid products for HII and SF regions:
  oxygen abundance, LogQ, derived log U, KDE uncertainties, raw pyqz flags,
  off-grid fractions, validity maps, and separately inferred integrated values.
* Added NebulaBayes 0.9.9 MAPPINGS-5.1 HII-grid products for HII and SF
  regions: oxygen abundance, log U, log(P/k), asymmetric 68-percent credible
  intervals, reduced chi-square, posterior/QC flags, validity maps, and
  separately inferred integrated values.
* Both methods use the existing dust-corrected six-line flux/error maps and
  the existing HII/SF masks, with one inference per unique spatial BINID and no
  additional S/N cut. Integrated model values use one common raw six-line
  aperture, one integrated Balmer-decrement correction, and are log-only.
* The adjacent `model_grid_compat.py` compatibility script is required when
  either model-grid method is enabled. It supplies exact-version runtime fixes
  for current Python/NumPy/SciPy/Matplotlib without modifying site-packages.

Changes (2026-08-24)
-----------------------
* Added the Peng et al. (2026) implementation of the Ji & Yan (2022; JY22)
  N2/S2/R3 Bayesian metallicity calibration for the existing HII and SF masks.
* JY22 starts from the raw six-line fluxes and independent raw errors, applies
  the configured Balmer-decrement correction internally, and propagates the
  complete first-order N2/S2/R3 covariance. No MaNGA-specific 1.25 error
  inflation, additional S/N cut, EW cut, or three-dimensional selection is used.
* Added posterior-mean O/H and log U maps, marginal equal-tailed 16th/84th
  percentiles, minimum chi-square, QC flags, validity maps, and log-only
  integrated HII/SF results from the released Peng2026 grid.
"""

# ------------------------------------------------------------------
# User Configuration Parameters
# ------------------------------------------------------------------

# Inclination correction toggle
# Set to True to apply cos(θ) inclination correction, False to disable
apply_inclination_correction = True

# Fixed distance scale adopted for consistency with previous MAUVE papers.
DISTANCE_MPC = 16.5
DISTANCE_REFERENCE = "Fixed MAUVE paper distance scale"

# Kennicutt & Evans (2012) Hα SFR coefficient on a Kroupa IMF, converted
# to Chabrier using the Salpeter-relative factors noted in that review.
SFR_HA_KROUPA_COEFF = 5.3e-42
SALPETER_TO_KROUPA = 0.67
SALPETER_TO_CHABRIER = 0.63
KROUPA_TO_CHABRIER = SALPETER_TO_CHABRIER / SALPETER_TO_KROUPA
SFR_HA_CHABRIER_COEFF = SFR_HA_KROUPA_COEFF * KROUPA_TO_CHABRIER

# Extinction-law configuration for k(λ)=A(λ)/E(B−V)
# Supported: "calzetti" (Calzetti 2000; default), "mw" (CCM89 Milky Way)
extinction_law = "calzetti"
mw_rv = 3.1  # Used only when extinction_law is "mw"/"ccm89".

# PyNeb-based Te products. Set False only if you deliberately want to run the
# legacy strong-line products without electron-density/direct-method outputs.
ENABLE_TE_METALLICITY_PRODUCTS = True

# HII-region photoionisation-grid products. model_grid_compat.py is required
# only by the legacy pyqz and NebulaBayes packages; JY22 is implemented by the
# adjacent pure diagnostic helper and the released Peng2026 FITS grid.
ENABLE_PYQZ_METALLICITY_PRODUCTS = True
ENABLE_NEBULABAYES_METALLICITY_PRODUCTS = True
ENABLE_JY22_METALLICITY_PRODUCTS = True

PYQZ_DIAGNOSTIC = "[NII]/[SII]+;[OIII]/[SII]+"
PYQZ_QZS = ["LogQ", "gas[O]+12"]
PYQZ_PK = 5.0
PYQZ_STRUCT = "pp"
PYQZ_KAPPA = float("inf")
PYQZ_SAMPLING = 2
PYQZ_ERROR_PDF = "normal"
PYQZ_SRS = 800
PYQZ_KDE_METHOD = "multiv"
PYQZ_KDE_QZ_SAMPLING = 101j
PYQZ_KDE_DO_SINGLES = True
PYQZ_FLAG_LEVEL = 2.0
PYQZ_RANDOM_SEED = 20260823
PYQZ_NPROC = 1

NB_INTERPD_GRID_SHAPE = (40, 20, 160)
NB_INTERP_ORDER = 1
NB_GRID_ERROR = 0.10
NB_NORM_LINE = "Hbeta"
NB_DEREDDEN = False
NB_PRIOR = "Uniform"
NB_LIKELIHOOD_LINES = None

JY22_GRID_RELATIVE_PATH = (
    "Peng2026/photoionization_models/photoionization_grid_interpolated.fits"
)
JY22_EXPECTED_GRID_SHA256 = (
    "d7a219b60c9a1ea8b29339988c84b1832028b28b3c417d1c3c58b420831eb38a"
)
JY22_LOG_U_MIN = -4.0
JY22_LOG_U_MAX = -0.5
JY22_SOLAR_O_H = 8.69
JY22_INTRINSIC_BALMER_RATIO = 2.86
JY22_ERROR_INFLATION = 1.0  # Do not transfer the MaNGA-specific factor 1.25.

SII_DENSITY_LOOKUP_TEMPERATURES = tuple(range(8000, 13001, 1000))
SII_DENSITY_LOOKUP_BASENAME = (
    "pyneb_sii_6716_6731_density_lookup_te8000_13000_step1000_ne20_100000_n4096"
)
SII_DENSITY_LOOKUP_NPZ = f"{SII_DENSITY_LOOKUP_BASENAME}.npz"
SII_DENSITY_LOOKUP_PNG = f"{SII_DENSITY_LOOKUP_BASENAME}.png"
SII_DENSITY_LOOKUP_CACHE = {}

# ------------------------------------------------------------------
# 0.  Command-line interface  (exactly as requested)
# ------------------------------------------------------------------

import argparse, logging, os, re, sys, time
from importlib import metadata as importlib_metadata
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

JY22_GRID_PATH = (
    Path(__file__).resolve().parent / JY22_GRID_RELATIVE_PATH
).resolve()

if (
    ENABLE_TE_METALLICITY_PRODUCTS
    or ENABLE_PYQZ_METALLICITY_PRODUCTS
    or ENABLE_NEBULABAYES_METALLICITY_PRODUCTS
):
    os.environ.setdefault(
        "MPLCONFIGDIR",
        str(Path(os.environ.get("TMPDIR", "/tmp")) / "matplotlib"),
    )

if (
    ENABLE_PYQZ_METALLICITY_PRODUCTS
    or ENABLE_NEBULABAYES_METALLICITY_PRODUCTS
    or ENABLE_JY22_METALLICITY_PRODUCTS
):
    try:
        from model_grid_diagnostics import (
            JY22_INTEGER_FIELDS,
            MODEL_LINE_BASES,
            NB_INTEGER_FIELDS,
            NB_LINE_LIST,
            NOT_EVALUATED_FLAG,
            PYQZ_INTEGER_FIELDS,
            BinSpectra,
            broadcast_bin_results,
            build_jy22_ratio_spectra,
            build_model_input_valid_mask,
            empty_jy22_result,
            empty_nebulabayes_result,
            empty_pyqz_result,
            extract_unique_bin_spectra,
            integrate_line_maps,
            load_jy22_grid,
            ordered_model_hdu_names,
            run_jy22_spectra,
            run_nebulabayes_spectra,
            run_pyqz_spectra,
        )
    except ImportError as exc:
        raise RuntimeError(
            "Enabled model-grid products require the adjacent "
            "model_grid_diagnostics.py file."
        ) from exc

if (
    ENABLE_PYQZ_METALLICITY_PRODUCTS
    or ENABLE_NEBULABAYES_METALLICITY_PRODUCTS
):
    try:
        from model_grid_compat import (
            load_nebulabayes_hii_grid,
            load_pyqz_compat,
            make_nebulabayes_hii_model,
        )
    except ImportError as exc:
        raise RuntimeError(
            "pyqz/NebulaBayes products require the adjacent "
            "model_grid_compat.py exact-version runtime layer."
        ) from exc

if ENABLE_TE_METALLICITY_PRODUCTS:
    try:
        import pyneb as pn
    except ImportError as exc:
        raise RuntimeError(
            "PyNeb-based Te metallicity products are enabled, but the active "
            f"Python interpreter cannot import pyneb: {sys.executable}. "
            "Install PyNeb into this interpreter, or set "
            "ENABLE_TE_METALLICITY_PRODUCTS = False to run only legacy "
            "strong-line products."
        ) from exc
else:
    pn = None

# ------------------------------------------------------------------
# Helper function for inclination correction
# ------------------------------------------------------------------

def _unique_paths(*paths: Path | None) -> list[Path]:
    unique: list[Path] = []
    seen: set[str] = set()

    for path in paths:
        if path is None:
            continue
        resolved = path.expanduser().resolve()
        if str(resolved) in seen:
            continue
        seen.add(str(resolved))
        unique.append(resolved)

    return unique


def _has_glob_chars(path: Path) -> bool:
    return any(char in str(path) for char in "*?[")


def _fits_path_alternates(path: Path) -> list[Path]:
    name = path.name
    lower_name = name.lower()
    if lower_name.endswith(".fits.gz"):
        return [path, path.with_name(name[:-3])]
    if lower_name.endswith(".fits"):
        return [path, path.with_name(f"{name}.gz")]
    return [path]


def _glob_candidate_patterns(path: Path) -> list[Path]:
    return _fits_path_alternates(path)


def _candidate_paths(*paths: Path | None) -> list[Path]:
    candidates: list[Path] = []
    seen: set[str] = set()

    for path in paths:
        if path is None:
            continue

        expanded = path.expanduser()
        if _has_glob_chars(expanded):
            path_candidates: list[Path] = []
            for pattern in _glob_candidate_patterns(expanded):
                path_candidates.extend(
                    sorted(
                        pattern.parent.glob(pattern.name),
                        key=lambda candidate: str(candidate),
                    )
                )
        else:
            path_candidates = _fits_path_alternates(expanded)

        for candidate in path_candidates:
            resolved = candidate.resolve()
            if str(resolved) in seen:
                continue
            seen.add(str(resolved))
            candidates.append(resolved)

    return candidates


def find_first_existing(*paths: Path | None) -> Path | None:
    for candidate in _candidate_paths(*paths):
        if candidate.exists():
            return candidate
    return None


def resolve_existing_path(label: str, *paths: Path | None) -> Path:
    resolved = find_first_existing(*paths)
    if resolved is not None:
        return resolved

    checked = "\n".join(f"  - {candidate}" for candidate in _unique_paths(*paths))
    raise FileNotFoundError(f"Could not find {label}. Checked:\n{checked}")


def build_input_candidates(
    root: Path | None, relative_path: Path, flat_name: str | None = None
) -> list[Path]:
    if root is None:
        return []

    candidates = [root / relative_path]
    if flat_name is not None:
        candidates.append(root / flat_name)

    return candidates


def build_named_input_candidates(
    root: Path | None, relative_dir: Path, names: list[str]
) -> list[Path]:
    candidates: list[Path] = []
    for name in names:
        candidates.extend(build_input_candidates(root, relative_dir / name, name))
    return candidates


def _input_search_dirs(root: Path | None, relative_dir: Path) -> list[Path]:
    if root is None:
        return []
    return _unique_paths(root / relative_dir, root)


def _keyword_sort_key(keyword: str) -> tuple[bool, str, str]:
    return (keyword != "", keyword.lower(), keyword)


def _format_keyword(keyword: str) -> str:
    return "<none>" if keyword == "" else keyword


def extract_shared_keyword(filename: str, galaxy_name: str, suffixes: list[str]) -> str | None:
    prefix = f"{galaxy_name}_"
    if not filename.startswith(prefix):
        return None

    remainder = filename[len(prefix):]
    for suffix in sorted(suffixes, key=len, reverse=True):
        suffix_variants = [suffix]
        lower_suffix = suffix.lower()
        if lower_suffix.endswith(".fits"):
            suffix_variants.append(f"{suffix}.gz")
        elif lower_suffix.endswith(".fits.gz"):
            suffix_variants.append(suffix[:-3])

        for suffix_variant in suffix_variants:
            if remainder == suffix_variant:
                return ""
            marker = f"_{suffix_variant}"
            if remainder.endswith(marker):
                return remainder[: -len(marker)]

    return None


def collect_keyworded_input_paths(
    root: Path | None,
    fallback_root: Path | None,
    relative_dir: Path,
    galaxy_name: str,
    suffixes: list[str],
) -> dict[str, Path]:
    matches: dict[str, Path] = {}

    for search_root in (root, fallback_root):
        for search_dir in _input_search_dirs(search_root, relative_dir):
            for suffix in suffixes:
                for pattern in (f"{galaxy_name}_{suffix}", f"{galaxy_name}_*_{suffix}"):
                    for candidate in _candidate_paths(search_dir / pattern):
                        if not candidate.exists() or candidate.is_dir():
                            continue

                        keyword = extract_shared_keyword(
                            candidate.name, galaxy_name, suffixes
                        )
                        if keyword is None or keyword in matches:
                            continue
                        matches[keyword] = candidate.resolve()

    return matches


def find_keyworded_input_path(
    root: Path | None,
    fallback_root: Path | None,
    relative_dir: Path,
    galaxy_name: str,
    suffixes: list[str],
    keyword: str,
) -> Path | None:
    return collect_keyworded_input_paths(
        root, fallback_root, relative_dir, galaxy_name, suffixes
    ).get(keyword)


def resolve_keyworded_input_group(
    group_label: str,
    root: Path | None,
    fallback_root: Path | None,
    relative_dir: Path,
    galaxy_name: str,
    specs: dict[str, list[str]],
) -> tuple[str, dict[str, Path]]:
    matches_by_label = {
        label: collect_keyworded_input_paths(
            root, fallback_root, relative_dir, galaxy_name, suffixes
        )
        for label, suffixes in specs.items()
    }

    keyword_sets = [set(matches) for matches in matches_by_label.values()]
    common_keywords = set.intersection(*keyword_sets) if keyword_sets else set()
    if not common_keywords:
        found = []
        for label, matches in matches_by_label.items():
            if matches:
                keywords = ", ".join(
                    _format_keyword(keyword)
                    for keyword in sorted(matches, key=_keyword_sort_key)
                )
            else:
                keywords = "none"
            found.append(f"  - {label}: {keywords}")

        raise FileNotFoundError(
            f"Could not find {group_label} with a shared optional keyword.\n"
            f"Each required product input must use the same text between "
            f"'{galaxy_name}_' and its product suffix.\n"
            "Found keyword groups:\n" + "\n".join(found)
        )

    chosen_keyword = sorted(common_keywords, key=_keyword_sort_key)[0]
    return chosen_keyword, {
        label: matches_by_label[label][chosen_keyword] for label in specs
    }


def build_further_output_path(input_path: Path) -> Path:
    lower_name = input_path.name.lower()
    if lower_name.endswith(".fits.gz"):
        suffix = ".fits"
        stem = input_path.name[: -len(".fits.gz")]
    else:
        suffix = input_path.suffix or ".fits"
        stem = input_path.name[: -len(input_path.suffix)] if input_path.suffix else input_path.name
    if stem.endswith("_further"):
        return input_path.parent / f"{stem}{suffix}"
    if stem.endswith("_extended"):
        stem = stem[: -len("_extended")]
    return input_path.parent / f"{stem}_further{suffix}"


KNOWN_LINE_WAVELENGTH_UM = {
    "HB4861": 0.4861,
    "OIII4958": 0.4958,
    "OIII5006": 0.5006,
    "NII5754": 0.5755,
    "OI6300": 0.6300,
    "SIII6312": 0.6312,
    "OI6363": 0.6363,
    "NII6548": 0.6548,
    "HA6562": 0.6562,
    "NII6583": 0.6583,
    "SII6716": 0.6716,
    "SII6730": 0.6730,
    "OII7318": 0.7319,
    "OII7319": 0.7320,
    "OII7329": 0.7329,
    "OII7330": 0.7331,
    "SIII9068": 0.9069,
}

REQUIRED_GAS_LINE_BASES = (
    "HB4861",
    "HA6562",
    "OIII5006",
    "NII6583",
    "SII6716",
    "SII6730",
)

TE_LINE_WAVELENGTH_UM = {
    "OIII4958": 0.4958,
    "NII5754": 0.5755,
    "SIII6312": 0.6312,
    "NII6548": 0.6548,
    "OII7318": 0.7319,
    "OII7319": 0.7320,
    "OII7329": 0.7329,
    "OII7330": 0.7331,
    "SIII9068": 0.9069,
}

TE_LINE_BASES = tuple(TE_LINE_WAVELENGTH_UM)

SPATIAL_BASE_MAP_NAMES = (
    "V_STARS2",
    "SIGMA_STARS2",
)


def infer_line_wavelength_um(line_base: str) -> float:
    """Return the vacuum/air wavelength in microns from the known registry or trailing digits."""
    if line_base in KNOWN_LINE_WAVELENGTH_UM:
        return KNOWN_LINE_WAVELENGTH_UM[line_base]

    match = re.search(r"(\d{4})$", line_base)
    if match is None:
        raise ValueError(
            f"Cannot infer wavelength for emission-line base '{line_base}'. "
            "Add it to KNOWN_LINE_WAVELENGTH_UM before dust correction."
        )
    return float(match.group(1)) / 10000.0


def optional_image_or_nan(
    hdul: fits.HDUList, extension_name: str, shape: tuple[int, int]
) -> tuple[np.ndarray, bool]:
    if extension_name not in hdul:
        return np.full(shape, np.nan, dtype=float), False
    data = np.asarray(hdul[extension_name].data, dtype=float)
    if data.shape != shape:
        raise ValueError(
            f"{extension_name} has shape {data.shape}, expected {shape}"
        )
    return data, True


def load_spatial_keep_mask(mask_path: Path, target_shape: tuple[int, int]) -> tuple[np.ndarray, str]:
    """
    Load the nGIST spatial mask as a boolean keep mask.

    nGIST `MASKFILE[MASK]` stores one row per spatial pixel in row-major order.
    In the current MAUVE v3tk products, MASK == 0 is the finite/unmasked region
    and MASK == 1 is the excluded region.
    """
    with fits.open(mask_path) as hdul:
        for hdu in hdul[1:]:
            data = hdu.data
            if isinstance(data, np.ndarray) and data.ndim == 2 and data.shape == target_shape:
                return np.asarray(data) == 0, f"{hdu.name} image; MASK == 0"

        if "MASKFILE" not in hdul:
            raise ValueError(f"{mask_path} has no MASKFILE table or 2D mask image")

        table = hdul["MASKFILE"].data
        if table is None or table.dtype.names is None or "MASK" not in table.dtype.names:
            raise ValueError(f"{mask_path} MASKFILE table has no MASK column")

        mask_values = np.asarray(table["MASK"])
        expected_size = int(np.prod(target_shape))
        if mask_values.size != expected_size:
            raise ValueError(
                f"{mask_path} MASK column has {mask_values.size} rows, "
                f"but the target map shape {target_shape} requires {expected_size}"
            )

        mask_image = mask_values.reshape(target_shape, order="C")
        return mask_image == 0, "MASKFILE[MASK] table; row-major MASK == 0"


def apply_spatial_keep_mask(
    keep_mask: np.ndarray, maps: dict[str, np.ndarray]
) -> dict[str, np.ndarray]:
    masked_maps: dict[str, np.ndarray] = {}
    for name, data in maps.items():
        if data.shape != keep_mask.shape:
            raise ValueError(
                f"Cannot apply spatial mask to {name}: map shape {data.shape} "
                f"does not match mask shape {keep_mask.shape}"
            )
        masked_maps[name] = np.where(keep_mask, data, np.nan)
    return masked_maps


def read_spatial_binid_map(bin_path: Path, expected_shape: tuple[int, ...]) -> np.ndarray:
    with fits.open(bin_path) as hdul:
        if "BINID" in hdul:
            raw_binid = np.asarray(hdul["BINID"].data, dtype=np.float64)
        elif len(hdul) > 1:
            raw_binid = np.asarray(hdul[1].data, dtype=np.float64)
        else:
            raise KeyError(f"Could not find BINID extension in {bin_path.name}.")

    if raw_binid.shape != expected_shape:
        raise ValueError(
            f"BINID shape mismatch for {bin_path.name}: "
            f"{raw_binid.shape} vs gas-map shape {expected_shape}."
        )

    finite_bin = np.isfinite(raw_binid) & (raw_binid >= 0)
    if not np.any(finite_bin):
        raise ValueError(f"BINID map in {bin_path.name} has no finite non-negative bins.")

    rounded_binid = np.rint(raw_binid[finite_bin])
    if not np.allclose(raw_binid[finite_bin], rounded_binid, rtol=0.0, atol=1.0e-6):
        raise ValueError(f"BINID map in {bin_path.name} contains non-integer bin values.")

    binid = np.full(raw_binid.shape, -1, dtype=np.int32)
    binid[finite_bin] = rounded_binid.astype(np.int32)
    return binid


def collapse_map_to_spatial_bins(
    data: np.ndarray, binid_map: np.ndarray, label: str
) -> np.ndarray:
    data = np.asarray(data, dtype=np.float64)
    if data.shape != binid_map.shape:
        raise ValueError(
            f"{label} shape mismatch: map shape {data.shape} vs BINID shape {binid_map.shape}."
        )

    bin_pixels = binid_map >= 0
    if not np.any(bin_pixels):
        raise ValueError(f"Cannot bin-collapse {label}: BINID map has no valid pixels.")

    valid = bin_pixels & np.isfinite(data)
    max_binid = int(np.max(binid_map[bin_pixels]))
    sums = np.bincount(
        binid_map[valid].ravel(),
        weights=data[valid].ravel(),
        minlength=max_binid + 1,
    )
    counts = np.bincount(binid_map[valid].ravel(), minlength=max_binid + 1)
    means = np.divide(
        sums,
        counts,
        out=np.full(max_binid + 1, np.nan, dtype=np.float64),
        where=counts > 0,
    )

    binned = np.full(data.shape, np.nan, dtype=np.float64)
    binned[valid] = means[binid_map[valid]]
    return binned


def collapse_named_maps_to_spatial_bins(
    binid_map: np.ndarray, maps: dict[str, np.ndarray]
) -> dict[str, np.ndarray]:
    return {
        name: collapse_map_to_spatial_bins(data, binid_map, name)
        for name, data in maps.items()
    }


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
p.add_argument(
    "--root",
    default=".",
    help="Root directory containing v3tk_v7.6.8/{galaxy} products and prior post-processing outputs",
)
p.add_argument(
    "--product-subdir",
    default="v3tk_v7.6.8",
    help="Product subdirectory under --root, e.g. v3tk_v7.6.8 or v3tk_v7.6.8_7000",
)
p.add_argument(
    "--fallback-root",
    default=".",
    help="Fallback directory searched when an input file is not found under --root",
)
p.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
args = p.parse_args()

loglvl = logging.INFO if args.verbose else logging.WARNING
logging.basicConfig(level=loglvl, format="%(levelname)s %(message)s", stream=sys.stdout)

t0   = time.perf_counter()
gal  = args.galaxy.upper()
root = Path(args.root).expanduser().resolve()
product_subdir = Path(args.product_subdir)
fallback_root = (
    Path(args.fallback_root).expanduser().resolve()
    if args.fallback_root is not None
    else None
)

gas_suffixes = [
    "gas_BIN_maps.fits",
    "gas_bin_maps.fits",
    "BIN_maps.fits",
    "bin_maps.fits",
]
bin_further_suffixes = [
    "spatial_binning_maps_further.fits",
    "SPATIAL_BINNING_maps_further.fits",
]
mask_suffixes = [
    "mask.fits",
    "MASK.fits",
]

product_dir = product_subdir / gal
input_keyword, gas_paths = resolve_keyworded_input_group(
    "SFR gas FITS input",
    root,
    fallback_root,
    product_dir,
    gal,
    {
        "gas-line FITS": gas_suffixes,
    },
)
src = gas_paths["gas-line FITS"]
out_path = root / product_dir / build_further_output_path(src).name
bin_further_path = (
    find_keyworded_input_path(
        root, fallback_root, product_dir, gal, bin_further_suffixes, input_keyword
    )
    or find_keyworded_input_path(
        root, fallback_root, Path("."), gal, bin_further_suffixes, input_keyword
    )
)
if bin_further_path is None:
    raise FileNotFoundError(
        "Could not find further spatial-binning FITS for SFR+Z.py. "
        "Expected the Mass.py output in the working directory, e.g. "
        f"{gal}_spatial_binning_maps_further.fits, or under {product_dir}."
    )
mask_path = (
    find_keyworded_input_path(
        root, fallback_root, product_dir, gal, mask_suffixes, input_keyword
    )
    or find_keyworded_input_path(
        root, fallback_root, Path("."), gal, mask_suffixes, input_keyword
    )
)
if mask_path is None and input_keyword != "":
    mask_path = (
        find_keyworded_input_path(root, fallback_root, product_dir, gal, mask_suffixes, "")
        or find_keyworded_input_path(root, fallback_root, Path("."), gal, mask_suffixes, "")
    )
inclination_path = Path("MAUVE_Inclination.dat")

# two key cut parameters
cut = 3 # FLUX/FLUX_ERR
noise = 20 # detection limit of FLUX, in the unit of 10^-20 erg/s

# ------------------------------------------------------------------
# 1.  Load raw gas-line maps and match the Mass.py further product
# ------------------------------------------------------------------

print(f"Shared input keyword ➜ {_format_keyword(input_keyword)}")
print(f"Reading gas-line FITS ➜ {src}")
print(f"Matched binning FITS ➜ {bin_further_path}")
if mask_path is not None:
    print(f"Spatial mask FITS ➜ {mask_path}")
else:
    print("Spatial mask FITS ➜ <not found; leaving gas maps unchanged>")
with fits.open(src) as hdul:
    V_STARS2 = np.asarray(hdul['V_STARS2'].data, dtype=float)
    SIGMA_STARS2 = np.asarray(hdul['SIGMA_STARS2'].data, dtype=float)
    gas_header = hdul['HA6562_FLUX'].header.copy()

    map_shape = np.asarray(hdul['HA6562_FLUX'].data).shape
    gas_line_bases: list[str] = []
    gas_line_maps: dict[str, np.ndarray] = {}
    available_hdu_names = {hdu.name for hdu in hdul}

    for hdu in hdul:
        ext_name = hdu.name
        if not ext_name.endswith("_FLUX"):
            continue

        line_base = ext_name[: -len("_FLUX")]
        flux_name = f"{line_base}_FLUX"
        err_name = f"{line_base}_FLUX_ERR"
        if err_name not in available_hdu_names:
            print(f"Warning: skipping {line_base}; {err_name} is missing")
            continue

        flux_data = np.asarray(hdul[flux_name].data, dtype=float)
        err_data = np.asarray(hdul[err_name].data, dtype=float)
        if flux_data.shape != map_shape or err_data.shape != map_shape:
            raise ValueError(
                f"{line_base} flux/error shape mismatch: "
                f"{flux_data.shape}/{err_data.shape}, expected {map_shape}"
            )

        gas_line_bases.append(line_base)
        gas_line_maps[flux_name] = flux_data
        gas_line_maps[err_name] = err_data

    GAS_LINE_BASES = tuple(gas_line_bases)
    missing_required = [
        line_base for line_base in REQUIRED_GAS_LINE_BASES
        if line_base not in GAS_LINE_BASES
    ]
    if missing_required:
        raise KeyError(
            "Missing required gas-line flux/error pairs for SFR+Z.py: "
            + ", ".join(missing_required)
        )

    globals().update(gas_line_maps)

    te_line_availability: dict[str, bool] = {}
    te_optional_maps: dict[str, np.ndarray] = {}
    for line_base in TE_LINE_BASES:
        te_line_availability[line_base] = line_base in GAS_LINE_BASES
        if te_line_availability[line_base]:
            continue

        te_optional_maps[f"{line_base}_FLUX"] = np.full(map_shape, np.nan, dtype=float)
        te_optional_maps[f"{line_base}_FLUX_ERR"] = np.full(map_shape, np.nan, dtype=float)

globals().update(te_optional_maps)

gas_binid_map = read_spatial_binid_map(bin_further_path, map_shape)
gas_bin_map_names = SPATIAL_BASE_MAP_NAMES + tuple(
    f"{line_base}_{suffix}"
    for line_base in GAS_LINE_BASES
    for suffix in ("FLUX", "FLUX_ERR")
)
globals().update(
    collapse_named_maps_to_spatial_bins(
        gas_binid_map,
        {name: globals()[name] for name in gas_bin_map_names},
    )
)
print(
    "Gas-line maps collapsed onto spatial BINID: "
    f"{np.unique(gas_binid_map[gas_binid_map >= 0]).size} bins"
)

print("Available gas-line flux/error pairs:")
for line_base in GAS_LINE_BASES:
    print(f"  - {line_base}")

if ENABLE_TE_METALLICITY_PRODUCTS:
    print("Te-method optional emission-line availability:")
    for line_base in TE_LINE_BASES:
        status = "available" if te_line_availability[line_base] else "missing"
        print(f"  - {line_base}: {status}")

gas_header

# ------------------------------------------------------------------
# 2.  Spatial-mask check
# ------------------------------------------------------------------

if mask_path is not None:
    spatial_keep_mask, spatial_mask_source = load_spatial_keep_mask(
        mask_path, HA6562_FLUX.shape
    )
    spatial_masked_map_names = SPATIAL_BASE_MAP_NAMES + tuple(
        f"{line_base}_{suffix}"
        for line_base in GAS_LINE_BASES
        for suffix in ("FLUX", "FLUX_ERR")
    )
    globals().update(
        apply_spatial_keep_mask(
            spatial_keep_mask,
            {name: globals()[name] for name in spatial_masked_map_names},
        )
    )
    keep_count = int(np.count_nonzero(spatial_keep_mask))
    masked_count = int(spatial_keep_mask.size - keep_count)
    finite_after = int(np.count_nonzero(np.isfinite(HA6562_FLUX)))
    print(
        "Spatial mask applied: "
        f"{spatial_mask_source}; keep={keep_count}, masked={masked_count}, "
        f"finite H-alpha after mask={finite_after}"
    )
else:
    print("Spatial mask not found; skipping additional spatial-mask application.")

# ------------------------------------------------------------------
# 3.  Roadmap
# ------------------------------------------------------------------

# Road map:
# 1. Calculate the Balmer Decrement (BD) from Hβ and Hα
# 2. Convert BD to gas E(B-V) using the selected extinction curve (default: Calzetti 2000)
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


def ccm89_k(w_um, Rv=3.1):
    """
    CCM89 (Cardelli, Clayton & Mathis 1989): k(λ)=A(λ)/E(B−V) with λ in microns.

    Uses CCM89 eqs. (2)–(5) exactly:
      IR:      0.3 <= x < 1.1
      Opt/NIR: 1.1 <= x < 3.3, y=x-1.82
      UV:      3.3 <= x <= 8.0 with Fa,Fb terms for x>5.9
      Far-UV:  8.0 < x <= 10.0
    where x = 1/λ (micron^-1).
    """
    w = np.asarray(w_um, dtype=float)
    x = 1.0 / w  # micron^-1

    a = np.full_like(x, np.nan, dtype=float)
    b = np.full_like(x, np.nan, dtype=float)

    # (2) Infrared: 0.3 <= x < 1.1
    ir = (x >= 0.3) & (x < 1.1)
    a[ir] = 0.574 * x[ir]**1.61
    b[ir] = -0.527 * x[ir]**1.61

    # (3) Optical/NIR: 1.1 <= x < 3.3, y = x - 1.82
    opt = (x >= 1.1) & (x < 3.3)
    y = x[opt] - 1.82
    a[opt] = (1.0
              + 0.17699*y
              - 0.50447*y**2
              - 0.02427*y**3
              + 0.72085*y**4
              + 0.01979*y**5
              - 0.77530*y**6
              + 0.32999*y**7)
    b[opt] = (1.41338*y
              + 2.28305*y**2
              + 1.07233*y**3
              - 5.38434*y**4
              - 0.62251*y**5
              + 5.30260*y**6
              - 2.09002*y**7)

    # (4) Ultraviolet: 3.3 <= x <= 8.0
    uv = (x >= 3.3) & (x <= 8.0)
    a[uv] = 1.752 - 0.316*x[uv] - 0.104/((x[uv] - 4.67)**2 + 0.341)
    b[uv] = -3.090 + 1.825*x[uv] + 1.206/((x[uv] - 4.62)**2 + 0.263)

    # Fa,Fb curvature terms for 5.9 <= x <= 8.0
    fuv = (x >= 5.9) & (x <= 8.0)
    y = x[fuv] - 5.9
    a[fuv] += -0.04473*y**2 - 0.009779*y**3
    b[fuv] +=  0.2130*y**2 + 0.1207*y**3

    # (5) Far-UV: 8 < x <= 10, use (x-8)
    faruv = (x > 8.0) & (x <= 10.0)
    y = x[faruv] - 8.0
    a[faruv] = -1.073 - 0.628*y + 0.137*y**2 - 0.070*y**3
    b[faruv] = 13.670 + 4.257*y - 0.420*y**2 + 0.374*y**3

    # Convert to k(λ)=A(λ)/E(B−V) = Rv*a + b (from CCM89 eq. 1)
    k = Rv * a + b
    return k.item() if k.ndim == 1 and k.size == 1 else k


def extinction_k(w_um, law=None, Rv=None):
    """Return k(λ)=A(λ)/E(B−V) for the selected extinction law; wavelengths in microns."""
    if law is None:
        law = extinction_law
    law = str(law).lower()

    if law in ("mw", "milkyway", "ccm", "ccm89"):
        if Rv is None:
            Rv = mw_rv
        return ccm89_k(w_um, Rv=Rv)
    if law in ("calzetti", "calz", "c00"):
        return calzetti_k(w_um)

    raise ValueError(f"Unknown extinction law: {law}")

GAS_LINE_WAVELENGTH_UM = {
    line_base: infer_line_wavelength_um(line_base)
    for line_base in GAS_LINE_BASES
}
GAS_LINE_K = {
    line_base: extinction_k(w_um)
    for line_base, w_um in GAS_LINE_WAVELENGTH_UM.items()
}

k_HB4861  = GAS_LINE_K["HB4861"]
k_HA6562  = GAS_LINE_K["HA6562"]
k_OIII5006= GAS_LINE_K["OIII5006"]
k_NII6583 = GAS_LINE_K["NII6583"]
k_SII6716 = GAS_LINE_K["SII6716"]
k_SII6730 = GAS_LINE_K["SII6730"]
TE_LINE_K = {
    line_base: GAS_LINE_K.get(line_base, extinction_k(w_um))
    for line_base, w_um in TE_LINE_WAVELENGTH_UM.items()
}

# Print k(λ) values used for dust correction (this will also appear in redirected *.log outputs)
print("--------------------------------------------------------------")
if str(extinction_law).lower() in ("mw", "milkyway", "ccm", "ccm89"):
    print(f"Extinction law for k(λ)=A(λ)/E(B−V): {extinction_law} (mw_rv={mw_rv})")
else:
    print(f"Extinction law for k(λ)=A(λ)/E(B−V): {extinction_law}")
print(f"k(Hβ 4861Å)  = {float(k_HB4861):.4f}   at λ=0.4861 µm")
print(f"k(Hα 6562Å)  = {float(k_HA6562):.4f}   at λ=0.6562 µm")
print(f"k([OIII]5006)= {float(k_OIII5006):.4f}   at λ=0.5006 µm")
print(f"k([NII]6583) = {float(k_NII6583):.4f}   at λ=0.6583 µm")
print(f"k([SII]6716) = {float(k_SII6716):.4f}   at λ=0.6716 µm")
print(f"k([SII]6730) = {float(k_SII6730):.4f}   at λ=0.6730 µm")
for line_base in GAS_LINE_BASES:
    if line_base in REQUIRED_GAS_LINE_BASES:
        continue
    print(
        f"k({line_base}) = {float(GAS_LINE_K[line_base]):.4f}   "
        f"at λ={GAS_LINE_WAVELENGTH_UM[line_base]:.4f} µm"
    )
if ENABLE_TE_METALLICITY_PRODUCTS:
    for line_base in TE_LINE_BASES:
        if line_base in GAS_LINE_K:
            continue
        print(
            f"k({line_base}) = {float(TE_LINE_K[line_base]):.4f}   "
            f"at λ={TE_LINE_WAVELENGTH_UM[line_base]:.4f} µm (line missing; NaN placeholder)"
        )
print("--------------------------------------------------------------")

R_int = 2.86

# Define a function to convert the BD to gas E(B-V) 
# using the formula E(B-V)_BD = 2.5/(k_HB4861 - k_HA6562) * np.log10(BD/R_int)
def convert_bd_to_ebv(BD, k_HB4861, k_HA6562, R_int=2.86):
    E_BV_BD = 2.5 / (k_HB4861 - k_HA6562) * np.log10(BD / R_int)
    return E_BV_BD


def convert_bd_to_ebv_error(ha_flux, hb_flux, ha_err, hb_err, k_HB4861, k_HA6562):
    """First-order uncertainty in gas E(B-V) from the Balmer decrement."""
    coeff = 2.5 / ((k_HB4861 - k_HA6562) * np.log(10.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        valid = (
            np.isfinite(ha_flux) & np.isfinite(hb_flux) &
            np.isfinite(ha_err) & np.isfinite(hb_err) &
            (ha_flux > 0) & (hb_flux > 0) &
            (ha_err >= 0) & (hb_err >= 0)
        )
        ebv_err = np.abs(coeff) * np.sqrt((ha_err / ha_flux)**2 + (hb_err / hb_flux)**2)
    return np.where(valid, ebv_err, np.nan)


# Calculate the gas E(B-V) from BD
E_BV_BD = convert_bd_to_ebv(BD, k_HB4861, k_HA6562, R_int)
E_BV_BD_ERR = convert_bd_to_ebv_error(
    HA6562_FLUX, HB4861_FLUX, HA6562_FLUX_ERR, HB4861_FLUX_ERR,
    k_HB4861, k_HA6562
)

# Use E(B-V)_BD to correct the fluxes
def correct_flux_with_ebv(flux, ebv, k):
    """Correct flux with gas E(B-V) and extinction coefficient k."""
    return flux * 10**(0.4 * k * ebv)


def correct_flux_error_with_ebv(flux, flux_err, ebv, k, ebv_err=None):
    """Propagate line-flux and Balmer-decrement uncertainty through dust correction."""
    scale = 10**(0.4 * k * ebv)
    scaled_flux_err = flux_err * scale
    if ebv_err is None:
        return scaled_flux_err

    ebv_term = np.abs(flux * scale) * np.log(10.0) * 0.4 * np.abs(k) * ebv_err
    with np.errstate(invalid="ignore"):
        return np.sqrt(scaled_flux_err**2 + ebv_term**2)


def ratio_range_mask(values, low=None, high=None):
    """Finite mask with optional inclusive lower/upper diagnostic-ratio bounds."""
    mask = np.isfinite(values)
    if low is not None:
        mask &= values >= low
    if high is not None:
        mask &= values <= high
    return mask


def mask_scalar_by_range(value, diagnostic, low=None, high=None):
    """Return a scalar value only when its diagnostic ratio is inside range."""
    return value if bool(np.asarray(ratio_range_mask(diagnostic, low, high)).all()) else np.nan


def _real_roots(roots, real_atol=1e-8):
    realish = roots[np.abs(roots.imag) <= real_atol].real
    if realish.size == 0 and roots.size:
        idx = np.argmin(np.abs(roots.imag))
        if np.abs(roots[idx].imag) <= 1e-6:
            realish = np.array([roots[idx].real])
    return realish


C20_X_RANGE = (-0.7, 0.3)


def select_c20_root(roots, oh_prior=None, x_range=C20_X_RANGE):
    """Choose a C20 polynomial root deterministically inside the adopted branch range."""
    realish = _real_roots(roots)
    if realish.size == 0:
        return np.nan

    low, high = x_range
    candidates = realish[(realish >= low) & (realish <= high)]
    if candidates.size == 0:
        return np.nan

    if oh_prior is not None and np.isfinite(oh_prior):
        target = float(oh_prior) - 8.69
    else:
        target = 0.0
    return candidates[np.argmin(np.abs(candidates - target))]


def c20_prior_value(oh_prior, iy=None, ix=None):
    """Return a scalar branch-selection prior from an array/scalar, or None."""
    if isinstance(oh_prior, np.ndarray):
        value = oh_prior[iy, ix]
    else:
        value = oh_prior
    return float(value) if value is not None and np.isfinite(value) else None


def apply_metallicity_range(values, errors=None, low=7.63, high=9.23):
    """Mask metallicities, and optional errors, outside the adopted valid range."""
    valid = np.isfinite(values) & (values >= low) & (values <= high)
    values_out = np.where(valid, values, np.nan)
    if errors is None:
        return values_out
    errors_out = np.where(valid & np.isfinite(errors), errors, np.nan)
    return values_out, errors_out


def integrated_flux_error(flux_err_map, mask=None, flux=None, ebv=None, k=None, ebv_err=None):
    """Quadrature-sum a line-flux error map and optionally dust-correct the result."""
    if mask is None:
        err = np.sqrt(np.nansum(flux_err_map**2))
    else:
        err = np.sqrt(np.nansum(np.where(mask, flux_err_map, np.nan)**2))
    if ebv is not None and k is not None:
        if flux is None:
            flux = np.nan
        err = correct_flux_error_with_ebv(flux, err, ebv, k, ebv_err)
    return err


def as_single_pixel(value):
    """Represent an integrated scalar as a 1x1 map for shared map solvers."""
    return np.array([[value]], dtype=float)


def single_pixel_value(value):
    """Extract the scalar from a 1x1 map-like result."""
    return float(np.asarray(value, dtype=float).ravel()[0])


C20_METHOD_NAMES = ("O3N2", "O3S2", "RS32", "R3", "N2", "S2")


def combine_c20_measurements(values, errors):
    """
    Combine C20 metallicities without letting the smallest fitting scatter
    silently select one calibration as the whole "combined" result.
    """
    values = np.asarray(values, dtype=float)
    errors = np.asarray(errors, dtype=float)
    valid = np.isfinite(values) & np.isfinite(errors) & (errors > 0)

    weights = np.where(valid, 1.0 / errors**2, 0.0)
    sum_weights = np.sum(weights, axis=0)
    has_value = sum_weights > 0

    with np.errstate(divide="ignore", invalid="ignore"):
        combined = np.sum(np.where(valid, values * weights, 0.0), axis=0) / sum_weights
        formal_error = np.sqrt(1.0 / sum_weights)

    combined = np.where(has_value, combined, np.nan)
    formal_error = np.where(has_value, formal_error, np.nan)

    residuals = np.where(valid, values - np.expand_dims(combined, axis=0), 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        method_scatter = np.sqrt(np.sum(weights * residuals**2, axis=0) / sum_weights)

    n_methods = np.sum(valid, axis=0).astype(int)
    method_scatter = np.where(has_value & (n_methods > 1), method_scatter, 0.0)
    combined_error = np.where(
        has_value, np.sqrt(formal_error**2 + method_scatter**2), np.nan
    )

    dominant_method = np.argmax(weights, axis=0).astype(int)
    dominant_method = np.where(has_value, dominant_method, -1)
    return combined, combined_error, dominant_method, n_methods


def combine_c20_scalar(methods):
    """Scalar wrapper for integrated C20 totals."""
    values = np.array([value for _, value, _ in methods], dtype=float)
    errors = np.array([error for _, _, error in methods], dtype=float)
    combined, combined_error, dominant_method, n_methods = combine_c20_measurements(
        values, errors
    )
    return (
        single_pixel_value(combined),
        single_pixel_value(combined_error),
        int(np.asarray(dominant_method).ravel()[0]),
        int(np.asarray(n_methods).ravel()[0]),
    )


def c20_method_label(method_index):
    """Human-readable C20 method label for diagnostics."""
    if 0 <= method_index < len(C20_METHOD_NAMES):
        return C20_METHOD_NAMES[method_index]
    return "none"

# Correct every available gas line with E(B-V)_BD.
for line_base in GAS_LINE_BASES:
    k_line = GAS_LINE_K[line_base]
    globals()[f"{line_base}_FLUX_corr"] = correct_flux_with_ebv(
        globals()[f"{line_base}_FLUX"], E_BV_BD, k_line
    )
    globals()[f"{line_base}_FLUX_ERR_corr"] = correct_flux_error_with_ebv(
        globals()[f"{line_base}_FLUX"],
        globals()[f"{line_base}_FLUX_ERR"],
        E_BV_BD,
        k_line,
        E_BV_BD_ERR,
    )

# Keep missing optional Te inputs as NaN corrected maps so method-specific
# detection gates can fail cleanly instead of becoming hard requirements.
for line_base in TE_LINE_BASES:
    if line_base in GAS_LINE_BASES:
        continue
    globals()[f"{line_base}_FLUX_corr"] = np.full(HA6562_FLUX.shape, np.nan, dtype=float)
    globals()[f"{line_base}_FLUX_ERR_corr"] = np.full(HA6562_FLUX.shape, np.nan, dtype=float)

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


def ratio_linear_error(numerator, denominator, numerator_err, denominator_err):
    """Propagate the 1-sigma error for a linear ratio numerator/denominator."""
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = numerator / denominator
        ratio_err = np.abs(ratio) * np.sqrt(
            (numerator_err / numerator)**2 + (denominator_err / denominator)**2
        )
    return ratio_err


def line_detection_mask(flux, flux_err, cut_value=cut, noise_value=noise):
    """Configured emission-line detection gate using raw line flux and error maps."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return (
            np.isfinite(flux) & np.isfinite(flux_err) &
            (flux_err > 0) & (flux / flux_err >= cut_value) &
            (flux >= noise_value)
        )


def quadrature_sum(*arrays):
    """Quadrature sum for independent uncertainty maps."""
    total = np.zeros_like(np.asarray(arrays[0], dtype=float), dtype=float)
    valid = np.zeros_like(total, dtype=bool)
    for array in arrays:
        arr = np.asarray(array, dtype=float)
        total += np.where(np.isfinite(arr), arr**2, 0.0)
        valid |= np.isfinite(arr)
    return np.where(valid, np.sqrt(total), np.nan)


def positive_ratio(numerator, denominator):
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = numerator / denominator
    return np.where(
        np.isfinite(ratio) & (numerator > 0) & (denominator > 0), ratio, np.nan
    )


def _values_for_mask(value, valid):
    if isinstance(value, np.ndarray):
        return value[valid]
    return value


def pyneb_temden_map(
    atom,
    int_ratio,
    valid_mask,
    *,
    tem=-1,
    den=-1,
    wave1=-1,
    wave2=-1,
    to_eval=None,
):
    """Evaluate PyNeb getTemDen on valid pixels, falling back to scalar calls if needed."""
    result = np.full_like(np.asarray(int_ratio, dtype=float), np.nan, dtype=float)
    valid = (
        np.asarray(valid_mask, dtype=bool) &
        np.isfinite(int_ratio) & (int_ratio > 0)
    )
    if isinstance(tem, np.ndarray):
        valid &= np.isfinite(tem) & (tem > 0)
    if isinstance(den, np.ndarray):
        valid &= np.isfinite(den) & (den > 0)
    if not np.any(valid):
        return result

    try:
        values = atom.getTemDen(
            int_ratio[valid],
            tem=_values_for_mask(tem, valid),
            den=_values_for_mask(den, valid),
            wave1=wave1,
            wave2=wave2,
            to_eval=to_eval,
        )
        result[valid] = np.asarray(values, dtype=float)
        return result
    except Exception:
        pass

    for iy, ix in zip(*np.where(valid)):
        try:
            tem_value = tem[iy, ix] if isinstance(tem, np.ndarray) else tem
            den_value = den[iy, ix] if isinstance(den, np.ndarray) else den
            result[iy, ix] = atom.getTemDen(
                int_ratio[iy, ix],
                tem=tem_value,
                den=den_value,
                wave1=wave1,
                wave2=wave2,
                to_eval=to_eval,
            )
        except Exception:
            result[iy, ix] = np.nan
    return result


def sii_density_lookup_paths(tem, min_density, max_density, n_grid, lookup_dir=None):
    lookup_dir = Path.cwd() if lookup_dir is None else Path(lookup_dir)
    if (
        int(round(float(tem))) in SII_DENSITY_LOOKUP_TEMPERATURES and
        abs(float(min_density) - 20.0) < 1e-9 and
        abs(float(max_density) - 1.0e5) < 1e-6 and
        int(n_grid) == 4096
    ):
        basename = SII_DENSITY_LOOKUP_BASENAME
    else:
        basename = (
            "pyneb_sii_6716_6731_density_lookup_"
            f"te{float(tem):.0f}_ne{float(min_density):.0f}_"
            f"{float(max_density):.0f}_n{int(n_grid)}"
        )
    return lookup_dir / f"{basename}.npz", lookup_dir / f"{basename}.png"


def select_sii_ratio_row(temperature_grid, ratio_grid_2d, tem):
    temperature_grid = np.asarray(temperature_grid, dtype=float)
    ratio_grid_2d = np.asarray(ratio_grid_2d, dtype=float)
    matches = np.where(np.isclose(temperature_grid, float(tem), rtol=0.0, atol=1e-9))[0]
    if matches.size == 0:
        raise ValueError(
            f"Requested Te={tem} K is not in the cached SII lookup temperatures: "
            f"{temperature_grid}"
        )
    return np.asarray(ratio_grid_2d[int(matches[0])], dtype=float)


def save_sii_density_lookup_png(density_grid, temperature_grid, ratio_grid_2d, png_path):
    """Write a line-plot view of the cached [S II] density-ratio grid."""
    try:
        import matplotlib
        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt

        density_grid = np.asarray(density_grid, dtype=float)
        temperature_grid = np.asarray(temperature_grid, dtype=float)
        ratio_grid_2d = np.asarray(ratio_grid_2d, dtype=float)
        colors = plt.get_cmap("plasma")(
            np.linspace(0.08, 0.92, len(temperature_grid))
        )
        fig, ax = plt.subplots(figsize=(8.0, 5.0), constrained_layout=True)
        for temp, color, ratio_row in zip(temperature_grid, colors, ratio_grid_2d):
            ax.plot(
                density_grid,
                ratio_row,
                color=color,
                lw=2.0,
                label=f"{temp:.0f} K",
            )
        ax.set_xscale("log")
        ax.set_xlabel(r"$n_e$ (cm$^{-3}$)")
        ax.set_ylabel(r"$I(6716) / I(6731)$")
        ax.set_title(r"PyNeb [S II] $\lambda6716/\lambda6731$ density lookup")
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(title=r"$T_e$", ncol=2, fontsize=8)
        fig.savefig(png_path, dpi=180)
        plt.close(fig)
    except Exception as exc:
        print(f"Warning: could not write SII density lookup PNG {png_path}: {exc}")


def load_or_create_sii_density_lookup(
    atom,
    *,
    tem,
    min_density,
    max_density,
    n_grid,
    lookup_dir=None,
):
    """Load the reusable PyNeb [S II] density lookup from pwd, or create it."""
    npz_path, png_path = sii_density_lookup_paths(
        tem, min_density, max_density, n_grid, lookup_dir=lookup_dir
    )
    requested_tem = int(round(float(tem)))
    cache_key = (
        requested_tem,
        float(min_density),
        float(max_density),
        int(n_grid),
        str(npz_path.resolve()),
    )
    if cache_key in SII_DENSITY_LOOKUP_CACHE:
        return SII_DENSITY_LOOKUP_CACHE[cache_key]

    if npz_path.exists():
        try:
            with np.load(npz_path) as cached:
                density_grid = np.asarray(cached["density_grid"], dtype=float)
                temperature_grid = np.asarray(cached["temperature_grid"], dtype=float)
                ratio_grid_2d = np.asarray(cached["ratio_grid_2d"], dtype=float)
                cached_min = float(cached["min_density"])
                cached_max = float(cached["max_density"])
                cached_n = int(cached["n_grid"])
            valid_cache = (
                ratio_grid_2d.shape == (temperature_grid.size, density_grid.size) and
                density_grid.size == int(n_grid) and
                np.all(np.isfinite(density_grid)) and
                np.all(np.isfinite(temperature_grid)) and
                np.all(np.isfinite(ratio_grid_2d)) and
                tuple(int(round(t)) for t in temperature_grid) == SII_DENSITY_LOOKUP_TEMPERATURES and
                requested_tem in tuple(int(round(t)) for t in temperature_grid) and
                abs(cached_min - float(min_density)) < 1e-9 and
                abs(cached_max - float(max_density)) < 1e-6 and
                cached_n == int(n_grid)
            )
            if valid_cache:
                if not png_path.exists():
                    save_sii_density_lookup_png(
                        density_grid, temperature_grid, ratio_grid_2d, png_path
                    )
                ratio_grid = select_sii_ratio_row(
                    temperature_grid, ratio_grid_2d, requested_tem
                )
                print(f"Loaded SII density lookup table: {npz_path}")
                SII_DENSITY_LOOKUP_CACHE[cache_key] = (density_grid, ratio_grid)
                return density_grid, ratio_grid
            print(f"Warning: ignoring incompatible SII density lookup table: {npz_path}")
        except Exception as exc:
            print(f"Warning: could not read SII density lookup table {npz_path}: {exc}")

    density_grid = np.geomspace(min_density, max_density, n_grid)
    temperature_grid = np.asarray(SII_DENSITY_LOOKUP_TEMPERATURES, dtype=float)
    ratio_rows = []
    for temp in temperature_grid:
        emiss_6716 = atom.getEmissivity(tem=temp, den=density_grid, wave=6716)
        emiss_6731 = atom.getEmissivity(tem=temp, den=density_grid, wave=6731)
        ratio_rows.append(np.asarray(emiss_6716 / emiss_6731, dtype=float))
    ratio_grid_2d = np.vstack(ratio_rows)
    ratio_grid = select_sii_ratio_row(temperature_grid, ratio_grid_2d, requested_tem)

    try:
        np.savez_compressed(
            npz_path,
            density_grid=density_grid,
            temperature_grid=temperature_grid,
            ratio_grid_2d=ratio_grid_2d,
            selected_tem=float(requested_tem),
            selected_ratio_grid=ratio_grid,
            min_density=float(min_density),
            max_density=float(max_density),
            n_grid=int(n_grid),
            wave1=6716,
            wave2=6731,
            ratio_name="I(6716)/I(6731)",
        )
        save_sii_density_lookup_png(
            density_grid, temperature_grid, ratio_grid_2d, png_path
        )
        print(f"Saved SII density lookup table: {npz_path}")
        print(f"Saved SII density lookup line plot: {png_path}")
    except Exception as exc:
        print(f"Warning: could not write SII density lookup cache in {Path.cwd()}: {exc}")

    SII_DENSITY_LOOKUP_CACHE[cache_key] = (density_grid, ratio_grid)
    return density_grid, ratio_grid


def sii_density_from_ratio_lookup(
    atom,
    ratio,
    valid_mask,
    *,
    tem=10000.0,
    min_density=20.0,
    max_density=1.0e5,
    n_grid=4096,
    exact_high_density=False,
    exact_high_density_threshold=1.0e4,
    return_exact_mask=False,
):
    """Invert [S II] 6716/6731 using a PyNeb emissivity lookup table."""
    ratio = np.asarray(ratio, dtype=float)
    result = np.full_like(ratio, np.nan, dtype=float)
    exact_mask = np.zeros_like(ratio, dtype=bool)
    valid = np.asarray(valid_mask, dtype=bool) & np.isfinite(ratio) & (ratio > 0)
    if not np.any(valid):
        return (result, exact_mask) if return_exact_mask else result

    density_grid, ratio_grid = load_or_create_sii_density_lookup(
        atom,
        tem=tem,
        min_density=min_density,
        max_density=max_density,
        n_grid=n_grid,
    )

    order = np.argsort(ratio_grid)
    sorted_ratio = ratio_grid[order]
    sorted_density = density_grid[order]
    low_ratio = sorted_ratio[0]
    high_ratio = sorted_ratio[-1]

    interp_density = np.interp(
        ratio[valid],
        sorted_ratio,
        sorted_density,
        left=np.nan,
        right=min_density,
    )
    interp_density = np.where(
        ratio[valid] >= high_ratio,
        min_density,
        np.where(ratio[valid] < low_ratio, np.nan, interp_density),
    )
    result[valid] = interp_density

    if exact_high_density:
        exact_mask = valid & (
            ~np.isfinite(result) |
            (np.isfinite(result) & (result >= exact_high_density_threshold))
        )
        for iy, ix in zip(*np.where(exact_mask)):
            try:
                exact_density = atom.getTemDen(
                    ratio[iy, ix], tem=tem, wave1=6716, wave2=6731
                )
                if np.isfinite(exact_density) and exact_density > 0:
                    result[iy, ix] = exact_density
            except Exception:
                pass

    return (result, exact_mask) if return_exact_mask else result


def pyneb_ion_abundance_map(atom, int_ratio, tem, den, valid_mask, *, to_eval):
    """Evaluate PyNeb getIonAbundance with Hbeta normalized to 100."""
    result = np.full_like(np.asarray(int_ratio, dtype=float), np.nan, dtype=float)
    valid = (
        np.asarray(valid_mask, dtype=bool) &
        np.isfinite(int_ratio) & (int_ratio > 0) &
        np.isfinite(tem) & (tem > 0) &
        np.isfinite(den) & (den > 0)
    )
    if not np.any(valid):
        return result

    try:
        values = atom.getIonAbundance(
            int_ratio=int_ratio[valid],
            tem=tem[valid],
            den=den[valid],
            to_eval=to_eval,
            Hbeta=100.0,
        )
        result[valid] = np.asarray(values, dtype=float)
        return result
    except Exception:
        pass

    for iy, ix in zip(*np.where(valid)):
        try:
            result[iy, ix] = atom.getIonAbundance(
                int_ratio=int_ratio[iy, ix],
                tem=tem[iy, ix],
                den=den[iy, ix],
                to_eval=to_eval,
                Hbeta=100.0,
            )
        except Exception:
            result[iy, ix] = np.nan
    return result


def pyneb_ion_abundance_error_map(
    atom,
    int_ratio,
    int_ratio_err,
    tem,
    tem_err,
    den,
    den_err,
    valid_mask,
    *,
    to_eval,
):
    central = pyneb_ion_abundance_map(atom, int_ratio, tem, den, valid_mask, to_eval=to_eval)
    valid = valid_mask & np.isfinite(central) & (central > 0)

    ratio_low = np.where(np.isfinite(int_ratio_err), np.maximum(int_ratio - int_ratio_err, 1e-30), np.nan)
    ratio_high = np.where(np.isfinite(int_ratio_err), int_ratio + int_ratio_err, np.nan)
    abund_low = pyneb_ion_abundance_map(atom, ratio_low, tem, den, valid, to_eval=to_eval)
    abund_high = pyneb_ion_abundance_map(atom, ratio_high, tem, den, valid, to_eval=to_eval)
    err_ratio = finite_difference_error(abund_low, abund_high, valid)

    tem_low = np.where(np.isfinite(tem_err), np.maximum(tem - tem_err, 1000.0), np.nan)
    tem_high = np.where(np.isfinite(tem_err), tem + tem_err, np.nan)
    abund_tem_low = pyneb_ion_abundance_map(atom, int_ratio, tem_low, den, valid, to_eval=to_eval)
    abund_tem_high = pyneb_ion_abundance_map(atom, int_ratio, tem_high, den, valid, to_eval=to_eval)
    err_tem = finite_difference_error(abund_tem_low, abund_tem_high, valid)

    den_low = np.where(np.isfinite(den_err), np.maximum(den - den_err, 20.0), np.nan)
    den_high = np.where(np.isfinite(den_err), den + den_err, np.nan)
    abund_den_low = pyneb_ion_abundance_map(atom, int_ratio, tem, den_low, valid, to_eval=to_eval)
    abund_den_high = pyneb_ion_abundance_map(atom, int_ratio, tem, den_high, valid, to_eval=to_eval)
    err_den = finite_difference_error(abund_den_low, abund_den_high, valid)

    total_err = quadrature_sum(err_ratio, err_tem, err_den)
    return central, np.where(valid & np.isfinite(total_err), total_err, np.nan)


def finite_difference_error(low_values, high_values, central_mask):
    error = 0.5 * np.abs(high_values - low_values)
    return np.where(central_mask & np.isfinite(error), error, np.nan)


def sanitize_temperature_map(temperature, valid_mask):
    return np.where(
        valid_mask & np.isfinite(temperature) &
        (temperature >= 1000.0) & (temperature <= 30000.0),
        temperature,
        np.nan,
    )


def brazzini_te_siii_from_nii(te_nii):
    """Brazzini+2024 Te(SIII)-Te(NII); input/output in K."""
    return 1.22 * te_nii - 2000.0


def brazzini_te_siii_from_nii_error(te_nii, te_nii_err):
    """Includes coefficient errors and 724 K intrinsic dispersion from Brazzini+2024."""
    return np.sqrt(
        (1.22 * te_nii_err)**2 +
        (0.01 * te_nii)**2 +
        (0.01 * 1.0e4)**2 +
        724.0**2
    )


def brazzini_te_oiii_from_siii(te_siii):
    """Brazzini+2024 Eq. Te(OIII)-Te(SIII); input/output in K."""
    return 0.80 * te_siii + 2000.0


def brazzini_te_oiii_from_siii_error(te_siii, te_siii_err):
    """Includes coefficient errors and 1270 K intrinsic dispersion from Brazzini+2024."""
    return np.sqrt(
        (0.80 * te_siii_err)**2 +
        (0.02 * te_siii)**2 +
        (0.02 * 1.0e4)**2 +
        1270.0**2
    )


def oxygen_abundance_log_map(o_plus, o_plus_err, o_double_plus, o_double_plus_err, valid_mask):
    total_oxygen = o_plus + o_double_plus
    valid = valid_mask & np.isfinite(total_oxygen) & (total_oxygen > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        abundance = 12.0 + np.log10(total_oxygen)
        abundance_err = (
            np.sqrt(o_plus_err**2 + o_double_plus_err**2) /
            (np.log(10.0) * total_oxygen)
        )
    return (
        np.where(valid, abundance, np.nan),
        np.where(valid & np.isfinite(abundance_err), abundance_err, np.nan),
    )


def kreckel_mendez_nii_oh_map(te_nii, te_nii_err, valid_mask):
    with np.errstate(invalid="ignore"):
        te_unit = te_nii / 1.0e4
        abundance = -1.19 * te_unit + 9.68
        abundance_err = np.sqrt(
            (1.19 * te_nii_err / 1.0e4)**2 +
            (0.14 * te_unit)**2 +
            0.15**2
        )
    valid = (
        valid_mask &
        np.isfinite(abundance) &
        np.isfinite(te_nii) &
        (te_nii >= 8000.0) &
        (te_nii <= 13000.0)
    )
    return (
        np.where(valid, abundance, np.nan),
        np.where(valid & np.isfinite(abundance_err), abundance_err, np.nan),
    )

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

# N2S2-N06 metallicity calculation function
def calculate_n2s2_n06_metallicity(nii6583_flux, ha6562_flux, sii6716_flux, sii6730_flux):
    """Calculate [O/H] using N2S2-N06 calibration:
    log(N2S2) = log([NII]λ6584 / ([SII]λ6716+λ6731)) = -0.25214 + 0.74100*x + 0.58181*x² + 0.17963*x³
    where x = 12+log(O/H) - 8.69 = log(Z/Z☉) 
    """
    # Use basic finite checks on emission lines
    good_mask = (np.isfinite(nii6583_flux) & np.isfinite(ha6562_flux) &
                np.isfinite(sii6716_flux) & np.isfinite(sii6730_flux) &
                (nii6583_flux > 0) & (ha6562_flux > 0) &
                (sii6716_flux > 0) & (sii6730_flux > 0))
    
    # Initialize output arrays
    oh_n2s2_n06 = np.full_like(nii6583_flux, np.nan)
    
    # Calculate N2S2 ratio
    sii_total = sii6716_flux + sii6730_flux  # [SII] λ6716+λ6731
    n2s2_ratio = np.log10(nii6583_flux / sii_total)
    
    # Apply N2S2-N06 calibration - solve cubic polynomial for x
    # log(N2S2) = -0.25214 + 0.74100*x + 0.58181*x² + 0.17963*x³
    # Rearrange to: 0.17963*x³ + 0.58181*x² + 0.74100*x + (-0.25214 - log(N2S2)) = 0
    c3 = 0.17963
    c2 = 0.58181
    c1 = 0.74100
    c0 = -0.25214
    
    if np.any(good_mask):
        valid_indices = np.where(good_mask)
        for i in range(len(valid_indices[0])):
            idx_y, idx_x = valid_indices[0][i], valid_indices[1][i]
            n2s2_val = n2s2_ratio[idx_y, idx_x]
            
            # Solve cubic equation: c3*x³ + c2*x² + c1*x + (c0 - n2s2_val) = 0
            poly_coeffs = [c3, c2, c1, (c0 - n2s2_val)]
            roots = np.roots(poly_coeffs)
            
            # Select the real root (use first real root found)
            real_roots = roots[np.isreal(roots)].real
            if len(real_roots) > 0:
                # Take the first real root without range restrictions
                x_final = real_roots[0]
                oh_n2s2_n06[idx_y, idx_x] = x_final + 8.69
    
    return oh_n2s2_n06, good_mask

# Calculate N2S2-N06 metallicity
O_H_N2S2_N06, n2s2_n06_good_mask = calculate_n2s2_n06_metallicity(
    NII6583_FLUX_corr, HA6562_FLUX_corr, SII6716_FLUX_corr, SII6730_FLUX_corr)
# N2S2-N06 metallicity calculated without range restrictions

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
    good_mask &= ratio_range_mask(o3n2_ratio, -1.1, 1.7)
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
    good_mask &= ratio_range_mask(n2_ratio, -1.6, -0.2)
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
    good_mask &= ratio_range_mask(o3n2_ratio, None, 1.9)
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
    good_mask &= ratio_range_mask(n2_ratio, -2.5, -0.3)
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
        x_solution1 = (-b + np.sqrt(discriminant[combined_mask])) / (2*a)
        x_solution2 = (-b - np.sqrt(discriminant[combined_mask])) / (2*a)

        idxs = np.argwhere(combined_mask)
        y_err_values = y_err[combined_mask]
        prior_values = (
            oh_d16_sf[combined_mask]
            if isinstance(oh_d16_sf, np.ndarray)
            else np.full_like(
                y_err_values,
                c20_prior_value(oh_d16_sf) if c20_prior_value(oh_d16_sf) is not None else np.nan,
                dtype=float,
            )
        )
        for idx, (iy, ix) in enumerate(idxs):
            roots = np.array([x_solution1[idx], x_solution2[idx]], dtype=complex)
            x_final = select_c20_root(roots, prior_values[idx])
            if not np.isfinite(x_final):
                continue

            derivative_x = np.abs(c1 + 2*c2*x_final)
            if derivative_x <= 0:
                continue
            x_err = y_err_values[idx] / derivative_x

            oh_o3n2_c20[iy, ix] = x_final + 8.69
            fitting_err = 0.09  # dex
            oh_o3n2_c20_err[iy, ix] = np.sqrt(x_err**2 + fitting_err**2)

    combined_mask = good_mask & np.isfinite(oh_o3n2_c20)

    return oh_o3n2_c20, oh_o3n2_c20_err, combined_mask

# Calculate O3N2-C20 metallicity with error propagation
O_H_O3N2_C20, O_H_O3N2_C20_ERR, o3n2_c20_good_mask = calculate_o3n2_c20_metallicity(
    HB4861_FLUX_corr, OIII5006_FLUX_corr, NII6583_FLUX_corr, HA6562_FLUX_corr,
    HB4861_FLUX_ERR_corr, OIII5006_FLUX_ERR_corr, NII6583_FLUX_ERR_corr, HA6562_FLUX_ERR_corr,
    O_H_D16)
# Set O_H_O3N2_C20 to be nan if outside the range of 7.63 and 9.23
O_H_O3N2_C20, O_H_O3N2_C20_ERR = apply_metallicity_range(
    O_H_O3N2_C20, O_H_O3N2_C20_ERR
)

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
            
            oh_prior = c20_prior_value(oh_d16_sf, idx_y, idx_x)
            x_final = select_c20_root(roots, oh_prior)
            if np.isfinite(x_final):
                derivative_x = (np.abs(c1 + 2*c2*x_final + 3*c3*x_final**2 + 4*c4*x_final**3))
                if derivative_x <= 0:
                    combined_mask[idx_y, idx_x] = False
                    continue
                x_err = y_err_val / derivative_x

                oh_o3s2_c20[idx_y, idx_x] = x_final + 8.69
                fitting_err = 0.11  # dex
                oh_o3s2_c20_err[idx_y, idx_x] = np.sqrt(x_err**2 + fitting_err**2)
            else:
                combined_mask[idx_y, idx_x] = False
    
    # Update the combined mask to only include spaxels where we found valid solutions
    combined_mask = combined_mask & np.isfinite(oh_o3s2_c20)
    
    return oh_o3s2_c20, oh_o3s2_c20_err, combined_mask

# Calculate O3S2-C20 metallicity with error propagation
O_H_O3S2_C20, O_H_O3S2_C20_ERR, o3s2_c20_good_mask = calculate_o3s2_c20_metallicity(
    HB4861_FLUX_corr, OIII5006_FLUX_corr, SII6716_FLUX_corr, SII6730_FLUX_corr,
    HB4861_FLUX_ERR_corr, OIII5006_FLUX_ERR_corr, SII6716_FLUX_ERR_corr, SII6730_FLUX_ERR_corr,
    O_H_D16)
# Set O_H_O3S2_C20 to be nan if outside the range of 7.63 and 9.23
O_H_O3S2_C20, O_H_O3S2_C20_ERR = apply_metallicity_range(
    O_H_O3S2_C20, O_H_O3S2_C20_ERR
)

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
        
        # Calculate linear-ratio errors before propagating through log10(A + B).
        oiii_hb_err = ratio_linear_error(
            oiii5006_flux[good_mask], hb4861_flux[good_mask],
            oiii5006_flux_err[good_mask], hb4861_flux_err[good_mask]
        )
        sii_ha_err = ratio_linear_error(
            sii_total, ha6563_flux[good_mask],
            sii_total_err, ha6563_flux_err[good_mask]
        )
        
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
            oh_prior = c20_prior_value(oh_d16_sf, iy, ix)
            x_final = select_c20_root(roots, oh_prior)
            if np.isfinite(x_final):
                derivative_x = (np.abs(c1 + 2*c2*x_final + 3*c3*x_final**2 + 4*c4*x_final**3))
                if derivative_x <= 0:
                    continue
                x_err = y_err_val / derivative_x

                oh_rs32_c20[iy, ix] = x_final + 8.69
                fitting_err = 0.08  # dex
                oh_rs32_c20_err[iy, ix] = np.sqrt(x_err**2 + fitting_err**2)

    combined_mask = good_mask & np.isfinite(oh_rs32_c20)
    return oh_rs32_c20, oh_rs32_c20_err, combined_mask

# Calculate RS32-C20 metallicity
O_H_RS32_C20, O_H_RS32_C20_ERR, rs32_c20_good_mask = calculate_rs32_c20_metallicity(HB4861_FLUX_corr, HA6562_FLUX_corr, 
                                                                  OIII5006_FLUX_corr, SII6716_FLUX_corr, SII6730_FLUX_corr,
                                                                  HB4861_FLUX_ERR_corr, HA6562_FLUX_ERR_corr,
                                                                  OIII5006_FLUX_ERR_corr, SII6716_FLUX_ERR_corr, SII6730_FLUX_ERR_corr,
                                                                  O_H_D16)
# Set O_H_RS32_C20 to be nan if outside the range of 7.63 and 9.23
O_H_RS32_C20, O_H_RS32_C20_ERR = apply_metallicity_range(
    O_H_RS32_C20, O_H_RS32_C20_ERR
)

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
            oh_prior = c20_prior_value(oh_d16_sf, iy, ix)
            x_final = select_c20_root(roots, oh_prior)
            if np.isfinite(x_final):
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
                    total_error = np.sqrt(obs_error**2 + fitting_error**2)
                    oh_r3_c20_err[iy, ix] = total_error

    combined_mask = good_mask & np.isfinite(oh_r3_c20)
    return oh_r3_c20, oh_r3_c20_err, combined_mask

# Calculate R3-C20 metallicity
O_H_R3_C20, O_H_R3_C20_ERR, r3_c20_good_mask = calculate_r3_c20_metallicity(HB4861_FLUX_corr, HB4861_FLUX_ERR_corr,
                                                                             OIII5006_FLUX_corr, OIII5006_FLUX_ERR_corr,
                                                                             O_H_D16)
# Set O_H_R3_C20 to be nan if outside the range of 7.63 and 9.23
O_H_R3_C20, O_H_R3_C20_ERR = apply_metallicity_range(
    O_H_R3_C20, O_H_R3_C20_ERR
)

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
      • If multiple such roots exist, choose the root closest to the D16 prior.
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

        for idx, ((iy, ix), y_val) in enumerate(zip(idxs, y)):
            if not np.isfinite(y_val):
                continue

            roots = np.roots([c4, c3, c2, c1, (c0 - y_val)])

            oh_prior = c20_prior_value(oh_d16_sf, iy, ix)
            x_final = select_c20_root(roots, oh_prior)
            if not np.isfinite(x_final):
                continue
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
                total_error = np.sqrt(obs_error**2 + fitting_error**2)
                oh_n2_c20_err[iy, ix] = total_error

    combined_mask = good_mask & np.isfinite(oh_n2_c20)
    return oh_n2_c20, oh_n2_c20_err, combined_mask

# Calculate N2-C20 metallicity
O_H_N2_C20, O_H_N2_C20_ERR, n2_c20_good_mask = calculate_n2_c20_metallicity(HA6562_FLUX_corr, HA6562_FLUX_ERR_corr,
                                                                             NII6583_FLUX_corr, NII6583_FLUX_ERR_corr,
                                                                             O_H_D16)
# Set O_H_N2_C20 to be nan if outside the range of 7.63 and 9.23
O_H_N2_C20, O_H_N2_C20_ERR = apply_metallicity_range(
    O_H_N2_C20, O_H_N2_C20_ERR
)

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
      • If multiple, choose the root closest to the D16 prior.
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

        for idx, ((iy, ix), y_val) in enumerate(zip(idxs, y)):
            if not np.isfinite(y_val):
                continue

            # Solve: c4*x^4 + c3*x^3 + c2*x^2 + c1*x + (c0 - y) = 0
            roots = np.roots([c4, c3, c2, c1, (c0 - y_val)])

            oh_prior = c20_prior_value(oh_d16_sf, iy, ix)
            x_final = select_c20_root(roots, oh_prior)
            if not np.isfinite(x_final):
                continue
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
                total_error = np.sqrt(obs_error**2 + fitting_error**2)
                oh_s2_c20_err[iy, ix] = total_error

    combined_mask = good_mask & np.isfinite(oh_s2_c20)
    return oh_s2_c20, oh_s2_c20_err, combined_mask

# Calculate S2-C20 metallicity
O_H_S2_C20, O_H_S2_C20_ERR, s2_c20_good_mask = calculate_s2_c20_metallicity(HA6562_FLUX_corr, HA6562_FLUX_ERR_corr,
                                                                             SII6716_FLUX_corr, SII6716_FLUX_ERR_corr,
                                                                             SII6730_FLUX_corr, SII6730_FLUX_ERR_corr,
                                                                             O_H_D16)
# Set O_H_S2_C20 to be nan if outside the range of 7.63 and 9.23
O_H_S2_C20, O_H_S2_C20_ERR = apply_metallicity_range(
    O_H_S2_C20, O_H_S2_C20_ERR
)

# Combined C20 metallicity calculation function

def calculate_combined_c20_metallicity(gal):
    """
    Calculate combined C20 metallicity from all finite C20 methods.
    
    For each spaxel, we:
    1. Calculate metallicity and error for all 6 C20 methods
    2. Compute an inverse-variance weighted mean
    3. Add method-to-method scatter to the formal combined error
    
    Returns:
        oh_combined_c20: Combined metallicity map
        oh_combined_c20_err: Combined error map
        method_map: Dominant-weight method for each spaxel (0-5)
        combined_mask: Combined valid spaxel mask
    """
    # Reuse the already loaded and corrected arrays from this run.
    hb4861_flux = HB4861_FLUX_corr
    oiii5006_flux = OIII5006_FLUX_corr
    sii6716_flux = SII6716_FLUX_corr
    sii6730_flux = SII6730_FLUX_corr

    hb4861_flux_err = HB4861_FLUX_ERR_corr
    oiii5006_flux_err = OIII5006_FLUX_ERR_corr
    sii6716_flux_err = SII6716_FLUX_ERR_corr
    sii6730_flux_err = SII6730_FLUX_ERR_corr

    ha6563_flux = HA6562_FLUX_corr
    ha6563_flux_err = HA6562_FLUX_ERR_corr
    nii6584_flux = NII6583_FLUX_corr
    nii6584_flux_err = NII6583_FLUX_ERR_corr
    oh_d16_sf = O_H_D16
    
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
    
    all_metallicities = np.where(all_masks, all_metallicities, np.nan)
    all_errors = np.where(all_masks, all_errors, np.nan)
    oh_combined_c20, oh_combined_c20_err, method_map, _ = combine_c20_measurements(
        all_metallicities, all_errors
    )
    combined_mask = np.isfinite(oh_combined_c20)
    
    return oh_combined_c20, oh_combined_c20_err, method_map, combined_mask

# Calculate combined C20 metallicity
O_H_COMBINED_C20, O_H_COMBINED_C20_ERR, combined_c20_method_map, combined_c20_good_mask = calculate_combined_c20_metallicity(gal)
# Set combined C20 to be nan if outside the range of 7.63 and 9.23
O_H_COMBINED_C20, O_H_COMBINED_C20_ERR = apply_metallicity_range(
    O_H_COMBINED_C20, O_H_COMBINED_C20_ERR
)
combined_c20_method_map = np.where(np.isfinite(O_H_COMBINED_C20), combined_c20_method_map, -1)

print(f"Combined C20 metallicity: median = {np.nanmedian(O_H_COMBINED_C20):.3f}, range = ({np.nanmin(O_H_COMBINED_C20):.3f}, {np.nanmax(O_H_COMBINED_C20):.3f})")
print(f"Combined C20 dominant-weight method usage: O3N2={np.sum(combined_c20_method_map==0)}, O3S2={np.sum(combined_c20_method_map==1)}, RS32={np.sum(combined_c20_method_map==2)}, R3={np.sum(combined_c20_method_map==3)}, N2={np.sum(combined_c20_method_map==4)}, S2={np.sum(combined_c20_method_map==5)}")

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
def flux_to_luminosity(flux, distance=DISTANCE_MPC):
    """
    Convert flux to luminosity.
    
    Parameters:
    flux : array-like
        Integrated line flux in 1e-20 erg/(s cm^2).
    distance : float
        Distance in Mpc.
        
    Returns:
    luminosity : array-like
        Luminosity in erg/s.
    """
    return (flux*1e-20*u.erg/u.s/u.cm**2 * 4*np.pi*(distance*u.Mpc)**2).cgs.value

# Calculate the luminosity of Halpha
HA6562_LUM = flux_to_luminosity(HA6562_FLUX_Corr)

# SFR map from Halpha luminosity, using Kennicutt & Evans (2012)
# Kroupa-to-Chabrier IMF conversion encoded in SFR_HA_CHABRIER_COEFF.
def calzetti_sfr(luminosity):
    """
    Convert Halpha luminosity to SFR with a Chabrier IMF coefficient.
    
    Parameters:
    luminosity : array-like
        Halpha luminosity in erg/s.
        
    Returns:
    sfr : array-like
        Star formation rate in solar masses per year.
    """
    return SFR_HA_CHABRIER_COEFF * luminosity  # 4.98e-42 for Chabrier IMF

# Calculate the SFR map from Halpha luminosity
SFR_map = calzetti_sfr(HA6562_LUM)

# Getting the SFR surface density
# Convert to surface density in M☉/yr/kpc²
# 1. Convert pixel area to physical area in kpc²
legacy_wcs2 = WCS(gas_header).celestial  # strip spectral axis
pixel_scale = (proj_plane_pixel_scales(legacy_wcs2) * u.deg).to(u.arcsec)
pixel_area_Mpc = ((pixel_scale[0]).to(u.rad).value*DISTANCE_MPC*u.Mpc)*(((pixel_scale[1]).to(u.rad).value*DISTANCE_MPC*u.Mpc))
pixel_area_kpc = pixel_area_Mpc.to(u.kpc**2)

# 2. Read galaxy inclination and calculate correction factor
if apply_inclination_correction:
    galaxy_inclination = read_galaxy_inclination(gal, inclination_path)
    if galaxy_inclination is not None:
        inclination_rad = np.deg2rad(galaxy_inclination)
        cos_inclination = np.cos(inclination_rad)
        # Calculate b/a factor: sqrt((1-q0^2)*cos^2(i) + q0^2) where q0=0.2
        ba_factor = np.abs(np.sqrt((1-0.2**2)*cos_inclination**2 + 0.2**2))
        print(f"Galaxy {gal} inclination: {galaxy_inclination}° (cos θ = {cos_inclination:.3f})")
        print(f"Inclination correction ENABLED: applying b/a = {ba_factor:.3f} (adopting intrinsic thickness q₀ = 0.2 for disc galaxy)")
    else:
        ba_factor = 1.0
        print(f"No inclination data found for {gal}, using ba_factor = 1.0")
else:
    ba_factor = 1.0
    print(f"Inclination correction DISABLED: using ba_factor = 1.0")

# 3. Convert SFR to surface density with inclination correction
SFR_surface_density_map = SFR_map / pixel_area_kpc.value
SFR_surface_density_map_corrected = SFR_surface_density_map * ba_factor  # Apply inclination correction
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
# Current classification uses the measured non-Balmer fluxes even when those
# lines fail QC; those cases are tracked as low-S/N/unclassified rather than
# handled with formal upper/lower-limit BPT censoring.
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
logN2_err = bpt_error_propagation(NII6583_FLUX_corr, HA6562_FLUX_corr,
                                   NII6583_FLUX_ERR_corr, HA6562_FLUX_ERR_corr)
logS2_err = bpt_error_propagation(SII6716_FLUX_corr + SII6730_FLUX_corr, HA6562_FLUX_corr,
                                    np.sqrt(SII6716_FLUX_ERR_corr**2 + SII6730_FLUX_ERR_corr**2), HA6562_FLUX_ERR_corr)
logO3_err = bpt_error_propagation(OIII5006_FLUX_corr, HB4861_FLUX_corr,
                                   OIII5006_FLUX_ERR_corr, HB4861_FLUX_ERR_corr)

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

# ====== Classification 1: SF = HII + Comp ======
mask_classified1_N2_HII_Comp = ((mask_N2_left_HII | mask_N2_left_Comp) & 
                               (mask_N2_right_HII | mask_N2_right_Comp) & 
                                 (mask_N2_down_HII | mask_N2_down_Comp) &
                                 (mask_N2_up_HII | mask_N2_up_Comp))
mask_classified1_N2_AGN = (mask_N2_left_AGN & mask_N2_right_AGN &
                          mask_N2_down_AGN & mask_N2_up_AGN)
mask_classified1_S2_HII = (mask_S2_left_HII & mask_S2_right_HII &
                          mask_S2_down_HII & mask_S2_up_HII)
mask_classified1_S2_LINER = (mask_S2_left_LINER & mask_S2_right_LINER &
                            mask_S2_down_LINER & mask_S2_up_LINER)
mask_classified1_S2_Seyfert = (mask_S2_left_Seyfert & mask_S2_right_Seyfert &
                             mask_S2_down_Seyfert & mask_S2_up_Seyfert)

# NII classified1 and unclassified1 masks
mask_N2_classified1 = (mask_classified1_N2_HII_Comp | mask_classified1_N2_AGN)
mask_N2_unclassified1 = ~mask_N2_classified1
# SII classified1 and unclassified1 masks
mask_S2_classified1 = (mask_classified1_S2_HII | mask_classified1_S2_LINER | mask_classified1_S2_Seyfert)
mask_S2_unclassified1 = ~mask_S2_classified1

# ====== Classification 2: SF = HII only ======
# For Classification 2, we need separate HII and Comp+AGN masks from Classification 1
mask_classified2_N2_HII = mask_classified1_N2_HII_Comp & mask_N2_HII  # Only HII part from HII_Comp
mask_classified2_N2_Comp_AGN = (mask_classified1_N2_HII_Comp & mask_N2_Comp) | mask_classified1_N2_AGN  # Comp + AGN

# NII classified2 and unclassified2 masks
mask_classified2_N2 = (mask_classified2_N2_HII | mask_classified2_N2_Comp_AGN)
mask_unclassified2_N2 = ~mask_classified2_N2  # Same as unclassified1 since detection criteria unchanged
mask_classified2_N2_SF = mask_classified2_N2_HII  # Only HII is SF in Classification 2
mask_classified2_N2_nonSF = mask_classified2_N2_Comp_AGN  # Comp + AGN are non-SF

# For S2 BPT, Classification 2 is same as Classification 1 since S2 only has HII as SF anyway
mask_classified2_S2_HII = mask_classified1_S2_HII
mask_classified2_S2_LINER_Seyfert = (mask_classified1_S2_LINER | mask_classified1_S2_Seyfert)

# SII classified2 and unclassified2 masks  
mask_classified2_S2 = (mask_classified2_S2_HII | mask_classified2_S2_LINER_Seyfert)
mask_unclassified2_S2 = ~mask_classified2_S2  # Same as unclassified1
mask_classified2_S2_SF = mask_classified2_S2_HII  # Only HII is SF
mask_classified2_S2_nonSF = mask_classified2_S2_LINER_Seyfert  # LINER + Seyfert are non-SF

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


# ------------------------------------------------------------------
# 7.  Final Masks: Track 4 cases for each classification
# ------------------------------------------------------------------

# ====== Classification 1: SF = HII + Comp (current default) ======
# definite SF spaxels: or HA_detected_HB_detected & mask_N2_classified1 & mask_N2_SF
mask_SF_N2 = ((HA_detected_HB_detected_NII_detected_OIII_detected & mask_N2_classified1 & mask_N2_SF) | 
                    (HA_detected_HB_detected_NII_not_detected_OIII_not_detected & mask_N2_classified1 & mask_N2_SF) | 
                    (HA_detected_HB_detected_NII_not_detected_OIII_detected & mask_N2_classified1 & mask_N2_SF) | 
                    (HA_detected_HB_detected_NII_detected_OIII_not_detected & mask_N2_classified1 & mask_N2_SF))
# get SFR as non-SF: : or HA_detected_HB_detected & mask_N2_classified1 & mask_N2_nonSF
mask_nonSF_N2 = ((HA_detected_HB_detected_NII_detected_OIII_detected & mask_N2_classified1 & mask_N2_nonSF) |
              (HA_detected_HB_detected_NII_not_detected_OIII_not_detected & mask_N2_classified1 & mask_N2_nonSF) |
              (HA_detected_HB_detected_NII_not_detected_OIII_detected & mask_N2_classified1 & mask_N2_nonSF) |
              (HA_detected_HB_detected_NII_detected_OIII_not_detected & mask_N2_classified1 & mask_N2_nonSF))
# all the unclassified1 spaxels: : or HA_detected_HB_detected & mask_N2_unclassified1
mask_unclassified1_N2_final = ((HA_detected_HB_detected_NII_not_detected_OIII_not_detected & mask_N2_unclassified1) | 
                       (HA_detected_HB_detected_NII_not_detected_OIII_detected & mask_N2_unclassified1) | 
                       (HA_detected_HB_detected_NII_detected_OIII_not_detected & mask_N2_unclassified1) | 
                       (HA_detected_HB_detected_NII_detected_OIII_detected & mask_N2_unclassified1))

# ====== Classification 2: SF = HII only (more conservative) ======
# HII spaxels only in N2 BPT
mask_SF_N2_class2 = ((HA_detected_HB_detected_NII_detected_OIII_detected & mask_classified2_N2 & mask_classified2_N2_SF) | 
                    (HA_detected_HB_detected_NII_not_detected_OIII_not_detected & mask_classified2_N2 & mask_classified2_N2_SF) | 
                    (HA_detected_HB_detected_NII_not_detected_OIII_detected & mask_classified2_N2 & mask_classified2_N2_SF) | 
                    (HA_detected_HB_detected_NII_detected_OIII_not_detected & mask_classified2_N2 & mask_classified2_N2_SF))
# Comp + AGN as non-SF in classification 2
mask_nonSF_N2_class2 = ((HA_detected_HB_detected_NII_detected_OIII_detected & mask_classified2_N2 & mask_classified2_N2_nonSF) |
                  (HA_detected_HB_detected_NII_not_detected_OIII_not_detected & mask_classified2_N2 & mask_classified2_N2_nonSF) |
                  (HA_detected_HB_detected_NII_not_detected_OIII_detected & mask_classified2_N2 & mask_classified2_N2_nonSF) |
                  (HA_detected_HB_detected_NII_detected_OIII_not_detected & mask_classified2_N2 & mask_classified2_N2_nonSF))
# all the unclassified2 spaxels (same as unclassified1)
mask_unclassified2_N2_final = ((HA_detected_HB_detected_NII_not_detected_OIII_not_detected & mask_unclassified2_N2) | 
                       (HA_detected_HB_detected_NII_not_detected_OIII_detected & mask_unclassified2_N2) | 
                       (HA_detected_HB_detected_NII_detected_OIII_not_detected & mask_unclassified2_N2) | 
                       (HA_detected_HB_detected_NII_detected_OIII_detected & mask_unclassified2_N2))
# the rest are upper-limit spaxels: 
mask_upper = (HA_not_detected | HA_detected_HB_not_detected)

# Something else might be useful

# all the classified1 spaxels
mask_classified1_N2 = mask_SF_N2 | mask_nonSF_N2

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


# ====== Classification 1: SF = HII (S2 BPT only has HII as SF anyway) ======
# definite SF spaxels: or HA_detected_HB_detected & mask_S2_classified1 & mask_S2_SF
mask_SF_S2 = ((HA_detected_HB_detected_SII_detected_OIII_detected & mask_S2_classified1 & mask_S2_SF) | 
              (HA_detected_HB_detected_SII_not_detected_OIII_not_detected & mask_S2_classified1 & mask_S2_SF) | 
              (HA_detected_HB_detected_SII_not_detected_OIII_detected & mask_S2_classified1 & mask_S2_SF) | 
              (HA_detected_HB_detected_SII_detected_OIII_not_detected & mask_S2_classified1 & mask_S2_SF))
# get SFR as non-SF: or HA_detected_HB_detected & mask_S2_classified1 & mask_S2_nonSF
mask_nonSF_S2 = ((HA_detected_HB_detected_SII_detected_OIII_detected & mask_S2_classified1 & mask_S2_nonSF) |
                  (HA_detected_HB_detected_SII_not_detected_OIII_not_detected & mask_S2_classified1 & mask_S2_nonSF) |
                  (HA_detected_HB_detected_SII_not_detected_OIII_detected & mask_S2_classified1 & mask_S2_nonSF) |
                  (HA_detected_HB_detected_SII_detected_OIII_not_detected & mask_S2_classified1 & mask_S2_nonSF))
# all the unconstrained spaxels: or HA_detected_HB_detected & mask_S2_unclassified1
mask_unclassified1_S2_final = ((HA_detected_HB_detected_SII_not_detected_OIII_not_detected & mask_S2_unclassified1) | 
                           (HA_detected_HB_detected_SII_not_detected_OIII_detected & mask_S2_unclassified1) | 
                           (HA_detected_HB_detected_SII_detected_OIII_not_detected & mask_S2_unclassified1) | 
                           (HA_detected_HB_detected_SII_detected_OIII_detected & mask_S2_unclassified1))

# ====== Classification 2: SF = HII only (same as Classification 1 for S2 BPT) ======
# HII spaxels only in S2 BPT 
mask_SF_S2_class2 = ((HA_detected_HB_detected_SII_detected_OIII_detected & mask_classified2_S2 & mask_classified2_S2_SF) | 
                    (HA_detected_HB_detected_SII_not_detected_OIII_not_detected & mask_classified2_S2 & mask_classified2_S2_SF) | 
                    (HA_detected_HB_detected_SII_not_detected_OIII_detected & mask_classified2_S2 & mask_classified2_S2_SF) | 
                    (HA_detected_HB_detected_SII_detected_OIII_not_detected & mask_classified2_S2 & mask_classified2_S2_SF))
# LINER + Seyfert as non-SF in classification 2
mask_nonSF_S2_class2 = ((HA_detected_HB_detected_SII_detected_OIII_detected & mask_classified2_S2 & mask_classified2_S2_nonSF) |
                       (HA_detected_HB_detected_SII_not_detected_OIII_not_detected & mask_classified2_S2 & mask_classified2_S2_nonSF) |
                       (HA_detected_HB_detected_SII_not_detected_OIII_detected & mask_classified2_S2 & mask_classified2_S2_nonSF) |
                       (HA_detected_HB_detected_SII_detected_OIII_not_detected & mask_classified2_S2 & mask_classified2_S2_nonSF))
# all the unclassified2 spaxels (same as unclassified1)
mask_unclassified2_S2_final = ((HA_detected_HB_detected_SII_not_detected_OIII_not_detected & mask_unclassified2_S2) | 
                           (HA_detected_HB_detected_SII_not_detected_OIII_detected & mask_unclassified2_S2) | 
                           (HA_detected_HB_detected_SII_detected_OIII_not_detected & mask_unclassified2_S2) | 
                           (HA_detected_HB_detected_SII_detected_OIII_detected & mask_unclassified2_S2))
# # the rest are upper spaxels: 
# mask_upper = (HA_not_detected | HA_detected_HB_not_detected)

# Something else might be useful

# all the constrained spaxels
mask_classified1_S2 = mask_SF_S2 | mask_nonSF_S2

# ------------------------------------------------------------------
# 8.  Masks: Combine two BPT
# ------------------------------------------------------------------

# ====== Classification 1: both (SF = HII + Comp in N2, HII in S2) ======

# SF: SF in both N2 and S2 BPT diagrams:
mask_SF_both = mask_SF_N2 & mask_SF_S2
# non-SF: constrained in both N2 and S2 BPT diagrams, but not SF in either or both:
mask_nonSF_both = ((mask_classified1_N2 & mask_classified1_S2) & ~mask_SF_both)
# Unclassified1: unconstrained in either N2 or S2 BPT diagrams:
mask_unclassified1_both = ((~(mask_classified1_N2 & mask_classified1_S2)) & HA_detected_HB_detected)
# all the constrained spaxels:
mask_classified1_both = (mask_classified1_N2 & mask_classified1_S2)

# ====== Classification 1: either (SF = HII + Comp in N2, HII in S2) ======

# SF: SF in either N2 or S2 BPT diagrams:
mask_SF_either = mask_SF_N2 | mask_SF_S2
# non-SF: constrained in either N2 or S2 BPT diagrams, but not SF in either or both:
mask_nonSF_either = ((mask_classified1_N2 | mask_classified1_S2) & ~mask_SF_either)
# Unclassified1: unconstrained in either N2 or S2 BPT diagrams:
mask_unclassified1_either = ((~(mask_classified1_N2 | mask_classified1_S2)) & HA_detected_HB_detected)
# all the constrained spaxels:
mask_classified1_either = (mask_classified1_N2 | mask_classified1_S2)

# ====== Classification 2: both (SF = HII only in both N2 and S2) ======

# SF: SF in both N2 and S2 BPT diagrams (HII only):
mask_SF_both_class2 = mask_SF_N2_class2 & mask_SF_S2_class2
# non-SF: constrained in both N2 and S2 BPT diagrams, but not SF:
mask_nonSF_both_class2 = ((mask_classified2_N2 & mask_classified2_S2) & ~mask_SF_both_class2)
# Unclassified2: unconstrained in either N2 or S2 BPT diagrams (same as unclassified1):
mask_unclassified2_both = ((~(mask_classified2_N2 & mask_classified2_S2)) & HA_detected_HB_detected)
# all the constrained spaxels:
mask_classified2_both = (mask_classified2_N2 & mask_classified2_S2)

# ====== Classification 2: either (SF = HII only in either N2 or S2) ======

# SF: SF in either N2 or S2 BPT diagrams (HII only):
mask_SF_either_class2 = mask_SF_N2_class2 | mask_SF_S2_class2
# non-SF: constrained in either N2 or S2 BPT diagrams, but not SF:
mask_nonSF_either_class2 = ((mask_classified2_N2 | mask_classified2_S2) & ~mask_SF_either_class2)
# Unclassified2: unconstrained in either N2 or S2 BPT diagrams (same as unclassified1):
mask_unclassified2_either = ((~(mask_classified2_N2 | mask_classified2_S2)) & HA_detected_HB_detected)
# all the constrained spaxels:
mask_classified2_either = (mask_classified2_N2 | mask_classified2_S2)

# Upper limit spaxels (same for both classifications):
mask_upper = (HA_not_detected | HA_detected_HB_not_detected)

# Independent exact-class BPT maps.
# Only full line detections whose central and +/-1 sigma BPT positions remain
# in the same class are assigned a positive class code. Full detections that
# are ambiguous remain 0; non-detections or invalid ratios remain -1.
NII_BPT = np.full_like(HA6562_FLUX, -1, dtype=np.int16)
SII_BPT = np.full_like(HA6562_FLUX, -1, dtype=np.int16)

mask_NII_BPT_valid = (
    HA_detected_HB_detected_NII_detected_OIII_detected
    & np.isfinite(logN2) & np.isfinite(logO3)
    & np.isfinite(logN2_err) & np.isfinite(logO3_err)
)
NII_BPT[mask_NII_BPT_valid] = 0
mask_NII_BPT_HII = (
    mask_NII_BPT_valid & mask_N2_HII
    & mask_N2_left_HII & mask_N2_right_HII
    & mask_N2_down_HII & mask_N2_up_HII
)
mask_NII_BPT_Comp = (
    mask_NII_BPT_valid & mask_N2_Comp
    & mask_N2_left_Comp & mask_N2_right_Comp
    & mask_N2_down_Comp & mask_N2_up_Comp
)
mask_NII_BPT_AGN = (
    mask_NII_BPT_valid & mask_N2_AGN
    & mask_N2_left_AGN & mask_N2_right_AGN
    & mask_N2_down_AGN & mask_N2_up_AGN
)
NII_BPT[mask_NII_BPT_HII] = 1
NII_BPT[mask_NII_BPT_Comp] = 2
NII_BPT[mask_NII_BPT_AGN] = 3

mask_SII_BPT_valid = (
    HA_detected_HB_detected_SII_detected_OIII_detected
    & np.isfinite(logS2) & np.isfinite(logO3)
    & np.isfinite(logS2_err) & np.isfinite(logO3_err)
)
SII_BPT[mask_SII_BPT_valid] = 0
mask_SII_BPT_HII = (
    mask_SII_BPT_valid & mask_S2_HII
    & mask_S2_left_HII & mask_S2_right_HII
    & mask_S2_down_HII & mask_S2_up_HII
)
mask_SII_BPT_LINER = (
    mask_SII_BPT_valid & mask_S2_LINER
    & mask_S2_left_LINER & mask_S2_right_LINER
    & mask_S2_down_LINER & mask_S2_up_LINER
)
mask_SII_BPT_Seyfert = (
    mask_SII_BPT_valid & mask_S2_Seyfert
    & mask_S2_left_Seyfert & mask_S2_right_Seyfert
    & mask_S2_down_Seyfert & mask_S2_up_Seyfert
)
SII_BPT[mask_SII_BPT_HII] = 1
SII_BPT[mask_SII_BPT_LINER] = 2
SII_BPT[mask_SII_BPT_Seyfert] = 3

# ------------------------------------------------------------------
# 9.  Append Σ_SFR layers (choose 'both' for now, but can be changed to 'either' or just fall back to 'N2' or 'S2')
# ------------------------------------------------------------------

def choose_BPT(choice='both', classification=1):
    """
    Choose BPT classification method and return SFR surface density maps, metallicity maps, line maps, and masks.
    
    Parameters:
    choice : str, optional
        BPT classification choice. Options: 'both', 'either', 'N2', 'S2' (default: 'both')
    classification : int, optional
        Classification scheme: 1 = SF includes HII+Comp, 2 = SF includes HII only (default: 1)
        
    Returns:
    tuple : (SFR_maps, metallicity_maps, line_maps, masks) where:
        - SFR_maps: tuple of four arrays (SF/HII, nonSF/nonHII, unconstrained, upper) for the chosen method
        - metallicity_maps: tuple of thirteen arrays for SF/HII regions only
        - line_maps: tuple of six arrays for SF/HII regions only
        - masks: tuple of four boolean arrays (mask_SF/HII, mask_nonSF/nonHII, mask_unclassified1, mask_upper)
    """
    # Get the appropriate masks based on choice and classification
    if classification == 1:
        # Classification 1: SF = HII + Comp in N2, HII in S2
        if choice == 'both':
            mask_SF = mask_SF_both
            mask_nonSF = mask_nonSF_both
            mask_unclassified1 = mask_unclassified1_both
        elif choice == 'either':
            mask_SF = mask_SF_either
            mask_nonSF = mask_nonSF_either
            mask_unclassified1 = mask_unclassified1_either
        elif choice == 'N2':
            mask_SF = mask_SF_N2
            mask_nonSF = mask_nonSF_N2
            mask_unclassified1 = mask_unclassified1_N2_final
        elif choice == 'S2':
            mask_SF = mask_SF_S2
            mask_nonSF = mask_nonSF_S2
            mask_unclassified1 = mask_unclassified1_S2_final
        else:
            raise ValueError(f"Invalid choice '{choice}'. Options: 'both', 'either', 'N2', 'S2'")
    elif classification == 2:
        # Classification 2: SF = HII only in both N2 and S2
        if choice == 'both':
            mask_SF = mask_SF_both_class2
            mask_nonSF = mask_nonSF_both_class2
            mask_unclassified1 = mask_unclassified2_both
        elif choice == 'either':
            mask_SF = mask_SF_either_class2
            mask_nonSF = mask_nonSF_either_class2
            mask_unclassified1 = mask_unclassified2_either
        elif choice == 'N2':
            mask_SF = mask_SF_N2_class2
            mask_nonSF = mask_nonSF_N2_class2
            mask_unclassified1 = mask_unclassified2_N2_final
        elif choice == 'S2':
            mask_SF = mask_SF_S2_class2
            mask_nonSF = mask_nonSF_S2_class2
            mask_unclassified1 = mask_unclassified2_S2_final
        else:
            raise ValueError(f"Invalid choice '{choice}'. Options: 'both', 'either', 'N2', 'S2'")
    else:
        raise ValueError(f"Invalid classification '{classification}'. Options: 1, 2")
    
    # Apply masks to create SFR surface density maps
    LOG_SFR_surface_density_map_SF = np.where(mask_SF, LOG_SFR_surface_density_map, np.nan)
    LOG_SFR_surface_density_map_nonSF = np.where(mask_nonSF, LOG_SFR_surface_density_map, np.nan)
    LOG_SFR_surface_density_map_unclassified1 = np.where(mask_unclassified1, LOG_SFR_surface_density_map, np.nan)
    LOG_SFR_surface_density_map_upper = np.where(mask_upper, LOG_SFR_surface_density_map, np.nan)
    
    # Apply masks to create non-log SFR surface-density maps.
    SFR_map_SF = np.where(mask_SF, SFR_surface_density_map_corrected, np.nan)
    SFR_map_nonSF = np.where(mask_nonSF, SFR_surface_density_map_corrected, np.nan)
    SFR_map_unclassified1 = np.where(mask_unclassified1, SFR_surface_density_map_corrected, np.nan)
    SFR_map_upper = np.where(mask_upper, SFR_surface_density_map_corrected, np.nan)
    
    # Apply SF mask to create metallicity maps (only for SF regions)
    O_H_D16_SF = np.where(mask_SF, O_H_D16, np.nan)
    O_H_PG16_SF = np.where(mask_SF, O_H_PG16, np.nan)
    O_H_N2S2_N06_SF = np.where(mask_SF, O_H_N2S2_N06, np.nan)
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
    HB4861_FLUX_ERR_SF = np.where(mask_SF, HB4861_FLUX_ERR_corr, np.nan)
    HA6562_FLUX_ERR_SF = np.where(mask_SF, HA6562_FLUX_ERR_corr, np.nan)
    OIII5006_FLUX_ERR_SF = np.where(mask_SF, OIII5006_FLUX_ERR_corr, np.nan)
    NII6583_FLUX_ERR_SF = np.where(mask_SF, NII6583_FLUX_ERR_corr, np.nan)
    SII6716_FLUX_ERR_SF = np.where(mask_SF, SII6716_FLUX_ERR_corr, np.nan)
    SII6730_FLUX_ERR_SF = np.where(mask_SF, SII6730_FLUX_ERR_corr, np.nan)

    # Return SFR maps, metallicity maps, metallicity error maps, line maps, and masks
    sfr_maps = (LOG_SFR_surface_density_map_SF, LOG_SFR_surface_density_map_nonSF, 
                LOG_SFR_surface_density_map_unclassified1, LOG_SFR_surface_density_map_upper)
    sfr_maps_regular = (SFR_map_SF, SFR_map_nonSF, SFR_map_unclassified1, SFR_map_upper)
    metallicity_maps = (O_H_D16_SF, O_H_PG16_SF, O_H_N2S2_N06_SF, O_H_O3N2_M13_SF, O_H_N2_M13_SF, O_H_O3N2_PP04_SF, O_H_N2_PP04_SF, O_H_O3N2_C20_SF, O_H_O3S2_C20_SF, O_H_RS32_C20_SF, O_H_R3_C20_SF, O_H_N2_C20_SF, O_H_S2_C20_SF, O_H_COMBINED_C20_SF)
    metallicity_error_maps = (O_H_O3N2_C20_SF_ERR, O_H_O3S2_C20_SF_ERR, O_H_RS32_C20_SF_ERR, O_H_R3_C20_SF_ERR, O_H_N2_C20_SF_ERR, O_H_S2_C20_SF_ERR, O_H_COMBINED_C20_SF_ERR)
    line_maps = (HB4861_FLUX_corr_SF, HA6562_FLUX_corr_SF, OIII5006_FLUX_corr_SF, 
                 NII6583_FLUX_corr_SF, SII6716_FLUX_corr_SF, SII6730_FLUX_corr_SF)
    masks = (mask_SF, mask_nonSF, mask_unclassified1, mask_upper)
    
    return sfr_maps, sfr_maps_regular, metallicity_maps, metallicity_error_maps, line_maps, masks

# Get the SFR surface density maps, metallicity maps, metallicity error maps, line maps, and masks using the default 'both' choice
# Classification 1: SF = HII + Comp
(LOG_SFR_surface_density_map_SF, LOG_SFR_surface_density_map_nonSF, 
 LOG_SFR_surface_density_map_unclassified1, LOG_SFR_surface_density_map_upper), (SFR_map_SF, SFR_map_nonSF, SFR_map_unclassified1, SFR_map_upper), (O_H_D16_SF, O_H_PG16_SF, O_H_N2S2_N06_SF, O_H_O3N2_M13_SF, O_H_N2_M13_SF, O_H_O3N2_PP04_SF, O_H_N2_PP04_SF, O_H_O3N2_C20_SF, O_H_O3S2_C20_SF, O_H_RS32_C20_SF, O_H_R3_C20_SF, O_H_N2_C20_SF, O_H_S2_C20_SF, O_H_COMBINED_C20_SF), (O_H_O3N2_C20_SF_ERR, O_H_O3S2_C20_SF_ERR, O_H_RS32_C20_SF_ERR, O_H_R3_C20_SF_ERR, O_H_N2_C20_SF_ERR, O_H_S2_C20_SF_ERR, O_H_COMBINED_C20_SF_ERR), (HB4861_FLUX_corr_SF, HA6562_FLUX_corr_SF, OIII5006_FLUX_corr_SF, NII6583_FLUX_corr_SF, SII6716_FLUX_corr_SF, SII6730_FLUX_corr_SF), (mask_SF, mask_nonSF, mask_unclassified1, mask_upper) = choose_BPT()

# Classification 2: HII only
(LOG_SFR_surface_density_map_HII, LOG_SFR_surface_density_map_nonHII, 
 LOG_SFR_surface_density_map_unclassified2, LOG_SFR_surface_density_map_upper_HII), (SFR_map_HII, SFR_map_nonHII, SFR_map_unclassified2, SFR_map_upper_HII), (O_H_D16_HII, O_H_PG16_HII, O_H_N2S2_N06_HII, O_H_O3N2_M13_HII, O_H_N2_M13_HII, O_H_O3N2_PP04_HII, O_H_N2_PP04_HII, O_H_O3N2_C20_HII, O_H_O3S2_C20_HII, O_H_RS32_C20_HII, O_H_R3_C20_HII, O_H_N2_C20_HII, O_H_S2_C20_HII, O_H_COMBINED_C20_HII), (O_H_O3N2_C20_HII_ERR, O_H_O3S2_C20_HII_ERR, O_H_RS32_C20_HII_ERR, O_H_R3_C20_HII_ERR, O_H_N2_C20_HII_ERR, O_H_S2_C20_HII_ERR, O_H_COMBINED_C20_HII_ERR), (HB4861_FLUX_corr_HII, HA6562_FLUX_corr_HII, OIII5006_FLUX_corr_HII, NII6583_FLUX_corr_HII, SII6716_FLUX_corr_HII, SII6730_FLUX_corr_HII), (mask_HII, mask_nonHII, mask_unclassified2, mask_upper_HII) = choose_BPT(classification=2)

# HII/SF views for every detected gas-line correction product.  The legacy
# choose_BPT return keeps the historical six-line names; this loop extends the
# same convention to optional lines and to corrected line-error maps.
for line_base in GAS_LINE_BASES:
    for product_suffix in ("FLUX_corr", "FLUX_ERR_corr"):
        source_name = f"{line_base}_{product_suffix}"
        source_map = globals()[source_name]
        globals()[f"{source_name}_HII"] = np.where(mask_HII, source_map, np.nan)
        globals()[f"{source_name}_SF"] = np.where(mask_SF, source_map, np.nan)

# ------------------------------------------------------------------
# 9b. pyqz, NebulaBayes, and JY22 HII-region model-grid products
# ------------------------------------------------------------------

MODEL_OUTPUT_MAPS = {}
MODEL_INTEGRATED_RESULTS = {}
PYQZ_PACKAGE_VERSION = None
PYQZ_GRID_ID = None
NB_PACKAGE_VERSION = None
NB_GRID_PATH = None
JY22_GRID = None


def model_progress_callback(label):
    """Build a concise progress reporter for a potentially long model pass."""
    started = time.perf_counter()

    def report(done, total, bin_id):
        interval = max(1, int(total) // 10)
        if done == 1 or done == total or done % interval == 0:
            elapsed = time.perf_counter() - started
            print(
                f"{label}: {done}/{total} spectra complete "
                f"(BINID {bin_id}, {elapsed:.1f} s)"
            )

    return report


def print_model_failures(label, failures):
    """Report every failed identifier compactly, grouped by exception text."""
    if not failures:
        print(f"{label}: no per-spectrum exceptions")
        return
    grouped = {}
    for bin_id, message in failures:
        grouped.setdefault(message, []).append(bin_id)
    print(f"{label}: {len(failures)} per-spectrum exception(s)")
    for message, bin_ids in grouped.items():
        for start in range(0, len(bin_ids), 25):
            identifiers = ",".join(
                str(bin_id) for bin_id in bin_ids[start : start + 25]
            )
            suffix = f": {message}" if start == 0 else " (same exception)"
            print(f"  BINID(s) {identifiers}{suffix}")


def selected_model_field(source_map, selection, *, integer=False):
    """Apply an HII/SF mask without changing the underlying shared fit."""
    if integer:
        output = np.full(source_map.shape, NOT_EVALUATED_FLAG, dtype=np.int16)
        output[selection] = np.asarray(source_map, dtype=np.int16)[selection]
        return output
    return np.where(selection, np.asarray(source_map, dtype=np.float64), np.nan)


def store_pyqz_region_maps(region, region_mask, broadcast_maps, input_valid):
    selection = np.asarray(region_mask, dtype=bool) & input_valid
    float_specs = {
        f"O_H_PYQZ_{region}": "o_h",
        f"O_H_PYQZ_{region}_ERR": "o_h_err",
        f"LOG_Q_PYQZ_{region}": "log_q",
        f"LOG_Q_PYQZ_{region}_ERR": "log_q_err",
        f"LOG_U_PYQZ_{region}": "log_u",
        f"LOG_U_PYQZ_{region}_ERR": "log_u_err",
        f"PYQZ_RS_OFFGRID_{region}": "rs_offgrid",
    }
    integer_specs = {
        f"PYQZ_FLAG_{region}": "flag",
        f"PYQZ_VALID_{region}": "valid",
    }
    for output_name, field in float_specs.items():
        MODEL_OUTPUT_MAPS[output_name] = selected_model_field(
            broadcast_maps[field], selection
        )
    for output_name, field in integer_specs.items():
        MODEL_OUTPUT_MAPS[output_name] = selected_model_field(
            broadcast_maps[field], selection, integer=True
        )


def store_nebulabayes_region_maps(region, region_mask, broadcast_maps, input_valid):
    selection = np.asarray(region_mask, dtype=bool) & input_valid
    float_specs = {
        f"O_H_NEBULABAYES_{region}": "o_h",
        f"O_H_NEBULABAYES_{region}_CI68_LOW": "o_h_ci68_low",
        f"O_H_NEBULABAYES_{region}_CI68_HIGH": "o_h_ci68_high",
        f"LOG_U_NEBULABAYES_{region}": "log_u",
        f"LOG_U_NEBULABAYES_{region}_CI68_LOW": "log_u_ci68_low",
        f"LOG_U_NEBULABAYES_{region}_CI68_HIGH": "log_u_ci68_high",
        f"LOG_PK_NEBULABAYES_{region}": "log_pk",
        f"LOG_PK_NEBULABAYES_{region}_CI68_LOW": "log_pk_ci68_low",
        f"LOG_PK_NEBULABAYES_{region}_CI68_HIGH": "log_pk_ci68_high",
        f"NB_CHI2_RED_{region}": "chi2_red",
    }
    integer_specs = {
        f"NB_NLOCALMAX_{region}": "n_localmax",
        f"NB_FLAG_{region}": "flag",
        f"NB_VALID_{region}": "valid",
    }
    for output_name, field in float_specs.items():
        MODEL_OUTPUT_MAPS[output_name] = selected_model_field(
            broadcast_maps[field], selection
        )
    for output_name, field in integer_specs.items():
        MODEL_OUTPUT_MAPS[output_name] = selected_model_field(
            broadcast_maps[field], selection, integer=True
        )


def store_jy22_region_maps(region, region_mask, broadcast_maps, input_valid):
    """Store JY22 posterior summaries and QC on one existing region mask."""
    selection = np.asarray(region_mask, dtype=bool) & input_valid
    float_specs = {
        f"O_H_JY22_{region}": "o_h",
        f"O_H_JY22_{region}_16": "o_h_16",
        f"O_H_JY22_{region}_84": "o_h_84",
        f"LOG_U_JY22_{region}": "log_u",
        f"LOG_U_JY22_{region}_16": "log_u_16",
        f"LOG_U_JY22_{region}_84": "log_u_84",
        f"JY22_CHI2_MIN_{region}": "chi2_min",
    }
    integer_specs = {
        f"JY22_FLAG_{region}": "flag",
        f"JY22_VALID_{region}": "valid",
    }
    for output_name, field in float_specs.items():
        MODEL_OUTPUT_MAPS[output_name] = selected_model_field(
            broadcast_maps[field], selection
        )
    for output_name, field in integer_specs.items():
        MODEL_OUTPUT_MAPS[output_name] = selected_model_field(
            broadcast_maps[field], selection, integer=True
        )


def corrected_integrated_model_spectrum(integrated):
    """Apply one integrated Balmer-decrement correction to six raw sums."""
    if integrated.n_pixels == 0:
        return None
    hb_flux = integrated.fluxes["HB4861"]
    ha_flux = integrated.fluxes["HA6562"]
    hb_error = integrated.errors["HB4861"]
    ha_error = integrated.errors["HA6562"]
    balmer_decrement = max(ha_flux / hb_flux, R_int)
    ebv = float(
        np.asarray(
            convert_bd_to_ebv(
                balmer_decrement, k_HB4861, k_HA6562, R_int
            )
        )
    )
    ebv_error = float(
        np.asarray(
            convert_bd_to_ebv_error(
                ha_flux,
                hb_flux,
                ha_error,
                hb_error,
                k_HB4861,
                k_HA6562,
            )
        )
    )
    corrected_fluxes = []
    corrected_errors = []
    for line_base in MODEL_LINE_BASES:
        raw_flux = integrated.fluxes[line_base]
        raw_error = integrated.errors[line_base]
        k_line = GAS_LINE_K[line_base]
        corrected_fluxes.append(
            float(np.asarray(correct_flux_with_ebv(raw_flux, ebv, k_line)))
        )
        corrected_errors.append(
            float(
                np.asarray(
                    correct_flux_error_with_ebv(
                        raw_flux,
                        raw_error,
                        ebv,
                        k_line,
                        ebv_error,
                    )
                )
            )
        )
    fluxes = np.asarray(corrected_fluxes, dtype=np.float64)
    errors = np.asarray(corrected_errors, dtype=np.float64)
    if not np.all(np.isfinite(fluxes) & (fluxes > 0.0)):
        return None
    if not np.all(np.isfinite(errors) & (errors > 0.0)):
        return None
    return {
        "fluxes": fluxes,
        "errors": errors,
        "balmer_decrement": balmer_decrement,
        "ebv": ebv,
        "ebv_error": ebv_error,
    }


def record_from_batch(batch_results, index):
    return {name: values[index].item() for name, values in batch_results.items()}


if (
    ENABLE_PYQZ_METALLICITY_PRODUCTS
    or ENABLE_NEBULABAYES_METALLICITY_PRODUCTS
    or ENABLE_JY22_METALLICITY_PRODUCTS
):
    model_corrected_flux_maps = {
        line_base: globals()[f"{line_base}_FLUX_corr"]
        for line_base in MODEL_LINE_BASES
    }
    model_corrected_error_maps = {
        line_base: globals()[f"{line_base}_FLUX_ERR_corr"]
        for line_base in MODEL_LINE_BASES
    }
    model_raw_flux_maps = {
        line_base: globals()[f"{line_base}_FLUX"]
        for line_base in MODEL_LINE_BASES
    }
    model_raw_error_maps = {
        line_base: globals()[f"{line_base}_FLUX_ERR"]
        for line_base in MODEL_LINE_BASES
    }
    model_region_union = np.asarray(mask_HII | mask_SF, dtype=bool)
    model_input_valid = build_model_input_valid_mask(
        model_corrected_flux_maps,
        model_corrected_error_maps,
        model_region_union,
    )
    model_bin_spectra = extract_unique_bin_spectra(
        gas_binid_map,
        model_input_valid,
        model_corrected_flux_maps,
        model_corrected_error_maps,
    )
    model_raw_bin_spectra = None
    if ENABLE_JY22_METALLICITY_PRODUCTS:
        model_raw_bin_spectra = extract_unique_bin_spectra(
            gas_binid_map,
            model_input_valid,
            model_raw_flux_maps,
            model_raw_error_maps,
        )
        if not (
            np.array_equal(model_raw_bin_spectra.bin_ids, model_bin_spectra.bin_ids)
            and np.array_equal(
                model_raw_bin_spectra.pixel_counts, model_bin_spectra.pixel_counts
            )
        ):
            raise RuntimeError(
                "Raw and corrected model spectra do not describe identical BINIDs."
            )
    print("--------------------------------------------------------------")
    print("HII-region model-grid inference")
    print(
        "Common six-line corrected-input aperture: "
        f"{np.count_nonzero(model_input_valid)} pixels in "
        f"{model_bin_spectra.bin_ids.size} unique spatial bins; no extra S/N cut"
    )
    print(f"Python interpreter: {sys.executable}")
    if ENABLE_PYQZ_METALLICITY_PRODUCTS or ENABLE_NEBULABAYES_METALLICITY_PRODUCTS:
        print("pyqz/NebulaBayes compatibility layer: model_grid_compat.py")

    pyqz_runtime = None
    pyqz_bin_run = None
    if ENABLE_PYQZ_METALLICITY_PRODUCTS:
        PYQZ_PACKAGE_VERSION = importlib_metadata.version("pyqz")
        pyqz_runtime = load_pyqz_compat()
        _, _, pyqz_grid_metadata = pyqz_runtime.get_grid(
            PYQZ_DIAGNOSTIC,
            Pk=PYQZ_PK,
            kappa=PYQZ_KAPPA,
            struct=PYQZ_STRUCT,
            sampling=PYQZ_SAMPLING,
        )
        PYQZ_GRID_ID = str(
            pyqz_grid_metadata.get("MV_id", "bundled MAPPINGS V grid")
        )
        print(
            f"pyqz {PYQZ_PACKAGE_VERSION}: diagnostic={PYQZ_DIAGNOSTIC}; "
            f"Pk={PYQZ_PK}; struct={PYQZ_STRUCT}; kappa={PYQZ_KAPPA}; "
            f"sampling={PYQZ_SAMPLING}; srs={PYQZ_SRS}; "
            f"KDE={PYQZ_KDE_METHOD}; seed={PYQZ_RANDOM_SEED}; "
            f"nproc={PYQZ_NPROC}"
        )
        print(f"pyqz runtime overlay: {pyqz_runtime.__file__}")
        print(f"pyqz grid identity: {PYQZ_GRID_ID}")
        pyqz_bin_run = run_pyqz_spectra(
            pyqz_runtime,
            model_bin_spectra,
            diagnostic=PYQZ_DIAGNOSTIC,
            qzs=PYQZ_QZS,
            pk=PYQZ_PK,
            kappa=PYQZ_KAPPA,
            struct=PYQZ_STRUCT,
            sampling=PYQZ_SAMPLING,
            error_pdf=PYQZ_ERROR_PDF,
            random_seed=PYQZ_RANDOM_SEED,
            srs=PYQZ_SRS,
            flag_level=PYQZ_FLAG_LEVEL,
            kde_method=PYQZ_KDE_METHOD,
            kde_qz_sampling=PYQZ_KDE_QZ_SAMPLING,
            kde_do_singles=PYQZ_KDE_DO_SINGLES,
            nproc=PYQZ_NPROC,
            progress=model_progress_callback("pyqz bins"),
        )
        print_model_failures("pyqz bins", pyqz_bin_run.failures)
        if (
            model_bin_spectra.bin_ids.size > 0
            and not np.any(pyqz_bin_run.results["valid"] == 1)
        ):
            raise RuntimeError("pyqz failed for every eligible spatial bin.")
        pyqz_broadcast = broadcast_bin_results(
            gas_binid_map,
            model_bin_spectra.bin_ids,
            pyqz_bin_run.results,
            integer_fields=set(PYQZ_INTEGER_FIELDS),
        )
        store_pyqz_region_maps(
            "HII", mask_HII, pyqz_broadcast, model_input_valid
        )
        store_pyqz_region_maps(
            "SF", mask_SF, pyqz_broadcast, model_input_valid
        )

    nebulabayes_model = None
    nebulabayes_bin_run = None
    if ENABLE_NEBULABAYES_METALLICITY_PRODUCTS:
        NB_PACKAGE_VERSION = importlib_metadata.version("NebulaBayes")
        _, nb_grid_parameters, NB_GRID_PATH = load_nebulabayes_hii_grid()
        nebulabayes_model = make_nebulabayes_hii_model(
            NB_LINE_LIST,
            interpd_grid_shape=NB_INTERPD_GRID_SHAPE,
            interp_order=NB_INTERP_ORDER,
            grid_error=NB_GRID_ERROR,
        )
        print(
            f"NebulaBayes {NB_PACKAGE_VERSION}: grid={NB_GRID_PATH}; "
            f"parameters={nb_grid_parameters}; shape={NB_INTERPD_GRID_SHAPE}; "
            f"order={NB_INTERP_ORDER}; grid_error={NB_GRID_ERROR}; "
            f"norm={NB_NORM_LINE}; deredden={NB_DEREDDEN}; prior={NB_PRIOR}"
        )
        nebulabayes_bin_run = run_nebulabayes_spectra(
            nebulabayes_model,
            model_bin_spectra,
            norm_line=NB_NORM_LINE,
            deredden=NB_DEREDDEN,
            prior=NB_PRIOR,
            likelihood_lines=NB_LIKELIHOOD_LINES,
            progress=model_progress_callback("NebulaBayes bins"),
        )
        print_model_failures(
            "NebulaBayes bins", nebulabayes_bin_run.failures
        )
        if (
            model_bin_spectra.bin_ids.size > 0
            and not np.any(nebulabayes_bin_run.results["valid"] == 1)
        ):
            raise RuntimeError(
                "NebulaBayes failed for every eligible spatial bin."
            )
        nebulabayes_broadcast = broadcast_bin_results(
            gas_binid_map,
            model_bin_spectra.bin_ids,
            nebulabayes_bin_run.results,
            integer_fields=set(NB_INTEGER_FIELDS),
        )
        store_nebulabayes_region_maps(
            "HII", mask_HII, nebulabayes_broadcast, model_input_valid
        )
        store_nebulabayes_region_maps(
            "SF", mask_SF, nebulabayes_broadcast, model_input_valid
        )

    jy22_bin_run = None
    if ENABLE_JY22_METALLICITY_PRODUCTS:
        if model_raw_bin_spectra is None:
            raise RuntimeError("JY22 raw spatial spectra were not prepared.")
        if not np.isclose(JY22_ERROR_INFLATION, 1.0, rtol=0.0, atol=0.0):
            raise RuntimeError(
                "JY22_ERROR_INFLATION is scientifically locked to 1.0; the "
                "MaNGA-specific factor 1.25 is not adopted."
            )
        JY22_GRID = load_jy22_grid(
            JY22_GRID_PATH,
            min_log_u=JY22_LOG_U_MIN,
            max_log_u=JY22_LOG_U_MAX,
            solar_o_h=JY22_SOLAR_O_H,
        )
        if JY22_GRID.sha256 != JY22_EXPECTED_GRID_SHA256:
            raise RuntimeError(
                "JY22 grid SHA-256 differs from the released grid locked by "
                f"this pipeline: {JY22_GRID.sha256}"
            )
        jy22_extinction_coefficients = np.asarray(
            [GAS_LINE_K[line_base] for line_base in MODEL_LINE_BASES],
            dtype=np.float64,
        )
        jy22_ratio_spectra = build_jy22_ratio_spectra(
            model_raw_bin_spectra,
            jy22_extinction_coefficients,
            intrinsic_balmer_ratio=JY22_INTRINSIC_BALMER_RATIO,
        )
        print(
            "JY22/Peng2026: ratios=N2,S2,R3; posterior means with marginal "
            "equal-tailed p16/p84; full raw-flux ratio covariance; "
            "error inflation=1.0"
        )
        print(
            f"JY22 grid: {JY22_GRID.source_path}; "
            f"sha256={JY22_GRID.sha256}; "
            f"requested logU={JY22_GRID.requested_log_u_bounds}; "
            f"retained logU={JY22_GRID.effective_log_u_bounds}; "
            f"solar O/H={JY22_GRID.solar_o_h}"
        )
        jy22_bin_run = run_jy22_spectra(
            JY22_GRID,
            jy22_ratio_spectra,
            progress=model_progress_callback("JY22 bins"),
        )
        print_model_failures("JY22 bins", jy22_bin_run.failures)
        if (
            model_raw_bin_spectra.bin_ids.size > 0
            and not np.any(jy22_bin_run.results["valid"] == 1)
        ):
            raise RuntimeError("JY22 failed for every eligible spatial bin.")
        jy22_broadcast = broadcast_bin_results(
            gas_binid_map,
            model_raw_bin_spectra.bin_ids,
            jy22_bin_run.results,
            integer_fields=set(JY22_INTEGER_FIELDS),
        )
        store_jy22_region_maps(
            "HII", mask_HII, jy22_broadcast, model_input_valid
        )
        store_jy22_region_maps(
            "SF", mask_SF, jy22_broadcast, model_input_valid
        )

    for region, region_mask in (("HII", mask_HII), ("SF", mask_SF)):
        integrated = integrate_line_maps(
            model_raw_flux_maps,
            model_raw_error_maps,
            region_mask,
            binid_map=gas_binid_map,
        )
        MODEL_INTEGRATED_RESULTS[region] = {
            "aperture": integrated,
            "spectrum": corrected_integrated_model_spectrum(integrated),
            "pyqz": empty_pyqz_result(),
            "nebulabayes": empty_nebulabayes_result(),
            "jy22": empty_jy22_result(),
            "jy22_ebv": np.nan,
        }

    integrated_regions = [
        region
        for region in ("HII", "SF")
        if MODEL_INTEGRATED_RESULTS[region]["spectrum"] is not None
    ]
    if integrated_regions:
        integrated_spectra = BinSpectra(
            bin_ids=np.arange(len(integrated_regions), dtype=np.int64),
            fluxes=np.asarray(
                [
                    MODEL_INTEGRATED_RESULTS[region]["spectrum"]["fluxes"]
                    for region in integrated_regions
                ],
                dtype=np.float64,
            ),
            errors=np.asarray(
                [
                    MODEL_INTEGRATED_RESULTS[region]["spectrum"]["errors"]
                    for region in integrated_regions
                ],
                dtype=np.float64,
            ),
            pixel_counts=np.asarray(
                [
                    MODEL_INTEGRATED_RESULTS[region]["aperture"].n_pixels
                    for region in integrated_regions
                ],
                dtype=np.int64,
            ),
        )
        if ENABLE_PYQZ_METALLICITY_PRODUCTS:
            for index, region in enumerate(integrated_regions):
                # Reset the documented seed for each integrated aperture.  If
                # HII and SF select exactly the same spectrum, their pyqz
                # result is then identical rather than differing by MC noise.
                region_spectrum = BinSpectra(
                    bin_ids=np.asarray([index], dtype=np.int64),
                    fluxes=integrated_spectra.fluxes[index : index + 1],
                    errors=integrated_spectra.errors[index : index + 1],
                    pixel_counts=integrated_spectra.pixel_counts[index : index + 1],
                )
                integrated_pyqz_run = run_pyqz_spectra(
                    pyqz_runtime,
                    region_spectrum,
                    diagnostic=PYQZ_DIAGNOSTIC,
                    qzs=PYQZ_QZS,
                    pk=PYQZ_PK,
                    kappa=PYQZ_KAPPA,
                    struct=PYQZ_STRUCT,
                    sampling=PYQZ_SAMPLING,
                    error_pdf=PYQZ_ERROR_PDF,
                    random_seed=PYQZ_RANDOM_SEED,
                    srs=PYQZ_SRS,
                    flag_level=PYQZ_FLAG_LEVEL,
                    kde_method=PYQZ_KDE_METHOD,
                    kde_qz_sampling=PYQZ_KDE_QZ_SAMPLING,
                    kde_do_singles=PYQZ_KDE_DO_SINGLES,
                    nproc=PYQZ_NPROC,
                )
                print_model_failures(
                    f"pyqz integrated {region}", integrated_pyqz_run.failures
                )
                MODEL_INTEGRATED_RESULTS[region]["pyqz"] = record_from_batch(
                    integrated_pyqz_run.results, 0
                )
        if ENABLE_NEBULABAYES_METALLICITY_PRODUCTS:
            integrated_nb_run = run_nebulabayes_spectra(
                nebulabayes_model,
                integrated_spectra,
                norm_line=NB_NORM_LINE,
                deredden=NB_DEREDDEN,
                prior=NB_PRIOR,
                likelihood_lines=NB_LIKELIHOOD_LINES,
            )
            print_model_failures(
                "NebulaBayes integrated", integrated_nb_run.failures
            )
            for index, region in enumerate(integrated_regions):
                MODEL_INTEGRATED_RESULTS[region][
                    "nebulabayes"
                ] = record_from_batch(integrated_nb_run.results, index)

    if ENABLE_JY22_METALLICITY_PRODUCTS:
        jy22_integrated_regions = [
            region
            for region in ("HII", "SF")
            if MODEL_INTEGRATED_RESULTS[region]["aperture"].n_pixels > 0
        ]
        if jy22_integrated_regions:
            jy22_integrated_raw_spectra = BinSpectra(
                bin_ids=np.arange(len(jy22_integrated_regions), dtype=np.int64),
                fluxes=np.asarray(
                    [
                        [
                            MODEL_INTEGRATED_RESULTS[region]["aperture"].fluxes[
                                line_base
                            ]
                            for line_base in MODEL_LINE_BASES
                        ]
                        for region in jy22_integrated_regions
                    ],
                    dtype=np.float64,
                ),
                errors=np.asarray(
                    [
                        [
                            MODEL_INTEGRATED_RESULTS[region]["aperture"].errors[
                                line_base
                            ]
                            for line_base in MODEL_LINE_BASES
                        ]
                        for region in jy22_integrated_regions
                    ],
                    dtype=np.float64,
                ),
                pixel_counts=np.asarray(
                    [
                        MODEL_INTEGRATED_RESULTS[region]["aperture"].n_pixels
                        for region in jy22_integrated_regions
                    ],
                    dtype=np.int64,
                ),
            )
            jy22_integrated_ratios = build_jy22_ratio_spectra(
                jy22_integrated_raw_spectra,
                jy22_extinction_coefficients,
                intrinsic_balmer_ratio=JY22_INTRINSIC_BALMER_RATIO,
            )
            jy22_integrated_run = run_jy22_spectra(
                JY22_GRID, jy22_integrated_ratios
            )
            print_model_failures(
                "JY22 integrated", jy22_integrated_run.failures
            )
            for index, region in enumerate(jy22_integrated_regions):
                MODEL_INTEGRATED_RESULTS[region]["jy22"] = record_from_batch(
                    jy22_integrated_run.results, index
                )
                MODEL_INTEGRATED_RESULTS[region]["jy22_ebv"] = float(
                    jy22_integrated_ratios.ebv[index]
                )

    print("Integrated HII/SF model-grid results (common six-line raw aperture)")
    for region in ("HII", "SF"):
        summary = MODEL_INTEGRATED_RESULTS[region]
        aperture = summary["aperture"]
        spectrum = summary["spectrum"]
        print(
            f"  {region}: {aperture.n_pixels} pixels, {aperture.n_bins} bins; "
            f"integrated E(B-V)="
            f"{spectrum['ebv']:.5f}+/-{spectrum['ebv_error']:.5f}"
            if spectrum is not None
            else f"  {region}: 0 valid common-aperture pixels"
        )
        if ENABLE_PYQZ_METALLICITY_PRODUCTS:
            result = summary["pyqz"]
            print(
                f"    pyqz: O/H={result['o_h']:.5f}+/-{result['o_h_err']:.5f}; "
                f"LogQ={result['log_q']:.5f}+/-{result['log_q_err']:.5f}; "
                f"logU={result['log_u']:.5f}+/-{result['log_u_err']:.5f}; "
                f"flag={result['flag']}; off-grid={result['rs_offgrid']:.2f}%; "
                f"valid={result['valid']} (uncertainties: 0.61 KDE contour)"
            )
        if ENABLE_NEBULABAYES_METALLICITY_PRODUCTS:
            result = summary["nebulabayes"]
            print(
                f"    NebulaBayes: O/H={result['o_h']:.5f} "
                f"[{result['o_h_ci68_low']:.5f}, {result['o_h_ci68_high']:.5f}]; "
                f"logU={result['log_u']:.5f} "
                f"[{result['log_u_ci68_low']:.5f}, {result['log_u_ci68_high']:.5f}]; "
                f"log(P/k)={result['log_pk']:.5f} "
                f"[{result['log_pk_ci68_low']:.5f}, {result['log_pk_ci68_high']:.5f}]; "
                f"chi2_red={result['chi2_red']:.4g}; "
                f"nlocalmax={result['n_localmax']}; flag={result['flag']}; "
                f"valid={result['valid']}"
            )
        if ENABLE_JY22_METALLICITY_PRODUCTS:
            result = summary["jy22"]
            print(
                f"    JY22: O/H={result['o_h']:.5f} "
                f"[{result['o_h_16']:.5f}, {result['o_h_84']:.5f}]; "
                f"logU={result['log_u']:.5f} "
                f"[{result['log_u_16']:.5f}, {result['log_u_84']:.5f}]; "
                f"chi2_min={result['chi2_min']:.4g}; "
                f"E(B-V)={summary['jy22_ebv']:.5f}; "
                f"flag={result['flag']}; valid={result['valid']} "
                "(posterior means; marginal equal-tailed p16/p84)"
            )
    print(
        "Model covariance note: JY22 uses the full first-order N2/S2/R3 "
        "covariance from raw line errors and the shared Balmer correction; "
        "pyqz/NebulaBayes independent-error APIs do not represent the analogous "
        "shared covariance."
    )
    print("--------------------------------------------------------------")

# ------------------------------------------------------------------
# 10.  PyNeb electron density and Te-based HII metallicity products
# ------------------------------------------------------------------

if ENABLE_TE_METALLICITY_PRODUCTS:
    PYNEB_S2 = pn.Atom("S", 2)
    PYNEB_N2 = pn.Atom("N", 2)
    PYNEB_S3 = pn.Atom("S", 3)
    PYNEB_O2 = pn.Atom("O", 2)
    PYNEB_O3 = pn.Atom("O", 3)

    TE_DETECTED = {
        "HB4861": line_detection_mask(HB4861_FLUX, HB4861_FLUX_ERR),
        "HA6562": line_detection_mask(HA6562_FLUX, HA6562_FLUX_ERR),
        "OIII5006": line_detection_mask(OIII5006_FLUX, OIII5006_FLUX_ERR),
        "NII6583": line_detection_mask(NII6583_FLUX, NII6583_FLUX_ERR),
        "SII6716": line_detection_mask(SII6716_FLUX, SII6716_FLUX_ERR),
        "SII6730": line_detection_mask(SII6730_FLUX, SII6730_FLUX_ERR),
    }
    for line_base in TE_LINE_BASES:
        TE_DETECTED[line_base] = line_detection_mask(
            globals()[f"{line_base}_FLUX"],
            globals()[f"{line_base}_FLUX_ERR"],
        )
    TE_DETECTED["NII6548_6583"] = line_detection_mask(
        NII6548_FLUX + NII6583_FLUX,
        quadrature_sum(NII6548_FLUX_ERR, NII6583_FLUX_ERR),
    )
    TE_DETECTED["OIII4959_5007"] = line_detection_mask(
        OIII4958_FLUX + OIII5006_FLUX,
        quadrature_sum(OIII4958_FLUX_ERR, OIII5006_FLUX_ERR),
    )
    TE_DETECTED["OII7325"] = line_detection_mask(
        OII7318_FLUX + OII7319_FLUX + OII7329_FLUX + OII7330_FLUX,
        quadrature_sum(
            OII7318_FLUX_ERR,
            OII7319_FLUX_ERR,
            OII7329_FLUX_ERR,
            OII7330_FLUX_ERR,
        ),
    )

    SII_DENSITY_VALID = TE_DETECTED["SII6716"] & TE_DETECTED["SII6730"]
    SII_DENSITY_RATIO = positive_ratio(SII6716_FLUX_corr, SII6730_FLUX_corr)
    SII_DENSITY_RATIO_ERR = ratio_linear_error(
        SII6716_FLUX_corr,
        SII6730_FLUX_corr,
        SII6716_FLUX_ERR_corr,
        SII6730_FLUX_ERR_corr,
    )

    NE_SII_RAW, NE_SII_EXACT_HIGH_DENSITY = sii_density_from_ratio_lookup(
        PYNEB_S2,
        SII_DENSITY_RATIO,
        SII_DENSITY_VALID,
        tem=10000.0,
        exact_high_density=True,
        return_exact_mask=True,
    )
    sii_low_density_ratio = PYNEB_S2.getLowDensRatio(wave1=6716, wave2=6731)
    NE_SII_FIXED20 = (
        SII_DENSITY_VALID &
        (
            (SII_DENSITY_RATIO >= sii_low_density_ratio) |
            (np.isfinite(NE_SII_RAW) & (NE_SII_RAW < 20.0))
        )
    )
    NE_SII = np.where(NE_SII_FIXED20, 20.0, NE_SII_RAW)
    NE_SII = np.where(
        SII_DENSITY_VALID & np.isfinite(NE_SII) & (NE_SII >= 20.0),
        NE_SII,
        np.nan,
    )
    NE_SII_ALL = np.where(
        np.isfinite(V_STARS2),
        np.where(np.isfinite(NE_SII), NE_SII, 20.0),
        np.nan,
    )

    ne_ratio_low = np.maximum(SII_DENSITY_RATIO - SII_DENSITY_RATIO_ERR, 1e-30)
    ne_ratio_high = SII_DENSITY_RATIO + SII_DENSITY_RATIO_ERR
    NE_SII_LOW = sii_density_from_ratio_lookup(
        PYNEB_S2,
        ne_ratio_low,
        SII_DENSITY_VALID,
        tem=10000.0,
        exact_high_density=True,
    )
    NE_SII_HIGH = sii_density_from_ratio_lookup(
        PYNEB_S2,
        ne_ratio_high,
        SII_DENSITY_VALID,
        tem=10000.0,
        exact_high_density=True,
    )
    NE_SII_LOW = np.where(np.isfinite(NE_SII_LOW) & (NE_SII_LOW >= 20.0), NE_SII_LOW, 20.0)
    NE_SII_HIGH = np.where(np.isfinite(NE_SII_HIGH) & (NE_SII_HIGH >= 20.0), NE_SII_HIGH, 20.0)
    NE_SII_ERR = finite_difference_error(NE_SII_LOW, NE_SII_HIGH, SII_DENSITY_VALID)
    NE_SII_ERR = np.where(NE_SII_FIXED20, 0.0, NE_SII_ERR)

    NII_TE_DENOM = NII6548_FLUX_corr + NII6583_FLUX_corr
    NII_TE_DENOM_ERR = quadrature_sum(NII6548_FLUX_ERR_corr, NII6583_FLUX_ERR_corr)
    NII_TE_RATIO = positive_ratio(NII5754_FLUX_corr, NII_TE_DENOM)
    NII_TE_RATIO_ERR = ratio_linear_error(
        NII5754_FLUX_corr,
        NII_TE_DENOM,
        NII5754_FLUX_ERR_corr,
        NII_TE_DENOM_ERR,
    )
    SIII_TE_DENOM = 3.5 * SIII9068_FLUX_corr
    SIII_TE_DENOM_ERR = 3.5 * SIII9068_FLUX_ERR_corr
    SIII_TE_RATIO = positive_ratio(SIII6312_FLUX_corr, SIII_TE_DENOM)
    SIII_TE_RATIO_ERR = ratio_linear_error(
        SIII6312_FLUX_corr,
        SIII_TE_DENOM,
        SIII6312_FLUX_ERR_corr,
        SIII_TE_DENOM_ERR,
    )
    OII7325_FLUX_CORR = (
        OII7318_FLUX_corr + OII7319_FLUX_corr +
        OII7329_FLUX_corr + OII7330_FLUX_corr
    )
    OII7325_FLUX_ERR_CORR = quadrature_sum(
        OII7318_FLUX_ERR_corr,
        OII7319_FLUX_ERR_corr,
        OII7329_FLUX_ERR_corr,
        OII7330_FLUX_ERR_corr,
    )
    OIII4959_5007_FLUX_CORR = OIII4958_FLUX_corr + OIII5006_FLUX_corr
    OIII4959_5007_FLUX_ERR_CORR = quadrature_sum(
        OIII4958_FLUX_ERR_corr,
        OIII5006_FLUX_ERR_corr,
    )
    OII7325_HBETA_RATIO = 100.0 * positive_ratio(OII7325_FLUX_CORR, HB4861_FLUX_corr)
    OII7325_HBETA_RATIO_ERR = 100.0 * ratio_linear_error(
        OII7325_FLUX_CORR,
        HB4861_FLUX_corr,
        OII7325_FLUX_ERR_CORR,
        HB4861_FLUX_ERR_corr,
    )
    OIII_HBETA_RATIO = 100.0 * positive_ratio(OIII4959_5007_FLUX_CORR, HB4861_FLUX_corr)
    OIII_HBETA_RATIO_ERR = 100.0 * ratio_linear_error(
        OIII4959_5007_FLUX_CORR,
        HB4861_FLUX_corr,
        OIII4959_5007_FLUX_ERR_CORR,
        HB4861_FLUX_ERR_corr,
    )

    def compute_te_region_products(region_mask):
        common_te = (
            region_mask &
            TE_DETECTED["HB4861"] &
            TE_DETECTED["HA6562"] &
            SII_DENSITY_VALID &
            np.isfinite(NE_SII)
        )
        nii_te_valid = (
            common_te &
            TE_DETECTED["NII5754"] &
            TE_DETECTED["NII6548_6583"]
        )
        siii_te_valid = (
            common_te &
            TE_DETECTED["SIII6312"] &
            TE_DETECTED["SIII9068"]
        )
        oii7325_valid = TE_DETECTED["OII7325"]
        oiii_valid = TE_DETECTED["OIII4959_5007"]
        br24_direct_valid = (
            nii_te_valid & siii_te_valid & oii7325_valid & oiii_valid
        )
        nii_oii7325_valid = nii_te_valid & oii7325_valid & oiii_valid
        nii_k25_valid = nii_te_valid

        te_nii = pyneb_temden_map(
            PYNEB_N2,
            NII_TE_RATIO,
            nii_te_valid,
            den=NE_SII,
            to_eval="L(5755)/(L(6548)+L(6584))",
        )
        te_nii = sanitize_temperature_map(te_nii, nii_te_valid)
        te_nii_ratio_low = pyneb_temden_map(
            PYNEB_N2,
            np.maximum(NII_TE_RATIO - NII_TE_RATIO_ERR, 1e-30),
            nii_te_valid,
            den=NE_SII,
            to_eval="L(5755)/(L(6548)+L(6584))",
        )
        te_nii_ratio_high = pyneb_temden_map(
            PYNEB_N2,
            NII_TE_RATIO + NII_TE_RATIO_ERR,
            nii_te_valid,
            den=NE_SII,
            to_eval="L(5755)/(L(6548)+L(6584))",
        )
        te_nii_den_low = pyneb_temden_map(
            PYNEB_N2,
            NII_TE_RATIO,
            nii_te_valid,
            den=np.maximum(NE_SII - NE_SII_ERR, 20.0),
            to_eval="L(5755)/(L(6548)+L(6584))",
        )
        te_nii_den_high = pyneb_temden_map(
            PYNEB_N2,
            NII_TE_RATIO,
            nii_te_valid,
            den=NE_SII + NE_SII_ERR,
            to_eval="L(5755)/(L(6548)+L(6584))",
        )
        te_nii_err = quadrature_sum(
            finite_difference_error(
                te_nii_ratio_low, te_nii_ratio_high, np.isfinite(te_nii)
            ),
            finite_difference_error(
                te_nii_den_low, te_nii_den_high, np.isfinite(te_nii)
            ),
        )
        nii_k25_valid = (
            nii_k25_valid &
            np.isfinite(te_nii) &
            (te_nii >= 8000.0) &
            (te_nii <= 13000.0)
        )

        te_siii = pyneb_temden_map(
            PYNEB_S3,
            SIII_TE_RATIO,
            siii_te_valid,
            den=NE_SII,
            to_eval="L(6312)/(L(9069)+2.5*L(9069))",
        )
        te_siii = sanitize_temperature_map(te_siii, siii_te_valid)
        te_siii_ratio_low = pyneb_temden_map(
            PYNEB_S3,
            np.maximum(SIII_TE_RATIO - SIII_TE_RATIO_ERR, 1e-30),
            siii_te_valid,
            den=NE_SII,
            to_eval="L(6312)/(L(9069)+2.5*L(9069))",
        )
        te_siii_ratio_high = pyneb_temden_map(
            PYNEB_S3,
            SIII_TE_RATIO + SIII_TE_RATIO_ERR,
            siii_te_valid,
            den=NE_SII,
            to_eval="L(6312)/(L(9069)+2.5*L(9069))",
        )
        te_siii_den_low = pyneb_temden_map(
            PYNEB_S3,
            SIII_TE_RATIO,
            siii_te_valid,
            den=np.maximum(NE_SII - NE_SII_ERR, 20.0),
            to_eval="L(6312)/(L(9069)+2.5*L(9069))",
        )
        te_siii_den_high = pyneb_temden_map(
            PYNEB_S3,
            SIII_TE_RATIO,
            siii_te_valid,
            den=NE_SII + NE_SII_ERR,
            to_eval="L(6312)/(L(9069)+2.5*L(9069))",
        )
        te_siii_err = quadrature_sum(
            finite_difference_error(
                te_siii_ratio_low, te_siii_ratio_high, np.isfinite(te_siii)
            ),
            finite_difference_error(
                te_siii_den_low, te_siii_den_high, np.isfinite(te_siii)
            ),
        )

        te_siii_nii_chain = sanitize_temperature_map(
            brazzini_te_siii_from_nii(te_nii),
            nii_te_valid & np.isfinite(te_nii),
        )
        te_siii_nii_chain_err = np.where(
            np.isfinite(te_siii_nii_chain),
            brazzini_te_siii_from_nii_error(te_nii, te_nii_err),
            np.nan,
        )
        te_oiii_br24 = sanitize_temperature_map(
            brazzini_te_oiii_from_siii(te_siii),
            br24_direct_valid & np.isfinite(te_siii),
        )
        te_oiii_br24_err = np.where(
            np.isfinite(te_oiii_br24),
            brazzini_te_oiii_from_siii_error(te_siii, te_siii_err),
            np.nan,
        )
        te_oiii_nii_chain = sanitize_temperature_map(
            brazzini_te_oiii_from_siii(te_siii_nii_chain),
            nii_oii7325_valid & np.isfinite(te_siii_nii_chain),
        )
        te_oiii_nii_chain_err = np.where(
            np.isfinite(te_oiii_nii_chain),
            brazzini_te_oiii_from_siii_error(
                te_siii_nii_chain, te_siii_nii_chain_err
            ),
            np.nan,
        )

        o_plus_br24, o_plus_br24_err = pyneb_ion_abundance_error_map(
            PYNEB_O2,
            OII7325_HBETA_RATIO,
            OII7325_HBETA_RATIO_ERR,
            te_nii,
            te_nii_err,
            NE_SII,
            NE_SII_ERR,
            br24_direct_valid & np.isfinite(te_nii),
            to_eval="L(7318)+L(7319)+L(7329)+L(7330)",
        )
        o_doubleplus_br24, o_doubleplus_br24_err = pyneb_ion_abundance_error_map(
            PYNEB_O3,
            OIII_HBETA_RATIO,
            OIII_HBETA_RATIO_ERR,
            te_oiii_br24,
            te_oiii_br24_err,
            NE_SII,
            NE_SII_ERR,
            br24_direct_valid & np.isfinite(te_oiii_br24),
            to_eval="L(4959)+L(5007)",
        )
        o_h_br24_direct, o_h_br24_direct_err = oxygen_abundance_log_map(
            o_plus_br24,
            o_plus_br24_err,
            o_doubleplus_br24,
            o_doubleplus_br24_err,
            br24_direct_valid,
        )

        o_plus_nii_oii7325, o_plus_nii_oii7325_err = pyneb_ion_abundance_error_map(
            PYNEB_O2,
            OII7325_HBETA_RATIO,
            OII7325_HBETA_RATIO_ERR,
            te_nii,
            te_nii_err,
            NE_SII,
            NE_SII_ERR,
            nii_oii7325_valid & np.isfinite(te_nii),
            to_eval="L(7318)+L(7319)+L(7329)+L(7330)",
        )
        (
            o_doubleplus_nii_oii7325,
            o_doubleplus_nii_oii7325_err,
        ) = pyneb_ion_abundance_error_map(
            PYNEB_O3,
            OIII_HBETA_RATIO,
            OIII_HBETA_RATIO_ERR,
            te_oiii_nii_chain,
            te_oiii_nii_chain_err,
            NE_SII,
            NE_SII_ERR,
            nii_oii7325_valid & np.isfinite(te_oiii_nii_chain),
            to_eval="L(4959)+L(5007)",
        )
        o_h_nii_oii7325, o_h_nii_oii7325_err = oxygen_abundance_log_map(
            o_plus_nii_oii7325,
            o_plus_nii_oii7325_err,
            o_doubleplus_nii_oii7325,
            o_doubleplus_nii_oii7325_err,
            nii_oii7325_valid,
        )

        o_h_nii_k25, o_h_nii_k25_err = kreckel_mendez_nii_oh_map(
            te_nii,
            te_nii_err,
            nii_k25_valid,
        )

        return {
            "COMMON_TE": common_te,
            "NII_TE_VALID": nii_te_valid,
            "SIII_TE_VALID": siii_te_valid,
            "BR24_DIRECT_VALID": br24_direct_valid,
            "NII_OII7325_VALID": nii_oii7325_valid,
            "NII_K25_VALID": nii_k25_valid,
            "TE_NII": te_nii,
            "TE_NII_ERR": te_nii_err,
            "TE_SIII_BR24": te_siii,
            "TE_SIII_BR24_ERR": te_siii_err,
            "TE_SIII_NII_CHAIN": te_siii_nii_chain,
            "TE_SIII_NII_CHAIN_ERR": te_siii_nii_chain_err,
            "TE_OIII_BR24": te_oiii_br24,
            "TE_OIII_BR24_ERR": te_oiii_br24_err,
            "TE_OIII_NII_CHAIN": te_oiii_nii_chain,
            "TE_OIII_NII_CHAIN_ERR": te_oiii_nii_chain_err,
            "O_PLUS_BR24": o_plus_br24,
            "O_PLUS_BR24_ERR": o_plus_br24_err,
            "O_DOUBLEPLUS_BR24": o_doubleplus_br24,
            "O_DOUBLEPLUS_BR24_ERR": o_doubleplus_br24_err,
            "O_H_BR24_DIRECT": o_h_br24_direct,
            "O_H_BR24_DIRECT_ERR": o_h_br24_direct_err,
            "O_PLUS_NII_OII7325": o_plus_nii_oii7325,
            "O_PLUS_NII_OII7325_ERR": o_plus_nii_oii7325_err,
            "O_DOUBLEPLUS_NII_OII7325": o_doubleplus_nii_oii7325,
            "O_DOUBLEPLUS_NII_OII7325_ERR": o_doubleplus_nii_oii7325_err,
            "O_H_NII_OII7325": o_h_nii_oii7325,
            "O_H_NII_OII7325_ERR": o_h_nii_oii7325_err,
            "O_H_NII_K25": o_h_nii_k25,
            "O_H_NII_K25_ERR": o_h_nii_k25_err,
        }

    te_hii_products = compute_te_region_products(mask_HII)
    te_sf_products = compute_te_region_products(mask_SF)

    NII_TE_VALID = te_hii_products["NII_TE_VALID"]
    SIII_TE_VALID = te_hii_products["SIII_TE_VALID"]
    BR24_DIRECT_VALID = te_hii_products["BR24_DIRECT_VALID"]
    NII_OII7325_VALID = te_hii_products["NII_OII7325_VALID"]
    NII_K25_VALID = te_hii_products["NII_K25_VALID"]
    TE_NII_HII = te_hii_products["TE_NII"]
    TE_NII_HII_ERR = te_hii_products["TE_NII_ERR"]
    TE_SIII_BR24_HII = te_hii_products["TE_SIII_BR24"]
    TE_SIII_BR24_HII_ERR = te_hii_products["TE_SIII_BR24_ERR"]
    TE_SIII_NII_CHAIN_HII = te_hii_products["TE_SIII_NII_CHAIN"]
    TE_SIII_NII_CHAIN_HII_ERR = te_hii_products["TE_SIII_NII_CHAIN_ERR"]
    TE_OIII_BR24_HII = te_hii_products["TE_OIII_BR24"]
    TE_OIII_BR24_HII_ERR = te_hii_products["TE_OIII_BR24_ERR"]
    TE_OIII_NII_CHAIN_HII = te_hii_products["TE_OIII_NII_CHAIN"]
    TE_OIII_NII_CHAIN_HII_ERR = te_hii_products["TE_OIII_NII_CHAIN_ERR"]
    O_PLUS_BR24_HII = te_hii_products["O_PLUS_BR24"]
    O_PLUS_BR24_HII_ERR = te_hii_products["O_PLUS_BR24_ERR"]
    O_DOUBLEPLUS_BR24_HII = te_hii_products["O_DOUBLEPLUS_BR24"]
    O_DOUBLEPLUS_BR24_HII_ERR = te_hii_products["O_DOUBLEPLUS_BR24_ERR"]
    O_H_BR24_DIRECT_HII = te_hii_products["O_H_BR24_DIRECT"]
    O_H_BR24_DIRECT_HII_ERR = te_hii_products["O_H_BR24_DIRECT_ERR"]
    O_PLUS_NII_OII7325_HII = te_hii_products["O_PLUS_NII_OII7325"]
    O_PLUS_NII_OII7325_HII_ERR = te_hii_products["O_PLUS_NII_OII7325_ERR"]
    O_DOUBLEPLUS_NII_OII7325_HII = te_hii_products["O_DOUBLEPLUS_NII_OII7325"]
    O_DOUBLEPLUS_NII_OII7325_HII_ERR = te_hii_products["O_DOUBLEPLUS_NII_OII7325_ERR"]
    O_H_NII_OII7325_HII = te_hii_products["O_H_NII_OII7325"]
    O_H_NII_OII7325_HII_ERR = te_hii_products["O_H_NII_OII7325_ERR"]
    O_H_NII_K25_HII = te_hii_products["O_H_NII_K25"]
    O_H_NII_K25_HII_ERR = te_hii_products["O_H_NII_K25_ERR"]

    BR24_DIRECT_VALID_SF = te_sf_products["BR24_DIRECT_VALID"]
    NII_OII7325_VALID_SF = te_sf_products["NII_OII7325_VALID"]
    NII_K25_VALID_SF = te_sf_products["NII_K25_VALID"]
    BR24_DIRECT_VALID_HII = BR24_DIRECT_VALID
    NII_OII7325_VALID_HII = NII_OII7325_VALID
    NII_K25_VALID_HII = NII_K25_VALID
    TE_NII_SF = te_sf_products["TE_NII"]
    TE_NII_SF_ERR = te_sf_products["TE_NII_ERR"]
    TE_SIII_BR24_SF = te_sf_products["TE_SIII_BR24"]
    TE_SIII_BR24_SF_ERR = te_sf_products["TE_SIII_BR24_ERR"]
    TE_SIII_NII_CHAIN_SF = te_sf_products["TE_SIII_NII_CHAIN"]
    TE_SIII_NII_CHAIN_SF_ERR = te_sf_products["TE_SIII_NII_CHAIN_ERR"]
    TE_OIII_BR24_SF = te_sf_products["TE_OIII_BR24"]
    TE_OIII_BR24_SF_ERR = te_sf_products["TE_OIII_BR24_ERR"]
    TE_OIII_NII_CHAIN_SF = te_sf_products["TE_OIII_NII_CHAIN"]
    TE_OIII_NII_CHAIN_SF_ERR = te_sf_products["TE_OIII_NII_CHAIN_ERR"]
    O_PLUS_BR24_SF = te_sf_products["O_PLUS_BR24"]
    O_PLUS_BR24_SF_ERR = te_sf_products["O_PLUS_BR24_ERR"]
    O_DOUBLEPLUS_BR24_SF = te_sf_products["O_DOUBLEPLUS_BR24"]
    O_DOUBLEPLUS_BR24_SF_ERR = te_sf_products["O_DOUBLEPLUS_BR24_ERR"]
    O_H_BR24_DIRECT_SF = te_sf_products["O_H_BR24_DIRECT"]
    O_H_BR24_DIRECT_SF_ERR = te_sf_products["O_H_BR24_DIRECT_ERR"]
    O_PLUS_NII_OII7325_SF = te_sf_products["O_PLUS_NII_OII7325"]
    O_PLUS_NII_OII7325_SF_ERR = te_sf_products["O_PLUS_NII_OII7325_ERR"]
    O_DOUBLEPLUS_NII_OII7325_SF = te_sf_products["O_DOUBLEPLUS_NII_OII7325"]
    O_DOUBLEPLUS_NII_OII7325_SF_ERR = te_sf_products["O_DOUBLEPLUS_NII_OII7325_ERR"]
    O_H_NII_OII7325_SF = te_sf_products["O_H_NII_OII7325"]
    O_H_NII_OII7325_SF_ERR = te_sf_products["O_H_NII_OII7325_ERR"]
    O_H_NII_K25_SF = te_sf_products["O_H_NII_K25"]
    O_H_NII_K25_SF_ERR = te_sf_products["O_H_NII_K25_ERR"]

    INTEGRATED_TE_LINE_K = {
        "HB4861": k_HB4861,
        "HA6562": k_HA6562,
        "OIII5006": k_OIII5006,
        "NII6583": k_NII6583,
        "SII6716": k_SII6716,
        "SII6730": k_SII6730,
    }
    INTEGRATED_TE_LINE_K.update(TE_LINE_K)

    def finite_scalar(value):
        value = float(np.asarray(value, dtype=float).ravel()[0])
        return value if np.isfinite(value) else np.nan

    def finite_positive(value):
        return np.isfinite(value) and value > 0

    def sum_raw_line(line_base, valid_mask):
        flux = np.asarray(globals()[f"{line_base}_FLUX"], dtype=float)
        flux_err = np.asarray(globals()[f"{line_base}_FLUX_ERR"], dtype=float)
        selected = np.asarray(valid_mask, dtype=bool) & np.isfinite(flux)
        if not np.any(selected):
            return np.nan, np.nan
        return (
            float(np.nansum(np.where(selected, flux, np.nan))),
            float(integrated_flux_error(flux_err, selected)),
        )

    def integrated_ebv_for_mask(valid_mask):
        hb, hb_err = sum_raw_line("HB4861", valid_mask)
        ha, ha_err = sum_raw_line("HA6562", valid_mask)
        if not (finite_positive(hb) and finite_positive(ha)):
            return np.nan, np.nan
        bd = ha / hb
        if bd < R_int:
            bd = R_int
        ebv = convert_bd_to_ebv(bd, k_HB4861, k_HA6562, R_int)
        ebv_err = convert_bd_to_ebv_error(
            ha, hb, ha_err, hb_err, k_HB4861, k_HA6562
        )
        return finite_scalar(ebv), finite_scalar(ebv_err)

    def corrected_integrated_lines(line_bases, valid_mask):
        ebv, ebv_err = integrated_ebv_for_mask(valid_mask)
        corrected = {}
        for line_base in line_bases:
            raw_flux, raw_err = sum_raw_line(line_base, valid_mask)
            if np.isfinite(ebv) and np.isfinite(raw_flux):
                k_line = INTEGRATED_TE_LINE_K[line_base]
                flux_corr = correct_flux_with_ebv(raw_flux, ebv, k_line)
                err_corr = correct_flux_error_with_ebv(
                    raw_flux, raw_err, ebv, k_line, ebv_err
                )
                corrected[line_base] = (
                    finite_scalar(flux_corr),
                    finite_scalar(err_corr),
                )
            else:
                corrected[line_base] = (np.nan, np.nan)
        return corrected, ebv, ebv_err

    def line_value(lines, line_base):
        return lines.get(line_base, (np.nan, np.nan))

    def compute_integrated_density(valid_mask):
        n_pix = int(np.count_nonzero(valid_mask))
        s6716, s6716_err = sum_raw_line("SII6716", valid_mask)
        s6730, s6730_err = sum_raw_line("SII6730", valid_mask)
        ratio = s6716 / s6730 if finite_positive(s6716) and finite_positive(s6730) else np.nan
        ratio_err = finite_scalar(
            ratio_linear_error(
                as_single_pixel(s6716),
                as_single_pixel(s6730),
                as_single_pixel(s6716_err),
                as_single_pixel(s6730_err),
            )
        )
        valid = np.array([[finite_positive(ratio)]])
        ne_raw_map, exact_mask = sii_density_from_ratio_lookup(
            PYNEB_S2,
            as_single_pixel(ratio),
            valid,
            tem=10000.0,
            exact_high_density=True,
            return_exact_mask=True,
        )
        ne_raw = finite_scalar(ne_raw_map)
        fixed20 = bool(
            valid[0, 0] and (
                ratio >= sii_low_density_ratio or
                (np.isfinite(ne_raw) and ne_raw < 20.0)
            )
        )
        ne = 20.0 if fixed20 else ne_raw
        ne = ne if np.isfinite(ne) and ne >= 20.0 else np.nan

        ratio_low = max(ratio - ratio_err, 1e-30) if np.isfinite(ratio_err) else np.nan
        ratio_high = ratio + ratio_err if np.isfinite(ratio_err) else np.nan
        ne_low = sii_density_from_ratio_lookup(
            PYNEB_S2,
            as_single_pixel(ratio_low),
            np.array([[finite_positive(ratio_low)]]),
            tem=10000.0,
            exact_high_density=True,
        )
        ne_high = sii_density_from_ratio_lookup(
            PYNEB_S2,
            as_single_pixel(ratio_high),
            np.array([[finite_positive(ratio_high)]]),
            tem=10000.0,
            exact_high_density=True,
        )
        ne_low = max(finite_scalar(ne_low), 20.0) if np.isfinite(finite_scalar(ne_low)) else 20.0
        ne_high = max(finite_scalar(ne_high), 20.0) if np.isfinite(finite_scalar(ne_high)) else 20.0
        ne_err = 0.0 if fixed20 else 0.5 * abs(ne_high - ne_low)
        ne_err = ne_err if np.isfinite(ne) and np.isfinite(ne_err) else np.nan
        return {
            "n_pix": n_pix,
            "value": ne,
            "err": ne_err,
            "fixed20": fixed20,
            "exact_high_density": bool(np.any(exact_mask)),
        }

    def compute_integrated_nii_te(valid_mask):
        line_bases = ["HB4861", "HA6562", "NII5754", "NII6548", "NII6583"]
        lines, ebv, ebv_err = corrected_integrated_lines(line_bases, valid_mask)
        density = compute_integrated_density(valid_mask)
        n5755, n5755_err = line_value(lines, "NII5754")
        n6548, n6548_err = line_value(lines, "NII6548")
        n6583, n6583_err = line_value(lines, "NII6583")
        denom = n6548 + n6583
        denom_err = np.sqrt(n6548_err**2 + n6583_err**2)
        ratio = positive_ratio(as_single_pixel(n5755), as_single_pixel(denom))
        ratio_err = ratio_linear_error(
            as_single_pixel(n5755),
            as_single_pixel(denom),
            as_single_pixel(n5755_err),
            as_single_pixel(denom_err),
        )
        valid = np.array([[finite_positive(single_pixel_value(ratio))]])
        ne = as_single_pixel(density["value"])
        ne_err = as_single_pixel(density["err"])
        te = pyneb_temden_map(
            PYNEB_N2,
            ratio,
            valid,
            den=ne,
            to_eval="L(5755)/(L(6548)+L(6584))",
        )
        te = sanitize_temperature_map(te, valid)
        te_ratio_low = pyneb_temden_map(
            PYNEB_N2,
            np.maximum(ratio - ratio_err, 1e-30),
            valid,
            den=ne,
            to_eval="L(5755)/(L(6548)+L(6584))",
        )
        te_ratio_high = pyneb_temden_map(
            PYNEB_N2,
            ratio + ratio_err,
            valid,
            den=ne,
            to_eval="L(5755)/(L(6548)+L(6584))",
        )
        te_den_low = pyneb_temden_map(
            PYNEB_N2,
            ratio,
            valid,
            den=np.maximum(ne - ne_err, 20.0),
            to_eval="L(5755)/(L(6548)+L(6584))",
        )
        te_den_high = pyneb_temden_map(
            PYNEB_N2,
            ratio,
            valid,
            den=ne + ne_err,
            to_eval="L(5755)/(L(6548)+L(6584))",
        )
        te_err = quadrature_sum(
            finite_difference_error(te_ratio_low, te_ratio_high, np.isfinite(te)),
            finite_difference_error(te_den_low, te_den_high, np.isfinite(te)),
        )
        return {
            "n_pix": int(np.count_nonzero(valid_mask)),
            "ebv": ebv,
            "ebv_err": ebv_err,
            "density": density,
            "lines": lines,
            "value": finite_scalar(te),
            "err": finite_scalar(te_err),
        }

    def compute_integrated_siii_te(valid_mask):
        line_bases = ["HB4861", "HA6562", "SIII6312", "SIII9068"]
        lines, ebv, ebv_err = corrected_integrated_lines(line_bases, valid_mask)
        density = compute_integrated_density(valid_mask)
        s6312, s6312_err = line_value(lines, "SIII6312")
        s9068, s9068_err = line_value(lines, "SIII9068")
        denom = 3.5 * s9068
        denom_err = 3.5 * s9068_err
        ratio = positive_ratio(as_single_pixel(s6312), as_single_pixel(denom))
        ratio_err = ratio_linear_error(
            as_single_pixel(s6312),
            as_single_pixel(denom),
            as_single_pixel(s6312_err),
            as_single_pixel(denom_err),
        )
        valid = np.array([[finite_positive(single_pixel_value(ratio))]])
        ne = as_single_pixel(density["value"])
        ne_err = as_single_pixel(density["err"])
        te = pyneb_temden_map(
            PYNEB_S3,
            ratio,
            valid,
            den=ne,
            to_eval="L(6312)/(L(9069)+2.5*L(9069))",
        )
        te = sanitize_temperature_map(te, valid)
        te_ratio_low = pyneb_temden_map(
            PYNEB_S3,
            np.maximum(ratio - ratio_err, 1e-30),
            valid,
            den=ne,
            to_eval="L(6312)/(L(9069)+2.5*L(9069))",
        )
        te_ratio_high = pyneb_temden_map(
            PYNEB_S3,
            ratio + ratio_err,
            valid,
            den=ne,
            to_eval="L(6312)/(L(9069)+2.5*L(9069))",
        )
        te_den_low = pyneb_temden_map(
            PYNEB_S3,
            ratio,
            valid,
            den=np.maximum(ne - ne_err, 20.0),
            to_eval="L(6312)/(L(9069)+2.5*L(9069))",
        )
        te_den_high = pyneb_temden_map(
            PYNEB_S3,
            ratio,
            valid,
            den=ne + ne_err,
            to_eval="L(6312)/(L(9069)+2.5*L(9069))",
        )
        te_err = quadrature_sum(
            finite_difference_error(te_ratio_low, te_ratio_high, np.isfinite(te)),
            finite_difference_error(te_den_low, te_den_high, np.isfinite(te)),
        )
        return {
            "n_pix": int(np.count_nonzero(valid_mask)),
            "ebv": ebv,
            "ebv_err": ebv_err,
            "density": density,
            "lines": lines,
            "value": finite_scalar(te),
            "err": finite_scalar(te_err),
        }

    def compute_integrated_oxygen(valid_mask, *, use_measured_siii):
        line_bases = [
            "HB4861", "HA6562", "NII5754", "NII6548", "NII6583",
            "OIII4958", "OIII5006", "OII7318", "OII7319", "OII7329",
            "OII7330",
        ]
        if use_measured_siii:
            line_bases += ["SIII6312", "SIII9068"]
        lines, ebv, ebv_err = corrected_integrated_lines(line_bases, valid_mask)
        density = compute_integrated_density(valid_mask)
        nii = compute_integrated_nii_te(valid_mask)
        te_nii = as_single_pixel(nii["value"])
        te_nii_err = as_single_pixel(nii["err"])
        if use_measured_siii:
            siii = compute_integrated_siii_te(valid_mask)
            te_oiii_value = finite_scalar(
                sanitize_temperature_map(
                    as_single_pixel(brazzini_te_oiii_from_siii(siii["value"])),
                    np.array([[np.isfinite(siii["value"])]]),
                )
            )
            te_oiii_err_value = finite_scalar(
                brazzini_te_oiii_from_siii_error(siii["value"], siii["err"])
            )
        else:
            te_siii_chain = finite_scalar(
                sanitize_temperature_map(
                    as_single_pixel(brazzini_te_siii_from_nii(nii["value"])),
                    np.array([[np.isfinite(nii["value"])]]),
                )
            )
            te_siii_chain_err = finite_scalar(
                brazzini_te_siii_from_nii_error(nii["value"], nii["err"])
            )
            te_oiii_value = finite_scalar(
                sanitize_temperature_map(
                    as_single_pixel(brazzini_te_oiii_from_siii(te_siii_chain)),
                    np.array([[np.isfinite(te_siii_chain)]]),
                )
            )
            te_oiii_err_value = finite_scalar(
                brazzini_te_oiii_from_siii_error(te_siii_chain, te_siii_chain_err)
            )

        hb, hb_err = line_value(lines, "HB4861")
        oii_flux = sum(line_value(lines, line)[0] for line in ["OII7318", "OII7319", "OII7329", "OII7330"])
        oii_err = np.sqrt(
            sum(line_value(lines, line)[1] ** 2 for line in ["OII7318", "OII7319", "OII7329", "OII7330"])
        )
        oiii_flux = line_value(lines, "OIII4958")[0] + line_value(lines, "OIII5006")[0]
        oiii_err = np.sqrt(line_value(lines, "OIII4958")[1] ** 2 + line_value(lines, "OIII5006")[1] ** 2)
        oii_ratio = 100.0 * positive_ratio(as_single_pixel(oii_flux), as_single_pixel(hb))
        oii_ratio_err = 100.0 * ratio_linear_error(
            as_single_pixel(oii_flux),
            as_single_pixel(hb),
            as_single_pixel(oii_err),
            as_single_pixel(hb_err),
        )
        oiii_ratio = 100.0 * positive_ratio(as_single_pixel(oiii_flux), as_single_pixel(hb))
        oiii_ratio_err = 100.0 * ratio_linear_error(
            as_single_pixel(oiii_flux),
            as_single_pixel(hb),
            as_single_pixel(oiii_err),
            as_single_pixel(hb_err),
        )
        valid = np.array([[np.isfinite(nii["value"]) and np.isfinite(te_oiii_value)]])
        ne = as_single_pixel(density["value"])
        ne_err = as_single_pixel(density["err"])
        o_plus, o_plus_err = pyneb_ion_abundance_error_map(
            PYNEB_O2,
            oii_ratio,
            oii_ratio_err,
            te_nii,
            te_nii_err,
            ne,
            ne_err,
            valid,
            to_eval="L(7318)+L(7319)+L(7329)+L(7330)",
        )
        o_doubleplus, o_doubleplus_err = pyneb_ion_abundance_error_map(
            PYNEB_O3,
            oiii_ratio,
            oiii_ratio_err,
            as_single_pixel(te_oiii_value),
            as_single_pixel(te_oiii_err_value),
            ne,
            ne_err,
            valid,
            to_eval="L(4959)+L(5007)",
        )
        oh, oh_err = oxygen_abundance_log_map(
            o_plus, o_plus_err, o_doubleplus, o_doubleplus_err, valid
        )
        return {
            "n_pix": int(np.count_nonzero(valid_mask)),
            "ebv": ebv,
            "ebv_err": ebv_err,
            "ne": density["value"],
            "ne_err": density["err"],
            "te_nii": nii["value"],
            "te_nii_err": nii["err"],
            "te_oiii": te_oiii_value,
            "te_oiii_err": te_oiii_err_value,
            "value": finite_scalar(oh),
            "err": finite_scalar(oh_err),
        }

    def compute_integrated_k25(valid_mask):
        nii = compute_integrated_nii_te(valid_mask)
        te = nii["value"]
        te_err = nii["err"]
        valid = np.isfinite(te) and 8000.0 <= te <= 13000.0
        if valid:
            te_unit = te / 1.0e4
            abundance = -1.19 * te_unit + 9.68
            abundance_err = np.sqrt(
                (1.19 * te_err / 1.0e4) ** 2 +
                (0.14 * te_unit) ** 2 +
                0.15 ** 2
            )
        else:
            abundance = np.nan
            abundance_err = np.nan
        return {
            "n_pix": nii["n_pix"],
            "ne": nii["density"]["value"],
            "ne_err": nii["density"]["err"],
            "te_nii": te,
            "te_nii_err": te_err,
            "value": abundance,
            "err": abundance_err,
        }

    def compute_integrated_te_summary(products):
        density = compute_integrated_density(products["COMMON_TE"])
        return {
            "NE_SII": density,
            "TE_NII": compute_integrated_nii_te(products["NII_TE_VALID"]),
            "BR24_DIRECT": compute_integrated_oxygen(
                products["BR24_DIRECT_VALID"], use_measured_siii=True
            ),
            "NII_OII7325": compute_integrated_oxygen(
                products["NII_OII7325_VALID"], use_measured_siii=False
            ),
            "NII_K25": compute_integrated_k25(products["NII_TE_VALID"]),
        }

    def format_measurement(value, err=None, precision=4):
        if not np.isfinite(value):
            return "nan"
        if err is not None and np.isfinite(err):
            return f"{value:.{precision}f} +/- {err:.{precision}f}"
        return f"{value:.{precision}f}"

    def print_integrated_te_summary(label, summary):
        ne = summary["NE_SII"]
        te = summary["TE_NII"]
        br24 = summary["BR24_DIRECT"]
        nii_oii = summary["NII_OII7325"]
        k25 = summary["NII_K25"]
        print(f"  {label}:")
        print(
            "    NE_SII: "
            f"{format_measurement(ne['value'], ne['err'], 2)} cm^-3 "
            f"(n_pix={ne['n_pix']}, fixed20={int(ne['fixed20'])}, "
            f"exact_high_density={int(ne['exact_high_density'])})"
        )
        print(
            "    Te(NII): "
            f"{format_measurement(te['value'], te['err'], 1)} K "
            f"(n_pix={te['n_pix']})"
        )
        print(
            "    O_H_BR24_DIRECT: "
            f"{format_measurement(br24['value'], br24['err'], 4)} dex "
            f"(n_pix={br24['n_pix']})"
        )
        print(
            "    O_H_NII_OII7325: "
            f"{format_measurement(nii_oii['value'], nii_oii['err'], 4)} dex "
            f"(n_pix={nii_oii['n_pix']})"
        )
        print(
            "    O_H_NII_K25: "
            f"{format_measurement(k25['value'], k25['err'], 4)} dex "
            f"(n_pix={k25['n_pix']})"
        )

    TE_INTEGRATED_SF = compute_integrated_te_summary(te_sf_products)
    TE_INTEGRATED_HII = compute_integrated_te_summary(te_hii_products)

    ne_sii_finite_count = int(np.count_nonzero(np.isfinite(NE_SII)))
    ne_sii_fixed_count = int(np.count_nonzero(NE_SII_FIXED20))
    ne_sii_fixed_fraction = (
        ne_sii_fixed_count / ne_sii_finite_count
        if ne_sii_finite_count > 0 else np.nan
    )
    print("PyNeb Te-method valid-pixel counts:")
    print(f"  NE_SII finite pixels: {ne_sii_finite_count}")
    print(
        "  NE_SII fixed at 20 cm^-3: "
        f"{ne_sii_fixed_count}/{ne_sii_finite_count} "
        f"({100.0 * ne_sii_fixed_fraction:.1f}%)"
    )
    print(f"  NE_SII exact high-density fallback pixels: {int(np.count_nonzero(NE_SII_EXACT_HIGH_DENSITY))}")
    print(f"  NE_SII_ALL finite pixels: {int(np.count_nonzero(np.isfinite(NE_SII_ALL)))}")
    print(f"  Te(NII) HII pixels: {int(np.count_nonzero(np.isfinite(TE_NII_HII)))}")
    print(f"  Te(NII) SF pixels: {int(np.count_nonzero(np.isfinite(TE_NII_SF)))}")
    print(f"  Brazzini+2024 direct HII pixels: {int(np.count_nonzero(np.isfinite(O_H_BR24_DIRECT_HII)))}")
    print(f"  Brazzini+2024 direct SF pixels: {int(np.count_nonzero(np.isfinite(O_H_BR24_DIRECT_SF)))}")
    print(f"  NII+OII7325 semi-direct HII pixels: {int(np.count_nonzero(np.isfinite(O_H_NII_OII7325_HII)))}")
    print(f"  NII+OII7325 semi-direct SF pixels: {int(np.count_nonzero(np.isfinite(O_H_NII_OII7325_SF)))}")
    print(f"  Kreckel+2025 NII-only HII pixels: {int(np.count_nonzero(np.isfinite(O_H_NII_K25_HII)))}")
    print(f"  Kreckel+2025 NII-only SF pixels: {int(np.count_nonzero(np.isfinite(O_H_NII_K25_SF)))}")
    print("Integrated PyNeb Te-method values:")
    print_integrated_te_summary("SF", TE_INTEGRATED_SF)
    print_integrated_te_summary("HII", TE_INTEGRATED_HII)

# ------------------------------------------------------------------
# 11.  Calculate the total Metallicity in SF regions (Classification 1)
# ------------------------------------------------------------------

# Sum raw line maps in SF regions first, then apply one integrated BD correction.
HB4861_FLUX_SF_total = np.nansum(np.where(mask_SF, HB4861_FLUX, np.nan))
HA6562_FLUX_SF_total = np.nansum(np.where(mask_SF, HA6562_FLUX, np.nan))
OIII5006_FLUX_SF_total = np.nansum(np.where(mask_SF, OIII5006_FLUX, np.nan))
NII6583_FLUX_SF_total = np.nansum(np.where(mask_SF, NII6583_FLUX, np.nan))
SII6716_FLUX_SF_total = np.nansum(np.where(mask_SF, SII6716_FLUX, np.nan))
SII6730_FLUX_SF_total = np.nansum(np.where(mask_SF, SII6730_FLUX, np.nan))
HB4861_FLUX_ERR_SF_total = integrated_flux_error(HB4861_FLUX_ERR, mask_SF)
HA6562_FLUX_ERR_SF_total = integrated_flux_error(HA6562_FLUX_ERR, mask_SF)
OIII5006_FLUX_ERR_SF_total = integrated_flux_error(OIII5006_FLUX_ERR, mask_SF)
NII6583_FLUX_ERR_SF_total = integrated_flux_error(NII6583_FLUX_ERR, mask_SF)
SII6716_FLUX_ERR_SF_total = integrated_flux_error(SII6716_FLUX_ERR, mask_SF)
SII6730_FLUX_ERR_SF_total = integrated_flux_error(SII6730_FLUX_ERR, mask_SF)

if (
    np.isfinite(HB4861_FLUX_SF_total)
    and np.isfinite(HA6562_FLUX_SF_total)
    and HB4861_FLUX_SF_total > 0
    and HA6562_FLUX_SF_total > 0
):
    BD_SF_total = HA6562_FLUX_SF_total / HB4861_FLUX_SF_total
    if BD_SF_total < R_int:
        BD_SF_total = R_int
    E_BV_BD_SF_total = convert_bd_to_ebv(BD_SF_total, k_HB4861, k_HA6562, R_int)
    E_BV_BD_SF_total_ERR = convert_bd_to_ebv_error(
        HA6562_FLUX_SF_total, HB4861_FLUX_SF_total,
        HA6562_FLUX_ERR_SF_total, HB4861_FLUX_ERR_SF_total,
        k_HB4861, k_HA6562
    )
    HB4861_FLUX_corr_SF_total = correct_flux_with_ebv(HB4861_FLUX_SF_total, E_BV_BD_SF_total, k_HB4861)
    HA6562_FLUX_corr_SF_total = correct_flux_with_ebv(HA6562_FLUX_SF_total, E_BV_BD_SF_total, k_HA6562)
    OIII5006_FLUX_corr_SF_total = correct_flux_with_ebv(OIII5006_FLUX_SF_total, E_BV_BD_SF_total, k_OIII5006)
    NII6583_FLUX_corr_SF_total = correct_flux_with_ebv(NII6583_FLUX_SF_total, E_BV_BD_SF_total, k_NII6583)
    SII6716_FLUX_corr_SF_total = correct_flux_with_ebv(SII6716_FLUX_SF_total, E_BV_BD_SF_total, k_SII6716)
    SII6730_FLUX_corr_SF_total = correct_flux_with_ebv(SII6730_FLUX_SF_total, E_BV_BD_SF_total, k_SII6730)
    HB4861_FLUX_ERR_corr_SF_total = correct_flux_error_with_ebv(
        HB4861_FLUX_SF_total, HB4861_FLUX_ERR_SF_total, E_BV_BD_SF_total,
        k_HB4861, E_BV_BD_SF_total_ERR
    )
    HA6562_FLUX_ERR_corr_SF_total = correct_flux_error_with_ebv(
        HA6562_FLUX_SF_total, HA6562_FLUX_ERR_SF_total, E_BV_BD_SF_total,
        k_HA6562, E_BV_BD_SF_total_ERR
    )
    OIII5006_FLUX_ERR_corr_SF_total = correct_flux_error_with_ebv(
        OIII5006_FLUX_SF_total, OIII5006_FLUX_ERR_SF_total, E_BV_BD_SF_total,
        k_OIII5006, E_BV_BD_SF_total_ERR
    )
    NII6583_FLUX_ERR_corr_SF_total = correct_flux_error_with_ebv(
        NII6583_FLUX_SF_total, NII6583_FLUX_ERR_SF_total, E_BV_BD_SF_total,
        k_NII6583, E_BV_BD_SF_total_ERR
    )
    SII6716_FLUX_ERR_corr_SF_total = correct_flux_error_with_ebv(
        SII6716_FLUX_SF_total, SII6716_FLUX_ERR_SF_total, E_BV_BD_SF_total,
        k_SII6716, E_BV_BD_SF_total_ERR
    )
    SII6730_FLUX_ERR_corr_SF_total = correct_flux_error_with_ebv(
        SII6730_FLUX_SF_total, SII6730_FLUX_ERR_SF_total, E_BV_BD_SF_total,
        k_SII6730, E_BV_BD_SF_total_ERR
    )
else:
    BD_SF_total = np.nan
    E_BV_BD_SF_total = np.nan
    E_BV_BD_SF_total_ERR = np.nan
    HB4861_FLUX_corr_SF_total = np.nan
    HA6562_FLUX_corr_SF_total = np.nan
    OIII5006_FLUX_corr_SF_total = np.nan
    NII6583_FLUX_corr_SF_total = np.nan
    SII6716_FLUX_corr_SF_total = np.nan
    SII6730_FLUX_corr_SF_total = np.nan
    HB4861_FLUX_ERR_corr_SF_total = np.nan
    HA6562_FLUX_ERR_corr_SF_total = np.nan
    OIII5006_FLUX_ERR_corr_SF_total = np.nan
    NII6583_FLUX_ERR_corr_SF_total = np.nan
    SII6716_FLUX_ERR_corr_SF_total = np.nan
    SII6730_FLUX_ERR_corr_SF_total = np.nan

# ------------------------------------------------------------------
# 10a. Calculate the total Metallicity in HII regions (Classification 2)
# ------------------------------------------------------------------

# Sum raw line maps in HII regions first, then apply one integrated BD correction.
HB4861_FLUX_HII_total = np.nansum(np.where(mask_HII, HB4861_FLUX, np.nan))
HA6562_FLUX_HII_total = np.nansum(np.where(mask_HII, HA6562_FLUX, np.nan))
OIII5006_FLUX_HII_total = np.nansum(np.where(mask_HII, OIII5006_FLUX, np.nan))
NII6583_FLUX_HII_total = np.nansum(np.where(mask_HII, NII6583_FLUX, np.nan))
SII6716_FLUX_HII_total = np.nansum(np.where(mask_HII, SII6716_FLUX, np.nan))
SII6730_FLUX_HII_total = np.nansum(np.where(mask_HII, SII6730_FLUX, np.nan))
HB4861_FLUX_ERR_HII_total = integrated_flux_error(HB4861_FLUX_ERR, mask_HII)
HA6562_FLUX_ERR_HII_total = integrated_flux_error(HA6562_FLUX_ERR, mask_HII)
OIII5006_FLUX_ERR_HII_total = integrated_flux_error(OIII5006_FLUX_ERR, mask_HII)
NII6583_FLUX_ERR_HII_total = integrated_flux_error(NII6583_FLUX_ERR, mask_HII)
SII6716_FLUX_ERR_HII_total = integrated_flux_error(SII6716_FLUX_ERR, mask_HII)
SII6730_FLUX_ERR_HII_total = integrated_flux_error(SII6730_FLUX_ERR, mask_HII)

if (
    np.isfinite(HB4861_FLUX_HII_total)
    and np.isfinite(HA6562_FLUX_HII_total)
    and HB4861_FLUX_HII_total > 0
    and HA6562_FLUX_HII_total > 0
):
    BD_HII_total = HA6562_FLUX_HII_total / HB4861_FLUX_HII_total
    if BD_HII_total < R_int:
        BD_HII_total = R_int
    E_BV_BD_HII_total = convert_bd_to_ebv(BD_HII_total, k_HB4861, k_HA6562, R_int)
    E_BV_BD_HII_total_ERR = convert_bd_to_ebv_error(
        HA6562_FLUX_HII_total, HB4861_FLUX_HII_total,
        HA6562_FLUX_ERR_HII_total, HB4861_FLUX_ERR_HII_total,
        k_HB4861, k_HA6562
    )
    HB4861_FLUX_corr_HII_total = correct_flux_with_ebv(HB4861_FLUX_HII_total, E_BV_BD_HII_total, k_HB4861)
    HA6562_FLUX_corr_HII_total = correct_flux_with_ebv(HA6562_FLUX_HII_total, E_BV_BD_HII_total, k_HA6562)
    OIII5006_FLUX_corr_HII_total = correct_flux_with_ebv(OIII5006_FLUX_HII_total, E_BV_BD_HII_total, k_OIII5006)
    NII6583_FLUX_corr_HII_total = correct_flux_with_ebv(NII6583_FLUX_HII_total, E_BV_BD_HII_total, k_NII6583)
    SII6716_FLUX_corr_HII_total = correct_flux_with_ebv(SII6716_FLUX_HII_total, E_BV_BD_HII_total, k_SII6716)
    SII6730_FLUX_corr_HII_total = correct_flux_with_ebv(SII6730_FLUX_HII_total, E_BV_BD_HII_total, k_SII6730)
else:
    BD_HII_total = np.nan
    E_BV_BD_HII_total = np.nan
    E_BV_BD_HII_total_ERR = np.nan
    HB4861_FLUX_corr_HII_total = np.nan
    HA6562_FLUX_corr_HII_total = np.nan
    OIII5006_FLUX_corr_HII_total = np.nan
    NII6583_FLUX_corr_HII_total = np.nan
    SII6716_FLUX_corr_HII_total = np.nan
    SII6730_FLUX_corr_HII_total = np.nan

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
else:
    O_H_PG16_SF_total = (a1_lower + a2_lower * log_R3_S2_SF_total + a3_lower * log_N2_SF_total + 
                      (a4_lower + a5_lower * log_R3_S2_SF_total + a6_lower * log_N2_SF_total) * log_S2_SF_total)

# N2S2-N06 metallicity calculation for SF total region
# Calculate N2S2 ratio for total SF region using N2S2-N06 calibration
if (np.isfinite(NII6583_FLUX_corr_SF_total) and np.isfinite(SII6716_FLUX_corr_SF_total) and
    np.isfinite(SII6730_FLUX_corr_SF_total) and 
    NII6583_FLUX_corr_SF_total > 0 and SII6716_FLUX_corr_SF_total > 0 and SII6730_FLUX_corr_SF_total > 0):
    
    sii_total_sf = SII6716_FLUX_corr_SF_total + SII6730_FLUX_corr_SF_total
    n2s2_ratio_sf_total = np.log10(NII6583_FLUX_corr_SF_total / sii_total_sf)
    
    # Solve cubic equation: 0.17963*x³ + 0.58181*x² + 0.74100*x + (-0.25214 - n2s2_ratio_sf_total) = 0
    c3 = 0.17963
    c2 = 0.58181
    c1 = 0.74100
    c0 = -0.25214
    
    poly_coeffs = [c3, c2, c1, (c0 - n2s2_ratio_sf_total)]
    roots = np.roots(poly_coeffs)
    
    # Select the real root (use first real root found)
    real_roots = roots[np.isreal(roots)].real
    if len(real_roots) > 0:
        # Take the first real root without range restrictions
        x_final = real_roots[0]
        O_H_N2S2_N06_SF_total = x_final + 8.69
    else:
        O_H_N2S2_N06_SF_total = np.nan
else:
    O_H_N2S2_N06_SF_total = np.nan

# O3N2-M13 (Marino et al. 2013) metallicity calculation (total)
# Calculate O3N2 ratio for total SF region using M13 calibration
oiii_hb_SF_total = OIII5006_FLUX_corr_SF_total / HB4861_FLUX_corr_SF_total
nii_ha_SF_total = NII6583_FLUX_corr_SF_total / HA6562_FLUX_corr_SF_total
o3n2_ratio_SF_total = np.log10(oiii_hb_SF_total / nii_ha_SF_total)
# Apply O3N2-M13 (Marino et al. 2013) calibration: [O/H] = 8.533 - 0.214 * O3N2
O_H_O3N2_M13_SF_total = 8.533 - 0.214 * o3n2_ratio_SF_total
O_H_O3N2_M13_SF_total = mask_scalar_by_range(
    O_H_O3N2_M13_SF_total, o3n2_ratio_SF_total, -1.1, 1.7
)

# N2-M13 (Marino et al. 2013) metallicity calculation (total)
# Calculate N2 ratio for total SF region using M13 calibration
n2_ratio_SF_total = np.log10(NII6583_FLUX_corr_SF_total / HA6562_FLUX_corr_SF_total)
# Apply N2-M13 (Marino et al. 2013) calibration: [O/H] = 8.743 + 0.462 * N2
O_H_N2_M13_SF_total = 8.743 + 0.462 * n2_ratio_SF_total
O_H_N2_M13_SF_total = mask_scalar_by_range(
    O_H_N2_M13_SF_total, n2_ratio_SF_total, -1.6, -0.2
)

# O3N2-PP04 (Pettini & Pagel 2004) metallicity calculation (total)
# Calculate O3N2 ratio for total SF region using PP04 calibration
# (reuse previously calculated ratios from O3N2-M13)
# Apply O3N2-PP04 (Pettini & Pagel 2004) calibration: [O/H] = 8.73 - 0.32 * O3N2
O_H_O3N2_PP04_SF_total = 8.73 - 0.32 * o3n2_ratio_SF_total
O_H_O3N2_PP04_SF_total = mask_scalar_by_range(
    O_H_O3N2_PP04_SF_total, o3n2_ratio_SF_total, None, 1.9
)

# N2-PP04 (Pettini & Pagel 2004) metallicity calculation (total)
# Calculate N2 ratio for total SF region using PP04 calibration
# (reuse previously calculated ratio from N2-M13)
# Apply N2-PP04 (Pettini & Pagel 2004) calibration: [O/H] = 9.37 + 2.03*N2 + 1.26*N2^2 + 0.32*N2^3
O_H_N2_PP04_SF_total = (9.37 + 2.03 * n2_ratio_SF_total + 
                       1.26 * n2_ratio_SF_total**2 + 
                       0.32 * n2_ratio_SF_total**3)
O_H_N2_PP04_SF_total = mask_scalar_by_range(
    O_H_N2_PP04_SF_total, n2_ratio_SF_total, -2.5, -0.3
)

# C20 metallicity calculations for SF total using the same polynomial solvers
# as the spaxel maps, with integrated line-flux errors propagated in quadrature.
_hb_sf = as_single_pixel(HB4861_FLUX_corr_SF_total)
_ha_sf = as_single_pixel(HA6562_FLUX_corr_SF_total)
_oiii_sf = as_single_pixel(OIII5006_FLUX_corr_SF_total)
_nii_sf = as_single_pixel(NII6583_FLUX_corr_SF_total)
_sii6716_sf = as_single_pixel(SII6716_FLUX_corr_SF_total)
_sii6730_sf = as_single_pixel(SII6730_FLUX_corr_SF_total)
_hb_sf_err = as_single_pixel(HB4861_FLUX_ERR_corr_SF_total)
_ha_sf_err = as_single_pixel(HA6562_FLUX_ERR_corr_SF_total)
_oiii_sf_err = as_single_pixel(OIII5006_FLUX_ERR_corr_SF_total)
_nii_sf_err = as_single_pixel(NII6583_FLUX_ERR_corr_SF_total)
_sii6716_sf_err = as_single_pixel(SII6716_FLUX_ERR_corr_SF_total)
_sii6730_sf_err = as_single_pixel(SII6730_FLUX_ERR_corr_SF_total)

_oh, _err, _ = calculate_o3n2_c20_metallicity(
    _hb_sf, _oiii_sf, _nii_sf, _ha_sf,
    _hb_sf_err, _oiii_sf_err, _nii_sf_err, _ha_sf_err, O_H_D16_SF_total
)
_oh, _err = apply_metallicity_range(_oh, _err)
O_H_O3N2_C20_SF_total = single_pixel_value(_oh)
O_H_O3N2_C20_SF_total_ERR = single_pixel_value(_err)

_oh, _err, _ = calculate_o3s2_c20_metallicity(
    _hb_sf, _oiii_sf, _sii6716_sf, _sii6730_sf,
    _hb_sf_err, _oiii_sf_err, _sii6716_sf_err, _sii6730_sf_err, O_H_D16_SF_total
)
_oh, _err = apply_metallicity_range(_oh, _err)
O_H_O3S2_C20_SF_total = single_pixel_value(_oh)
O_H_O3S2_C20_SF_total_ERR = single_pixel_value(_err)

_oh, _err, _ = calculate_rs32_c20_metallicity(
    _hb_sf, _ha_sf, _oiii_sf, _sii6716_sf, _sii6730_sf,
    _hb_sf_err, _ha_sf_err, _oiii_sf_err, _sii6716_sf_err, _sii6730_sf_err,
    O_H_D16_SF_total
)
_oh, _err = apply_metallicity_range(_oh, _err)
O_H_RS32_C20_SF_total = single_pixel_value(_oh)
O_H_RS32_C20_SF_total_ERR = single_pixel_value(_err)

_oh, _err, _ = calculate_r3_c20_metallicity(
    _hb_sf, _hb_sf_err, _oiii_sf, _oiii_sf_err, O_H_D16_SF_total
)
_oh, _err = apply_metallicity_range(_oh, _err)
O_H_R3_C20_SF_total = single_pixel_value(_oh)
O_H_R3_C20_SF_total_ERR = single_pixel_value(_err)

_oh, _err, _ = calculate_n2_c20_metallicity(
    _ha_sf, _ha_sf_err, _nii_sf, _nii_sf_err, O_H_D16_SF_total
)
_oh, _err = apply_metallicity_range(_oh, _err)
O_H_N2_C20_SF_total = single_pixel_value(_oh)
O_H_N2_C20_SF_total_ERR = single_pixel_value(_err)

_oh, _err, _ = calculate_s2_c20_metallicity(
    _ha_sf, _ha_sf_err, _sii6716_sf, _sii6716_sf_err, _sii6730_sf, _sii6730_sf_err,
    O_H_D16_SF_total
)
_oh, _err = apply_metallicity_range(_oh, _err)
O_H_S2_C20_SF_total = single_pixel_value(_oh)
O_H_S2_C20_SF_total_ERR = single_pixel_value(_err)

(
    O_H_COMBINED_C20_SF_total,
    O_H_COMBINED_C20_SF_total_ERR,
    O_H_COMBINED_C20_SF_total_METHOD,
    O_H_COMBINED_C20_SF_total_NMETHODS,
) = combine_c20_scalar((
    ("O3N2", O_H_O3N2_C20_SF_total, O_H_O3N2_C20_SF_total_ERR),
    ("O3S2", O_H_O3S2_C20_SF_total, O_H_O3S2_C20_SF_total_ERR),
    ("RS32", O_H_RS32_C20_SF_total, O_H_RS32_C20_SF_total_ERR),
    ("R3", O_H_R3_C20_SF_total, O_H_R3_C20_SF_total_ERR),
    ("N2", O_H_N2_C20_SF_total, O_H_N2_C20_SF_total_ERR),
    ("S2", O_H_S2_C20_SF_total, O_H_S2_C20_SF_total_ERR),
))

# Dopita et al. (2016) metallicity calculation (HII total)
y_HII_total = np.log10(NII6583_FLUX_corr_HII_total / (SII6716_FLUX_corr_HII_total + SII6730_FLUX_corr_HII_total)) + 0.264*np.log10(NII6583_FLUX_corr_HII_total / HA6562_FLUX_corr_HII_total)
O_H_D16_HII_total = 8.77 + y_HII_total + 0.45*(y_HII_total + 0.3)**5

# Pilyugin & Grebel (2016) metallicity calculation (HII total)
OIII_scaled_HII_total = 1.33 * OIII5006_FLUX_corr_HII_total
NII_scaled_HII_total = 1.34 * NII6583_FLUX_corr_HII_total
N2_HII_total = NII_scaled_HII_total / HB4861_FLUX_corr_HII_total
S2_HII_total = (SII6716_FLUX_corr_HII_total + SII6730_FLUX_corr_HII_total) / HB4861_FLUX_corr_HII_total
R3_HII_total = OIII_scaled_HII_total / HB4861_FLUX_corr_HII_total
log_R3_S2_HII_total = np.log10(R3_HII_total/S2_HII_total)
log_N2_HII_total = np.log10(N2_HII_total)
log_S2_HII_total = np.log10(S2_HII_total)
O_H_PG16_HII_total = []
if log_N2_HII_total >= -0.6:
    O_H_PG16_HII_total = (a1_upper + a2_upper * log_R3_S2_HII_total + a3_upper * log_N2_HII_total + 
                      (a4_upper + a5_upper * log_R3_S2_HII_total + a6_upper * log_N2_HII_total) * log_S2_HII_total)
else:
    O_H_PG16_HII_total = (a1_lower + a2_lower * log_R3_S2_HII_total + a3_lower * log_N2_HII_total +
                        (a4_lower + a5_lower * log_R3_S2_HII_total + a6_lower * log_N2_HII_total) * log_S2_HII_total)

# N2S2-N06 metallicity calculation for HII total region
if (np.isfinite(NII6583_FLUX_corr_HII_total) and np.isfinite(SII6716_FLUX_corr_HII_total) and
    np.isfinite(SII6730_FLUX_corr_HII_total) and 
    NII6583_FLUX_corr_HII_total > 0 and SII6716_FLUX_corr_HII_total > 0 and SII6730_FLUX_corr_HII_total > 0):
    
    sii_total_hii = SII6716_FLUX_corr_HII_total + SII6730_FLUX_corr_HII_total
    n2s2_ratio_hii_total = np.log10(NII6583_FLUX_corr_HII_total / sii_total_hii)
    
    # Solve cubic equation: 0.17963*x³ + 0.58181*x² + 0.74100*x + (-0.25214 - n2s2_ratio_hii_total) = 0
    c3 = 0.17963
    c2 = 0.58181
    c1 = 0.74100
    c0 = -0.25214
    
    poly_coeffs = [c3, c2, c1, (c0 - n2s2_ratio_hii_total)]
    roots = np.roots(poly_coeffs)
    
    # Select the real root (use first real root found)
    real_roots = roots[np.isreal(roots)].real
    if len(real_roots) > 0:
        # Take the first real root without range restrictions
        x_final = real_roots[0]
        O_H_N2S2_N06_HII_total = x_final + 8.69
    else:
        O_H_N2S2_N06_HII_total = np.nan
else:
    O_H_N2S2_N06_HII_total = np.nan

# Other metallicity calculations for HII total
oiii_hb_HII_total = OIII5006_FLUX_corr_HII_total / HB4861_FLUX_corr_HII_total
nii_ha_HII_total = NII6583_FLUX_corr_HII_total / HA6562_FLUX_corr_HII_total
o3n2_ratio_HII_total = np.log10(oiii_hb_HII_total / nii_ha_HII_total)

# O3N2-M13 (Marino et al. 2013) metallicity calculation (HII total)
O_H_O3N2_M13_HII_total = 8.533 - 0.214 * o3n2_ratio_HII_total
O_H_O3N2_M13_HII_total = mask_scalar_by_range(
    O_H_O3N2_M13_HII_total, o3n2_ratio_HII_total, -1.1, 1.7
)

# N2-M13 (Marino et al. 2013) metallicity calculation (HII total)
n2_ratio_HII_total = np.log10(NII6583_FLUX_corr_HII_total / HA6562_FLUX_corr_HII_total)
O_H_N2_M13_HII_total = 8.743 + 0.462 * n2_ratio_HII_total
O_H_N2_M13_HII_total = mask_scalar_by_range(
    O_H_N2_M13_HII_total, n2_ratio_HII_total, -1.6, -0.2
)

# O3N2-PP04 (Pettini & Pagel 2004) metallicity calculation (HII total)
O_H_O3N2_PP04_HII_total = 8.73 - 0.32 * o3n2_ratio_HII_total
O_H_O3N2_PP04_HII_total = mask_scalar_by_range(
    O_H_O3N2_PP04_HII_total, o3n2_ratio_HII_total, None, 1.9
)

# N2-PP04 (Pettini & Pagel 2004) metallicity calculation (HII total)
O_H_N2_PP04_HII_total = (9.37 + 2.03 * n2_ratio_HII_total + 
                       1.26 * n2_ratio_HII_total**2 + 
                       0.32 * n2_ratio_HII_total**3)
O_H_N2_PP04_HII_total = mask_scalar_by_range(
    O_H_N2_PP04_HII_total, n2_ratio_HII_total, -2.5, -0.3
)

# C20 metallicity calculations for HII total using the same polynomial solvers
# as the spaxel maps. The old linear approximation is intentionally avoided.
HB4861_FLUX_ERR_corr_HII_total = integrated_flux_error(
    HB4861_FLUX_ERR, mask_HII, flux=HB4861_FLUX_HII_total,
    ebv=E_BV_BD_HII_total, k=k_HB4861, ebv_err=E_BV_BD_HII_total_ERR
)
HA6562_FLUX_ERR_corr_HII_total = integrated_flux_error(
    HA6562_FLUX_ERR, mask_HII, flux=HA6562_FLUX_HII_total,
    ebv=E_BV_BD_HII_total, k=k_HA6562, ebv_err=E_BV_BD_HII_total_ERR
)
OIII5006_FLUX_ERR_corr_HII_total = integrated_flux_error(
    OIII5006_FLUX_ERR, mask_HII, flux=OIII5006_FLUX_HII_total,
    ebv=E_BV_BD_HII_total, k=k_OIII5006, ebv_err=E_BV_BD_HII_total_ERR
)
NII6583_FLUX_ERR_corr_HII_total = integrated_flux_error(
    NII6583_FLUX_ERR, mask_HII, flux=NII6583_FLUX_HII_total,
    ebv=E_BV_BD_HII_total, k=k_NII6583, ebv_err=E_BV_BD_HII_total_ERR
)
SII6716_FLUX_ERR_corr_HII_total = integrated_flux_error(
    SII6716_FLUX_ERR, mask_HII, flux=SII6716_FLUX_HII_total,
    ebv=E_BV_BD_HII_total, k=k_SII6716, ebv_err=E_BV_BD_HII_total_ERR
)
SII6730_FLUX_ERR_corr_HII_total = integrated_flux_error(
    SII6730_FLUX_ERR, mask_HII, flux=SII6730_FLUX_HII_total,
    ebv=E_BV_BD_HII_total, k=k_SII6730, ebv_err=E_BV_BD_HII_total_ERR
)

_hb_hii = as_single_pixel(HB4861_FLUX_corr_HII_total)
_ha_hii = as_single_pixel(HA6562_FLUX_corr_HII_total)
_oiii_hii = as_single_pixel(OIII5006_FLUX_corr_HII_total)
_nii_hii = as_single_pixel(NII6583_FLUX_corr_HII_total)
_sii6716_hii = as_single_pixel(SII6716_FLUX_corr_HII_total)
_sii6730_hii = as_single_pixel(SII6730_FLUX_corr_HII_total)
_hb_hii_err = as_single_pixel(HB4861_FLUX_ERR_corr_HII_total)
_ha_hii_err = as_single_pixel(HA6562_FLUX_ERR_corr_HII_total)
_oiii_hii_err = as_single_pixel(OIII5006_FLUX_ERR_corr_HII_total)
_nii_hii_err = as_single_pixel(NII6583_FLUX_ERR_corr_HII_total)
_sii6716_hii_err = as_single_pixel(SII6716_FLUX_ERR_corr_HII_total)
_sii6730_hii_err = as_single_pixel(SII6730_FLUX_ERR_corr_HII_total)

_oh, _err, _ = calculate_o3n2_c20_metallicity(
    _hb_hii, _oiii_hii, _nii_hii, _ha_hii,
    _hb_hii_err, _oiii_hii_err, _nii_hii_err, _ha_hii_err, O_H_D16_HII_total
)
_oh, _err = apply_metallicity_range(_oh, _err)
O_H_O3N2_C20_HII_total = single_pixel_value(_oh)
O_H_O3N2_C20_HII_total_ERR = single_pixel_value(_err)

_oh, _err, _ = calculate_o3s2_c20_metallicity(
    _hb_hii, _oiii_hii, _sii6716_hii, _sii6730_hii,
    _hb_hii_err, _oiii_hii_err, _sii6716_hii_err, _sii6730_hii_err,
    O_H_D16_HII_total
)
_oh, _err = apply_metallicity_range(_oh, _err)
O_H_O3S2_C20_HII_total = single_pixel_value(_oh)
O_H_O3S2_C20_HII_total_ERR = single_pixel_value(_err)

_oh, _err, _ = calculate_rs32_c20_metallicity(
    _hb_hii, _ha_hii, _oiii_hii, _sii6716_hii, _sii6730_hii,
    _hb_hii_err, _ha_hii_err, _oiii_hii_err, _sii6716_hii_err, _sii6730_hii_err,
    O_H_D16_HII_total
)
_oh, _err = apply_metallicity_range(_oh, _err)
O_H_RS32_C20_HII_total = single_pixel_value(_oh)
O_H_RS32_C20_HII_total_ERR = single_pixel_value(_err)

_oh, _err, _ = calculate_r3_c20_metallicity(
    _hb_hii, _hb_hii_err, _oiii_hii, _oiii_hii_err, O_H_D16_HII_total
)
_oh, _err = apply_metallicity_range(_oh, _err)
O_H_R3_C20_HII_total = single_pixel_value(_oh)
O_H_R3_C20_HII_total_ERR = single_pixel_value(_err)

_oh, _err, _ = calculate_n2_c20_metallicity(
    _ha_hii, _ha_hii_err, _nii_hii, _nii_hii_err, O_H_D16_HII_total
)
_oh, _err = apply_metallicity_range(_oh, _err)
O_H_N2_C20_HII_total = single_pixel_value(_oh)
O_H_N2_C20_HII_total_ERR = single_pixel_value(_err)

_oh, _err, _ = calculate_s2_c20_metallicity(
    _ha_hii, _ha_hii_err, _sii6716_hii, _sii6716_hii_err, _sii6730_hii,
    _sii6730_hii_err, O_H_D16_HII_total
)
_oh, _err = apply_metallicity_range(_oh, _err)
O_H_S2_C20_HII_total = single_pixel_value(_oh)
O_H_S2_C20_HII_total_ERR = single_pixel_value(_err)

(
    O_H_COMBINED_C20_HII_total,
    O_H_COMBINED_C20_HII_total_ERR,
    O_H_COMBINED_C20_HII_total_METHOD,
    O_H_COMBINED_C20_HII_total_NMETHODS,
) = combine_c20_scalar((
    ("O3N2", O_H_O3N2_C20_HII_total, O_H_O3N2_C20_HII_total_ERR),
    ("O3S2", O_H_O3S2_C20_HII_total, O_H_O3S2_C20_HII_total_ERR),
    ("RS32", O_H_RS32_C20_HII_total, O_H_RS32_C20_HII_total_ERR),
    ("R3", O_H_R3_C20_HII_total, O_H_R3_C20_HII_total_ERR),
    ("N2", O_H_N2_C20_HII_total, O_H_N2_C20_HII_total_ERR),
    ("S2", O_H_S2_C20_HII_total, O_H_S2_C20_HII_total_ERR),
))

# ------------------------------------------------------------------
# 11.  Calculate the total Metallicity in total available regions
# ------------------------------------------------------------------

# Sum raw line maps in the total region first, then apply one integrated BD correction.
HB4861_FLUX_total = np.nansum(HB4861_FLUX)
HA6562_FLUX_total = np.nansum(HA6562_FLUX)
OIII5006_FLUX_total = np.nansum(OIII5006_FLUX)
NII6583_FLUX_total = np.nansum(NII6583_FLUX)
SII6716_FLUX_total = np.nansum(SII6716_FLUX)
SII6730_FLUX_total = np.nansum(SII6730_FLUX)

# Calculate total raw flux errors by error propagation (sqrt of sum of squares).
HB4861_FLUX_ERR_total = integrated_flux_error(HB4861_FLUX_ERR)
HA6562_FLUX_ERR_total = integrated_flux_error(HA6562_FLUX_ERR)
OIII5006_FLUX_ERR_total = integrated_flux_error(OIII5006_FLUX_ERR)
NII6583_FLUX_ERR_total = integrated_flux_error(NII6583_FLUX_ERR)
SII6716_FLUX_ERR_total = integrated_flux_error(SII6716_FLUX_ERR)
SII6730_FLUX_ERR_total = integrated_flux_error(SII6730_FLUX_ERR)

if (
    np.isfinite(HB4861_FLUX_total)
    and np.isfinite(HA6562_FLUX_total)
    and HB4861_FLUX_total > 0
    and HA6562_FLUX_total > 0
):
    BD_total = HA6562_FLUX_total / HB4861_FLUX_total
    if BD_total < R_int:
        BD_total = R_int
    E_BV_BD_total = convert_bd_to_ebv(BD_total, k_HB4861, k_HA6562, R_int)
    E_BV_BD_total_ERR = convert_bd_to_ebv_error(
        HA6562_FLUX_total, HB4861_FLUX_total,
        HA6562_FLUX_ERR_total, HB4861_FLUX_ERR_total,
        k_HB4861, k_HA6562
    )

    # Correct integrated total fluxes with the uniform E(B-V)
    HB4861_FLUX_corr_total = correct_flux_with_ebv(HB4861_FLUX_total, E_BV_BD_total, k_HB4861)
    HA6562_FLUX_corr_total = correct_flux_with_ebv(HA6562_FLUX_total, E_BV_BD_total, k_HA6562)
    OIII5006_FLUX_corr_total = correct_flux_with_ebv(OIII5006_FLUX_total, E_BV_BD_total, k_OIII5006)
    NII6583_FLUX_corr_total = correct_flux_with_ebv(NII6583_FLUX_total, E_BV_BD_total, k_NII6583)
    SII6716_FLUX_corr_total = correct_flux_with_ebv(SII6716_FLUX_total, E_BV_BD_total, k_SII6716)
    SII6730_FLUX_corr_total = correct_flux_with_ebv(SII6730_FLUX_total, E_BV_BD_total, k_SII6730)
    HB4861_FLUX_ERR_total = correct_flux_error_with_ebv(
        HB4861_FLUX_total, HB4861_FLUX_ERR_total, E_BV_BD_total,
        k_HB4861, E_BV_BD_total_ERR
    )
    HA6562_FLUX_ERR_total = correct_flux_error_with_ebv(
        HA6562_FLUX_total, HA6562_FLUX_ERR_total, E_BV_BD_total,
        k_HA6562, E_BV_BD_total_ERR
    )
    OIII5006_FLUX_ERR_total = correct_flux_error_with_ebv(
        OIII5006_FLUX_total, OIII5006_FLUX_ERR_total, E_BV_BD_total,
        k_OIII5006, E_BV_BD_total_ERR
    )
    NII6583_FLUX_ERR_total = correct_flux_error_with_ebv(
        NII6583_FLUX_total, NII6583_FLUX_ERR_total, E_BV_BD_total,
        k_NII6583, E_BV_BD_total_ERR
    )
    SII6716_FLUX_ERR_total = correct_flux_error_with_ebv(
        SII6716_FLUX_total, SII6716_FLUX_ERR_total, E_BV_BD_total,
        k_SII6716, E_BV_BD_total_ERR
    )
    SII6730_FLUX_ERR_total = correct_flux_error_with_ebv(
        SII6730_FLUX_total, SII6730_FLUX_ERR_total, E_BV_BD_total,
        k_SII6730, E_BV_BD_total_ERR
    )
else:
    BD_total = np.nan
    E_BV_BD_total = np.nan
    E_BV_BD_total_ERR = np.nan
    HB4861_FLUX_corr_total = np.nan
    HA6562_FLUX_corr_total = np.nan
    OIII5006_FLUX_corr_total = np.nan
    NII6583_FLUX_corr_total = np.nan
    SII6716_FLUX_corr_total = np.nan
    SII6730_FLUX_corr_total = np.nan
    HB4861_FLUX_ERR_total = np.nan
    HA6562_FLUX_ERR_total = np.nan
    OIII5006_FLUX_ERR_total = np.nan
    NII6583_FLUX_ERR_total = np.nan
    SII6716_FLUX_ERR_total = np.nan
    SII6730_FLUX_ERR_total = np.nan

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

# N2S2-N06 metallicity calculation for total available regions
if (np.isfinite(NII6583_FLUX_corr_total) and np.isfinite(SII6716_FLUX_corr_total) and
    np.isfinite(SII6730_FLUX_corr_total) and 
    NII6583_FLUX_corr_total > 0 and SII6716_FLUX_corr_total > 0 and SII6730_FLUX_corr_total > 0):
    
    sii_total_region = SII6716_FLUX_corr_total + SII6730_FLUX_corr_total
    n2s2_ratio_total = np.log10(NII6583_FLUX_corr_total / sii_total_region)
    
    # Solve cubic equation: 0.17963*x³ + 0.58181*x² + 0.74100*x + (-0.25214 - n2s2_ratio_total) = 0
    c3 = 0.17963
    c2 = 0.58181
    c1 = 0.74100
    c0 = -0.25214
    
    poly_coeffs = [c3, c2, c1, (c0 - n2s2_ratio_total)]
    roots = np.roots(poly_coeffs)
    
    # Select the real root (use first real root found)
    real_roots = roots[np.isreal(roots)].real
    if len(real_roots) > 0:
        # Take the first real root without range restrictions
        x_final = real_roots[0]
        O_H_N2S2_N06_total = x_final + 8.69
    else:
        O_H_N2S2_N06_total = np.nan
else:
    O_H_N2S2_N06_total = np.nan

# O3N2-M13 (Marino et al. 2013) metallicity calculation (total)
# Calculate O3N2 ratio for total SF region using M13 calibration
oiii_hb_total = OIII5006_FLUX_corr_total / HB4861_FLUX_corr_total
nii_ha_total = NII6583_FLUX_corr_total / HA6562_FLUX_corr_total
o3n2_ratio_total = np.log10(oiii_hb_total / nii_ha_total)
# Apply O3N2-M13 (Marino et al. 2013) calibration: [O/H] = 8.533 - 0.214 * O3N2
O_H_O3N2_M13_total = 8.533 - 0.214 * o3n2_ratio_total
O_H_O3N2_M13_total = mask_scalar_by_range(
    O_H_O3N2_M13_total, o3n2_ratio_total, -1.1, 1.7
)

# N2-M13 (Marino et al. 2013) metallicity calculation (total)
# Calculate N2 ratio for total SF region using M13 calibration
n2_ratio_total = np.log10(NII6583_FLUX_corr_total / HA6562_FLUX_corr_total)
# Apply N2-M13 (Marino et al. 2013) calibration: [O/H] = 8.743 + 0.462 * N2
O_H_N2_M13_total = 8.743 + 0.462 * n2_ratio_total
O_H_N2_M13_total = mask_scalar_by_range(
    O_H_N2_M13_total, n2_ratio_total, -1.6, -0.2
)

# O3N2-PP04 (Pettini & Pagel 2004) metallicity calculation (total)
O_H_O3N2_PP04_total = 8.73 - 0.32 * o3n2_ratio_total
O_H_O3N2_PP04_total = mask_scalar_by_range(
    O_H_O3N2_PP04_total, o3n2_ratio_total, None, 1.9
)

# N2-PP04 (Pettini & Pagel 2004) metallicity calculation (total)
O_H_N2_PP04_total = (9.37 + 2.03 * n2_ratio_total + 
                       1.26 * n2_ratio_total**2 + 
                       0.32 * n2_ratio_total**3)
O_H_N2_PP04_total = mask_scalar_by_range(
    O_H_N2_PP04_total, n2_ratio_total, -2.5, -0.3
)

# C20 metallicity calculations for the total region, using the same solvers
# as the spaxel maps and the integrated-BD-corrected total flux errors.
_hb_total = as_single_pixel(HB4861_FLUX_corr_total)
_ha_total = as_single_pixel(HA6562_FLUX_corr_total)
_oiii_total = as_single_pixel(OIII5006_FLUX_corr_total)
_nii_total = as_single_pixel(NII6583_FLUX_corr_total)
_sii6716_total = as_single_pixel(SII6716_FLUX_corr_total)
_sii6730_total = as_single_pixel(SII6730_FLUX_corr_total)
_hb_total_err = as_single_pixel(HB4861_FLUX_ERR_total)
_ha_total_err = as_single_pixel(HA6562_FLUX_ERR_total)
_oiii_total_err = as_single_pixel(OIII5006_FLUX_ERR_total)
_nii_total_err = as_single_pixel(NII6583_FLUX_ERR_total)
_sii6716_total_err = as_single_pixel(SII6716_FLUX_ERR_total)
_sii6730_total_err = as_single_pixel(SII6730_FLUX_ERR_total)

_oh, _err, _ = calculate_o3n2_c20_metallicity(
    _hb_total, _oiii_total, _nii_total, _ha_total,
    _hb_total_err, _oiii_total_err, _nii_total_err, _ha_total_err,
    O_H_D16_total
)
_oh, _err = apply_metallicity_range(_oh, _err)
O_H_O3N2_C20_total = single_pixel_value(_oh)
O_H_O3N2_C20_total_ERR = single_pixel_value(_err)

_oh, _err, _ = calculate_o3s2_c20_metallicity(
    _hb_total, _oiii_total, _sii6716_total, _sii6730_total,
    _hb_total_err, _oiii_total_err, _sii6716_total_err, _sii6730_total_err,
    O_H_D16_total
)
_oh, _err = apply_metallicity_range(_oh, _err)
O_H_O3S2_C20_total = single_pixel_value(_oh)
O_H_O3S2_C20_total_ERR = single_pixel_value(_err)

_oh, _err, _ = calculate_rs32_c20_metallicity(
    _hb_total, _ha_total, _oiii_total, _sii6716_total, _sii6730_total,
    _hb_total_err, _ha_total_err, _oiii_total_err, _sii6716_total_err,
    _sii6730_total_err, O_H_D16_total
)
_oh, _err = apply_metallicity_range(_oh, _err)
O_H_RS32_C20_total = single_pixel_value(_oh)
O_H_RS32_C20_total_ERR = single_pixel_value(_err)

_oh, _err, _ = calculate_r3_c20_metallicity(
    _hb_total, _hb_total_err, _oiii_total, _oiii_total_err, O_H_D16_total
)
_oh, _err = apply_metallicity_range(_oh, _err)
O_H_R3_C20_total = single_pixel_value(_oh)
O_H_R3_C20_total_ERR = single_pixel_value(_err)

_oh, _err, _ = calculate_n2_c20_metallicity(
    _ha_total, _ha_total_err, _nii_total, _nii_total_err, O_H_D16_total
)
_oh, _err = apply_metallicity_range(_oh, _err)
O_H_N2_C20_total = single_pixel_value(_oh)
O_H_N2_C20_total_ERR = single_pixel_value(_err)

_oh, _err, _ = calculate_s2_c20_metallicity(
    _ha_total, _ha_total_err, _sii6716_total, _sii6716_total_err,
    _sii6730_total, _sii6730_total_err, O_H_D16_total
)
_oh, _err = apply_metallicity_range(_oh, _err)
O_H_S2_C20_total = single_pixel_value(_oh)
O_H_S2_C20_total_ERR = single_pixel_value(_err)

(
    O_H_COMBINED_C20_total,
    O_H_COMBINED_C20_total_ERR,
    O_H_COMBINED_C20_total_METHOD,
    O_H_COMBINED_C20_total_NMETHODS,
) = combine_c20_scalar((
    ("O3N2", O_H_O3N2_C20_total, O_H_O3N2_C20_total_ERR),
    ("O3S2", O_H_O3S2_C20_total, O_H_O3S2_C20_total_ERR),
    ("RS32", O_H_RS32_C20_total, O_H_RS32_C20_total_ERR),
    ("R3", O_H_R3_C20_total, O_H_R3_C20_total_ERR),
    ("N2", O_H_N2_C20_total, O_H_N2_C20_total_ERR),
    ("S2", O_H_S2_C20_total, O_H_S2_C20_total_ERR),
))

# ------------------------------------------------------------------
# 12.  Output the results
# ------------------------------------------------------------------

with fits.open(src) as hdul:
    new_hdul = fits.HDUList([hdu.copy() for hdu in hdul])
ORIGINAL_GAS_HDU_COUNT = len(new_hdul)

# Add provenance information to primary header
new_hdul[0].header['BPTMODE'] = 'both'
new_hdul[0].header['CUT_SN'] = cut
new_hdul[0].header['NOISE'] = noise  # in 1e-20 erg s-1 cm-2
new_hdul[0].header['DIST_MPC'] = DISTANCE_MPC
new_hdul[0].header['DISTREF'] = DISTANCE_REFERENCE
new_hdul[0].header['SFRIMF'] = 'Chabrier'
new_hdul[0].header['SFRCOEF'] = (SFR_HA_CHABRIER_COEFF, 'Halpha SFR coefficient, Msun/yr per erg/s')
new_hdul[0].header['BPTLIMIT'] = 'Low-S/N non-Balmer lines use measured fluxes, not limit-aware BPT'
new_hdul[0].header['SFRNOTE'] = 'All-spaxel SFR includes upper-limit substitutions where Balmer QC fails'
new_hdul[0].header['C20COMB'] = 'ivar+scatter'
new_hdul[0].header['BPTMAPS'] = '-1 unknown, 0 unclassified, positives are stable classes'
new_hdul[0].header['DUSTLAW'] = str(extinction_law)
if ENABLE_PYQZ_METALLICITY_PRODUCTS or ENABLE_NEBULABAYES_METALLICITY_PRODUCTS:
    new_hdul[0].header['MGRIDCMP'] = (
        'model_grid_compat.py',
        'Required exact-version compat layer',
    )
    new_hdul[0].header.add_history(
        'pyqz/NebulaBayes require model_grid_compat.py; site-packages unchanged.'
    )
    new_hdul[0].header.add_history(
        'pyqz/NebulaBayes use existing dust-corrected six-line HII/SF products.'
    )
    new_hdul[0].header.add_history(
        'Their independent-error APIs omit shared Hbeta/Balmer covariance.'
    )
if ENABLE_PYQZ_METALLICITY_PRODUCTS:
    new_hdul[0].header['PYQZVER'] = PYQZ_PACKAGE_VERSION
    new_hdul[0].header['PYQZGRD'] = 'MAPPINGS V 5.0.16'
    new_hdul[0].header['PYQZDG'] = PYQZ_DIAGNOSTIC
    new_hdul[0].header['PYQZPK'] = (PYQZ_PK, 'Fixed log10(P/k / (K cm-3))')
    new_hdul[0].header['PYQZSTR'] = PYQZ_STRUCT
    new_hdul[0].header['PYQZKAP'] = ('inf', 'Kappa electron distribution')
    new_hdul[0].header['PYQZSMP'] = PYQZ_SAMPLING
    new_hdul[0].header['PYQZSRS'] = PYQZ_SRS
    new_hdul[0].header['PYQZKDE'] = PYQZ_KDE_METHOD
    new_hdul[0].header['PYQZSED'] = PYQZ_RANDOM_SEED
    new_hdul[0].header.add_history(
        'pyqz pressure is fixed at log(P/k)=5.0; no pressure map is inferred.'
    )
if ENABLE_NEBULABAYES_METALLICITY_PRODUCTS:
    new_hdul[0].header['NBAYVER'] = NB_PACKAGE_VERSION
    new_hdul[0].header['NBGRID'] = 'MAPPINGS 5.1 HII'
    new_hdul[0].header['NBSHAPE'] = ','.join(
        str(value) for value in NB_INTERPD_GRID_SHAPE
    )
    new_hdul[0].header['NBORDER'] = NB_INTERP_ORDER
    new_hdul[0].header['NBGRERR'] = NB_GRID_ERROR
    new_hdul[0].header['NBPRIOR'] = NB_PRIOR
    new_hdul[0].header['NBDERED'] = NB_DEREDDEN
    new_hdul[0].header['NBNORM'] = NB_NORM_LINE
    new_hdul[0].header.add_history(
        'NebulaBayes grid_error=0.10 enters likelihood; reported chi2 excludes it.'
    )
if ENABLE_JY22_METALLICITY_PRODUCTS:
    if JY22_GRID is None:
        raise RuntimeError("JY22 is enabled but its validated grid metadata is absent.")
    new_hdul[0].header['JY22EN'] = (True, 'JY22/Peng2026 products enabled')
    new_hdul[0].header['JY22GRID'] = (
        JY22_GRID.source_path.name,
        'Peng2026 grid',
    )
    new_hdul[0].header['JY22SHA'] = JY22_GRID.sha256
    new_hdul[0].header['JY22REF'] = ('Zenodo 21717332', 'Released grid archive')
    new_hdul[0].header['JY22PRIR'] = (
        'flat-grid',
        'Flat prior over retained grid nodes',
    )
    new_hdul[0].header['JY22UMIN'] = (
        JY22_GRID.requested_log_u_bounds[0],
        'Requested minimum log U',
    )
    new_hdul[0].header['JY22UMAX'] = (
        JY22_GRID.requested_log_u_bounds[1],
        'Requested maximum log U',
    )
    new_hdul[0].header['JY22UELO'] = (
        JY22_GRID.effective_log_u_bounds[0],
        'Retained minimum grid log U',
    )
    new_hdul[0].header['JY22UEHI'] = (
        JY22_GRID.effective_log_u_bounds[1],
        'Retained maximum grid log U',
    )
    new_hdul[0].header['JY22SOL'] = (
        JY22_GRID.solar_o_h,
        'Solar 12+log(O/H) added to log Z/Zsun',
    )
    new_hdul[0].header['JY22RAT'] = ('N2,S2,R3', 'Internal ratio order')
    new_hdul[0].header['JY22COV'] = ('full', 'Full first-order ratio covariance')
    new_hdul[0].header['JY22EINF'] = (
        JY22_ERROR_INFLATION,
        'Flux-error inflation; MaNGA 1.25 not used',
    )
    new_hdul[0].header['JY22BD'] = (
        JY22_INTRINSIC_BALMER_RATIO,
        'Intrinsic Halpha/Hbeta floor',
    )
    new_hdul[0].header.add_history(
        'JY22 uses raw six-line flux/errors, existing HII/SF masks, no extra cut.'
    )
    new_hdul[0].header.add_history(
        'JY22 corrects SII6716 and SII6730 separately before summing S2.'
    )
    new_hdul[0].header.add_history(
        'JY22 assumes independent raw line errors before Jacobian propagation.'
    )
    new_hdul[0].header.add_history(
        'JY22 uses a flat discrete prior over the retained grid nodes.'
    )
    new_hdul[0].header.add_history(
        'Released grid spans logU=-4 to +1; inference is cut at requested -0.5.'
    )
    new_hdul[0].header.add_history(
        'No -0.5 node exists; the retained upper logU node is -0.5384615385.'
    )
    new_hdul[0].header.add_history(
        'JY22 means are central; p16/p84 interpolate cumulative marginals.'
    )

# Gas E(B-V)
hdu_E_BV_BD = fits.ImageHDU(E_BV_BD.astype(np.float64),
                           header=gas_header, name="Gas_E_BV_BD")
hdu_E_BV_BD.header['BUNIT'] = 'mag'
new_hdul.append(hdu_E_BV_BD)
hdu_E_BV_BD_ERR = fits.ImageHDU(E_BV_BD_ERR.astype(np.float64),
                               header=gas_header, name="Gas_E_BV_BD_ERR")
hdu_E_BV_BD_ERR.header['BUNIT'] = 'mag'
hdu_E_BV_BD_ERR.header['COMMENT'] = '1-sigma uncertainty from Halpha/Hbeta flux errors'
new_hdul.append(hdu_E_BV_BD_ERR)

# Independent BPT classification maps
hdu_NII_BPT = fits.ImageHDU(NII_BPT.astype(np.int16),
                            header=gas_header, name="NII_BPT")
hdu_NII_BPT.header['BUNIT'] = 'class'
hdu_NII_BPT.header['COMMENT'] = '-1=unknown/non-detection, 0=unclassified, 1=HII, 2=Comp, 3=AGN'
hdu_NII_BPT.header['COMMENT'] = 'Uses dust-corrected fluxes; positive classes require stable central +/-1sigma class'
new_hdul.append(hdu_NII_BPT)
hdu_SII_BPT = fits.ImageHDU(SII_BPT.astype(np.int16),
                            header=gas_header, name="SII_BPT")
hdu_SII_BPT.header['BUNIT'] = 'class'
hdu_SII_BPT.header['COMMENT'] = '-1=unknown/non-detection, 0=unclassified, 1=HII, 2=LINER, 3=Seyfert'
hdu_SII_BPT.header['COMMENT'] = 'Uses dust-corrected fluxes; positive classes require stable central +/-1sigma class'
new_hdul.append(hdu_SII_BPT)
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

# Line flux maps for SF regions (Classification 1)
hdu_HB4861_FLUX_corr_SF = fits.ImageHDU(HB4861_FLUX_corr_SF.astype(np.float64),
                                        header=gas_header, name="HB4861_FLUX_corr_SF")
hdu_HB4861_FLUX_corr_SF.header['BUNIT'] = '1e-20 erg s-1 cm-2'
hdu_HB4861_FLUX_corr_SF.header['COMMENT'] = 'H-beta flux in SF regions (Classification 1)'
new_hdul.append(hdu_HB4861_FLUX_corr_SF)
hdu_HA6562_FLUX_corr_SF = fits.ImageHDU(HA6562_FLUX_corr_SF.astype(np.float64),
                                       header=gas_header, name="HA6562_FLUX_corr_SF")
hdu_HA6562_FLUX_corr_SF.header['BUNIT'] = '1e-20 erg s-1 cm-2'
hdu_HA6562_FLUX_corr_SF.header['COMMENT'] = 'H-alpha flux in SF regions (Classification 1)'
new_hdul.append(hdu_HA6562_FLUX_corr_SF)
hdu_OIII5006_FLUX_corr_SF = fits.ImageHDU(OIII5006_FLUX_corr_SF.astype(np.float64),
                                         header=gas_header, name="OIII5006_FLUX_corr_SF")
hdu_OIII5006_FLUX_corr_SF.header['BUNIT'] = '1e-20 erg s-1 cm-2'
hdu_OIII5006_FLUX_corr_SF.header['COMMENT'] = '[OIII]5007 flux in SF regions (Classification 1)'
new_hdul.append(hdu_OIII5006_FLUX_corr_SF)
hdu_NII6583_FLUX_corr_SF = fits.ImageHDU(NII6583_FLUX_corr_SF.astype(np.float64),
                                        header=gas_header, name="NII6583_FLUX_corr_SF")
hdu_NII6583_FLUX_corr_SF.header['BUNIT'] = '1e-20 erg s-1 cm-2'
hdu_NII6583_FLUX_corr_SF.header['COMMENT'] = '[NII]6583 flux in SF regions (Classification 1)'
new_hdul.append(hdu_NII6583_FLUX_corr_SF)
hdu_SII6716_FLUX_corr_SF = fits.ImageHDU(SII6716_FLUX_corr_SF.astype(np.float64),
                                        header=gas_header, name="SII6716_FLUX_corr_SF")
hdu_SII6716_FLUX_corr_SF.header['BUNIT'] = '1e-20 erg s-1 cm-2'
hdu_SII6716_FLUX_corr_SF.header['COMMENT'] = '[SII]6716 flux in SF regions (Classification 1)'
new_hdul.append(hdu_SII6716_FLUX_corr_SF)
hdu_SII6730_FLUX_corr_SF = fits.ImageHDU(SII6730_FLUX_corr_SF.astype(np.float64),
                                        header=gas_header, name="SII6730_FLUX_corr_SF")
hdu_SII6730_FLUX_corr_SF.header['BUNIT'] = '1e-20 erg s-1 cm-2'
hdu_SII6730_FLUX_corr_SF.header['COMMENT'] = '[SII]6730 flux in SF regions (Classification 1)'
new_hdul.append(hdu_SII6730_FLUX_corr_SF)

# Line flux maps for HII regions (Classification 2)
hdu_HB4861_FLUX_corr_HII = fits.ImageHDU(HB4861_FLUX_corr_HII.astype(np.float64),
                                         header=gas_header, name="HB4861_FLUX_corr_HII")
hdu_HB4861_FLUX_corr_HII.header['BUNIT'] = '1e-20 erg s-1 cm-2'
hdu_HB4861_FLUX_corr_HII.header['COMMENT'] = 'H-beta flux in HII regions (Classification 2)'
new_hdul.append(hdu_HB4861_FLUX_corr_HII)
hdu_HA6562_FLUX_corr_HII = fits.ImageHDU(HA6562_FLUX_corr_HII.astype(np.float64),
                                        header=gas_header, name="HA6562_FLUX_corr_HII")
hdu_HA6562_FLUX_corr_HII.header['BUNIT'] = '1e-20 erg s-1 cm-2'
hdu_HA6562_FLUX_corr_HII.header['COMMENT'] = 'H-alpha flux in HII regions (Classification 2)'
new_hdul.append(hdu_HA6562_FLUX_corr_HII)
hdu_OIII5006_FLUX_corr_HII = fits.ImageHDU(OIII5006_FLUX_corr_HII.astype(np.float64),
                                          header=gas_header, name="OIII5006_FLUX_corr_HII")
hdu_OIII5006_FLUX_corr_HII.header['BUNIT'] = '1e-20 erg s-1 cm-2'
hdu_OIII5006_FLUX_corr_HII.header['COMMENT'] = '[OIII]5007 flux in HII regions (Classification 2)'
new_hdul.append(hdu_OIII5006_FLUX_corr_HII)
hdu_NII6583_FLUX_corr_HII = fits.ImageHDU(NII6583_FLUX_corr_HII.astype(np.float64),
                                         header=gas_header, name="NII6583_FLUX_corr_HII")
hdu_NII6583_FLUX_corr_HII.header['BUNIT'] = '1e-20 erg s-1 cm-2'
hdu_NII6583_FLUX_corr_HII.header['COMMENT'] = '[NII]6583 flux in HII regions (Classification 2)'
new_hdul.append(hdu_NII6583_FLUX_corr_HII)
hdu_SII6716_FLUX_corr_HII = fits.ImageHDU(SII6716_FLUX_corr_HII.astype(np.float64),
                                         header=gas_header, name="SII6716_FLUX_corr_HII")
hdu_SII6716_FLUX_corr_HII.header['BUNIT'] = '1e-20 erg s-1 cm-2'
hdu_SII6716_FLUX_corr_HII.header['COMMENT'] = '[SII]6716 flux in HII regions (Classification 2)'
new_hdul.append(hdu_SII6716_FLUX_corr_HII)
hdu_SII6730_FLUX_corr_HII = fits.ImageHDU(SII6730_FLUX_corr_HII.astype(np.float64),
                                         header=gas_header, name="SII6730_FLUX_corr_HII")
hdu_SII6730_FLUX_corr_HII.header['BUNIT'] = '1e-20 erg s-1 cm-2'
hdu_SII6730_FLUX_corr_HII.header['COMMENT'] = '[SII]6730 flux in HII regions (Classification 2)'
new_hdul.append(hdu_SII6730_FLUX_corr_HII)
# Corrected Hα luminosity
hdu_halpha_lum = fits.ImageHDU(HA6562_LUM.astype(np.float64),
                               header=gas_header, name="Halpha_Luminosity_corr")
hdu_halpha_lum.header['BUNIT'] = 'erg/s'
hdu_halpha_lum.header['COMMENT'] = 'Includes Balmer-corrected detections plus upper-limit substitutions where Balmer QC fails'
new_hdul.append(hdu_halpha_lum)
# Corrected SFR
hdu_sfr = fits.ImageHDU(SFR_map.astype(np.float64),
                        header=gas_header, name="Halpha_SFR_corr")
hdu_sfr.header['BUNIT'] = 'M_sun/yr'
hdu_sfr.header['COMMENT'] = 'Uses Chabrier coefficient; all-spaxel map includes upper-limit substitutions'
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
hdu_logsfr_unclassified1 = fits.ImageHDU(LOG_SFR_surface_density_map_unclassified1.astype(np.float64),
                                           header=gas_header, name="LOGSFR_SURFACE_DENSITY_UNCLASSIFIED1")
hdu_logsfr_unclassified1.header['BUNIT'] = 'log(M_sun/yr/kpc2)'
new_hdul.append(hdu_logsfr_unclassified1)

# HII-specific SFR maps (Classification 2)
hdu_logsfr_hii = fits.ImageHDU(LOG_SFR_surface_density_map_HII.astype(np.float64),
                              header=gas_header, name="LOGSFR_SURFACE_DENSITY_HII")
hdu_logsfr_hii.header['BUNIT'] = 'log(M_sun/yr/kpc2)'
hdu_logsfr_hii.header['COMMENT'] = 'SFR surface density in HII regions only (Classification 2)'
new_hdul.append(hdu_logsfr_hii)
hdu_logsfr_nonhii = fits.ImageHDU(LOG_SFR_surface_density_map_nonHII.astype(np.float64),
                                 header=gas_header, name="LOGSFR_SURFACE_DENSITY_NONHII")
hdu_logsfr_nonhii.header['BUNIT'] = 'log(M_sun/yr/kpc2)'
hdu_logsfr_nonhii.header['COMMENT'] = 'SFR surface density in non-HII regions (Classification 2)'
new_hdul.append(hdu_logsfr_nonhii)
hdu_logsfr_unclassified2 = fits.ImageHDU(LOG_SFR_surface_density_map_unclassified2.astype(np.float64),
                                        header=gas_header, name="LOGSFR_SURFACE_DENSITY_UNCLASSIFIED2")
hdu_logsfr_unclassified2.header['BUNIT'] = 'log(M_sun/yr/kpc2)'
hdu_logsfr_unclassified2.header['COMMENT'] = 'SFR surface density in unclassified regions (Classification 2)'
new_hdul.append(hdu_logsfr_unclassified2)
hdu_logsfr_upper = fits.ImageHDU(LOG_SFR_surface_density_map_upper.astype(np.float64),
                                   header=gas_header, name="LOGSFR_SURFACE_DENSITY_UPPER")
hdu_logsfr_upper.header['BUNIT'] = 'log(M_sun/yr/kpc2)'
new_hdul.append(hdu_logsfr_upper)

# Regular SFR maps (not log) - SF regions (Classification 1)
hdu_sfr_sf = fits.ImageHDU(SFR_map_SF.astype(np.float64),
                          header=gas_header, name="Halpha_SFR_corr_SF")
hdu_sfr_sf.header['BUNIT'] = 'M_sun/yr/kpc2'
hdu_sfr_sf.header['COMMENT'] = 'SFR surface density in SF regions (Classification 1)'
new_hdul.append(hdu_sfr_sf)
hdu_sfr_nonsf = fits.ImageHDU(SFR_map_nonSF.astype(np.float64),
                             header=gas_header, name="Halpha_SFR_corr_nonSF")
hdu_sfr_nonsf.header['BUNIT'] = 'M_sun/yr/kpc2'
hdu_sfr_nonsf.header['COMMENT'] = 'SFR surface density in non-SF regions (Classification 1)'
new_hdul.append(hdu_sfr_nonsf)
hdu_sfr_unclassified1 = fits.ImageHDU(SFR_map_unclassified1.astype(np.float64),
                                     header=gas_header, name="Halpha_SFR_corr_unclassified1")
hdu_sfr_unclassified1.header['BUNIT'] = 'M_sun/yr/kpc2'
hdu_sfr_unclassified1.header['COMMENT'] = 'SFR surface density in unclassified regions (Classification 1)'
new_hdul.append(hdu_sfr_unclassified1)

# Regular SFR maps (not log) - HII regions (Classification 2)
hdu_sfr_hii = fits.ImageHDU(SFR_map_HII.astype(np.float64),
                           header=gas_header, name="Halpha_SFR_corr_HII")
hdu_sfr_hii.header['BUNIT'] = 'M_sun/yr/kpc2'
hdu_sfr_hii.header['COMMENT'] = 'SFR surface density in HII regions only (Classification 2)'
new_hdul.append(hdu_sfr_hii)
hdu_sfr_nonhii = fits.ImageHDU(SFR_map_nonHII.astype(np.float64),
                              header=gas_header, name="Halpha_SFR_corr_nonHII")
hdu_sfr_nonhii.header['BUNIT'] = 'M_sun/yr/kpc2'
hdu_sfr_nonhii.header['COMMENT'] = 'SFR surface density in non-HII regions (Classification 2)'
new_hdul.append(hdu_sfr_nonhii)
hdu_sfr_unclassified2 = fits.ImageHDU(SFR_map_unclassified2.astype(np.float64),
                                     header=gas_header, name="Halpha_SFR_corr_unclassified2")
hdu_sfr_unclassified2.header['BUNIT'] = 'M_sun/yr/kpc2'
hdu_sfr_unclassified2.header['COMMENT'] = 'SFR surface density in unclassified regions (Classification 2)'
new_hdul.append(hdu_sfr_unclassified2)

# [O/H]
hdu_O_H_D16_SF = fits.ImageHDU(O_H_D16_SF.astype(np.float64),
                             header=gas_header, name="O_H_D16_SF")
hdu_O_H_D16_SF.header['BUNIT'] = '12+log(O/H)'
new_hdul.append(hdu_O_H_D16_SF)
hdu_O_H_PG16_SF = fits.ImageHDU(O_H_PG16_SF.astype(np.float64),
                             header=gas_header, name="O_H_PG16_SF")
hdu_O_H_PG16_SF.header['BUNIT'] = '12+log(O/H)'
new_hdul.append(hdu_O_H_PG16_SF)
hdu_O_H_N2S2_N06_SF = fits.ImageHDU(O_H_N2S2_N06_SF.astype(np.float64),
                             header=gas_header, name="O_H_N2S2_N06_SF")
hdu_O_H_N2S2_N06_SF.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_N2S2_N06_SF.header['COMMENT'] = 'N2S2-N06 metallicity calibration in SF regions'
new_hdul.append(hdu_O_H_N2S2_N06_SF)
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
hdu_O_H_COMBINED_C20_SF.header['COMMENT'] = 'Inverse-variance combined C20 metallicity in SF regions'
new_hdul.append(hdu_O_H_COMBINED_C20_SF)

hdu_COMBINED_C20_METHOD = fits.ImageHDU(combined_c20_method_map.astype(np.int16),
                             header=gas_header, name="COMBINED_C20_METHOD")
hdu_COMBINED_C20_METHOD.header['BUNIT'] = 'method_index'
hdu_COMBINED_C20_METHOD.header['COMMENT'] = 'Dominant C20 weight: -1=none, 0=O3N2, 1=O3S2, 2=RS32, 3=R3, 4=N2, 5=S2'
new_hdul.append(hdu_COMBINED_C20_METHOD)

# Metallicity error maps for SF regions (Classification 1)
hdu_O_H_O3N2_C20_SF_ERR = fits.ImageHDU(O_H_O3N2_C20_SF_ERR.astype(np.float64),
                             header=gas_header, name="O_H_O3N2_C20_SF_ERR")
hdu_O_H_O3N2_C20_SF_ERR.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_O3N2_C20_SF_ERR.header['COMMENT'] = 'O3N2-C20 metallicity error in SF regions'
new_hdul.append(hdu_O_H_O3N2_C20_SF_ERR)
hdu_O_H_O3S2_C20_SF_ERR = fits.ImageHDU(O_H_O3S2_C20_SF_ERR.astype(np.float64),
                             header=gas_header, name="O_H_O3S2_C20_SF_ERR")
hdu_O_H_O3S2_C20_SF_ERR.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_O3S2_C20_SF_ERR.header['COMMENT'] = 'O3S2-C20 metallicity error in SF regions'
new_hdul.append(hdu_O_H_O3S2_C20_SF_ERR)
hdu_O_H_RS32_C20_SF_ERR = fits.ImageHDU(O_H_RS32_C20_SF_ERR.astype(np.float64),
                             header=gas_header, name="O_H_RS32_C20_SF_ERR")
hdu_O_H_RS32_C20_SF_ERR.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_RS32_C20_SF_ERR.header['COMMENT'] = 'RS32-C20 metallicity error in SF regions'
new_hdul.append(hdu_O_H_RS32_C20_SF_ERR)
hdu_O_H_R3_C20_SF_ERR = fits.ImageHDU(O_H_R3_C20_SF_ERR.astype(np.float64),
                             header=gas_header, name="O_H_R3_C20_SF_ERR")
hdu_O_H_R3_C20_SF_ERR.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_R3_C20_SF_ERR.header['COMMENT'] = 'R3-C20 metallicity error in SF regions'
new_hdul.append(hdu_O_H_R3_C20_SF_ERR)
hdu_O_H_N2_C20_SF_ERR = fits.ImageHDU(O_H_N2_C20_SF_ERR.astype(np.float64),
                             header=gas_header, name="O_H_N2_C20_SF_ERR")
hdu_O_H_N2_C20_SF_ERR.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_N2_C20_SF_ERR.header['COMMENT'] = 'N2-C20 metallicity error in SF regions'
new_hdul.append(hdu_O_H_N2_C20_SF_ERR)
hdu_O_H_S2_C20_SF_ERR = fits.ImageHDU(O_H_S2_C20_SF_ERR.astype(np.float64),
                             header=gas_header, name="O_H_S2_C20_SF_ERR")
hdu_O_H_S2_C20_SF_ERR.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_S2_C20_SF_ERR.header['COMMENT'] = 'S2-C20 metallicity error in SF regions'
new_hdul.append(hdu_O_H_S2_C20_SF_ERR)
hdu_O_H_COMBINED_C20_SF_ERR = fits.ImageHDU(O_H_COMBINED_C20_SF_ERR.astype(np.float64),
                             header=gas_header, name="O_H_COMBINED_C20_SF_ERR")
hdu_O_H_COMBINED_C20_SF_ERR.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_COMBINED_C20_SF_ERR.header['COMMENT'] = 'Combined C20 error including formal weight and method scatter'
new_hdul.append(hdu_O_H_COMBINED_C20_SF_ERR)

# HII-specific metallicity maps (Classification 2)
hdu_O_H_D16_HII = fits.ImageHDU(O_H_D16_HII.astype(np.float64),
                             header=gas_header, name="O_H_D16_HII")
hdu_O_H_D16_HII.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_D16_HII.header['COMMENT'] = 'D16 metallicity in HII regions only (Classification 2)'
new_hdul.append(hdu_O_H_D16_HII)
hdu_O_H_PG16_HII = fits.ImageHDU(O_H_PG16_HII.astype(np.float64),
                             header=gas_header, name="O_H_PG16_HII")
hdu_O_H_PG16_HII.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_PG16_HII.header['COMMENT'] = 'PG16 metallicity in HII regions only (Classification 2)'
new_hdul.append(hdu_O_H_PG16_HII)
hdu_O_H_N2S2_N06_HII = fits.ImageHDU(O_H_N2S2_N06_HII.astype(np.float64),
                             header=gas_header, name="O_H_N2S2_N06_HII")
hdu_O_H_N2S2_N06_HII.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_N2S2_N06_HII.header['COMMENT'] = 'N2S2-N06 metallicity in HII regions only (Classification 2)'
new_hdul.append(hdu_O_H_N2S2_N06_HII)
hdu_O_H_O3N2_M13_HII = fits.ImageHDU(O_H_O3N2_M13_HII.astype(np.float64),
                             header=gas_header, name="O_H_O3N2_M13_HII")
hdu_O_H_O3N2_M13_HII.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_O3N2_M13_HII.header['COMMENT'] = 'O3N2-M13 metallicity in HII regions only (Classification 2)'
new_hdul.append(hdu_O_H_O3N2_M13_HII)
hdu_O_H_N2_M13_HII = fits.ImageHDU(O_H_N2_M13_HII.astype(np.float64),
                             header=gas_header, name="O_H_N2_M13_HII")
hdu_O_H_N2_M13_HII.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_N2_M13_HII.header['COMMENT'] = 'N2-M13 metallicity in HII regions only (Classification 2)'
new_hdul.append(hdu_O_H_N2_M13_HII)
hdu_O_H_O3N2_PP04_HII = fits.ImageHDU(O_H_O3N2_PP04_HII.astype(np.float64),
                             header=gas_header, name="O_H_O3N2_PP04_HII")
hdu_O_H_O3N2_PP04_HII.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_O3N2_PP04_HII.header['COMMENT'] = 'O3N2-PP04 metallicity in HII regions only (Classification 2)'
new_hdul.append(hdu_O_H_O3N2_PP04_HII)
hdu_O_H_N2_PP04_HII = fits.ImageHDU(O_H_N2_PP04_HII.astype(np.float64),
                             header=gas_header, name="O_H_N2_PP04_HII")
hdu_O_H_N2_PP04_HII.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_N2_PP04_HII.header['COMMENT'] = 'N2-PP04 metallicity in HII regions only (Classification 2)'
new_hdul.append(hdu_O_H_N2_PP04_HII)
hdu_O_H_O3N2_C20_HII = fits.ImageHDU(O_H_O3N2_C20_HII.astype(np.float64),
                             header=gas_header, name="O_H_O3N2_C20_HII")
hdu_O_H_O3N2_C20_HII.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_O3N2_C20_HII.header['COMMENT'] = 'O3N2-C20 metallicity in HII regions only (Classification 2)'
new_hdul.append(hdu_O_H_O3N2_C20_HII)
hdu_O_H_O3S2_C20_HII = fits.ImageHDU(O_H_O3S2_C20_HII.astype(np.float64),
                             header=gas_header, name="O_H_O3S2_C20_HII")
hdu_O_H_O3S2_C20_HII.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_O3S2_C20_HII.header['COMMENT'] = 'O3S2-C20 metallicity in HII regions only (Classification 2)'
new_hdul.append(hdu_O_H_O3S2_C20_HII)
hdu_O_H_RS32_C20_HII = fits.ImageHDU(O_H_RS32_C20_HII.astype(np.float64),
                             header=gas_header, name="O_H_RS32_C20_HII")
hdu_O_H_RS32_C20_HII.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_RS32_C20_HII.header['COMMENT'] = 'RS32-C20 metallicity in HII regions only (Classification 2)'
new_hdul.append(hdu_O_H_RS32_C20_HII)
hdu_O_H_R3_C20_HII = fits.ImageHDU(O_H_R3_C20_HII.astype(np.float64),
                             header=gas_header, name="O_H_R3_C20_HII")
hdu_O_H_R3_C20_HII.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_R3_C20_HII.header['COMMENT'] = 'R3-C20 metallicity in HII regions only (Classification 2)'
new_hdul.append(hdu_O_H_R3_C20_HII)
hdu_O_H_N2_C20_HII = fits.ImageHDU(O_H_N2_C20_HII.astype(np.float64),
                             header=gas_header, name="O_H_N2_C20_HII")
hdu_O_H_N2_C20_HII.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_N2_C20_HII.header['COMMENT'] = 'N2-C20 metallicity in HII regions only (Classification 2)'
new_hdul.append(hdu_O_H_N2_C20_HII)
hdu_O_H_S2_C20_HII = fits.ImageHDU(O_H_S2_C20_HII.astype(np.float64),
                             header=gas_header, name="O_H_S2_C20_HII")
hdu_O_H_S2_C20_HII.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_S2_C20_HII.header['COMMENT'] = 'S2-C20 metallicity in HII regions only (Classification 2)'
new_hdul.append(hdu_O_H_S2_C20_HII)
hdu_O_H_COMBINED_C20_HII = fits.ImageHDU(O_H_COMBINED_C20_HII.astype(np.float64),
                             header=gas_header, name="O_H_COMBINED_C20_HII")
hdu_O_H_COMBINED_C20_HII.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_COMBINED_C20_HII.header['COMMENT'] = 'Inverse-variance combined C20 metallicity in HII regions only (Classification 2)'
new_hdul.append(hdu_O_H_COMBINED_C20_HII)

# Metallicity error maps for HII regions (Classification 2)
hdu_O_H_O3N2_C20_HII_ERR = fits.ImageHDU(O_H_O3N2_C20_HII_ERR.astype(np.float64),
                             header=gas_header, name="O_H_O3N2_C20_HII_ERR")
hdu_O_H_O3N2_C20_HII_ERR.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_O3N2_C20_HII_ERR.header['COMMENT'] = 'O3N2-C20 metallicity error in HII regions (Classification 2)'
new_hdul.append(hdu_O_H_O3N2_C20_HII_ERR)
hdu_O_H_O3S2_C20_HII_ERR = fits.ImageHDU(O_H_O3S2_C20_HII_ERR.astype(np.float64),
                             header=gas_header, name="O_H_O3S2_C20_HII_ERR")
hdu_O_H_O3S2_C20_HII_ERR.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_O3S2_C20_HII_ERR.header['COMMENT'] = 'O3S2-C20 metallicity error in HII regions (Classification 2)'
new_hdul.append(hdu_O_H_O3S2_C20_HII_ERR)
hdu_O_H_RS32_C20_HII_ERR = fits.ImageHDU(O_H_RS32_C20_HII_ERR.astype(np.float64),
                             header=gas_header, name="O_H_RS32_C20_HII_ERR")
hdu_O_H_RS32_C20_HII_ERR.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_RS32_C20_HII_ERR.header['COMMENT'] = 'RS32-C20 metallicity error in HII regions (Classification 2)'
new_hdul.append(hdu_O_H_RS32_C20_HII_ERR)
hdu_O_H_R3_C20_HII_ERR = fits.ImageHDU(O_H_R3_C20_HII_ERR.astype(np.float64),
                             header=gas_header, name="O_H_R3_C20_HII_ERR")
hdu_O_H_R3_C20_HII_ERR.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_R3_C20_HII_ERR.header['COMMENT'] = 'R3-C20 metallicity error in HII regions (Classification 2)'
new_hdul.append(hdu_O_H_R3_C20_HII_ERR)
hdu_O_H_N2_C20_HII_ERR = fits.ImageHDU(O_H_N2_C20_HII_ERR.astype(np.float64),
                             header=gas_header, name="O_H_N2_C20_HII_ERR")
hdu_O_H_N2_C20_HII_ERR.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_N2_C20_HII_ERR.header['COMMENT'] = 'N2-C20 metallicity error in HII regions (Classification 2)'
new_hdul.append(hdu_O_H_N2_C20_HII_ERR)
hdu_O_H_S2_C20_HII_ERR = fits.ImageHDU(O_H_S2_C20_HII_ERR.astype(np.float64),
                             header=gas_header, name="O_H_S2_C20_HII_ERR")
hdu_O_H_S2_C20_HII_ERR.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_S2_C20_HII_ERR.header['COMMENT'] = 'S2-C20 metallicity error in HII regions (Classification 2)'
new_hdul.append(hdu_O_H_S2_C20_HII_ERR)
hdu_O_H_COMBINED_C20_HII_ERR = fits.ImageHDU(O_H_COMBINED_C20_HII_ERR.astype(np.float64),
                             header=gas_header, name="O_H_COMBINED_C20_HII_ERR")
hdu_O_H_COMBINED_C20_HII_ERR.header['BUNIT'] = '12+log(O/H)'
hdu_O_H_COMBINED_C20_HII_ERR.header['COMMENT'] = 'Combined C20 error including formal weight and method scatter'
new_hdul.append(hdu_O_H_COMBINED_C20_HII_ERR)

if ENABLE_TE_METALLICITY_PRODUCTS:
    te_hdu_specs = [
        ("NE_SII", NE_SII, "cm-3", "PyNeb [S II] 6716/6730 density; finite values below 20 cm-3 are clamped"),
        ("NE_SII_ALL", NE_SII_ALL, "cm-3", "NE_SII where measured; otherwise 20 cm-3 for spatially usable pixels"),
        ("NE_SII_ERR", NE_SII_ERR, "cm-3", "First-order [S II] density uncertainty from line-ratio error"),
        ("TE_NII_HII", TE_NII_HII, "K", "HII-gated Te([N II]) from 5755/(6548+6584) using PyNeb"),
        ("TE_NII_HII_ERR", TE_NII_HII_ERR, "K", "First-order Te([N II]) uncertainty from line-ratio and density errors"),
        ("TE_SIII_BR24_HII", TE_SIII_BR24_HII, "K", "HII-gated Te([S III]) from 6312/(9069+9532); 9532=2.5*9069"),
        ("TE_SIII_BR24_HII_ERR", TE_SIII_BR24_HII_ERR, "K", "First-order Te([S III]) uncertainty from line-ratio and density errors"),
        ("TE_OIII_BR24_HII", TE_OIII_BR24_HII, "K", "Brazzini+2024 Te([O III]) inferred from measured Te([S III])"),
        ("TE_OIII_BR24_HII_ERR", TE_OIII_BR24_HII_ERR, "K", "Includes Brazzini+2024 relation coefficient errors and intrinsic scatter"),
        ("TE_OIII_NII_CHAIN_HII", TE_OIII_NII_CHAIN_HII, "K", "MAUVE semi-direct Te([O III]) from Te([N II]) via Brazzini+2024 Te chains"),
        ("TE_OIII_NII_CHAIN_HII_ERR", TE_OIII_NII_CHAIN_HII_ERR, "K", "Includes propagated Te([N II]) error and Brazzini+2024 relation scatters"),
        ("O_PLUS_BR24_HII", O_PLUS_BR24_HII, "O+/H+", "Brazzini-style O+ from [O II] 7318+7319+7329+7330"),
        ("O_PLUS_BR24_HII_ERR", O_PLUS_BR24_HII_ERR, "O+/H+", "Propagated O+ uncertainty"),
        ("O_DPLUS_BR24_HII", O_DOUBLEPLUS_BR24_HII, "O++/H+", "Brazzini-style O++ from [O III] 4959+5007"),
        ("O_DPLUS_BR24_HII_ERR", O_DOUBLEPLUS_BR24_HII_ERR, "O++/H+", "Propagated O++ uncertainty"),
        ("O_H_BR24_DIRECT_HII", O_H_BR24_DIRECT_HII, "12+log(O/H)", "Brazzini+2024 multi-zone direct-style O/H in HII regions"),
        ("O_H_BR24_DIRECT_HII_ERR", O_H_BR24_DIRECT_HII_ERR, "dex", "Propagated Brazzini+2024 direct-style O/H uncertainty"),
        ("O_PLUS_NII_OII7325_HII", O_PLUS_NII_OII7325_HII, "O+/H+", "Semi-direct O+ from [O II] 7325 using Te([N II])"),
        ("O_PLUS_NII_OII7325_HII_ERR", O_PLUS_NII_OII7325_HII_ERR, "O+/H+", "Propagated semi-direct O+ uncertainty"),
        ("O_DPLUS_NII_OII7325_HII", O_DOUBLEPLUS_NII_OII7325_HII, "O++/H+", "Semi-direct O++ from [O III] with Te([O III]) inferred from Te([N II])"),
        ("O_DPLUS_NII_OII7325_HII_ERR", O_DOUBLEPLUS_NII_OII7325_HII_ERR, "O++/H+", "Propagated semi-direct O++ uncertainty"),
        ("O_H_NII_OII7325_HII", O_H_NII_OII7325_HII, "12+log(O/H)", "MAUVE NII+OII7325 semi-direct ionic-sum O/H in HII regions"),
        ("O_H_NII_OII7325_HII_ERR", O_H_NII_OII7325_HII_ERR, "dex", "Propagated NII+OII7325 semi-direct O/H uncertainty"),
        ("O_H_NII_K25_HII", O_H_NII_K25_HII, "12+log(O/H)", "Kreckel+2025/Mendez-Delgado+2023 Te([N II]) metallicity proxy"),
        ("O_H_NII_K25_HII_ERR", O_H_NII_K25_HII_ERR, "dex", "Includes Te([N II]) and calibration coefficient uncertainties"),
        ("TE_NII_SF", TE_NII_SF, "K", "SF-gated Te([N II]) from 5755/(6548+6584) using PyNeb"),
        ("TE_NII_SF_ERR", TE_NII_SF_ERR, "K", "First-order Te([N II]) uncertainty from line-ratio and density errors"),
        ("TE_SIII_BR24_SF", TE_SIII_BR24_SF, "K", "SF-gated Te([S III]) from 6312/(9069+9532); 9532=2.5*9069"),
        ("TE_SIII_BR24_SF_ERR", TE_SIII_BR24_SF_ERR, "K", "First-order Te([S III]) uncertainty from line-ratio and density errors"),
        ("TE_OIII_BR24_SF", TE_OIII_BR24_SF, "K", "Brazzini+2024 Te([O III]) inferred from measured Te([S III]) in SF regions"),
        ("TE_OIII_BR24_SF_ERR", TE_OIII_BR24_SF_ERR, "K", "Includes Brazzini+2024 relation coefficient errors and intrinsic scatter"),
        ("TE_OIII_NII_CHAIN_SF", TE_OIII_NII_CHAIN_SF, "K", "MAUVE semi-direct Te([O III]) from Te([N II]) via Brazzini+2024 Te chains in SF regions"),
        ("TE_OIII_NII_CHAIN_SF_ERR", TE_OIII_NII_CHAIN_SF_ERR, "K", "Includes propagated Te([N II]) error and Brazzini+2024 relation scatters"),
        ("O_PLUS_BR24_SF", O_PLUS_BR24_SF, "O+/H+", "Brazzini-style O+ from [O II] 7318+7319+7329+7330 in SF regions"),
        ("O_PLUS_BR24_SF_ERR", O_PLUS_BR24_SF_ERR, "O+/H+", "Propagated O+ uncertainty"),
        ("O_DPLUS_BR24_SF", O_DOUBLEPLUS_BR24_SF, "O++/H+", "Brazzini-style O++ from [O III] 4959+5007 in SF regions"),
        ("O_DPLUS_BR24_SF_ERR", O_DOUBLEPLUS_BR24_SF_ERR, "O++/H+", "Propagated O++ uncertainty"),
        ("O_H_BR24_DIRECT_SF", O_H_BR24_DIRECT_SF, "12+log(O/H)", "Brazzini+2024 multi-zone direct-style O/H in SF regions"),
        ("O_H_BR24_DIRECT_SF_ERR", O_H_BR24_DIRECT_SF_ERR, "dex", "Propagated Brazzini+2024 direct-style O/H uncertainty"),
        ("O_PLUS_NII_OII7325_SF", O_PLUS_NII_OII7325_SF, "O+/H+", "Semi-direct O+ from [O II] 7325 using Te([N II]) in SF regions"),
        ("O_PLUS_NII_OII7325_SF_ERR", O_PLUS_NII_OII7325_SF_ERR, "O+/H+", "Propagated semi-direct O+ uncertainty"),
        ("O_DPLUS_NII_OII7325_SF", O_DOUBLEPLUS_NII_OII7325_SF, "O++/H+", "Semi-direct O++ from [O III] with Te([O III]) inferred from Te([N II]) in SF regions"),
        ("O_DPLUS_NII_OII7325_SF_ERR", O_DOUBLEPLUS_NII_OII7325_SF_ERR, "O++/H+", "Propagated semi-direct O++ uncertainty"),
        ("O_H_NII_OII7325_SF", O_H_NII_OII7325_SF, "12+log(O/H)", "MAUVE NII+OII7325 semi-direct ionic-sum O/H in SF regions"),
        ("O_H_NII_OII7325_SF_ERR", O_H_NII_OII7325_SF_ERR, "dex", "Propagated NII+OII7325 semi-direct O/H uncertainty"),
        ("O_H_NII_K25_SF", O_H_NII_K25_SF, "12+log(O/H)", "Kreckel+2025/Mendez-Delgado+2023 Te([N II]) metallicity proxy in SF regions"),
        ("O_H_NII_K25_SF_ERR", O_H_NII_K25_SF_ERR, "dex", "Includes Te([N II]) and calibration coefficient uncertainties"),
    ]
    for name, data, bunit, comment in te_hdu_specs:
        hdu = fits.ImageHDU(np.asarray(data, dtype=np.float64), header=gas_header, name=name)
        hdu.header["BUNIT"] = bunit
        hdu.header["REF1"] = "Brazzini+2024 A&A 691 A173"
        hdu.header["REF2"] = "Kreckel+2025 A&A 703 A42"
        hdu.header["REF3"] = "Mendez-Delgado+2023 Nature"
        hdu.header["COMMENT"] = comment
        new_hdul.append(hdu)

    te_mask_specs = [
        ("NE_SII_FIXED20", NE_SII_FIXED20, "1 where NE_SII was fixed to 20 cm-3"),
        ("NE_SII_EXACT_HIDEN", NE_SII_EXACT_HIGH_DENSITY, "1 where NE_SII used exact PyNeb getTemDen high-density fallback"),
        ("BR24_DIRECT_VALID_HII", BR24_DIRECT_VALID_HII, "1 where all HII Brazzini+2024 direct-style inputs passed detection gates"),
        ("BR24_DIRECT_VALID_SF", BR24_DIRECT_VALID_SF, "1 where all SF Brazzini+2024 direct-style inputs passed detection gates"),
        ("NII_OII7325_VALID_HII", NII_OII7325_VALID_HII, "1 where all HII NII+OII7325 semi-direct inputs passed detection gates"),
        ("NII_OII7325_VALID_SF", NII_OII7325_VALID_SF, "1 where all SF NII+OII7325 semi-direct inputs passed detection gates"),
        ("NII_K25_VALID_HII", NII_K25_VALID_HII, "1 where all HII Kreckel+2025 NII-only inputs passed detection gates"),
        ("NII_K25_VALID_SF", NII_K25_VALID_SF, "1 where all SF Kreckel+2025 NII-only inputs passed detection gates"),
    ]
    for name, data, comment in te_mask_specs:
        hdu = fits.ImageHDU(np.asarray(data, dtype=np.int16), header=gas_header, name=name)
        hdu.header["BUNIT"] = "0/1"
        hdu.header["COMMENT"] = comment
        new_hdul.append(hdu)


CARTA_DIMENSIONLESS_BUNIT = "1"
FLUX_BUNIT = "1e-20 erg s-1 cm-2"


def append_ordered_image(
    hdul: fits.HDUList,
    name: str,
    data,
    bunit: str,
    comment: str | list[str] | tuple[str, ...] | None = None,
    dtype=np.float64,
    refs: bool = False,
    references: tuple[str, ...] | list[str] | None = None,
) -> None:
    hdu = fits.ImageHDU(np.asarray(data, dtype=dtype), header=gas_header, name=name)
    hdu.header["BUNIT"] = bunit
    if refs and references is not None:
        raise ValueError("Use either refs=True or method-specific references, not both.")
    if refs:
        references = (
            "Brazzini+2024 A&A 691 A173",
            "Kreckel+2025 A&A 703 A42",
            "Mendez-Delgado+2023 Nature",
        )
    if references is not None:
        for index, reference in enumerate(references, start=1):
            hdu.header[f"REF{index}"] = reference
    if comment is not None:
        comments = comment if isinstance(comment, (list, tuple)) else (comment,)
        for item in comments:
            hdu.header["COMMENT"] = item
    hdul.append(hdu)


def append_region_pair(
    hdul: fits.HDUList,
    stem: str,
    hii_data,
    sf_data,
    bunit: str,
    hii_comment: str,
    sf_comment: str,
    dtype=np.float64,
    refs: bool = False,
    references: tuple[str, ...] | list[str] | None = None,
) -> None:
    append_ordered_image(
        hdul, f"{stem}_HII", hii_data, bunit, hii_comment, dtype, refs, references
    )
    append_ordered_image(
        hdul, f"{stem}_SF", sf_data, bunit, sf_comment, dtype, refs, references
    )


def append_named_hii_sf_pair(
    hdul: fits.HDUList,
    hii_name: str,
    sf_name: str,
    hii_data,
    sf_data,
    bunit: str,
    hii_comment: str | list[str] | tuple[str, ...],
    sf_comment: str | list[str] | tuple[str, ...],
    dtype=np.float64,
    refs: bool = False,
    references: tuple[str, ...] | list[str] | None = None,
) -> None:
    append_ordered_image(
        hdul, hii_name, hii_data, bunit, hii_comment, dtype, refs, references
    )
    append_ordered_image(
        hdul, sf_name, sf_data, bunit, sf_comment, dtype, refs, references
    )


def build_ordered_output_hdul(base_hdus) -> fits.HDUList:
    ordered_hdul = fits.HDUList([hdu.copy() for hdu in base_hdus])

    # 1. Gas E(B-V).
    append_ordered_image(ordered_hdul, "Gas_E_BV_BD", E_BV_BD, "mag")
    append_ordered_image(
        ordered_hdul,
        "Gas_E_BV_BD_ERR",
        E_BV_BD_ERR,
        "mag",
        "1-sigma uncertainty from Halpha/Hbeta flux errors",
    )

    # 2. Dust-corrected line fluxes and errors for every fitted line.
    for line_base in GAS_LINE_BASES:
        append_ordered_image(
            ordered_hdul,
            f"{line_base}_FLUX_corr",
            globals()[f"{line_base}_FLUX_corr"],
            FLUX_BUNIT,
            f"{line_base} dust-corrected flux from Balmer-decrement E(B-V)",
        )
        append_ordered_image(
            ordered_hdul,
            f"{line_base}_FLUX_ERR_corr",
            globals()[f"{line_base}_FLUX_ERR_corr"],
            FLUX_BUNIT,
            f"{line_base} dust-corrected 1-sigma flux uncertainty",
        )
        append_region_pair(
            ordered_hdul,
            f"{line_base}_FLUX_corr",
            globals()[f"{line_base}_FLUX_corr_HII"],
            globals()[f"{line_base}_FLUX_corr_SF"],
            FLUX_BUNIT,
            f"{line_base} dust-corrected flux in HII regions (Classification 2)",
            f"{line_base} dust-corrected flux in SF regions (Classification 1)",
        )
        append_region_pair(
            ordered_hdul,
            f"{line_base}_FLUX_ERR_corr",
            globals()[f"{line_base}_FLUX_ERR_corr_HII"],
            globals()[f"{line_base}_FLUX_ERR_corr_SF"],
            FLUX_BUNIT,
            f"{line_base} dust-corrected flux uncertainty in HII regions",
            f"{line_base} dust-corrected flux uncertainty in SF regions",
        )

    # 3. BPT and SFR products.
    append_ordered_image(
        ordered_hdul,
        "NII_BPT",
        NII_BPT,
        CARTA_DIMENSIONLESS_BUNIT,
        [
            "Class map: -1=unknown/non-detection, 0=unclassified, 1=HII, 2=Comp, 3=AGN",
            "Uses dust-corrected fluxes; positive classes require stable central +/-1sigma class",
        ],
        dtype=np.int16,
    )
    append_ordered_image(
        ordered_hdul,
        "SII_BPT",
        SII_BPT,
        CARTA_DIMENSIONLESS_BUNIT,
        [
            "Class map: -1=unknown/non-detection, 0=unclassified, 1=HII, 2=LINER, 3=Seyfert",
            "Uses dust-corrected fluxes; positive classes require stable central +/-1sigma class",
        ],
        dtype=np.int16,
    )
    append_ordered_image(
        ordered_hdul,
        "Halpha_Luminosity_corr",
        HA6562_LUM,
        "erg/s",
        "Includes Balmer-corrected detections plus upper-limit substitutions where Balmer QC fails",
    )
    append_ordered_image(
        ordered_hdul,
        "Halpha_SFR_corr",
        SFR_map,
        "M_sun/yr",
        "Uses Chabrier coefficient; all-spaxel map includes upper-limit substitutions",
    )
    append_ordered_image(
        ordered_hdul,
        "LOGSFR_SURFACE_DENSITY",
        LOG_SFR_surface_density_map,
        CARTA_DIMENSIONLESS_BUNIT,
        "log10 SFR surface density in M_sun/yr/kpc2",
    )
    append_region_pair(
        ordered_hdul,
        "LOGSFR_SURFACE_DENSITY",
        LOG_SFR_surface_density_map_HII,
        LOG_SFR_surface_density_map_SF,
        CARTA_DIMENSIONLESS_BUNIT,
        "log10 SFR surface density in HII regions (M_sun/yr/kpc2)",
        "log10 SFR surface density in SF regions (M_sun/yr/kpc2)",
    )
    append_ordered_image(
        ordered_hdul,
        "LOGSFR_SURFACE_DENSITY_NONHII",
        LOG_SFR_surface_density_map_nonHII,
        CARTA_DIMENSIONLESS_BUNIT,
        "log10 SFR surface density in non-HII regions (M_sun/yr/kpc2)",
    )
    append_ordered_image(
        ordered_hdul,
        "LOGSFR_SURFACE_DENSITY_NONSF",
        LOG_SFR_surface_density_map_nonSF,
        CARTA_DIMENSIONLESS_BUNIT,
        "log10 SFR surface density in non-SF regions (M_sun/yr/kpc2)",
    )
    append_ordered_image(
        ordered_hdul,
        "LOGSFR_SURFACE_DENSITY_UNCLASSIFIED2",
        LOG_SFR_surface_density_map_unclassified2,
        CARTA_DIMENSIONLESS_BUNIT,
        "log10 SFR surface density in Classification-2 unclassified regions",
    )
    append_ordered_image(
        ordered_hdul,
        "LOGSFR_SURFACE_DENSITY_UNCLASSIFIED1",
        LOG_SFR_surface_density_map_unclassified1,
        CARTA_DIMENSIONLESS_BUNIT,
        "log10 SFR surface density in Classification-1 unclassified regions",
    )
    append_ordered_image(
        ordered_hdul,
        "LOGSFR_SURFACE_DENSITY_UPPER",
        LOG_SFR_surface_density_map_upper,
        CARTA_DIMENSIONLESS_BUNIT,
        "log10 SFR surface density upper-limit map",
    )
    append_region_pair(
        ordered_hdul,
        "Halpha_SFR_corr",
        SFR_map_HII,
        SFR_map_SF,
        "M_sun/yr/kpc2",
        "SFR surface density in HII regions (Classification 2)",
        "SFR surface density in SF regions (Classification 1)",
    )
    append_ordered_image(
        ordered_hdul,
        "Halpha_SFR_corr_nonHII",
        SFR_map_nonHII,
        "M_sun/yr/kpc2",
        "SFR surface density in non-HII regions (Classification 2)",
    )
    append_ordered_image(
        ordered_hdul,
        "Halpha_SFR_corr_nonSF",
        SFR_map_nonSF,
        "M_sun/yr/kpc2",
        "SFR surface density in non-SF regions (Classification 1)",
    )
    append_ordered_image(
        ordered_hdul,
        "Halpha_SFR_corr_unclassified2",
        SFR_map_unclassified2,
        "M_sun/yr/kpc2",
        "SFR surface density in Classification-2 unclassified regions",
    )
    append_ordered_image(
        ordered_hdul,
        "Halpha_SFR_corr_unclassified1",
        SFR_map_unclassified1,
        "M_sun/yr/kpc2",
        "SFR surface density in Classification-1 unclassified regions",
    )

    # 4. Strong-line metallicities, always as immediate HII/SF pairs.
    strong_metallicity_specs = [
        ("O_H_D16", O_H_D16_HII, O_H_D16_SF, "Dopita+2016 strong-line oxygen abundance"),
        ("O_H_PG16", O_H_PG16_HII, O_H_PG16_SF, "Pilyugin & Grebel 2016 strong-line oxygen abundance"),
        ("O_H_N2S2_N06", O_H_N2S2_N06_HII, O_H_N2S2_N06_SF, "N2S2-N06 oxygen abundance"),
        ("O_H_O3N2_M13", O_H_O3N2_M13_HII, O_H_O3N2_M13_SF, "O3N2-M13 oxygen abundance"),
        ("O_H_N2_M13", O_H_N2_M13_HII, O_H_N2_M13_SF, "N2-M13 oxygen abundance"),
        ("O_H_O3N2_PP04", O_H_O3N2_PP04_HII, O_H_O3N2_PP04_SF, "O3N2-PP04 oxygen abundance"),
        ("O_H_N2_PP04", O_H_N2_PP04_HII, O_H_N2_PP04_SF, "N2-PP04 oxygen abundance"),
        ("O_H_O3N2_C20", O_H_O3N2_C20_HII, O_H_O3N2_C20_SF, "O3N2-C20 oxygen abundance"),
        ("O_H_O3S2_C20", O_H_O3S2_C20_HII, O_H_O3S2_C20_SF, "O3S2-C20 oxygen abundance"),
        ("O_H_RS32_C20", O_H_RS32_C20_HII, O_H_RS32_C20_SF, "RS32-C20 oxygen abundance"),
        ("O_H_R3_C20", O_H_R3_C20_HII, O_H_R3_C20_SF, "R3-C20 oxygen abundance"),
        ("O_H_N2_C20", O_H_N2_C20_HII, O_H_N2_C20_SF, "N2-C20 oxygen abundance"),
        ("O_H_S2_C20", O_H_S2_C20_HII, O_H_S2_C20_SF, "S2-C20 oxygen abundance"),
        ("O_H_COMBINED_C20", O_H_COMBINED_C20_HII, O_H_COMBINED_C20_SF, "Inverse-variance combined C20 oxygen abundance"),
    ]
    for stem, hii_data, sf_data, description in strong_metallicity_specs:
        append_region_pair(
            ordered_hdul,
            stem,
            hii_data,
            sf_data,
            CARTA_DIMENSIONLESS_BUNIT,
            f"12+log(O/H); {description} in HII regions",
            f"12+log(O/H); {description} in SF regions",
        )

    c20_error_specs = [
        ("O_H_O3N2_C20", O_H_O3N2_C20_HII_ERR, O_H_O3N2_C20_SF_ERR, "O3N2-C20"),
        ("O_H_O3S2_C20", O_H_O3S2_C20_HII_ERR, O_H_O3S2_C20_SF_ERR, "O3S2-C20"),
        ("O_H_RS32_C20", O_H_RS32_C20_HII_ERR, O_H_RS32_C20_SF_ERR, "RS32-C20"),
        ("O_H_R3_C20", O_H_R3_C20_HII_ERR, O_H_R3_C20_SF_ERR, "R3-C20"),
        ("O_H_N2_C20", O_H_N2_C20_HII_ERR, O_H_N2_C20_SF_ERR, "N2-C20"),
        ("O_H_S2_C20", O_H_S2_C20_HII_ERR, O_H_S2_C20_SF_ERR, "S2-C20"),
        ("O_H_COMBINED_C20", O_H_COMBINED_C20_HII_ERR, O_H_COMBINED_C20_SF_ERR, "Combined C20"),
    ]
    for stem, hii_data, sf_data, description in c20_error_specs:
        append_named_hii_sf_pair(
            ordered_hdul,
            f"{stem}_HII_ERR",
            f"{stem}_SF_ERR",
            hii_data,
            sf_data,
            CARTA_DIMENSIONLESS_BUNIT,
            f"1-sigma dex uncertainty for {description} 12+log(O/H) in HII regions",
            f"1-sigma dex uncertainty for {description} 12+log(O/H) in SF regions",
        )

    append_ordered_image(
        ordered_hdul,
        "COMBINED_C20_METHOD",
        combined_c20_method_map,
        CARTA_DIMENSIONLESS_BUNIT,
        "Dominant C20 weight: -1=none, 0=O3N2, 1=O3S2, 2=RS32, 3=R3, 4=N2, 5=S2",
        dtype=np.int16,
    )

    # 5. HII-region photoionisation-grid diagnostics.
    pyqz_references = (
        "Dopita+2013 ApJS 208 10",
        "Vogt+2015 ApJ 799 54 (pyqz)",
    )
    nebulabayes_references = (
        "Thomas+2018 ApJ 856 89 (NebulaBayes)",
    )
    jy22_references = (
        "Ji & Yan 2022 A&A 659 A112",
        "Peng+2026 arXiv:2608.20239",
        "Zenodo 21717332 (released model grid)",
    )
    if ENABLE_PYQZ_METALLICITY_PRODUCTS:
        pyqz_parameter_specs = [
            (
                "O_H_PYQZ",
                "12+log(O/H), pyqz combined multivariate-KDE estimate",
            ),
            (
                "LOG_Q_PYQZ",
                "log10(q / (cm s-1)), pyqz combined multivariate-KDE estimate",
            ),
            (
                "LOG_U_PYQZ",
                "log10 dimensionless ionisation parameter; LogQ-10.47712125472",
            ),
        ]
        for stem, description in pyqz_parameter_specs:
            append_named_hii_sf_pair(
                ordered_hdul,
                f"{stem}_HII",
                f"{stem}_SF",
                MODEL_OUTPUT_MAPS[f"{stem}_HII"],
                MODEL_OUTPUT_MAPS[f"{stem}_SF"],
                CARTA_DIMENSIONLESS_BUNIT,
                f"{description}; existing HII mask",
                f"{description}; existing SF mask",
                references=pyqz_references,
            )
            append_named_hii_sf_pair(
                ordered_hdul,
                f"{stem}_HII_ERR",
                f"{stem}_SF_ERR",
                MODEL_OUTPUT_MAPS[f"{stem}_HII_ERR"],
                MODEL_OUTPUT_MAPS[f"{stem}_SF_ERR"],
                CARTA_DIMENSIONLESS_BUNIT,
                "pyqz uncertainty: maximum half-extent of peak-normalised "
                "0.61 KDE contour; HII mask",
                "pyqz uncertainty: maximum half-extent of peak-normalised "
                "0.61 KDE contour; SF mask",
                references=pyqz_references,
            )
        append_named_hii_sf_pair(
            ordered_hdul,
            "PYQZ_FLAG_HII",
            "PYQZ_FLAG_SF",
            MODEL_OUTPUT_MAPS["PYQZ_FLAG_HII"],
            MODEL_OUTPUT_MAPS["PYQZ_FLAG_SF"],
            CARTA_DIMENSIONLESS_BUNIT,
            "Raw pyqz flag in HII regions; -99 means not evaluated",
            "Raw pyqz flag in SF regions; -99 means not evaluated",
            dtype=np.int16,
            references=pyqz_references,
        )
        append_named_hii_sf_pair(
            ordered_hdul,
            "PYQZ_RS_OFFGRID_HII",
            "PYQZ_RS_OFFGRID_SF",
            MODEL_OUTPUT_MAPS["PYQZ_RS_OFFGRID_HII"],
            MODEL_OUTPUT_MAPS["PYQZ_RS_OFFGRID_SF"],
            "%",
            "Percentage of pyqz Monte Carlo realisations off-grid; HII mask",
            "Percentage of pyqz Monte Carlo realisations off-grid; SF mask",
            references=pyqz_references,
        )
        append_named_hii_sf_pair(
            ordered_hdul,
            "PYQZ_VALID_HII",
            "PYQZ_VALID_SF",
            MODEL_OUTPUT_MAPS["PYQZ_VALID_HII"],
            MODEL_OUTPUT_MAPS["PYQZ_VALID_SF"],
            CARTA_DIMENSIONLESS_BUNIT,
            "1=finite pyqz O/H and LogQ, 0=fit failed, -99=not evaluated; HII",
            "1=finite pyqz O/H and LogQ, 0=fit failed, -99=not evaluated; SF",
            dtype=np.int16,
            references=pyqz_references,
        )

    if ENABLE_NEBULABAYES_METALLICITY_PRODUCTS:
        nb_parameter_specs = [
            (
                "O_H_NEBULABAYES",
                "12+log(O/H), NebulaBayes marginal-posterior peak",
            ),
            (
                "LOG_U_NEBULABAYES",
                "log10 dimensionless ionisation parameter, posterior peak",
            ),
            (
                "LOG_PK_NEBULABAYES",
                "log10(P/k / (K cm-3)), NebulaBayes posterior peak",
            ),
        ]
        for stem, description in nb_parameter_specs:
            append_named_hii_sf_pair(
                ordered_hdul,
                f"{stem}_HII",
                f"{stem}_SF",
                MODEL_OUTPUT_MAPS[f"{stem}_HII"],
                MODEL_OUTPUT_MAPS[f"{stem}_SF"],
                CARTA_DIMENSIONLESS_BUNIT,
                f"{description}; existing HII mask",
                f"{description}; existing SF mask",
                references=nebulabayes_references,
            )
            for bound, bound_description in (
                ("LOW", "lower"),
                ("HIGH", "upper"),
            ):
                append_named_hii_sf_pair(
                    ordered_hdul,
                    f"{stem}_HII_CI68_{bound}",
                    f"{stem}_SF_CI68_{bound}",
                    MODEL_OUTPUT_MAPS[f"{stem}_HII_CI68_{bound}"],
                    MODEL_OUTPUT_MAPS[f"{stem}_SF_CI68_{bound}"],
                    CARTA_DIMENSIONLESS_BUNIT,
                    f"NebulaBayes equal-tailed 68-percent {bound_description} "
                    "credible bound; HII",
                    f"NebulaBayes equal-tailed 68-percent {bound_description} "
                    "credible bound; SF",
                    references=nebulabayes_references,
                )
        append_named_hii_sf_pair(
            ordered_hdul,
            "NB_CHI2_RED_HII",
            "NB_CHI2_RED_SF",
            MODEL_OUTPUT_MAPS["NB_CHI2_RED_HII"],
            MODEL_OUTPUT_MAPS["NB_CHI2_RED_SF"],
            CARTA_DIMENSIONLESS_BUNIT,
            "NebulaBayes best-model reduced chi-square from observational errors; HII",
            "NebulaBayes best-model reduced chi-square from observational errors; SF",
            references=nebulabayes_references,
        )
        append_named_hii_sf_pair(
            ordered_hdul,
            "NB_NLOCALMAX_HII",
            "NB_NLOCALMAX_SF",
            MODEL_OUTPUT_MAPS["NB_NLOCALMAX_HII"],
            MODEL_OUTPUT_MAPS["NB_NLOCALMAX_SF"],
            CARTA_DIMENSIONLESS_BUNIT,
            "Maximum marginal-PDF local-maximum count; -99 not evaluated; HII",
            "Maximum marginal-PDF local-maximum count; -99 not evaluated; SF",
            dtype=np.int16,
            references=nebulabayes_references,
        )
        append_named_hii_sf_pair(
            ordered_hdul,
            "NB_FLAG_HII",
            "NB_FLAG_SF",
            MODEL_OUTPUT_MAPS["NB_FLAG_HII"],
            MODEL_OUTPUT_MAPS["NB_FLAG_SF"],
            CARTA_DIMENSIONLESS_BUNIT,
            [
                "NebulaBayes QC bitmask: 1=O/H edge, 2=logU edge, 4=logP/k edge",
                "8=open/non-finite CI68, 16=fit exception, -99=not evaluated; HII",
            ],
            [
                "NebulaBayes QC bitmask: 1=O/H edge, 2=logU edge, 4=logP/k edge",
                "8=open/non-finite CI68, 16=fit exception, -99=not evaluated; SF",
            ],
            dtype=np.int16,
            references=nebulabayes_references,
        )
        append_named_hii_sf_pair(
            ordered_hdul,
            "NB_VALID_HII",
            "NB_VALID_SF",
            MODEL_OUTPUT_MAPS["NB_VALID_HII"],
            MODEL_OUTPUT_MAPS["NB_VALID_SF"],
            CARTA_DIMENSIONLESS_BUNIT,
            "1=finite O/H, logU, logP/k; 0=fit failed; -99=not evaluated; HII",
            "1=finite O/H, logU, logP/k; 0=fit failed; -99=not evaluated; SF",
            dtype=np.int16,
            references=nebulabayes_references,
        )

    if ENABLE_JY22_METALLICITY_PRODUCTS:
        jy22_parameter_specs = (
            ("O_H_JY22", "12+log(O/H), JY22 posterior mean"),
            ("LOG_U_JY22", "log10 dimensionless ionisation parameter, JY22 posterior mean"),
        )
        for stem, description in jy22_parameter_specs:
            append_named_hii_sf_pair(
                ordered_hdul,
                f"{stem}_HII",
                f"{stem}_SF",
                MODEL_OUTPUT_MAPS[f"{stem}_HII"],
                MODEL_OUTPUT_MAPS[f"{stem}_SF"],
                CARTA_DIMENSIONLESS_BUNIT,
                [
                    f"{description}; existing HII mask",
                    "N2/S2/R3 likelihood uses full propagated ratio covariance",
                ],
                [
                    f"{description}; existing SF mask",
                    "N2/S2/R3 likelihood uses full propagated ratio covariance",
                ],
                references=jy22_references,
            )
            for percentile, percentile_description in (
                ("16", "16th"),
                ("84", "84th"),
            ):
                append_named_hii_sf_pair(
                    ordered_hdul,
                    f"{stem}_HII_{percentile}",
                    f"{stem}_SF_{percentile}",
                    MODEL_OUTPUT_MAPS[f"{stem}_HII_{percentile}"],
                    MODEL_OUTPUT_MAPS[f"{stem}_SF_{percentile}"],
                    CARTA_DIMENSIONLESS_BUNIT,
                    f"JY22 marginal equal-tailed {percentile_description} posterior percentile; HII",
                    f"JY22 marginal equal-tailed {percentile_description} posterior percentile; SF",
                    references=jy22_references,
                )
        append_named_hii_sf_pair(
            ordered_hdul,
            "JY22_CHI2_MIN_HII",
            "JY22_CHI2_MIN_SF",
            MODEL_OUTPUT_MAPS["JY22_CHI2_MIN_HII"],
            MODEL_OUTPUT_MAPS["JY22_CHI2_MIN_SF"],
            CARTA_DIMENSIONLESS_BUNIT,
            "Minimum JY22 chi-square using full N2/S2/R3 covariance; HII",
            "Minimum JY22 chi-square using full N2/S2/R3 covariance; SF",
            references=jy22_references,
        )
        append_named_hii_sf_pair(
            ordered_hdul,
            "JY22_FLAG_HII",
            "JY22_FLAG_SF",
            MODEL_OUTPUT_MAPS["JY22_FLAG_HII"],
            MODEL_OUTPUT_MAPS["JY22_FLAG_SF"],
            CARTA_DIMENSIONLESS_BUNIT,
            [
                "JY22 QC bits: 1=O/H edge, 2=logU edge, 4=posterior invalid",
                "8=covariance invalid, 16=exception, -99=not evaluated; HII",
            ],
            [
                "JY22 QC bits: 1=O/H edge, 2=logU edge, 4=posterior invalid",
                "8=covariance invalid, 16=exception, -99=not evaluated; SF",
            ],
            dtype=np.int16,
            references=jy22_references,
        )
        append_named_hii_sf_pair(
            ordered_hdul,
            "JY22_VALID_HII",
            "JY22_VALID_SF",
            MODEL_OUTPUT_MAPS["JY22_VALID_HII"],
            MODEL_OUTPUT_MAPS["JY22_VALID_SF"],
            CARTA_DIMENSIONLESS_BUNIT,
            "1=all JY22 summaries and chi2 finite; 0=invalid; -99=not evaluated; HII",
            "1=all JY22 summaries and chi2 finite; 0=invalid; -99=not evaluated; SF",
            dtype=np.int16,
            references=jy22_references,
        )

    # 6. PyNeb density, temperature, ionic-abundance, and Te-metallicity products.
    if ENABLE_TE_METALLICITY_PRODUCTS:
        append_ordered_image(
            ordered_hdul,
            "NE_SII",
            NE_SII,
            "cm-3",
            "PyNeb [S II] 6716/6730 density; finite values below 20 cm-3 are clamped",
            refs=True,
        )
        append_ordered_image(
            ordered_hdul,
            "NE_SII_ALL",
            NE_SII_ALL,
            "cm-3",
            "NE_SII where measured; otherwise 20 cm-3 for spatially usable pixels",
            refs=True,
        )
        append_ordered_image(
            ordered_hdul,
            "NE_SII_ERR",
            NE_SII_ERR,
            "cm-3",
            "First-order [S II] density uncertainty from line-ratio error",
            refs=True,
        )

        te_pair_specs = [
            ("TE_NII_HII", "TE_NII_SF", TE_NII_HII, TE_NII_SF, "K", "Te([N II]) from 5755/(6548+6584) using PyNeb"),
            ("TE_NII_HII_ERR", "TE_NII_SF_ERR", TE_NII_HII_ERR, TE_NII_SF_ERR, "K", "First-order Te([N II]) uncertainty"),
            ("TE_SIII_BR24_HII", "TE_SIII_BR24_SF", TE_SIII_BR24_HII, TE_SIII_BR24_SF, "K", "Te([S III]) from 6312/(9069+9532); 9532=2.5*9069"),
            ("TE_SIII_BR24_HII_ERR", "TE_SIII_BR24_SF_ERR", TE_SIII_BR24_HII_ERR, TE_SIII_BR24_SF_ERR, "K", "First-order Te([S III]) uncertainty"),
            ("TE_OIII_BR24_HII", "TE_OIII_BR24_SF", TE_OIII_BR24_HII, TE_OIII_BR24_SF, "K", "Brazzini+2024 Te([O III]) inferred from Te([S III])"),
            ("TE_OIII_BR24_HII_ERR", "TE_OIII_BR24_SF_ERR", TE_OIII_BR24_HII_ERR, TE_OIII_BR24_SF_ERR, "K", "Brazzini+2024 Te([O III]) uncertainty"),
            ("TE_OIII_NII_CHAIN_HII", "TE_OIII_NII_CHAIN_SF", TE_OIII_NII_CHAIN_HII, TE_OIII_NII_CHAIN_SF, "K", "Te([O III]) inferred from Te([N II]) via Brazzini+2024 Te chains"),
            ("TE_OIII_NII_CHAIN_HII_ERR", "TE_OIII_NII_CHAIN_SF_ERR", TE_OIII_NII_CHAIN_HII_ERR, TE_OIII_NII_CHAIN_SF_ERR, "K", "Te([O III]) from Te([N II]) uncertainty"),
            ("O_PLUS_BR24_HII", "O_PLUS_BR24_SF", O_PLUS_BR24_HII, O_PLUS_BR24_SF, CARTA_DIMENSIONLESS_BUNIT, "O+/H+ from [O II] 7325 in the Brazzini-style method"),
            ("O_PLUS_BR24_HII_ERR", "O_PLUS_BR24_SF_ERR", O_PLUS_BR24_HII_ERR, O_PLUS_BR24_SF_ERR, CARTA_DIMENSIONLESS_BUNIT, "1-sigma O+/H+ uncertainty for the Brazzini-style method"),
            ("O_DPLUS_BR24_HII", "O_DPLUS_BR24_SF", O_DOUBLEPLUS_BR24_HII, O_DOUBLEPLUS_BR24_SF, CARTA_DIMENSIONLESS_BUNIT, "O++/H+ from [O III] 4959+5007 in the Brazzini-style method"),
            ("O_DPLUS_BR24_HII_ERR", "O_DPLUS_BR24_SF_ERR", O_DOUBLEPLUS_BR24_HII_ERR, O_DOUBLEPLUS_BR24_SF_ERR, CARTA_DIMENSIONLESS_BUNIT, "1-sigma O++/H+ uncertainty for the Brazzini-style method"),
            ("O_H_BR24_DIRECT_HII", "O_H_BR24_DIRECT_SF", O_H_BR24_DIRECT_HII, O_H_BR24_DIRECT_SF, CARTA_DIMENSIONLESS_BUNIT, "12+log(O/H) from Brazzini+2024 multi-zone direct-style method"),
            ("O_H_BR24_DIRECT_HII_ERR", "O_H_BR24_DIRECT_SF_ERR", O_H_BR24_DIRECT_HII_ERR, O_H_BR24_DIRECT_SF_ERR, CARTA_DIMENSIONLESS_BUNIT, "1-sigma dex uncertainty for Brazzini+2024 12+log(O/H)"),
            ("O_PLUS_NII_OII7325_HII", "O_PLUS_NII_OII7325_SF", O_PLUS_NII_OII7325_HII, O_PLUS_NII_OII7325_SF, CARTA_DIMENSIONLESS_BUNIT, "O+/H+ from [O II] 7325 using Te([N II])"),
            ("O_PLUS_NII_OII7325_HII_ERR", "O_PLUS_NII_OII7325_SF_ERR", O_PLUS_NII_OII7325_HII_ERR, O_PLUS_NII_OII7325_SF_ERR, CARTA_DIMENSIONLESS_BUNIT, "1-sigma O+/H+ uncertainty for the NII+OII7325 method"),
            ("O_DPLUS_NII_OII7325_HII", "O_DPLUS_NII_OII7325_SF", O_DOUBLEPLUS_NII_OII7325_HII, O_DOUBLEPLUS_NII_OII7325_SF, CARTA_DIMENSIONLESS_BUNIT, "O++/H+ from [O III] with Te([O III]) inferred from Te([N II])"),
            ("O_DPLUS_NII_OII7325_HII_ERR", "O_DPLUS_NII_OII7325_SF_ERR", O_DOUBLEPLUS_NII_OII7325_HII_ERR, O_DOUBLEPLUS_NII_OII7325_SF_ERR, CARTA_DIMENSIONLESS_BUNIT, "1-sigma O++/H+ uncertainty for the NII+OII7325 method"),
            ("O_H_NII_OII7325_HII", "O_H_NII_OII7325_SF", O_H_NII_OII7325_HII, O_H_NII_OII7325_SF, CARTA_DIMENSIONLESS_BUNIT, "12+log(O/H) from MAUVE NII+OII7325 semi-direct ionic sum"),
            ("O_H_NII_OII7325_HII_ERR", "O_H_NII_OII7325_SF_ERR", O_H_NII_OII7325_HII_ERR, O_H_NII_OII7325_SF_ERR, CARTA_DIMENSIONLESS_BUNIT, "1-sigma dex uncertainty for NII+OII7325 12+log(O/H)"),
            ("O_H_NII_K25_HII", "O_H_NII_K25_SF", O_H_NII_K25_HII, O_H_NII_K25_SF, CARTA_DIMENSIONLESS_BUNIT, "12+log(O/H) from Kreckel+2025/Mendez-Delgado+2023 Te([N II]) proxy"),
            ("O_H_NII_K25_HII_ERR", "O_H_NII_K25_SF_ERR", O_H_NII_K25_HII_ERR, O_H_NII_K25_SF_ERR, CARTA_DIMENSIONLESS_BUNIT, "1-sigma dex uncertainty for Kreckel+2025/Mendez-Delgado+2023 12+log(O/H)"),
        ]
        for hii_name, sf_name, hii_data, sf_data, bunit, description in te_pair_specs:
            append_named_hii_sf_pair(
                ordered_hdul,
                hii_name,
                sf_name,
                hii_data,
                sf_data,
                bunit,
                f"{description} in HII regions",
                f"{description} in SF regions",
                refs=True,
            )

        append_ordered_image(
            ordered_hdul,
            "NE_SII_FIXED20",
            NE_SII_FIXED20,
            CARTA_DIMENSIONLESS_BUNIT,
            "1 where NE_SII was fixed to 20 cm-3",
            dtype=np.int16,
            refs=True,
        )
        append_ordered_image(
            ordered_hdul,
            "NE_SII_EXACT_HIDEN",
            NE_SII_EXACT_HIGH_DENSITY,
            CARTA_DIMENSIONLESS_BUNIT,
            "1 where NE_SII used exact PyNeb getTemDen high-density fallback",
            dtype=np.int16,
            refs=True,
        )
        te_mask_pairs = [
            ("BR24_DIRECT_VALID", BR24_DIRECT_VALID_HII, BR24_DIRECT_VALID_SF, "all Brazzini+2024 direct-style inputs passed detection gates"),
            ("NII_OII7325_VALID", NII_OII7325_VALID_HII, NII_OII7325_VALID_SF, "all NII+OII7325 semi-direct inputs passed detection gates"),
            ("NII_K25_VALID", NII_K25_VALID_HII, NII_K25_VALID_SF, "all Kreckel+2025 NII-only inputs passed detection gates"),
        ]
        for stem, hii_data, sf_data, description in te_mask_pairs:
            append_region_pair(
                ordered_hdul,
                stem,
                hii_data,
                sf_data,
                CARTA_DIMENSIONLESS_BUNIT,
                f"1 where HII {description}",
                f"1 where SF {description}",
                dtype=np.int16,
                refs=True,
            )

    extension_names = [hdu.name for hdu in ordered_hdul[1:]]
    duplicate_names = sorted(
        name for name in set(extension_names) if extension_names.count(name) > 1
    )
    if duplicate_names:
        raise RuntimeError(
            "Ordered SFR+Z FITS schema contains duplicate EXTNAME values: "
            + ", ".join(duplicate_names)
        )
    if (
        ENABLE_PYQZ_METALLICITY_PRODUCTS
        or ENABLE_NEBULABAYES_METALLICITY_PRODUCTS
        or ENABLE_JY22_METALLICITY_PRODUCTS
    ):
        expected_model_names = [
            name
            for name in ordered_model_hdu_names()
            if (
                (ENABLE_PYQZ_METALLICITY_PRODUCTS and "PYQZ" in name)
                or (
                    ENABLE_NEBULABAYES_METALLICITY_PRODUCTS
                    and ("NEBULABAYES" in name or name.startswith("NB_"))
                )
                or (ENABLE_JY22_METALLICITY_PRODUCTS and "JY22" in name)
            )
        ]
    else:
        expected_model_names = []
    if expected_model_names:
        model_start = extension_names.index("COMBINED_C20_METHOD") + 1
        actual_model_names = extension_names[
            model_start : model_start + len(expected_model_names)
        ]
        if actual_model_names != expected_model_names:
            raise RuntimeError(
                "Ordered model-grid FITS schema differs from its tested contract."
            )
        integer_model_names = {
            name
            for name in expected_model_names
            if "FLAG" in name or "VALID" in name or "NLOCALMAX" in name
        }
        for name in expected_model_names:
            hdu = ordered_hdul[name]
            if hdu.data.shape != HA6562_FLUX.shape:
                raise RuntimeError(
                    f"Model-grid HDU {name} shape {hdu.data.shape} differs from "
                    f"gas-map shape {HA6562_FLUX.shape}."
                )
            expected_dtype = np.dtype(
                np.int16 if name in integer_model_names else np.float64
            )
            if (
                hdu.data.dtype.kind != expected_dtype.kind
                or hdu.data.dtype.itemsize != expected_dtype.itemsize
            ):
                raise RuntimeError(
                    f"Model-grid HDU {name} has dtype {hdu.data.dtype}, expected "
                    f"{expected_dtype}."
                )
            expected_bunit = (
                "%" if "RS_OFFGRID" in name else CARTA_DIMENSIONLESS_BUNIT
            )
            if hdu.header.get("BUNIT") != expected_bunit:
                raise RuntimeError(
                    f"Model-grid HDU {name} has BUNIT={hdu.header.get('BUNIT')!r}, "
                    f"expected {expected_bunit!r}."
                )
            if "COMMENT" not in hdu.header:
                raise RuntimeError(f"Model-grid HDU {name} has no science comment.")
            for wcs_key in (
                "CTYPE1", "CTYPE2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2",
                "CDELT1", "CDELT2", "CD1_1", "CD1_2", "CD2_1", "CD2_2",
            ):
                if (
                    wcs_key in gas_header
                    and hdu.header.get(wcs_key) != gas_header.get(wcs_key)
                ):
                    raise RuntimeError(
                        f"Model-grid HDU {name} changed WCS card {wcs_key}."
                    )
    return ordered_hdul


new_hdul = build_ordered_output_hdul(new_hdul[:ORIGINAL_GAS_HDU_COUNT])
new_hdul.writeto(out_path, overwrite=True)
print("Further file written ➜", out_path.resolve())

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
print("Number of pixels with Halpha detected, Hbeta detected, NII not detected, OIII not detected and unclassified1:", 
      np.sum(HA_detected_HB_detected_NII_not_detected_OIII_not_detected & mask_N2_unclassified1))
print("Number of pixels with Halpha detected, Hbeta detected, NII not detected, OIII detected and unclassified1:",
      np.sum(HA_detected_HB_detected_NII_not_detected_OIII_detected & mask_N2_unclassified1))
print("Number of pixels with Halpha detected, Hbeta detected, NII detected, OIII not detected and unclassified1:",
      np.sum(HA_detected_HB_detected_NII_detected_OIII_not_detected & mask_N2_unclassified1))
print("Number of pixels with Halpha detected, Hbeta detected, NII detected, OIII detected and unclassified1:",
      np.sum(HA_detected_HB_detected_NII_detected_OIII_detected & mask_N2_unclassified1))
print("Number of pixels with Halpha detected, Hbeta detected, SII not detected, OIII not detected and unclassified1:", 
      np.sum(HA_detected_HB_detected_SII_not_detected_OIII_not_detected & mask_N2_unclassified1))
print("Number of pixels with Halpha detected, Hbeta detected, SII not detected, OIII detected and unclassified1:",
      np.sum(HA_detected_HB_detected_SII_not_detected_OIII_detected & mask_N2_unclassified1))
print("Number of pixels with Halpha detected, Hbeta detected, SII detected, OIII not detected and unclassified1:",
      np.sum(HA_detected_HB_detected_SII_detected_OIII_not_detected & mask_N2_unclassified1))
print("Number of pixels with Halpha detected, Hbeta detected, SII detected, OIII detected and unclassified1:",
      np.sum(HA_detected_HB_detected_SII_detected_OIII_detected & mask_N2_unclassified1))
print("--------------------------------------------------------------")
print(f"Total Halpha luminosity map sum (detections + upper-limit substitutions): {np.nansum(HA6562_LUM):.2e} erg/s")
print(f"Total corrected Halpha luminosity from SF region: {np.nansum(HA6562_LUM[mask_SF]):.2e} erg/s")
print(f"Total Halpha SFR map sum (detections + upper-limit substitutions): {np.nansum(SFR_map):.2f} M☉/yr or in log10 scale: {np.log10(np.nansum(SFR_map)):.2f} log(M☉/yr)")
print(f"Total Halpha SFR from SF region: {np.nansum(SFR_map[mask_SF]):.2f} M☉/yr or in log10 scale: {np.log10(np.nansum(SFR_map[mask_SF])):.2f} log(M☉/yr)")
print(f"Total Halpha SFR from HII region: {np.nansum(SFR_map[mask_HII]):.2f} M☉/yr or in log10 scale: {np.log10(np.nansum(SFR_map[mask_HII])):.2f} log(M☉/yr)")
print("--------------------------------------------------------------")
print("[O/H] D16 SF: Total metallicity in SF region: ", O_H_D16_SF_total)
print("[O/H] PG16 SF: Total metallicity in SF region: ", O_H_PG16_SF_total)
print("[O/H] N2S2-N06 SF: Total metallicity in SF region: ", O_H_N2S2_N06_SF_total)
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
print("[O/H] Combined-C20 SF: Total metallicity in SF region: ", O_H_COMBINED_C20_SF_total,
      f"(weighted {O_H_COMBINED_C20_SF_total_NMETHODS} methods; dominant {c20_method_label(O_H_COMBINED_C20_SF_total_METHOD)})")
print("--------------------------------------------------------------")
print("[O/H] D16 HII: Total metallicity in HII region: ", O_H_D16_HII_total)
print("[O/H] PG16 HII: Total metallicity in HII region: ", O_H_PG16_HII_total)
print("[O/H] N2S2-N06 HII: Total metallicity in HII region: ", O_H_N2S2_N06_HII_total)
print("[O/H] O3N2-M13 HII: Total metallicity in HII region: ", O_H_O3N2_M13_HII_total)
print("[O/H] N2-M13 HII: Total metallicity in HII region: ", O_H_N2_M13_HII_total)
print("[O/H] O3N2-PP04 HII: Total metallicity in HII region: ", O_H_O3N2_PP04_HII_total)
print("[O/H] N2-PP04 HII: Total metallicity in HII region: ", O_H_N2_PP04_HII_total)
print("[O/H] O3N2-C20 HII: Total metallicity in HII region: ", O_H_O3N2_C20_HII_total)
print("[O/H] O3S2-C20 HII: Total metallicity in HII region: ", O_H_O3S2_C20_HII_total)
print("[O/H] RS32-C20 HII: Total metallicity in HII region: ", O_H_RS32_C20_HII_total)
print("[O/H] R3-C20 HII: Total metallicity in HII region: ", O_H_R3_C20_HII_total)
print("[O/H] N2-C20 HII: Total metallicity in HII region: ", O_H_N2_C20_HII_total)
print("[O/H] S2-C20 HII: Total metallicity in HII region: ", O_H_S2_C20_HII_total)
print("[O/H] Combined-C20 HII: Total metallicity in HII region: ", O_H_COMBINED_C20_HII_total,
      f"(weighted {O_H_COMBINED_C20_HII_total_NMETHODS} methods; dominant {c20_method_label(O_H_COMBINED_C20_HII_total_METHOD)})")
print("--------------------------------------------------------------")
print("[O/H] D16: Total metallicity in total region: ", O_H_D16_total)
print("[O/H] PG16: Total metallicity in total region: ", O_H_PG16_total)
print("[O/H] N2S2-N06: Total metallicity in total region: ", O_H_N2S2_N06_total)
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
print("[O/H] Combined-C20: Total metallicity in total region: ", O_H_COMBINED_C20_total,
      f"(weighted {O_H_COMBINED_C20_total_NMETHODS} methods; dominant {c20_method_label(O_H_COMBINED_C20_total_METHOD)})")
print("--------------------------------------------------------------")


print(f"Total runtime: {time.perf_counter() - t0:.1f} s")

# End of file
