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
* Added Milky Way (CCM89; Cardelli, Clayton & Mathis 1989) extinction curve as k(λ)=A(λ)/E(B−V).
* Set Milky Way (CCM89, R_V=3.1) as the default extinction law for k-values (Balmer-decrement dust correction).
* Kept Calzetti (2000) curve intact and available via `extinction_k(..., law="calzetti")`.

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
# Supported: "mw" (CCM89 Milky Way; default), "calzetti" (Calzetti 2000)
extinction_law = "mw"
mw_rv = 3.1

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


def _candidate_paths(*paths: Path | None) -> list[Path]:
    candidates: list[Path] = []
    seen: set[str] = set()

    for path in paths:
        if path is None:
            continue

        expanded = path.expanduser()
        if _has_glob_chars(expanded):
            path_candidates = sorted(
                expanded.parent.glob(expanded.name), key=lambda candidate: str(candidate)
            )
        else:
            path_candidates = [expanded]

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
        if remainder == suffix:
            return ""
        marker = f"_{suffix}"
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
                    for candidate in sorted(search_dir.glob(pattern), key=lambda path: str(path)):
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


def build_extended_output_path(input_path: Path) -> Path:
    suffix = input_path.suffix or ".fits"
    stem = input_path.name[: -len(input_path.suffix)] if input_path.suffix else input_path.name
    if stem.endswith("_extended"):
        return Path(f"{stem}{suffix}")
    return Path(f"{stem}_extended{suffix}")

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
    "gas_BIN_maps_extended.fits",
    "gas_bin_maps_extended.fits",
    "BIN_maps_extended.fits",
    "bin_maps_extended.fits",
]
kin_suffixes = [
    "KIN_maps_extended.fits",
    "kin_maps_extended.fits",
]
bin_extended_suffixes = [
    "SPATIAL_BINNING_maps_extended.fits",
    "spatial_binning_maps_extended.fits",
]

product_dir = Path("products/v0.6") / gal
input_keyword, product_paths = resolve_keyworded_input_group(
    "SFR product FITS inputs",
    root,
    fallback_root,
    product_dir,
    gal,
    {
        "gas-line FITS": gas_suffixes,
        "extended binning FITS": bin_extended_suffixes,
    },
)
src = product_paths["gas-line FITS"]
out_path = build_extended_output_path(src)
bin_extended_path = product_paths["extended binning FITS"]
kin_path = find_keyworded_input_path(
    root,
    fallback_root,
    product_dir,
    gal,
    kin_suffixes,
    input_keyword,
)
inclination_path = Path("MAUVE_Inclination.dat")

# two key cut parameters
cut = 3 # FLUX/FLUX_ERR
noise = 20 # detection limit of FLUX, in the unit of 10^-20 erg/s

# ------------------------------------------------------------------
# 1.  Load gas-line maps (extended version file if it already exists)
# ------------------------------------------------------------------

print(f"Shared input keyword ➜ {_format_keyword(input_keyword)}")
print(f"Reading gas-line FITS ➜ {src}")
print(f"Matched binning FITS ➜ {bin_extended_path}")
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
    gas_header = hdul['HA6562_FLUX'].header.copy()
    hdul.close()

gas_header

# ------------------------------------------------------------------
# 2.  Foreground-star removal  (verbatim from notebook snippet)
# ------------------------------------------------------------------

if kin_path is not None:
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
# 2. Convert BD to gas E(B-V) using the selected extinction curve (default: Milky Way CCM89)
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

k_HB4861  = extinction_k(0.4861)   # CCM89(MW, Rv=3.1) ≈ 3.609
k_HA6562  = extinction_k(0.6562)   # CCM89(MW, Rv=3.1) ≈ 2.535
k_OIII5006= extinction_k(0.5006)
k_NII6583 = extinction_k(0.6583)
k_SII6716 = extinction_k(0.6716)
k_SII6730 = extinction_k(0.6730)

# Print k(λ) values used for dust correction (this will also appear in redirected *.log outputs)
print("--------------------------------------------------------------")
print(f"Extinction law for k(λ)=A(λ)/E(B−V): {extinction_law} (mw_rv={mw_rv})")
print(f"k(Hβ 4861Å)  = {float(k_HB4861):.4f}   at λ=0.4861 µm")
print(f"k(Hα 6562Å)  = {float(k_HA6562):.4f}   at λ=0.6562 µm")
print(f"k([OIII]5006)= {float(k_OIII5006):.4f}   at λ=0.5006 µm")
print(f"k([NII]6583) = {float(k_NII6583):.4f}   at λ=0.6583 µm")
print(f"k([SII]6716) = {float(k_SII6716):.4f}   at λ=0.6716 µm")
print(f"k([SII]6730) = {float(k_SII6730):.4f}   at λ=0.6730 µm")
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

# Correct the fluxes with E(B-V)_BD
HB4861_FLUX_corr = correct_flux_with_ebv(HB4861_FLUX, E_BV_BD, k_HB4861)
HA6562_FLUX_corr = correct_flux_with_ebv(HA6562_FLUX, E_BV_BD, k_HA6562)
OIII5006_FLUX_corr = correct_flux_with_ebv(OIII5006_FLUX, E_BV_BD, k_OIII5006)
NII6583_FLUX_corr = correct_flux_with_ebv(NII6583_FLUX, E_BV_BD, k_NII6583)
SII6716_FLUX_corr = correct_flux_with_ebv(SII6716_FLUX, E_BV_BD, k_SII6716)
SII6730_FLUX_corr = correct_flux_with_ebv(SII6730_FLUX, E_BV_BD, k_SII6730)
HB4861_FLUX_ERR_corr = correct_flux_error_with_ebv(HB4861_FLUX, HB4861_FLUX_ERR, E_BV_BD, k_HB4861, E_BV_BD_ERR)
HA6562_FLUX_ERR_corr = correct_flux_error_with_ebv(HA6562_FLUX, HA6562_FLUX_ERR, E_BV_BD, k_HA6562, E_BV_BD_ERR)
OIII5006_FLUX_ERR_corr = correct_flux_error_with_ebv(OIII5006_FLUX, OIII5006_FLUX_ERR, E_BV_BD, k_OIII5006, E_BV_BD_ERR)
NII6583_FLUX_ERR_corr = correct_flux_error_with_ebv(NII6583_FLUX, NII6583_FLUX_ERR, E_BV_BD, k_NII6583, E_BV_BD_ERR)
SII6716_FLUX_ERR_corr = correct_flux_error_with_ebv(SII6716_FLUX, SII6716_FLUX_ERR, E_BV_BD, k_SII6716, E_BV_BD_ERR)
SII6730_FLUX_ERR_corr = correct_flux_error_with_ebv(SII6730_FLUX, SII6730_FLUX_ERR, E_BV_BD, k_SII6730, E_BV_BD_ERR)

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

# ------------------------------------------------------------------
# 10.  Calculate the total Metallicity in SF regions (Classification 1)
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
