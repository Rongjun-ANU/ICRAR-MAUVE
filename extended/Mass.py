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

Changes (2025-11-26)
-----------------------
* Updated solar absolute magnitude values from Willmer (2018, ApJS, 236, 47):
  - Bessell R-band: M_r_sun = 4.61 (AB magnitude)
  - SDSS r-band: M_r_sun = 4.65 (AB magnitude)
* Reference: https://iopscience.iop.org/article/10.3847/1538-4365/aabfdf (Table 3)

Changes (2025-12-02)
-----------------------
* Switched extinction correction to use the Calzetti (2000) attenuation curve.
* Added `calzetti_k(w_um)` helper function to calculate k(λ).
* Updated internal attenuation calculation: A_r = k_r_calz * EBV_star.
* Effective wavelengths used: Bessell R ~ 0.64 µm, SDSS r ~ 0.623 µm.

Changes (2026-03-31)
-----------------------
* Added `--fallback-root` for secondary input lookup.
* Datacube, binning, SFH, and weights inputs are now resolved independently,
  with the primary root checked before the fallback root.
* This allows mixed-root runs where different required inputs come from CANFAR
  and local storage in the same execution.

Changes (2026-04-03)
-----------------------
* Added first-order uncertainty propagation for newly derived mass products.
* Implemented uncertainty maps for:
    - FLUX_R_corr_ERR
    - magnitude_r_ERR
    - magnitude_r_uncorrected_ERR
    - LOGMSTAR_ERR
    - LOGMASS_SURFACE_DENSITY_ERR
* R-band photometric uncertainties are propagated from cube STAT variance when
    available; otherwise a bin-wise scatter fallback is used.
* EBV uncertainty from SFH maps is included when an EBV error extension exists.
* Clarified strict photometric propagation:
        - AB maggies are represented as a linear functional of the spectrum using
            coefficients extracted from the same speclite filter operator used for
            magnitude calculations.
        - STAT variance is propagated with the exact discrete first-order Jacobian
            of that linear operator, assuming independent spectral bins (no covariance).
        - Summary equations (per spaxel, on the discrete wavelength grid):
            maggies = sum_i (a_i * f_i)
            Var(maggies) = sum_i (a_i^2 * Var(f_i))
            sigma_mag = (2.5/ln 10) * sigma_maggies / maggies
* Added performance controls:
        - `--ncpus` to control BLAS/OpenMP thread count.
        - `--row-block-size` for chunked vectorized row processing.
* Performance/logging refinement for large cubes:
    - The vectorized photometric pass now operates only on wavelength bins
        where the filter linear coefficients are non-zero.
    - Added more frequent elapsed-time progress prints during row-block processing.

Changes (2026-04-23)
-----------------------
* Fixed the composite `M/L_R` calculation to normalize SFH template weights per bin
  before combining SSP properties, and added an explicit interpretation switch for
  how nGIST `SFH` weights are converted into an R-band `M/L`.
* Replaced the approximate `V-R`-based LIGHT-mode conversion with an exact
  normalization-window conversion using the SSP template spectra themselves to
  measure the pre-normalization mean flux within the fitted wavelength range
  (`LMIN`–`LMAX`), while keeping the R-band luminosity term from the BaSTI
  photometry table.
* Added an optional `--template-library-dir` override and an automatic fallback
  from the exact `light_norm_to_r` mode to the approximate `light_v_to_r` mode
  when the original SSP template library is not available on the machine.
* Added diagnostics for raw SFH weight sums and for the resulting candidate `M/L_R`
  distributions (simple-light, exact light norm→R, approximate light V→R, and
  mass-harmonic), so different runs can be compared directly in the logs.
* Relaxed the R-band photometry validity criterion from "all filter-support wavelength
  bins must be finite" to a support-fraction test with renormalization for small masked
  gaps inside the filter bandpass.
* Added diagnostics for strict-vs-relaxed filter-support acceptance counts and the
  retained filter-support fraction across spaxels.
* Separated cube-photometry validity from the nGIST binning-product `FLUX` mask so
  integrated R-band totals and reconstructed binned maps are not biased by binning-map
  footprint differences.
* Clarified integrated summary logging by printing both uncorrected and extinction-
  corrected R-band totals, plus the integrated stellar `M/L_R`.

Changes (2026-04-26)
-----------------------
* Added flexible input filename resolution for datacube version suffixes and
  upper/lower-case MAUVE product names while preserving normalized `_extended`
  output naming.
* Fixed the `M/L_R` BINID lookup to support weights tables with more rows than
  the currently populated BINID range.
* Clarified stellar-mass surface-density comments to match the stored
  `log(Msol/kpc2)` unit.

Changes (2026-04-27)
-----------------------
* Non-cube product inputs are now resolved as a matched optional-keyword group,
  e.g. `{galaxy}_{keyword}_SPATIAL_BINNING_maps.fits`,
  `{galaxy}_{keyword}_SFH_maps.fits`, and
  `{galaxy}_{keyword}_sfh_weights.fits`.

"""

# ------------------------------------------------------------------
# User Configuration Parameters
# ------------------------------------------------------------------

# Inclination correction toggle
# Set to True to apply cos(θ) inclination correction, False to disable
apply_inclination_correction = True

# Interpretation of nGIST SFH template weights when converting them into a
# composite R-band mass-to-light ratio.
# Available options:
#   "light_norm_to_r": exact LIGHT-mode conversion using the SSP template mean
#                      flux in the fitted wavelength window and the BaSTI
#                      R-band luminosity proxy (default).
#   "light_r":      weights are already R-band light fractions.
#   "light_v_to_r": approximate LIGHT-mode fallback using SSP V-R colours.
#   "mass":         weights are mass fractions, giving the harmonic-mean M/L_R.
sfh_weight_interpretation = "light_norm_to_r"

# If the exact LIGHT-mode conversion is requested but the original SSP template
# library is unavailable, fall back to this approximate interpretation.
sfh_weight_interpretation_fallback = "light_v_to_r"

# Minimum retained fraction of the R-band filter support required to accept
# a spaxel photometry measurement when some wavelength bins are masked/NaN.
min_filter_support_fraction = 0.99

# ------------------------------------------------------------------
# 0.  Command‑line interface
# ------------------------------------------------------------------
import argparse
import os
import re
import sys
from pathlib import Path

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


def resolve_existing_path(label: str, *paths: Path | None) -> Path:
    candidates = _candidate_paths(*paths)
    for candidate in candidates:
        if candidate.exists():
            return candidate

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
    p.add_argument(
        "--fallback-root",
        default=".",
        help="Fallback directory searched when a required input is not found under --root",
    )
    p.add_argument(
        "--disable-stat-propagation",
        action="store_true",
        help="Disable STAT-based photometric uncertainty propagation and use bin-scatter fallback only",
    )
    p.add_argument(
        "--ncpus",
        type=int,
        default=0,
        help="Number of BLAS/OpenMP threads (0 = library default)",
    )
    p.add_argument(
        "--row-block-size",
        type=int,
        default=16,
        help="Number of image rows processed per vectorized block (default: 16)",
    )
    p.add_argument(
        "--template-library-dir",
        default=None,
        help=(
            "Optional path to the SSP template library used by the SFH fit "
            "(needed for the exact light_norm_to_r M/L conversion)"
        ),
    )
    return p.parse_args()

args          = parse_args()
galaxy        = args.galaxy.upper()   # ensure consistent capitalisation
rootdir       = Path(args.root).expanduser().resolve()
fallback_root = (
    Path(args.fallback_root).expanduser().resolve()
    if args.fallback_root is not None
    else None
)
disable_stat_propagation = args.disable_stat_propagation
ncpus = max(int(args.ncpus), 0)
row_block_size = max(int(args.row_block_size), 1)
template_library_dir_override = (
    Path(args.template_library_dir).expanduser().resolve()
    if args.template_library_dir is not None
    else None
)

# Ensure line-buffered logs so progress appears promptly when stdout is piped.
try:
    sys.stdout.reconfigure(line_buffering=True)
    sys.stderr.reconfigure(line_buffering=True)
except Exception:
    pass

if ncpus > 0:
    for env_var in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "BLIS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[env_var] = str(ncpus)

# ------------------------------------------------------------------
# 1.  File paths derived from CLI args
# ------------------------------------------------------------------
cube_names = [
    f"{galaxy}_DATACUBE_FINAL_WCS_Pall_mad_red_v3.fits",
    f"{galaxy}_DATACUBE_FINAL_WCS_Pall_mad_red_v*.fits",
]
bin_suffixes = [
    "SPATIAL_BINNING_maps.fits",
    "spatial_binning_maps.fits",
]
sfh_suffixes = [
    "SFH_maps.fits",
    "sfh_maps.fits",
]
weight_suffixes = [
    "sfh-weights.fits",
    "sfh_weights.fits",
]

cube_path = resolve_existing_path(
    "datacube FITS",
    *build_named_input_candidates(rootdir, Path("cubes/v3.0"), cube_names),
    *build_named_input_candidates(fallback_root, Path("cubes/v3.0"), cube_names),
)

input_keyword, product_paths = resolve_keyworded_input_group(
    "stellar-population product FITS inputs",
    rootdir,
    fallback_root,
    Path("products/v0.6") / galaxy,
    galaxy,
    {
        "spatial binning FITS": bin_suffixes,
        "SFH FITS": sfh_suffixes,
        "SFH weights FITS": weight_suffixes,
    },
)
bin_path = product_paths["spatial binning FITS"]
sfh_path = product_paths["SFH FITS"]
weight_path = product_paths["SFH weights FITS"]
script_dir = Path(__file__).resolve().parent
phot_path = resolve_existing_path(
    "photometry table",
    Path("BaSTI+Chabrier.dat"),
    script_dir / "BaSTI+Chabrier.dat",
    script_dir.parent / "data" / "IC3392" / "BaSTI+Chabrier.dat",
)
out_path = build_extended_output_path(bin_path)   # output in CWD
inclination_path = resolve_existing_path(
    "inclination table",
    Path("MAUVE_Inclination.dat"),
    script_dir / "MAUVE_Inclination.dat",
)

# For backwards compatibility with variable names in the original notebook
vor_path     = bin_path
binning_path = bin_path

print("\n=== Using the following files ===")
print("Primary root :", rootdir)
if fallback_root is not None:
    print("Fallback root:", fallback_root)
print("Input keyword:", _format_keyword(input_keyword))
print("Cube         :", cube_path)
print("Binning map  :", bin_path)
print("SFH map      :", sfh_path)
print("Weights      :", weight_path)
print("Photometry   :", phot_path)
if template_library_dir_override is not None:
    print("Template lib :", template_library_dir_override)
print("Output       :", out_path, "\n")
if disable_stat_propagation:
    print("STAT propagation: DISABLED by --disable-stat-propagation")
if ncpus > 0:
    print(f"Thread control: forcing {ncpus} BLAS/OpenMP threads")
else:
    print("Thread control: using library default thread settings")
print(f"Vectorized processing row-block size: {row_block_size}")

# ------------------------------------------------------------------
# 2.  Imports (same as original)
# ------------------------------------------------------------------
import glob, warnings, gc, time
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


def trapz_linear_weights(x):
    """Return coefficients c such that trapz(y, x) == sum(c * y)."""
    x = np.asarray(x, dtype=np.float64)
    if x.ndim != 1 or x.size < 2:
        raise ValueError("x must be a 1D array with at least two points")

    dx = np.diff(x)
    ctrap = np.empty_like(x, dtype=np.float64)
    ctrap[0] = 0.5 * dx[0]
    ctrap[-1] = 0.5 * dx[-1]
    if x.size > 2:
        ctrap[1:-1] = 0.5 * (dx[:-1] + dx[1:])
    return ctrap


def quadrature_sum(*terms):
    """Quadrature-sum independent 1-sigma terms while tolerating NaNs."""
    if len(terms) == 0:
        raise ValueError("quadrature_sum requires at least one term")

    total = np.zeros_like(np.asarray(terms[0], dtype=np.float64), dtype=np.float64)
    any_finite = np.zeros_like(total, dtype=bool)

    for term in terms:
        arr = np.asarray(term, dtype=np.float64)
        finite = np.isfinite(arr)
        total[finite] += arr[finite] ** 2
        any_finite |= finite

    out = np.sqrt(total)
    out[~any_finite] = np.nan
    return out


def finite_median(arr):
    """Return nanmedian if finite values exist, else NaN without warning."""
    arr = np.asarray(arr)
    finite = np.isfinite(arr)
    if not np.any(finite):
        return np.nan
    return np.nanmedian(arr[finite])


def finite_percentiles(arr, q=(16, 50, 84)):
    """Return finite-value percentiles; NaN if no finite values."""
    arr = np.asarray(arr, dtype=np.float64)
    finite = np.isfinite(arr)
    if not np.any(finite):
        return np.full(len(q), np.nan, dtype=np.float64)
    return np.nanpercentile(arr[finite], q)


def finite_minmax(arr):
    """Return (min, max) over finite values; (NaN, NaN) if none."""
    arr = np.asarray(arr, dtype=np.float64)
    finite = np.isfinite(arr)
    if not np.any(finite):
        return (np.nan, np.nan)
    return (float(np.nanmin(arr[finite])), float(np.nanmax(arr[finite])))


def extract_filter_values(filter_output, filter_name):
    """Extract numeric array for a single filter from speclite output."""
    arr = np.asarray(filter_output)
    if arr.dtype.names is not None:
        if filter_name in arr.dtype.names:
            return np.asarray(arr[filter_name], dtype=np.float64)
        if len(arr.dtype.names) == 1:
            return np.asarray(arr[arr.dtype.names[0]], dtype=np.float64)
        raise ValueError(
            f"Could not identify filter '{filter_name}' in output fields {arr.dtype.names}"
        )
    return np.asarray(arr, dtype=np.float64)


def build_ab_maggies_linear_coeff(
    filter_response,
    wave,
    flux_unit,
    basis_chunk=128,
    verbose=False,
):
    """Build exact discrete coefficients for speclite AB maggies on this wavelength grid."""
    n_wave = len(wave)
    coeff = np.empty(n_wave, dtype=np.float64)
    n_chunks = (n_wave + basis_chunk - 1) // basis_chunk
    t0 = time.perf_counter()

    for chunk_i, start in enumerate(range(0, n_wave, basis_chunk), start=1):
        end = min(start + basis_chunk, n_wave)
        block = np.zeros((end - start, n_wave), dtype=np.float64)
        block[np.arange(end - start), np.arange(start, end)] = 1.0
        block_maggies = filter_response.get_ab_maggies(block * flux_unit, wavelength=wave)
        coeff[start:end] = extract_filter_values(block_maggies, filter_response.name)

        if verbose and (
            chunk_i == 1
            or chunk_i == n_chunks
            or chunk_i % max(1, n_chunks // 8) == 0
        ):
            elapsed_s = time.perf_counter() - t0
            print(
                f"  Linear-operator chunk {chunk_i}/{n_chunks} "
                f"(wave bins {start}:{end}, elapsed {elapsed_s:.1f}s)"
            )

    return coeff


_MILES_TEMPLATE_RE = re.compile(
    r"Mch(?P<slope>\d+\.\d+)"
    r"Z(?P<mh_sign>[mp])(?P<mh>\d+\.\d+)"
    r"T(?P<age>\d+\.\d+)"
    r"_iT(?P<alpha_sign>[mp])(?P<alpha>\d+\.\d+)"
)


def parse_template_filename_key(template_path: Path) -> tuple[float, float, float] | None:
    """Parse a MILES-style SSP filename into rounded (age[Gyr], [M/H], [alpha/Fe])."""
    match = _MILES_TEMPLATE_RE.search(template_path.name)
    if match is None:
        return None

    mh = float(match.group("mh"))
    if match.group("mh_sign") == "m":
        mh *= -1.0

    alpha = float(match.group("alpha"))
    if match.group("alpha_sign") == "m":
        alpha *= -1.0

    age = float(match.group("age"))
    return (round(age, 2), round(mh, 2), round(alpha, 2))


def resolve_template_library_dir(
    library_header_value: str | None,
    script_dir: Path,
    override_dir: Path | None = None,
    rootdir: Path | None = None,
    fallback_root: Path | None = None,
) -> Path:
    """Resolve the SSP template-library directory used to create the SFH weights."""
    library_raw = str(library_header_value).strip() if library_header_value is not None else ""
    library_path = Path(library_raw.rstrip("/")) if library_raw else None
    library_name = library_path.name if library_path is not None and library_path.name else "MILES_safe"

    return resolve_existing_path(
        f"template library directory for {library_name}",
        override_dir,
        library_path,
        Path(library_name),
        script_dir / library_name,
        script_dir.parent / library_name,
        script_dir.parent / "data" / "IC3392" / library_name,
        script_dir.parent.parent / "data" / "IC3392" / library_name,
        (rootdir / library_name) if rootdir is not None else None,
        (fallback_root / library_name) if fallback_root is not None else None,
    )


def build_template_norm_flux_lookup(
    template_dir: Path,
    lmin: float,
    lmax: float,
) -> dict[tuple[float, float, float], float]:
    """
    Return the pre-normalization mean flux in the fitted spectral window for each SSP.

    This mirrors pPXF's LIGHT-mode normalization, where the output weights
    represent mean-flux fractions within norm_range / [LMIN, LMAX].
    """
    lookup: dict[tuple[float, float, float], float] = {}

    for template_path in sorted(template_dir.glob("*.fits")):
        key = parse_template_filename_key(template_path)
        if key is None:
            continue

        with fits.open(template_path, memmap=True) as hdul:
            template_flux = np.asarray(hdul[0].data, dtype=np.float64)
            template_hdr = hdul[0].header

        wave = template_hdr["CRVAL1"] + template_hdr["CDELT1"] * np.arange(
            template_hdr["NAXIS1"], dtype=np.float64
        )
        in_window = (wave >= lmin) & (wave <= lmax)
        if not np.any(in_window):
            continue

        mean_flux = float(np.nanmean(template_flux[in_window]))
        if np.isfinite(mean_flux) and mean_flux > 0:
            lookup[key] = mean_flux

    if len(lookup) == 0:
        raise RuntimeError(
            f"No valid SSP template mean-flux measurements were found in {template_dir} "
            f"for the normalization window [{lmin}, {lmax}] Angstrom."
        )

    return lookup

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

fname = phot_path

# ---------- 2.2  Load data, skip the two header lines ----------
tbl = np.genfromtxt(
    fname, dtype=None, encoding=None, names=names,
    comments='#', skip_header=2, autostrip=True)

# ---------- 2.3  Keep only Chabrier rows ----------
mask = (tbl['IMF'] == 'Ch')
phot = tbl[mask]

print(f"Loaded {len(phot)} out of {len(tbl)} rows of data from {fname.name}")


# --- 3.1  Build a lookup dict keyed by (Age[Gyr], [M/H]) rounded to 2 dec ----
key_phot = {}

for row in phot:
    age_gyr = round(float(row['Age']), 2)
    mh_dex = round(float(row['MH']), 2)
    key_phot[(age_gyr, mh_dex)] = {
        "ml_r": float(row['ML_R']),
        "v_minus_r": float(row['VminusR']),
        # Any fixed zeropoint cancels, so an R-band luminosity proxy is enough.
        "lum_r_proxy": float(10.0 ** (-0.4 * float(row['R']))),
    }

grid = grid_data  # the FITS_rec you already loaded

# 1) convert the opaque FITS_rec into plain ndarrays
w = weights_data['WEIGHTS'].astype(np.float64)             # (n_bin, n_ssp)
ml_r_ssp = np.full(len(grid), np.nan, dtype=np.float64)
v_minus_r_ssp = np.full(len(grid), np.nan, dtype=np.float64)
lum_r_proxy_ssp = np.full(len(grid), np.nan, dtype=np.float64)
grid_template_keys: list[tuple[float, float, float]] = []

for i, (logage, mh, alpha) in enumerate(grid):
    age_gyr = round(10 ** float(logage), 2)
    mh_dex = round(float(mh), 2)
    alpha_dex = round(float(alpha), 2)
    grid_template_keys.append((age_gyr, mh_dex, alpha_dex))

    phot_entry = key_phot.get((age_gyr, mh_dex))
    if phot_entry is None:
        continue

    ml_r_ssp[i] = phot_entry["ml_r"]
    v_minus_r_ssp[i] = phot_entry["v_minus_r"]
    lum_r_proxy_ssp[i] = phot_entry["lum_r_proxy"]

if not np.all(np.isfinite(ml_r_ssp)) or np.any(ml_r_ssp <= 0):
    raise ValueError("Encountered invalid SSP ML_R values from the photometry table.")
if not np.all(np.isfinite(v_minus_r_ssp)):
    raise ValueError("Encountered invalid SSP V-R colours from the photometry table.")
if not np.all(np.isfinite(lum_r_proxy_ssp)) or np.any(lum_r_proxy_ssp <= 0):
    raise ValueError("Encountered invalid SSP R-band luminosity proxies from the photometry table.")

# 2) normalize template weights per Voronoi bin before combining SSP properties.
weight_sum = w.sum(axis=1, dtype=np.float64)
if not np.all(np.isfinite(weight_sum)) or np.any(weight_sum <= 0):
    raise ValueError("Encountered non-finite or non-positive SFH weight sums.")

w_norm = w / weight_sum[:, None]
norm_check = w_norm.sum(axis=1, dtype=np.float64)

# Candidate interpretations of the saved SFH weights.
# - light_r: old behaviour, valid only if the weights are already R-band light fractions.
# - light_norm_to_r: exact LIGHT-mode conversion using the template mean flux
#   in the fitted wavelength window (ppxf norm_range / nGIST LMIN-LMAX).
# - light_v_to_r: approximate LIGHT-mode conversion using SSP V-R colours.
# - mass: use the harmonic-mean M/L_R that would apply if the weights were mass fractions.
light_r_simple = np.sum(w_norm * ml_r_ssp[None, :], axis=1, dtype=np.float64)

norm_flux_ssp = np.full(len(grid), np.nan, dtype=np.float64)
light_norm_to_r = np.full(weight_sum.shape, np.nan, dtype=np.float64)
template_library_dir: Path | None = None
template_library_status = "not requested"
template_library_error: Exception | None = None

try:
    template_library_dir = resolve_template_library_dir(
        weights_hdr.get("LIBRARY"),
        script_dir=script_dir,
        override_dir=template_library_dir_override,
        rootdir=rootdir,
        fallback_root=fallback_root,
    )
    template_norm_flux_lookup = build_template_norm_flux_lookup(
        template_library_dir,
        lmin=float(weights_hdr["LMIN"]),
        lmax=float(weights_hdr["LMAX"]),
    )

    missing_template_keys: list[tuple[float, float, float]] = []
    for i, template_key in enumerate(grid_template_keys):
        flux_value = template_norm_flux_lookup.get(template_key)
        if flux_value is None:
            missing_template_keys.append(template_key)
            continue
        norm_flux_ssp[i] = flux_value

    if missing_template_keys:
        sample_missing = ", ".join(str(key) for key in missing_template_keys[:5])
        raise KeyError(
            "Could not match all SFH grid SSPs to template spectra in "
            f"{template_library_dir}. Missing {len(missing_template_keys)} keys; "
            f"first few: {sample_missing}"
        )
    if not np.all(np.isfinite(norm_flux_ssp)) or np.any(norm_flux_ssp <= 0):
        raise ValueError("Encountered invalid SSP mean-flux values from the template library.")

    lr_over_lnorm = lum_r_proxy_ssp / norm_flux_ssp
    light_norm_to_r_num = np.sum(
        w_norm * (ml_r_ssp * lr_over_lnorm)[None, :],
        axis=1,
        dtype=np.float64,
    )
    light_norm_to_r_den = np.sum(
        w_norm * lr_over_lnorm[None, :],
        axis=1,
        dtype=np.float64,
    )
    light_norm_to_r = np.divide(
        light_norm_to_r_num,
        light_norm_to_r_den,
        out=np.full_like(light_norm_to_r_num, np.nan, dtype=np.float64),
        where=light_norm_to_r_den > 0,
    )
    template_library_status = "available"
except (FileNotFoundError, KeyError, RuntimeError, ValueError) as exc:
    template_library_error = exc
    template_library_status = f"unavailable ({exc})"

lr_over_lv = 10.0 ** (0.4 * v_minus_r_ssp)
light_v_to_r_num = np.sum(w_norm * (ml_r_ssp * lr_over_lv)[None, :], axis=1, dtype=np.float64)
light_v_to_r_den = np.sum(w_norm * lr_over_lv[None, :], axis=1, dtype=np.float64)
light_v_to_r = np.divide(
    light_v_to_r_num,
    light_v_to_r_den,
    out=np.full_like(light_v_to_r_num, np.nan, dtype=np.float64),
    where=light_v_to_r_den > 0,
)
mass_harmonic = np.divide(
    1.0,
    np.sum(w_norm / ml_r_ssp[None, :], axis=1, dtype=np.float64),
    out=np.full(weight_sum.shape, np.nan, dtype=np.float64),
    where=np.all(np.isfinite(w_norm / ml_r_ssp[None, :]), axis=1),
)

ml_bin_candidates = {
    "light_r": light_r_simple,
    "light_norm_to_r": light_norm_to_r,
    "light_v_to_r": light_v_to_r,
    "mass": mass_harmonic,
}
effective_sfh_weight_interpretation = sfh_weight_interpretation
if sfh_weight_interpretation not in ml_bin_candidates:
    raise ValueError(
        "Unsupported sfh_weight_interpretation="
        f"{sfh_weight_interpretation!r}. Expected one of {sorted(ml_bin_candidates)}."
    )
if sfh_weight_interpretation == "light_norm_to_r" and not np.all(np.isfinite(light_norm_to_r)):
    if sfh_weight_interpretation_fallback not in ml_bin_candidates:
        raise ValueError(
            "Unsupported sfh_weight_interpretation_fallback="
            f"{sfh_weight_interpretation_fallback!r}. Expected one of {sorted(ml_bin_candidates)}."
        )
    effective_sfh_weight_interpretation = sfh_weight_interpretation_fallback
    print(
        "Warning: Exact light_norm_to_r conversion is unavailable because the "
        "original SSP template library could not be resolved. "
        f"Falling back to {effective_sfh_weight_interpretation!r}."
    )
    if template_library_error is not None:
        print(f"Template-library resolution detail: {template_library_error}")
ml_bin = ml_bin_candidates[effective_sfh_weight_interpretation]

print(
    "Raw SFH weight-sum [min, p16, p50, p84, max] = "
    f"[{np.min(weight_sum):.4e}, "
    f"{np.percentile(weight_sum, 16):.4e}, "
    f"{np.percentile(weight_sum, 50):.4e}, "
    f"{np.percentile(weight_sum, 84):.4e}, "
    f"{np.max(weight_sum):.4e}]"
)
print(
    "Normalized SFH weight-sum check [min, max] = "
    f"[{np.min(norm_check):.6f}, {np.max(norm_check):.6f}]"
)
print(
    "Composite ML_R candidates [median] "
    f"(light_r, light_norm_to_r, light_v_to_r, mass_harmonic) = "
    f"[{np.nanmedian(light_r_simple):.4f}, "
    f"{np.nanmedian(light_norm_to_r):.4f}, "
    f"{np.nanmedian(light_v_to_r):.4f}, "
    f"{np.nanmedian(mass_harmonic):.4f}]"
)
if template_library_dir is not None:
    print(
        "Template-library LIGHT normalization for ML_R: "
        f"{template_library_dir} (LMIN={float(weights_hdr['LMIN']):.1f}, "
        f"LMAX={float(weights_hdr['LMAX']):.1f} Angstrom)"
    )
else:
    print(
        "Template-library LIGHT normalization for ML_R: "
        f"{template_library_status}"
    )
print(
    "Selected SFH-weight interpretation for ML_R: "
    f"{effective_sfh_weight_interpretation}"
)

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
binid_max = int(np.nanmax(binning_BINID[valid])) if np.any(valid) else -1
if binid_max >= len(ml_bin):
    raise ValueError(
        "BINID values exceed the number of rows in the SFH weights file: "
        f"max BINID={binid_max}, n_weight_rows={len(ml_bin)}"
    )
print(f"BINID range check for ML_R mapping: max BINID={binid_max}, n_weight_rows={len(ml_bin)}")
ml_lut = np.full(binid_max + 1, np.nan, dtype=np.float32)
ml_lut[:] = ml_bin[: binid_max + 1].astype(np.float32)
binning_MLR[valid] = ml_lut[binning_BINID[valid].astype(int)]


# ---------------------------------------------------------------------
# 1.  R-band magnitude map per spaxel using speclite
# ---------------------------------------------------------------------
print("Loading:", cube_path.resolve())
has_stat = False
with fits.open(cube_path, memmap=True) as cube:
    data = cube["DATA"].data                                 # (nz, ny, nx) memmap
    hdr  = cube["DATA"].header
    nz, ny, nx = data.shape
    stat_data = None

    if "STAT" in cube and not disable_stat_propagation:
        stat_data = cube["STAT"].data
        if stat_data.shape != data.shape:
            print("Warning: STAT extension shape does not match DATA; ignoring STAT.")
            stat_data = None
        else:
            has_stat = True
            print("Using STAT variance extension for photometric error propagation.")
    elif "STAT" in cube and disable_stat_propagation:
        print("STAT extension found but STAT-based propagation is disabled; using bin-wise scatter fallback.")
    else:
        print("Warning: No STAT extension found; using bin-wise scatter fallback for photometric uncertainty.")

    data_gib = data.nbytes / (1024.0 ** 3)
    print(f"Cube dimensions (nz, ny, nx): {nz}, {ny}, {nx}")
    print(f"Approx. DATA footprint: {data_gib:.2f} GiB")
    if stat_data is not None:
        stat_gib = stat_data.nbytes / (1024.0 ** 3)
        print(f"Approx. STAT footprint: {stat_gib:.2f} GiB")

    # wavelength grid (native header unit → Å)
    spec_wcs = WCS(hdr).sub(["spectral"])                  # 1-axis WCS
    wave_native = spec_wcs.all_pix2world(
        np.arange(nz)[:, None], 0
    )[:, 0]                                                 # numeric values
    wave = (wave_native * spec_wcs.wcs.cunit[0]).to(u.AA)   # use cunit[0]

    # Load filter using speclite (more robust than synphot for this purpose)
    f_r = filters.load_filter('bessell-R')
    # decamDR1noatm-r, bessell-R or 'sdss2010noatm-r'/'sdss2010-r' for SDSS r-band

    flux_unit = 1e-20 * u.erg / (u.s * u.cm**2 * u.AA)
    F0_ref = 3631e-23  # Reference flux in erg s⁻¹ cm⁻² Hz⁻¹ for AB magnitude zero point

    print("Building exact speclite linear operator for AB maggies...")
    basis_chunk = max(64, min(512, row_block_size * 64))
    maggies_coeff = build_ab_maggies_linear_coeff(
        f_r,
        wave,
        flux_unit,
        basis_chunk=basis_chunk,
        verbose=True,
    )
    maggies_coeff_sq = maggies_coeff ** 2
    support = maggies_coeff != 0
    support_count = int(np.sum(support))
    support_idx = np.flatnonzero(support)
    if support_idx.size == 0:
        raise RuntimeError("Filter support is empty on this wavelength grid; cannot compute AB maggies.")

    support_contiguous = bool(np.all(np.diff(support_idx) == 1))
    if support_contiguous:
        support_sel = slice(int(support_idx[0]), int(support_idx[-1]) + 1)
        coeff_eff = maggies_coeff[support_sel]
        coeff_sq_eff = maggies_coeff_sq[support_sel]
        print(
            "Using contiguous filter-support slice: "
            f"[{support_sel.start}:{support_sel.stop}] ({coeff_eff.size} bins)."
        )
    else:
        support_sel = support_idx
        coeff_eff = maggies_coeff[support_sel]
        coeff_sq_eff = maggies_coeff_sq[support_sel]
        print(f"Using non-contiguous filter-support index set ({coeff_eff.size} bins).")

    n_wave_eff = coeff_eff.size
    coeff_support_total = float(np.sum(np.abs(coeff_eff)))
    if coeff_support_total <= 0:
        raise RuntimeError("Effective filter-support coefficient sum is not positive.")

    print(
        "Photometric propagation method: exact discrete first-order Jacobian through "
        f"speclite AB-maggies operator (support bins={support_count}/{nz})."
    )
    print(
        "STAT assumption: independent spectral-bin variances (no covariance provided by cube)."
    )

    maggies_map = np.full((ny, nx), np.nan, dtype=np.float64)
    filter_support_fraction_map = np.full((ny, nx), np.nan, dtype=np.float32)
    if has_stat:
        maggies_var_map = np.full((ny, nx), np.nan, dtype=np.float64)

    strict_allfinite_pixels = 0
    accepted_support_pixels = 0

    est_io_factor = 2 if has_stat else 1
    est_io_gib = (
        float(n_wave_eff) * float(ny) * float(nx) * float(data.dtype.itemsize) * float(est_io_factor)
    ) / (1024.0 ** 3)
    stream_label = "DATA+STAT" if has_stat else "DATA only"
    print(
        "Computing AB maggies with vectorized blocks... "
        f"(estimated support-slice stream volume: {est_io_gib:.2f} GiB, {stream_label})"
    )
    loop_t0 = time.perf_counter()
    n_blocks = (ny + row_block_size - 1) // row_block_size
    progress_every_blocks = max(1, 16 // row_block_size)
    last_print_t = loop_t0

    for block_i, j0 in enumerate(range(0, ny, row_block_size), start=1):
        j1 = min(j0 + row_block_size, ny)
        block_t0 = time.perf_counter()

        raw_block = data[support_sel, j0:j1, :]                            # (nz_eff, b, nx)
        finite_flux = np.isfinite(raw_block)
        finite_flux_support = np.all(finite_flux, axis=0)
        finite_flux_weight = np.einsum(
            "z,zpn->pn",
            np.abs(coeff_eff),
            finite_flux.astype(np.float32),
            optimize=True,
        )
        support_fraction = finite_flux_weight / coeff_support_total
        good_support = support_fraction >= min_filter_support_fraction
        strict_allfinite_pixels += int(np.count_nonzero(finite_flux_support))
        accepted_support_pixels += int(np.count_nonzero(good_support))

        zero = np.array(0.0, dtype=raw_block.dtype)
        block_flux = np.where(finite_flux, raw_block, zero).reshape(n_wave_eff, -1)
        block_maggies = np.einsum("z,zp->p", coeff_eff, block_flux, optimize=True)
        block_maggies = block_maggies.reshape(j1 - j0, nx)
        with np.errstate(divide="ignore", invalid="ignore"):
            block_maggies = np.where(good_support, block_maggies / support_fraction, np.nan)
        maggies_map[j0:j1, :] = block_maggies
        filter_support_fraction_map[j0:j1, :] = support_fraction.astype(np.float32)

        if has_stat:
            assert stat_data is not None
            raw_stat = stat_data[support_sel, j0:j1, :]
            finite_stat = np.isfinite(raw_stat) & (raw_stat >= 0)
            finite_stat_weight = np.einsum(
                "z,zpn->pn",
                np.abs(coeff_eff),
                finite_stat.astype(np.float32),
                optimize=True,
            )
            stat_support_fraction = finite_stat_weight / coeff_support_total
            good_stat_support = stat_support_fraction >= min_filter_support_fraction
            block_stat = np.where(finite_stat, raw_stat, zero).reshape(n_wave_eff, -1)
            block_var = np.einsum("z,zp->p", coeff_sq_eff, block_stat, optimize=True)
            block_var = block_var.reshape(j1 - j0, nx)
            with np.errstate(divide="ignore", invalid="ignore"):
                block_var = np.where(
                    good_support & good_stat_support,
                    block_var / support_fraction**2,
                    np.nan,
                )
            maggies_var_map[j0:j1, :] = block_var

        block_t1 = time.perf_counter()
        block_sec = block_t1 - block_t0
        elapsed_sec = block_t1 - loop_t0
        eta_sec = (elapsed_sec / block_i) * (n_blocks - block_i) if block_i > 0 else np.nan
        frac = 100.0 * j1 / ny
        if (
            block_i == 1
            or block_i % progress_every_blocks == 0
            or j1 == ny
            or (block_t1 - last_print_t) >= 20.0
        ):
            print(
                f"Processed rows {j0}:{j1}/{ny} "
                f"({frac:5.1f}%, block {block_i}/{n_blocks}, "
                f"block {block_sec:.2f}s, elapsed {elapsed_sec/60.0:.2f} min, ETA {eta_sec/60.0:.2f} min)"
            )
            last_print_t = block_t1

    total_min = (time.perf_counter() - loop_t0) / 60.0
    throughput_mib_s = (est_io_gib * 1024.0) / (total_min * 60.0) if total_min > 0 else np.nan
    print(f"AB maggies vectorized pass completed in {total_min:.2f} min")
    if np.isfinite(throughput_mib_s):
        print(
            "Estimated effective streaming throughput: "
            f"{throughput_mib_s:.1f} MiB/s over support-slice read volume"
        )

    finite_support_fraction = filter_support_fraction_map[np.isfinite(filter_support_fraction_map)]
    if finite_support_fraction.size > 0:
        p16_sup, p50_sup, p84_sup = np.percentile(finite_support_fraction, [16, 50, 84])
        print(
            "Filter-support fraction [p16, p50, p84] = "
            f"[{p16_sup:.6f}, {p50_sup:.6f}, {p84_sup:.6f}]"
        )
    print(
        "Strict all-finite support pixels vs accepted relaxed-support pixels = "
        f"{strict_allfinite_pixels:,} vs {accepted_support_pixels:,} "
        f"(threshold={min_filter_support_fraction:.4f})"
    )

    m_r_map = np.where(maggies_map > 0, -2.5 * np.log10(maggies_map), np.nan).astype(np.float32)

    flux_map = F0_ref * maggies_map
    flux_err_map = np.full_like(flux_map, np.nan, dtype=np.float64)
    if has_stat:
        with np.errstate(invalid="ignore"):
            flux_err_map = F0_ref * np.sqrt(np.where(maggies_var_map >= 0, maggies_var_map, np.nan))

    # ---------------------------------------------------------------------
    # 2.  Collapse to Voronoi bins
    # ---------------------------------------------------------------------
    with fits.open(vor_path) as hd_v:
        # ---------- PATCH BEGIN ----------
        # keep BINID as float so NaNs survive the read
        BINID_f32 = hd_v["BINID"].data.astype(np.float32)
        muse_hdr2 = hd_v["BINID"].header

        # build an *integer* copy with sentinel -1 for “no bin”
        bad_pix = (~np.isfinite(BINID_f32)) | (BINID_f32 < 0)
        BINID = np.full(BINID_f32.shape, -1, dtype=np.int32)
        good_bin = ~bad_pix
        BINID[good_bin] = BINID_f32[good_bin].astype(np.int32)
        # ---------- PATCH END ----------

    uniq = np.unique(BINID[BINID >= 0])                 # keep this line

    # Explicitly release cube-backed arrays before downstream processing.
    del data
    if has_stat:
        del stat_data

gc.collect()

# Keep the cube-photometry validity separate from the binning-product mask.
# Otherwise the "total" R-band flux changes when the binning map changes.
phot_valid_mask = np.isfinite(m_r_map)
bin_member_mask = BINID >= 0
pixel_valid_for_binned_maps = bin_member_mask & phot_valid_mask

# Average flux in each bin using only pixels with valid cube photometry.
flux_map_clean = np.where(pixel_valid_for_binned_maps, flux_map, np.nan)
valid_pixels = sum_labels((~np.isnan(flux_map_clean)).astype(np.int32), BINID, uniq)
sum_flux = sum_labels(np.nan_to_num(flux_map_clean), BINID, uniq)
mean_flux_err_stat = np.full_like(sum_flux, np.nan, dtype=np.float64)
mean_flux = np.divide(sum_flux, valid_pixels, 
                     out=np.full_like(sum_flux, np.nan), 
                     where=valid_pixels > 0)

# Always compute a sample-based fallback uncertainty on the bin mean flux.
flux_sq = np.where(np.isnan(flux_map_clean), np.nan, flux_map_clean**2)
sum_flux_sq = sum_labels(np.nan_to_num(flux_sq), BINID, uniq)
sample_var = np.full_like(sum_flux, np.nan, dtype=np.float64)
multi_pix = valid_pixels > 1
sample_var[multi_pix] = (
    sum_flux_sq[multi_pix]
    - (sum_flux[multi_pix] ** 2) / valid_pixels[multi_pix]
) / (valid_pixels[multi_pix] - 1)
# For singleton bins, scatter is undefined; use zero-variance fallback.
sample_var[valid_pixels == 1] = 0.0
sample_var = np.where(sample_var >= 0, sample_var, 0.0)
mean_flux_err_sample = np.divide(
    np.sqrt(sample_var),
    np.sqrt(valid_pixels),
    out=np.full_like(sum_flux, np.nan, dtype=np.float64),
    where=valid_pixels > 0,
)

if has_stat:
    err_valid = pixel_valid_for_binned_maps & np.isfinite(flux_map) & np.isfinite(flux_err_map)
    valid_err_pixels = sum_labels(err_valid.astype(np.int32), BINID, uniq)
    sum_var_flux = sum_labels(np.where(err_valid, flux_err_map**2, 0.0), BINID, uniq)
    mean_flux_err_stat = np.divide(
        np.sqrt(sum_var_flux),
        valid_err_pixels,
        out=np.full_like(sum_flux, np.nan, dtype=np.float64),
        where=valid_err_pixels > 0,
    )

# Prefer STAT-propagated errors where valid; otherwise fall back to sample-based.
mean_flux_err = np.where(
    np.isfinite(mean_flux_err_stat) & (mean_flux_err_stat >= 0),
    mean_flux_err_stat,
    mean_flux_err_sample,
)

print(
    f"Finite binned flux-error bins: {np.sum(np.isfinite(mean_flux_err))}/{mean_flux_err.size} "
    f"(STAT={np.sum(np.isfinite(mean_flux_err_stat))}, sample={np.sum(np.isfinite(mean_flux_err_sample))})"
)

# Convert back to magnitudes
mean_mag = -2.5 * np.log10(mean_flux / F0_ref)
with np.errstate(divide="ignore", invalid="ignore"):
    mean_mag_err = np.where(
        (mean_flux > 0) & np.isfinite(mean_flux_err),
        (2.5 / np.log(10.0)) * (mean_flux_err / mean_flux),
        np.nan,
    )

with np.errstate(divide="ignore", invalid="ignore"):
    mean_rel_flux_err = np.where(
        (mean_flux > 0) & np.isfinite(mean_flux_err),
        mean_flux_err / mean_flux,
        np.nan,
    )

print(f"Finite binned magnitude-error bins: {np.sum(np.isfinite(mean_mag_err))}/{mean_mag_err.size}")
median_rel_flux_err = finite_median(mean_rel_flux_err)
if np.isfinite(median_rel_flux_err):
    print(f"Median binned relative flux error : {median_rel_flux_err:.4e}")
else:
    print("Median binned relative flux error : N/A (no finite values)")

min_mag_err, max_mag_err = finite_minmax(mean_mag_err)
if np.isfinite(min_mag_err):
    print(f"Binned magnitude error min/max   : {min_mag_err:.4e} .. {max_mag_err:.4e} mag")

p16_mag, p50_mag, p84_mag = finite_percentiles(mean_mag_err)
if np.isfinite(p50_mag):
    print(
        "Binned magnitude error [p16, p50, p84] "
        f"= [{p16_mag:.4e}, {p50_mag:.4e}, {p84_mag:.4e}] mag"
    )
else:
    print("Binned magnitude error [p16, p50, p84] = N/A (no finite values)")

# Create lookup table for bin-averaged magnitudes
lut = np.full(int(BINID.max()) + 1, np.nan, dtype=np.float32)
lut[uniq] = mean_mag
m_r_binned = np.where(pixel_valid_for_binned_maps, lut[BINID], np.nan)     # (ny, nx)

lut_err = np.full(int(BINID.max()) + 1, np.nan, dtype=np.float32)
lut_err[uniq] = mean_mag_err.astype(np.float32)
m_r_binned_err = np.where(pixel_valid_for_binned_maps, lut_err[BINID], np.nan)

flux_map_binned = F0_ref * 10**(-0.4 * m_r_binned)  # Convert mag to flux

# ---------------------------------------------------------------------
# 3.  Galactic-extinction correction (Calzetti 2000)
# ---------------------------------------------------------------------
# 1) Read the stellar EBV map from nGIST (internal dust for stars)
with fits.open(sfh_path) as sfh_hdul:
    EBV = sfh_hdul["EBV"].data.astype(np.float32)
    EBV_ERR = None
    EBV_ERR_extname = None
    for extname in ("EBV_ERR", "EBV_ERROR", "E_BV_ERR", "EBV_STD", "EBV_SIGMA"):
        if extname in sfh_hdul:
            EBV_ERR = sfh_hdul[extname].data.astype(np.float32)
            EBV_ERR_extname = extname
            break

if EBV_ERR is not None:
    print(f"Using {EBV_ERR_extname} extension for EBV error propagation.")
else:
    EBV_ERR = np.full_like(EBV, np.nan, dtype=np.float32)
    print("Warning: No EBV error extension found; EBV uncertainty term will be omitted.")

# 2) Effective wavelength of your R-band filter, in microns
#    (rough numbers: Bessell R ~ 0.64 µm, SDSS r ~ 0.623 µm)
if 'bessell' in f_r.name.lower():
    lam_eff_um = 0.64
    M_r_sun = 4.61    # Solar absolute magnitude in Bessell R (AB magnitude)
elif 'sdss' in f_r.name.lower():
    lam_eff_um = 0.623
    M_r_sun = 4.65    # Solar absolute magnitude in SDSS r (AB magnitude)
else:
    print(f"Warning: Unknown filter {f_r.name}, using SDSS coefficients")
    lam_eff_um = 0.623
    M_r_sun = 4.65

# 3) Calzetti k(λ) for that band
k_r_calz = calzetti_k(lam_eff_um)

# 4) Internal attenuation of the stellar continuum
A_r = k_r_calz * EBV
A_r_err = np.where(np.isfinite(EBV_ERR), np.abs(k_r_calz) * EBV_ERR, np.nan)

# 5) Correct magnitudes for *internal* attenuation
m_r_corr = m_r_binned - A_r
m_r_corr_err = quadrature_sum(m_r_binned_err, A_r_err)
m_r_corr[~pixel_valid_for_binned_maps] = np.nan
m_r_corr_err[~pixel_valid_for_binned_maps] = np.nan
m_r_corr_err[~np.isfinite(m_r_corr)] = np.nan

# magnitude back to nanomaggies in Legacy survey format
def magnitude_to_nanomaggies(magnitude):
    return 10**((22.5 - magnitude) / 2.5)
FLUX_R_corr = magnitude_to_nanomaggies(m_r_corr)  # Convert to flux
FLUX_R_corr_err = np.abs(FLUX_R_corr) * (np.log(10.0) / 2.5) * m_r_corr_err
FLUX_R_corr_err[~np.isfinite(FLUX_R_corr)] = np.nan

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
logM_star_err = np.where(M_star > 0, 0.4 * m_r_corr_err, np.nan)
M_star_err = np.abs(M_star) * np.log(10.0) * logM_star_err
logM_star[~pixel_valid_for_binned_maps] = np.nan
logM_star_err[~pixel_valid_for_binned_maps] = np.nan

print("R-band magnitude calculation completed!")
print(f"Filter used: {f_r.name}")
print(f"Magnitude range: {np.nanmin(m_r_corr):.2f} to {np.nanmax(m_r_corr):.2f}")
print(f"Distance modulus: {distmod:.2f} mag")

# ──────────────────────────────────────────────────────────────────
#  Integrated summary lines
# ──────────────────────────────────────────────────────────────────
total_flux_uncorrected = np.nansum(np.where(phot_valid_mask, flux_map, np.nan))
total_flux_corrected = np.nansum(np.where(np.isfinite(FLUX_R_corr), FLUX_R_corr, np.nan))
tot_mag_uncorrected = -2.5 * np.log10(total_flux_uncorrected / F0_ref)
tot_mag_corrected = 22.5 - 2.5 * np.log10(total_flux_corrected)
total_L_R_linear = np.nansum(L_Lsun)
total_M_R_linear = np.nansum(M_star)
tot_L_R  = np.log10(total_L_R_linear)          # log₁₀(L/L☉)
tot_M_R  = np.log10(total_M_R_linear)          # log₁₀(M/M☉)
integrated_ml_r = total_M_R_linear / total_L_R_linear

print(f"Total R-band magnitude (uncorrected) : {tot_mag_uncorrected:.3f} mag (AB)")
print(f"Total R-band magnitude (corrected)   : {tot_mag_corrected:.3f} mag (AB)")
print(f"Total R-band luminosity              : {tot_L_R:.3f} log10(L☉)")
print(f"Total stellar mass (R)               : {tot_M_R:.3f} log10(M☉)")
print(f"Integrated stellar M/L_R             : {integrated_ml_r:.4f} M☉/L☉")
median_logmstar_err = finite_median(logM_star_err)
if np.isfinite(median_logmstar_err):
    p16_logm, p50_logm, p84_logm = finite_percentiles(logM_star_err)
    print(f"Median LOGMSTAR_ERR      : {median_logmstar_err:.4e} dex")
    print(
        "LOGMSTAR_ERR [p16, p50, p84] "
        f"= [{p16_logm:.4e}, {p50_logm:.4e}, {p84_logm:.4e}] dex"
    )
    min_logm_err, max_logm_err = finite_minmax(logM_star_err)
    if np.isfinite(min_logm_err):
        print(f"LOGMSTAR_ERR min/max     : {min_logm_err:.4e} .. {max_logm_err:.4e} dex")
else:
    print("Median LOGMSTAR_ERR      : N/A (no finite values)")

print(
    "Error budget note (LOGMSTAR): includes photometric+extinction propagation; "
    "excludes M/L model, distance, and inclination/area systematic uncertainties."
)


# 1) read the whole HDUList from disk
with fits.open(binning_path) as hdul:
    # 2) clone all existing HDUs into a new list
    new_hdul = fits.HDUList([hdu.copy() for hdu in hdul])

# 3) build the new FLUX_R_corr image HDU
FLUX_R_corr_hdu = fits.ImageHDU(
    data=FLUX_R_corr.astype(np.float64), name="FLUX_R_corr")

FLUX_R_corr_err_hdu = fits.ImageHDU(
    data=FLUX_R_corr_err.astype(np.float64), name="FLUX_R_corr_ERR")

# 4) keep WCS and pixel-scale info by copying the original BINID header
FLUX_R_corr_hdu.header.update(binning_hdr)         # you created binning_hdr earlier
FLUX_R_corr_hdu.header["EXTNAME"] = "FLUX_R_corr"
FLUX_R_corr_hdu.header["BUNIT"]   = "nanomaggies"  # physical units

FLUX_R_corr_err_hdu.header.update(binning_hdr)
FLUX_R_corr_err_hdu.header["EXTNAME"] = "FLUX_R_corr_ERR"
FLUX_R_corr_err_hdu.header["BUNIT"]   = "nanomaggies"
FLUX_R_corr_err_hdu.header["COMMENT"] = "1-sigma uncertainty in FLUX_R_corr"

# 5) append and write to disk
new_hdul.append(FLUX_R_corr_hdu)
new_hdul.append(FLUX_R_corr_err_hdu)
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

    mass_err_hdu = fits.ImageHDU(data=logM_star_err.astype(np.float64),
                                 header=binning_hdr, name="LOGMSTAR_ERR")
    mass_err_hdu.header["BUNIT"] = "dex"
    mass_err_hdu.header["COMMENT"] = "1-sigma uncertainty in LOGMSTAR"
    hdul.append(mass_err_hdu)

    hdul.flush()                                             # write in-place

print("Stellar-mass layer saved ➜", out_path.resolve())


# Getting the stellar mass surface density
# Convert to surface density in M☉/kpc²
# 1. Convert pixel area to physical area in kpc²
legacy_wcs2 = WCS(binning_hdr).celestial  # strip spectral axis
pixel_scale = (proj_plane_pixel_scales(legacy_wcs2) * u.deg).to(u.arcsec)
pixel_area_Mpc = ((pixel_scale[0]).to(u.rad).value*16.5*u.Mpc)*(((pixel_scale[1]).to(u.rad).value*16.5*u.Mpc))
pixel_area_kpc = pixel_area_Mpc.to(u.kpc**2)

# 2. Read galaxy inclination and calculate correction factor
if apply_inclination_correction:
    galaxy_inclination = read_galaxy_inclination(galaxy, str(inclination_path))
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
stellar_mass_surface_density_err = M_star_err / pixel_area_kpc
stellar_mass_surface_density_corrected_err = stellar_mass_surface_density_err * ba_factor
# stellar_mass_surface_density_corrected = stellar_mass_surface_density 

# 4. Convert to log10 scale
log_stellar_mass_surface_density = np.where(
    stellar_mass_surface_density_corrected.value > 0,
    np.log10(stellar_mass_surface_density_corrected.value),
    np.nan,
)
with np.errstate(divide="ignore", invalid="ignore"):
    log_stellar_mass_surface_density_err = np.where(
        stellar_mass_surface_density_corrected.value > 0,
        stellar_mass_surface_density_corrected_err.value
        / (stellar_mass_surface_density_corrected.value * np.log(10.0)),
        np.nan,
    )
log_stellar_mass_surface_density_err[~pixel_valid_for_binned_maps] = np.nan
median_logsigma_err = finite_median(log_stellar_mass_surface_density_err)
if np.isfinite(median_logsigma_err):
    p16_logsig, p50_logsig, p84_logsig = finite_percentiles(log_stellar_mass_surface_density_err)
    print(f"Median LOGMASS_SURFACE_DENSITY_ERR: {median_logsigma_err:.4e} dex")
    print(
        "LOGMASS_SURFACE_DENSITY_ERR [p16, p50, p84] "
        f"= [{p16_logsig:.4e}, {p50_logsig:.4e}, {p84_logsig:.4e}] dex"
    )
    min_logsig_err, max_logsig_err = finite_minmax(log_stellar_mass_surface_density_err)
    if np.isfinite(min_logsig_err):
        print(
            "LOGMASS_SURFACE_DENSITY_ERR min/max: "
            f"{min_logsig_err:.4e} .. {max_logsig_err:.4e} dex"
        )
else:
    print("Median LOGMASS_SURFACE_DENSITY_ERR: N/A (no finite values)")

mass_err_common = np.isfinite(logM_star_err) & np.isfinite(log_stellar_mass_surface_density_err)
if np.any(mass_err_common):
    delta_mass_err = np.abs(
        logM_star_err[mass_err_common] - log_stellar_mass_surface_density_err[mass_err_common]
    )
    print(
        "Consistency check |LOGMSTAR_ERR - LOGMASS_SURFACE_DENSITY_ERR| "
        f"median={np.nanmedian(delta_mass_err):.4e} dex, max={np.nanmax(delta_mass_err):.4e} dex"
    )
else:
    print(
        "Consistency check |LOGMSTAR_ERR - LOGMASS_SURFACE_DENSITY_ERR|: "
        "N/A (no common finite values)"
    )


with fits.open(out_path, mode="append") as hdul:             # open existing file
    mass_density_hdu = fits.ImageHDU(
        data=log_stellar_mass_surface_density.astype(np.float64),  # like others
        header=binning_hdr, name="LOGMASS_SURFACE_DENSITY")
    mass_density_hdu.header["BUNIT"] = "log(Msol/kpc2)"  # units keyword
    hdul.append(mass_density_hdu)                                    # add as 10-th HDU 

    mass_density_err_hdu = fits.ImageHDU(
        data=log_stellar_mass_surface_density_err.astype(np.float64),
        header=binning_hdr, name="LOGMASS_SURFACE_DENSITY_ERR")
    mass_density_err_hdu.header["BUNIT"] = "dex"
    mass_density_err_hdu.header["COMMENT"] = "1-sigma uncertainty in LOGMASS_SURFACE_DENSITY"
    hdul.append(mass_density_err_hdu)

    hdul.flush()                                             # write in-place

print("Stellar mass surface density layer saved ➜", out_path.resolve())


with fits.open(out_path, mode="append") as hdul:             # open existing file
    m_r_hdu = fits.ImageHDU(data=m_r_corr.astype(np.float64),  # like others
                             header=binning_hdr, name="magnitude_r")
    m_r_hdu.header["BUNIT"] = "mag_AB"                   # units keyword
    hdul.append(m_r_hdu)                                    # add as 9-th HDU

    m_r_err_hdu = fits.ImageHDU(data=m_r_corr_err.astype(np.float64),
                                header=binning_hdr, name="magnitude_r_ERR")
    m_r_err_hdu.header["BUNIT"] = "mag_AB"
    m_r_err_hdu.header["COMMENT"] = "1-sigma uncertainty in magnitude_r"
    hdul.append(m_r_err_hdu)

    hdul.flush()                                             # write in-place

print("r-band magnitude layer saved ➜", out_path.resolve())


with fits.open(out_path, mode="append") as hdul:             # open existing file
    m_r_uncorrected_hdu = fits.ImageHDU(data=m_r_binned.astype(np.float64),  # like others
                             header=binning_hdr, name="magnitude_r_uncorrected")
    m_r_uncorrected_hdu.header["BUNIT"] = "mag_AB"                   # units keyword
    hdul.append(m_r_uncorrected_hdu)                                    # add as 9-th HDU

    m_r_uncorrected_err_hdu = fits.ImageHDU(
        data=m_r_binned_err.astype(np.float64),
        header=binning_hdr, name="magnitude_r_uncorrected_ERR")
    m_r_uncorrected_err_hdu.header["BUNIT"] = "mag_AB"
    m_r_uncorrected_err_hdu.header["COMMENT"] = "1-sigma uncertainty in magnitude_r_uncorrected"
    hdul.append(m_r_uncorrected_err_hdu)

    hdul.flush()                                             # write in-place

print("r-band uncorrected magnitude layer saved ➜", out_path.resolve())
