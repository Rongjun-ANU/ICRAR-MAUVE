#!/usr/bin/env python3
"""Shared PHANGS STAT-HDU variance-gap checks and fixes.

The workflow targets the three PHANGS-native cubes staged in
/scratch/pawsey1308/mauve/cubes/v3tk. Detection is restricted to the useful
4700-9350 A open wavelength interval, excludes the AO/LGS gap, and ignores
spatial pixels masked by {GALID}_mask.fits. A bad STAT sample with finite DATA
is filled spectrally from one or two immediate positive STAT bounds. DATA may be
positive, zero, or negative in this spectral case and is preserved. Joint
DATA/STAT NaNs are grouped per wavelength into 8-connected spatial components;
only components with a complete one-pixel boundary containing finite DATA and
positive finite STAT are filled. Each accepted voxel uses original neighbors
inside a 1 arcsec-diameter circle (2.5-pixel radius).
"""

import argparse
import concurrent.futures
import os
import shutil
import sys
from datetime import datetime
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple, Union

import numpy as np
from astropy.io import fits
from scipy import ndimage


DEFAULT_GALAXIES = ("NGC4254", "NGC4321", "NGC4535")
DEFAULT_CUBE_DIR = Path("/scratch/pawsey1308/mauve/cubes/v3tk")
LOG_NAME = "check_phangs_variance.log"
WAVE_MIN = 4700.0
WAVE_MAX = 9350.0
AO_GAP_RANGES = ((5800.0, 5970.0),)
SPATIAL_DIAMETER_ARCSEC = 1.0
PIXEL_SCALE_ARCSEC = 0.2
SPATIAL_RADIUS_PIXELS = SPATIAL_DIAMETER_ARCSEC / (2.0 * PIXEL_SCALE_ARCSEC)
EIGHT_CONNECTED_STRUCTURE = np.ones((3, 3), dtype=np.uint8)


class StatGap(object):
    __slots__ = (
        "z_start",
        "z_end",
        "y",
        "x",
        "fillable",
        "reason",
        "wave_start",
        "wave_end",
        "lower_z",
        "upper_z",
    )

    def __init__(
        self,
        z_start,
        z_end,
        y,
        x,
        fillable,
        reason,
        wave_start=None,
        wave_end=None,
        lower_z=None,
        upper_z=None,
    ):
        self.z_start = z_start
        self.z_end = z_end
        self.y = y
        self.x = x
        self.fillable = fillable
        self.reason = reason
        self.wave_start = wave_start
        self.wave_end = wave_end
        self.lower_z = lower_z
        self.upper_z = upper_z

    @property
    def length(self):
        return self.z_end - self.z_start + 1


class SpatialGap(object):
    __slots__ = (
        "z",
        "y",
        "x",
        "fillable",
        "reason",
        "wave",
        "data_fill",
        "stat_fill",
    )

    def __init__(
        self,
        z,
        y,
        x,
        fillable,
        reason,
        wave=None,
        data_fill=None,
        stat_fill=None,
    ):
        self.z = z
        self.y = y
        self.x = x
        self.fillable = fillable
        self.reason = reason
        self.wave = wave
        self.data_fill = data_fill
        self.stat_fill = stat_fill


class ResidualStatSample(object):
    __slots__ = ("z", "y", "x", "wave", "data_value", "stat_value")

    def __init__(self, z, y, x, wave, data_value, stat_value):
        self.z = z
        self.y = y
        self.x = x
        self.wave = wave
        self.data_value = data_value
        self.stat_value = stat_value


class ResidualStatAudit(object):
    __slots__ = ("count", "samples")

    def __init__(self, count, samples):
        self.count = count
        self.samples = samples


class GapReport(object):
    __slots__ = (
        "galid",
        "cube_path",
        "shape",
        "gaps",
        "spatial_gaps",
        "residual_audit",
    )

    def __init__(
        self,
        galid,
        cube_path,
        shape,
        gaps,
        spatial_gaps=(),
        residual_audit=None,
    ):
        self.galid = galid
        self.cube_path = cube_path
        self.shape = shape
        self.gaps = gaps
        self.spatial_gaps = spatial_gaps
        self.residual_audit = residual_audit

    @property
    def spectral_fillable_count(self):
        return sum(1 for gap in self.gaps if gap.fillable)

    @property
    def spatial_fillable_count(self):
        return sum(1 for gap in self.spatial_gaps if gap.fillable)

    @property
    def fillable_count(self):
        return self.spectral_fillable_count + self.spatial_fillable_count

    @property
    def unfillable_count(self):
        return len(self.gaps) + len(self.spatial_gaps) - self.fillable_count


def cube_path_for(galid: str, cube_dir: Path = DEFAULT_CUBE_DIR) -> Path:
    return cube_dir / f"{galid}_PHANGS_DATACUBE_native.fits"


def fixed_path_for(galid: str, cube_dir: Path = DEFAULT_CUBE_DIR) -> Path:
    return cube_dir / f"{galid}_PHANGS_DATACUBE_native_fixed.fits"


def mask_path_for(galid: str, cube_dir: Path = DEFAULT_CUBE_DIR) -> Path:
    return cube_dir / f"{galid}_mask.fits"


def parse_galaxies(argv=None):
    parser = argparse.ArgumentParser(
        description="Check or fix non-positive STAT gaps in PHANGS native cubes."
    )
    parser.add_argument(
        "galaxies",
        nargs="*",
        help="Galaxy IDs to process. Defaults to NGC4254 NGC4321 NGC4535.",
    )
    parser.add_argument(
        "--cube-dir",
        default=str(DEFAULT_CUBE_DIR),
        help=f"Directory containing PHANGS native cubes. Default: {DEFAULT_CUBE_DIR}",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1)),
        help="Number of worker processes for scanning. Default: SLURM_CPUS_PER_TASK or CPU count.",
    )
    parser.add_argument(
        "--log",
        default=LOG_NAME,
        help=f"Log path for check output. Default: {LOG_NAME}",
    )
    return parser.parse_args(argv)


def selected_galaxies(names):
    if not names:
        return DEFAULT_GALAXIES
    return tuple(name.strip().upper() for name in names if name.strip())


def positive_finite(values: np.ndarray) -> np.ndarray:
    return np.isfinite(values) & (values > 0)


def wavelength_axis_from_header(header, nz):
    crval = header.get("CRVAL3")
    crpix = header.get("CRPIX3", 1.0)
    cdelt = header.get("CDELT3", header.get("CD3_3"))
    if crval is None or cdelt is None:
        return None
    wave = float(crval) + (np.arange(nz, dtype=float) + 1.0 - float(crpix)) * float(cdelt)
    unit = str(header.get("CUNIT3", "")).lower()
    if unit in ("m", "meter", "metre"):
        wave = wave * 1.0e10
    elif unit in ("nm", "nanometer", "nanometre"):
        wave = wave * 10.0
    return wave


def wavelength_axis_from_hdul(hdul, nz):
    for hdu_name in ("STAT", "DATA", "FLUX"):
        if hdu_name in hdul:
            wave = wavelength_axis_from_header(hdul[hdu_name].header, nz)
            if wave is not None:
                return wave
    return wavelength_axis_from_header(hdul[0].header, nz)


def checked_wavelength_mask(wavelengths):
    mask = (wavelengths > WAVE_MIN) & (wavelengths < WAVE_MAX)
    for start, end in AO_GAP_RANGES:
        mask &= ~((wavelengths >= start) & (wavelengths <= end))
    return mask


def mask_data_from_path(mask_path):
    if not mask_path.exists():
        return None
    with fits.open(mask_path, memmap=True) as hdul:
        for hdu in hdul:
            if hdu.data is not None:
                data = np.asarray(hdu.data)
                if data.ndim == 2:
                    return np.array(data, copy=True)
                if data.ndim == 3 and data.shape[0] == 1:
                    return np.array(data[0], copy=True)
    return None


def is_masked_spaxel(mask_data, y, x):
    if mask_data is None:
        return False
    if y >= mask_data.shape[0] or x >= mask_data.shape[1]:
        return True
    value = mask_data[y, x]
    return (not np.isfinite(value)) or value != 0


def checked_data_mask(data_cube, checked_mask, y, x):
    if data_cube is None:
        return checked_mask
    spectrum = np.asarray(data_cube[:, y, x])
    return checked_mask & np.isfinite(spectrum)


def unmasked_spatial_pixels(mask_data, shape):
    ny, nx = shape
    if mask_data is None:
        return np.ones((ny, nx), dtype=bool)
    unmasked = np.zeros((ny, nx), dtype=bool)
    overlap_y = min(ny, mask_data.shape[0])
    overlap_x = min(nx, mask_data.shape[1])
    overlap = np.asarray(mask_data[:overlap_y, :overlap_x])
    unmasked[:overlap_y, :overlap_x] = np.isfinite(overlap) & (overlap == 0)
    return unmasked


def audit_residual_stat(
    data_cube,
    stat_data,
    checked_mask,
    wavelengths=None,
    mask_data=None,
    sample_limit=20,
):
    if data_cube.shape != stat_data.shape:
        raise ValueError(
            f"DATA and STAT cubes must have the same shape, got "
            f"{data_cube.shape} and {stat_data.shape}"
        )
    if checked_mask.shape != (stat_data.shape[0],):
        raise ValueError(
            f"Checked wavelength mask must have shape {(stat_data.shape[0],)}, "
            f"got {checked_mask.shape}"
        )
    if wavelengths is not None and wavelengths.shape != (stat_data.shape[0],):
        raise ValueError(
            f"Wavelength array must have shape {(stat_data.shape[0],)}, "
            f"got {wavelengths.shape}"
        )
    sample_limit = max(0, int(sample_limit))
    unmasked = unmasked_spatial_pixels(mask_data, stat_data.shape[1:])
    count = 0
    samples = []
    for z in np.flatnonzero(checked_mask):
        data_plane = np.asarray(data_cube[z])
        stat_plane = np.asarray(stat_data[z])
        residual = unmasked & np.isfinite(data_plane) & ~positive_finite(stat_plane)
        plane_count = int(np.count_nonzero(residual))
        count += plane_count
        remaining = sample_limit - len(samples)
        if plane_count and remaining > 0:
            for y, x in np.argwhere(residual)[:remaining]:
                wave = None if wavelengths is None else float(wavelengths[z])
                samples.append(
                    ResidualStatSample(
                        int(z),
                        int(y),
                        int(x),
                        wave,
                        float(data_plane[y, x]),
                        float(stat_plane[y, x]),
                    )
                )
    return ResidualStatAudit(count, tuple(samples))


def spatial_circle_offsets():
    limit = int(np.floor(SPATIAL_RADIUS_PIXELS))
    offsets = []
    for dy in range(-limit, limit + 1):
        for dx in range(-limit, limit + 1):
            if dx == 0 and dy == 0:
                continue
            if dx * dx + dy * dy <= SPATIAL_RADIUS_PIXELS**2:
                offsets.append((dy, dx))
    return tuple(offsets)


SPATIAL_CIRCLE_OFFSETS = spatial_circle_offsets()


def spatial_mean_at(cube, z, y, x, mask_data=None, positive_only=False):
    ny, nx = cube.shape[1:]
    values = []
    for dy, dx in SPATIAL_CIRCLE_OFFSETS:
        neighbor_y = y + dy
        neighbor_x = x + dx
        if neighbor_y < 0 or neighbor_y >= ny or neighbor_x < 0 or neighbor_x >= nx:
            continue
        if is_masked_spaxel(mask_data, neighbor_y, neighbor_x):
            continue
        value = cube[z, neighbor_y, neighbor_x]
        if not np.isfinite(value):
            continue
        if positive_only and value <= 0:
            continue
        values.append(float(value))
    if not values:
        return None
    return float(np.mean(values))


def make_spatial_gap(data_cube, stat_data, z, y, x, wavelengths=None, mask_data=None):
    data_fill = spatial_mean_at(data_cube, z, y, x, mask_data)
    stat_fill = spatial_mean_at(stat_data, z, y, x, mask_data, positive_only=True)
    fillable = data_fill is not None and stat_fill is not None
    reason = "bounded_component_spatial_mean_1arcsec" if fillable else "no_valid_spatial_neighbors"
    wave = None if wavelengths is None else float(wavelengths[z])
    return SpatialGap(z, y, x, fillable, reason, wave, data_fill, stat_fill)


def bounded_joint_nan_components(
    data_plane, stat_plane, mask_data=None, unmasked_pixels=None
):
    if data_plane.shape != stat_plane.shape:
        raise ValueError(
            f"DATA and STAT planes must have the same shape, got "
            f"{data_plane.shape} and {stat_plane.shape}"
        )
    ny, nx = data_plane.shape
    if unmasked_pixels is None:
        unmasked = unmasked_spatial_pixels(mask_data, (ny, nx))
    else:
        unmasked = unmasked_pixels
        if unmasked.shape != (ny, nx):
            raise ValueError(
                f"Unmasked-pixel array must have shape {(ny, nx)}, got "
                f"{unmasked.shape}"
            )
    targets = np.isnan(data_plane) & np.isnan(stat_plane) & unmasked
    if not np.any(targets):
        return ()
    labels, component_count = ndimage.label(targets, structure=EIGHT_CONNECTED_STRUCTURE)
    if component_count == 0:
        return ()

    accepted = []
    for label_id, component_slice in enumerate(ndimage.find_objects(labels), start=1):
        if component_slice is None:
            continue
        y_slice, x_slice = component_slice
        if y_slice.start == 0 or x_slice.start == 0 or y_slice.stop == ny or x_slice.stop == nx:
            continue

        padded_y = slice(y_slice.start - 1, y_slice.stop + 1)
        padded_x = slice(x_slice.start - 1, x_slice.stop + 1)
        local_component = labels[padded_y, padded_x] == label_id
        boundary = ndimage.binary_dilation(
            local_component, structure=EIGHT_CONNECTED_STRUCTURE
        ) & ~local_component

        boundary_data = np.asarray(data_plane[padded_y, padded_x])[boundary]
        boundary_stat = np.asarray(stat_plane[padded_y, padded_x])[boundary]
        boundary_unmasked = unmasked[padded_y, padded_x][boundary]
        boundary_valid = (
            boundary_unmasked
            & np.isfinite(boundary_data)
            & np.isfinite(boundary_stat)
            & (boundary_stat > 0)
        )
        if boundary_valid.size == 0 or not np.all(boundary_valid):
            continue

        coordinates = np.argwhere(local_component)
        coordinates[:, 0] += padded_y.start
        coordinates[:, 1] += padded_x.start
        accepted.append(coordinates)
    return tuple(accepted)


def find_spatial_gaps(
    data_cube,
    stat_data,
    checked_mask,
    wavelengths=None,
    mask_data=None,
    z_start=0,
    z_end=None,
):
    if data_cube.shape != stat_data.shape:
        raise ValueError(
            f"DATA and STAT cubes must have the same shape, got "
            f"{data_cube.shape} and {stat_data.shape}"
        )
    if z_end is None:
        z_end = stat_data.shape[0]

    spatial_gaps = []  # type: List[SpatialGap]
    unmasked = unmasked_spatial_pixels(mask_data, stat_data.shape[1:])
    for z in range(z_start, z_end):
        if not checked_mask[z]:
            continue
        components = bounded_joint_nan_components(
            data_cube[z], stat_data[z], mask_data, unmasked
        )
        for coordinates in components:
            component_gaps = [
                make_spatial_gap(
                    data_cube,
                    stat_data,
                    z,
                    int(y),
                    int(x),
                    wavelengths,
                    mask_data,
                )
                for y, x in coordinates
            ]
            if component_gaps and all(gap.fillable for gap in component_gaps):
                spatial_gaps.extend(component_gaps)
    return tuple(spatial_gaps)


def find_runs(mask: np.ndarray) -> Iterable[Tuple[int, int]]:
    indices = np.flatnonzero(mask)
    if indices.size == 0:
        return
    breaks = np.flatnonzero(np.diff(indices) > 1)
    starts = np.r_[indices[0], indices[breaks + 1]]
    ends = np.r_[indices[breaks], indices[-1]]
    for start, end in zip(starts, ends):
        yield int(start), int(end)


def positive_bound_indices(spectrum: np.ndarray, start: int, end: int, checked_mask=None):
    before_idx = start - 1
    after_idx = end + 1
    lower_z = None
    upper_z = None
    if before_idx >= 0 and (checked_mask is None or checked_mask[before_idx]):
        if positive_finite(spectrum[before_idx]):
            lower_z = before_idx
    if after_idx < spectrum.size and (checked_mask is None or checked_mask[after_idx]):
        if positive_finite(spectrum[after_idx]):
            upper_z = after_idx
    return lower_z, upper_z


def classify_gap(spectrum: np.ndarray, start: int, end: int, checked_mask=None) -> Tuple[bool, str]:
    lower_z, upper_z = positive_bound_indices(spectrum, start, end, checked_mask)
    if lower_z is not None and upper_z is not None:
        return True, "bounded_by_positive_values"
    if lower_z is not None:
        return True, "bounded_below_by_positive_value"
    if upper_z is not None:
        return True, "bounded_above_by_positive_value"
    return False, "no_positive_finite_bound"


def find_stat_gaps(
    stat_data: np.ndarray,
    galid: str = "",
    cube_path: Union[Path, str] = "",
    wavelengths=None,
    data_cube=None,
    mask_data=None,
) -> GapReport:
    if stat_data.ndim != 3:
        raise ValueError(f"STAT HDU must be a 3-D cube, got shape {stat_data.shape}")

    if wavelengths is None:
        checked_mask = np.ones(stat_data.shape[0], dtype=bool)
    else:
        checked_mask = checked_wavelength_mask(wavelengths)

    gaps = []  # type: List[StatGap]
    nz, ny, nx = stat_data.shape
    for y in range(ny):
        for x in range(nx):
            if is_masked_spaxel(mask_data, y, x):
                continue
            eligible_mask = checked_data_mask(data_cube, checked_mask, y, x)
            spectrum = np.asarray(stat_data[:, y, x])
            bad = (~positive_finite(spectrum)) & eligible_mask
            for start, end in find_runs(bad):
                fillable, reason = classify_gap(spectrum, start, end, eligible_mask)
                lower_z, upper_z = positive_bound_indices(spectrum, start, end, eligible_mask)
                wave_start = None if wavelengths is None else float(wavelengths[start])
                wave_end = None if wavelengths is None else float(wavelengths[end])
                gaps.append(
                    StatGap(
                        start,
                        end,
                        y,
                        x,
                        fillable,
                        reason,
                        wave_start,
                        wave_end,
                        lower_z,
                        upper_z,
                    )
                )

    spatial_gaps = ()
    if data_cube is not None:
        spatial_gaps = find_spatial_gaps(
            data_cube,
            stat_data,
            checked_mask,
            wavelengths,
            mask_data,
        )

    return GapReport(
        galid,
        Path(cube_path),
        tuple(stat_data.shape),
        tuple(gaps),
        spatial_gaps,
    )


def _scan_y_range(
    args: Tuple[str, Path, int, int]
) -> Tuple[Tuple[int, ...], Tuple[StatGap, ...]]:
    galid, cube_path, y_start, y_end = args
    with fits.open(cube_path, memmap=True) as hdul:
        stat = hdul["STAT"].data
        data_cube = hdul["DATA"].data if "DATA" in hdul else None
        wavelengths = wavelength_axis_from_hdul(hdul, stat.shape[0])
        checked_mask = np.ones(stat.shape[0], dtype=bool) if wavelengths is None else checked_wavelength_mask(wavelengths)
        mask_data = mask_data_from_path(mask_path_for(galid, cube_path.parent))
        gaps = []  # type: List[StatGap]
        nz, ny, nx = stat.shape
        if y_start < 0 or y_end > ny:
            raise ValueError(f"Bad y range {y_start}:{y_end} for STAT shape {stat.shape}")
        for y in range(y_start, y_end):
            for x in range(nx):
                if is_masked_spaxel(mask_data, y, x):
                    continue
                eligible_mask = checked_data_mask(data_cube, checked_mask, y, x)
                spectrum = np.asarray(stat[:, y, x])
                bad = (~positive_finite(spectrum)) & eligible_mask
                for start, end in find_runs(bad):
                    fillable, reason = classify_gap(spectrum, start, end, eligible_mask)
                    lower_z, upper_z = positive_bound_indices(spectrum, start, end, eligible_mask)
                    wave_start = None if wavelengths is None else float(wavelengths[start])
                    wave_end = None if wavelengths is None else float(wavelengths[end])
                    gaps.append(
                        StatGap(
                            start,
                            end,
                            y,
                            x,
                            fillable,
                            reason,
                            wave_start,
                            wave_end,
                            lower_z,
                            upper_z,
                        )
                    )

        return (nz, ny, nx), tuple(gaps)


def _scan_z_range(args):
    galid, cube_path, z_start, z_end = args
    with fits.open(cube_path, memmap=True) as hdul:
        if "DATA" not in hdul:
            return tuple(hdul["STAT"].shape), ()
        stat = hdul["STAT"].data
        data_cube = hdul["DATA"].data
        wavelengths = wavelength_axis_from_hdul(hdul, stat.shape[0])
        checked_mask = (
            np.ones(stat.shape[0], dtype=bool)
            if wavelengths is None
            else checked_wavelength_mask(wavelengths)
        )
        mask_data = mask_data_from_path(mask_path_for(galid, cube_path.parent))
        spatial_gaps = find_spatial_gaps(
            data_cube,
            stat,
            checked_mask,
            wavelengths,
            mask_data,
            z_start,
            z_end,
        )
        return tuple(stat.shape), spatial_gaps


def y_ranges(ny: int, workers: int) -> List[Tuple[int, int]]:
    workers = max(1, min(workers, ny))
    edges = np.linspace(0, ny, workers + 1, dtype=int)
    return [(int(edges[i]), int(edges[i + 1])) for i in range(workers) if edges[i] < edges[i + 1]]


def z_ranges(nz: int, workers: int) -> List[Tuple[int, int]]:
    return y_ranges(nz, workers)


def scan_cube(cube_path: Path, galid: str, workers: int) -> GapReport:
    if workers <= 1:
        return _scan_single(cube_path, galid)

    with fits.open(cube_path, memmap=True) as hdul:
        if "STAT" not in hdul:
            raise KeyError(f"{cube_path} has no STAT HDU")
        shape = tuple(hdul["STAT"].shape)
        has_data = "DATA" in hdul
    if len(shape) != 3:
        raise ValueError(f"STAT HDU must be a 3-D cube, got shape {shape}")

    spatial_ranges = z_ranges(shape[0], workers) if has_data else []
    spectral_ranges = y_ranges(shape[1], workers)
    spectral_tasks = [
        (galid, cube_path, start, end) for start, end in spectral_ranges
    ]
    spatial_tasks = [
        (galid, cube_path, start, end) for start, end in spatial_ranges
    ]
    gaps = []  # type: List[StatGap]
    spatial_gaps = []  # type: List[SpatialGap]
    max_workers = max(len(spectral_tasks), len(spatial_tasks))
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        for chunk_shape, chunk_gaps in executor.map(_scan_y_range, spectral_tasks):
            if tuple(chunk_shape) != shape:
                raise RuntimeError(f"STAT shape changed while scanning {cube_path}")
            gaps.extend(chunk_gaps)
        for chunk_shape, chunk_spatial_gaps in executor.map(_scan_z_range, spatial_tasks):
            if tuple(chunk_shape) != shape:
                raise RuntimeError(f"STAT shape changed while scanning {cube_path}")
            spatial_gaps.extend(chunk_spatial_gaps)
    gaps.sort(key=lambda gap: (gap.y, gap.x, gap.z_start, gap.z_end))
    spatial_gaps.sort(key=lambda gap: (gap.y, gap.x, gap.z))
    return GapReport(galid, cube_path, shape, tuple(gaps), tuple(spatial_gaps))


def fill_value_for_gap(stat_data: np.ndarray, gap: StatGap) -> float:
    bound_indices = [index for index in (gap.lower_z, gap.upper_z) if index is not None]
    if not bound_indices:
        raise ValueError("Cannot fill a gap without a finite positive bound")
    values = np.array([stat_data[index, gap.y, gap.x] for index in bound_indices], dtype=float)
    return float(np.mean(values))


def fill_stat_gaps(stat_data: np.ndarray, gaps: Sequence[StatGap]) -> np.ndarray:
    fixed = np.array(stat_data, copy=True)
    for gap in gaps:
        if not gap.fillable:
            continue
        fixed[gap.z_start : gap.z_end + 1, gap.y, gap.x] = fill_value_for_gap(fixed, gap)
    return fixed


def fill_spatial_gaps(data_cube, stat_data, spatial_gaps):
    fixed_data = np.array(data_cube, copy=True)
    fixed_stat = np.array(stat_data, copy=True)
    for gap in spatial_gaps:
        if not gap.fillable:
            continue
        fixed_data[gap.z, gap.y, gap.x] = gap.data_fill
        fixed_stat[gap.z, gap.y, gap.x] = gap.stat_fill
    return fixed_data, fixed_stat


def fix_cube(input_path: Path, output_path: Path, workers: int = 1, galid: str = "") -> GapReport:
    report = scan_cube(input_path, galid or input_path.stem, workers) if workers > 1 else _scan_single(input_path, galid)
    if report.fillable_count == 0:
        return report

    shutil.copy2(input_path, output_path)
    with fits.open(output_path, mode="update", memmap=True) as hdul:
        stat = hdul["STAT"].data
        data_cube = hdul["DATA"].data if "DATA" in hdul else None
        for gap in report.gaps:
            if not gap.fillable:
                continue
            stat[gap.z_start : gap.z_end + 1, gap.y, gap.x] = fill_value_for_gap(stat, gap)
        if data_cube is not None:
            for gap in report.spatial_gaps:
                if not gap.fillable:
                    continue
                data_cube[gap.z, gap.y, gap.x] = gap.data_fill
                stat[gap.z, gap.y, gap.x] = gap.stat_fill
            wavelengths = wavelength_axis_from_hdul(hdul, stat.shape[0])
            checked_mask = (
                np.ones(stat.shape[0], dtype=bool)
                if wavelengths is None
                else checked_wavelength_mask(wavelengths)
            )
            mask_data = mask_data_from_path(
                mask_path_for(galid or report.galid, input_path.parent)
            )
            report.residual_audit = audit_residual_stat(
                data_cube,
                stat,
                checked_mask,
                wavelengths,
                mask_data,
            )
        hdul.flush()
    return report


def _scan_single(cube_path: Path, galid: str) -> GapReport:
    with fits.open(cube_path, memmap=True) as hdul:
        if "STAT" not in hdul:
            raise KeyError(f"{cube_path} has no STAT HDU")
        stat = hdul["STAT"].data
        data_cube = hdul["DATA"].data if "DATA" in hdul else None
        wavelengths = wavelength_axis_from_hdul(hdul, stat.shape[0])
        mask_data = mask_data_from_path(mask_path_for(galid, cube_path.parent))
        return find_stat_gaps(stat, galid, cube_path, wavelengths, data_cube, mask_data)


def format_gap(galid: str, gap: StatGap) -> str:
    wave_text = "wave=unknown"
    if gap.wave_start is not None and gap.wave_end is not None:
        wave_text = f"wave={gap.wave_start:.2f}:{gap.wave_end:.2f}"
    return (
        f"{galid} (x,y)=({gap.x},{gap.y}) z={gap.z_start}:{gap.z_end} "
        f"{wave_text} length={gap.length} reason={gap.reason}"
    )


def format_spatial_gap(galid, gap):
    wave_text = "wave=unknown" if gap.wave is None else f"wave={gap.wave:.2f}"
    return (
        f"{galid} (x,y)=({gap.x},{gap.y}) z={gap.z} {wave_text} "
        f"method=spatial_1arcsec data_fill={gap.data_fill:.6g} "
        f"stat_fill={gap.stat_fill:.6g} reason={gap.reason}"
    )


def format_residual_audit(galid, audit):
    lines = [
        f"[{galid}] post_fix_residual_finite_DATA_invalid_STAT={audit.count} "
        f"sample_count={len(audit.samples)}"
    ]
    for sample in audit.samples:
        wave_text = "wave=unknown" if sample.wave is None else f"wave={sample.wave:.2f}"
        lines.append(
            f"{galid} post_fix_residual (x,y)=({sample.x},{sample.y}) "
            f"z={sample.z} {wave_text} DATA={sample.data_value:.6g} "
            f"STAT={sample.stat_value:.6g}"
        )
    omitted = audit.count - len(audit.samples)
    if omitted:
        lines.append(f"[{galid}] post_fix_residual_samples_omitted={omitted}")
    return lines


def format_report(report: GapReport) -> List[str]:
    target_gaps = [gap for gap in report.gaps if gap.fillable]
    spatial_targets = [gap for gap in report.spatial_gaps if gap.fillable]
    lines = [
        f"[{report.galid}] cube={report.cube_path}",
        f"[{report.galid}] STAT shape={report.shape} "
        f"spectral_fillable_gaps={len(target_gaps)} "
        f"spatial_fillable_voxels={len(spatial_targets)}",
    ]
    lines.extend(format_gap(report.galid, gap) for gap in target_gaps)
    lines.extend(format_spatial_gap(report.galid, gap) for gap in spatial_targets)
    if not target_gaps and not spatial_targets:
        lines.append(f"[{report.galid}] no fillable spectral or spatial gaps found")
    if report.residual_audit is not None:
        lines.extend(format_residual_audit(report.galid, report.residual_audit))
    return lines


def append_log(log_path: Path, lines: Sequence[str]) -> None:
    with log_path.open("a", encoding="utf-8") as handle:
        for line in lines:
            handle.write(line + "\n")


def start_log(log_path: Path, lines: Sequence[str]):
    backup_path = None
    if log_path.exists() and log_path.stat().st_size > 0:
        timestamp = datetime.now().strftime("%Y%m%dT%H%M%S%f")
        backup_path = log_path.with_name(
            f"{log_path.stem}.{timestamp}{log_path.suffix}"
        )
        log_path.replace(backup_path)
    with log_path.open("w", encoding="utf-8") as handle:
        for line in lines:
            handle.write(line + "\n")
    return backup_path


def check_galaxy(galid: str, cube_dir: Path, workers: int) -> GapReport:
    cube_path = cube_path_for(galid, cube_dir)
    if not cube_path.exists():
        raise FileNotFoundError(cube_path)
    return scan_cube(cube_path, galid, workers)


def check_main(argv=None) -> int:
    args = parse_galaxies(argv)
    cube_dir = Path(args.cube_dir)
    log_path = Path(args.log)
    galaxies = selected_galaxies(args.galaxies)
    header = [
        "",
        f"=== PHANGS STAT check started {datetime.now().isoformat(timespec='seconds')} ===",
        f"cube_dir={cube_dir} workers={args.workers} galaxies={' '.join(galaxies)}",
    ]
    print("\n".join(header))
    backup_log = start_log(log_path, header)
    if backup_log is not None:
        print(f"Rotated previous log to {backup_log}")
    status = 0
    for galid in galaxies:
        try:
            report = check_galaxy(galid, cube_dir, args.workers)
            lines = format_report(report)
        except Exception as exc:
            status = 1
            lines = [f"[{galid}] ERROR {type(exc).__name__}: {exc}"]
        print("\n".join(lines))
        append_log(log_path, lines)
    return status


def fix_main(argv=None) -> int:
    args = parse_galaxies(argv)
    cube_dir = Path(args.cube_dir)
    galaxies = selected_galaxies(args.galaxies)
    status = 0
    for galid in galaxies:
        input_path = cube_path_for(galid, cube_dir)
        output_path = fixed_path_for(galid, cube_dir)
        try:
            if not input_path.exists():
                raise FileNotFoundError(input_path)
            report = fix_cube(input_path, output_path, args.workers, galid)
            report_lines = format_report(report)
            for line in report_lines:
                print(line)
            if report.residual_audit is not None:
                append_log(
                    Path(args.log),
                    format_residual_audit(galid, report.residual_audit),
                )
            if report.fillable_count:
                print(f"[{galid}] wrote {output_path}")
            else:
                print(f"[{galid}] no fillable gaps; fixed cube was not written")
        except Exception as exc:
            status = 1
            print(f"[{galid}] ERROR {type(exc).__name__}: {exc}", file=sys.stderr)
    return status
