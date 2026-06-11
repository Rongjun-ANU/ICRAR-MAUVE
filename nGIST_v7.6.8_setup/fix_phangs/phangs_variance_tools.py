#!/usr/bin/env python3
"""Shared PHANGS STAT-HDU variance-gap checks and fixes.

The workflow targets the three PHANGS-native cubes staged in
/scratch/pawsey1308/mauve/cubes/v3tk. Detection is restricted to the useful
4750-9100 A wavelength interval, excludes the AO/LGS gap, ignores spatial
pixels masked by {GALID}_mask.fits, and ignores spaxels whose DATA HDU is
non-finite in the checked non-AO wavelength range. A STAT gap is fillable only
when it is bracketed by finite positive STAT values inside the checked mask.
Fixing copies the input cube and replaces each fillable bad STAT run with the
nanmean of the immediate bracketing STAT values.
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


DEFAULT_GALAXIES = ("NGC4254", "NGC4321", "NGC4535")
DEFAULT_CUBE_DIR = Path("/scratch/pawsey1308/mauve/cubes/v3tk")
LOG_NAME = "check_phangs_variance.log"
WAVE_MIN = 4750.0
WAVE_MAX = 9100.0
AO_GAP_RANGES = ((5800.0, 5970.0),)


class StatGap(object):
    __slots__ = ("z_start", "z_end", "y", "x", "fillable", "reason", "wave_start", "wave_end")

    def __init__(self, z_start, z_end, y, x, fillable, reason, wave_start=None, wave_end=None):
        self.z_start = z_start
        self.z_end = z_end
        self.y = y
        self.x = x
        self.fillable = fillable
        self.reason = reason
        self.wave_start = wave_start
        self.wave_end = wave_end

    @property
    def length(self):
        return self.z_end - self.z_start + 1


class GapReport(object):
    __slots__ = ("galid", "cube_path", "shape", "gaps")

    def __init__(self, galid, cube_path, shape, gaps):
        self.galid = galid
        self.cube_path = cube_path
        self.shape = shape
        self.gaps = gaps

    @property
    def fillable_count(self):
        return sum(1 for gap in self.gaps if gap.fillable)

    @property
    def unfillable_count(self):
        return len(self.gaps) - self.fillable_count


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
    mask = (wavelengths >= WAVE_MIN) & (wavelengths <= WAVE_MAX)
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


def data_has_nonfinite_in_checked_range(data_cube, checked_mask, y, x):
    if data_cube is None:
        return False
    spectrum = np.asarray(data_cube[:, y, x])
    return np.any(~np.isfinite(spectrum[checked_mask]))


def find_runs(mask: np.ndarray) -> Iterable[Tuple[int, int]]:
    indices = np.flatnonzero(mask)
    if indices.size == 0:
        return
    breaks = np.flatnonzero(np.diff(indices) > 1)
    starts = np.r_[indices[0], indices[breaks + 1]]
    ends = np.r_[indices[breaks], indices[-1]]
    for start, end in zip(starts, ends):
        yield int(start), int(end)


def classify_gap(spectrum: np.ndarray, start: int, end: int, checked_mask=None) -> Tuple[bool, str]:
    before_idx = start - 1
    after_idx = end + 1
    if before_idx < 0 or after_idx >= spectrum.size:
        return False, "edge_or_checked_window"
    if checked_mask is not None and (not checked_mask[before_idx] or not checked_mask[after_idx]):
        return False, "edge_or_checked_window"
    if positive_finite(spectrum[before_idx]) and positive_finite(spectrum[after_idx]):
        return True, "bounded_by_positive_values"
    return False, "neighbor_value_not_positive_finite"


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
            if data_has_nonfinite_in_checked_range(data_cube, checked_mask, y, x):
                continue
            spectrum = np.asarray(stat_data[:, y, x])
            bad = (~positive_finite(spectrum)) & checked_mask
            for start, end in find_runs(bad):
                fillable, reason = classify_gap(spectrum, start, end, checked_mask)
                wave_start = None if wavelengths is None else float(wavelengths[start])
                wave_end = None if wavelengths is None else float(wavelengths[end])
                gaps.append(StatGap(start, end, y, x, fillable, reason, wave_start, wave_end))

    return GapReport(galid, Path(cube_path), tuple(stat_data.shape), tuple(gaps))


def _scan_y_range(args: Tuple[str, Path, int, int]) -> Tuple[Tuple[int, ...], Tuple[StatGap, ...]]:
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
                if data_has_nonfinite_in_checked_range(data_cube, checked_mask, y, x):
                    continue
                spectrum = np.asarray(stat[:, y, x])
                bad = (~positive_finite(spectrum)) & checked_mask
                for start, end in find_runs(bad):
                    fillable, reason = classify_gap(spectrum, start, end, checked_mask)
                    wave_start = None if wavelengths is None else float(wavelengths[start])
                    wave_end = None if wavelengths is None else float(wavelengths[end])
                    gaps.append(StatGap(start, end, y, x, fillable, reason, wave_start, wave_end))
        return (nz, ny, nx), tuple(gaps)


def y_ranges(ny: int, workers: int) -> List[Tuple[int, int]]:
    workers = max(1, min(workers, ny))
    edges = np.linspace(0, ny, workers + 1, dtype=int)
    return [(int(edges[i]), int(edges[i + 1])) for i in range(workers) if edges[i] < edges[i + 1]]


def scan_cube(cube_path: Path, galid: str, workers: int) -> GapReport:
    if workers <= 1:
        return _scan_single(cube_path, galid)

    with fits.open(cube_path, memmap=True) as hdul:
        if "STAT" not in hdul:
            raise KeyError(f"{cube_path} has no STAT HDU")
        shape = tuple(hdul["STAT"].shape)
    if len(shape) != 3:
        raise ValueError(f"STAT HDU must be a 3-D cube, got shape {shape}")

    ranges = y_ranges(shape[1], workers)
    tasks = [(galid, cube_path, start, end) for start, end in ranges]
    gaps = []  # type: List[StatGap]
    with concurrent.futures.ProcessPoolExecutor(max_workers=len(tasks)) as executor:
        for chunk_shape, chunk_gaps in executor.map(_scan_y_range, tasks):
            if tuple(chunk_shape) != shape:
                raise RuntimeError(f"STAT shape changed while scanning {cube_path}")
            gaps.extend(chunk_gaps)
    gaps.sort(key=lambda gap: (gap.y, gap.x, gap.z_start, gap.z_end))
    return GapReport(galid, cube_path, shape, tuple(gaps))


def fill_value_for_gap(stat_data: np.ndarray, gap: StatGap) -> float:
    values = np.array(
        [
            stat_data[gap.z_start - 1, gap.y, gap.x],
            stat_data[gap.z_end + 1, gap.y, gap.x],
        ],
        dtype=float,
    )
    return float(np.nanmean(values))


def fill_stat_gaps(stat_data: np.ndarray, gaps: Sequence[StatGap]) -> np.ndarray:
    fixed = np.array(stat_data, copy=True)
    for gap in gaps:
        if not gap.fillable:
            continue
        fixed[gap.z_start : gap.z_end + 1, gap.y, gap.x] = fill_value_for_gap(fixed, gap)
    return fixed


def fix_cube(input_path: Path, output_path: Path, workers: int = 1, galid: str = "") -> GapReport:
    report = scan_cube(input_path, galid or input_path.stem, workers) if workers > 1 else _scan_single(input_path, galid)
    if report.fillable_count == 0:
        return report

    shutil.copy2(input_path, output_path)
    with fits.open(output_path, mode="update", memmap=True) as hdul:
        stat = hdul["STAT"].data
        for gap in report.gaps:
            if not gap.fillable:
                continue
            stat[gap.z_start : gap.z_end + 1, gap.y, gap.x] = fill_value_for_gap(stat, gap)
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


def format_report(report: GapReport) -> List[str]:
    target_gaps = [gap for gap in report.gaps if gap.fillable]
    lines = [
        f"[{report.galid}] cube={report.cube_path}",
        f"[{report.galid}] STAT shape={report.shape} fillable_gaps={len(target_gaps)}",
    ]
    lines.extend(format_gap(report.galid, gap) for gap in target_gaps)
    if not target_gaps:
        lines.append(f"[{report.galid}] no fillable non-positive or non-finite STAT gaps found")
    return lines


def append_log(log_path: Path, lines: Sequence[str]) -> None:
    with log_path.open("a", encoding="utf-8") as handle:
        for line in lines:
            handle.write(line + "\n")


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
    append_log(log_path, header)
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
            for line in format_report(report):
                print(line)
            if report.fillable_count:
                print(f"[{galid}] wrote {output_path}")
            else:
                print(f"[{galid}] no fillable gaps; fixed cube was not written")
        except Exception as exc:
            status = 1
            print(f"[{galid}] ERROR {type(exc).__name__}: {exc}", file=sys.stderr)
    return status
