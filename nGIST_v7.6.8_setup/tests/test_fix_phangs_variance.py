import importlib.util
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits


ROOT = Path(__file__).resolve().parents[1]
FIX_DIR = ROOT / "fix_phangs"


def load_module():
    spec = importlib.util.spec_from_file_location(
        "phangs_variance_tools", FIX_DIR / "phangs_variance_tools.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def write_cube(path: Path, stat_data: np.ndarray):
    fits.HDUList(
        [
            fits.PrimaryHDU(),
            fits.ImageHDU(data=np.ones_like(stat_data, dtype=np.float32), name="DATA"),
            fits.ImageHDU(data=stat_data.astype(np.float32), name="STAT"),
        ]
    ).writeto(path)


def test_detects_only_non_positive_runs_bounded_by_positive_values():
    tools = load_module()
    stat = np.ones((8, 1, 1), dtype=np.float32)
    stat[3:5, 0, 0] = 0.0
    stat[0, 0, 0] = np.nan
    stat[7, 0, 0] = -1.0

    report = tools.find_stat_gaps(stat)

    assert [(gap.z_start, gap.z_end, gap.y, gap.x, gap.fillable) for gap in report.gaps] == [
        (0, 0, 0, 0, False),
        (3, 4, 0, 0, True),
        (7, 7, 0, 0, False),
    ]


def test_fill_uses_equal_length_windows_before_and_after_gap():
    tools = load_module()
    stat = np.array([10, 20, 0, 0, 50, 70], dtype=np.float32).reshape(6, 1, 1)
    gaps = tools.find_stat_gaps(stat).gaps

    fixed = tools.fill_stat_gaps(stat, gaps)

    assert fixed[:, 0, 0].tolist() == [10, 20, 37.5, 37.5, 50, 70]


def test_fix_cube_writes_fixed_output_and_preserves_original():
    tools = load_module()
    input_path = Path.cwd() / "synthetic_PHANGS_DATACUBE_native.fits"
    output_path = Path.cwd() / "synthetic_PHANGS_DATACUBE_native_fixed.fits"
    stat = np.array([4, 8, np.nan, 16], dtype=np.float32).reshape(4, 1, 1)
    write_cube(input_path, stat)

    try:
        report = tools.fix_cube(input_path, output_path)
        with fits.open(input_path) as original, fits.open(output_path) as fixed:
            assert np.isnan(original["STAT"].data[2, 0, 0])
            assert fixed["STAT"].data[:, 0, 0].tolist() == [4.0, 8.0, 12.0, 16.0]
        assert report.fillable_count == 1
    finally:
        input_path.unlink(missing_ok=True)
        output_path.unlink(missing_ok=True)
