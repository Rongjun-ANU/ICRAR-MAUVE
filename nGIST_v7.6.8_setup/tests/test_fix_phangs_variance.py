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


def write_cube(path: Path, stat_data: np.ndarray, data_data=None):
    if data_data is None:
        data_data = np.ones_like(stat_data, dtype=np.float32)
    fits.HDUList(
        [
            fits.PrimaryHDU(),
            fits.ImageHDU(data=np.asarray(data_data, dtype=np.float32), name="DATA"),
            fits.ImageHDU(data=stat_data.astype(np.float32), name="STAT"),
        ]
    ).writeto(path)


def test_detects_non_positive_runs_with_at_least_one_positive_bound():
    tools = load_module()
    stat = np.ones((8, 1, 1), dtype=np.float32)
    stat[3:5, 0, 0] = 0.0
    stat[0, 0, 0] = np.nan
    stat[7, 0, 0] = -1.0

    report = tools.find_stat_gaps(stat)

    assert [(gap.z_start, gap.z_end, gap.y, gap.x, gap.fillable) for gap in report.gaps] == [
        (0, 0, 0, 0, True),
        (3, 4, 0, 0, True),
        (7, 7, 0, 0, True),
    ]


def test_fill_uses_immediate_positive_bounds():
    tools = load_module()
    stat = np.array([10, 20, 0, 0, 50, 70], dtype=np.float32).reshape(6, 1, 1)
    gaps = tools.find_stat_gaps(stat).gaps

    fixed = tools.fill_stat_gaps(stat, gaps)

    assert fixed[:, 0, 0].tolist() == [10, 20, 35.0, 35.0, 50, 70]


def test_fills_multiple_gaps_at_one_spaxel_including_one_sided_gap():
    tools = load_module()
    stat = np.array([10, 0, 20, 0, 0], dtype=np.float32).reshape(5, 1, 1)

    gaps = tools.find_stat_gaps(stat).gaps
    fixed = tools.fill_stat_gaps(stat, gaps)

    assert [(gap.z_start, gap.z_end, gap.fillable) for gap in gaps] == [
        (1, 1, True),
        (3, 4, True),
    ]
    assert fixed[:, 0, 0].tolist() == [10, 15, 20, 20, 20]


def test_checked_wavelength_range_is_open_4700_to_9350_and_excludes_ao_gap():
    tools = load_module()
    wavelengths = np.array(
        [4699.0, 4700.0, 4701.25, 5799.0, 5800.0, 5900.0, 5970.0, 5971.0, 9348.75, 9350.0]
    )

    checked = tools.checked_wavelength_mask(wavelengths)

    assert checked.tolist() == [False, False, True, True, False, False, False, True, True, False]


def test_spatial_components_respect_wavelength_and_ao_exclusions():
    tools = load_module()
    wavelengths = np.array([4700.0, 4701.25, 5900.0, 5971.0, 9350.0])
    data = np.full((5, 7, 7), 2.0, dtype=np.float32)
    stat = np.full((5, 7, 7), 20.0, dtype=np.float32)
    data[:, 3, 3] = np.nan
    stat[:, 3, 3] = np.nan

    report = tools.find_stat_gaps(stat, wavelengths=wavelengths, data_cube=data)

    assert [gap.z for gap in report.spatial_gaps] == [1, 3]


def test_one_sided_gap_does_not_use_excluded_window_neighbor():
    tools = load_module()
    stat = np.array([999, 0, 10], dtype=np.float32).reshape(3, 1, 1)
    wavelengths = np.array([4700.0, 4701.25, 4702.5])

    gaps = tools.find_stat_gaps(stat, wavelengths=wavelengths).gaps
    fixed = tools.fill_stat_gaps(stat, gaps)

    assert [(gap.z_start, gap.z_end, gap.fillable) for gap in gaps] == [(1, 1, True)]
    assert fixed[:, 0, 0].tolist() == [999, 10, 10]


def test_data_nan_at_other_wavelength_does_not_skip_stat_gap():
    tools = load_module()
    stat = np.ones((6, 1, 1), dtype=np.float32)
    stat[1, 0, 0] = np.nan
    data = np.ones_like(stat)
    data[4, 0, 0] = np.nan

    gaps = tools.find_stat_gaps(stat, data_cube=data).gaps

    assert [(gap.z_start, gap.z_end, gap.fillable) for gap in gaps] == [(1, 1, True)]


def test_stat_gap_with_nonpositive_finite_data_is_a_spectral_target():
    tools = load_module()
    for data_value in (-2.0, 0.0):
        stat = np.ones((4, 1, 1), dtype=np.float32)
        stat[1, 0, 0] = np.nan
        data = np.ones_like(stat)
        data[1, 0, 0] = data_value

        report = tools.find_stat_gaps(stat, data_cube=data)

        assert [(gap.z_start, gap.z_end, gap.fillable) for gap in report.gaps] == [
            (1, 1, True)
        ]
        assert report.spatial_gaps == ()


def test_residual_audit_counts_only_checked_unmasked_finite_data_bad_stat():
    tools = load_module()
    wavelengths = np.array([4700.0, 4701.25, 5900.0, 5971.25])
    checked = tools.checked_wavelength_mask(wavelengths)
    data = np.ones((4, 3, 3), dtype=np.float32)
    stat = np.ones_like(data)
    mask = np.zeros((3, 3), dtype=np.uint8)

    stat[0, 0, 0] = np.nan  # Excluded wavelength endpoint.
    data[1, 0, 0] = -2.0
    stat[1, 0, 0] = np.nan  # Count: finite DATA, bad STAT.
    data[1, 0, 1] = np.nan
    stat[1, 0, 1] = np.nan  # Exclude: joint DATA/STAT NaN.
    mask[1, 1] = 1
    stat[1, 1, 1] = np.nan  # Exclude: masked spaxel.
    stat[2, 2, 1] = np.nan  # Exclude: AO gap.
    data[3, 2, 2] = 0.0
    stat[3, 2, 2] = -1.0  # Count: finite DATA, non-positive STAT.

    audit = tools.audit_residual_stat(
        data,
        stat,
        checked,
        wavelengths=wavelengths,
        mask_data=mask,
        sample_limit=10,
    )

    assert audit.count == 2
    assert [(sample.z, sample.y, sample.x) for sample in audit.samples] == [
        (1, 0, 0),
        (3, 2, 2),
    ]


def test_residual_audit_format_is_compact_and_reports_omitted_samples():
    tools = load_module()
    sample = tools.ResidualStatSample(7, 2, 3, 5000.0, -2.0, np.nan)
    audit = tools.ResidualStatAudit(25, (sample,))

    lines = tools.format_residual_audit("NGC4321", audit)

    assert lines == [
        "[NGC4321] post_fix_residual_finite_DATA_invalid_STAT=25 sample_count=1",
        "NGC4321 post_fix_residual (x,y)=(3,2) z=7 wave=5000.00 "
        "DATA=-2 STAT=nan",
        "[NGC4321] post_fix_residual_samples_omitted=24",
    ]


def test_joint_nan_is_excluded_spectrally_but_considered_spatially():
    tools = load_module()
    stat = np.ones((4, 1, 1), dtype=np.float32)
    stat[1, 0, 0] = np.nan
    data = np.ones_like(stat)
    data[1, 0, 0] = np.nan

    report = tools.find_stat_gaps(stat, data_cube=data)

    assert report.gaps == ()
    assert report.spatial_gaps == ()


def test_joint_data_stat_nan_uses_one_arcsec_spatial_means():
    tools = load_module()
    data = np.full((1, 7, 7), 1000.0, dtype=np.float32)
    stat = np.full((1, 7, 7), 2000.0, dtype=np.float32)
    for y in range(7):
        for x in range(7):
            if (x - 3) ** 2 + (y - 3) ** 2 <= 2.5**2:
                data[0, y, x] = 2.0
                stat[0, y, x] = 20.0
    data[0, 3, 5] = 6.0
    stat[0, 3, 5] = 60.0
    data[0, 3, 3] = np.nan
    stat[0, 3, 3] = np.nan

    report = tools.find_stat_gaps(stat, data_cube=data)
    fixed_data, fixed_stat = tools.fill_spatial_gaps(data, stat, report.spatial_gaps)

    assert report.gaps == ()
    assert len(report.spatial_gaps) == 1
    assert report.spatial_gaps[0].fillable
    assert np.isclose(report.spatial_gaps[0].data_fill, 2.2)
    assert np.isclose(report.spatial_gaps[0].stat_fill, 22.0)
    assert np.isclose(fixed_data[0, 3, 3], 2.2)
    assert np.isclose(fixed_stat[0, 3, 3], 22.0)


def test_spatial_component_touching_mask_is_rejected():
    tools = load_module()
    data = np.full((1, 7, 7), 2.0, dtype=np.float32)
    stat = np.full((1, 7, 7), 20.0, dtype=np.float32)
    data[0, 3, 4] = 1000.0
    stat[0, 3, 4] = 1000.0
    data[0, 3, 3] = np.nan
    stat[0, 3, 3] = np.nan
    mask = np.zeros((7, 7), dtype=np.uint8)
    mask[3, 4] = 1

    report = tools.find_stat_gaps(stat, data_cube=data, mask_data=mask)

    assert report.spatial_gaps == ()


def test_spatial_component_with_finite_nonpositive_data_boundary_is_accepted():
    tools = load_module()
    for boundary_value in (-2.0, 0.0):
        data = np.full((1, 7, 7), 2.0, dtype=np.float32)
        stat = np.full((1, 7, 7), 20.0, dtype=np.float32)
        data[0, 3, 2:5] = np.nan
        stat[0, 3, 2:5] = np.nan
        data[0, 3, 5] = boundary_value

        report = tools.find_stat_gaps(stat, data_cube=data)

        assert len(report.spatial_gaps) == 3
        assert {(gap.x, gap.y) for gap in report.spatial_gaps} == {
            (2, 3),
            (3, 3),
            (4, 3),
        }
        assert all(gap.fillable for gap in report.spatial_gaps)


def test_spatial_component_with_nonpositive_stat_boundary_is_rejected():
    tools = load_module()
    data = np.full((1, 7, 7), 2.0, dtype=np.float32)
    stat = np.full((1, 7, 7), 20.0, dtype=np.float32)
    data[0, 3, 3] = np.nan
    stat[0, 3, 3] = np.nan
    stat[0, 2, 2] = 0.0

    report = tools.find_stat_gaps(stat, data_cube=data)

    assert report.spatial_gaps == ()


def test_spatial_component_touching_image_edge_is_rejected():
    tools = load_module()
    data = np.full((1, 7, 7), 2.0, dtype=np.float32)
    stat = np.full((1, 7, 7), 20.0, dtype=np.float32)
    data[0, :, :3] = np.nan
    stat[0, :, :3] = np.nan

    report = tools.find_stat_gaps(stat, data_cube=data)

    assert report.spatial_gaps == ()


def test_diagonally_connected_component_with_invalid_stat_boundary_is_rejected():
    tools = load_module()
    data = np.full((1, 8, 8), 2.0, dtype=np.float32)
    stat = np.full((1, 8, 8), 20.0, dtype=np.float32)
    data[0, 3, 3] = np.nan
    stat[0, 3, 3] = np.nan
    data[0, 4, 4] = np.nan
    stat[0, 4, 4] = np.nan
    stat[0, 2, 2] = 0.0

    report = tools.find_stat_gaps(stat, data_cube=data)

    assert report.spatial_gaps == ()


def test_diagonally_connected_component_with_positive_boundary_is_accepted():
    tools = load_module()
    data = np.full((1, 8, 8), 2.0, dtype=np.float32)
    stat = np.full((1, 8, 8), 20.0, dtype=np.float32)
    data[0, 3, 3] = np.nan
    stat[0, 3, 3] = np.nan
    data[0, 4, 4] = np.nan
    stat[0, 4, 4] = np.nan

    report = tools.find_stat_gaps(stat, data_cube=data)

    assert [(gap.y, gap.x) for gap in report.spatial_gaps] == [(3, 3), (4, 4)]


def test_thick_bounded_component_is_rejected_all_or_none_without_inner_means():
    tools = load_module()
    data = np.full((1, 11, 11), 2.0, dtype=np.float32)
    stat = np.full((1, 11, 11), 20.0, dtype=np.float32)
    data[0, 2:9, 2:9] = np.nan
    stat[0, 2:9, 2:9] = np.nan

    report = tools.find_stat_gaps(stat, data_cube=data)

    assert report.spatial_gaps == ()


def test_fix_cube_fills_joint_nan_in_data_and_stat_spatially():
    tools = load_module()
    input_path = Path.cwd() / "synthetic_spatial_PHANGS_DATACUBE_native.fits"
    output_path = Path.cwd() / "synthetic_spatial_PHANGS_DATACUBE_native_fixed.fits"
    data = np.full((1, 7, 7), 3.0, dtype=np.float32)
    stat = np.full((1, 7, 7), 30.0, dtype=np.float32)
    data[0, 3, 3] = np.nan
    stat[0, 3, 3] = np.nan
    write_cube(input_path, stat, data)

    try:
        report = tools.fix_cube(input_path, output_path)
        with fits.open(input_path) as original, fits.open(output_path) as fixed:
            assert np.isnan(original["DATA"].data[0, 3, 3])
            assert np.isnan(original["STAT"].data[0, 3, 3])
            assert fixed["DATA"].data[0, 3, 3] == 3.0
            assert fixed["STAT"].data[0, 3, 3] == 30.0
        assert report.spatial_fillable_count == 1
    finally:
        input_path.unlink(missing_ok=True)
        output_path.unlink(missing_ok=True)


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


def test_fix_cube_preserves_negative_data_and_reports_no_case1_residual():
    tools = load_module()
    input_path = Path.cwd() / "synthetic_negative_DATA_native.fits"
    output_path = Path.cwd() / "synthetic_negative_DATA_native_fixed.fits"
    stat = np.ones((4, 1, 1), dtype=np.float32)
    stat[1, 0, 0] = np.nan
    data = np.ones_like(stat)
    data[1, 0, 0] = -2.0
    write_cube(input_path, stat, data)

    try:
        report = tools.fix_cube(input_path, output_path)
        with fits.open(output_path) as fixed:
            assert fixed["DATA"].data[1, 0, 0] == -2.0
            assert fixed["STAT"].data[1, 0, 0] == 1.0
        assert report.residual_audit is not None
        assert report.residual_audit.count == 0
    finally:
        input_path.unlink(missing_ok=True)
        output_path.unlink(missing_ok=True)


def test_fix_main_appends_post_fix_audit_to_checker_log():
    tools = load_module()
    cube_dir = Path.cwd() / "synthetic_fix_main_cubes"
    cube_dir.mkdir(exist_ok=True)
    input_path = cube_dir / "NGC4321_PHANGS_DATACUBE_native.fits"
    output_path = cube_dir / "NGC4321_PHANGS_DATACUBE_native_fixed.fits"
    log_path = cube_dir / "check_phangs_variance.log"
    stat = np.ones((4, 1, 1), dtype=np.float32)
    stat[1, 0, 0] = np.nan
    data = np.ones_like(stat)
    data[1, 0, 0] = -2.0
    write_cube(input_path, stat, data)

    try:
        status = tools.fix_main(
            [
                "NGC4321",
                "--cube-dir",
                str(cube_dir),
                "--workers",
                "1",
                "--log",
                str(log_path),
            ]
        )

        assert status == 0
        assert (
            "[NGC4321] post_fix_residual_finite_DATA_invalid_STAT=0 sample_count=0"
            in log_path.read_text(encoding="utf-8")
        )
    finally:
        input_path.unlink(missing_ok=True)
        output_path.unlink(missing_ok=True)
        log_path.unlink(missing_ok=True)
        cube_dir.rmdir()


def test_start_log_replaces_previous_run():
    tools = load_module()
    log_path = Path.cwd() / "synthetic_check_phangs_variance.log"

    try:
        log_path.write_text("stale 12 GB run\n", encoding="utf-8")
        backup_path = tools.start_log(log_path, ["new run"])

        assert log_path.read_text(encoding="utf-8") == "new run\n"
        assert backup_path is not None
        assert backup_path.read_text(encoding="utf-8") == "stale 12 GB run\n"
    finally:
        log_path.unlink(missing_ok=True)
        if "backup_path" in locals() and backup_path is not None:
            backup_path.unlink(missing_ok=True)
