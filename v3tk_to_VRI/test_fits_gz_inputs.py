import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def test_fits_stem_strips_double_extension():
    import fits_path_utils as utils

    assert utils.fits_stem(Path("NGC4064_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk.fits.gz")) == (
        "NGC4064_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk"
    )


def test_fits_glob_matches_compressed_counterpart():
    import fits_path_utils as utils

    with tempfile.TemporaryDirectory() as tmpdir:
        root = Path(tmpdir)
        compressed = root / "NGC4064_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk_VRI.fits.gz"
        compressed.touch()

        matches = utils.expand_fits_glob(root, "*_DATACUBE*_VRI.fits")

        assert matches == [compressed]


def test_v3tk_to_vri_default_output_from_gzip_input():
    import v3tk_to_VRI as converter

    source = Path("/tmp/NGC4064_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk.fits.gz")

    assert converter.default_output_path(source) == Path(
        "/tmp/NGC4064_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk_VRI.fits"
    )


def test_v3tk_to_vri_default_output_from_phangs_native_input():
    import v3tk_to_VRI as converter

    source = Path("/tmp/NGC4254_PHANGS_DATACUBE_native.fits")

    assert converter.default_output_path(source) == Path(
        "/tmp/NGC4254_PHANGS_DATACUBE_native_VRI.fits"
    )


def test_vri_stage_defaults_discover_phangs_native_products():
    import v3tk_combined_VRI_image as combined
    import v3tk_get_legacy as legacy
    import v3tk_observed_VRI_image as observed

    with tempfile.TemporaryDirectory() as tmpdir:
        root = Path(tmpdir)
        product = root / "NGC4254_PHANGS_DATACUBE_native_VRI.fits"
        product.touch()

        observed_args = observed._parse_args([])
        legacy_args = legacy._parse_args([])
        combined_args = combined._parse_args([])

        observed_jobs = observed._discover_jobs(root, observed_args.pattern)
        legacy_jobs = legacy._discover_jobs(root, legacy_args.pattern)
        combined_jobs = combined._discover_jobs(root, combined_args.pattern, combined_args.suffix)

        assert [job.galaxy_id for job in observed_jobs] == ["NGC4254"]
        assert [job.galaxy_id for job in legacy_jobs] == ["NGC4254"]
        assert [job.galaxy_id for job in combined_jobs] == ["NGC4254"]
        assert observed_jobs[0].output_png.name == "NGC4254_observed_VRI.png"
        assert legacy_jobs[0].legacy_reprojected_jpg.name == "NGC4254_legacy_reprojected.jpg"
        assert combined_jobs[0].output_png.name == "NGC4254_combined_VRI.png"


def test_batch_wrapper_does_not_unzip_inputs():
	script = (ROOT / "v3tk_to_VRI.sh").read_text()

	assert "gunzip" not in script
	assert "dest_fits" not in script
	assert 'v3tk_to_VRI.py "$dest_input"' in script


def test_batch_wrapper_stages_public_phangs_native_inputs():
    script = (ROOT / "v3tk_to_VRI.sh").read_text()

    for galid in ("NGC4254", "NGC4321", "NGC4535"):
        assert f"{galid}_PHANGS_DATACUBE_native.fits" in script
    assert "vos:phangs/RELEASES/PHANGS-MUSE/DR1.0/DATACUBES" in script
    assert "copt" not in script.lower()


if __name__ == "__main__":
    for name, func in sorted(globals().items()):
        if name.startswith("test_") and callable(func):
            func()
            print(f"PASS {name}")
