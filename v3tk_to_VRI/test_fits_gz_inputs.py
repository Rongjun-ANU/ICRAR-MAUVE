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


def test_batch_wrapper_does_not_unzip_inputs():
	script = (ROOT / "v3tk_to_VRI.sh").read_text()

	assert "gunzip" not in script
	assert "dest_fits" not in script
	assert 'v3tk_to_VRI.py "$dest_input"' in script


if __name__ == "__main__":
    for name, func in sorted(globals().items()):
        if name.startswith("test_") and callable(func):
            func()
            print(f"PASS {name}")
