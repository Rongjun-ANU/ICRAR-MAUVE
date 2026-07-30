import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from astropy.io import fits


ROOT = Path(__file__).resolve().parent


def write_source(path: Path) -> None:
    primary = fits.PrimaryHDU()
    primary.header["GALAXY"] = "NGC4254"

    sfr = fits.ImageHDU(
        np.arange(6, dtype=np.float64).reshape(2, 3),
        name="LOGSFR_SURFACE_DENSITY",
    )
    sfr.header["BUNIT"] = "dex(solMass yr-1 kpc-2)"
    sfr.header["CTYPE1"] = "RA---TAN"

    sfr_hii = fits.ImageHDU(
        np.full((2, 3), 2.5, dtype=np.float64),
        name="LOGSFR_SURFACE_DENSITY_HII",
    )
    sfr_hii.header["BUNIT"] = "dex(solMass yr-1 kpc-2)"
    unrelated = fits.ImageHDU(np.ones((2, 3)), name="O_H_D16_HII")

    path.parent.mkdir(parents=True)
    fits.HDUList([primary, unrelated, sfr, sfr_hii]).writeto(path)


def assert_cutout(source: Path, output: Path) -> None:
    with fits.open(source) as source_hdul, fits.open(output) as output_hdul:
        assert [hdu.name for hdu in output_hdul] == [
            "PRIMARY",
            "LOGSFR_SURFACE_DENSITY",
            "LOGSFR_SURFACE_DENSITY_HII",
        ]
        assert output_hdul[0].header == source_hdul[0].header
        for name in ("LOGSFR_SURFACE_DENSITY", "LOGSFR_SURFACE_DENSITY_HII"):
            assert output_hdul[name].header == source_hdul[name].header
            np.testing.assert_array_equal(output_hdul[name].data, source_hdul[name].data)


class SfrCutoutTests(unittest.TestCase):
    def test_python_extracts_only_logsfr_surface_density_hdus_with_headers(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            product_root = Path(tmpdir)
            galaxy_dir = product_root / "v3tk_v7.6.8" / "NGC4254"
            source = galaxy_dir / "NGC4254_gas_bin_maps_further.fits"
            output = galaxy_dir / "NGC4254_SFR_maps_further.fits"
            write_source(source)

            result = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "SFR_cutout.py"),
                    "-g",
                    "NGC4254",
                    "--root",
                    str(product_root),
                    "--product-subdir",
                    "v3tk_v7.6.8",
                ],
                capture_output=True,
                text=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            assert_cutout(source, output)

    def test_shell_accepts_7000_run_selector_and_galaxy(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            product_root = Path(tmpdir)
            galaxy_dir = product_root / "v3tk_v7.6.8_7000" / "NGC4254"
            source = galaxy_dir / "NGC4254_gas_bin_maps_further.fits"
            output = galaxy_dir / "NGC4254_SFR_maps_further.fits"
            write_source(source)

            env = os.environ.copy()
            env["PRODUCT_PARENT"] = str(product_root)
            env["PYTHON_BIN"] = sys.executable
            result = subprocess.run(
                ["bash", str(ROOT / "SFR_cutout.sh"), "7000", "NGC4254"],
                cwd=ROOT,
                env=env,
                capture_output=True,
                text=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            assert_cutout(source, output)


if __name__ == "__main__":
    unittest.main()
