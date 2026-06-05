import importlib
import os
import unittest

import numpy as np
from astropy.io import fits

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/mplconfig_qc_tests")
os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/xdg_cache_qc_tests")


class QCMassCorrectionTests(unittest.TestCase):
    def test_mass_correction_map_masks_low_snr_and_invalid_ml(self):
        mod = importlib.import_module("QC_ngist_v3tk_v768")

        ebv = np.array([[0.10, 0.20], [0.30, 0.40]])
        snr = np.array([[25.0, 24.9], [30.0, 40.0]])
        ml_r = np.array([[2.0, 2.0], [np.nan, -1.0]])

        sfh_hdul = fits.HDUList(
            [
                fits.PrimaryHDU(),
                fits.ImageHDU(data=ebv, name="EBV"),
                fits.ImageHDU(data=snr, name="SNR_POSTFIT"),
            ]
        )
        mass_hdul = fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(data=ml_r, name="ML_R")])

        result = mod.mass_correction_map(sfh_hdul, mass_hdul, snr_threshold=25.0)

        expected_valid = np.log10(2.0) + mod.calzetti_k(0.623) * 0.10 / 2.5
        self.assertTrue(np.isclose(result[0, 0], expected_valid))
        self.assertTrue(np.isnan(result[0, 1]))
        self.assertTrue(np.isnan(result[1, 0]))
        self.assertTrue(np.isnan(result[1, 1]))


if __name__ == "__main__":
    unittest.main()
