import importlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw


SCRIPT = Path(__file__).with_name("auto_arrange_and_combine_named.py")


class NamedMosaicTests(unittest.TestCase):
    def module(self):
        self.assertTrue(SCRIPT.is_file(), "The separate named observed-image arranger is missing")
        return importlib.import_module(SCRIPT.stem)

    def test_black_margins_can_overlap_without_covering_faint_pixels(self):
        m = self.module()
        # Complementary L shapes occupy overlapping bounding boxes.
        a = np.zeros((30, 30, 3), dtype=np.uint8)
        a[:5, :, :] = (120, 50, 30)
        a[:, :5, :] = (120, 50, 30)
        a[10, 10] = (1, 0, 0)
        b = np.zeros((30, 30, 3), dtype=np.uint8)
        b[25:, 5:, :] = (30, 50, 120)
        b[5:, 25:, :] = (30, 50, 120)
        canvas = Image.new("RGB", (30, 30), "black")
        occupied = np.zeros((30, 30), dtype=bool)
        for data in (a, b):
            mask = data.max(axis=2) > 0
            m.paste_checked(canvas, occupied, Image.fromarray(data), mask, (0, 0))
        rendered = np.asarray(canvas)
        self.assertEqual(tuple(rendered[10, 10]), (1, 0, 0))
        for data in (a, b):
            mask = data.max(axis=2) > 0
            np.testing.assert_array_equal(rendered[mask], data[mask])
        with self.assertRaisesRegex(ValueError, "overlap"):
            m.paste_checked(canvas, occupied, Image.fromarray(a), a.max(axis=2) > 0, (0, 0))

    def test_conservative_grid_keeps_one_pixel_at_cell_edge(self):
        m = self.module()
        mask = np.zeros((17, 19), dtype=bool)
        mask[7, 7] = mask[16, 18] = True
        coarse = m.coarsen_mask(mask, 8)
        self.assertEqual(coarse.shape, (3, 3))
        self.assertTrue(coarse[0, 0])
        self.assertTrue(coarse[2, 2])
        self.assertEqual(int(coarse.sum()), 2)

    def test_rotation_preserves_quarter_turn_pixels(self):
        m = self.module()
        a = np.arange(7 * 11 * 3, dtype=np.uint8).reshape(7, 11, 3)
        for angle in (0, 90, 180, 270):
            np.testing.assert_array_equal(
                np.asarray(m.rotate_image(Image.fromarray(a), angle)),
                np.rot90(a, angle // 90),
            )

    def test_general_rotations_stay_within_signed_quarter_turn(self):
        m = self.module()
        self.assertEqual(m.general_rotation_angles(30), [-90, -60, -30, 0, 30, 60, 90])
        self.assertEqual(m.general_rotation_angles(90), [-90, 0, 90])
        self.assertTrue(all(-90 <= angle <= 90 for angle in m.general_rotation_angles(15)))

    def test_equivalent_placements_prefer_least_rotation(self):
        m = self.module()
        mask = np.ones((2, 2), dtype=bool)
        variants = [[
            m.Variant(angle, (0, 0, 2, 2), (0, 0), (0, 0, 0, 0), (2, 2), mask)
            for angle in (-90, 0, 90)
        ]]
        result = m.pack_at_size(variants, [0], 4, 4, 1)
        self.assertEqual(result[0].variant.angle, 0)

    def test_smaller_rotation_is_preferred_when_both_variants_fit(self):
        m = self.module()
        variants = [[
            m.Variant(0, (0, 0, 2, 2), (0, 0), (0, 0, 0, 0), (2, 2),
                      np.ones((2, 2), dtype=bool)),
            m.Variant(90, (0, 0, 2, 1), (0, 0), (0, 0, 0, 0), (2, 1),
                      np.ones((1, 2), dtype=bool)),
        ]]
        result = m.pack_at_size(variants, [0], 4, 4, 1)
        self.assertEqual(result[0].variant.angle, 0)

    def test_cli_observed_only_yellow_names_and_verified_layout(self):
        self.module()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name, color in (("NGC1", (110, 40, 20)), ("NGC4567_8", (20, 40, 110))):
                im = Image.new("RGB", (180, 140), "black")
                ImageDraw.Draw(im).polygon([(90, 5), (165, 70), (90, 135), (15, 70)], fill=color)
                im.putpixel((91, 4), (1, 0, 0))
                im.save(root / f"{name}_observed_VRI.png")
            Image.new("RGB", (40, 40), "white").save(root / "NGC999_combined_VRI.png")
            Image.new("RGB", (40, 40), "white").save(root / "All_observed_VRI.png")
            command = [sys.executable, str(SCRIPT), "*_VRI.png", "16", "9",
                       "--font-size", "14", "--grid-size", "4", "--attempts", "2",
                       "--rotation-step", "90"]
            run = subprocess.run(command, cwd=root, text=True, capture_output=True)
            self.assertEqual(run.returncode, 0, run.stderr)
            output = root / "All_observed_VRI_named_16_9.png"
            report = json.loads(output.with_suffix(".layout.json").read_text())
            self.assertEqual({item["label"] for item in report["placements"]}, {"NGC1", "NGC4567_8"})
            self.assertEqual(report["status"], "HEURISTIC_ONLY")
            self.assertEqual(report["validation"]["overlapping_reserved_pixels"], 0)
            self.assertEqual(report["canvas"][0] * 9, report["canvas"][1] * 16)
            with Image.open(output) as im:
                pixels = np.asarray(im)
                self.assertTrue(np.any(np.all(pixels == (255, 255, 0), axis=2)))
                self.assertEqual(int(np.all(pixels == (1, 0, 0), axis=2).sum()), 2)
                for item in report["placements"]:
                    with Image.open(item["source"]) as source:
                        rotated = self.module().rotate_image(source.convert("RGB"), item["angle_degrees"])
                    patch = np.asarray(rotated.crop(item["rotated_crop_box"]))
                    x, y = item["image_position"]
                    region = pixels[y:y + patch.shape[0], x:x + patch.shape[1]]
                    mask = patch.max(axis=2) > 0
                    np.testing.assert_array_equal(region[mask], patch[mask])
            first_report = output.with_suffix(".layout.json").read_bytes()
            run = subprocess.run(command, cwd=root, text=True, capture_output=True)
            self.assertEqual(run.returncode, 0, run.stderr)
            self.assertEqual(first_report, output.with_suffix(".layout.json").read_bytes())

    def test_default_all_and_extended_mode_use_exact_names_and_membership(self):
        m = self.module()
        expected = (
            "IC3392", "NGC4064", "NGC4192", "NGC4293", "NGC4298",
            "NGC4330", "NGC4383", "NGC4396", "NGC4419", "NGC4457",
            "NGC4501", "NGC4522", "NGC4694", "NGC4698",
        )
        self.assertEqual(m.EXTENDED_GALAXY_IDS, expected)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in (*expected, "NGC9999"):
                Image.new("RGB", (30, 24), (40, 20, 10)).save(
                    root / f"{name}_observed_VRI.png"
                )
            common = [
                "--no-rotate", "--no-pa-alignment", "--font-size", "8",
                "--grid-size", "2", "--gap", "1", "--attempts", "1",
            ]
            run = subprocess.run([sys.executable, str(SCRIPT), *common], cwd=root,
                                 text=True, capture_output=True)
            self.assertEqual(run.returncode, 0, run.stderr)
            all_output = root / "All_observed_VRI_named__16_9.png"
            self.assertTrue(all_output.is_file())
            all_report = json.loads(all_output.with_suffix(".layout.json").read_text())
            self.assertEqual(all_report["settings"]["selection_mode"], "all")
            self.assertEqual(len(all_report["placements"]), 15)
            self.assertEqual(all_report["canvas"][0] * 9, all_report["canvas"][1] * 16)

            run = subprocess.run([sys.executable, str(SCRIPT), "extended", *common],
                                 cwd=root, text=True, capture_output=True)
            self.assertEqual(run.returncode, 0, run.stderr)
            extended_output = root / "14_observed_VRI_named_16_9.png"
            self.assertTrue(extended_output.is_file())
            extended_report = json.loads(
                extended_output.with_suffix(".layout.json").read_text()
            )
            self.assertEqual(extended_report["settings"]["selection_mode"], "extended")
            self.assertEqual(
                tuple(item["label"] for item in extended_report["placements"]),
                expected,
            )
            self.assertNotIn("NGC9999", {item["label"] for item in extended_report["placements"]})

    def test_extended_mode_reports_missing_required_galaxies(self):
        self.module()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            Image.new("RGB", (10, 10), "white").save(root / "IC3392_observed_VRI.png")
            run = subprocess.run(
                [sys.executable, str(SCRIPT), "extended", "--no-rotate"],
                cwd=root, text=True, capture_output=True,
            )
            self.assertNotEqual(run.returncode, 0)
            self.assertIn("missing 13 required", run.stderr)
            self.assertIn("NGC4064", run.stderr)

    def test_cli_rejects_no_observed_inputs(self):
        self.module()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            Image.new("RGB", (10, 10), "white").save(root / "NGC1_combined_VRI.png")
            run = subprocess.run([sys.executable, str(SCRIPT), "*_combined_VRI.png"],
                                 cwd=root, text=True, capture_output=True)
            self.assertNotEqual(run.returncode, 0)
            self.assertIn("observed_VRI", run.stderr)
            self.assertFalse(list(root.glob("All*")))

    def test_catalogue_pa_uses_fits_orientation_and_keeps_inclination_distinct(self):
        m = self.module()
        from astropy.io import fits
        from astropy.wcs import WCS

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            table = root / "Brown2021Table1.txt"
            table.write_text("Galaxy\tR.A.\tDecl.\tv_opt\ti\tP.A.\n"
                             "NGC 4192 ^a\t12h\t14deg\t-118\t83\t333\n")
            catalogue = m.read_pa_table(table)
            self.assertEqual(catalogue["NGC4192"], {"pa_deg": 333.0, "inclination_deg": 83.0})
            image_path = root / "NGC4192_observed_VRI.png"
            Image.new("RGB", (40, 30), (100, 50, 20)).save(image_path)
            wcs = WCS(naxis=2)
            wcs.wcs.crpix = [20, 15]
            wcs.wcs.cdelt = [-1 / 3600, 1 / 3600]
            wcs.wcs.crval = [180, 14]
            wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
            fits_path = root / "NGC4192_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk_VRI.fits"
            fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(np.ones((30, 40)), wcs.to_header(), name="V_FLUX")]).writeto(fits_path)
            angles, metadata = m.pa_alignment(m.load_source(image_path), catalogue)
            np.testing.assert_allclose(sorted(angles), [-63, 27], atol=0.002)
            self.assertTrue(all(-90 <= angle <= 90 for angle in angles))
            self.assertEqual(metadata["inclination_deg"], 83)
            # Rotated WCS: do not blindly assume that every PNG has north up.
            wcs.wcs.pc = np.array([[0, -1], [1, 0]])
            fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(np.ones((30, 40)), wcs.to_header(), name="V_FLUX")]).writeto(fits_path, overwrite=True)
            _, metadata = m.pa_alignment(m.load_source(image_path), catalogue)
            self.assertAlmostEqual(metadata["major_axis_angle_from_right_deg"], -27, delta=0.002)

    def test_output_paths_cannot_replace_combined_png_with_report(self):
        self.module()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "NGC1_observed_VRI.png"
            Image.new("RGB", (10, 10), "white").save(source)
            other = root / "NGC1_combined_VRI.png"
            other.write_bytes(b"keep existing product")
            for option in ("--report-file", "--output"):
                run = subprocess.run([sys.executable, str(SCRIPT), str(source), option, str(other)],
                                     cwd=root, text=True, capture_output=True)
                self.assertNotEqual(run.returncode, 0)
                self.assertEqual(other.read_bytes(), b"keep existing product")


if __name__ == "__main__":
    unittest.main()
