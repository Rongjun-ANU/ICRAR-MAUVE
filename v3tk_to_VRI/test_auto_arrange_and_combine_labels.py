import tempfile
import subprocess
import sys
import unittest
from pathlib import Path

import auto_arrange_and_combine as mosaic
from PIL import Image


class LabelHelpersTest(unittest.TestCase):
	def test_parse_label_mode_removes_only_trailing_label(self):
		items, add_labels = mosaic._parse_label_mode(["*_combined_VRI.png", "16", "9", "label"])
		self.assertEqual(items, ["*_combined_VRI.png", "16", "9"])
		self.assertTrue(add_labels)

	def test_parse_label_mode_preserves_existing_command(self):
		items, add_labels = mosaic._parse_label_mode(["*_combined_VRI.png", "16", "9"])
		self.assertEqual(items, ["*_combined_VRI.png", "16", "9"])
		self.assertFalse(add_labels)

	def test_galaxy_id_removes_combined_vri_suffix(self):
		self.assertEqual(
			mosaic._galaxy_id_from_path(Path("NGC4567_8_combined_VRI.png")),
			"NGC4567_8",
		)


class LabelRenderingTest(unittest.TestCase):
	def test_save_canvas_without_labels_preserves_existing_rendering(self):
		image = Image.new("RGBA", (30, 20), (10, 20, 30, 255))
		with tempfile.TemporaryDirectory() as tmpdir:
			output = Path(tmpdir) / "mosaic.png"
			mosaic.save_canvas(output, [image], 40, 30, [(5, 5)])
			with Image.open(output) as rendered:
				self.assertEqual(rendered.getpixel((5, 5)), (10, 20, 30, 255))
				self.assertEqual(rendered.getpixel((0, 0)), (0, 0, 0, 255))

	def test_save_canvas_draws_red_label_inside_image_top_left(self):
		image = Image.new("RGBA", (120, 60), (255, 255, 255, 255))
		with tempfile.TemporaryDirectory() as tmpdir:
			output = Path(tmpdir) / "mosaic.png"
			mosaic.save_canvas(output, [image], 140, 80, [(10, 10)], labels=["NGC4380"])
			with Image.open(output) as rendered:
				red_pixels = [
					(x, y)
					for y in range(10, 40)
					for x in range(10, 100)
					if rendered.getpixel((x, y))[0] > 180
					and rendered.getpixel((x, y))[1] < 100
					and rendered.getpixel((x, y))[2] < 100
				]
		self.assertTrue(red_pixels)


class LabelCliTest(unittest.TestCase):
	def test_trailing_label_renders_galaxy_id_with_ratio(self):
		with tempfile.TemporaryDirectory() as tmpdir:
			workdir = Path(tmpdir)
			for galaxy_id in ("NGC4380", "NGC4501"):
				Image.new("RGBA", (120, 60), (255, 255, 255, 255)).save(
					workdir / f"{galaxy_id}_combined_VRI.png"
				)
			command = [
				sys.executable,
				str(Path(mosaic.__file__).resolve()),
				"--fast",
				"--no-reuse-existing",
				"--no-sync-compatible",
				"*_combined_VRI.png",
				"16",
				"9",
				"label",
			]
			completed = subprocess.run(command, cwd=workdir, text=True, capture_output=True)
			self.assertEqual(completed.returncode, 0, completed.stderr)
			self.assertFalse((workdir / "All_combined_VRI_16_9.png").exists())
			with Image.open(workdir / "All_combined_VRI_label_16_9.png") as rendered:
				self.assertTrue(
					any(
						r > 180 and g < 100 and b < 100
						for r, g, b, _ in rendered.convert("RGBA").getdata()
					)
				)
			self.assertTrue((workdir / "All_combined_VRI_label_16_9.proof.txt").is_file())


if __name__ == "__main__":
	unittest.main()
