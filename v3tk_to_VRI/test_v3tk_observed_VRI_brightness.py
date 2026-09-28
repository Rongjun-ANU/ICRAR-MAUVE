import inspect
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

import v3tk_combined_VRI_image as combined
import v3tk_get_legacy as legacy
import v3tk_observed_VRI_image as observed


class ColorPreservingBrightnessTest(unittest.TestCase):
	def _match_without_raw_legacy_seam(
		self,
		vri,
		legacy_rgb,
		mask,
		*,
		transition_width,
		target_max_luminance,
		normalization_mask=None,
	):
		match_edge = observed._match_vri_edge_to_legacy
		parameters = inspect.signature(match_edge).parameters
		self.assertNotIn(
			"seam_width",
			parameters,
			"raw Legacy seam blending must not be available inside the VRI footprint",
		)
		kwargs = dict(
			transition_width=transition_width,
			target_max_luminance=target_max_luminance,
		)
		if normalization_mask is not None:
			self.assertIn(
				"normalization_mask",
				parameters,
				"the brightness reference must be restricted to the galaxy-centre aperture",
			)
			kwargs["normalization_mask"] = normalization_mask
		return match_edge(
			vri,
			legacy_rgb,
			mask,
			**kwargs,
		)

	def test_default_gamma_brightens_faint_and_mid_tones(self):
		args = observed._parse_args([])

		self.assertEqual(args.gamma, 0.65)
		self.assertEqual(args.post_boost, 2.75)
		self.assertEqual(args.target_max_luminance, 255.0)
		self.assertEqual(args.center_radius_fraction, 0.22)
		self.assertEqual(args.legacy_transition_width, 20)
		self.assertFalse(hasattr(args, "legacy_seam_width"))

	def test_center_reference_aperture_includes_offset_nucleus_and_excludes_foreground_star(self):
		self.assertTrue(hasattr(observed, "_central_reference_mask"))
		mask = np.ones((101, 101), dtype=bool)

		reference = observed._central_reference_mask(mask, radius_fraction=0.22)

		self.assertTrue(reference[50, 70])  # nucleus at d/min = 0.198
		self.assertFalse(reference[50, 75])  # foreground star at d/min = 0.248

	def test_foreground_star_does_not_set_galaxy_center_brightness_scale(self):
		self.assertTrue(hasattr(observed, "_central_reference_mask"))
		mask = np.ones((101, 101), dtype=bool)
		reference = observed._central_reference_mask(mask, radius_fraction=0.22)
		legacy_rgb = np.zeros((101, 101, 3), dtype=np.float32)
		base_vri = np.full((101, 101, 3), 0.10, dtype=np.float32)
		base_vri[50, 70] = 0.60  # galaxy-centre peak inside the aperture
		base_vri[50, 75] = 0.80  # brighter foreground star outside the aperture
		brighter_star_vri = base_vri.copy()
		brighter_star_vri[50, 75] = 0.95

		base_result = self._match_without_raw_legacy_seam(
			base_vri,
			legacy_rgb,
			mask,
			transition_width=0,
			target_max_luminance=255.0,
			normalization_mask=reference,
		)
		brighter_star_result = self._match_without_raw_legacy_seam(
			brighter_star_vri,
			legacy_rgb,
			mask,
			transition_width=0,
			target_max_luminance=255.0,
			normalization_mask=reference,
		)

		np.testing.assert_allclose(base_result[50, 70], 1.0, atol=1e-7)
		np.testing.assert_array_equal(
			brighter_star_result[50, 70],
			base_result[50, 70],
		)

	def test_tone_curve_lifts_shadows_while_protecting_highlights(self):
		self.assertTrue(hasattr(observed, "_apply_color_preserving_brightness"))
		apply_brightness = observed._apply_color_preserving_brightness
		rgb = np.array(
			[
				[[0.04, 0.02, 0.01], [0.80, 0.40, 0.20]],
				[[0.00, 0.00, 0.00], [0.20, 0.10, 0.05]],
			],
			dtype=np.float32,
		)

		bright = apply_brightness(rgb, gamma=0.65, boost=2.5)
		original_max = rgb.max(axis=-1)
		bright_max = bright.max(axis=-1)
		expected_max = np.clip(
			original_max
			+ 2.5
			* (np.power(original_max, 0.65) - original_max)
			* (1.0 - original_max),
			0.0,
			1.0,
		)

		np.testing.assert_allclose(bright_max, expected_max, rtol=1e-6, atol=1e-7)
		positive = original_max > 0
		np.testing.assert_allclose(
			bright[positive] / bright_max[positive, None],
			rgb[positive] / original_max[positive, None],
			rtol=1e-6,
			atol=1e-7,
		)
		self.assertGreater(bright[0, 0].max(), rgb[0, 0].max())
		self.assertLess(bright[0, 1].max(), 1.0)
		self.assertGreater(
			bright[0, 0].max() / rgb[0, 0].max(),
			bright[1, 1].max() / rgb[1, 1].max(),
		)
		self.assertGreater(
			bright[1, 1].max() / rgb[1, 1].max(),
			bright[0, 1].max() / rgb[0, 1].max(),
		)

	def test_common_luminance_target_is_not_globally_limited_by_one_colored_pixel(self):
		mask = np.ones((3, 3), dtype=bool)
		vri = np.full((3, 3, 3), 0.60, dtype=np.float32)
		vri[0, 0] = np.array([0.80, 0.00, 0.00], dtype=np.float32)
		legacy_rgb = np.zeros_like(vri)

		result = self._match_without_raw_legacy_seam(
			vri,
			legacy_rgb,
			mask,
			transition_width=0,
			target_max_luminance=255.0,
		)

		weights = np.array([0.2126, 0.7152, 0.0722], dtype=np.float32)
		luma = np.tensordot(result, weights, axes=([-1], [0]))
		self.assertAlmostEqual(float(luma[mask].max()), 1.0, places=6)
		self.assertLessEqual(float(result.max()), 1.0)

	def test_legacy_transition_preserves_vri_boundary_and_common_luminance(self):
		self.assertTrue(hasattr(observed, "_match_vri_edge_to_legacy"))

		mask = np.zeros((15, 15), dtype=bool)
		mask[1:14, 1:14] = True
		vri = np.zeros((15, 15, 3), dtype=np.float32)
		vri[mask] = np.array([0.36, 0.32, 0.28], dtype=np.float32)
		boundary_detail = np.array([0.06, 0.03, 0.015], dtype=np.float32)
		vri[1, 7] += boundary_detail
		vri[7, 7] = np.array([0.50, 0.45, 0.40], dtype=np.float32)
		legacy = np.full((15, 15, 3), [0.10, 0.20, 0.30], dtype=np.float32)

		base_vri = vri.copy()
		base_vri[1, 7] -= boundary_detail
		base_result = self._match_without_raw_legacy_seam(
			base_vri,
			legacy,
			mask,
			transition_width=4,
			target_max_luminance=225.0,
		)
		result = self._match_without_raw_legacy_seam(
			vri,
			legacy,
			mask,
			transition_width=4,
			target_max_luminance=225.0,
		)

		# The standalone image remains black outside the VRI footprint.
		np.testing.assert_array_equal(result[~mask], 0.0)
		# VRI structure remains present even at the exact footprint boundary; raw
		# Legacy pixels are never copied or cross-faded into the valid footprint.
		preserved_boundary_detail = result[1, 7] - base_result[1, 7]
		self.assertGreater(
			np.linalg.norm(preserved_boundary_detail),
			0.95 * np.linalg.norm(boundary_detail),
		)
		# Beyond the transition, the interior VRI RGB ratios remain coupled.
		center = result[7, 7]
		np.testing.assert_allclose(center / center[0], vri[7, 7] / vri[7, 7, 0], atol=1e-7)
		# The common display maximum is Rec.709 luma 225/255.
		luma = np.tensordot(result, np.array([0.2126, 0.7152, 0.0722]), axes=([-1], [0]))
		self.assertAlmostEqual(float(luma[mask].max()), 225.0 / 255.0, places=6)

	def test_legacy_transition_preserves_vri_detail_and_ignores_covered_legacy(self):
		self.assertTrue(hasattr(observed, "_match_vri_edge_to_legacy"))
		match_edge = observed._match_vri_edge_to_legacy

		mask = np.zeros((21, 21), dtype=bool)
		mask[1:20, 1:20] = True
		base_vri = np.zeros((21, 21, 3), dtype=np.float32)
		base_vri[mask] = np.array([0.15, 0.12, 0.09], dtype=np.float32)
		base_vri[10, 10] = np.array([0.45, 0.40, 0.35], dtype=np.float32)
		detail_vri = base_vri.copy()
		detail = np.array([0.06, 0.03, 0.015], dtype=np.float32)
		probe = (4, 10)  # Inside the brightness transition.
		detail_vri[probe] += detail
		legacy = np.full((21, 21, 3), [0.25, 0.24, 0.23], dtype=np.float32)
		legacy_changed_under_vri = legacy.copy()
		legacy_changed_under_vri[mask] = np.array([0.95, 0.05, 0.80], dtype=np.float32)
		weights = np.array([0.2126, 0.7152, 0.0722], dtype=np.float32)
		target = float(np.dot(base_vri[10, 10], weights) * 255.0)

		base_result = self._match_without_raw_legacy_seam(
			base_vri,
			legacy,
			mask,
			transition_width=6,
			target_max_luminance=target,
		)
		detail_result = self._match_without_raw_legacy_seam(
			detail_vri,
			legacy,
			mask,
			transition_width=6,
			target_max_luminance=target,
		)
		covered_legacy_result = self._match_without_raw_legacy_seam(
			base_vri,
			legacy_changed_under_vri,
			mask,
			transition_width=6,
			target_max_luminance=target,
		)

		# A sharp VRI feature keeps nearly all of its contrast: it is not alpha-blurred.
		preserved_detail = detail_result[probe] - base_result[probe]
		self.assertGreater(np.linalg.norm(preserved_detail), 0.95 * np.linalg.norm(detail))
		# Legacy pixels covered by the VRI footprint are removed before edge matching,
		# so changing them cannot alter any VRI-footprint output pixel.
		np.testing.assert_array_equal(covered_legacy_result[mask], base_result[mask])

	def test_combiner_replaces_legacy_footprint_instead_of_blending(self):
		self.assertTrue(hasattr(combined, "_replace_legacy_footprint_with_vri"))
		replace_footprint = combined._replace_legacy_footprint_with_vri
		legacy_rgb = np.full((3, 4, 3), [20, 100, 200], dtype=np.uint8)
		vri_rgb = np.full((3, 4, 3), [240, 20, 10], dtype=np.uint8)
		mask = np.zeros((3, 4), dtype=bool)
		mask[1, 1:3] = True

		result = replace_footprint(vri_rgb, legacy_rgb, mask)

		np.testing.assert_array_equal(result[mask], vri_rgb[mask])
		np.testing.assert_array_equal(result[~mask], legacy_rgb[~mask])

	def test_discovery_requires_the_matching_legacy_image(self):
		with TemporaryDirectory() as tmp:
			root = Path(tmp)
			fits_path = root / "NGC4321_PHANGS_DATACUBE_native_VRI.fits"
			fits_path.touch()

			jobs = observed._discover_jobs(root, "*_DATACUBE*_VRI.fits")

		self.assertEqual(len(jobs), 1)
		self.assertEqual(jobs[0].legacy_reprojected_jpg.name, "NGC4321_legacy_reprojected.jpg")

	def test_full_wrapper_runs_legacy_before_observed_before_combined(self):
		wrapper = Path("v3tk_VRI_image.sh").read_text(encoding="utf-8")
		legacy_at = wrapper.index('"${PY_RUN[@]}" v3tk_get_legacy.py')
		observed_at = wrapper.index('"${PY_RUN[@]}" v3tk_observed_VRI_image.py')
		combined_at = wrapper.index('"${PY_RUN[@]}" v3tk_combined_VRI_image.py')

		self.assertLess(legacy_at, observed_at)
		self.assertLess(observed_at, combined_at)

	def test_legacy_stage_does_not_depend_on_an_observed_png(self):
		with TemporaryDirectory() as tmp:
			root = Path(tmp)
			(root / "NGC4321_PHANGS_DATACUBE_native_VRI.fits").touch()

			jobs = legacy._discover_jobs(root, "*_DATACUBE*_VRI.fits")

		self.assertEqual(len(jobs), 1)
		self.assertFalse(hasattr(jobs[0], "observed_png"))


if __name__ == "__main__":
	unittest.main()
