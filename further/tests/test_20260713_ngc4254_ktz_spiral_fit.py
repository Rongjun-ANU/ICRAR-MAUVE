import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "20260713_NGC4254_KTZ_spiral_fit.ipynb"


def load_notebook():
    assert NOTEBOOK.exists(), f"Notebook has not been created: {NOTEBOOK}"
    return json.loads(NOTEBOOK.read_text(encoding="utf-8"))


class NotebookContractTests(unittest.TestCase):
    def test_notebook_has_explanation_before_every_code_cell(self):
        cells = load_notebook()["cells"]
        self.assertEqual(cells[0]["cell_type"], "markdown")
        for index, cell in enumerate(cells):
            if cell["cell_type"] == "code":
                self.assertGreater(index, 0)
                self.assertEqual(cells[index - 1]["cell_type"], "markdown")
                explanation = "".join(cells[index - 1]["source"])
                self.assertGreaterEqual(len(explanation.split()), 35)

    def test_notebook_declares_required_sections_and_parameters(self):
        nb = load_notebook()
        text = "\n".join("".join(cell["source"]) for cell in nb["cells"])
        required = [
            "LOGSFR_SURFACE_DENSITY_HII",
            "BIN_ID",
            "m_arms",
            "pitch_angle",
            "Theta0",
            "h_R",
            "eta",
            "harmonic_n",
            "harmonic_g",
            "harmonic_alpha",
            "Synthetic injection-recovery",
            "Spiral skeleton",
            "Arm profile",
            "Interpretation limits",
        ]
        for item in required:
            self.assertIn(item, text)

    def test_executed_notebook_has_saved_result_contract(self):
        nb = load_notebook()
        code_cells = [cell for cell in nb["cells"] if cell["cell_type"] == "code"]
        self.assertTrue(code_cells)
        self.assertTrue(all(cell["execution_count"] is not None for cell in code_cells))
        text_outputs = []
        image_count = 0
        for cell in code_cells:
            for output in cell.get("outputs", []):
                if output.get("output_type") == "stream":
                    text = output.get("text", [])
                    text_outputs.extend([text] if isinstance(text, str) else text)
                data = output.get("data", {})
                image_count += int("image/png" in data)
        joined = "".join(text_outputs)
        self.assertIn("SYNTHETIC_RECOVERY_PASS", joined)
        self.assertIn("REAL_FIT_COMPLETE", joined)
        self.assertIn("POSITIVITY_CHECK_PASS", joined)
        self.assertGreaterEqual(image_count, 6)


if __name__ == "__main__":
    unittest.main()
