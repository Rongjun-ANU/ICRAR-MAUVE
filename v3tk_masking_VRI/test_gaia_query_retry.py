import ast
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch


ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "make_ngist_masks_from_catalogs_VRI.py"


def load_query_gaia_sources():
    tree = ast.parse(SCRIPT.read_text(), filename=str(SCRIPT))
    node = next(
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == "query_gaia_sources"
    )
    units = SimpleNamespace(Quantity=object, deg=object())
    namespace = {"SkyCoord": object, "u": units, "Config": object, "Gaia": None}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(SCRIPT), "exec"), namespace)
    return namespace


class _FakeRadius:
    def to(self, _unit):
        return SimpleNamespace(value=0.1)


class _FakeJob:
    def __init__(self, result):
        self._result = result

    def get_results(self):
        return self._result


class TestGaiaQueryRetry(unittest.TestCase):
    def test_connection_resets_are_retried_until_the_query_succeeds(self):
        namespace = load_query_gaia_sources()
        expected_result = object()

        class EventuallySuccessfulGaia:
            attempts = 0

            @classmethod
            def launch_job_async(cls, _query, dump_to_file=False):
                self.assertFalse(dump_to_file)
                cls.attempts += 1
                if cls.attempts <= 4:
                    raise ConnectionResetError(54, "Connection reset by peer")
                return _FakeJob(expected_result)

        center = SimpleNamespace(
            ra=SimpleNamespace(deg=185.388),
            dec=SimpleNamespace(deg=14.609),
        )
        cfg = SimpleNamespace(gaia_star_mode="foreground", gaia_gmag_max=21.0)
        namespace["Gaia"] = EventuallySuccessfulGaia

        with patch("time.sleep", return_value=None):
            result = namespace["query_gaia_sources"](center, _FakeRadius(), cfg)

        self.assertIs(result, expected_result)
        self.assertEqual(EventuallySuccessfulGaia.attempts, 5)


if __name__ == "__main__":
    unittest.main()
