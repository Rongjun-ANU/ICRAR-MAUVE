import ast
import unittest
from collections import defaultdict
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from astropy.table import Table


ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "make_ngist_masks_from_catalogs_VRI.py"


def load_query_legacy_sources():
    tree = ast.parse(SCRIPT.read_text(), filename=str(SCRIPT))
    node = next(
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef)
        and node.name == "query_legacy_dr9_tractor_and_photoz"
    )
    units = SimpleNamespace(Quantity=object, deg=object())
    namespace = {
        "SkyCoord": object,
        "u": units,
        "Config": object,
        "np": np,
        "pyvo": None,
        "_as_str": str,
    }
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(SCRIPT), "exec"), namespace)
    return namespace


class _FakeRadius:
    def to_value(self, _unit):
        return 0.1


class _FakeTapResult:
    def __init__(self, table):
        self._table = table

    def to_table(self):
        return self._table


class TestLegacyQueryRetry(unittest.TestCase):
    def test_each_legacy_query_retries_until_it_succeeds(self):
        namespace = load_query_legacy_sources()
        ls_ids = np.arange(1, 502, dtype=int)
        tractor = Table(
            {
                "ls_id": ls_ids,
                "ra": np.full(len(ls_ids), 185.388),
                "dec": np.full(len(ls_ids), 14.609),
                "type": np.full(len(ls_ids), "EXP"),
                "release": np.full(len(ls_ids), 9001),
                "brickid": np.full(len(ls_ids), 42),
                "objid": ls_ids,
                "shape_r": np.full(len(ls_ids), 1.5),
                "shape_e1": np.full(len(ls_ids), 0.1),
                "shape_e2": np.full(len(ls_ids), 0.2),
            }
        )
        photoz = Table(
            {
                "ls_id": ls_ids,
                "z_phot_l95": np.full(len(ls_ids), 0.4),
                "z_phot_u95": np.full(len(ls_ids), 0.8),
                "z_phot_mean": np.full(len(ls_ids), 0.6),
            }
        )

        class EventuallySuccessfulTapService:
            attempts = defaultdict(int)

            def run_sync(self, query, maxrec=None):
                if "FROM ls_dr9.photo_z" in query:
                    if "WHERE ls_id IN (501)" in query:
                        key = "photoz_chunk_2"
                        result = photoz[500:]
                    else:
                        key = "photoz_chunk_1"
                        result = photoz[:500]
                elif "shape_r" in query:
                    key = "tractor_with_shape"
                    result = tractor
                else:
                    key = "tractor_minimal"
                    result = tractor

                self.attempts[key] += 1
                if self.attempts[key] <= 4:
                    raise ConnectionResetError(54, "Connection reset by peer")
                return _FakeTapResult(result)

        service = EventuallySuccessfulTapService()
        namespace["pyvo"] = SimpleNamespace(
            dal=SimpleNamespace(TAPService=lambda _url: service)
        )
        center = SimpleNamespace(
            ra=SimpleNamespace(deg=185.388),
            dec=SimpleNamespace(deg=14.609),
        )
        cfg = SimpleNamespace(enable_legacy=True)

        with patch("time.sleep", return_value=None):
            result = namespace["query_legacy_dr9_tractor_and_photoz"](
                center, _FakeRadius(), cfg
            )

        self.assertIsNotNone(result)
        self.assertEqual(len(result), 501)
        self.assertEqual(service.attempts["tractor_with_shape"], 5)
        self.assertEqual(service.attempts["tractor_minimal"], 0)
        self.assertEqual(service.attempts["photoz_chunk_1"], 5)
        self.assertEqual(service.attempts["photoz_chunk_2"], 5)


if __name__ == "__main__":
    unittest.main()
