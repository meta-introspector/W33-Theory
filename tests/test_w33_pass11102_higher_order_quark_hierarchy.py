"""Regression for Pass 11102: in the 12 survivors the non-top generations stay parametrically degenerate at all orders."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11102_higher_order_quark_hierarchy as P  # noqa: E402

C = json.loads(P.OUT.read_text())


def test_reflection_symmetry_at_all_orders():
    assert C["summary"]["higgs_cases"] == 72 and C["summary"]["reflection_symmetric"] == 72
    for v in C["models"].values():
        for s in ("up", "down"):
            for h in v[s]:
                assert P.reflection_symmetric(h)
                e = h["exponents"]["1.0"]
                assert e[0] == 0 and e[1] == e[2] and e[1] > 0      # one heavy state, the other two at the same order
