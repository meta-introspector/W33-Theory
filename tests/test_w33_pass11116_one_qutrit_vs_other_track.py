"""Regression for Pass 11116: twist = gauge element 104/104; FI Z3 in 3; no SU(3) gluing (243) in any model."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11116_one_qutrit_vs_other_track as P  # noqa: E402


def test_classification():
    r = P.summarize()
    assert r['models'] == 104 and r['twist_gauge_with_centres'] == 104 and r['twist_gauge_u1_only'] == 97
    assert r['twist_is_fi_z3'] == ['4', '25', '29'] and r['glueable_su3_models'] == []
    assert r['survivors_twist_gauge']
