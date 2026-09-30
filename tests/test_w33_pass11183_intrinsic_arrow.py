"""Regression for Pass 11183: intrinsic arrow -- two qutrits exact (0 or 2 trits), formula check, frozen 3-qutrit data."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11183_intrinsic_arrow as P  # noqa: E402


def test_arrow():
    assert P.check_formula(20)
    t2 = P.two_qutrit()
    A = {}
    for k, v in t2.items():
        for (a, b), c in v.items():
            A[a] = A.get(a, 0) + c
    assert A == {0: 19152, 2: 32688}
    d = json.loads((ROOT / "data" / "w33_pass11183_intrinsic_arrow.json").read_text())
    vals = {int(k.split('|')[0][1:]) for v in d['three_qutrit_sample'].values() for k in v}
    assert vals <= {0, 2, 3} and 1 not in vals
    assert d['paper_ticks']['clock_K']['A'] == 0 and d['paper_ticks']['perfect_tick_VKV']['A'] == 2
    assert d['paper_ticks']['f9_gate']['A'] == 3
