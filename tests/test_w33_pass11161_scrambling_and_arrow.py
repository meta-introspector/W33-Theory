"""Regression for Pass 11161: I3 <= 0 for every two-qutrit Clifford gate; perfect = maximal scramblers; the arrow."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11161_scrambling_and_arrow as P  # noqa: E402


def test_scrambling():
    r = P.summarize()
    assert r['I3_census'] == {'-2': 13824, '-1': 36864, '0': 1152} and r['I3_never_positive']
    assert r['formulas_checked'] and r['depolarising_cases'] > 0 and r['depolarising_max_dev'] < 1e-12
