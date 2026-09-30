"""Regression for Pass 11175: five admissible patterns; F9-linear census 2 642 411 520 realising two of them."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11175_five_qutrit_patterns as P  # noqa: E402
import w33_pass11175_scan_five_qutrit as S  # noqa: E402


def test_patterns():
    r = P.summarize()
    assert r['admissible_counts'] == {'double cross': 100, 'C4+C6': 600, 'cross+perm': 600, 'C10': 1440, 'all': 1}
    assert r['sp10_samples'] == 7000000 and r['sp10_perfect'] == 0
    assert r['f9_unitaries'] == 2642411520 and r['f9_patterns'] == {'cross+perm': 15728640, 'C10': 6291456}
    assert r['realised'] == ['C10', 'cross+perm'] and r['circulant_stabiliser']['stabiliser'] == 4


def test_f9_controls():
    assert S.f9_perfect(3)['row_sets'] == 2048 and S.f9_perfect(4)['row_sets'] == 0
