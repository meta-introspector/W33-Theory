"""Regression for Pass 11163: exactly 286654464 = 2^17 3^7 perfect three-qutrit Clifford gates (128/4095 of Sp(6,3))."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11163_perfect_three_qutrit_count as P  # noqa: E402


def test_count():
    r = P.summarize()
    assert r['perfect_count'] == 286654464 == 2 ** 17 * 3 ** 7 and r['fraction_exact'] == '128/4095'
    assert r['count_is_24_times'] and r['within_sample_estimate'] and r['two_qutrit_control'] == 13824
