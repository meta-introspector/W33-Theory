"""Regression for Pass 11170: no F9-linear perfect four-qutrit gate; explicit perfect five-qutrit gate (AME(10,3))."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11170_perfect_gate_ladder as P  # noqa: E402


def test_ladder():
    r = P.summarize()
    assert r['f9_perfect_row_sets'] == {'2': 32, '3': 2048, '4': 0}
    assert r['circulant_perfect'] == {'3': 24, '4': 0, '5': 80}
    assert r['five_qutrit_symplectic'] and r['five_qutrit_cut_entropies'] == [5] and r['five_qutrit_cuts'] == 252
    assert r['three_qutrit_control_cut_entropies'] == [3] and r['ladder']['4'] is False
