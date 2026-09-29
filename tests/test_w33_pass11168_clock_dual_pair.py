"""Regression for Pass 11168: no Lorentz-invariant clock dynamics is perfect; perfect kick-tick-kick steps all have pi = v0<->v2."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11168_clock_dual_pair as P  # noqa: E402


def test_clock_dual_pair():
    r = P.summarize()
    assert r['clock_is_kinetic_with_gram_M'] and r['order_K_V'] == 24 and r['perfect_in_K_V'] == 0
    assert r['transverse_blocks_vanish_in_K_V'] and r['lorentz_order'] == 48
    assert r['order_K_V_Lorentz'] == 576 and r['perfect_in_K_V_Lorentz'] == 0
    assert r['K_V_perfect'] == 0 and r['K_V_K_perfect'] == 0 and r['V_K_V_perfect'] == 108
    assert r['V_K_V_pis'] == [[2, 1, 0]] and r['V_K_V_all_lorentz_breaking'] and r['determinant_formulas_hold']
    assert r['K_V_K_V_perfect'] == 30 and r['V1_K_V2_perfect_pairs'] == 23328
