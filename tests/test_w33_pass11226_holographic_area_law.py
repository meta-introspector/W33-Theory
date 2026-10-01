"""Regression for Pass 11226: exact Ryu-Takayanagi in perfect-tensor networks of the substrate's tick."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def test_perfect_tick_network_recomputed():
    import w33_pass11226_holographic_area_law as H
    Lp, p = H.perfect_tick_choi()
    assert H.is_perfect(Lp, 4)
    T, bonds, net = H.square_hyperbolic_patch()
    r = H.rt_check(T, bonds, net, Lp, H.hyperbolic_square_order)
    assert r["intervals"] == 200 and r["rt_exact"] == 200


def test_frozen():
    d = json.loads((ROOT / "data" / "w33_pass11226_holographic_area_law.json").read_text())
    assert d["rt_exact_everywhere"] and d["ame6_perfect"] and d["perfect_tick_choi_is_AME43"]
    assert d["pentagon_{5,4}"]["rt_exact"] == 300 and d["square_{4,4}"]["rt_exact"] == 210
    c = d["control_SUM_{4,5}"]
    assert not d["control_SUM_choi_is_AME43"] and c["entropy_below_cut"] > 0
