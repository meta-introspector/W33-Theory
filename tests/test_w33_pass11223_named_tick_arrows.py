"""Regression for Pass 11223: arrows of the named ticks (recomputed, ~10 s)."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def test_recomputed():
    import w33_pass11223_named_tick_arrows as N
    r = N.run()
    assert r["group_order"] == 51840 and r["arrow_free_total"] == "133/360"
    b = r["by_relation"]
    assert (b["local"]["A0"], b["perfect"]["A0"], b["partial"]["A0"]) == (816, 4896, 13440)
    assert b["perfect"]["p_arrow_free"] == "17/48"
    assert r["named_two_qutrit"]["p"]["A"] == 0 and r["named_two_qutrit"]["h_star"]["A"] == 2
    assert r["transvections"]["A"] == {0: 80}
    t3 = r["named_three_qutrit"]
    assert t3["perfect_tick_VKV"]["A"] == 2 and t3["f9_gate"]["A"] == 3 and t3["clock_K"]["A"] == 0
