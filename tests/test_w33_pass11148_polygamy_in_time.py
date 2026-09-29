"""Regression for Pass 11148: temporal entanglement is polygamous; spatial N_AB + N_AC = 4/sqrt15."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11148_polygamy_in_time as P  # noqa: E402


def test_polygamy():
    r = P.summarize()
    assert r['temporal_AB_AC'] == 2.0 and r['temporal_all_pairs'] == 3.0
    assert r['ab_ac_equals_4_over_sqrt15'] and r['spatial_all_pairs'] < 1.2
