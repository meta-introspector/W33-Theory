"""Regression for Pass 11132: single-Higgs bottom lowered to eps^2 in exactly the 14 down-unlocked models."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11132_bottom_tau_across_21 as P  # noqa: E402


def test_bottom_tau():
    r = P.summarize()
    assert r['n_models'] == 21 and r['conjugation_bijection'] and r['bottom_eps3_without_T']
    assert len(r['lowered_models']) == 14 and r['lowered_equals_down_unlocked']
    assert all(v >= 2 for v in r['tau_at_eps2_higgs_count'].values())
