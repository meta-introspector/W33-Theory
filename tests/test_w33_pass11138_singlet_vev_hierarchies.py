"""Regression for Pass 11138: random per-singlet VEV hierarchies never split the light generations."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11138_singlet_vev_hierarchies as P  # noqa: E402


def test_vevs():
    r = P.summarize()
    assert r['draws'] == 108 and not r['hierarchy_found']
    assert r['max_down_spread'] < 0.2 and r['max_light_lepton_spread'] < 0.2
