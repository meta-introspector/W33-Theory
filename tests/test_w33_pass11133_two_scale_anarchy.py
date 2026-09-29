"""Regression for Pass 11133: no ratio of VEV scales splits the conjugate-sector generations (anarchy)."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11132_bottom_tau_across_21 as P  # noqa: E402


def test_anarchy():
    r = P.summarize()
    assert r['three_distinct_exponent_cases'] == 0 and r['down_anarchic_every_r']
