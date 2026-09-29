"""Regression for Pass 11140: zero-radius gauge algebra su(16) in all 21 neutral-exit survivors."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11140_zero_radius_u16 as P  # noqa: E402


def test_u16():
    r = P.summarize()
    assert r['n_models'] == 21 and r['su16_in_all'] and r['phase_rules_agree']
