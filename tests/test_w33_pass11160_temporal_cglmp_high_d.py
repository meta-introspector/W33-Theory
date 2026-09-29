"""Regression for Pass 11160: temporal CGLMP up to d = 10; gap monotone; method check at d = 3."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11160_temporal_cglmp_high_d as P  # noqa: E402


def test_high_d():
    r = P.summarize()
    assert abs(r['method_check_d3'] - 3.1628065077657) < 1e-9 and r['gaps_monotone']
    assert abs(r['table']['10']['spatial'] - 3.1396) < 1e-3 and r['table']['10']['temporal'] > 3.65
    assert r['temporal_deficit_slope'] < -0.6 and r['spatial_deficit_slope'] > -0.3
