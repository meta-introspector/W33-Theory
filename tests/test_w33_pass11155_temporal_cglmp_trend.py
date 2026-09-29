"""Regression for Pass 11155: temporal CGLMP deficit shrinks with d much faster than the spatial one."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11155_temporal_cglmp_trend as P  # noqa: E402


def test_trend():
    r = P.summarize()
    assert r['gaps_monotone'] and abs(r['table']['6']['spatial'] - 3.0497) < 1e-3 and abs(r['table']['7']['spatial'] - 3.0776) < 1e-3
    assert r['temporal_deficit_loglog_slope'] < -0.5 and r['temporal_deficit_loglog_slope'] < 3 * r['spatial_deficit_loglog_slope']
