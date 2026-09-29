"""Regression for Pass 11147: temporal CGLMP optimum exceeds the spatial one for d >= 3, gap growing."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11147_temporal_cglmp_by_dimension as P  # noqa: E402


def test_by_dimension():
    r = P.summarize()
    t = r['table']
    assert abs(t['2']['spatial'] - 2 * 2 ** 0.5) < 1e-8 and r['qubit_coincidence']
    assert abs(t['3']['spatial'] - r['d3_spatial_closed_form']) < 1e-8
    assert abs(t['4']['spatial'] - 2.9727) < 1e-3 and abs(t['5']['spatial'] - 3.0157) < 1e-3
    assert r['gap_monotone'] and 3.1627 < t['3']['temporal'] < 3.1629
