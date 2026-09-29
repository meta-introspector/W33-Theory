"""Regression for Pass 11115: no tachyon-free Wilson-line minimum on the B = 0 axis."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11115_axis_duality_and_window as P  # noqa: E402


def test_axis_window():
    r = P.summarize()
    assert r['fricke']['10']['symmetric'] and not r['fricke']['2']['symmetric']
    assert r['model77_window_start'] == 0.87 and r['model77_first_tachyonic'] == 0.85 and r['model77_monotone']
