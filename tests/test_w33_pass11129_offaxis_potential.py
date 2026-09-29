"""Regression for Pass 11129: B is a near-flat direction; B = 0 a minimum just above the onset."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11129_offaxis_potential as P  # noqa: E402


def test_offaxis():
    r = P.summarize()
    assert r['cosine_shape'] and r['B0_minimum_near_onset'] and r['B0_maximum_at_2']
    assert r['max_A_over_radial'] < 1e-4
