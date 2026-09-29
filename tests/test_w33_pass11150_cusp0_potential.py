"""Regression for Pass 11150: near cusp 0 the potential tracks its Fricke image and rolls into the disk."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11150_cusp0_potential as P  # noqa: E402


def test_cusp0():
    r = P.summarize()
    assert all(v < 0.01 for v in r['fricke_relative_difference'].values())
    assert r['rolls_into_disk'] and r['slope_toward_disk'] < -1000
