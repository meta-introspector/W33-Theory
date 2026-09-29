"""Regression for Pass 11122: the six neutral-exit survivors roll monotonically to the onset along a volume law."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11122_six_potentials as P  # noqa: E402


def test_six_potentials():
    r = P.summarize()
    assert len(r['models']) == 6
    assert r['all_monotone'] and r['volume_law_all'] and r['subvolume_negative_all']
    assert all(v['slope_near_onset'] > 500 for v in r['models'].values())
