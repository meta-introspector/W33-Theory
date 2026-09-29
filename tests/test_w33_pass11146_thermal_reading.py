"""Regression for Pass 11146: no Witten-sector tachyon without holonomy; TFD partial transpose = Euclidean Choi."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11146_thermal_reading as P  # noqa: E402


def test_thermal():
    r = P.summarize()
    assert r['n_models'] == 21 and r['no_tachyon_without_holonomy'] and r['tachyonic_with_holonomy_at_small_radius']
    for row in r['tfd']:
        assert row['pt_equals_jamiolkowski_half_evolution'] < 1e-12 and abs(row['negativity'] - row['formula']) < 1e-12
    assert abs(r['tfd'][0]['negativity'] - 1) < 1e-12
