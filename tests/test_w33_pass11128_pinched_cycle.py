"""Regression for Pass 11128: linear winding law down to y = 0.2; golden root = phi^2/sqrt3; no tachyon-free island."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11128_pinched_cycle as P  # noqa: E402


def test_pinched_cycle():
    r = P.summarize()
    assert r['neutral_only_down_to_0p2'] and r['golden_root_is_phi2_over_sqrt3']
    assert len(r['golden_models']) == 4 and len(r['sqrt3_models']) == 2
    assert abs(r['quadratic_residual']) < 1e-12 and r['deepest'] < -0.44
    for c in r['tachyon_counts'].values():
        assert c['0.5'] == 4 and c['0.3'] == 6 and c['0.2'] == 8
