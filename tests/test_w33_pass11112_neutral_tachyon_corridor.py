"""Regression for Pass 11112: neutral/charged tachyons on different tori; potential shrinks the neutral torus first."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11112_neutral_tachyon_corridor as P  # noqa: E402


def test_corridor():
    r = P.summarize()
    assert r['bfield_degenerate'] and r['bfield_points'] == 36
    assert r['torus_separation'] and r['torus_points'] == 12
    for m in ('model57', 'model53'):
        assert r[m]['neutral_torus_shrinks_first']
        assert all(c['dLambda_strip_cap'] < -0.5 and abs(c['dLambda_tail']) < 0.3 for c in r[m]['comparisons'])
