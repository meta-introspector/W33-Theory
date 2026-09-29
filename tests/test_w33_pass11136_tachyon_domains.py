"""Regression for Pass 11136: tachyon domains are hyperbolic disks at the Gamma_0(3) special points."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11136_tachyon_domains as P  # noqa: E402


def test_domains():
    r = P.summarize()
    assert all(v['neutral'] < 1e-8 and v['charged'] < 1e-8 for v in r['max_deviation'].values())
    assert r['charged_absent_consistent'] and r['golden_top_is_phi2_over_sqrt3'] and r['sqrt3_horodisk_tangent']
    assert abs(r['golden_hyperbolic_centre'] - r['fricke_point']) < 1e-12
    assert abs(r['golden_eR'] - r['phi_squared']) < 1e-12
    assert abs(r['charged_hyperbolic_centre'] - r['elliptic_point_im']) < 1e-12
    assert abs(r['charged_eR'] - r['two_plus_sqrt3']) < 1e-12
    assert r['min_gap_neutral_over_charged'] > 0.19
