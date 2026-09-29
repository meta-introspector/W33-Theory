"""Regression for Pass 11121: Hesse parameter at the fixed points; the bridge table."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11121_hesse_gauge_bridge as P  # noqa: E402


def test_fixed_points():
    r = P.main()
    assert abs(r['lambda_rho_cubed'][0] - 1) < 1e-12 and r['lambda_i_is_1_plus_sqrt3']
    assert abs(r['jcurve_i'][0] - 1728) < 1e-6
    assert len(r['bridge']) == 7
