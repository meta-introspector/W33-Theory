"""Regression for Pass 11109: A2 theta-function Yukawas; m_c/m_t = 1/2 at the stabilised T* = rho."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11109_instanton_yukawas_at_rho as P  # noqa: E402


def test_duality_and_ratios():
    r = P.main()
    assert r['duality_closes_double_triangle'] and not r['duality_closes_single_triangle']
    assert abs(r['ratio_at_rho'] - 0.5) < 1e-12 and abs(r['ratio_on_axis']['1.0'] - (3 ** 0.5 - 1) / 2) < 1e-12
    assert abs(r['minimal_triangle_area_fraction'] - 1 / 6) < 1e-12
    assert 3.1 < r['ImT_for_mc_over_mt']['0.0036'] < 3.3
