"""Regression for Pass 11162: 4/sqrt15 is a strict local maximum modulo local unitaries; 96 restarts find nothing higher."""
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11162_polygamy_local_maximum as P  # noqa: E402


def test_local_maximum():
    r = P.summarize()
    assert abs(r['f_at_psi'] - 4 / math.sqrt(15)) < 1e-12 and r['grad_max'] < 1e-6
    assert abs(r['pt_gap_from_zero'] - 2 / 15) < 1e-9
    assert r['hessian_negative'] == 32 and r['hessian_positive'] == 0
    assert r['hessian_zero'] == r['lu_orbit_dim'] == 20 and r['hessian_least_negative'] < -0.7
    assert r['restarts'] == 96 and r['restarts_above_target'] == 0 and r['restarts_at_target'] == 96
