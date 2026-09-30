"""Regression for Pass 11187: certified branch and bound on the torus-symmetric class terminates."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11187_torus_class_bnb as P  # noqa: E402


def test_bnb():
    r = P.summarize()
    assert r['closed_form_checked'] and r['complete']
    assert r['negative_control_not_concave'] and r['negative_control_rejected'] and r['optimum_neighbourhood_accepted']
    assert r['hessian_fd_max_error'] < 1e-4 and r['closed_by_concavity'] > 0
