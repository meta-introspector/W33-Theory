"""Regression for Pass 11159: N_AB(x) = sqrt((1+15x)(1-x))/4, maximum 4/sqrt15 at x = 7/15."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11159_explicit_polygamy_state as P  # noqa: E402


def test_explicit():
    r = P.summarize()
    assert r['closed_form_matches_symbolic'] and r['critical_point'] == ['7/15'] and r['max_value'] == '4*sqrt(15)/15'
    assert abs(r['explicit_state_NAB'] - 2 / 15 ** 0.5) < 1e-12 and abs(r['explicit_state_NBC'] - 1 / 15) < 1e-12
