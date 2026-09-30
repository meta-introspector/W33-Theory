"""Regression for Pass 11173: explicit strategy with I_d = 4 - 2 L(d)/(d-1), L(d) = 3 + 2(s-2)/(s(s-1)), s = sqrt d."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11173_temporal_cglmp_to_four as P  # noqa: E402


def test_to_four():
    r = P.summarize()
    assert r['exact_matches_formula'] and r['formula_matches_float'] and r['ramp_matches_direct']
    assert abs(r['I_1024'] - 3.9940166) < 1e-6
    assert r['exact_L_perfect_squares']['16'] == '10/3'
