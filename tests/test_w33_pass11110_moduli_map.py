"""Regression for Pass 11110: model-2 one-loop potential over the reachable moduli."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11110_moduli_map as P  # noqa: E402


def test_no_minimum_inside_tachyon_free_region():
    r = P.summarize()
    assert r['lambda_positive_everywhere'] and r['family_min_at_rho'] and r['hierarchy_point_above_rho']
    assert r['lowest_is_at_smallest_radius'] and r['lowest_wilson_line_point']['T_WL'] == ['1.75j', '1.75j']
    assert len(r['points']) == 20
