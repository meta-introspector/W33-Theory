"""Regression for Pass 11181: exact intrinsically-entangling fractions and the prime rule (GAP class data)."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11181_intrinsic_entanglement as P  # noqa: E402


def test_fractions():
    r = P.summarize()
    assert r['n2']['intrinsically_entangling'] == '19/45' and r['n2']['unique_split'] == '19/48'
    assert r['n3']['intrinsically_entangling'] == '7922/12285' and r['n3']['unique_split'] == '7/48'
    assert r['n2']['orders_fixed_point_free'] == [5, 9]
    for n in ('n2', 'n3'):
        assert r[n]['prime_rule_holds'] and r[n]['local_orders_are_2_3_numbers'] and r[n]['burnside_average_fixed'] == '1'
