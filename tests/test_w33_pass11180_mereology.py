"""Regression for Pass 11180: two-qutrit quantum mereology census over the 45 factorisations."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11180_mereology as P  # noqa: E402


def test_census():
    r = P.summarize()
    assert r['factorisations'] == 45 and r['total'] == 51840
    assert r['intrinsically_entangling'] == 21888 and r['intrinsically_entangling_orders'] == ['5', '9']
    assert r['unique_subsystems'] == 20520 and r['max_perfect'] == 27
    prof = {(x['local'], x['perfect'], x['partial']): x['count'] for x in r['rows']}
    assert prof[(6, 27, 12)] == 480 and prof[(0, 15, 30)] == 10368 and prof[(0, 9, 36)] == 11520
