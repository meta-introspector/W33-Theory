"""Regression for Pass 11179: best four-qutrit near-perfect ticks -- 2 deficient cuts of rank 2, paired structure."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11179_four_qutrit_near_miss as P  # noqa: E402


def test_near_miss():
    r = P.summarize()
    assert r['count_minimum'] == 2 and r['deficit_minimum'] == 4
    assert r['all_rank_two'] and r['all_three_uniform'] and r['all_third_pairing_maximal']
    assert r['all_ticks_symplectic'] and r['all_single_blocks_invertible']
    assert set(r['example']['two_party_ranks'].values()) == {2}
