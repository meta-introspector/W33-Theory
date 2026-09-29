"""Regression for Pass 11166: the spectral clock tick is a Clifford propagator; translation-invariant clocks are never perfect."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11166_clock_tick_clifford as P  # noqa: E402


def test_clock():
    r = P.summarize()
    assert r['U_cubed_minus_I'] < 1e-10 and r['is_clifford'] and r['symplectic_ok'] and r['shifts_conserved']
    assert r['block_ranks'] == [[2, 0, 1], [0, 2, 0], [1, 0, 2]] and r['column_law'] == [1, 1, 1] and not r['perfect']
    assert r['infos_match_stabilizer_formula']
    assert r['translation_invariant_max_offdiag_rank'] == 1 and r['translation_invariant_perfect'] == 0
