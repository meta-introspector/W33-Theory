"""Regression for Pass 11164: F9-linear perfect three-qutrit gates = unitaries over F9 with nonzero entries; K = P^2 + i 1."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11164_three_qutrit_tetracode as P  # noqa: E402


def test_f9():
    r = P.summarize()
    assert r['U2_order'] == 96 and r['U2_perfect'] == 64
    assert r['U3_order'] == 24192 and r['U3_perfect'] == 12288 and r['embedding_symplectic']
    assert r['K_shift_plus_i_ones_unitary'] and r['K_column_norms'] == [1, 1, 2]
    assert r['K_superregular_graph_is_MDS_6_3_4_over_F9']
