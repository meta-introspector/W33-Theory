"""Regression for Pass 11186: AME(10,3) graph gates realise the same two patterns; determinant law of every order."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11186_five_qutrit_beyond_f9 as P  # noqa: E402


def test_beyond_f9():
    r = P.summarize()
    assert r['determinant_law_every_order'] and r['jacobi_complement_symmetry']
    assert r['graphs_found'] == 71 and r['gates_all_perfect_symplectic']
    assert r['patterns'] == {'C10': 5112, 'cross+perm': 12780} and r['new_beyond_f9'] == []
