"""Regression for Pass 11207: A(S) = n - c(S).  Frozen checks plus live constructions for each indecomposable type."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11207_arrow_theorem as T  # noqa: E402

D = json.loads((ROOT / "data" / "w33_pass11207_arrow_theorem.json").read_text())


def test_frozen():
    assert D['constructions_all_optimal']
    for n, k in (('n2', 20), ('n3', 74), ('n4', 278)):
        assert D['classes'][n]['classes'] == k and D['classes'][n]['A_equals_n_minus_c'] == k
        assert D['classes'][n]['jordan_formula_ok'] == k
    for n in ('n5', 'n6'):
        assert D['random'][n]['all_exact'] and D['random'][n]['jordan_formula_ok']


def test_live_constructions():
    assert T.check_alpha(3) == 3 and T.check_alpha(4, other_class=True) == 4
    assert T.check_beta(3) == 3 and T.check_beta(5) == 5
    assert T.check_gamma(1, 3) == (3, 3) and T.check_gamma(2, 1) == (2, 2)
    assert T.check_delta([1, 1, 2]) == (2, 2) and T.check_delta(T.polymul([1, 1, 2], [1, 1, 2])) == (4, 4)


def test_live_random_exact():
    rows = T.exact_arrow_random(5, 6, 7)
    assert all(r['exact'] for r in rows)


def test_live_closed_form():
    import numpy as np
    assert T.arrow(np.eye(8, dtype=np.int64), 4) == 0
    assert T.arrow(T.regular_unipotent(5), 5) == 5
    assert T.arrow(T.W_matrix(3), 3) == 3
