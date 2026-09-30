"""Regression for Pass 11207 (part 2): arrow bound beyond four qutrits (frozen search results and semisimplicity census; live
half-moving split for the five-qutrit regular unipotent and a live census check at n = 2)."""
import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11207_half_moving_splits as P  # noqa: E402

D = json.loads((ROOT / "data" / "w33_pass11207_half_moving_splits.json").read_text())


def test_frozen_search():
    for n in ('n5', 'n6'):
        d = D[n]
        assert d['all_symplectic'] and d['all_without_have_half_moving_split'] and d['all_bilagrangian']
        assert d['without'] > 0 and d['rows'][0]['kind'] == 'regular unipotent' and not d['rows'][0]['invariant_plane']
    assert D['n5']['without'] == 136 and D['n5']['irreducible_charpoly'] == 27


def test_frozen_census():
    c = D['semisimple_census']
    assert Fraction(c['n2']['T_semisimple_fraction']) == Fraction(25, 36)
    for n in ('n2', 'n3', 'n4'):
        assert c[n]['S_semisimple_implies_T_semisimple']
        assert 0.67 < float(Fraction(c[n]['T_semisimple_fraction'])) < 0.70


def test_live_unipotent_and_census():
    n = 5
    U = P.regular_unipotent(n)
    pts = P.points(2 * n)
    assert not P.invariant_plane(U, n, pts)
    xs, _ = P.half_moving_split(U, n, pts)
    ok, ex = P.verify_split(U, n, xs)
    assert ok and ex == [1] * n and P.bilagrangian(U, n, xs)
    assert not P.is_semisimple(P.T_of(U, n))            # the regular unipotent is outside the theorem's hypothesis
    assert P.is_semisimple(np.eye(4, dtype=np.int64)) and not P.is_semisimple(P.regular_unipotent(2))
