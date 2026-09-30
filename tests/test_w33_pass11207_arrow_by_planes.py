"""Regression for Pass 11207 (part 1): exact intrinsic arrow per conjugacy class (frozen data; live check of the plane formula on
the n = 2 classes and a few n = 3 classes against direct enumeration)."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11207_arrow_by_planes as P  # noqa: E402


def test_frozen():
    d = json.loads((ROOT / "data" / "w33_pass11207_arrow_by_planes.json").read_text())
    for n, vals in ((2, [0, 2]), (3, [0, 2, 3]), (4, [0, 2, 3, 4])):
        s = d[f"n{n}"]
        assert s['total_ok'] and s['all_checks'] and s['plane_formula_matches'] and s['planes_ok']
        assert s['A_values'] == vals and s['A_max'] == n
        assert s['A_below_n_iff_invariant_plane']
    assert d['n2']['A_distribution_psp'] == {'0': 9576, '2': 16344}
    assert d['n3']['A_fractions'] == {'0': '55241/884520', '2': '581/1080', '3': '8836/22113'}
    assert d['n3']['classes'] == 74 and d['n4']['classes'] == 278


def test_live_two_qutrits():
    cls = P.load_classes(2)
    Bs = np.array(P.M.factorisations(2))
    Pl = P.planes_from_splits(Bs)
    assert len(Pl) == 90 and len(cls) == 20
    for c in cls:
        E = P.A.exports(c['S'], Bs, 2).sum(1)
        a, inv, _ = P.arrow_by_planes(c['S'], Pl, 2)
        assert a == int(E.min())
        assert (a < 2) == (inv > 0)
        assert P.is_symplectic(c['S'], 2) and P.M.porder(c['S'], 2) == c['order']
