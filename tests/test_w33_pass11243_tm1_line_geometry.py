"""Regression for Pass 11243: TM1 as a configuration on a line of W33."""
import itertools
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def test_columns_and_stabilisers():
    import w33_pass11243_tm1_line_geometry as G
    S4 = sorted(itertools.permutations(range(4)))
    cols = G.columns_relative_to_point(S4, 0, True)
    for c in cols:
        target = {"moves X": [2 / 3, 1 / 6, 1 / 6], "fixes X": [0.5, 0.5, 0], "double transposition": [1 / 3] * 3}
        assert np.allclose(c["column"], target[c["kind"]], atol=1e-5)
    dX, dQ = G.D4[0], G.D4[1]
    assert len(G.stabiliser_of(dX, S4)) == 3
    assert [list(p) for p in G.stabiliser_of(dX - dQ, S4)] == [[0, 1, 2, 3], [1, 0, 2, 3]]


def test_frozen():
    d = json.loads((ROOT / "data" / "w33_pass11243_tm1_line_geometry.json").read_text())
    assert d["line_stabiliser"] == dict(sp43_order=51840, line_stabiliser_order=1296, image_order=24, kernel_order=54)
    assert d["tm1_iff_moves_X"]
    assert d["transverse_rank_TM1_pair"] == 2 and d["transverse_rank_aligned_pair"] == 0
    assert d["sequestered_outcomes"] == {"TM1": 48, "theta13 = 0 column": 48}
    assert d["cross_invariants_by_outcome"]["TM1"] == [[16.0, 8.0]]
    assert d["cross_invariants_by_outcome"]["theta13 = 0 column"] == [[0.0, 8.0]]
