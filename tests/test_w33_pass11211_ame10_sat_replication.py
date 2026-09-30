"""Regression for Pass 11211 (replication of Pass 11190): sign structure of AME(10,3) states.  Frozen SAT verdicts and positive control; live checks
of the local Klein classification and of the sign/pattern dictionary on one AME(10,3) graph."""
import itertools
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11211_ame10_sat_replication as A  # noqa: E402

D = json.loads((ROOT / "data" / "w33_pass11211_ame10_sat_replication.json").read_text())


def test_frozen():
    assert D['local']['labelled_graphs'] == 1052 and D['local']['realisable_classes'] == ['K4+3K1', 'prism+K1']
    pc = D['positive_control']
    assert pc['graphs'] == 71 and pc['zero_sum'] and pc['all_prism'] and pc['duality'] and pc['steiner']
    assert pc['dictionary_ok'] and pc['pattern_counts'] == {'C10': 5112, 'cross+perm': 12780}
    s = D['sat']
    assert s['base'] and not s['Q1_K4_3K1_somewhere'] and not s['Q2_forbidden_pattern']
    assert s['control_C10'] and s['control_cross_perm']
    assert not s['glucose_Q1'] and not s['glucose_Q2']
    assert D['theorem_steiner'] and D['theorem_patterns']


def test_live_local_classification():
    g = A.sign_graphs()
    assert len(g) == 1052
    c = Counter((x['name'], x['realisable']) for x in g)
    assert c == {('prism+K1', True): 840, ('K4+3K1', True): 70, ('K33+K1', False): 140, ('7K1', False): 2}


def test_live_signs_one_graph():
    G = A.SC.to_mat(np.array(json.loads(A.TABU.read_text())['graphs'][5]))
    sg = A.signs(G)
    vals = {U: A.to_value(U, s) for U, s in sg.items()}
    assert all(sum(s) % 3 == 0 for s in sg.values())
    blocks = [tuple(p for p in range(10) if p not in U) for U in A.SIX if vals[U] is None]
    cover = Counter(t for B in blocks for t in itertools.combinations(B, 3))
    assert len(blocks) == 30 and set(cover.values()) == {1} and len(cover) == 120
    assert all(A.same(vals[U], i, j) == A.same(vals[Up], i, j) for U, Up, i, j in A.duality_pairs())
    for I in list(itertools.combinations(range(10), 5))[::21]:
        O = [p for p in range(10) if p not in I]
        P = A.pattern_from_values(I, [vals[tuple(sorted([i] + O))] for i in I])
        S = A.SC.gate_from_graph(G, list(I))
        B = S.reshape(5, 2, 5, 2).transpose(0, 2, 1, 3)
        det = (B[..., 0, 0] * B[..., 1, 1] - B[..., 0, 1] * B[..., 1, 0]) % 3
        assert ((det == 2) == (P == 1)).all()
        assert A.shape(P) == ('cross+perm' if any(tuple(sorted(set(I) - {i})) in blocks for i in I) else 'C10')
