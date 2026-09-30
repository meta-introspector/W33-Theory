#!/usr/bin/env python3
"""Pass 11175: perfect five-qutrit Clifford gates -- their possible orientation patterns, their rarity, and which patterns
occur.

PATTERN THEOREM.  For a perfect S in Sp(10,3) every block determinant is +-1 and every column and row sums to 1 mod 3
(Pass 11158), so each column holds k entries +1 and 5 - k entries -1 with 10 - k = 1 mod 3: k = 0 or 3 -- five or two
orientation reversals -- and likewise each row.  Up to relabelling inputs and outputs there are exactly five such
patterns of -1's (checked by exhaustive enumeration of the 11^5 row choices):
    'C10'   every row and column has two -1's, forming one 10-cycle of the bipartite graph;
    'C4+C6' two -1's per row and column, forming a 4-cycle and a 6-cycle;
    'cross+perm'  one full row and one full column of -1's plus a permutation on the remaining 4 x 4;
    'double cross' two full rows and two full columns;
    'all'   every block orientation-reversing.
RARITY.  Perfect <=> all 25 blocks and all 100 two-party 4 x 4 minors invertible (the 3- and 4-party minors follow from
det S[I,J] = +-det S[I^c,J^c], Jacobi's identity for symplectic S).  7 000 000 samples of Sp(10,3) (products of 200
random transvections; analysis/w33_pass11175_scan_sp10_3.py) contain NO perfect gate: the fraction is below 4.3e-7
(95%), against 128/4095 = 0.031 for three qutrits and 4/15 for two.
One local orbit is far smaller than that: the explicit circulant of Pass 11170 has local stabiliser of order 4, so its
local orbit (24^10/4 = 1.59e13 gates) is a fraction 1.0e-13 of Sp(10,3) -- invisible to sampling.
OCCURRENCE.  analysis/w33_pass11175_scan_five_qutrit.py:
  (a) F9-linear, exhaustively (orbit-stabiliser over the column monomial unitaries, whose orbits on unit rows are
      labelled by the number of norm-1 entries; controls: 2048 row sets for 3 x 3, 0 for 4 x 4):
      22 020 096 row sets = 2 642 411 520 superregular unitary 5 x 5 matrices over F9 = 2^23 3^2 5 7, a fraction
      0.00256 of U(5,F9) (compare 2/3, 32/63, 0 for n = 2, 3, 4).  They realise only TWO of the five patterns:
      'cross+perm' (15 728 640 row sets) and 'C10' (6 291 456).
  (b) AME(10,3) weighted graph states by hill climbing (each would give a perfect gate for every choice of five input
      legs): NONE found in 6 x 25 minutes -- they exist (Danielsen, arXiv:1106.2428) but are not reachable this way.
  So the F9-linear gates realise 'C10' and 'cross+perm'; whether 'C4+C6', 'double cross' or 'all' occur is OPEN.
"""
from __future__ import annotations

import itertools
import json
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SCAN = ROOT / "data" / "w33_pass11175_five_qutrit_scan.json"
SAMPLE = ROOT / "data" / "w33_pass11175_sp10_3_sample.json"
OUT = ROOT / "data" / "w33_pass11175_five_qutrit_patterns.json"


def canon(P):
    n = P.shape[0]
    best = None
    for r in itertools.permutations(range(n)):
        Q = P[list(r)]
        cols = sorted(tuple(Q[:, j]) for j in range(n))
        key = tuple(itertools.chain(*zip(*cols)))
        if best is None or key < best:
            best = key
    return ''.join(map(str, best))


def name(key):
    P = np.array([int(c) for c in key]).reshape(5, 5)
    full_rows = int((P.sum(1) == 5).sum())
    if P.sum() == 25:
        return 'all'
    if full_rows == 2:
        return 'double cross'
    if full_rows == 1:
        return 'cross+perm'
    # 2-regular bipartite: count cycles
    seen, cycles = set(), []
    for r0 in range(5):
        if ('r', r0) in seen:
            continue
        length, node, prev = 0, ('r', r0), None
        while True:
            seen.add(node)
            kind, i = node
            nbrs = [('c', j) for j in range(5) if P[i, j]] if kind == 'r' else [('r', j) for j in range(5) if P[j, i]]
            nxt = [x for x in nbrs if x != prev]
            prev, node = node, nxt[0]
            length += 1
            if node == ('r', r0):
                break
        cycles.append(length)
    return 'C10' if cycles == [10] else 'C4+C6'


def admissible_patterns():
    rows = [tuple(1 if k in s else 0 for k in range(5)) for s in itertools.combinations(range(5), 2)] + [(1,) * 5]
    classes = Counter()
    for combo in itertools.product(rows, repeat=5):
        P = np.array(combo)
        if all(P[:, j].sum() in (2, 5) for j in range(5)):
            classes[canon(P)] += 1
    return classes


def summarize():
    adm = admissible_patterns()
    names = {k: name(k) for k in adm}
    res = dict(pass_id=11175, admissible_classes=len(adm), admissible_names=sorted(set(names.values())),
               admissible_counts={names[k]: v for k, v in adm.items()})
    if SAMPLE.exists():
        s = json.loads(SAMPLE.read_text())
        res.update(sp10_samples=s['samples'], sp10_perfect=s['perfect'],
                   sp10_fraction_upper_95=3.0 / s['samples'] if s['perfect'] == 0 else None)
    if SCAN.exists():
        d = json.loads(SCAN.read_text())
        f9 = d['f9']
        res.update(f9_rows=f9['rows'], f9_row_sets=f9['row_sets'], f9_unitaries=f9['unitary_matrices'],
                   f9_patterns={name(k): v for k, v in f9['patterns'].items()},
                   graphs_found=len(d['graphs']), graph_gates_all_perfect=d['graph_gates_all_perfect_symplectic'],
                   graph_gate_patterns=dict(Counter({name(k): 0 for k in d['graph_gate_patterns']}) +
                                            Counter({name(k): v for k, v in d['graph_gate_patterns'].items()})))
        res['circulant_stabiliser'] = d.get('circulant_stabiliser')
        realised = set(res['f9_patterns']) | set(res['graph_gate_patterns'])
        res['realised'] = sorted(realised)
        res['not_realised'] = sorted(set(names.values()) - realised)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
