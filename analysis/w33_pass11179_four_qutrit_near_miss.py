#!/usr/bin/env python3
"""Pass 11179: how close four qutrits come to a perfect tick.

AME(8,3) does not exist (shadow inequalities; Huber, Eltschka, Siewert, Guhne, arXiv:1708.06298), so every 8-qutrit
stabilizer state -- up to local Cliffords, every weighted qutrit graph on 8 vertices -- has a balanced cut S|S^c
(|S| = 4) with Gamma[S, S^c] singular over F3.  Tabu search (analysis/w33_pass11179_scan_four_qutrit_nearmiss.py):
  * minimising the NUMBER of deficient balanced cuts: 120 of 120 restarts reach exactly 2, never 1 (0 is excluded);
  * minimising the TOTAL deficit sum (4 - rank): 36 of 36 restarts reach exactly 4.
Every optimum found has the same rigid structure (checked here on each recorded optimum):
  * the two deficient cuts each have rank 2 -- two trits short of the maximum four, not one;
  * the state is 3-uniform: every 3|5 cut is maximal (all 56);
  * the 8 qutrits split into four pairs P0, P1, P2, P3 with the deficient cuts P0+P1 | P2+P3 and P0+P2 | P1+P3, while
    the third pairing P0+P3 | P1+P2 is maximal.
As a tick.  Taking the maximal cut P0+P3 | P1+P2 as input | output gives a four-qutrit Clifford gate (symplectic,
checked) whose Choi state is maximal on every cut except the two above: every single-qutrit block is invertible, and
the only defect is that the four ALIGNED pair-to-pair blocks (input pair P0 or P3 to output pair P1 or P2) have rank 2
instead of 4 -- each carries two trits, not four (Jacobi pairs them: det S[I,J] = +-det S[I^c,J^c]).
Scope: the minimum 2 (and total deficit 4) is what the search reaches every time -- a numerical statement, not a proof
that one deficient cut is impossible (the shadow bound excludes only 0).
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11179_scan_four_qutrit_nearmiss as N  # noqa: E402

COUNT = ROOT / "data" / "w33_pass11179_nearmiss_scan.json"
DEFICIT = ROOT / "data" / "w33_pass11179_nearmiss_deficit_scan.json"
OUT = ROOT / "data" / "w33_pass11179_four_qutrit_near_miss.json"


def inv_mod3(A):
    n = A.shape[0]
    M = np.concatenate([A % 3, np.eye(n, dtype=np.int64)], 1)
    r = 0
    for c in range(n):
        p = next((i for i in range(r, n) if M[i, c] % 3), None)
        if p is None:
            return None
        M[[r, p]] = M[[p, r]]
        M[r] = (M[r] * pow(int(M[r, c]), -1, 3)) % 3
        for i in range(n):
            if i != r and M[i, c]:
                M[i] = (M[i] - M[i, c] * M[r]) % 3
        r += 1
    return M[:, n:]


def gate(G, inputs):
    nq = G.shape[0]
    outputs = [v for v in range(nq) if v not in inputs]
    L = np.zeros((nq, 2 * nq), np.int64)
    for v in range(nq):
        L[v, 2 * v] = 1
        L[v, 1::2] = G[v] % 3
    ci = [2 * q + t for q in inputs for t in range(2)]
    co = [2 * q + t for q in outputs for t in range(2)]
    Ai = inv_mod3(L[:, ci].T)
    M = (L[:, co].T @ Ai) % 3
    T = np.diag([1, -1] * len(inputs)) % 3
    return (M @ T) % 3, outputs


def structure(w):
    G = N.to_mats(np.array(w)[None])[0]
    cuts, u3 = N.profile(w)
    (S1, r1), (S2, r2) = cuts
    P0 = sorted(set(S1) & set(S2))
    P1, P2 = sorted(set(S1) - set(P0)), sorted(set(S2) - set(P0))
    P3 = sorted(set(range(8)) - set(S1) - set(S2))
    third = P0 + P3
    Sc = [v for v in range(8) if v not in third]
    third_rank = N.rank3(G[np.ix_(third, Sc)])
    # the tick: inputs P0 + P3, outputs P1 + P2
    S, outs = gate(G, P0 + P3)
    n = 4
    J = np.zeros((8, 8), np.int64)
    for k in range(n):
        J[2 * k, 2 * k + 1], J[2 * k + 1, 2 * k] = 1, -1
    symp = bool(np.array_equal((S.T @ J @ S) % 3, J % 3))
    B = S.reshape(n, 2, n, 2).transpose(0, 2, 1, 3)
    det1 = (B[..., 0, 0] * B[..., 1, 1] - B[..., 0, 1] * B[..., 1, 0]) % 3
    ins = P0 + P3
    in_pos = {q: k for k, q in enumerate(ins)}
    out_pos = {q: k for k, q in enumerate(outs)}

    def two_party_rank(out_pair, in_pair):
        rows = [2 * out_pos[q] + t for q in out_pair for t in range(2)]
        cols = [2 * in_pos[q] + t for q in in_pair for t in range(2)]
        return N.rank3(S[np.ix_(rows, cols)])
    ranks = {f"out{op}|in{ip}": two_party_rank(op, ip) for op in (P1, P2) for ip in (P0, P3)}
    return dict(pairs=[P0, P1, P2, P3], deficient_ranks=[r1, r2], third_pairing_rank=third_rank, three_uniform=u3,
                tick_symplectic=symp, single_blocks_invertible=bool((det1 != 0).all()), two_party_ranks=ranks)


def summarize():
    c = json.loads(COUNT.read_text())
    d = json.loads(DEFICIT.read_text())
    structs = [structure(b['weights']) for b in c['best'][:10]] + [structure(b['weights']) for b in d['best'][:10]]
    res = dict(pass_id=11179, count_distribution=c['distribution'], count_minimum=c['minimum'],
               deficit_distribution=d['distribution'], deficit_minimum=d['minimum'],
               checked=len(structs),
               all_rank_two=all(s['deficient_ranks'] == [2, 2] for s in structs),
               all_three_uniform=all(s['three_uniform'] for s in structs),
               all_third_pairing_maximal=all(s['third_pairing_rank'] == 4 for s in structs),
               all_ticks_symplectic=all(s['tick_symplectic'] for s in structs),
               all_single_blocks_invertible=all(s['single_blocks_invertible'] for s in structs),
               example=structs[0])
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
