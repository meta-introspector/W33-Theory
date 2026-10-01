#!/usr/bin/env python3
"""Pass 11223: the arrow of the substrate's own named ticks, and how perfect entanglement relates to the arrow.

A(S) = n - c(S) (Passes 11207, 11217) is conjugation-invariant: A = 0 exactly when S keeps n mutually orthogonal planes
invariant, i.e. when S is a product of single-qutrit gates in SOME tensor factorisation.  Pass 11193 classified
two-qutrit Cliffords relative to the STANDARD factorisation (local 1152 / perfect 13824 / partial 36864).  This pass
crosses the two.

Two qutrits (all 51840 elements of Sp(4,3)):
                 A = 0    A = 2     P(A = 0)
    local          816      336      17/24   (the A = 2 ones swap the qutrits: an exchange is an arrow)
    perfect       4896     8928      17/48   (a third of the maximally entangling gates are local in another split)
    partial      13440    23424      35/96
    total        19152    32688     133/360  (= Pass 11204's arrow-free fraction)
Named ticks: the fixed perfect gate p of Pass 11193 has A = 0 (perfect in the standard split, a product in another);
the local word h_* (order 12, with a qutrit swap) has A = 2; SUM has A = 0; every one of the 80 transvections (the
W33 point symmetries) has A = 0.
Three qutrits (Pass 11182's paper ticks): the clock K, K^2 and the dual kick V have A = 0 (products of single-qutrit
gates in 108 splits).  The perfect interacting tick VKV has A = 2, c = 1: one qutrit is dynamically protected even though
VKV is local in no split.  The F9 gate has the maximal A = 3, c = 0: in the 27 splits where Pass 11182 found it
'local' it must permute the three planes with no fixed plane, i.e. it cycles the qutrits.
Reading: 'maximally entangling' is a property of a gate AND a factorisation; the arrow is the factorisation-free
residue.  The paper's perfect two-qutrit tick carries no intrinsic arrow; exchanges (swaps, qutrit cycles) do.
"""
from __future__ import annotations

import json
import sys
from collections import Counter, deque
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11177_perfect_gates_tritangent as TG  # noqa: E402
import w33_pass11182_paper_ticks_mereology as PT  # noqa: E402
import w33_pass11193_optimal_perfect_gate_compiler as C  # noqa: E402
import w33_pass11207_arrow_theorem as AT  # noqa: E402

OUT = ROOT / "data" / "w33_pass11223_named_tick_arrows.json"
J = TG.J % 3
H_STAR = np.array([[0, 0, 0, 1], [0, 0, 2, 0], [0, 1, 0, 0], [2, 1, 0, 0]], np.int64)
SUM = np.array([[1, 0, 0, 0], [0, 1, 0, 2], [1, 0, 1, 0], [0, 0, 0, 1]], np.int64).T


def order(S, cap=100):
    P = S.copy()
    for k in range(1, cap):
        if np.array_equal(P % 3, np.eye(len(S), dtype=np.int64)):
            return k
        P = (P @ S) % 3
    return None


def sp43():
    vecs = [np.array(v) for v in np.ndindex(3, 3, 3, 3) if any(v)]
    gens = [(np.eye(4, dtype=np.int64) + np.outer(v, v @ J)) % 3 for v in vecs]
    gens = [g for g in gens if np.array_equal((g.T @ J @ g) % 3, J)]
    seen = {C.key(np.eye(4, dtype=np.int64)): np.eye(4, dtype=np.int64)}
    dq = deque([np.eye(4, dtype=np.int64)])
    while dq:
        g = dq.popleft()
        for t in gens:
            h = (g @ t) % 3
            k = C.key(h)
            if k not in seen:
                seen[k] = h
                dq.append(h)
    return list(seen.values()), gens


def run():
    G, transvections = sp43()
    tab = Counter((C.relation(S), AT.arrow(S, 2)) for S in G)
    rows = {}
    for rel in ("local", "perfect", "partial"):
        tot = sum(v for (r, a), v in tab.items() if r == rel)
        z = tab.get((rel, 0), 0)
        rows[rel] = dict(total=tot, A0=z, A2=tot - z, p_arrow_free=str(Fraction(z, tot)))
    # local gates with A = 2 are exactly the qutrit-swapping ones?
    swap_A2 = Counter()
    for S in G:
        if C.relation(S) == "local":
            swapped = bool(S[0:2, 2:4].any())
            swap_A2[(swapped, AT.arrow(S, 2))] += 1
    p = C.fixed_perfect_gate()
    named2 = {name: dict(relation=C.relation(S), A=AT.arrow(S, 2), order=order(S))
              for name, S in (("p", p), ("p^2", (p @ p) % 3), ("h_star", H_STAR), ("SUM", SUM),
                              ("p h_star p", (p @ H_STAR @ p) % 3))}
    trans_A = Counter(AT.arrow(T, 2) for T in transvections)
    named3 = {name: dict(A=AT.arrow(S % 3, 3), c=3 - AT.arrow(S % 3, 3), order=order(S % 3))
              for name, S in PT.ticks().items()}
    res = dict(pass_id=11223, group_order=len(G), by_relation=rows,
               local_swap_vs_arrow={f"swap={k[0]},A={k[1]}": v for k, v in sorted(swap_A2.items())},
               arrow_free_total=str(Fraction(sum(r["A0"] for r in rows.values()), len(G))),
               named_two_qutrit=named2, transvections=dict(count=len(transvections), A=dict(trans_A)),
               named_three_qutrit=named3)
    return res


def main():
    res = run()
    OUT.write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
