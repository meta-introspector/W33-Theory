#!/usr/bin/env python3
"""Pass 11183: an arrow of time that does not depend on the choice of subsystems.

Theorem 4.7 / Pass 11161: a tick followed by forgetting the partner loses the information about a qutrit's past that
the tick exported to the partner.  For an n-qutrit Clifford tick S and a tensor factorisation F (Pass 11180), the
stabilizer entropy formula gives, for each qutrit q of F,
        I(q_in : q_out) = 2 - rank S_F[others, q]        (S_F = B_F^-1 S B_F; the column block of q, other rows),
so forgetting everything but q destroys  export_q(S, F) = rank S_F[others, q]  in {0, 1, 2} trits (checked against the
entropy formula of Pass 11165).  The ARROW STRENGTH of S in the split F is E(S, F) = sum_q export_q, and the
INTRINSIC ARROW of S is
        A(S) = min over all splits F of E(S, F)      -- the information loss no choice of subsystems can avoid.
A(S) = 0 iff S is local AND keeps every subsystem in place in some split (a tick that swaps subsystems moves a whole
past into the partner, so forgetting the partner loses it); so A > 0 for all intrinsically entangling ticks of Pass 11181
and also for ticks whose only local splits permute the subsystems.
RESULTS.
  * Two qutrits (exact, all of Sp(4,3)): A is 0 (19152 ticks) or exactly 2 trits (32688 = 227/360 of all ticks) --
    never 1: symplecticity forces rank S_BA = rank S_AB, so the exported information is even.  A = 2 for every tick of
    order 5 or 9, and for the order-4 and order-6 ticks whose only local splits swap the two qutrits.
  * Three qutrits (sample of 240 plus the paper's ticks): A takes only the values 0, 2, 3 -- never 1: if two qutrits'
    columns are local and non-permuting, symplectic orthogonality makes the third plane invariant too, so information
    can never leak out of one subsystem alone; the minimal nonzero intrinsic arrow is two trits.
  * The paper's ticks: the clock has A = 0 (its 108 local splits are all arrow-free); the perfect interacting tick has
    A = 2 (attained in 162 splits); the F9 gate K = P^2 + i 1 has A = 3 -- it is local in 27 splits but cycles the
    qutrits in all of them, so every observer loses at least one trit per qutrit per tick.
Reading: an arrow of time that no choice of subsystems removes exists exactly for these ticks, and it comes in units of
at least two trits -- one full qutrit record, the unit of Theorem 4.7.
"""
from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11180_mereology as M  # noqa: E402

OUT = ROOT / "data" / "w33_pass11183_intrinsic_arrow.json"


def rank_cols(X):
    """rank over F3 of a batch of (m, r, 2) matrices"""
    zero = (X % 3 == 0).all(axis=(-1, -2))
    r = X.shape[1]
    full = np.zeros(X.shape[0], bool)
    for i in range(r):
        for j in range(i + 1, r):
            full |= (X[:, i, 0] * X[:, j, 1] - X[:, i, 1] * X[:, j, 0]) % 3 != 0
    return np.where(zero, 0, np.where(full, 2, 1))


def exports(S, Bs, n):
    """(num_splits, n) array of export_q(S, F)"""
    J = M.form(n)
    Binv = np.einsum('ij,fkj,kl->fil', -J, Bs, J) % 3
    SF = np.einsum('fij,jk,fkl->fil', Binv, S % 3, Bs) % 3
    out = np.zeros((len(Bs), n), np.int64)
    for q in range(n):
        rows = [2 * p + t for p in range(n) if p != q for t in range(2)]
        out[:, q] = rank_cols(SF[:, rows][:, :, [2 * q, 2 * q + 1]])
    return out


def check_formula(n_checks=40, seed=1):
    """export_q in the standard split = 2 - I(q_in : q_out) from the entropy formula (Pass 11165)"""
    import w33_pass11165_three_qutrit_arrow as A3
    rng = np.random.default_rng(seed)
    ok = True
    for _ in range(n_checks):
        S, _ = A3.random_symplectic(rng)
        ex = exports(S, np.eye(6, dtype=np.int64)[None], 3)[0]
        I = A3.infos(S)
        ok &= all(ex[q] == 2 - I[f"in{q}|out({q},)"] for q in range(3))
    return bool(ok)


def two_qutrit():
    import w33_pass11156_perfect_spacetime_gates as G
    Bs = np.array(M.factorisations(2))
    table = defaultdict(Counter)
    for S in G.sp43():
        E = exports(S, Bs, 2).sum(1)
        table[(M.porder(S % 3, 2), int((E == 0).sum()))][(int(E.min()), int(E.max()))] += 1
    return table


def three_qutrit(n=240, seed=11183):
    import w33_pass11165_three_qutrit_arrow as A3
    import w33_pass11182_paper_ticks_mereology as P
    Bs = P.load_facs()
    rng = np.random.default_rng(seed)
    table = defaultdict(Counter)
    for _ in range(n):
        S, _ = A3.random_symplectic(rng)
        E = exports(S, Bs, 3).sum(1)
        table[M.porder(S % 3, 3)][(int(E.min()), int(E.max()), int((E == 0).sum()))] += 1
    paper = {}
    for name, S in P.ticks().items():
        E = exports(S, Bs, 3).sum(1)
        paper[name] = dict(A=int(E.min()), max=int(E.max()), arrow_free_splits=int((E == 0).sum()),
                           splits_at_minimum=int((E == E.min()).sum()), order=M.porder(S % 3, 3))
    return table, paper


def summarize():
    t2 = two_qutrit()
    rows2 = [dict(order=k[0], local_splits=k[1], A_max={f"{a}|{b}": c for (a, b), c in sorted(v.items())})
             for k, v in sorted(t2.items())]
    A2 = Counter()
    for k, v in t2.items():
        for (a, b), c in v.items():
            A2[a] += c
    t3, paper = three_qutrit()
    rows3 = {str(o): {f"A{a}|max{b}|free{f}": c for (a, b, f), c in sorted(v.items())} for o, v in sorted(t3.items())}
    res = dict(pass_id=11183, formula_checked=check_formula(), two_qutrit=rows2,
               two_qutrit_A_distribution={str(k): v for k, v in sorted(A2.items())},
               three_qutrit_sample=rows3, paper_ticks=paper)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    r = summarize()
    print('formula', r['formula_checked'])
    for row in r['two_qutrit']:
        print(row)
    print('A distribution (two qutrits)', r['two_qutrit_A_distribution'])
    for o, v in r['three_qutrit_sample'].items():
        print(o, v)
    for k, v in r['paper_ticks'].items():
        print(k, v)
