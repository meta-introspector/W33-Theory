"""Pass 11192: time reversal on the two-qutrit substrate -- the anti-symplectic half of W(E6) and the 45 splits.

The automorphism group of W(3,3) is PSp(4,3).2 = W(E6); its outer coset consists of the anti-symplectic maps
T (T^t J T = -J), realised on Hilbert space by anti-unitaries (Appleby 2005; the multi-qutrit case via complex
conjugation (x,z) -> (x,-z) = tau): Wigner's time reversals.  Pass 11180 measured, for unitary ticks, in which of the 45
tensor factorisations a tick is local; Pass 11183 defined the intrinsic arrow A.  Here the same two measurements are
made on all 51840 anti-symplectic matrices tau*S, S in Sp(4,3):

  * local(T)  = number of splits in which T is block-monomial (it maps every factor to a factor),
  * A(T)      = min over the 45 splits of the total export sum_q rank T_F[others, q].

and the 36 reflections of W(E6) are located among them (anti-symplectic involutions with class size 36 in
PSp(4,3).2), with their split profile.
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11156_perfect_spacetime_gates as G  # noqa: E402
import w33_pass11180_mereology as M  # noqa: E402
import w33_pass11183_intrinsic_arrow as AR  # noqa: E402

OUT = ROOT / "data" / "w33_pass11192_time_reversal_mereology.json"
TAU = np.diag([1, 2, 1, 2]).astype(np.int64)


def key(S):
    """projective key: S and -S identified"""
    a, b = tuple((S % 3).ravel()), tuple(((-S) % 3).ravel())
    return min(a, b)


def local_counts(T, Bs):
    J = M.form(2)
    Binv = np.einsum('ij,fkj,kl->fil', -J, Bs, J) % 3
    TF = np.einsum('fij,jk,fkl->fil', Binv, T % 3, Bs) % 3
    Bk = TF.reshape(-1, 2, 2, 2, 2).transpose(0, 1, 3, 2, 4)
    nz = ~(Bk == 0).all(axis=(-1, -2))
    return int(((nz.sum(-1) == 1).all(-1)).sum())


def run():
    Bs = np.array(M.factorisations(2))
    J = M.form(2)
    sp = [S % 3 for S in G.sp43()]
    assert len(sp) == 51840
    anti = [(TAU @ S) % 3 for S in sp]
    for T in anti[:50]:
        assert np.array_equal((T.T @ J @ T) % 3, (-J) % 3)
    # projective elements of the outer coset, and their conjugacy classes inside PSp(4,3).2 (conjugating by Sp only
    # is enough for sizes up to a factor; we compute the orbit under conjugation by the whole extended group)
    reps = {}
    for T in anti:
        reps.setdefault(key(T), T)
    assert len(reps) == 25920
    table = Counter()
    rows = {}
    for k, T in reps.items():
        E = AR.exports(T, Bs, 2).sum(1)
        loc = local_counts(T, Bs)
        o = M.porder(T, 2)
        rows[k] = (o, loc, int(E.min()))
        table[(o, loc, int(E.min()))] += 1
    # conjugacy classes of involutions in the outer coset (projective order 2)
    ext = sp + anti
    inv = [k for k, v in rows.items() if v[0] == 2]
    seen, classes = set(), []
    for k in inv:
        if k in seen:
            continue
        T = np.array(k).reshape(4, 4)
        cls = set()
        for g in ext:
            # g^-1 over F3: g symplectic -> -J g^T J ; g anti-symplectic -> J g^T J
            s = 1 if np.array_equal((g.T @ J @ g) % 3, J % 3) else -1
            gi = (-s * (J @ g.T @ J)) % 3
            cls.add(key((g @ T @ gi) % 3))
        seen |= cls
        classes.append(dict(size=len(cls), order=2, local=rows[k][1], A=rows[k][2]))
    arrow = Counter()
    for (o, loc, a), n in table.items():
        arrow[a] += n
    res =dict(pass_id=11192, outer_coset_projective_elements=25920,
               arrow_distribution={int(a): n for a, n in sorted(arrow.items())},
               never_local=int(sum(n for (o, loc, a), n in table.items() if loc == 0)),
               profile={f"order{o}|local{loc}|A{a}": n for (o, loc, a), n in sorted(table.items())},
               involution_classes=classes)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1)
    for k, v in res.items():
        print(k, ":", v)


if __name__ == "__main__":
    main()
