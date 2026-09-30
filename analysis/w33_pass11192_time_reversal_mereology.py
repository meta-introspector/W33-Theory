#!/usr/bin/env python3
"""Pass 11192: time reversal on the substrate -- the anti-symplectic half of W(E6) and the subsystems it respects.

A time reversal of n qutrits (an antiunitary Clifford operation, mod Paulis and phases) acts on the phase space F_3^{2n}
by an ANTI-SYMPLECTIC map theta: omega(theta x, theta y) = -omega(x, y); complex conjugation is tau = (x, z) -> (x, -z).
Projectively they form the outer coset of PGSp(2n,3) = PSp(2n,3).2; for two qutrits PGSp(4,3) = W(E6), whose outer
(det -1) coset contains the 36 reflections.  GAP (analysis/gap/w33_pass11192_antisymplectic_classes.g) gives the outer
classes (10 for n = 2, 26 for n = 3).  A time reversal is LOCAL in a split F if it maps every plane of F onto a plane of
F (a product of single-qutrit antiunitaries, possibly composed with a relabelling of the qutrits).
THEOREMS (every n).
  (1) theta^2 = +1  =>  theta is a LOCAL COMPLEX CONJUGATION in some split: its eigenspaces E+, E- are transverse
      Lagrangians (omega(x, y) = -omega(x, y) on each), and dual bases e_i in E+, f_i in E- with omega(e_i, f_j) = d_ij
      give theta-invariant planes <e_i, f_i> on which theta = diag(1, -1) = tau.
  (2) theta^2 = -1 (a 'Kramers' time reversal: T^2 = the parity operator)  =>  n is EVEN: h(x, y) = omega(theta x, y)
      + i omega(x, y) (i acting as theta) is an F9-bilinear alternating nondegenerate form on the F9-space V.  Moreover no
      single qutrit admits an anti-symplectic map with square -1 (a 2x2 map of determinant -1 squares to a scalar only
      when its trace vanishes, and then to +1), so such a theta is local only in splits where it PAIRS the qutrits
      (a fixed-point-free involution) -- for odd n it is local in no split at all.
RESULTS.
  * n = 2 (exact, all 51840 anti-symplectic matrices and the 10 outer classes of W(E6)):
      theta^2 = +1 : 1080 matrices = 540 projective = the W(E6) class of products of three orthogonal reflections;
                     each is local in 7 of the 45 splits -- 6 times as local complex conjugation, once with a swap.
      theta^2 = -1 : 72 matrices = 36 projective = THE 36 REFLECTIONS OF E6; each is local in exactly 15 splits, ALL of
                     them swapping the two qutrits.  The E6 reflections are exactly the Kramers-type time reversals.
      local in no split: 21888 matrices (42.2%) -- intrinsically entangling time reversals.
  * n = 3 (26 outer classes of PGSp(6,3)): no class with theta^2 = -1 (as theorem (2) requires); every involution class
    is a local complex conjugation somewhere; the fraction of time reversals local in no split is frozen in the JSON.
"""
from __future__ import annotations

import json
import re
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11180_mereology as M  # noqa: E402

GAPOUT = ROOT / "data" / "w33_pass11192_gap_outer_classes.txt"
OUT = ROOT / "data" / "w33_pass11192_time_reversal_mereology.json"
PSP = {2: 25920, 3: 4585351680}


def load_outer(n):
    text = re.sub(r"\s+", " ", GAPOUT.read_text().replace("\\\n", ""))
    out = []
    for nn, order, size, mat in re.findall(r"OUTER n=(\d) order=(\d+) size=(\d+) M=\s*(\[ \[.*?\] \])", text):
        if int(nn) == n:
            rows = re.findall(r"\[([^\[\]]*)\]", mat)
            Mg = np.array([[int(x) for x in r.split(",")] for r in rows], np.int64)
            out.append(dict(order=int(order), size=int(size), T=Mg.T % 3))
    return out


def anti_symplectic(T, n):
    J = M.form(n)
    return bool(((T.T @ J @ T + J) % 3 == 0).all())


def square_sign(T, n):
    I = np.eye(2 * n, dtype=np.int64)
    T2 = (T @ T) % 3
    return 1 if np.array_equal(T2, I) else (-1 if np.array_equal(T2, (2 * I) % 3) else 0)


def locality(T, Bs, n):
    """(#splits where T is local, #of those where it fixes every plane, Counter of the induced qutrit permutations'
    cycle types)"""
    J = M.form(n)
    Binv = np.einsum('ij,fkj,kl->fil', -J, Bs, J) % 3
    TF = np.einsum('fij,jk,fkl->fil', Binv, T % 3, Bs) % 3
    Bk = TF.reshape(-1, n, 2, n, 2).transpose(0, 1, 3, 2, 4)
    nz = ~(Bk == 0).all((-1, -2))
    loc = (nz.sum(-1) == 1).all(-1)
    idx = np.flatnonzero(loc)
    cyc = Counter()
    fixed = 0
    for f in idx:
        perm = [int(np.argmax(nz[f, i])) for i in range(n)]
        seen, lens = set(), []
        for s in range(n):
            if s not in seen:
                L, x = 0, s
                while x not in seen:
                    seen.add(x)
                    x = perm[x]
                    L += 1
                lens.append(L)
        cyc[tuple(sorted(lens))] += 1
        fixed += all(perm[i] == i for i in range(n))
    return int(loc.sum()), int(fixed), {'+'.join(map(str, k)): v for k, v in sorted(cyc.items())}


def is_local_conjugation(T, B, n):
    """in the split with adapted basis B, T fixes every plane and acts on it as an anti-symplectic involution"""
    J = M.form(n)
    TF = ((-J @ B.T @ J) @ T @ B) % 3
    for q in range(n):
        blk = TF[2 * q:2 * q + 2, 2 * q:2 * q + 2]
        if not (np.delete(TF[:, 2 * q:2 * q + 2], [2 * q, 2 * q + 1], 0) == 0).all():
            return False
        if not np.array_equal((blk @ blk) % 3, np.eye(2, dtype=np.int64)):
            return False
    return True


def one_qutrit_lemma():
    """every 2x2 anti-symplectic (det = -1) map squares to +1 or to a non-scalar; never to -1"""
    import itertools
    sq = Counter()
    for a, b, c, d in itertools.product(range(3), repeat=4):
        if (a * d - b * c) % 3 == 2:
            T = np.array([[a, b], [c, d]])
            sq[square_sign(T, 1)] += 1
    return {str(k): v for k, v in sq.items()}


def two_qutrit_exhaustive():
    import w33_pass11156_perfect_spacetime_gates as G
    Sp = np.array(list(G.sp43())) % 3
    tau = np.diag([1, 2, 1, 2]).astype(np.int64)
    Th = np.einsum('ij,njk->nik', tau, Sp) % 3
    Bs = np.array(M.factorisations(2))
    tab = Counter()
    for T in Th:
        l, f, cyc = locality(T, Bs, 2)
        tab[(square_sign(T, 2), l, f)] += 1
    return {f"sq{s}|local{l}|fixing{f}": c for (s, l, f), c in sorted(tab.items())}, len(Th)


def summarize():
    import w33_pass11182_paper_ticks_mereology as P
    out = dict(pass_id=11192, one_qutrit=one_qutrit_lemma())
    ex, total = two_qutrit_exhaustive()
    out['n2_exhaustive'] = dict(total=total, table=ex)
    for n, Bs in ((2, np.array(M.factorisations(2))), (3, P.load_facs())):
        rows = []
        for c in load_outer(n):
            T = c['T']
            l, f, cyc = locality(T, Bs, n)
            conj = any(is_local_conjugation(T, Bs[i], n) for i in range(len(Bs))) if square_sign(T, n) == 1 else None
            rows.append(dict(order=c['order'], size=c['size'], anti_symplectic=anti_symplectic(T, n),
                              square=square_sign(T, n), local_splits=l, plane_fixing_splits=f, cycle_types=cyc,
                              local_conjugation_somewhere=conj))
        total = sum(r['size'] for r in rows)
        nowhere = sum(r['size'] for r in rows if r['local_splits'] == 0)
        out[f"n{n}"] = dict(
            outer_classes=len(rows), total=total, total_ok=total == PSP[n],
            all_anti_symplectic=all(r['anti_symplectic'] for r in rows),
            involution_classes=[dict(size=r['size'], square=r['square'], local=r['local_splits'],
                                     fixing=r['plane_fixing_splits'], cycles=r['cycle_types'])
                                for r in rows if r['order'] == 2],
            kramers_classes=sum(1 for r in rows if r['square'] == -1),
            every_involution_is_local_conjugation=all(r['local_conjugation_somewhere'] for r in rows
                                                      if r['square'] == 1),
            kramers_only_pairing=all(r['plane_fixing_splits'] == 0 and all('1' not in k.split('+')
                                                                            for k in r['cycle_types'])
                                     for r in rows if r['square'] == -1),
            local_nowhere=nowhere, local_nowhere_fraction=str(Fraction(nowhere, total)),
            local_nowhere_orders=sorted({r['order'] for r in rows if r['local_splits'] == 0}),
            rows=rows)
    OUT.write_text(json.dumps(out, indent=1, sort_keys=True))
    return out


if __name__ == "__main__":
    r = summarize()
    print('one qutrit squares:', r['one_qutrit'])
    print('two qutrits exhaustive:', r['n2_exhaustive'])
    for n in (2, 3):
        d = r[f"n{n}"]
        print(f"n={n}", {k: v for k, v in d.items() if k != 'rows'})
        for row in d['rows']:
            print('   ', {k: row[k] for k in ('order', 'size', 'square', 'local_splits', 'plane_fixing_splits',
                                              'cycle_types')})
