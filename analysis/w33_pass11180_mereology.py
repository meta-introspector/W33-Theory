#!/usr/bin/env python3
"""Pass 11180: quantum mereology on W(3,3) -- the dynamics decides the subsystems, and some ticks entangle in EVERY
subsystem decomposition.

A tensor factorisation of n qutrits (stabilizer-compatible) is a decomposition F_3^{2n} = P_1 + ... + P_n into mutually
orthogonal nondegenerate planes (45 for n = 2 -- the tritangent planes, Pass 11177; 110565 for n = 3, Pass 11178).  For a
Clifford tick S and a factorisation F with adapted symplectic basis B_F, S_F = B_F^{-1} S B_F is 'S as seen by the
subsystems of F':
    S is LOCAL in F      <=>  S_F is block-monomial (S maps each plane of F onto a plane of F);
    S is PERFECT in F    <=>  every block of S_F is invertible (maximal scrambling relative to F);
    otherwise S is partially entangling in F.
Entanglement generation is relative to the factorisation (Zanardi 2001; Zanardi-Lidar-Lloyd 2004; Carroll-Singh 2021 --
'quantum mereology').  Here the whole census is exact.
TWO QUTRITS (all 51840 elements of Sp(4,3); cross-checked against the octet computation of Pass 11177):
the profile (#local, #perfect, #partial) over the 45 factorisations depends only on the projective order and class:
    order 1:  (45, 0, 0)                 order 2: (13, 0, 32) x90, (5, 24, 16) x540
    order 3:  (9, 0, 36) x160, (3, 0, 42) x960, (6, 27, 12) x480
    order 4:  (1, 12, 32) x7560          order 12: (1, 12, 32) x8640
    order 5:  (0, 15, 30) x10368         order 9:  (0, 9, 36) x11520
    order 6:  (1, 0, 44) x1440, (1, 18, 26) x2880, (2, 15, 28) x4320, (4, 9, 32) x2880
  * INTRINSICALLY ENTANGLING: a tick is local in NO factorisation exactly when its projective order is 5 or 9
    (21888 of 51840, 42.2%) -- no choice of subsystems makes it non-entangling;
  * UNIQUE EMERGENT SUBSYSTEMS: 20520 ticks (39.6%; orders 4, 6, 12) are local in exactly one factorisation, which the
    dynamics therefore singles out;
  * RELATIVITY OF PERFECTION: the 480 order-3 ticks of profile (6, 27, 12) are product gates for 6 subsystem splits and
    maximal scramblers for 27.
The class names (U4(2) = PSp(4,3), ATLAS) are attached by GAP (analysis/gap/w33_pass11180_mereology_classes.g).
"""
from __future__ import annotations

import itertools
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
OUT = ROOT / "data" / "w33_pass11180_mereology.json"
GAPOUT = ROOT / "data" / "w33_pass11180_gap_classes.txt"


def form(n):
    J = np.zeros((2 * n, 2 * n), np.int64)
    for k in range(n):
        J[2 * k, 2 * k + 1], J[2 * k + 1, 2 * k] = 1, -1
    return J


def factorisations(n):
    """all tensor factorisations of n qutrits, each as an adapted symplectic basis matrix (columns u1,v1,...,un,vn)"""
    J = form(n)
    D = 2 * n
    vecs = [np.array(v) for v in itertools.product(range(3), repeat=D) if any(v)]
    pts = {}
    for v in vecs:
        i = next(k for k in range(D) if v[k])
        key = tuple((v * pow(int(v[i]), -1, 3)) % 3)
        pts[key] = np.array(key)
    P = np.array(list(pts.values()))                                   # projective points
    W = (P @ J @ P.T) % 3
    planes = {}
    for a in range(len(P)):
        for b in range(a + 1, len(P)):
            if W[a, b]:
                u, v = P[a], (P[b] * pow(int(W[a, b]), -1, 3)) % 3     # omega(u, v) = 1
                span = frozenset(tuple(((x * u + y * v) % 3)) for x in range(3) for y in range(3))
                if span not in planes:
                    planes[span] = (u, v)
    plist = list(planes.items())
    if n == 2:
        facs = {}
        for span, (u, v) in plist:
            rest = [(s2, b2) for s2, b2 in plist if all((b2[0] @ J @ w) % 3 == 0 and (b2[1] @ J @ w) % 3 == 0 for w in (u, v))]
            for s2, (u2, v2) in rest:
                key = frozenset([span, s2])
                if key not in facs:
                    facs[key] = np.stack([u, v, u2, v2], 1)
        return list(facs.values())
    # n = 3: P, then an orthogonal factorisation of the 4-dim complement
    facs = {}
    orth = defaultdict(list)
    for i, (s1, (u1, v1)) in enumerate(plist):
        for j, (s2, (u2, v2)) in enumerate(plist):
            if j > i and all((a @ J @ b) % 3 == 0 for a in (u1, v1) for b in (u2, v2)):
                orth[i].append(j)
    for i in range(len(plist)):
        oi = set(orth[i])
        for j in orth[i]:
            for k in orth[j]:
                if k in oi:
                    key = frozenset([plist[i][0], plist[j][0], plist[k][0]])
                    if key not in facs:
                        cols = [plist[i][1][0], plist[i][1][1], plist[j][1][0], plist[j][1][1], plist[k][1][0], plist[k][1][1]]
                        facs[key] = np.stack(cols, 1)
    return list(facs.values())


def profile(S, Bs, n, signatures=False):
    """counts of factorisations in which S is local / perfect / partial (and optionally the block-rank signatures)"""
    J = form(n)
    Binv = np.einsum('ij,fkj,kl->fil', -J, Bs, J) % 3                    # -J B^T J
    SF = np.einsum('fij,jk,fkl->fil', Binv, S % 3, Bs) % 3
    Bk = SF.reshape(-1, n, 2, n, 2).transpose(0, 1, 3, 2, 4)
    det = (Bk[..., 0, 0] * Bk[..., 1, 1] - Bk[..., 0, 1] * Bk[..., 1, 0]) % 3
    zero = (Bk == 0).all(axis=(-1, -2))
    nonzero_per_row = (~zero).sum(-1)
    local = (nonzero_per_row == 1).all(-1)
    perfect = (det != 0).all(axis=(-1, -2))
    out = dict(local=int(local.sum()), perfect=int(perfect.sum()), partial=int((~local & ~perfect).sum()))
    if signatures:
        code = np.where(zero, 0, np.where(det == 0, 1, np.where(det == 1, 2, 3)))
        out['codes'] = code
    return out


def porder(S, n):
    P = np.eye(2 * n, dtype=np.int64)
    I = np.eye(2 * n, dtype=np.int64)
    for k in range(1, 400):
        P = (P @ S) % 3
        if np.array_equal(P, I) or np.array_equal(P, (2 * I) % 3):
            return k
    return None


def two_qutrit_census():
    import w33_pass11156_perfect_spacetime_gates as G
    Bs = np.array(factorisations(2))
    table = defaultdict(Counter)
    for S in G.sp43():
        p = profile(S, Bs, 2)
        table[(p['local'], p['perfect'], p['partial'])][porder(S % 3, 2)] += 1
    return len(Bs), table


def summarize():
    nF, table = two_qutrit_census()
    rows = [dict(local=k[0], perfect=k[1], partial=k[2], orders={str(o): c for o, c in sorted(v.items())},
                 count=sum(v.values())) for k, v in sorted(table.items())]
    total = sum(r['count'] for r in rows)
    intrinsic = [r for r in rows if r['local'] == 0]
    unique = [r for r in rows if r['local'] == 1]
    res = dict(pass_id=11180, factorisations=nF, rows=rows, total=total,
               intrinsically_entangling=sum(r['count'] for r in intrinsic),
               intrinsically_entangling_orders=sorted({o for r in intrinsic for o in r['orders']}),
               unique_subsystems=sum(r['count'] for r in unique),
               max_perfect=max(r['perfect'] for r in rows),
               gap_classes=GAPOUT.read_text() if GAPOUT.exists() else None)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    r = summarize()
    for row in r['rows']:
        print(row)
    print({k: v for k, v in r.items() if k not in ('rows', 'gap_classes')})
