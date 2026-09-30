#!/usr/bin/env python3
"""Pass 11207 (part 2): the arrow bound A(S) <= n beyond four qutrits -- half-moving splits at five and six qutrits.

Part 1 (w33_pass11207_arrow_by_planes; n = 2, 3 also in Pass 11188) proved, exactly for n <= 4, that the intrinsic arrow satisfies A(S) <= n with equality iff S fixes no
nondegenerate plane, and showed for every n that
    (i)  A(S) <= n - 1  =>  S fixes a nondegenerate plane           (so ticks with no invariant plane have A >= n);
    (ii) S fixes a plane P  =>  A(S) <= A(S restricted to P^perp)     (induction on n).
So for n = 5 the bound A <= 5 follows from the exact n = 4 theorem for every tick WITH an invariant plane, and for a
tick WITHOUT one it is equivalent to a HALF-MOVING SPLIT: an orthogonal decomposition into planes P_i each meeting its
image, P_i cap S P_i != 0.  Planes span(x, Sx) with omega(x, Sx) != 0 meet their image in S x, so a half-moving split
is found by a depth-first search for n mutually orthogonal such planes.
Invariant nondegenerate planes are detected exactly: such a plane either lies in a +-1 eigenspace on which omega is
nonzero, or is span(x, Sx) for a non-eigenvector x with S^2 x in span(x, Sx).
RESULTS.
  * Every conjugacy class of PSp(10,3) (GAP, if data/w33_pass11207_gap_class_reps_n5.txt is complete) and every sampled
    tick of Sp(10,3), Sp(12,3) -- including regular unipotent ticks (one Jordan block) and ticks whose characteristic
    polynomial is irreducible -- either fixes a nondegenerate plane or has a half-moving split.  The search finds one
    within a few thousand nodes: half-moving splits are abundant, not rare.
  * Hence A(S) <= 5 on all of Sp(10,3) when the GAP list is complete (with A = 5 exactly on the classes with no invariant
    plane), and A(S) <= n on every sampled tick for n = 6.
"""
from __future__ import annotations

import itertools
import json
import re
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11165_three_qutrit_arrow as A3  # noqa: E402
import w33_pass11180_mereology as M  # noqa: E402

GAPOUT = ROOT / "data" / "w33_pass11207_gap_class_reps_n5.txt"
OUT = ROOT / "data" / "w33_pass11207_half_moving_splits.json"
PSP5 = 3 ** 25 * (3 ** 2 - 1) * (3 ** 4 - 1) * (3 ** 6 - 1) * (3 ** 8 - 1) * (3 ** 10 - 1) // 2


def points(D):
    return np.array([v for v in itertools.product(range(3), repeat=D)
                     if any(v) and v[next(i for i in range(D) if v[i])] == 1], np.int64)


def invariant_plane(S, n, pts):
    """does S fix a nondegenerate plane?  (exact test, see module docstring)"""
    J = M.form(n)
    SX = (pts @ S.T) % 3
    eig_p = ((SX - pts) % 3 == 0).all(1)
    eig_m = ((SX + pts) % 3 == 0).all(1)
    for E in (pts[eig_p], pts[eig_m]):
        if len(E) > 1 and ((E @ J @ E.T) % 3 != 0).any():
            return True
    w = np.einsum('pi,ij,pj->p', pts, J, SX) % 3
    sel = (w != 0) & ~eig_p & ~eig_m
    U, V = pts[sel], SX[sel]
    S2 = (V @ S.T) % 3
    inv = np.zeros(len(U), bool)
    for a in range(3):
        for b in range(3):
            inv |= (((a * U + b * V - S2) % 3) == 0).all(1)
    return bool(inv.any())


def half_moving_split(S, n, pts, restarts=200, node_budget=50_000, seed=0):
    """n mutually orthogonal nondegenerate planes span(x, Sx): depth-first search in a seeded random order with restarts
    (a fixed point order is pathological for structured ticks such as the regular unipotent at n = 6).
    Returns (list of x's or None, nodes used)."""
    J = M.form(n)
    SX = (pts @ S.T) % 3
    w = np.einsum('pi,ij,pj->p', pts, J, SX) % 3
    eig = ((SX - pts) % 3 == 0).all(1) | ((SX + pts) % 3 == 0).all(1)
    good = np.flatnonzero((w != 0) & ~eig)
    X = np.stack([pts[good], SX[good]], 2)
    XJ = np.einsum('mai,ab->mib', X, J)
    rng = np.random.default_rng(seed)
    total = [0]

    def orth(i, cands):
        if len(cands) == 0:
            return cands
        W = np.einsum('ib,pbj->pij', XJ[i], X[cands]) % 3
        return cands[~W.reshape(len(cands), 4).any(1)]

    def dfs(cands, chosen, budget):
        if len(chosen) == n:
            return chosen
        for t, i in enumerate(cands):
            budget[0] -= 1
            total[0] += 1
            if budget[0] < 0:
                return None
            rest = cands[t + 1:]
            if len(rest) < n - len(chosen) - 1:
                return None
            r = dfs(orth(i, rest), chosen + [int(i)], budget)
            if r:
                return r
        return None
    for _ in range(restarts):
        r = dfs(rng.permutation(len(X)), [], [node_budget])
        if r:
            return [pts[good[i]] for i in r], total[0]
    return None, total[0]


def verify_split(S, n, xs):
    """the planes span(x, Sx) are nondegenerate, mutually orthogonal and each meets its image: E(S, F) = n"""
    J = M.form(n)
    P = [np.stack([x, (S @ x) % 3], 1) for x in xs]
    ok = all(int(p[:, 0] @ J @ p[:, 1]) % 3 != 0 for p in P)
    ok &= all(((p.T @ J @ q) % 3 == 0).all() for p, q in itertools.combinations(P, 2))
    B = np.concatenate([p @ np.diag([1, pow(int(p[:, 0] @ J @ p[:, 1]) % 3, -1, 3)]) % 3 for p in P], 1) % 3
    exports = A3_exports(S, B, n)
    return bool(ok), exports


def bilagrangian(S, n, xs):
    """L = span(x_i) is Lagrangian for omega AND for omega_T(x, y) = omega(x, (S + S^-1) y), and L cap SL = 0"""
    J = M.form(n)
    L = np.stack(xs, 1) % 3
    Sinv = (-J @ S.T @ J) % 3
    T = (S + Sinv) % 3
    lag = _rank3(L) == n and ((L.T @ J @ L) % 3 == 0).all()
    lag_T = ((L.T @ J @ T @ L) % 3 == 0).all()
    transverse = _rank3(np.concatenate([L, (S @ L) % 3], 1)) == 2 * n
    return bool(lag and lag_T and transverse)


def A3_exports(S, B, n):
    import w33_pass11183_intrinsic_arrow as AR
    return [int(x) for x in AR.exports(S, B[None], n)[0]]


def regular_unipotent(n):
    """[[A, A E_nn], [0, A^-T]] in the (e_1..e_n, f_1..f_n) basis, A a unipotent Jordan block, moved to (e1,f1,...)"""
    A = np.eye(n, dtype=np.int64) + np.eye(n, k=1, dtype=np.int64)
    Ainv = np.round(np.linalg.inv(A)).astype(np.int64) % 3
    E = np.zeros((n, n), np.int64)
    E[n - 1, n - 1] = 1
    U = np.block([[A, (A @ E) % 3], [np.zeros((n, n), np.int64), Ainv.T % 3]]) % 3
    perm = [k for i in range(n) for k in (i, n + i)]
    return U[np.ix_(perm, perm)] % 3


def _rank3(A):
    A = [list(map(int, r)) for r in (np.asarray(A) % 3)]
    r, rows, cols = 0, len(A), len(A[0])
    for c in range(cols):
        p = next((i for i in range(r, rows) if A[i][c] % 3), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        inv = pow(A[r][c], -1, 3)
        A[r] = [(x * inv) % 3 for x in A[r]]
        for i in range(rows):
            if i != r and A[i][c] % 3:
                f = A[i][c]
                A[i] = [(a - f * b) % 3 for a, b in zip(A[i], A[r])]
        r += 1
    return r


def charpoly_degrees(S):
    from sympy import GF, Matrix, Poly, symbols
    x = symbols('x')
    p = Poly(Matrix(S.tolist()).charpoly(x).as_expr(), x, modulus=3)
    return sorted(Poly(f, x, modulus=3).degree() for f, e in p.factor_list()[1] for _ in range(e))


def is_semisimple(A):
    """A (mod 3) is semisimple iff the product of the distinct irreducible factors of its characteristic polynomial
    annihilates it"""
    from sympy import Matrix, Poly, symbols
    x = symbols('x')
    p = Poly(Matrix((A % 3).tolist()).charpoly(x).as_expr(), x, modulus=3)
    r = Poly(1, x, modulus=3)
    for f, _e in p.factor_list()[1]:
        r = r * Poly(f, x, modulus=3)
    coeffs = [int(c) % 3 for c in r.all_coeffs()]
    D = len(A)
    R = np.zeros((D, D), np.int64)
    for c in coeffs:                                   # Horner
        R = (R @ (A % 3) + c * np.eye(D, dtype=np.int64)) % 3
    return not R.any()


def T_of(S, n):
    J = M.form(n)
    return (S + (-J @ S.T @ J)) % 3


def semisimple_census():
    """fraction of PSp(2n,3), n = 2, 3, 4, whose T = S + S^-1 is semisimple (the classes covered by the theorem)"""
    import w33_pass11207_arrow_by_planes as P88
    out = {}
    for n in (2, 3, 4):
        cls = P88.load_classes(n)
        tot = sum(c['size'] for c in cls)
        cov = [c for c in cls if is_semisimple(T_of(c['S'], n))]
        ss = [c for c in cls if is_semisimple(c['S'])]
        out[f"n{n}"] = dict(classes=len(cls), T_semisimple_classes=len(cov),
                            T_semisimple_fraction=f"{sum(c['size'] for c in cov)}/{tot}",
                            S_semisimple_fraction=f"{sum(c['size'] for c in ss)}/{tot}",
                            S_semisimple_implies_T_semisimple=all(is_semisimple(T_of(c['S'], n)) for c in ss))
    return out


def load_classes_n5():
    if not GAPOUT.exists():
        return None
    text = re.sub(r"\s+", " ", GAPOUT.read_text().replace("\\\n", ""))
    head = re.search(r"n=5 \|PSp\|=(\d+) classes=(\d+)", text)
    if head is None:
        return None
    out = []
    for order, size, mat in re.findall(r"CLASS n=5 order=(\d+) size=(\d+) fixed=-1 real=\w+ M=\s*(\[ \[.*?\] \])",
                                       text):
        rows = re.findall(r"\[([^\[\]]*)\]", mat)
        Mg = np.array([[int(x) for x in r.split(",")] for r in rows], np.int64)
        out.append(dict(order=int(order), size=int(size), S=Mg.T % 3))
    if len(out) != int(head.group(2)):
        return None
    return out


def porder(S, n, limit=20000):
    """projective order (S^k = +-1), with a limit large enough for Sp(12,3)"""
    P, I = np.eye(2 * n, dtype=np.int64), np.eye(2 * n, dtype=np.int64)
    for k in range(1, limit):
        P = (P @ S) % 3
        if np.array_equal(P, I) or np.array_equal(P, (2 * I) % 3):
            return k
    return None


def study(S, n, pts):
    S = S % 3
    row = dict(order=porder(S, n), T_semisimple=is_semisimple(T_of(S, n)), symplectic=bool(((S.T @ M.form(n) @ S - M.form(n)) % 3 == 0).all()))
    if invariant_plane(S, n, pts):
        row.update(invariant_plane=True, A_bound=n - 1)
        return row
    xs, nodes = half_moving_split(S, n, pts)
    row.update(invariant_plane=False, nodes=nodes, half_moving=xs is not None)
    if xs is not None:
        ok, ex = verify_split(S, n, xs)
        row.update(split_verified=ok and sum(ex) == n and all(e == 1 for e in ex), A=n,
                   bilagrangian=bilagrangian(S, n, xs))
    return row


def summarize(n5_random=300, n6_random=40, seed=11191):
    t0 = time.time()
    rng = np.random.default_rng(seed)
    out = dict(pass_id=11207)
    for n, N in ((5, n5_random), (6, n6_random)):
        pts = points(2 * n)
        rows = []
        U = regular_unipotent(n)
        r = study(U, n, pts)
        r.update(kind='regular unipotent', single_jordan_block=_rank3(U - np.eye(2 * n, dtype=np.int64)) == 2 * n - 1)
        rows.append(r)
        for _ in range(N):
            S, _J = A3.random_symplectic(rng, n=n, steps=40 * n)
            S = S.astype(np.int64) % 3
            r = study(S, n, pts)
            r['kind'] = 'random'
            if not r['invariant_plane'] and n == 5:
                r['charpoly_degrees'] = charpoly_degrees(S)
            rows.append(r)
        no_inv = [r for r in rows if not r['invariant_plane']]
        out[f"n{n}"] = dict(
            sampled=len(rows), with_invariant_plane=len(rows) - len(no_inv), without=len(no_inv),
            all_symplectic=all(r['symplectic'] for r in rows),
            all_without_have_half_moving_split=all(r.get('half_moving') and r.get('split_verified') for r in no_inv),
            all_bilagrangian=all(r.get('bilagrangian') for r in no_inv),
            max_nodes=max((r['nodes'] for r in no_inv), default=0),
            irreducible_charpoly=sum(1 for r in no_inv if r.get('charpoly_degrees') == [2 * n]),
            orders_without=sorted({r['order'] for r in no_inv}), rows=rows)
        OUT.write_text(json.dumps(out, indent=1, sort_keys=True, default=str))
    cls = load_classes_n5()
    if cls is not None:
        pts = points(10)
        rows = []
        for c in cls:
            r = study(c['S'], 5, pts)
            r.update(gap_order=c['order'], size=c['size'], order_ok=r['order'] == c['order'])
            rows.append(r)
        no_inv = [r for r in rows if not r['invariant_plane']]
        total = sum(r['size'] for r in rows)
        out['n5_classes'] = dict(
            classes=len(rows), total_ok=total == PSP5, all_checks=all(r['symplectic'] and r['order_ok'] for r in rows),
            without_invariant_plane=len(no_inv),
            all_without_have_half_moving_split=all(r.get('half_moving') and r.get('split_verified') for r in no_inv),
            all_bilagrangian=all(r.get('bilagrangian') for r in no_inv),
            A_equals_5_fraction_of_psp=f"{sum(r['size'] for r in no_inv)}/{total}", rows=rows)
    out['semisimple_census'] = semisimple_census()
    for n in (5, 6):
        rows = out[f"n{n}"]['rows']
        out[f"n{n}"]['T_semisimple_sampled'] = sum(1 for r in rows if r['T_semisimple'])
    out['seconds'] = round(time.time() - t0, 1)
    OUT.write_text(json.dumps(out, indent=1, sort_keys=True, default=str))
    return out


if __name__ == "__main__":
    r = summarize()
    print(r['semisimple_census'])
    for k in ('n5', 'n6', 'n5_classes'):
        if k in r:
            print(k, {kk: v for kk, v in r[k].items() if kk != 'rows'})
