#!/usr/bin/env python3
"""Pass 11224: AME(10,3) stabilizer states are F9-linear up to local Cliffords -- the Hermitian GRS code, not Glynn's.

An AME(10,3) stabilizer state is a Lagrangian C of F_3^20 (10 qutrits) in which every 5 qutrits carry no stabilizer
element: a self-dual additive [[10,0,6]]_3 code.  It is F9-LINEAR up to local Cliffords if each qutrit's F_3^2 can be
given an F9-structure -- an element M_v of order 8 in GL(2,3) acting as multiplication by a primitive xi in F9 -- with
(+) M_v mapping C to C.  Local Cliffords conjugate the M_v, so allowing all 12 order-8 elements per qutrit tests
F9-linearity up to LC.  Because every 5 qutrits form an information set, C = {(x, E x)} over S = {0..4}, and C is
(+) M_v-invariant iff E M_S E^-1 is block diagonal with order-8 blocks: a batch test over 12^5 choices of M_S.

An F9-linear AME(10,3) code is a [10,5,6]_9 MDS code.  Length q+1 = 10 over F9 admits two kinds: generalised
Reed-Solomon (from the normal rational curve) and Glynn's non-classical 10-arc in PG(4,9) (Glynn 1986).  They are told
apart by the Schur square: dim(C * C) = 2k - 1 = 9 for GRS, larger otherwise.

Results.
  * All 71 AME(10,3) graph states of master's Pass 11186 census (data/w33_pass11186_ame10_tabu.json) are F9-linear up
    to LC, each with exactly 4 F9 structures (+-xi, Frobenius) of uniform trace; as [10,5,6]_9 codes they have Schur
    square of dimension 10, so none is GRS: by Glynn's classification of 10-arcs in PG(4,9) they are all GLYNN'S code.
  * No GRS [10,5]_9 code is Hermitian self-dual for any coordinate weights (0 of 1024), so there is no Reed-Solomon
    AME(10,3) state; Glynn's code is Hermitian self-dual for exactly the uniform weights (2 of 1024, = +-1).
  * Local-Clifford automorphisms of the Glynn state: the arc's stabiliser in PGL(5,9) (360) and its Frobenius-twisted
    part (360) all have coordinate scalings of constant norm, giving |Aut_LC| = 4 (360 + 360) = 2880, acting on the 10
    qutrits as PGL(2,9) (sharply 3-transitive; linear part PSL(2,9) = A6) -- the symmetry master's Pass 11190 found for
    the sign structure.
  * Orbit count: 24^10 10! / 2880 = 79 888 260 016 373 760 = 6^10 x 2 642 411 520 / 2, master's EXHAUSTIVE count of
    F9-linear perfect five-qutrit gates (Pass 11175) times the 6^10 local F9 structures over the 2 structures per code.
    So the F9-linear AME(10,3) states form ONE orbit: the AME(10,3) stabilizer state of F9 type is unique up to local
    Cliffords and qutrit relabelling, and there are exactly 79 888 260 016 373 760 perfect five-qutrit Clifford gates of
    that type.
  * Open: whether non-F9-linear AME(10,3) stabilizer states exist (none in 300 tabu restarts).
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11175_scan_five_qutrit as F  # noqa: E402
import w33_pass11190_ame10_sign_structure as P  # noqa: E402

OUT = ROOT / "data" / "w33_pass11224_ame10_f9_structure.json"
GL23 = [np.array([[a, b], [c, d]]) for a in range(3) for b in range(3) for c in range(3) for d in range(3)
        if (a * d - b * c) % 3]


def mpow(M, k):
    R = np.eye(len(M), dtype=np.int64)
    for _ in range(k):
        R = (R @ M) % 3
    return R


ORD8 = [M for M in GL23 if np.array_equal(mpow(M, 4), (2 * np.eye(2, dtype=np.int64)) % 3)]   # M^4 = -1


def inv3(A):
    A = A.copy() % 3
    n = len(A)
    M = np.hstack([A, np.eye(n, dtype=np.int64)])
    r = 0
    for c in range(n):
        p = next(i for i in range(r, n) if M[i, c])
        M[[r, p]] = M[[p, r]]
        M[r] = (M[r] * (1 if M[r, c] == 1 else 2)) % 3
        for i in range(n):
            if i != r and M[i, c]:
                M[i] = (M[i] - M[i, c] * M[r]) % 3
        r += 1
    return M[:, n:]


def encoder(G):
    """E with C = {(x_S, E x_S)}, S = qutrits 0..4 (coordinates 0..9), S^c = qutrits 5..9 (coordinates 10..19)"""
    L = P.stabiliser(G) % 3                     # rows generate C, columns (x_v, z_v) per qutrit
    A = L[:, :10]
    B = L[:, 10:]
    Ai = inv3(A)
    return (B.T @ Ai.T) % 3                     # c_{S^c} = E c_S  for c = row combination


def f9_structures(G):
    E = encoder(G)
    Ei = inv3(E)
    combos = np.array(list(itertools.product(range(len(ORD8)), repeat=5)))
    Ms = np.array(ORD8)
    out = []
    for chunk in np.array_split(combos, 16):
        MS = np.zeros((len(chunk), 10, 10), np.int64)
        for q in range(5):
            MS[:, 2 * q:2 * q + 2, 2 * q:2 * q + 2] = Ms[chunk[:, q]]
        T = np.einsum('ij,njk,kl->nil', E, MS, Ei) % 3
        mask = np.ones(10, bool)
        off = np.ones((10, 10), bool)
        for q in range(5):
            off[2 * q:2 * q + 2, 2 * q:2 * q + 2] = False
        ok = ~(T[:, off].any(1))
        for idx in np.flatnonzero(ok):
            blocks = [T[idx, 2 * q:2 * q + 2, 2 * q:2 * q + 2] for q in range(5)]
            if all(any(np.array_equal(b, m) for m in ORD8) for b in blocks):
                out.append([int(x) for x in chunk[idx]] + [next(i for i, m in enumerate(ORD8) if np.array_equal(b, m))
                                                            for b in blocks])
    return out


# ---------------------------------------------------------------- F9 arithmetic: xi with xi^2 = xi + 1 (x^2-x-1 irreducible over F3)
def f9_mul(a, b):
    # a = (a0, a1) = a0 + a1 xi ; xi^2 = xi + 1
    a0, a1 = a
    b0, b1 = b
    c0 = a0 * b0 + a1 * b1
    c1 = a0 * b1 + a1 * b0 + a1 * b1
    return (c0 % 3, c1 % 3)


def to_f9_code(G, struct):
    """express C as an F9-subspace of F9^10: coordinate v of a codeword = (a, b) with c_v = a e + b M_v e, e = (1, 0).
    Requires every M_v of trace 1, so that M_v^2 = M_v + 1 matches xi^2 = xi + 1 (multiplication by the same xi)."""
    assert all(int(np.trace(ORD8[k]) % 3) == 1 for k in struct)
    L = P.stabiliser(G) % 3
    Mv = [ORD8[k] for k in struct]
    rows = []
    for r in L:
        coords = []
        for v in range(10):
            c = r[2 * v:2 * v + 2]
            e = np.array([1, 0])
            Me = Mv[v] @ e % 3
            basis = np.stack([e, Me], 1)
            ab = (inv3(basis) @ c) % 3
            coords.append((int(ab[0]), int(ab[1])))
        rows.append(coords)
    return rows, Mv


def f9_rank(vectors):
    """rank over F9 of a list of vectors (each a list of (a0, a1)) -- Gaussian elimination with F9 arithmetic"""
    elems = [(a, b) for a in range(3) for b in range(3)]
    inv = {}
    for x in elems:
        for y in elems:
            if f9_mul(x, y) == (1, 0):
                inv[x] = y
    M = [list(v) for v in vectors]
    rank, ncol = 0, len(M[0]) if M else 0
    for c in range(ncol):
        piv = next((i for i in range(rank, len(M)) if M[i][c] != (0, 0)), None)
        if piv is None:
            continue
        M[rank], M[piv] = M[piv], M[rank]
        iv = inv[M[rank][c]]
        M[rank] = [f9_mul(iv, x) for x in M[rank]]
        for i in range(len(M)):
            if i != rank and M[i][c] != (0, 0):
                f = M[i][c]
                M[i] = [((x[0] - y[0]) % 3, (x[1] - y[1]) % 3) for x, y in zip(M[i], [f9_mul(f, z) for z in M[rank]])]
        rank += 1
    return rank


def schur_square_dim(rows):
    """the F9-span of the code (closure of the F3-rows under xi) and of its componentwise square"""
    xi = (0, 1)
    span = rows + [[f9_mul(xi, x) for x in r] for r in rows]
    k = f9_rank(span)
    prods = [[f9_mul(a, b) for a, b in zip(r, s)] for r, s in itertools.combinations_with_replacement(span, 2)]
    return k, f9_rank(prods)


# F9 = F3[xi]/(xi^2 - xi - 1); element a + b xi encoded as a + 3b
def enc(a, b): return (a % 3) + 3 * (b % 3)
def dec(x): return x % 3, x // 3
ADD = np.zeros((9, 9), int); MUL = np.zeros((9, 9), int)
for x in range(9):
    for y in range(9):
        a0, a1 = dec(x); b0, b1 = dec(y)
        ADD[x, y] = enc(a0 + b0, a1 + b1)
        MUL[x, y] = enc(a0 * b0 + a1 * b1, a0 * b1 + a1 * b0 + a1 * b1)
NEG = [enc(-dec(x)[0], -dec(x)[1]) for x in range(9)]
INV = [next(y for y in range(9) if MUL[x, y] == 1) if x else None for x in range(9)]
def pw(x, k):
    r = 1
    for _ in range(k): r = MUL[r, x]
    return r
FROB = [pw(x, 3) for x in range(9)]
NORM = [pw(x, 4) for x in range(9)]   # in {1, 2} for x != 0 (2 = -1)
xi = enc(0, 1)
def vadd(u, v): return [ADD[a, b] for a, b in zip(u, v)]
def smul(s, v): return [MUL[s, a] for a in v]
def matvec(A, v):
    out = []
    for row in A:
        acc = 0
        for a, b in zip(row, v): acc = ADD[acc, MUL[a, b]]
        out.append(acc)
    return out
def g_solve(Mcols, rhs):
    """solve sum_i x_i Mcols[i] = rhs over F9 (Mcols: 5 independent vectors)"""
    n = len(rhs)
    A = [[Mcols[j][i] for j in range(n)] + [rhs[i]] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c])
        A[c], A[p] = A[p], A[c]
        iv = INV[A[c][c]]; A[c] = [MUL[iv, x] for x in A[c]]
        for r in range(n):
            if r != c and A[r][c]:
                f = A[r][c]; A[r] = [ADD[x, NEG[MUL[f, y]]] for x, y in zip(A[r], A[c])]
    return [A[i][n] for i in range(n)]
def matinv(M):
    n = len(M)
    A = [list(M[i]) + [1 if i == j else 0 for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c])
        A[c], A[p] = A[p], A[c]
        iv = INV[A[c][c]]; A[c] = [MUL[iv, x] for x in A[c]]
        for r in range(n):
            if r != c and A[r][c]:
                f = A[r][c]; A[r] = [ADD[x, NEG[MUL[f, y]]] for x, y in zip(A[r], A[c])]
    return [row[n:] for row in A]
def matmul(A, B):
    n, m, k = len(A), len(B[0]), len(B)
    out = []
    for i in range(n):
        row = []
        for j in range(m):
            acc = 0
            for t in range(k): acc = ADD[acc, MUL[A[i][t], B[t][j]]]
            row.append(acc)
        out.append(row)
    return out
def normalize(v):
    i = next(k for k in range(len(v)) if v[k]); s = INV[v[i]]
    return tuple(MUL[s, a] for a in v), v[i]
eta = xi
pts = [[1, t, ADD[pw(t, 2), MUL[eta, pw(t, 6)]], pw(t, 3), pw(t, 4)] for t in range(9)] + [[0, 0, 0, 0, 1]]
key = {normalize(p)[0]: i for i, p in enumerate(pts)}
frame = list(range(6))
def frame_map(src, dst):
    """A with A src_i = lambda_i dst_i (i < 6), src/dst lists of 6 vectors"""
    a = g_solve(src[:5], src[5]); b = g_solve(dst[:5], dst[5])
    if 0 in a or 0 in b: return None
    S = [[smul(a[j], src[j])[i] for j in range(5)] for i in range(5)]   # columns a_j src_j
    D = [[smul(b[j], dst[j])[i] for j in range(5)] for i in range(5)]
    return matmul(D, matinv(S))


def arc_automorphisms():
    """projective (and Frobenius-twisted) maps of PG(4,9) preserving Glynn's arc: permutation and coordinate norms"""
    res = {"linear": [], "semilinear": []}
    for kind in res:
        src_pts = pts if kind == "linear" else [[FROB[x] for x in p_] for p_ in pts]
        src = [src_pts[i] for i in frame]
        for img in itertools.permutations(range(10), 6):
            A = frame_map(src, [pts[i] for i in img])
            if A is None:
                continue
            perm, lam, ok = [], [], True
            for v in range(10):
                nk, lead = normalize(matvec(A, src_pts[v]))
                if nk not in key:
                    ok = False
                    break
                u = key[nk]
                pn, plead = normalize(pts[u])
                lam.append(MUL[lead, INV[plead]])
                perm.append(u)
            if ok and len(set(perm)) == 10:
                res[kind].append((tuple(perm), tuple(NORM[l_] for l_ in lam)))
    return res


def porder(p_):
    q = list(range(10))
    k = 0
    while True:
        q = [p_[i] for i in q]
        k += 1
        if q == list(range(10)):
            return k


def hermitian_weightings(rows9):
    """weights u in {1,2}^10 with sum_i u_i x_i y_i^3 = 0 on the code (rows over F9, ints 0..8)"""
    out = []
    for u in itertools.product((1, 2), repeat=10):
        ok = True
        for r in rows9:
            for s_ in rows9:
                acc = 0
                for i in range(10):
                    acc = ADD[acc, MUL[u[i], MUL[r[i], FROB[s_[i]]]]]
                if acc:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            out.append(u)
    return out


def glynn_and_grs():
    import math
    rs = [[1 if j == 0 else pw(t, j) for t in range(9)] + [1 if j == 4 else 0] for j in range(5)]
    gl = [[c[i] for c in pts] for i in range(5)]
    aut = arc_automorphisms()
    lin = sum(1 for _, n_ in aut["linear"] if len(set(n_)) == 1)
    semi = sum(1 for _, n_ in aut["semilinear"] if len(set(n_)) == 1)
    perms = {p_ for k in aut for p_, _ in aut[k]}
    from collections import Counter
    orders = Counter(porder(p_) for p_ in perms)
    trip = {(p_[0], p_[1], p_[2]) for p_ in perms}
    order_lc = 4 * (lin + semi)
    n0 = 2642411520
    return dict(grs_hermitian_weightings=len(hermitian_weightings(rs)),
                glynn_hermitian_weightings=len(hermitian_weightings(gl)),
                arc_linear=len(aut["linear"]), arc_semilinear=len(aut["semilinear"]),
                aut_lc_order=order_lc, permutation_group_order=len(perms),
                permutation_element_orders={str(k): v for k, v in sorted(orders.items())},
                sharply_3_transitive=len(trip) == len(perms) == 720,
                orbit_size=24 ** 10 * math.factorial(10) // order_lc,
                master_f9_count_times_structures=6 ** 10 * n0 // 2,
                single_orbit=24 ** 10 * math.factorial(10) // order_lc == 6 ** 10 * n0 // 2)


def run():
    d = json.loads(P.TABU.read_text())
    res = dict(pass_id=11224, order8_elements=len(ORD8), states=[])
    for gi, w in enumerate(d["graphs"]):
        G = F.to_mat(np.array(w))
        st = f9_structures(G)
        traces = [sorted({int(np.trace(ORD8[k]) % 3) for k in s_}) for s_ in st]
        row = dict(graph=gi, f9_structures=len(st), uniform_trace=all(len(t) == 1 for t in traces))
        tr1 = [s_ for s_ in st if all(int(np.trace(ORD8[k]) % 3) == 1 for k in s_)]
        if tr1:
            rows, _ = to_f9_code(G, tr1[0])
            k, sq = schur_square_dim(rows)
            row.update(f9_dimension=k, schur_square_dim=sq, grs=bool(sq == 2 * k - 1))
        res["states"].append(row)
        print(row, flush=True)
    res["all_f9_linear"] = all(r["f9_structures"] > 0 for r in res["states"])
    res["all_grs"] = all(r.get("grs", False) for r in res["states"])
    res["structures_per_state"] = sorted({r["f9_structures"] for r in res["states"]})
    res["glynn"] = glynn_and_grs()
    return res


def main():
    res = run()
    OUT.write_text(json.dumps(res, indent=1))
    print(json.dumps({k: v for k, v in res.items() if k != "states"}, indent=1))


if __name__ == "__main__":
    main()
