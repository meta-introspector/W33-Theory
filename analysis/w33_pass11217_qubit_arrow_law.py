#!/usr/bin/env python3
"""Pass 11217: the arrow law A(S) = n - c(S) holds for qubits -- characteristic 2, every number of qubits.

Pass 11207 proved A(S) = n - c(S) for every odd local dimension and left qubits (Sp(2n,2)) as a conjecture verified
class by class through n = 5 (Pass 11208).  The odd-characteristic proof breaks in characteristic 2 because the form
b(x, y) = omega(x, Sy) on a bi-Lagrangian can be alternating, and then it has no orthogonal basis.

THE CRITERION (characteristic 2; proved in the note, checked here).  Write q(x) = omega(x, Sx).  Its polar form is
omega(x, (S + S^-1) y), so q is additive on every bi-Lagrangian L (Lagrangian for omega and for omega(., (S+S^-1).)).
  A split of V into n planes, each meeting its image, exists  <=>  some bi-Lagrangian L has
  (a) L n SL spanned by eigenvectors of S (over F_2: fixed vectors), and (b) q|L != 0 or L = L n SL.
(<=: diagonalise b on a complement of the radical -- possible because b is non-alternating -- and pair each radical
eigenvector with a dual vector.  =>: one vector per plane, chosen in P with SP's line, spans such an L.)  The
criterion is additive over orthogonal S-invariant sums.  Since any split has at most c invariant planes (lower bound
A >= n - c, characteristic-free), the law reduces to the orthogonally indecomposable S-modules of dimension >= 4.
These are:
  (i)   U + U* with U = F[x]/(f^m), f != f*: the graph of multiplication by a unit g with l(tau g) != 0 is transverse;
  (ii)  F[x]/(f^m), f = f* != x+1 (the extension R/R_0 is etale): the R_0-line R_0 eta is transverse, and q != 0 for a
        suitable unit eta because the norm is onto and units span R_0;
  (iii) V(2k), one Jordan block (Hesselink): in R = F[t]/(t^2k), sigma(t) = t/(1+t), s = t sigma(t) = t^2/(1+t),
        P = F[s]; l vanishes on P, L = P eta is transverse, and q|L = lambda(a^2 N(eta)) is nonzero for some eta when
        k >= 3 because the span of {a^2 N(eta)} contains every s^j with j != 1 and lambda(s^(k-1)) != 0;
  (iv)  W(k) = J_k + J_k*: L_g = P(1, 0) + P(t, g) is bi-Lagrangian, is Lagrangian iff sigma(g) in P-perp, has radical
        ann(g), and has q != 0 iff sigma(g) is not in Q-perp (Q = span{c^2 s t}); a g of valuation 0 (k even) or 1
        (k odd) exists for k >= 4 because a vector space is never the union of two proper subspaces;
  (v)   the three finite exceptions V(4), W(2), W(3), found by exhaustive search.
So A(S) = n - c(S) for every S in Sp(2n, 2) and every n.  Qubit closed form (checked on every class, n <= 5):
        c(S) = m_1(1)/2 + chi_2 * m_2(1) + m_1(x^2+x+1),
m_s(f) the number of Jordan blocks of size s at f and chi_2 = 1 iff some v in ker(S+1)^2 has omega(v, (S+1)v) = 1
(Hesselink's index: whether the size-2 blocks can be taken as transvection planes V(2)).
"""
from __future__ import annotations

import itertools
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11208_arrow_universality as U  # noqa: E402

OUT = ROOT / "data" / "w33_pass11217_qubit_arrow_law.json"


# ---------------------------------------------------------------- F_2 linear algebra
def rref(M):
    M = np.array(M, np.int64) % 2
    if M.ndim == 1:
        M = M[None]
    r = 0
    for c in range(M.shape[1]):
        piv = [i for i in range(r, M.shape[0]) if M[i, c]]
        if not piv:
            continue
        M[[r, piv[0]]] = M[[piv[0], r]]
        for i in range(M.shape[0]):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        r += 1
        if r == M.shape[0]:
            break
    return M[:r]


def rank(M):
    return len(rref(M)) if len(M) else 0


def nullspace(M):
    """basis (rows) of {x : M x = 0}"""
    M = np.array(M, np.int64) % 2
    R = rref(M)
    piv = [int(np.argmax(row)) for row in R]
    free = [c for c in range(M.shape[1]) if c not in piv]
    out = []
    for f in free:
        x = np.zeros(M.shape[1], np.int64)
        x[f] = 1
        for row, p in zip(R, piv):
            x[p] = row[f]
        out.append(x)
    return np.array(out, np.int64).reshape(len(out), M.shape[1])


def intersect(A, B):
    K = nullspace(np.vstack([A, B]).T)
    if len(K) == 0:
        return np.zeros((0, A.shape[1]), np.int64)
    return rref((K[:, :len(A)] @ A) % 2)


def inverse(M):
    D = len(M)
    R = rref(np.hstack([M % 2, np.eye(D, dtype=np.int64)]))
    assert (R[:, :D] == np.eye(D, dtype=np.int64)).all()
    return R[:, D:]


def mpow(M, e):
    R = np.eye(len(M), dtype=np.int64)
    for _ in range(e):
        R = (R @ M) % 2
    return R


def dimker(M):
    return len(M) - rank(M)


# ---------------------------------------------------------------- the criterion and the split it builds
def criterion(S, O, L):
    """(bi-Lagrangian, radical dim, radical fixed, q nonzero) for a Lagrangian L (rows) of (F_2^D, O)"""
    Sinv = (inverse(O) @ S.T @ O) % 2
    T = (S + Sinv) % 2
    bi = not ((L @ O @ T @ L.T) % 2).any()
    R = intersect(L, (L @ S.T) % 2)
    fixed = len(R) == 0 or not ((R @ S.T + R) % 2).any()
    q = any(int(x @ O @ S @ x) % 2 for x in L)
    return bi, len(R), fixed, q


def good(S, O, L):
    bi, r, fixed, q = criterion(S, O, L)
    return bi and fixed and (q or r == len(L))


def build_split(S, O, L):
    """the split of the criterion: cyclic planes span(l, Sl) for a b-orthonormal basis of a complement of the radical,
    and planes span(e, x) through the fixed radical vectors e"""
    D = len(S)
    R = intersect(L, (L @ S.T) % 2)
    basis = [r for r in R]
    comp = []
    for v in L:
        if rank(np.array(basis + comp + [v])) > len(basis) + len(comp):
            comp.append(v)

    def b(x, y):
        return int(x @ O @ S @ y) % 2
    # b-orthonormal basis of span(comp): q(x) = b(x, x) is additive in characteristic 2
    ortho, W = [], [np.array(c) for c in comp]
    while W:
        w = next((x for x in W if b(x, x)), None)
        if w is not None:
            ortho.append(w)
            W = [x for x in W if x is not w]
            W = [(x + b(x, w) * w) % 2 for x in W]          # b(w, w) = 1
            W = list(rref(np.array(W))) if W and rank(np.array(W)) else []
            continue
        # alternating remainder: trade a hyperbolic pair (u, v) against an earlier x0 with b(x0, x0) = 1
        u = W[0]
        v = next(x for x in W[1:] if b(u, x))
        x0 = ortho.pop()
        ortho += [(x0 + u) % 2, (x0 + v) % 2, (x0 + u + v) % 2]
        W = [x for x in W if x is not u and x is not v]
        W = [(x + b(x, v) * u + b(x, u) * v) % 2 for x in W]
        W = list(rref(np.array(W))) if W and rank(np.array(W)) else []
    planes = [np.stack([l, (S @ l) % 2]) for l in ortho]
    # fixed-vector planes inside the omega-complement of the cyclic planes
    if planes:
        C = np.vstack(planes)
        Wp = nullspace((C @ O) % 2)
    else:
        Wp = np.eye(D, dtype=np.int64)
    xs = []
    for j, e in enumerate(basis):
        A = np.array([(Wp @ O @ ei) % 2 for ei in basis]).T   # columns: omega(w_row, e_i)
        rhs = np.array([1 if i == j else 0 for i in range(len(basis))])
        sol = solve(A, rhs)
        x = (sol @ Wp) % 2
        for i, xi in enumerate(xs):
            x = (x + (int(xi @ O @ x) % 2) * basis[i]) % 2
        xs.append(x)
    planes += [np.stack([e, x]) for e, x in zip(basis, xs)]
    return planes


def solve(A, rhs):
    """one solution y of A^T-style system: find y with y @ Wp-images... A (m x k), rhs (k): y (m) with y @ A = rhs"""
    m, k = A.shape
    M = np.hstack([A.T % 2, rhs[:, None] % 2])
    R = rref(M)
    y = np.zeros(m, np.int64)
    for row in R:
        p = int(np.argmax(row[:m])) if row[:m].any() else None
        if p is None:
            assert not row[m], "inconsistent"
            continue
        y[p] = row[m]
    assert ((y @ A) % 2 == rhs % 2).all()
    return y


def verify_split(S, O, planes):
    D = len(S)
    ok = len(planes) == D // 2 and rank(np.vstack(planes)) == D
    for i, P in enumerate(planes):
        ok &= int(P[0] @ O @ P[1]) % 2 == 1
        SP = (P @ S.T) % 2
        ok &= rank(np.vstack([P, SP])) <= 3                  # P meets its image
        for Q in planes[i + 1:]:
            ok &= not ((P @ O @ Q.T) % 2).any()
    return bool(ok)


# ---------------------------------------------------------------- truncated polynomial rings over F_2
def pmul(a, b, mod):
    """product of coefficient vectors modulo the monic polynomial mod (low degree first)"""
    d = len(mod) - 1
    c = np.zeros(2 * d, np.int64)
    for i in np.flatnonzero(a):
        c[i:i + d] ^= b
    for i in range(2 * d - 1, d - 1, -1):
        if c[i]:
            c[i - d:i + 1] ^= mod
    return c[:d]


def mult_matrix(a, mod):
    d = len(mod) - 1
    M = np.zeros((d, d), np.int64)
    for i in range(d):
        e = np.zeros(d, np.int64)
        e[i] = 1
        M[:, i] = pmul(a, e, mod)
    return M


def poly(coeffs_high_to_low):
    return np.array(coeffs_high_to_low[::-1], np.int64)


def ring_inverse(a, mod):
    M = mult_matrix(a, mod)
    one = np.zeros(len(mod) - 1, np.int64)
    one[0] = 1
    return (inverse(M) @ one) % 2


def span_ok(rows):
    return rref(np.array(rows))


# ---------------------------------------------------------------- (iii) V(2k) and (iv) W(k): the t-adic models
def t_model(K):
    """R = F_2[t]/(t^K): mult matrices of t, x = 1 + t, sigma (t -> t/(1+t)), and s = t^2/(1+t)"""
    mod = np.zeros(K + 1, np.int64)
    mod[K] = 1
    e = [np.eye(K, dtype=np.int64)[i] for i in range(K)]
    one, t = e[0], (e[1] if K > 1 else np.zeros(K, np.int64))
    inv1pt = np.ones(K, np.int64)                              # 1/(1+t) = sum t^i
    sigt = pmul(t, inv1pt, mod)
    Sig = np.zeros((K, K), np.int64)
    p = one.copy()
    for i in range(K):
        Sig[:, i] = p
        p = pmul(p, sigt, mod)
    s = pmul(pmul(t, t, mod), inv1pt, mod)
    return dict(mod=mod, one=one, t=t, x=(one + t) % 2, Sig=Sig, s=s,
                X=mult_matrix((one + t) % 2, mod))


def spow(m, j):
    out = m["one"].copy()
    for _ in range(j):
        out = pmul(out, m["s"], m["mod"])
    return out


def V_piece(k, ell):
    """V(2k): R = F_2[t]/(t^2k), omega(a, b) = ell(a sigma(b))"""
    m = t_model(2 * k)
    D = 2 * k
    O = np.zeros((D, D), np.int64)
    for i in range(D):
        for j in range(D):
            ei = np.eye(D, dtype=np.int64)[i]
            O[i, j] = int(ell @ pmul(ei, m["Sig"][:, j], m["mod"])) % 2
    return m, m["X"], O


def V_ells(k):
    """all ell with ell(P) = 0 and ell(t^(2k-1)) = 1 (alternating, nondegenerate)"""
    m = t_model(2 * k)
    P = np.array([spow(m, j) for j in range(k)])
    K = nullspace(P)                                           # ell with P ell = 0
    out = []
    for c in itertools.product(range(2), repeat=len(K)):
        ell = (np.array(c) @ K) % 2 if len(K) else np.zeros(2 * k, np.int64)
        if ell[2 * k - 1]:
            out.append(ell)
    return out


def check_V(k):
    res = []
    for ell in V_ells(k):
        m, S, O = V_piece(k, ell)
        assert (O == O.T).all() and not np.diag(O).any() and rank(O) == 2 * k
        assert ((S.T @ O @ S) % 2 == O).all()
        P = [spow(m, j) for j in range(k)]
        # the proof's candidates: eta = 1 and eta = 1 + t s^j
        cands = [m["one"]] + [(m["one"] + pmul(m["t"], spow(m, j), m["mod"])) % 2 for j in range(k)]
        hit = None
        for eta in cands:
            L = rref(np.array([pmul(p, eta, m["mod"]) for p in P]))
            if len(L) != k:
                continue
            bi, r, fixed, q = criterion(S, O, L)
            assert bi and r == 0, "P eta must be a transverse bi-Lagrangian"
            if q:
                hit = L
                break
        ok = hit is not None and verify_split(S, O, build_split(S, O, hit))
        res.append(ok)
    return dict(k=k, forms=len(res), all_split=all(res))


def W_piece(k):
    m = t_model(k)
    D = 2 * k
    O = np.zeros((D, D), np.int64)
    for i in range(k):
        for j in range(k):
            ei = np.eye(k, dtype=np.int64)[i]
            v = pmul(ei, m["Sig"][:, j], m["mod"])[k - 1]
            O[i, k + j] ^= v
            O[k + j, i] ^= v
    S = np.zeros((D, D), np.int64)
    S[:k, :k] = m["X"]
    S[k:, k:] = m["X"]
    return m, S, O


def check_W(k):
    m, S, O = W_piece(k)
    assert (O == O.T).all() and not np.diag(O).any() and rank(O) == 2 * k and ((S.T @ O @ S) % 2 == O).all()
    np_ = (k + 1) // 2
    P = np.array([spow(m, j) for j in range(np_)])
    # P-perp (for <u, v> = l(uv)) and Q = span{c^2 s t}
    Pp = nullspace(np.array([[pmul(a, np.eye(k, dtype=np.int64)[i], m["mod"])[k - 1] for i in range(k)] for a in P]))
    st = pmul(m["s"], m["t"], m["mod"])
    Q = np.array([pmul(pmul(c, c, m["mod"]), st, m["mod"]) for c in P])
    val_needed = 0 if k % 2 == 0 else 1
    # containments that the proof shows fail
    Pp_in_t = all(not y[:1].any() for y in Pp)
    Pp_in_t2 = all(not y[:2].any() for y in Pp)
    Q_in_P = rank(np.vstack([P, Q])) == rank(P)
    found = None
    for c in itertools.product(range(2), repeat=len(Pp)):
        y = (np.array(c) @ Pp) % 2
        if not y.any():
            continue
        val = int(np.argmax(y))
        if val != val_needed:
            continue
        if not any(pmul(qq, y, m["mod"])[k - 1] for qq in Q):      # y in Q-perp: q vanishes on L
            continue
        found = y
        break
    ok = False
    if found is not None:
        g = (m["Sig"] @ found) % 2                               # sigma is an involution
        zero = np.zeros(k, np.int64)
        rows = [np.concatenate([a, zero]) for a in P] + \
               [np.concatenate([pmul(a, m["t"], m["mod"]), pmul(a, g, m["mod"])]) for a in P]
        L = rref(np.array(rows))
        bi, r, fixed, q = criterion(S, O, L)
        ok = len(L) == k and not ((L @ O @ L.T) % 2).any() and bi and fixed and q and r == (k % 2) \
            and verify_split(S, O, build_split(S, O, L))
    return dict(k=k, Pperp_contains_unit=not Pp_in_t, Pperp_not_in_t2=not Pp_in_t2, Q_not_in_P=not Q_in_P,
                construction_split=bool(ok))


# ---------------------------------------------------------------- (i) f != f* and (ii) f = f* != x+1
def check_hyperbolic(f_high, m_exp):
    """U + U*, U = F_2[x]/(f^m): V = R + R, S = (x, x^-1), omega((a,b),(a',b')) = l(ab') + l(a'b), l = top coeff"""
    f = poly(f_high)
    mod = np.array([1], np.int64)
    for _ in range(m_exp):
        mod = np.convolve(mod, f) % 2
    d = len(mod) - 1
    e = np.eye(d, dtype=np.int64)
    x = e[1] if d > 1 else np.zeros(d, np.int64)
    X = mult_matrix(x, mod)
    Xi = inverse(X)
    D = 2 * d
    ell = e[d - 1]
    O = np.zeros((D, D), np.int64)
    for i in range(d):
        for j in range(d):
            v = int(ell @ pmul(e[i], e[j], mod)) % 2
            O[i, d + j] ^= v
            O[d + j, i] ^= v
    S = np.zeros((D, D), np.int64)
    S[:d, :d], S[d:, d:] = X, Xi
    assert rank(O) == D and ((S.T @ O @ S) % 2 == O).all()
    tau = (x + (Xi @ e[0]) % 2) % 2
    for gc in itertools.product(range(2), repeat=d):
        g = np.array(gc)
        if rank(mult_matrix(g, mod)) < d or not int(ell @ pmul(tau, g, mod)) % 2:
            continue
        G = mult_matrix(g, mod)
        L = rref(np.hstack([np.eye(d, dtype=np.int64), G.T % 2]))
        bi, r, fixed, q = criterion(S, O, L)
        return dict(f="".join(map(str, f_high)), m=m_exp, dim=D, transverse=r == 0, bi=bi, q=q,
                    split=verify_split(S, O, build_split(S, O, L)))
    return dict(f="".join(map(str, f_high)), m=m_exp, dim=D, split=False)


def check_selfdual(f_high, m_exp):
    """R = F_2[x]/(f^m), f = f*, omega(a, b) = l(a sigma(b)), sigma(x) = x^-1, l(R_0) = 0; L = R_0 eta"""
    f = poly(f_high)
    mod = np.array([1], np.int64)
    for _ in range(m_exp):
        mod = np.convolve(mod, f) % 2
    d = len(mod) - 1
    e = np.eye(d, dtype=np.int64)
    X = mult_matrix(e[1], mod)
    Xi = inverse(X)
    xinv = (Xi @ e[0]) % 2
    Sig = np.zeros((d, d), np.int64)
    p = e[0].copy()
    for i in range(d):
        Sig[:, i] = p
        p = pmul(p, xinv, mod)
    tau = (e[1] + xinv) % 2
    R0 = [e[0]]
    for _ in range(d):
        R0.append(pmul(R0[-1], tau, mod))
    R0 = rref(np.array(R0))
    assert 2 * len(R0) == d
    # choose l with l(R_0) = 0 and omega nondegenerate
    K = nullspace(R0)
    for c in itertools.product(range(2), repeat=len(K)):
        ell = (np.array(c) @ K) % 2
        if not ell.any():
            continue
        O = np.array([[int(ell @ pmul(e[i], Sig[:, j], mod)) % 2 for j in range(d)] for i in range(d)])
        if rank(O) == d:
            break
    assert (O == O.T).all() and not np.diag(O).any() and ((X.T @ O @ X) % 2 == O).all()
    for ec in itertools.product(range(2), repeat=d):
        eta = np.array(ec)
        if rank(mult_matrix(eta, mod)) < d:
            continue
        L = rref(np.array([pmul(r, eta, mod) for r in R0]))
        bi, r, fixed, q = criterion(X, O, L)
        if bi and r == 0 and q:
            return dict(f="".join(map(str, f_high)), m=m_exp, dim=d, transverse=True, q=True,
                        split=verify_split(X, O, build_split(X, O, L)))
    return dict(f="".join(map(str, f_high)), m=m_exp, dim=d, split=False)


# ---------------------------------------------------------------- (v) exceptions and the class-level checks
def lagrangians(O):
    D = len(O)
    n = D // 2
    L = []
    for c in itertools.product(range(2), repeat=D):
        v = np.array(c)
        if v.any() and all((v @ O @ w) % 2 == 0 for w in L) and rank(np.array(L + [v])) == len(L) + 1:
            L.append(v)
            if len(L) == n:
                break
    L0 = rref(np.array(L))
    tv = [np.array(x) for x in itertools.product(range(2), repeat=D) if any(x)]
    seen = {L0.tobytes(): L0}
    fr = [L0]
    while fr:
        nf = []
        for L in fr:
            for v in tv:
                L2 = rref((L + np.outer((L @ O @ v) % 2, v)) % 2)
                if L2.tobytes() not in seen:
                    seen[L2.tobytes()] = L2
                    nf.append(L2)
        fr = nf
    return list(seen.values())


def exception(name, S, O):
    Ls = lagrangians(O)
    goods = [L for L in Ls if good(S, O, L)]
    ok = bool(goods) and verify_split(S, O, build_split(S, O, goods[0]))
    return dict(piece=name, lagrangians=len(Ls), criterion_lagrangians=len(goods), split=ok)


def c_formula(S, n):
    D = 2 * n
    I = np.eye(D, dtype=np.int64)
    J = U.form(n, 2)
    N = (S + I) % 2
    k1, k2, k3 = dimker(N), dimker(mpow(N, 2)), dimker(mpow(N, 3))
    m1, m2 = 2 * k1 - k2, (k2 - k1) - (k3 - k2)
    K2 = nullspace(mpow(N, 2))
    chi2 = 0
    for c in itertools.product(range(2), repeat=len(K2)):
        v = (np.array(c) @ K2) % 2 if len(K2) else None
        if v is not None and v.any() and (v @ J @ N @ v) % 2:
            chi2 = 1
            break
    F = (S @ S + S + I) % 2
    f1, f2 = dimker(F), dimker(mpow(F, 2))
    return m1 // 2 + chi2 * m2 + (f1 // 2 - (f2 - f1) // 2)


def class_checks(n):
    cls, head = U.load_classes(U.GAP_SMALL if n < 5 else U.GAP_BIG, 2, n)
    Pl = U.all_planes(n, 2)
    J = U.form(n, 2)
    Ls = lagrangians(J) if n <= 4 else None
    rows = []
    for c in cls:
        S = c["S"] % 2
        A, cc, _ = U.arrow_and_c(S, Pl, n, 2) if n <= 4 else (None, None, None)
        if n == 5:
            cc = c_brute(S, Pl, n)
        row = dict(order=c["order"], size=c["size"], c=cc, c_formula=c_formula(S, n))
        if n <= 4:
            row["A"] = A
            if cc == 0:
                gl = [L for L in Ls if good(S, J, L)]
                row["criterion"] = bool(gl)
                row["split_built"] = bool(gl) and verify_split(S, J, build_split(S, J, gl[0]))
        rows.append(row)
    return rows


def c_brute(S, Pl, n):
    J = U.form(n, 2)
    d = U.dvals(S, Pl, 2)
    G2 = Pl[d == 2]
    G2J = np.einsum('pai,ab->pib', G2, J)
    best = [0]

    def orth(i, c):
        if len(c) == 0:
            return c
        W = np.einsum('ib,pbj->pij', G2J[i], G2[c]) % 2
        return c[~W.reshape(len(c), 4).any(1)]

    def dfs(c, k):
        best[0] = max(best[0], k)
        if best[0] == n or k + len(c) <= best[0]:
            return
        for t, i in enumerate(c):
            if k + len(c) - t <= best[0]:
                return
            dfs(orth(i, c[t + 1:]), k + 1)
            if best[0] == n:
                return
    dfs(np.arange(len(G2)), 0)
    return best[0]


def verifier_control(trials=200, seed=11217):
    """random splits of W(4) must be rejected by the same verifier that accepts the constructed ones"""
    _, S, O = W_piece(4)
    rng = np.random.default_rng(seed)
    accepted = 0
    for _ in range(trials):
        planes, W = [], np.eye(8, dtype=np.int64)
        while len(planes) < 4:
            B = rref(W)
            u = (rng.integers(0, 2, len(B)) @ B) % 2
            if not u.any():
                continue
            v = next((c for c in ((rng.integers(0, 2, len(B)) @ B) % 2 for _ in range(30)) if (u @ O @ c) % 2), None)
            if v is None:
                continue
            planes.append(np.stack([u, v]))
            W = nullspace((np.vstack(planes) @ O) % 2)
        accepted += verify_split(S, O, planes)
    return dict(trials=trials, accepted=int(accepted))


def run(include_n5=True):
    t0 = time.time()
    res = dict(pass_id=11217)
    res["V2k"] = [check_V(k) for k in range(3, 8)]
    res["Wk"] = [check_W(k) for k in range(4, 13)]
    res["hyperbolic_f_ne_fstar"] = [check_hyperbolic(f, m) for f, m in
                                    [([1, 0, 1, 1], 1), ([1, 0, 1, 1], 2), ([1, 0, 0, 1, 1], 1),
                                     ([1, 0, 0, 1, 0, 1], 1), ([1, 1, 0, 0, 1], 2)]]
    res["selfdual_f_eq_fstar"] = [check_selfdual(f, m) for f, m in
                                  [([1, 1, 1], 2), ([1, 1, 1], 3), ([1, 1, 1, 1, 1], 1), ([1, 1, 1, 1, 1], 2),
                                   ([1, 0, 0, 1, 0, 0, 1], 1), ([1, 1, 1], 4)]]
    mV2, SV2, OV2 = V_piece(2, V_ells(2)[0])
    exc = [exception("V(4)", SV2, OV2)]
    for k in (2, 3):
        _, S, O = W_piece(k)
        exc.append(exception(f"W({k})", S, O))
    res["exceptions"] = exc
    res["verifier_control"] = verifier_control()
    res["classes"] = {}
    for n in ((2, 3, 4, 5) if include_n5 else (2, 3, 4)):
        rows = class_checks(n)
        summ = dict(classes=len(rows), c_formula_exact=all(r["c"] == r["c_formula"] for r in rows))
        if n <= 4:
            summ["law"] = all(r["A"] == n - r["c"] for r in rows)
            summ["c0_classes"] = sum(r["c"] == 0 for r in rows)
            summ["criterion_on_every_c0_class"] = all(r.get("criterion", True) for r in rows)
            summ["split_built_on_every_c0_class"] = all(r.get("split_built", True) for r in rows)
        res["classes"][str(n)] = summ
    res["all_constructions_split"] = (all(r["all_split"] for r in res["V2k"])
                                      and all(r["construction_split"] for r in res["Wk"])
                                      and all(r["split"] for r in res["hyperbolic_f_ne_fstar"])
                                      and all(r["split"] for r in res["selfdual_f_eq_fstar"])
                                      and all(r["split"] for r in res["exceptions"]))
    res["W_containments_fail"] = all(r["Q_not_in_P"] and (r["Pperp_contains_unit"] if r["k"] % 2 == 0
                                                          else r["Pperp_not_in_t2"]) for r in res["Wk"])
    res["seconds"] = round(time.time() - t0, 1)
    return res


def main():
    res = run()
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
