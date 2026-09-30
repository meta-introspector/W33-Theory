#!/usr/bin/env python3
"""Pass 11199: the intrinsic arrow of time for every number of qutrits --  A(S) = n - c(S).

A(S) (Pass 11183) is the information a Clifford tick S in Sp(2n,3) exports per tick in the best subsystem split,
minimised over all splits; Pass 11188: A(S) = 2n - max sum_P dim(P cap SP) over orthogonal families of nondegenerate
planes.  Let c(S) be the maximal number of mutually orthogonal S-INVARIANT nondegenerate planes (qutrits the dynamics
can keep to themselves).
THEOREM (every n).   A(S) = n - c(S).   In particular A(S) <= n, A(S) = n iff S fixes no nondegenerate plane, and A
never takes the value 1 (c = n - 1 forces c = n).
LOWER BOUND (elementary).  In any split, a planes are invariant (d = 2, a <= c) and the others have d <= 1, so
E = 2n - sum d >= 2n - 2a - (n - a) = n - a >= n - c.
UPPER BOUND.  Remove a maximum orthogonal family of invariant planes; the S-invariant complement W has none.  A is
subadditive over orthogonal S-invariant pieces, so it suffices that every tick with no invariant nondegenerate plane has
a HALF-MOVING SPLIT (all planes meeting their image).  Split V orthogonally and S-invariantly:
  * by the characteristic polynomial into components V_f + V_f* (f != f* irreducible) and V_f (f = f*);
  * a component with f = x -+ 1 is +-(a unipotent), which splits orthogonally into V(2m) (one Jordan block of even size)
    and W(k) (two dual Jordan blocks of odd size k); a component with f = f* of degree 2d is a Hermitian space over
    E = F_{3^{2d}} on which S = mu U, U unipotent unitary, and splits into single unitary Jordan blocks.
Each indecomposable type has an explicit construction (checked here over a range of sizes):
  (delta) f != f*:  S = A + A^-T on U + U*; L = {(x, Phi x)} with Phi symmetric, Phi A = A^T Phi (Taussky-Zassenhaus);
          L cap SL = 0 because A^2 - 1 is invertible.
  (alpha) V(2m): N = (u - 1)(u + 1)^-1 is skew and nilpotent; L = <v, N^2 v, ..., N^{2m-2} v> (v cyclic) is Lagrangian,
          T-stable (T = u + u^-1 is even in N), and L cap uL = 0 by parity.
  (beta)  W(k), k odd: reduce by the eigenvector e (u e = e): e^perp / e = J_{k-1} + J_{k-1}^-T, where the parity
          Lagrangian <N^{even} v, N^{even} w> works; its lift L has L cap uL = <e>, and the RADICAL LEMMA below gives the
          split with one plane through e.
  (gamma) mu J_k unitary: the real form V_R = F-span{eps_j N^j v} (eps_j = 1 for even j, kappa (kappa-bar = -kappa) for
          odd j) is Lagrangian, T-stable, and V_R cap S V_R = 0 because mu is not in F.
BI-LAGRANGIAN / RADICAL LEMMA.  If L is Lagrangian for omega and for omega_T (T = S + S^-1), then b(x, y) = omega(x, Sy)
is symmetric on L with radical L cap SL.  If the radical is 0 (or spanned by an eigenvector e), an orthogonal basis x_i
of b (modulo e) gives planes <x_i, S x_i>, mutually orthogonal and nondegenerate (plus the plane through e orthogonal to
all of them), each meeting its image.
CLOSED FORM.  Jordan types add over orthogonal S-invariant sums, and an invariant nondegenerate plane is one of three
2-dimensional indecomposables (a pair of size-1 blocks at lambda = +-1, a size-2 block at lambda = +-1, or a size-1
F9-block of x^2 + 1), so
        c(S) = sum_{lambda = +-1} ( m_1(lambda)/2 + m_2(lambda) ) + m_1^{F9}(x^2 + 1),
with m_s the number of Jordan blocks of size s (from ranks of (S - lambda)^s and (S^2 + 1)^s).  So A(S) is computable
in polynomial time for any n (function `arrow`).
CHECKS HERE.
  * A = n - c(S) on all 20 + 74 + 278 conjugacy classes of PSp(4,3), PSp(6,3), PSp(8,3) (A from Pass 11188).
  * Every construction verified by building the split and checking it: V(2m) m <= 8 (both form classes), W(k) k <= 11,
    unitary blocks over F9 (k <= 5) and F81 (k <= 3), non-self-reciprocal pairs up to (x^3+2x+1)^2.
  * Random ticks at n = 5, 6: A computed EXACTLY (a split attaining the lower bound n - c) -- it equals n - c every time.
"""
from __future__ import annotations

import itertools
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
from sympy import Poly, symbols

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

OUT = ROOT / "data" / "w33_pass11199_arrow_theorem.json"
_X = symbols('x')


# ---------------- F3 linear algebra ----------------
def rref(A):
    A = np.array(A, np.int64) % 3
    m, n = A.shape
    piv, r = [], 0
    for c in range(n):
        p = next((i for i in range(r, m) if A[i, c]), None)
        if p is None:
            continue
        A[[r, p]] = A[[p, r]]
        A[r] = (A[r] * pow(int(A[r, c]), -1, 3)) % 3
        for i in range(m):
            if i != r and A[i, c]:
                A[i] = (A[i] - A[i, c] * A[r]) % 3
        piv.append(c)
        r += 1
    return A[:r], piv


def rank(A):
    return len(rref(A)[1]) if np.size(A) else 0


def nullspace(A):
    A = np.array(A, np.int64) % 3
    R, piv = rref(A)
    n = A.shape[1]
    out = []
    for f in (c for c in range(n) if c not in piv):
        v = np.zeros(n, np.int64)
        v[f] = 1
        for i, c in enumerate(piv):
            v[c] = (-R[i, f]) % 3
        out.append(v)
    return np.array(out, np.int64).reshape(-1, n)


def inv(A):
    n = len(A)
    R, piv = rref(np.concatenate([np.asarray(A) % 3, np.eye(n, dtype=np.int64)], 1))
    assert piv[:n] == list(range(n)), 'singular'
    return R[:, n:] % 3


def intersect_dim(A, B):
    return rank(A) + rank(B) - rank(np.concatenate([A, B], 1))


def form_std(n):
    J = np.zeros((2 * n, 2 * n), np.int64)
    for k in range(n):
        J[2 * k, 2 * k + 1], J[2 * k + 1, 2 * k] = 1, 2
    return J


def orthogonal_basis(b):
    """invertible C with C^T b C diagonal (b symmetric mod 3); returns C and the diagonal"""
    B = np.array(b, np.int64) % 3
    r = len(B)
    C = np.eye(r, dtype=np.int64)
    for k in range(r):
        idx = list(range(k, r))
        if not any(B[i, i] for i in idx):
            pq = next(((p, q) for p in idx for q in idx if p < q and B[p, q]), None)
            if pq is None:
                break
            p, q = pq
            B[p] = (B[p] + B[q]) % 3
            B[:, p] = (B[:, p] + B[:, q]) % 3
            C[:, p] = (C[:, p] + C[:, q]) % 3
        j = next(i for i in idx if B[i, i])
        B[[k, j]] = B[[j, k]]
        B[:, [k, j]] = B[:, [j, k]]
        C[:, [k, j]] = C[:, [j, k]]
        piv = pow(int(B[k, k]), -1, 3)
        for i in range(k + 1, r):
            f = (B[i, k] * piv) % 3
            B[i] = (B[i] - f * B[k]) % 3
            B[:, i] = (B[:, i] - f * B[:, k]) % 3
            C[:, i] = (C[:, i] - f * C[:, k]) % 3
    return C % 3, np.diag(B) % 3


# ---------------- the radical lemma and the split check ----------------
def split_from_lagrangian(S, L, Om):
    D = len(Om)
    n = D // 2
    assert rank(L) == n and not ((L.T @ Om @ L) % 3).any(), 'not Lagrangian'
    T = (S + inv(S)) % 3
    assert not ((L.T @ Om @ T @ L) % 3).any(), 'not omega_T-isotropic'
    b = (L.T @ Om @ S @ L) % 3
    assert ((b - b.T) % 3 == 0).all()
    C, dg = orthogonal_basis(b)
    X = (L @ C) % 3
    planes = [np.stack([X[:, i], (S @ X[:, i]) % 3], 1) for i in range(n) if dg[i]]
    rad = [X[:, i] for i in range(n) if not dg[i]]
    if rad:
        assert len(rad) == 1, 'radical of dimension > 1'
        e = rad[0]
        assert not (((S @ e) - e) % 3).any() or not (((S @ e) + e) % 3).any(), 'radical is not an eigenvector'
        Pm = np.concatenate(planes, 1)
        W = nullspace((Pm.T @ Om) % 3).T
        assert W.shape[1] == 2
        planes.append(W)
    return planes


def check_split(S, planes, Om):
    n = len(Om) // 2
    assert len(planes) == n
    for P in planes:
        assert rank(P) == 2 and ((P.T @ Om @ P) % 3).any(), 'degenerate plane'
    for P, Q in itertools.combinations(planes, 2):
        assert not ((P.T @ Om @ Q) % 3).any(), 'not orthogonal'
    ds = [intersect_dim(P, (S @ P) % 3) for P in planes]
    return sum(2 - d for d in ds)


def cayley_inv(u):
    I = np.eye(len(u), dtype=np.int64)
    return ((u - I) @ inv((u + I) % 3)) % 3


def mpow(N, k):
    out = np.eye(len(N), dtype=np.int64)
    for _ in range(k):
        out = (out @ N) % 3
    return out


# ---------------- (alpha) V(2m) ----------------
def regular_unipotent(m):
    import w33_pass11191_arrow_five_six_qutrits as P91
    return P91.regular_unipotent(m)


def check_alpha(m, other_class=False):
    u = regular_unipotent(m)
    Om = form_std(m)
    if other_class:
        tau = np.diag([1, 2] * m).astype(np.int64)
        u = (tau @ u @ tau) % 3
    N = cayley_inv(u)
    top = mpow(N, 2 * m - 1)
    v = next(np.array(x) for x in itertools.product(range(3), repeat=2 * m) if (top @ np.array(x) % 3).any())
    L = np.stack([(mpow(N, 2 * a) @ v) % 3 for a in range(m)], 1)
    return check_split(u, split_from_lagrangian(u, L, Om), Om)


# ---------------- (beta) W(k) ----------------
def W_matrix(k):
    J = (np.eye(k, dtype=np.int64) + np.eye(k, k=1, dtype=np.int64)) % 3
    blk = np.zeros((2 * k, 2 * k), np.int64)
    blk[:k, :k] = J
    blk[k:, k:] = inv(J).T % 3
    perm = [x for i in range(k) for x in (i, k + i)]
    return blk[np.ix_(perm, perm)] % 3


def check_beta(k):
    u = W_matrix(k)
    Om = form_std(k)
    D = 2 * k
    N = cayley_inv(u)
    e = np.zeros(D, np.int64)
    e[0] = 1
    v = np.zeros(D, np.int64)
    v[2 * (k - 1)] = 1
    ustar = [2 * j + 1 for j in range(k)]
    rows = [((mpow(N, i) @ v) % 3) @ Om % 3 for i in range(k - 1)] + [e @ Om % 3]
    rhs = [1 if i == k - 2 else 0 for i in range(k - 1)] + [0]
    A = np.array(rows)[:, ustar]
    aug = np.concatenate([A, np.array(rhs)[:, None]], 1)
    R, piv = rref(aug)
    w = np.zeros(D, np.int64)
    sol = np.zeros(k, np.int64)
    for i, c in enumerate(piv):
        sol[c] = R[i, -1]
    w[ustar] = sol
    vecs = [e] + [(mpow(N, 2 * a) @ v) % 3 for a in range((k - 1) // 2)] \
        + [(mpow(N, 2 * b) @ w) % 3 for b in range((k - 1) // 2)]
    L = np.stack(vecs, 1)
    return check_split(u, split_from_lagrangian(u, L, Om), Om)


# ---------------- (gamma) unitary Jordan blocks ----------------
class Field:
    """F_{3^r} = F3[y]/(g), elements as coefficient arrays (low degree first)"""

    def __init__(self, r):
        self.r = r
        for coeffs in itertools.product(range(3), repeat=r):
            g = Poly([1] + list(coeffs), _X, modulus=3)
            if g.is_irreducible:
                self.g = [int(c) % 3 for c in g.all_coeffs()]
                break
        self.els = [np.array(c) for c in itertools.product(range(3), repeat=r)]

    def mul(self, a, b):
        r = self.r
        prod = np.zeros(2 * r - 1, np.int64)
        for i in range(r):
            for j in range(r):
                prod[i + j] += a[i] * b[j]
        for p in range(2 * r - 2, r - 1, -1):
            c = prod[p] % 3
            if c:
                prod[p] = 0
                for t in range(1, r + 1):
                    prod[p - t] -= c * self.g[t]
        return prod[:r] % 3

    def pow(self, a, e):
        out = np.zeros(self.r, np.int64)
        out[0] = 1
        for _ in range(e):
            out = self.mul(out, a)
        return out

    def mulmat(self, a):
        cols = []
        for j in range(self.r):
            b = np.zeros(self.r, np.int64)
            b[j] = 1
            cols.append(self.mul(a, b))
        return np.stack(cols, 1) % 3

    def trace(self, a):
        s, x = np.zeros(self.r, np.int64), a.copy()
        for _ in range(self.r):
            s = (s + x) % 3
            x = self.pow(x, 3)
        assert not s[1:].any()
        return int(s[0])


def check_gamma(d, k):
    E = Field(2 * d)
    q = 3 ** d

    def conj(a):
        return E.pow(a, q)
    one = np.zeros(2 * d, np.int64)
    one[0] = 1
    mu = next(a for a in E.els if (E.pow(a, q + 1) == one).all() and not (conj(a) == a).all())
    kappa = next(a for a in E.els if a.any() and not ((conj(a) + a) % 3).any())
    fbasis = []
    for a in E.els:
        if (conj(a) == a).all() and a.any() and rank(np.stack(fbasis + [a], 1)) == len(fbasis) + 1:
            fbasis.append(a)
        if len(fbasis) == d:
            break
    lam = kappa if k % 2 == 1 else one
    D = 2 * d * k
    unit = [np.eye(2 * d, dtype=np.int64)[a] for a in range(2 * d)]
    Om = np.zeros((D, D), np.int64)
    for i in range(k):
        j = k - 1 - i
        h = lam if i % 2 == 0 else (-lam) % 3
        for a in range(2 * d):
            for b in range(2 * d):
                Om[i * 2 * d + a, j * 2 * d + b] = E.trace(E.mul(E.mul(unit[a], conj(unit[b])), h))
    assert not ((Om + Om.T) % 3).any() and rank(Om) == D
    Nm = np.zeros((k, k), np.int64)
    for j in range(k - 1):
        Nm[j + 1, j] = 1
    I = np.eye(k, dtype=np.int64)
    Um = ((I + Nm) @ inv((I - Nm) % 3)) % 3
    S = (np.kron(Um, np.eye(2 * d, dtype=np.int64)) @ np.kron(np.eye(k, dtype=np.int64), E.mulmat(mu))) % 3
    assert not ((S.T @ Om @ S - Om) % 3).any()
    cols = []
    for j in range(k):
        eps = one if j % 2 == 0 else kappa
        for f in fbasis:
            vec = np.zeros(D, np.int64)
            vec[j * 2 * d:(j + 1) * 2 * d] = E.mul(eps, f)
            cols.append(vec)
    L = np.stack(cols, 1) % 3
    return check_split(S, split_from_lagrangian(S, L, Om), Om), D // 2


# ---------------- (delta) f != f* ----------------
def companion(p):
    r = len(p) - 1
    C = np.zeros((r, r), np.int64)
    for i in range(1, r):
        C[i, i - 1] = 1
    for i, c in enumerate(reversed(p[1:])):
        C[i, r - 1] = (-c) % 3
    return C


def check_delta(p, seed=0):
    A = companion(p)
    r = len(A)
    S = np.zeros((2 * r, 2 * r), np.int64)
    S[:r, :r] = A
    S[r:, r:] = inv(A).T % 3
    Om = np.zeros((2 * r, 2 * r), np.int64)
    Om[:r, r:] = np.eye(r, dtype=np.int64)
    Om[r:, :r] = 2 * np.eye(r, dtype=np.int64)
    idx = [(i, j) for i in range(r) for j in range(i, r)]
    mats = []
    for (i, j) in idx:
        B = np.zeros((r, r), np.int64)
        B[i, j] = B[j, i] = 1
        mats.append(B)
    Nsp = nullspace(np.array([((B @ A - A.T @ B) % 3).flatten() for B in mats]).T)
    rng = np.random.default_rng(seed)
    for _ in range(1000):
        v = (rng.integers(0, 3, len(Nsp)) @ Nsp) % 3
        Phi = sum(int(v[t]) * mats[t] for t in range(len(idx))) % 3
        if rank(Phi) == r:
            break
    L = np.concatenate([np.eye(r, dtype=np.int64), Phi], 0)
    return check_split(S, split_from_lagrangian(S, L, Om), Om), r


def polymul(*ps):
    out = Poly(1, _X, modulus=3)
    for p in ps:
        out = out * Poly(p, _X, modulus=3)
    return [int(c) % 3 for c in out.all_coeffs()]


# ---------------- random ticks: exact A at n = 5, 6 ----------------
def invariant_planes(S, n, pts):
    J = form_std(n)
    SX = (pts @ S.T) % 3
    out = []
    for sgn in (1, 2):
        E = pts[((SX - sgn * pts) % 3 == 0).all(1)]
        W = (E @ J @ E.T) % 3
        for a in range(len(E)):
            for b in np.flatnonzero(W[a, a + 1:]) + a + 1:
                out.append(np.stack([E[a], E[b]], 1))
    w = np.einsum('pi,ij,pj->p', pts, J, SX) % 3
    eig = ((SX - pts) % 3 == 0).all(1) | ((SX + pts) % 3 == 0).all(1)
    sel = np.flatnonzero((w != 0) & ~eig)
    U, V = pts[sel], SX[sel]
    S2 = (V @ S.T) % 3
    inv_ = np.zeros(len(U), bool)
    for a in range(3):
        for b in range(3):
            inv_ |= (((a * U + b * V - S2) % 3) == 0).all(1)
    out += [np.stack([U[i], V[i]], 1) for i in np.flatnonzero(inv_)]
    uniq = {}
    for P in out:
        key = frozenset(tuple((a * P[:, 0] + b * P[:, 1]) % 3) for a in range(3) for b in range(3))
        uniq.setdefault(key, P)
    return list(uniq.values())


def max_orthogonal_family(planes, n):
    J = form_std(n)
    m = len(planes)
    adj = np.array([[not ((P.T @ J @ Q) % 3).any() for Q in planes] for P in planes]) if m else np.zeros((0, 0), bool)
    best = [[]]

    def dfs(cands, chosen):
        if len(chosen) > len(best[0]):
            best[0] = chosen
        if len(chosen) + len(cands) <= len(best[0]) or len(best[0]) == n:
            return
        for t, i in enumerate(cands):
            rest = [j for j in cands[t + 1:] if adj[i, j]]
            dfs(rest, chosen + [i])
    dfs(list(range(m)), [])
    return [planes[i] for i in best[0]]


def exact_arrow_random(n, N, seed):
    import w33_pass11165_three_qutrit_arrow as A3
    import w33_pass11191_arrow_five_six_qutrits as P91
    rng = np.random.default_rng(seed)
    pts = P91.points(2 * n)
    J = form_std(n)
    rows = []
    for _ in range(N):
        S, _J = A3.random_symplectic(rng, n=n, steps=40 * n)
        S = S.astype(np.int64) % 3
        fam = max_orthogonal_family(invariant_planes(S, n, pts), n)
        c = len(fam)
        if c == n - 1:
            c = n
        if c < n:
            Pm = np.concatenate(fam, 1) if fam else np.zeros((2 * n, 0), np.int64)
            inW = ~(((pts @ J @ Pm) % 3) != 0).any(1) if fam else np.ones(len(pts), bool)
            xs, _ = _hm_in(S, n, pts[inW], n - len(fam))
            ok = xs is not None
            if ok:
                planes = list(fam) + [np.stack([x, (S @ x) % 3], 1) for x in xs]
                E = check_split(S, planes, J)
            else:
                E = None
        else:
            E = 0
            ok = True
        rows.append(dict(c=c, n_minus_c=n - c, split_found=ok, E=E, exact=ok and E == n - c,
                         jordan_ok=c_jordan(S, n) == c))
    return rows


def _hm_in(S, n, pts, need, restarts=200, budget=50_000, seed=0):
    """need mutually orthogonal planes span(x, Sx) (x among pts, all inside an S-invariant subspace)"""
    J = form_std(n)
    SX = (pts @ S.T) % 3
    w = np.einsum('pi,ij,pj->p', pts, J, SX) % 3
    eig = ((SX - pts) % 3 == 0).all(1) | ((SX + pts) % 3 == 0).all(1)
    good = np.flatnonzero((w != 0) & ~eig)
    X = np.stack([pts[good], SX[good]], 2)
    XJ = np.einsum('mai,ab->mib', X, J)
    rng = np.random.default_rng(seed)

    def orth(i, c):
        if len(c) == 0:
            return c
        W = np.einsum('ib,pbj->pij', XJ[i], X[c]) % 3
        return c[~W.reshape(len(c), 4).any(1)]

    def dfs(c, chosen, bud):
        if len(chosen) == need:
            return chosen
        for t, i in enumerate(c):
            bud[0] -= 1
            if bud[0] < 0 or len(c) - t < need - len(chosen):
                return None
            r = dfs(orth(i, c[t + 1:]), chosen + [int(i)], bud)
            if r:
                return r
        return None
    for _ in range(restarts):
        r = dfs(rng.permutation(len(X)), [], [budget])
        if r:
            return [pts[good[i]] for i in r], None
    return None, None


# ---------------- closed form: c(S) from Jordan data ----------------
def jordan_counts(M, D, smax):
    """number of Jordan blocks of size s = 1..smax of the operator M on its nilpotent part: r_{s-1} - 2 r_s + r_{s+1}"""
    r = [D]
    Pw = np.eye(D, dtype=np.int64)
    for _s in range(1, smax + 2):
        Pw = (Pw @ M) % 3
        r.append(rank(Pw))
    return [r[s - 1] - 2 * r[s] + r[s + 1] for s in range(1, smax + 1)]


def c_jordan(S, n):
    """c(S) = sum over lambda = +-1 of (m_1/2 + m_2) + (number of F9-blocks of size 1 for x^2 + 1)"""
    D = 2 * n
    I = np.eye(D, dtype=np.int64)
    c = 0
    for lam in (1, 2):
        m1, m2 = jordan_counts((S - lam * I) % 3, D, 2)
        c += m1 // 2 + m2
    (m1,) = jordan_counts((S @ S + I) % 3, D, 1)
    return c + m1 // 2


def arrow(S, n):
    """the intrinsic arrow in closed form: A(S) = n - c(S)"""
    return n - c_jordan(S % 3, n)


# ---------------- all classes n <= 4 ----------------
def classes_check():
    import w33_pass11182_paper_ticks_mereology as P
    import w33_pass11188_exact_arrow_by_class as P88
    D = json.loads(P88.OUT.read_text())
    out = {}
    for n in (2, 3, 4):
        cls = P88.load_classes(n)
        Pl = P88.all_planes(n) if n == 4 else P88.planes_from_splits(
            np.array(P88.M.factorisations(2)) if n == 2 else P.load_facs())
        J = P88.M.form(n)
        ok, jok, dist = 0, 0, Counter()
        for c, r in zip(cls, D[f"n{n}"]['rows']):
            d = P88.dvals(c['S'], Pl)
            G = Pl[d == 2].astype(np.int64)
            GJ = np.einsum('pai,ab->pib', G, J)
            best = [0]

            def orth(i, cc):
                if len(cc) == 0:
                    return cc
                W = np.einsum('ib,pbj->pij', GJ[i], G[cc]) % 3
                return cc[~W.reshape(len(cc), 4).any(1)]

            def dfs(cc, k):
                best[0] = max(best[0], k)
                if best[0] == n or k + len(cc) <= best[0]:
                    return
                for t, i in enumerate(cc):
                    if k + len(cc) - t <= best[0]:
                        return
                    dfs(orth(i, cc[t + 1:]), k + 1)
                    if best[0] == n:
                        return
            dfs(np.arange(len(G)), 0)
            ok += r['A'] == n - best[0]
            jok += c_jordan(c['S'], n) == best[0]
            dist[best[0]] += c['size']
        out[f"n{n}"] = dict(classes=len(cls), A_equals_n_minus_c=ok, jordan_formula_ok=jok,
                            c_distribution_psp={str(k): v for k, v in sorted(dist.items())})
    return out


def summarize():
    t0 = time.time()
    res = dict(pass_id=11199)
    res['alpha'] = {f"V({2*m}){'b' if oc else 'a'}": check_alpha(m, oc) for m in range(1, 9) for oc in (False, True)}
    res['beta'] = {f"W({k})": check_beta(k) for k in (3, 5, 7, 9, 11)}
    res['gamma'] = {f"F{3**(2*d)} J{k}": check_gamma(d, k) for d, ks in ((1, range(1, 6)), (2, range(1, 4))) for k in ks}
    f1, f2 = [1, 1, 2], [1, 0, 2, 1]
    res['delta'] = {name: check_delta(p) for name, p in (('x2+x+2', f1), ('(x2+x+2)^2', polymul(f1, f1)),
                                                        ('(x2+x+2)^3', polymul(f1, f1, f1)), ('x3+2x+1', f2),
                                                        ('(x3+2x+1)^2', polymul(f2, f2)))}
    res['constructions_all_optimal'] = (
        all(E == (m if m > 1 else 0) for k, E in res['alpha'].items() for m in [int(k[2:k.index(')')]) // 2])
        and all(E == int(k[2:-1]) for k, E in res['beta'].items())
        and all(E == n or (n == 1 and E == 0) for E, n in res['gamma'].values())
        and all(E == n for E, n in res['delta'].values()))
    res['classes'] = classes_check()
    res['random'] = {}
    for n, N in ((5, 200), (6, 40)):
        rows = exact_arrow_random(n, N, 11199 + n)
        res['random'][f"n{n}"] = dict(sampled=N, all_exact=all(r['exact'] for r in rows),
                                      jordan_formula_ok=all(r['jordan_ok'] for r in rows),
                                      A_distribution=dict(Counter(str(r['n_minus_c']) for r in rows)))
    res['seconds'] = round(time.time() - t0, 1)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
