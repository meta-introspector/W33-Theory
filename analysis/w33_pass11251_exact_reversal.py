"""Pass 11251: the exact best time reversal of (I(x)T) SUM (I(x)T^2) -- F_T = F_min^2, proved.

Pass 11238 found that the cheapest single-sector T-violator U = (I(x)T) SUM (I(x)T^2) has best time-reversal fidelity
0.7123860142010856, and Pass 11237 matched it to F_min^2 = ((1 + 2cos 2pi/9)/3)^2 only to double precision.

Proof strategy (computer-assisted, exact):
  * EXACT ARITHMETIC.  Every qutrit Clifford, after removing a global phase, has entries in 3^{-r/2} mu_18, and the cubic
    gate T = diag(zeta^{x^3}) has entries in mu_9 (zeta = exp(2 pi i/9)).  So V U^* V^dag U has entries in
    3^{-r} Z[zeta] and the trace t lies in Q(zeta).  The global phase of V cancels in V U^* V^dag.  All arithmetic below
    is in Z[zeta] = Z[x]/(x^6 + x^3 + 1) with exact integers.
  * ATTAINED VALUE.  The float search over all 51840 x 81 anti-unitary two-qutrit Cliffords locates the maximisers;
    for one maximiser the trace is recomputed exactly and |t|^2/81 is compared with F_min^4 as elements of Q(zeta).
  * NO LARGER VALUE (Galois gap).  For every V, Y = 81 |t|^2 is an algebraic integer of the cubic field Q(zeta)^+ (r <= 2
    for two qutrits, so 9t is in Z[zeta]).  Every Galois conjugate sigma(V) is again unitary up to 3^{r/2} (sigma
    commutes with complex conjugation in this abelian field), so all conjugates of Y lie in [0, 81^2].  Two different
    values Y1 != Y2 therefore differ by at least 1/(2 * 81^2)^2 (their difference has integer norm >= 1 and the other
    two conjugates are bounded by 2 * 81^2).  The float evaluation errs by far less than the resulting gap in F, so the
    float maximum identifies the exact maximum.
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11227_two_qutrit_cp_mixing as P2  # noqa: E402

OUT = ROOT / "data" / "w33_pass11251_exact_reversal.json"
ZETA = np.exp(2j * np.pi / 9)


# ---------------------------------------------------------------- Z[zeta9] arithmetic (coefficient vectors, length 6)
def reduce(c):
    """reduce a coefficient vector of any length modulo x^6 + x^3 + 1"""
    c = list(c)
    for k in range(len(c) - 1, 5, -1):           # x^k = -x^{k-3} - x^{k-6}
        a = c[k]
        if a:
            c[k - 3] -= a
            c[k - 6] -= a
        c[k] = 0
    return tuple(c[:6]) + (0,) * (6 - min(6, len(c)))


ZERO = (0,) * 6
ONE = (1, 0, 0, 0, 0, 0)


def zpow(k):
    """zeta9^k"""
    c = [0] * 9
    c[k % 9] = 1
    return reduce(c)


def z18pow(k):
    """zeta18^k = (-1)^k zeta9^{5k}  (zeta18 = -zeta9^5)"""
    return tuple((-1) ** (k % 2) * v for v in zpow(5 * k))


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def mul(a, b):
    c = [0] * 11
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    c[i + j] += x * y
    return reduce(c)


def conj(a):
    """complex conjugation zeta -> zeta^{-1} = zeta^8"""
    out = ZERO
    for i, x in enumerate(a):
        if x:
            out = add(out, tuple(x * v for v in zpow(-i)))
    return out


def galois(a, k):
    """the automorphism zeta -> zeta^k (k coprime to 9)"""
    out = ZERO
    for i, x in enumerate(a):
        if x:
            out = add(out, tuple(x * v for v in zpow(i * k)))
    return out


def to_complex(a):
    return complex(sum(x * ZETA ** i for i, x in enumerate(a)))


def matmul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[_dot([A[i][j] for j in range(m)], [B[j][k] for j in range(m)]) for k in range(p)] for i in range(n)]


def _dot(u, v):
    s = ZERO
    for x, y in zip(u, v):
        if x != ZERO and y != ZERO:
            s = add(s, mul(x, y))
    return s


def mconj(A):
    return [[conj(x) for x in row] for row in A]


def dagger(A):
    return [[conj(A[j][i]) for j in range(len(A))] for i in range(len(A[0]))]


def trace(A):
    s = ZERO
    for i in range(len(A)):
        s = add(s, A[i][i])
    return s


def exactify(V, tol=1e-9):
    """V (float, unitary, all nonzero entries of equal modulus 3^{-r/2} up to a global phase) -> (W, r) with
    W over Z[zeta] (entries in {0} u mu_18) and V = v0 * W, |v0|^2 = 3^{-r}"""
    flat = V.ravel()
    k0 = int(np.argmax(np.abs(flat) > tol))
    v0 = flat[k0]
    r = round(-2 * np.log(abs(v0)) / np.log(3))
    assert abs(abs(v0) ** 2 - 3.0 ** -r) < 1e-9
    W = []
    for row in V / v0:
        out = []
        for z in row:
            if abs(z) < tol:
                out.append(ZERO)
                continue
            assert abs(abs(z) - 1) < 1e-9, "entries of unequal modulus"
            k = np.angle(z) / (2 * np.pi / 18)
            kr = round(k)
            assert abs(k - kr) < 1e-7, "phase not an 18th root of unity"
            out.append(z18pow(kr))
        W.append(out)
    return W, r


# ---------------------------------------------------------------- Q(zeta) helpers with Fractions
def qmul(a, b):
    c = [Fraction(0)] * 11
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    c = list(c)
    for k in range(10, 5, -1):
        a_ = c[k]
        c[k - 3] -= a_
        c[k - 6] -= a_
        c[k] = 0
    return tuple(c[:6])


def minpoly_real(a):
    """minimal polynomial over Q of a real element a of Q(zeta)^+ (degree 1 or 3), via its conjugates under
    zeta -> zeta^k, k in {1, 2, 4}"""
    conjs = [tuple(Fraction(x) for x in galois(tuple(int(v) for v in a), k)) if all(Fraction(v).denominator == 1 for v in a)
             else None for k in (1, 2, 4)]
    if conjs[0] is None:
        raise ValueError("pass integral elements")
    one = tuple(Fraction(v) for v in ONE)
    e1 = tuple(sum(c[i] for c in conjs) for i in range(6))
    e2 = tuple(qmul(conjs[0], conjs[1])[i] + qmul(conjs[0], conjs[2])[i] + qmul(conjs[1], conjs[2])[i] for i in range(6))
    e3 = qmul(qmul(conjs[0], conjs[1]), conjs[2])
    for e in (e1, e2, e3):
        assert all(v == 0 for v in e[1:]), "symmetric function not rational"
    del one
    return [1, -e1[0], e2[0], -e3[0]]


def float_maximisers(U, reps, tol=1e-9):
    Pa = P2.PA
    Y = np.einsum('aji,jk,akl->ail', Pa.conj(), U, Pa)
    vals = np.empty((len(reps), len(Pa)))
    for s in range(0, len(reps), 4000):
        C = reps[s:s + 4000]
        B = np.einsum('rij,jk,rlk->ril', C, U.conj(), C.conj())
        vals[s:s + 4000] = np.abs(np.einsum('rjk,akj->ra', B, Y)) / 9
    return vals


def run():
    reps = np.load(P2.CACHE)
    T, I3 = P2.T, P2.I3
    U = np.kron(I3, T) @ P2.SUM @ np.kron(I3, T @ T)
    vals = float_maximisers(U, reps)
    A = np.abs(reps)
    min_nonzero_modulus = float(A[A > 1e-9].min())
    assert min_nonzero_modulus > 1 / 3 - 1e-9          # every V has r <= 2, so 9 t_true lies in Z[zeta]
    fmax = float(vals.max())
    levels = np.unique(np.round(vals.ravel(), 9))
    arg = np.argwhere(vals > fmax - 1e-9)
    r_i, a_i = arg[0]
    V = P2.PA[a_i] @ reps[r_i]
    W, r = exactify(V)
    Ue, rU = exactify(U)
    assert rU == 0
    M = matmul(matmul(matmul(W, mconj(Ue)), dagger(W)), Ue)
    t = trace(M)                                     # true trace = 3^{-r} t (times |v0|^2 = 3^{-r})
    # |t_true|^2 / 81 = t conj(t) / (3^{2r} 81)
    tt = mul(t, conj(t))
    # F_min = (1 + zeta + zeta^8)/3, so F_min^4 = (1 + zeta + zeta^8)^4 / 81
    g = add(add(ONE, zpow(1)), zpow(8))
    g4 = mul(mul(g, g), mul(g, g))
    lhs = tt                                         # = 3^{2r} * 81 * F_T^2 ... compare F_T^2 with F_min^4
    rhs = tuple(v * 3 ** (2 * r) for v in g4)        # F_T^2 = tt / (3^{2r} 81), F_min^4 = g4 / 81
    equal = (lhs == rhs)
    # F_T itself = F_min^2 (both positive): exact
    Y = tt                                           # 81 |t_true|^2 * 3^{2r} / 81 ... integral element
    mp = minpoly_real(Y)
    # Galois gap: Y_true = 81 |t_true|^2 integral with conjugates in [0, 81^2]
    gapY = 1.0 / (2 * 81 ** 2) ** 2
    F2 = (1 + 2 * np.cos(2 * np.pi / 9)) ** 2 / 9
    gapF = gapY / (2 * 81 * 81 * F2)                 # dF >= dY / (2 * 81^2 * F) at F ~ F_min^2 (Y = 81^2 F^2)
    above = sorted(x for x in levels if x > fmax - 1e-3)
    second = float(sorted(levels)[-2]) if len(levels) > 1 else None
    # structure of the maximisers: how many, and are they products?
    prods = 0
    for r2, a2 in arg:
        Vm = P2.PA[a2] @ reps[r2]
        R = Vm.reshape(3, 3, 3, 3).transpose(0, 2, 1, 3).reshape(9, 9)
        if np.linalg.matrix_rank(R, tol=1e-8) == 1:
            prods += 1
    return dict(
        pass_id=11251,
        tick="(I(x)T) SUM (I(x)T^2)",
        float_max=fmax,
        F_min_squared=F2,
        float_minus_Fmin2=fmax - F2,
        n_maximisers=int(len(arg)),
        n_product_maximisers=prods,
        maximiser_r=r,
        exact_trace_coeffs=list(t),
        exact_t_conj_t=list(tt),
        F_T_squared_equals_F_min_fourth_exactly=bool(equal),
        minpoly_of_t_conj_t=[str(x) for x in mp],
        galois_gap_in_Y=gapY,
        galois_gap_in_F=gapF,
        float_error_bound_used=1e-13,
        min_nonzero_clifford_modulus=min_nonzero_modulus,
        second_level=second,
        proved_F_T_equals_F_min_squared=bool(equal and gapF > 1e-13),
        levels_near_top=[float(x) for x in above],
    )


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
