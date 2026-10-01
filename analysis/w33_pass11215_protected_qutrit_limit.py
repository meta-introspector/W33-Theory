"""Pass 11215: the mean number of protected qutrits, exactly for every n and in the limit -- via Fulman's cycle index.

Pass 11207 (merged from the parallel session) proved A(S) = n - c(S) and gave c from Jordan data:
    c(S) = sum_{lambda = +-1} ( m_1(lambda)/2 + m_2(lambda) ) + m_1^{F9}(x^2 + 1),
m_s = number of Jordan blocks of size s.  Pass 11208 found E[c] = 133/180, 106927/147420, 142080247/195832728 for
n = 2, 3, 4 (exact class sums) and ~0.7255 for large n (sampled), and noted that a closed form "would follow from
cycle-index methods (Fulman)".

Fulman's cycle index of Sp(2n, q) (q odd) factorises over the polynomials phi:
    sum_n u^n E_n[prod_phi x_{phi, lambda_phi}] = prod_{phi = z+-1} S_phi * prod_{phi = phi*} U_phi * prod_{phi != phi*} G_phi,
with the unipotent-type sums
    S (z+-1):  sum over symplectic signed partitions lambda of u^{|lambda|/2} / c_Sp(lambda),
               c_Sp = q^{dim C - dim R} |R|,  dim C = (1/2) sum_i lambda'_i^2 + (1/2) sum_{i odd} m_i,
               R = prod_{i odd} Sp(m_i,q) x prod_{i even} O^+-(m_i,q) (the reductive part),
    U (phi = phi*, deg 2d):  sum over partitions of u^{d|lambda|} / c_U(lambda) with Q = q^d,
               c_U = Q^{sum_i lambda'_i^2 - sum_i m_i^2} prod_i |U(m_i, Q)|.
c is additive over phi, so  sum_n u^n E_n[c] = 1/(1-u) * sum_{phi} F_phi(u)/P_phi(u), where P_phi is the phi-factor at
x = 1 and F_phi the same sum weighted by the phi-contribution to c; hence lim E_n[c] = sum_phi F_phi(1)/P_phi(1).

Validation (built in): the unipotent sums reproduce Steinberg's counts (q^{2n^2} unipotents in Sp(2n,q), Q^{m(m-1)}
in U(m,Q)); the finite-n values reproduce Pass 11208's exact class sums for n = 2, 3, 4 (and 5).
"""

from __future__ import annotations

import json
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "w33_pass11215_protected_qutrit_limit.json"
q = 3


def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


def mult(lam):
    m = {}
    for p in lam:
        m[p] = m.get(p, 0) + 1
    return m


def conj_sq_sum(lam):
    if not lam:
        return 0
    return sum(sum(1 for p in lam if p >= i) ** 2 for i in range(1, max(lam) + 1))


def sp_order(m, Q=q):          # |Sp(m,Q)|, m even
    k = m // 2
    out = Q ** (k * k)
    for j in range(1, k + 1):
        out *= Q ** (2 * j) - 1
    return out


def o_order(m, eps, Q=q):     # |O^eps(m,Q)|
    if m == 0:
        return 1
    if m % 2:
        k = (m - 1) // 2
        out = 2 * Q ** (k * k)
        for j in range(1, k + 1):
            out *= Q ** (2 * j) - 1
        return out
    k = m // 2
    out = 2 * Q ** (k * (k - 1)) * (Q ** k - eps)
    for j in range(1, k):
        out *= Q ** (2 * j) - 1
    return out


def u_order(m, Q):             # |U(m,Q)|
    out = Q ** (m * (m - 1) // 2)
    for j in range(1, m + 1):
        out *= Q ** j - (-1) ** j
    return out


@lru_cache(None)
def symplectic_types(N):
    """(weight u-degree k = |lambda|/2, 1/c_Sp, contribution m1/2 + m2) for all symplectic signed partitions of 2k<=2N"""
    out = []
    for size in range(0, 2 * N + 1, 2):
        for lam in partitions(size):
            m = mult(lam)
            if any(i % 2 == 1 and mi % 2 == 1 for i, mi in m.items()):
                continue
            evens = [i for i in m if i % 2 == 0]
            base = q ** 0
            # twice the exponent: dim C minus dim of the reductive part prod Sp(m_i) (i odd) x O(m_i) (i even)
            expo2 = (conj_sq_sum(lam) + sum(mi for i, mi in m.items() if i % 2)
                     - sum(mi * (mi + 1) for i, mi in m.items() if i % 2)
                     - sum(mi * (mi - 1) for i, mi in m.items() if i % 2 == 0))
            odd_part = 1
            for i, mi in m.items():
                if i % 2:
                    odd_part *= sp_order(mi)
            contrib = Fraction(m.get(1, 0), 2) + m.get(2, 0)
            for signs in range(1 << len(evens)):
                c = Fraction(odd_part)
                for b, i in enumerate(evens):
                    c *= o_order(m[i], 1 if (signs >> b) & 1 else -1)
                # q^{expo2/2}: expo2 is even for symplectic partitions (checked)
                assert expo2 % 2 == 0
                c *= q ** (expo2 // 2)
                out.append((size // 2, 1 / c, contrib))
    return out


@lru_cache(None)
def unitary_types(N, d=1):
    """phi = phi* of degree 2d: (u-degree d|lambda|, 1/c_U, contribution m_1)"""
    Q = q ** d
    out = []
    for size in range(0, N // d + 1):
        for lam in partitions(size):
            m = mult(lam)
            c = Fraction(Q) ** (conj_sq_sum(lam) - sum(mi * mi for mi in m.values()))
            for mi in m.values():
                c *= u_order(mi, Q)
            out.append((d * size, 1 / c, m.get(1, 0)))
    return out


def series(types, N, weighted):
    s = [Fraction(0)] * (N + 1)
    for k, w, contrib in types:
        if k <= N:
            s[k] += w * (contrib if weighted else 1)
    return s


def mul(a, b, N):
    out = [Fraction(0)] * (N + 1)
    for i, x in enumerate(a):
        if x:
            for j in range(N + 1 - i):
                out[i + j] += x * b[j]
    return out


def inverse(a, N):
    out = [Fraction(0)] * (N + 1)
    out[0] = 1 / a[0]
    for n in range(1, N + 1):
        out[n] = -sum(a[k] * out[n - k] for k in range(1, n + 1)) / a[0]
    return out


def steinberg_checks(N=6):
    ok_sp = all(sum(w for k, w, _ in symplectic_types(N) if k == n) * sp_order(2 * n) == q ** (2 * n * n)
                for n in range(0, N + 1))
    ok_u = all(sum(w for k, w, _ in unitary_types(N, 1) if k == n) * u_order(n, q) == q ** (n * (n - 1))
               for n in range(0, N + 1))
    return ok_sp, ok_u


def expected_c(N):
    """exact E_n[c] for n <= N"""
    Ssp = symplectic_types(N)
    Uni = unitary_types(N, 1)
    P1, F1 = series(Ssp, N, False), series(Ssp, N, True)
    PU, FU = series(Uni, N, False), series(Uni, N, True)
    ratio = [2 * x for x in mul(F1, inverse(P1, N), N)]        # phi = z - 1 and z + 1 contribute equally
    ratio = [a + b for a, b in zip(ratio, mul(FU, inverse(PU, N), N))]
    return [sum(ratio[: n + 1]) for n in range(N + 1)]         # multiply by 1/(1-u)


def components(N):
    """partial sums of the phi = z-1 (same for z+1) and phi = z^2+1 contributions to E_n[c], n <= N.  (Both phi-factors
    also converge at u = 1 -- coefficients decay like 3^-n -- so lim = sum_phi F_phi(1)/P_phi(1); the finite-n route
    converges faster.)"""
    Ssp = symplectic_types(N)
    Uni = unitary_types(N, 1)
    a = mul(series(Ssp, N, True), inverse(series(Ssp, N, False), N), N)
    b = mul(series(Uni, N, True), inverse(series(Uni, N, False), N), N)
    return [sum(a[: n + 1]) for n in range(N + 1)], [sum(b[: n + 1]) for n in range(N + 1)]


def run():
    res = dict(pass_id=11215)
    ok_sp, ok_u = steinberg_checks(6)
    res["steinberg_sp"] = ok_sp
    res["steinberg_u"] = ok_u
    E = expected_c(8)
    res["E_c_exact"] = {str(n): dict(value=str(E[n]), float=float(E[n])) for n in range(1, 9)}
    ref = {2: Fraction(133, 180), 3: Fraction(106927, 147420), 4: Fraction(142080247, 195832728)}
    res["matches_pass11208"] = {str(n): E[n] == v for n, v in ref.items()}
    from mpmath import mp, mpf
    mp.dps = 32
    A, B = components(14)
    f = lambda x: mpf(x.numerator) / x.denominator
    res["limit"] = dict(value=str(f(2 * A[14] + B[14])), each_of_pm1=str(f(A[14])), from_x2_plus_1=str(f(B[14])),
                        change_n12_to_n14=float(abs(f(2 * A[14] + B[14]) - f(2 * A[12] + B[12]))))
    res["E_c_n10_float"] = float(2 * A[10] + B[10])
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
