"""Pass 11253: every time-reversal fidelity of qutrit Clifford+T dynamics has F^2 in the cubic field Q(cos 2pi/9).

THEOREM (any number n of qutrits).  Let U be a word in qutrit Cliffords and the cubic gate T = diag(zeta^{x^3}),
zeta = exp(2 pi i/9), and let V be a Clifford.  Up to global phases (which cancel in V U^* V^dag), every Clifford has
entries in 3^{-r/2} mu_18 and T has entries in mu_9 -- all in Q(zeta) up to a real power of sqrt 3.  Hence
t = tr(V U^* V^dag U) lies in Q(zeta) (the sqrt 3 powers pair up into rational powers of 3, since sqrt 3 is real), and
|t|^2 = t conj(t) is fixed by complex conjugation, so it lies in the maximal real subfield Q(zeta)^+ = Q(cos 2pi/9),
which is CUBIC.  The best time-reversal fidelity F_T = max_V |t|/3^n is attained, so F_T^2 is in Q(cos 2pi/9): it is
rational or a cubic irrationality.  (Pass 11237's levels 0.939 and 0.7258 were 'not identified at degree <= 12' only
because their minimal polynomials have coefficients beyond the 10^5 search bound -- denominators are powers of 3.)

COMPUTATION.  The random word streams of Pass 11235 (two qutrits) and Pass 11237 (one qutrit) are replayed with the
same seeds.  For every violating word: the float search finds a maximising reversal; U and V are rebuilt exactly over
Z[zeta] (Clifford factors exactified from their float matrices, T exact), divided by the largest power of sqrt(-3);
t is computed exactly, and F^2 is written as a + b c + e c^2, c = 2cos(2pi/9), with its minimal polynomial over Q.
The float and exact values must agree to 1e-12 for every word.
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11213_cubic_t_violation as P1  # noqa: E402
import w33_pass11227_two_qutrit_cp_mixing as P2  # noqa: E402
import w33_pass11235_two_qutrit_t_spectrum as S35  # noqa: E402
import w33_pass11251_exact_reversal as E  # noqa: E402

OUT = ROOT / "data" / "w33_pass11253_cubic_field_levels.json"
S = E.add(E.zpow(3), tuple(-v for v in E.zpow(6)))          # sqrt(-3) = zeta^3 - zeta^6
C = E.add(E.zpow(1), E.zpow(8))                               # c = 2 cos(2 pi/9)
C2 = E.mul(C, C)


def divide_by_s(a):
    """a / sqrt(-3) if it lies in Z[zeta], else None  (a/s = a * conj(s) / 3 = -a s / 3)"""
    p = E.mul(a, S)
    if all(v % 3 == 0 for v in p):
        return tuple(-v // 3 for v in p)
    return None


def reduce_scale(W, r):
    """W / s^r  ->  divide out common factors of s = sqrt(-3); returns (W', r')"""
    while r > 0:
        D = [[divide_by_s(x) for x in row] for row in W]
        if any(x is None for row in D for x in row):
            break
        W, r = D, r - 1
    return W, r


def exact_product(factors):
    """factors: list of float unitaries applied right-to-left as given (factors[0] is leftmost);
    each is either an equal-modulus Clifford or a diagonal T-type gate with mu_9 entries"""
    W, r = None, 0
    for F in factors:
        Wf, rf = E.exactify(F) if not _is_diag_mu9(F) else (_diag_exact(F), 0)
        W = Wf if W is None else E.matmul(W, Wf)
        r += rf
    return reduce_scale(W, r)


def _is_diag_mu9(F):
    return np.allclose(F, np.diag(np.diag(F))) and np.allclose(np.abs(np.diag(F)), 1) and \
        all(abs(np.angle(z) / (2 * np.pi / 9) - round(np.angle(z) / (2 * np.pi / 9))) < 1e-7 for z in np.diag(F))


def _diag_exact(F):
    n = len(F)
    W = [[E.ZERO] * n for _ in range(n)]
    for i, z in enumerate(np.diag(F)):
        W[i][i] = E.zpow(round(np.angle(z) / (2 * np.pi / 9)))
    return W


def in_basis(y):
    """write a real element y of Q(zeta) (Fraction coefficients) as a + b c + e c^2"""
    import sympy
    a, b, e = sympy.symbols("a b e")
    eqs = []
    for i in range(6):
        eqs.append(a * E.ONE[i] + b * C[i] + e * C2[i] - y[i])
    sol = sympy.solve(eqs, [a, b, e], dict=True)
    assert len(sol) == 1
    return [Fraction(str(sol[0][v])) for v in (a, b, e)]


def minpoly(y):
    """minimal polynomial (monic, rational) of a real element y of Q(zeta), from its conjugates zeta -> zeta^k,
    k = 1, 2, 4"""
    den = 1
    for v in y:
        den = den * Fraction(v).denominator // np.gcd(den, Fraction(v).denominator)
    yi = tuple(int(Fraction(v) * den) for v in y)
    conjs = [E.galois(yi, k) for k in (1, 2, 4)]
    e1 = [sum(cj[i] for cj in conjs) for i in range(6)]
    p01, p02, p12 = E.mul(conjs[0], conjs[1]), E.mul(conjs[0], conjs[2]), E.mul(conjs[1], conjs[2])
    e2 = [p01[i] + p02[i] + p12[i] for i in range(6)]
    e3 = E.mul(p01, conjs[2])
    for e in (e1, e2, e3):
        assert all(v == 0 for v in e[1:])
    coeffs = [Fraction(1), Fraction(-e1[0], den), Fraction(e2[0], den ** 2), Fraction(-e3[0], den ** 3)]
    if coeffs[1] ** 2 == 3 * coeffs[2] and coeffs[1] ** 3 == 27 * coeffs[3]:     # (x - y)^3 = x^3 - 3y x^2 + 3y^2 x - y^3
        return [Fraction(1), coeffs[1] / 3]
    return coeffs


def exact_level(U_factors, V, n_qutrits):
    """exact F^2 for the word U (list of float factors, leftmost first) and the reversal V (float Clifford)"""
    WU, rU = exact_product(U_factors)
    WV, rV = E.exactify(V)
    WV, rV = reduce_scale(WV, rV)
    sgn = (-1) ** rU                                   # conj(W / s^r) = (-1)^r conj(W) / s^r
    M = E.matmul(E.matmul(E.matmul(WV, E.mconj(WU)), E.dagger(WV)), WU)
    t = E.trace(M)
    # V U^* V^dag U = (V-scale) * W..., |V-scale|^2 = 3^{-rV}; U^* U scales: |s|^{-2 rU} = 3^{-rU}
    tt = E.mul(t, E.conj(t))
    d2 = 9 ** n_qutrits
    F2 = tuple(Fraction(v, 3 ** (2 * rV + 2 * rU) * d2) for v in tt)
    del sgn
    return F2, rU, rV


def sqrt_in_field(abc):
    """is F^2 = a + b c + e c^2 (c = 2cos 2pi/9, c^3 = 3c - 1) the square of an element of Q(c)?  Candidates from the
    real conjugates (sigma_k c = 2cos(2 pi k/9), k = 1, 2, 4) with every sign pattern, verified exactly."""
    import itertools as it
    a, b, e = abc
    cs = [2 * np.cos(2 * np.pi * k / 9) for k in (1, 2, 4)]
    vals = [float(a + b * x + e * x * x) for x in cs]
    if min(vals) < 0:
        return None
    Vm = np.array([[1, x, x * x] for x in cs])
    for signs in it.product((1, -1), repeat=3):
        y = np.linalg.solve(Vm, [s * np.sqrt(v) for s, v in zip(signs, vals)])
        cand = [Fraction(float(v)).limit_denominator(3 ** 12) for v in y]
        p, q, r = cand
        sq = (p * p - 2 * q * r, 2 * p * q + 6 * q * r - r * r, q * q + 2 * p * r + 3 * r * r)
        if sq == (a, b, e):
            return [str(x) for x in (cand if cand[0] + cand[1] * cs[0] + cand[2] * cs[0] ** 2 > 0 else [-x for x in cand])]
    return None


def value(y):
    return sum(float(v) * np.real(np.exp(2j * np.pi * i / 9)) for i, v in enumerate(y))


def one_qutrit_stream(per_depth=400, seed=11228):
    """Pass 11237's one-qutrit stream: U = Cf[i0], then U = Cf[j] T U"""
    Cf = P1.clifford1()
    rng = np.random.default_rng(seed)
    for d in range(1, 7):
        for _ in range(per_depth):
            idx = [rng.integers(len(Cf))] + [int(rng.integers(len(Cf))) for _ in range(d)]
            factors = []
            for j in reversed(idx[1:]):
                factors += [Cf[j], P1.T1]
            factors.append(Cf[idx[0]])
            yield d, factors


def two_qutrit_stream(per_depth=60, seed=11235):
    """Pass 11235's stream (same rng calls as S35.random_word)"""
    reps = np.load(P2.CACHE)
    rng = np.random.default_rng(seed)
    T1 = np.kron(P2.T, P2.I3)
    T2 = np.kron(P2.I3, P2.T)
    for d in range(1, 5):
        for _ in range(per_depth):
            fac = [P2.PA[rng.integers(81)] @ reps[rng.integers(len(reps))]]
            for _ in range(d):
                c = P2.PA[rng.integers(81)] @ reps[rng.integers(len(reps))]
                g = T1 if rng.integers(2) == 0 else T2
                fac = [c, g] + fac
            yield d, fac


def prod(factors):
    U = factors[0]
    for F in factors[1:]:
        U = U @ F
    return U


def best_one(U, Cf):
    vals = [abs(np.trace(V @ U.conj() @ V.conj().T @ U)) / 3 for V in Cf]
    k = int(np.argmax(vals))
    return vals[k], Cf[k]


def best_two(U, reps):
    vals = E.float_maximisers(U, reps)
    r, a = np.unravel_index(int(np.argmax(vals)), vals.shape)
    return float(vals[r, a]), P2.PA[a] @ reps[r]


def run():
    res = dict(pass_id=11253)
    levels = {}
    worst = 0.0
    Cf = P1.clifford1()
    n1 = 0
    for d, fac in one_qutrit_stream():
        U = prod(fac)
        f, V = best_one(U, Cf)
        if f > 1 - 1e-9:
            continue
        key = ("1q", round(f, 9))
        if key in levels:
            continue
        F2, rU, rV = exact_level(fac, V, 1)
        worst = max(worst, abs(value(F2) - f * f))
        abc = in_basis(F2)
        levels[key] = dict(qutrits=1, depth=d, F=f, F2_abc=[str(x) for x in abc], F_abc=sqrt_in_field(abc),
                           F2_minpoly=[str(x) for x in minpoly(F2)], scale_rU=rU, scale_rV=rV)
        n1 += 1
    reps = np.load(P2.CACHE)
    n2 = 0
    for d, fac in two_qutrit_stream():
        U = prod(fac)
        f, V = best_two(U, reps)
        if f > 1 - 1e-9:
            continue
        key = ("2q", round(f, 9))
        if key in levels:
            continue
        F2, rU, rV = exact_level(fac, V, 2)
        worst = max(worst, abs(value(F2) - f * f))
        abc = in_basis(F2)
        levels[key] = dict(qutrits=2, depth=d, F=f, F2_abc=[str(x) for x in abc], F_abc=sqrt_in_field(abc),
                           F2_minpoly=[str(x) for x in minpoly(F2)], scale_rU=rU, scale_rV=rV)
        n2 += 1
        print(key, levels[key]["F2_minpoly"], flush=True)
    res["levels"] = sorted(levels.values(), key=lambda x: (x["qutrits"], -x["F"]))
    res["distinct_levels_one_qutrit"] = n1
    res["distinct_levels_two_qutrit"] = n2
    res["max_float_vs_exact_F2_deviation"] = worst
    res["all_F2_in_Q_cos_2pi_9"] = True                 # by construction: in_basis succeeded for every level
    res["cubic_levels"] = sum(1 for v in levels.values() if len(v["F2_minpoly"]) == 4)
    res["rational_levels"] = sum(1 for v in levels.values() if len(v["F2_minpoly"]) == 2)
    res["levels_with_F_itself_in_Q_cos_2pi_9"] = sum(1 for v in levels.values() if v["F_abc"] is not None)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "levels"}, indent=1))
    for lv in res["levels"]:
        print(lv["qutrits"], round(lv["F"], 9), lv["F2_minpoly"], lv["F2_abc"])


if __name__ == "__main__":
    main()
