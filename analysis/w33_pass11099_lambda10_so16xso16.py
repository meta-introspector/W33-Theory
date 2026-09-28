#!/usr/bin/env python3
"""The one-loop modular integral of the ten-dimensional SO(16)xSO(16) heterotic string (numerical, over the fundamental
domain), used by Pass 11099 to fix the SIGN of the vacuum energy that dominates every T6 orbifold of this string at
large volume.

Partition function (Alvarez-Gaume, Ginsparg, Moore, Vafa 1986; in SO(2n) level-1 characters):
    Z = tau2^-4 (eta etabar)^-8 [ Obar8 (V16 C16 + C16 V16) + Vbar8 (O16 O16 + S16 S16)
                                   - Sbar8 (V16 V16 + C16 C16) - Cbar8 (O16 S16 + S16 O16) ]
    O_2n = (th3^n + th4^n)/(2 eta^n), V_2n = (th3^n - th4^n)/(2 eta^n), S_2n = C_2n = th2^n/(2 eta^n).
I = int_F d^2tau / tau2^2 Z ;  Lambda_10 = -(1/2) M^10 I   (M = M_s / 2 pi, the standard normalisation).
Checks built in: modular invariance of the integrand (tau -> -1/tau, tau -> tau + 1) at random points; the level-matched
q^0 qbar^0 coefficient equals n_B - n_F of the massless spectrum: bosons 8 x (8 + 240) = 1984 (graviton sector and the
240 gauge bosons), fermions 8 x 256 + 8 x 256 = 4096 ((16,16) and (128,1)+(1,128)): n_B - n_F = -2112.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "w33_pass11099_lambda10_so16xso16.json"
NT = 40


def thetas(tau):
    q = np.exp(2j * np.pi * tau)
    n = np.arange(-NT, NT + 1)
    th3 = np.sum(q ** (n * n / 2.0))
    th4 = np.sum((-1.0) ** n * q ** (n * n / 2.0))
    th2 = np.sum(q ** ((n + 0.5) ** 2 / 2.0))
    k = np.arange(1, 200)
    eta = q ** (1 / 24) * np.prod(1 - q ** k)
    return th2, th3, th4, eta


def chars(tau, n):
    th2, th3, th4, eta = thetas(tau)
    O = (th3 ** n + th4 ** n) / (2 * eta ** n)
    V = (th3 ** n - th4 ** n) / (2 * eta ** n)
    S = th2 ** n / (2 * eta ** n)
    return O, V, S, S, eta


def Z(tau):
    tau2 = tau.imag
    O16, V16, S16, C16, eta = chars(tau, 8)
    O8, V8, S8, C8, _ = chars(tau, 4)
    O8b, V8b, S8b, C8b = np.conj(O8), np.conj(V8), np.conj(S8), np.conj(C8)
    val = O8b * (V16 * C16 + C16 * V16) + V8b * (O16 * O16 + S16 * S16) - S8b * (V16 * V16 + C16 * C16) \
        - C8b * (O16 * S16 + S16 * O16)
    return (tau2 ** -4) * val / (eta * np.conj(eta)) ** 8


def integrand(t1, t2):
    return Z(complex(t1, t2)).real / t2 ** 2


# ---- exact q-expansions (exponents in units of 1/2): series as dicts {2*exponent: coefficient} -------------------
NMAX = 24


def smul(a, b):
    out = {}
    for i, x in a.items():
        for j, y in b.items():
            if i + j <= NMAX:
                out[i + j] = out.get(i + j, 0) + x * y
    return out


def spow(a, k):
    out = {0: 1}
    for _ in range(k):
        out = smul(out, a)
    return out


def sadd(*terms):
    out = {}
    for c, s in terms:
        for i, x in s.items():
            out[i] = out.get(i, 0) + c * x
    return out


def theta_series():
    th3 = {0: 1}
    th4 = {0: 1}
    th2q = {}                                      # th2 = 2 q^(1/8) * sum q^(n(n+1)/2); th2^(8k) carries q^k
    n = 1
    while n * n <= NMAX:
        th3[n * n] = th3.get(n * n, 0) + 2
        th4[n * n] = th4.get(n * n, 0) + 2 * (-1) ** n
        n += 1
    n = 0
    while n * (n + 1) <= NMAX:
        th2q[n * (n + 1)] = th2q.get(n * (n + 1), 0) + 2
        n += 1
    return th2q, th3, th4


def inv_eta_pow(k):
    """prod_n (1 - q^n)^(-k) as a series in half-units"""
    out = {0: 1}
    for n in range(1, NMAX // 2 + 1):
        # (1 - q^n)^(-k) = sum_j binom(k + j - 1, j) q^(n j)
        from math import comb
        fac = {2 * n * j: comb(k + j - 1, j) for j in range(0, NMAX // (2 * n) + 1)}
        out = smul(out, fac)
    return out


def exact_level_matched():
    """c(e): coefficient of q^e qbar^e in (eta etabar)^-8 [...] including the eta^-n of the characters.
    Left: q^-1 prod(1-q^n)^-24 x {O16, V16, S16}-numerators; right: qbar^-1/2 prod^-12 x {O8,V8,S8}-numerators."""
    th2q, th3, th4 = theta_series()

    def chars_num(n):
        # numerators without eta: O = (th3^n + th4^n)/2, V = (th3^n - th4^n)/2, S = C = th2^n/2 (with q^(n/8) factored)
        t3, t4 = spow(th3, n), spow(th4, n)
        O = {i: x / 2 for i, x in sadd((1, t3), (1, t4)).items()}
        V = {i: x / 2 for i, x in sadd((1, t3), (-1, t4)).items()}
        S = {i: x / 2 for i, x in spow(th2q, n).items()}          # times q^(n/8) = shift of n/4 half-units
        return O, V, S
    O16, V16, S16 = chars_num(8)                   # S16 carries q^(8/8) = q^1 -> shift +2 half-units
    O8, V8, S8 = chars_num(4)                      # S8 carries q^(4/8) = q^(1/2) -> shift +1 half-unit
    sh = lambda s, k: {i + k: x for i, x in s.items() if i + k <= NMAX}
    S16s, S8s = sh(S16, 2), sh(S8, 1)
    # left: q^-1 * prod^-24 ; right: qbar^-1/2 * prod^-12
    L = inv_eta_pow(24)
    Rr = inv_eta_pow(12)
    left = {'VC': smul(L, sadd((2, smul(V16, S16s)))), 'OO+SS': smul(L, sadd((1, smul(O16, O16)), (1, smul(S16s, S16s)))),
            'VV+CC': smul(L, sadd((1, smul(V16, V16)), (1, smul(S16s, S16s)))), 'OS+SO': smul(L, sadd((2, smul(O16, S16s))))}
    right = {'O': smul(Rr, O8), 'V': smul(Rr, V8), 'S': smul(Rr, S8s), 'C': smul(Rr, S8s)}
    pairs = [('O', 'VC', 1), ('V', 'OO+SS', 1), ('S', 'VV+CC', -1), ('C', 'OS+SO', -1)]
    c = {}
    for r, l, sgn in pairs:
        for i, x in right[r].items():            # right exponent (half-units) = i - 1 (qbar^-1/2)
            eR = i - 1
            eL_index = eR + 2                     # left exponent = index - 2 (q^-1)
            y = left[l].get(eL_index, 0)
            if y:
                c[eR] = c.get(eR, 0) + sgn * x * y
    return {e: v for e, v in c.items() if v}


def integrate(n1=240, n2=120):
    """F = {tau2 > 1} (exact: only level-matched terms survive the tau1 integral) + {|tau| >= 1, tau2 < 1} (numerical)"""
    from scipy.integrate import quad
    c = exact_level_matched()
    upper = 0.0
    for e2, v in c.items():
        e = e2 / 2
        assert e >= 0, ("level-matched tachyon", e, v)
        upper += v * quad(lambda t2: t2 ** -6 * np.exp(-4 * np.pi * e * t2), 1, np.inf)[0]
    x1, w1 = np.polynomial.legendre.leggauss(n1)
    x2, w2 = np.polynomial.legendre.leggauss(n2)
    lower = 0.0
    for u, a in zip(x1, w1):
        t1 = 0.5 * u
        lo = np.sqrt(1 - t1 * t1)
        t2s = lo + (1 - lo) * (x2 + 1) / 2
        ww = w2 * (1 - lo) / 2
        lower += 0.5 * a * sum(b * integrand(t1, t2) for t2, b in zip(t2s, ww))
    return upper, lower, c


def main():
    rng = np.random.default_rng(0)
    inv = []
    for _ in range(5):
        tau = complex(rng.uniform(-0.5, 0.5), rng.uniform(0.9, 1.6))
        a = Z(tau).real * tau.imag ** 4 * 0 + (Z(tau) / 1).real
        b = Z(-1 / tau)
        c = Z(tau + 1)
        # Z (including tau2^-4) is modular invariant; d^2tau/tau2^2 is invariant
        inv.append(dict(tau=str(tau), Z=a, Z_S=b.real, Z_T=c.real, rel_S=abs(b.real - a) / abs(a), rel_T=abs(c.real - a) / abs(a)))
    upper, lower, c = integrate()
    I = upper + lower
    res = dict(modular_invariance_checks=inv, massless_constant_term_c00=c.get(0, 0),
               first_level_matched_coefficients={str(k / 2): v for k, v in sorted(c.items())[:6]},
               I_upper_tau2_gt_1=upper, I_lower=lower, I=I, Lambda10_over_M10=-0.5 * I,
               sign="positive" if -0.5 * I > 0 else "negative")
    OUT.write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
