#!/usr/bin/env python3
"""Pass 11173: the temporal CGLMP value of one qudit tends to 4 -- an explicit strategy reaches
        I_d = 4 - 2 L(d)/(d - 1),   L(d) = 3 + 2(sqrt d - 2)/(sqrt d (sqrt d - 1)),
so I_d >= 4 - 6/(d-1) - O(d^{-3/2}); with Pass 11171 (I_d < 4 for every finite d), sup_d I_d = 4, never attained.

THE RAMP.  For CGLMP with d outcomes every term collapses to one linear ramp: with the cyclic lag
delta_xy in {0,...,d-1} of Bob's outcome from the target relation (delta_00 = a - b, delta_10 = b - a - 1,
delta_11 = a - b, delta_01 = b - a, all mod d),
        I_d = 4 - (2/(d-1)) * L,     L = sum_{x,y} E[delta_xy].
A local model has delta_00 + delta_10 + delta_11 + delta_01 = -1 mod d, so L >= d - 1 and I_d <= 2.  Approaching 4
means keeping the TOTAL EXPECTED LAG BOUNDED as d grows.
THE STRATEGY (read off the numerical optimisers of Pass 11171).  psi = |0> is a 'rest vector':
  A_0 = B_0 = standard basis (psi = A_0 vector 0 = B_0 vector 0);
  B_1 vector 0 = |0> too, and on the complement B_1 is B_0 shifted by one: B_1[a] = |a+1> (a = 1..d-2), B_1[d-1] = |1>;
  A_1 = Q B_1 with Q the Householder reflection exchanging |0> and the uniform vector g = sum_j |j>/sqrt d.
Setting x = 0 leaves the qudit on the rest vector, which neither later measurement disturbs (lags 0).  Setting x = 1
writes a uniformly random outcome into the shift register (B_0 label = B_1 label + 1); the only costs are the leakage of
the rest vector (each A_1 vector has overlap 1/sqrt d with it).  Exactly, with s = sqrt d, the Q entries are 1/s,
1 - 1/(s(s-1)), -1/(s(s-1)), and the lag splits as ~1 (the spread vector A_1[0] = g) + ~1 (landing on the rest vector,
cost d - 1 with probability 1/d^2 per outcome) + ~1 (Householder spreading):  L(d) = 3 + 2(s-2)/(s(s-1)).
Checked exactly (sympy, perfect squares d = 4..64) and in floating point to d = 1024 (I_1024 = 3.99402).
The numerical optima (Pass 11171) do slightly better, L ~ 1.7 (deficit ~3.4/d): same architecture -- p(a|0) concentrated
on one vector common to B_0 and B_1 (0.978 at d = 8), B_1 = B_0 shifted by one on its complement -- with the spread tuned.
Reading: the qudit's only memory between the two times is its post-measurement vector; it can still carry the SETTING
(rest vector or shift register) as well as the outcome, at a cost O(1/d).  That is why time can approach the algebraic
maximum while space (CGLMP -> 2.9696 for d -> infinity) cannot.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "w33_pass11173_temporal_cglmp_to_four.json"
CONV = {(0, 0): lambda a, b, d: (a - b) % d, (1, 0): lambda a, b, d: (b - a - 1) % d,
        (1, 1): lambda a, b, d: (a - b) % d, (0, 1): lambda a, b, d: (b - a) % d}


def construct(d):
    I = np.eye(d)
    perm = [0] + [a + 1 for a in range(1, d - 1)] + [1]
    B1 = I[:, perm]
    g = np.ones(d) / np.sqrt(d)
    v = I[:, 0] - g
    Q = I - 2 * np.outer(v, v) / (v @ v)
    return [I, Q @ B1], [I, B1], I[:, 0]


def lag_and_value(A, B, psi):
    d = len(psi)
    aa, bb = np.meshgrid(np.arange(d), np.arange(d), indexing='ij')
    L = 0.0
    for (x, y), fn in CONV.items():
        p = np.abs(A[x].conj().T @ psi) ** 2
        T = (np.abs(B[y].conj().T @ A[x]) ** 2).T
        L += float((p[:, None] * T * fn(aa, bb, d)).sum())
    return L, 4 - 2 * L / (d - 1)


def cglmp_direct(A, B, psi):
    """independent evaluation with the standard CGLMP coefficients (Pass 11171's coeff)"""
    import sys
    sys.path.insert(0, str(ROOT / "analysis"))
    import w33_pass11171_scan_temporal_cglmp_large_d as S
    d = len(psi)
    C = S.coeff(d)
    tot = 0.0
    for x in range(2):
        p = np.abs(A[x].conj().T @ psi) ** 2
        for y in range(2):
            T = (np.abs(B[y].conj().T @ A[x]) ** 2).T
            tot += float((C[x, y] * p[:, None] * T).sum())
    return tot


def L_formula(d):
    s = np.sqrt(d)
    return 3 + 2 * (s - 2) / (s * (s - 1))


def L_exact_square(s):
    s = sp.Integer(s)
    d = int(s * s)
    Qv = lambda b, j: (1 / s) if (j == 0 or b == 0) else ((1 - 1 / (s * (s - 1))) if b == j else (-1 / (s * (s - 1))))
    perm = lambda a: 0 if a == 0 else (a + 1 if a <= d - 2 else 1)
    tot = sp.Integer(0)
    for a in range(d):
        j = perm(a)
        for b in range(d):
            binv = 0 if b == 0 else (b - 1 if 2 <= b <= d - 1 else d - 1)
            tot += sp.Rational(1, d) * Qv(b, j) ** 2 * (((b - a - 1) % d) + ((a - binv) % d))
    return sp.nsimplify(tot)


def summarize():
    rows = {}
    for d in (3, 4, 6, 8, 10, 16, 24, 64, 256, 1024):
        A, B, psi = construct(d)
        L, v = lag_and_value(A, B, psi)
        rows[str(d)] = dict(L=L, I=v, L_formula=L_formula(d), direct=cglmp_direct(A, B, psi) if d <= 256 else None)
    exact = {str(s * s): str(L_exact_square(s)) for s in (2, 3, 4, 5, 6, 8)}
    exact_ok = all(sp.nsimplify(sp.Integer(3) + sp.Rational(2 * (s - 2), s * (s - 1))) == sp.nsimplify(exact[str(s * s)])
                   for s in (2, 3, 4, 5, 6, 8))
    res = dict(pass_id=11173, family=rows, exact_L_perfect_squares=exact, exact_matches_formula=bool(exact_ok),
               formula_matches_float=all(abs(r['L'] - r['L_formula']) < 1e-9 for r in rows.values()),
               ramp_matches_direct=all(r['direct'] is None or abs(r['direct'] - r['I']) < 1e-9 for r in rows.values()),
               I_1024=rows['1024']['I'])
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    r = summarize()
    print(json.dumps({k: v for k, v in r.items() if k != 'family'}, indent=1))
    for d, row in r['family'].items():
        print(d, round(row['I'], 6), round(row['L'], 6), row['direct'])
