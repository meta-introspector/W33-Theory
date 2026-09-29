#!/usr/bin/env python3
"""Pass 11159: an explicit state achieving 4/sqrt15, and an analytic derivation of the value.

In the eigenbases of its one-party reduced states, the numerical optimum of N_AB + N_AC (Passes 11148, 11154) is
W-like: weight 7/15 on |000>, and every excitation of A accompanied by an excitation of exactly one of B, C.  The family
    |psi(x)> = sqrt(x)|000> + sqrt((1-x)/4) sum_{a=1,2} |a> (|0 a> + |a 0>),     x in [0, 1],
has (exact sympy) partial-transpose spectrum for rho_AB
    x;  (1 - x +- sqrt((1 + 15x)(1 - x)))/8  (each twice);  (1 - x)/4 (three times);  -(1 - x)/4,
so its negativity is
    N_AB(x) = N_AC(x) = sqrt((1 + 15x)(1 - x)) / 4.
The maximum is at d/dx[(1 + 15x)(1 - x)] = 14 - 30x = 0, i.e. x = 7/15, where (1 + 15x)(1 - x) = 64/15 and
N_AB = N_AC = 2/sqrt15:  N_AB + N_AC = 4/sqrt15 -- derived, not only identified.  The explicit state
    |psi*> = sqrt(7/15)|000> + sqrt(2/15) sum_{a=1,2} |a>(|0a> + |a0>)
reproduces every invariant of the numerical optimum (spectra (7,4,4)/15, (11,2,2)/15; N_BC = 1/15) to machine precision.
Status: 4/sqrt15 is the maximum of this family exactly and coincides with the numerical global optimum of 46 random
restarts; a proof that no state outside the family does better remains open.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "w33_pass11159_explicit_polygamy_state.json"


def neg(rho):
    pt = rho.reshape(3, 3, 3, 3).transpose(0, 3, 2, 1).reshape(9, 9)
    return float((np.abs(np.linalg.eigvalsh((pt + pt.conj().T) / 2)).sum() - 1) / 2)


def pair(v, keep):
    T = v.reshape(3, 3, 3)
    M = {(0, 1): T.reshape(9, 3), (0, 2): T.transpose(0, 2, 1).reshape(9, 3), (1, 2): T.transpose(1, 2, 0).reshape(9, 3)}[keep]
    return M @ M.conj().T


def state(x):
    v = np.zeros(27)
    v[0] = np.sqrt(x)
    b = np.sqrt((1 - x) / 4)
    for a in (1, 2):
        v[a * 9 + a] = b
        v[a * 9 + a * 3] = b
    return v


def symbolic_spectrum():
    x = sp.symbols('x', positive=True)
    a, b = sp.sqrt(x), sp.sqrt((1 - x) / 4)
    V = [0] * 27
    V[0] = a
    for k in (1, 2):
        V[k * 9 + k] = b
        V[k * 9 + k * 3] = b
    T = lambda i, j, l: V[i * 9 + j * 3 + l]
    pt = sp.zeros(9, 9)
    for i in range(3):
        for j in range(3):
            for ip in range(3):
                for jp in range(3):
                    pt[i * 3 + j, ip * 3 + jp] = sum(T(i, jp, l) * T(ip, j, l) for l in range(3))
    ev = pt.eigenvals()
    negsum = sum(-m * e for e, m in ev.items() if e.subs(x, sp.Rational(7, 15)) < 0)
    nx = sp.simplify(negsum)
    return x, ev, nx


def summarize():
    x, ev, nx = symbolic_spectrum()
    closed = sp.sqrt((1 + 15 * x) * (1 - x)) / 4
    same = all(abs(float((nx - closed).subs(x, t))) < 1e-12 for t in (0.3, 0.5, 7 / 15, 0.7))
    crit = sp.solve(sp.diff((1 + 15 * x) * (1 - x), x), x)
    v = state(7 / 15)
    T = v.reshape(3, 3, 3)
    sA = sorted(np.linalg.svd(T.reshape(3, 9), compute_uv=False) ** 2 * 15)
    sB = sorted(np.linalg.svd(T.transpose(1, 0, 2).reshape(3, 9), compute_uv=False) ** 2 * 15)
    res = dict(pass_id=11159, negativity_closed_form='sqrt((1+15x)(1-x))/4', closed_form_matches_symbolic=same,
               critical_point=[str(c) for c in crit], max_value=str(sp.simplify(2 * closed.subs(x, sp.Rational(7, 15)))),
               explicit_state_NAB=neg(pair(v, (0, 1))), explicit_state_NAC=neg(pair(v, (0, 2))), explicit_state_NBC=neg(pair(v, (1, 2))),
               spectra_times_15={'A': [round(s, 9) for s in sA], 'B': [round(s, 9) for s in sB]})
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
