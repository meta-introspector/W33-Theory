"""Pass 11237: what the time-reversal fidelity levels are -- exact algebraic numbers (120-digit recomputation).

Passes 11228, 11235 and 11238 found a discrete spectrum of best time-reversal fidelities F_T for Clifford+cubic
dynamics: 0.8440296 (every minimal violator), 0.9392626, 0.8100955, 0.7123860, 0.7643761, 0.7257877, ...

Method (sound identification -- a first attempt with findpoly at double precision and coefficients up to 1e7 produced
spurious degree-2 'fits' for every number, because such fits have more free digits than 15-digit data; it is replaced):
  * the 216 single-qutrit Cliffords are rebuilt EXACTLY in mpmath (120 digits) by breadth-first search from H, S, X, Z;
  * random one-qutrit words are regenerated with the same random stream as the float pass and screened in floats;
  * for each target level, a word realising it is recomputed exactly: F_T = max_V |tr(V U^* V^dag U)|/3 over all 216 V;
  * the minimal polynomial of F_T^2 is sought with small coefficients (|c| <= 10^5, degree <= 12) and accepted only if it
    vanishes to 10^-105 at the 120-digit value.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11213_cubic_t_violation as P1  # noqa: E402
import w33_pass11228_t_witness as W1  # noqa: E402

OUT = ROOT / "data" / "w33_pass11237_fidelity_levels.json"
TARGETS = [0.8440296287, 0.9392625771, 0.8100954849, 0.7123860140, 0.7643761380, 0.7257876540]


def mp_setup():
    import mpmath as mp
    mp.mp.dps = 120
    w = mp.exp(2j * mp.pi / 3)
    z9 = mp.exp(2j * mp.pi / 9)
    H = mp.matrix([[w ** (j * k) / mp.sqrt(3) for k in range(3)] for j in range(3)])
    S = mp.diag([1, 1, w])
    X = mp.matrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    Z = mp.diag([1, w, w ** 2])
    T = mp.diag([1, z9, z9 ** 8])
    return mp, H, S, X, Z, T


def exact_cliffords(mp, gens):
    def tofloat(M):
        return np.array([[complex(M[i, j]) for j in range(3)] for i in range(3)])
    start = mp.eye(3)
    seen = {P1.key(tofloat(start)): start}
    frontier = [start]
    while frontier:
        new = []
        for U in frontier:
            for g in gens:
                V = g * U
                k = P1.key(tofloat(V))
                if k not in seen:
                    seen[k] = V
                    new.append(V)
        frontier = new
    return seen


def exact_fidelity(mp, U, Cl):
    best = mp.mpf(0)
    Uc = U.conjugate()
    for V in Cl:
        M = V * Uc * V.H * U
        best = max(best, abs(M[0, 0] + M[1, 1] + M[2, 2]) / 3)
    return best


def identify(mp, x):
    """minimal polynomial with |c| <= 10^5, degree <= 12, verified to 10^-105 at 120 digits (far more digits than
    free parameters, so a hit is not a fit)"""
    for deg in range(1, 13):
        c = mp.findpoly(x, deg, maxcoeff=10 ** 5, tol=mp.mpf(10) ** -100)
        if c:
            val = sum(ci * x ** (len(c) - 1 - i) for i, ci in enumerate(c))
            if abs(val) < mp.mpf(10) ** -105:
                return deg, [int(v) for v in c]
    return None, None


def in_cubic_field(mp, x):
    """integer relation between x and 1, cos(2pi/9), cos(2pi/9)^2"""
    c = mp.cos(2 * mp.pi / 9)
    return mp.pslq([x, 1, c, c ** 2], maxcoeff=10 ** 8, maxsteps=10 ** 6)


def run(per_depth=400, seed=11228):
    mp, H, S, X, Z, T = mp_setup()
    exact = exact_cliffords(mp, [H, S, X, Z])
    assert len(exact) == 216
    Cf = P1.clifford1()
    keys = [P1.key(V) for V in Cf]
    Cl = [exact[k] for k in keys]                       # same order as the float list
    rng = np.random.default_rng(seed)
    found = {}
    for d in range(1, 7):
        for _ in range(per_depth):
            idx = [rng.integers(len(Cf))] + [int(rng.integers(len(Cf))) for _ in range(d)]
            U = Cf[idx[0]]
            for j in idx[1:]:
                U = Cf[j] @ P1.T1 @ U
            f = W1.fidelity(U, Cf)
            for t in TARGETS:
                if abs(f - t) < 1e-8 and t not in found:
                    Ue = Cl[idx[0]]
                    for j in idx[1:]:
                        Ue = Cl[j] * T * Ue
                    found[t] = Ue
    res = dict(pass_id=11237, levels=[])
    for t in TARGETS:
        if t not in found:
            res["levels"].append(dict(level=t, found=False))
            continue
        F = exact_fidelity(mp, found[t], Cl)
        deg, c = identify(mp, F ** 2)
        rel = in_cubic_field(mp, F ** 2) if deg == 3 else None
        res["levels"].append(dict(level=t, F_digits=mp.nstr(F, 60), F2_minpoly_degree=deg, F2_minpoly=c,
                                  F2_in_Q_cos_2pi_9=None if rel is None else [int(v) for v in rel]))
        print(t, mp.nstr(F, 25), deg, c, flush=True)
    q = (1 + 2 * mp.cos(2 * mp.pi / 9)) / 3
    res["F_min"] = mp.nstr(q, 40)
    res["F_min_squared"] = mp.nstr(q ** 2, 40)
    res["F_min_squared_is_the_two_qutrit_level_0_712386014"] = bool(abs(q ** 2 - mp.mpf("0.712386014")) < 1e-9)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
