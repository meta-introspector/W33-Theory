"""Pass 11191: polygamy across local dimension -- prior art, and the qutrit bulge is not special.

The optimum psi* = sqrt(7/15)|000> + sqrt(2/15) sum_{a=1,2} |a>(|0a> + |a0>) of N_AB + N_AC (Passes 11148-11187) lies in
the family  d|000> + a sum_{j>=1} |j0j> + b sum_{j>=1} |jj0>  of Allen & Meyer, 'Polynomial monogamy relations for
entanglement negativity', PRL 118, 080402 (2017), arXiv:1502.04807, who conjecture from numerics that this family traces
the whole achievable negativity region for D > 2 (proved for qubits).  The repo's 'global optimality open' is that
conjecture specialised to the sum; the value 4/sqrt15 itself was not found in the literature (Pass 11188-11192 sweep).

This pass puts the qutrit number in context.  With the negativity normalised so a maximally entangled pair has N = 1,
N = (||rho^T||_1 - 1)/(D - 1) (for D = 3 this is the repo's N), it maximises N_AB + N_AC on the Allen-Meyer family for
D = 2..7 and, as a check that the family contains the optimum, over the full state space for D = 2 (random restarts).
The three-qubit maximum (8 sqrt2 - 4)/7 = 1.0448 exceeds the qutrit 4/sqrt15 = 1.0328: exceeding 1 is not a qutrit
effect, and the bulge shrinks with D.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "w33_pass11191_polygamy_across_dimension.json"


def negativity(psi, D, pair):
    T = psi.reshape(D, D, D)
    if pair == "AB":
        rho = np.einsum('abc,xyc->abxy', T, T.conj())
    else:
        rho = np.einsum('abc,xbz->acxz', T, T.conj())
    rt = rho.transpose(0, 3, 2, 1).reshape(D * D, D * D)     # partial transpose on the second party
    ev = np.linalg.eigvalsh((rt + rt.conj().T) / 2)
    return (np.abs(ev).sum() - 1) / (D - 1)


def family_state(p, D):
    d, a, b = p
    psi = np.zeros((D, D, D))
    psi[0, 0, 0] = d
    for j in range(1, D):
        psi[j, 0, j] = a
        psi[j, j, 0] = b
    return (psi / np.linalg.norm(psi)).ravel()


def f_family(p, D):
    psi = family_state(p, D)
    return negativity(psi, D, "AB") + negativity(psi, D, "AC")


def best_family(D, restarts=60, seed=0):
    rng = np.random.default_rng(seed)
    best = (-1, None)
    for _ in range(restarts):
        r = minimize(lambda p: -f_family(p, D), rng.normal(size=3), method="Nelder-Mead",
                     options=dict(xatol=1e-12, fatol=1e-14, maxiter=20000))
        if -r.fun > best[0]:
            best = (-r.fun, r.x)
    p = best[1]
    psi = family_state(p, D).reshape(D, D, D)
    return float(best[0]), dict(d2=float(psi[0, 0, 0] ** 2), a2=float(psi[1, 0, 1] ** 2), b2=float(psi[1, 1, 0] ** 2))


def best_full(D, restarts=40, seed=1):
    rng = np.random.default_rng(seed)
    n = D ** 3

    def f(x):
        psi = x[:n] + 1j * x[n:]
        psi = psi / np.linalg.norm(psi)
        return -(negativity(psi, D, "AB") + negativity(psi, D, "AC"))
    best = -1
    for _ in range(restarts):
        r = minimize(f, rng.normal(size=2 * n), method="BFGS", options=dict(maxiter=4000))
        r = minimize(f, r.x, method="Nelder-Mead", options=dict(maxiter=40000, xatol=1e-12, fatol=1e-14))
        best = max(best, -r.fun)
    return float(best)


def run():
    res = dict(pass_id=11191, normalisation="N = (||rho^T||_1 - 1)/(D - 1)")
    fam = {}
    for D in range(2, 8):
        v, w = best_family(D)
        fam[D] = dict(max_sum=v, weights=w)
    res["family_max"] = fam
    res["qubit_closed_form"] = (8 * np.sqrt(2) - 4) / 7
    res["qutrit_closed_form"] = 4 / np.sqrt(15)
    assert abs(fam[2]["max_sum"] - res["qubit_closed_form"]) < 1e-6
    assert abs(fam[3]["max_sum"] - res["qutrit_closed_form"]) < 1e-6
    res["full_space_D2"] = best_full(2)
    res["full_space_D2_matches_family"] = abs(res["full_space_D2"] - fam[2]["max_sum"]) < 1e-5
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=float)
    for k, v in res.items():
        print(k, ":", v)


if __name__ == "__main__":
    main()
