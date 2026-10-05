"""Pass 11504: is the cubed-phase test exact on all of PU(3), or only on Clifford+T words?

Pass 11500: on every one-qutrit Clifford+T word to depth 9, "for some anti-symplectic L of F3^2, |c(Lp)| = |c(p)| and
(c(Lp)/c(p))^3 is constant" coincides with reversibility.  Pass 11252's exact criterion needs more:
c(Lp) = mu omega^f(p) c(p) with f AFFINE.  A counterexample on PU(3) would be a unitary whose coefficients obey the relation
for a NON-affine f, and which is not reversible.

The relation is LINEAR in c: c lies in the mu-eigenspace of the 9x9 monomial matrix P_(L,f): (P c)(Lp) = omega^f(p) c(p).
So the candidates are unitaries U = sum_p c(p) W(p) inside such eigenspaces.  For random (L, non-affine f) and every
eigenvalue mu, minimise ||U U^dag - 1||^2 over the eigenspace (BFGS, 3 starts); every unitary found is tested with
Theorem 1's exact overlap decider.
"""

from __future__ import annotations

import itertools
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11350_linear_decider as L  # noqa: E402
import w33_pass11422_depth5_exact as E  # noqa: E402
import w33_pass11500_exact_depth9 as P9  # noqa: E402

OUT = ROOT / "data" / "w33_pass11504_cubed_phase_test.json"


def run(trials=4000, seed=11504):
    E._init()
    w = np.exp(2j * np.pi / 3)
    labs = [(a, b) for a in range(3) for b in range(3)]
    anti = [np.array([[a, b], [c, d]]) for a, b, c, d in itertools.product(range(3), repeat=4) if (a * d - b * c) % 3 == 2]
    Lperm = [[labs.index(tuple(int(x) for x in (S @ np.array(p)) % 3)) for p in labs] for S in anti]
    Aff = np.concatenate([np.ones((9, 1), np.int64), np.array(labs)], 1) % 3
    rng = np.random.default_rng(seed)
    st = Counter()
    found = []
    for _ in range(trials):
        li = int(rng.integers(24))
        f = rng.integers(3, size=9)
        if L.solve_affine(Aff, f % 3) is not None:
            continue
        st["non-affine (L, f) tried"] += 1
        Pm = np.zeros((9, 9), complex)
        for p in range(9):
            Pm[Lperm[li][p], p] = w ** f[p]
        vals, vecs = np.linalg.eig(Pm)
        for mu in np.unique(np.round(vals, 8)):
            Q, _ = np.linalg.qr(vecs[:, np.abs(vals - mu) < 1e-6])
            d = Q.shape[1]

            def loss(x):
                U = np.einsum('p,pij->ij', Q @ (x[:d] + 1j * x[d:]), P9.WP)
                return np.linalg.norm(U @ U.conj().T - np.eye(3)) ** 2

            best = min((minimize(loss, rng.normal(size=2 * d), method='BFGS') for _ in range(3)), key=lambda r: r.fun)
            st["eigenspaces searched"] += 1
            if best.fun < 1e-20:
                U = np.einsum('p,pij->ij', Q @ (best.x[:d] + 1j * best.x[d:]), P9.WP)
                ov = float(P9.X.overlaps(U[None])[0])
                st["unitaries found"] += 1
                st["... reversible (overlap = 3)" if ov > 3 - 1e-7 else "... NOT reversible (counterexample)"] += 1
                found.append(dict(L=li, f=[int(x) for x in f], mu=[float(mu.real), float(mu.imag)], overlap=ov,
                                  prefilter=bool(P9.prefilter(U[None])[0])))
    res = dict(pass_id=11504, trials=trials, counts=dict(st), unitaries=found)
    ce = st["... NOT reversible (counterexample)"]
    res["paper_sentence"] = (
        r"on $PU(3)$, a search of %d non-affine phase relations found %d unitaries, %s." % (
            st["eigenspaces searched"], st["unitaries found"],
            "all reversible, so no counterexample" if ce == 0 else "%d of them not reversible" % ce))
    res["ledger_phrase"] = ("%d unitaries on non-affine relations, %d counterexamples" % (st["unitaries found"], ce))
    print(json.dumps(res, indent=1), flush=True)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1)


if __name__ == "__main__":
    main()
