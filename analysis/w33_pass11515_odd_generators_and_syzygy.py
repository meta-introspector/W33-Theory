"""Pass 11515: the five time-odd generators of a qutrit state written down, and their first relation.

Pass 11502 counted one new odd generator in each degree 6..10 (none in 11..13) over the even ring.  Here:
  (1) EXPLICIT GENERATORS.  For each degree d = 6..10 the simplest stabiliser-probability monomial m whose odd Reynolds
      average R_m = (1/432) sum_g sgn(g) m(g psi) lies outside the span of {even x lower generators}; monomials are tried in
      order of the number of distinct stabiliser states, then lexicographically.
  (2) THE FIRST SYZYGY.  If the odd module were free on generators of degrees 6..10 over the even ring E, degree k would
      hold sum_i E_(k - d_i) products.  With Pass 11491's exact E_k, O_k this count equals O_k for k = 6..12 and exceeds it by
      exactly one at k = 13: the module is free through degree 12 and has its first relation in degree 13.  Verified: the
      products e * g_i (e over an even basis) have rank O_k, with full column rank for k <= 12 and one dependency at 13.
Numerical ranks on 160 random unit rays (a bidegree-(k,k) form is fixed by its values on the sphere); gaps recorded.
"""

from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11502_time_odd_module as T  # noqa: E402

OUT = ROOT / "data" / "w33_pass11515_odd_generators_and_syzygy.json"
NAMES = ["Z0", "Z1", "Z2", "X0", "X1", "X2", "Y0", "Y1", "Y2", "W0", "W1", "W2"]   # MUBs Z, X, Y = XZ, W = XZ^2 (Pass 11434)


def monomials(deg):
    """candidate monomials ordered by number of distinct states, then exponent pattern, then states"""
    out = []
    for nstates in range(1, 6):
        for states in itertools.combinations(range(12), nstates):
            for exps in itertools.product(range(1, deg + 1), repeat=nstates):
                if sum(exps) == deg:
                    out.append(tuple(s for s, e in zip(states, exps) for _ in range(e)))
    return out


def colrank(cols):
    if not cols:
        return 0
    return T.rank(np.stack(cols, 1))


def run(kmax=13):
    I = T.Invariants()
    E = {k: I.basis(k, False, T.MOLIEN_E[k]) for k in range(0, kmax + 1)}
    gens = {}
    rows = {}
    for d in range(6, 11):
        span = [E[k][:, a] * gens[g] for g in gens for k in [d - g] if k >= 0 for a in range(E[k].shape[1])]
        r0 = colrank(span)
        for m in monomials(d):
            v = I.reynolds(m, True)
            if np.linalg.norm(v) < 1e-12:
                continue
            v = v / np.linalg.norm(v)
            if colrank(span + [v]) > r0:
                gens[d] = v
                rows[str(d)] = dict(monomial=[NAMES[i] for i in m], mubs=sorted({NAMES[i][0] for i in m}),
                                    exponents=sorted([m.count(s) for s in set(m)], reverse=True),
                                    rank_before=r0, O_d=T.MOLIEN_O[d], gap=T.GAPS[-1])
                break
        print(d, rows[str(d)], flush=True)
    syz = {}
    for k in range(6, kmax + 1):
        cols = [E[k - g][:, a] * gens[g] for g in gens if k - g >= 0 for a in range(E[k - g].shape[1])]
        r = colrank(cols)
        syz[str(k)] = dict(products=len(cols), rank=r, O_k=T.MOLIEN_O[k], relations=len(cols) - r, gap=T.GAPS[-1])
        print(k, syz[str(k)], flush=True)
    predicted = {str(k): sum(T.MOLIEN_E[k - g] for g in range(6, 11) if k - g >= 0) - T.MOLIEN_O[k]
                 for k in range(6, kmax + 1)}
    return dict(pass_id=11515, generators=rows, syzygies=syz, predicted_relations_from_molien=predicted)


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
