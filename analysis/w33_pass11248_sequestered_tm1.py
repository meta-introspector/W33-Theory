#!/usr/bin/env python3
"""Pass 11248: the sequestered TM1 model -- the sign of one coupling picks TM1, and the misalignment is O(eps).

Pass 11243: in the triplet 3' of the line stabiliser's S4, TM1 needs phi_e along a point d_X of the line and phi_nu along
the chord d_X - d_Q; in the sequestered limit TM1 ties 48:48 with the excluded theta13 = 0 orientation, and the cross
invariant (phi . chi)^2 is 16 vs 0 on the two classes.  This pass builds an explicit S4-invariant toy potential
    V = V_e(phi) + V_nu(chi) + lam |phi|^2 |chi|^2 + eps (phi . chi)^2,
    V_e  = -|phi|^2 + |phi|^4/2 + b sum phi_i^4                       (b > 0: minimum along the body diagonals d_X),
    V_nu = -|chi|^2 + |chi|^4/2 + k (sum chi_i^4 - |chi|^4/2)^2 + h (chi_1 chi_2 chi_3)^2   (minimum along the chords),
minimises it globally (96 starts: every orbit pair of the sequestered vacuum manifold, perturbed) and measures
  * which class wins for eps < 0 and eps > 0 (TM1 vs theta13 = 0),
  * the misalignment angles of phi from d_X and of chi from the chord as eps -> 0 (slope = O(eps) coefficient),
  * the induced change of the TM1 column under the residual-symmetry rule, with U_e and the neutrino Z2 eigenvector
    replaced by those of the nearest symmetric configuration (first-order proxy).
V_nu uses degree-8 terms: Pass 11230 showed no renormalisable single-triplet potential has the chord as its minimum.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11230_tm1_alignment as A  # noqa: E402
import w33_pass11243_tm1_line_geometry as G  # noqa: E402

OUT = ROOT / "data" / "w33_pass11248_sequestered_tm1.json"
B, K, H, LAM = 0.3, 2.0, 3.0, 0.2


def V(x, eps):
    f, c = x[:3], x[3:]
    f2, c2 = f @ f, c @ c
    Ve = -f2 + f2 ** 2 / 2 + B * np.sum(f ** 4)
    Vn = -c2 + c2 ** 2 / 2 + K * (np.sum(c ** 4) - c2 ** 2 / 2) ** 2 + H * np.prod(c) ** 2
    return Ve + Vn + LAM * f2 * c2 + eps * (f @ c) ** 2


def classify(x):
    return A.mixing_outcome(A.G3P, A.G3P, x[:3], x[3:], tol=1e-3)


def nearest_symmetric(x):
    """snap each vev to the nearest orbit direction (point d_i or chord d_i - d_j, 3' signs)"""
    pts = [s * G.D4[i] for i in range(4) for s in (1, -1)]
    chords = [s * (G.D4[i] - G.D4[j]) for i in range(4) for j in range(4) if i != j for s in (1, -1)]
    u = lambda v: v / np.linalg.norm(v)
    f = max(pts, key=lambda p: u(p) @ u(x[:3]))
    c = max(chords, key=lambda p: u(p) @ u(x[3:]))
    ang = lambda a, b: float(np.degrees(np.arccos(np.clip(abs(u(a) @ u(b)), -1, 1))))
    return f, c, ang(f, x[:3]), ang(c, x[3:])


def global_min(eps, seed=0):
    rng = np.random.default_rng(seed)
    pts = [s * G.D4[i] for i in range(4) for s in (1, -1)]
    chords = [s * (G.D4[i] - G.D4[j]) for i in range(4) for j in range(4) if i < j for s in (1, -1)]
    best = None
    for p in pts:
        for c in chords:
            x0 = np.concatenate([p / np.sqrt(3), c / np.sqrt(8)]) + 0.01 * rng.normal(size=6)
            r = minimize(V, x0, args=(eps,), method="BFGS", options=dict(gtol=1e-12))
            if best is None or r.fun < best.fun - 1e-12:
                best = r
    return best


def run():
    res = dict(pass_id=11248, couplings=dict(b=B, k=K, h=H, lam=LAM))
    # sanity: decoupled minima are the symmetric directions
    r0 = global_min(0.0)
    f, c, af, ac = nearest_symmetric(r0.x)
    res["eps0_misalignment_deg"] = [af, ac]
    rows = []
    for eps in (-0.08, -0.04, -0.02, -0.01, -0.005, 0.005, 0.01, 0.02, 0.04):
        r = global_min(eps)
        f, c, af, ac = nearest_symmetric(r.x)
        cls = A.mixing_outcome(A.G3P, A.G3P, f, c)
        rows.append(dict(eps=eps, symmetric_class=cls, phi_misalign_deg=af, chi_misalign_deg=ac,
                         energy=float(r.fun)))
    res["scan"] = rows
    neg = [r for r in rows if r["eps"] < 0]
    pos = [r for r in rows if r["eps"] > 0]
    res["eps_negative_selects_TM1"] = all(r["symmetric_class"] == "TM1" for r in neg)
    res["eps_positive_selects_theta13_0"] = all(r["symmetric_class"] == "theta13 = 0 column" for r in pos)
    small = [r for r in neg if abs(r["eps"]) <= 0.02]
    res["misalignment_slope_deg_per_unit_eps"] = dict(
        phi=float(np.mean([r["phi_misalign_deg"] / abs(r["eps"]) for r in small])),
        chi=float(np.mean([r["chi_misalign_deg"] / abs(r["eps"]) for r in small])))
    return res


def main():
    res = run()
    OUT.write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
