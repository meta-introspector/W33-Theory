#!/usr/bin/env python3
"""Pass 11172: the whole negativity polygamy frontier of three qutrits -- spatial polygamy is a thin bulge, and it is
carried (numerically, for every weight) by one closed-form family.

Pass 11167 found the closed form on the symmetry class of the 4/sqrt15 optimum,
    psi = alpha|000> + beta sum_{a=1,2} |aa0> + gamma sum_{a=1,2} |a0a>,   s = beta^2, t = gamma^2, alpha^2 = 1 - 2s - 2t,
    N_AB = sqrt(t^2 + 4 alpha^2 s) - t + s,     N_AC = sqrt(s^2 + 4 alpha^2 t) - s + t.
Question: the whole achievable region of (N_AB, N_AC) over ALL three-qutrit pure states.  For each weight lambda we compare
the family maximum of N_AB + lambda N_AC with a GLOBAL weighted dual see-saw (psi <- top eigenvector of
-P1^{T_B} (x) I - lambda P2^{T_C} (x) I, P_i <- negative-eigenspace projectors; monotone), 150 random starts per lambda.
Results:
  * family = global to < 1e-12 at every lambda tested (0, 0.5, 0.9, 0.93, ..., 1);
  * for lambda <= lambda_c = 0.92143 the maximum is exactly 1, at the corner (N_AB, N_AC) = (1, 0) (a maximally entangled
    AB pair, C decoupled): no polygamy at all;
  * for lambda_c < lambda <= 1 the maximiser moves along the family arc, from (0.7206, 0.3033) at lambda_c through the
    symmetric point (2/sqrt15, 2/sqrt15) at lambda = 1 (and symmetrically beyond).
So the upper-right boundary of the convex hull of the polygamy region is: the segment from (1, 0) to the arc endpoint,
the family arc, and the mirror segment to (0, 1).  Beyond the monogamous line N_AB + N_AC = 1, space allows only a thin
bulge (at most 0.0328), and only for nearly balanced sharing.  In time the same sum is 2 (Pass 11148).
Scope: exact on the family; numerical (dual and primal searches) for the claim that the family carries the boundary.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11167_polygamy_global as G  # noqa: E402

OUT = ROOT / "data" / "w33_pass11172_polygamy_frontier.json"


def fam(s, t):
    al = max(0.0, 1 - 2 * s - 2 * t)
    return np.sqrt(t * t + 4 * al * s) - t + s, np.sqrt(s * s + 4 * al * t) - s + t


def family_max(lam):
    S, T = np.meshgrid(np.linspace(0, 0.5, 801), np.linspace(0, 0.5, 801))
    ok = 2 * S + 2 * T <= 1
    Al = np.clip(1 - 2 * S - 2 * T, 0, None)
    F = np.where(ok, (np.sqrt(T ** 2 + 4 * Al * S) - T + S) + lam * (np.sqrt(S ** 2 + 4 * Al * T) - S + T), -9)
    k = np.unravel_index(F.argmax(), F.shape)
    x0 = np.array([S[k], T[k]])

    def neg(x):
        s, t = x
        if s < 0 or t < 0 or 2 * s + 2 * t > 1:
            return 9
        a, b = fam(s, t)
        return -(a + lam * b)
    r = minimize(neg, x0, method='Nelder-Mead', options=dict(xatol=1e-13, fatol=1e-15, maxiter=4000))
    cands = [(-r.fun, r.x), (F.max(), x0), (sum(w * v for w, v in zip((1, lam), fam(1 / 3, 0))), np.array([1 / 3, 0.0])),
             (sum(w * v for w, v in zip((1, lam), fam(0, 1 / 3))), np.array([0.0, 1 / 3]))]
    best, x = max(cands, key=lambda c: c[0])
    s, t = x
    return float(best), [float(s), float(t)], [float(v) for v in fam(s, t)]


def interior_max(lam):
    """family maximum of N_AB + lam N_AC away from the corners (t >= 0.02)"""
    best = -9
    for s0, t0 in [(0.2, 0.1), (0.25, 0.08), (0.15, 0.13), (0.3, 0.05)]:
        def neg(x):
            s, t = x
            if s < 0 or t < 0.02 or 2 * s + 2 * t > 1:
                return 9
            a, b = fam(s, t)
            return -(a + lam * b)
        r = minimize(neg, [s0, t0], method='Nelder-Mead', options=dict(xatol=1e-14, fatol=1e-16, maxiter=6000))
        best = max(best, -r.fun)
    return best


def critical_lambda():
    lo, hi = 0.85, 0.95
    for _ in range(60):
        mid = (lo + hi) / 2
        if interior_max(mid) > 1:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def seesaw_weighted(rng, lam, iters=500):
    psi = rng.normal(size=27) + 1j * rng.normal(size=27)
    psi /= np.linalg.norm(psi)
    val = -1
    for _ in range(iters):
        T = psi.reshape(3, 3, 3)
        Ps = []
        for M in (T.reshape(9, 3), T.transpose(0, 2, 1).reshape(9, 3)):
            w, V = np.linalg.eigh(G.pt(M @ M.conj().T))
            Vn = V[:, w < 0]
            Ps.append(Vn @ Vn.conj().T)
        H1 = np.kron(G.pt(Ps[0]), np.eye(3))
        H2 = np.kron(G.pt(Ps[1]), np.eye(3)).reshape(3, 3, 3, 3, 3, 3).transpose(0, 2, 1, 3, 5, 4).reshape(27, 27)
        H = -(H1 + lam * H2)
        w, V = np.linalg.eigh((H + H.conj().T) / 2)
        new = float(w[-1])
        psi = V[:, -1]
        if abs(new - val) < 1e-13:
            break
        val = new
    return val


def summarize(n_lam=11, starts=150):
    rng = np.random.default_rng(11172)
    rows = []
    lams = [0.0, 0.5, 0.9, 0.93, 0.935, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0]
    for lam in lams:
        fm, st, pair = family_max(lam)
        gl = max(seesaw_weighted(rng, lam) for _ in range(starts))
        rows.append(dict(lam=float(lam), family=fm, family_point=st, family_pair=pair, global_seesaw=gl, gap=gl - fm))
    res = dict(pass_id=11172, rows=rows, max_gap=max(r['gap'] for r in rows),
               family_equals_global=all(abs(r['gap']) < 1e-7 for r in rows),
               symmetric_point=[2 / np.sqrt(15), 2 / np.sqrt(15)], critical_lambda=critical_lambda())
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    r = summarize()
    for row in r['rows']:
        print(round(row['lam'], 2), round(row['family'], 9), round(row['global_seesaw'], 9), '%.1e' % row['gap'], [round(x, 4) for x in row['family_pair']])
    print('max gap', r['max_gap'], r['family_equals_global'], 'critical lambda', r['critical_lambda'])
