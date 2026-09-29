#!/usr/bin/env python3
"""Pass 11118: the gradient flow of the one-loop potential over the two Wilson-line radii of model 57 -- which tachyon
boundary does the vacuum reach first?

Pass 11112: model 57's Standard-Model-neutral winding tachyons live on Wilson-line torus 1, the charged ones on torus 2,
both with critical radius Im T = 1.512 (Pass 11107), and pairwise comparisons showed the potential pulls torus 1 down
faster.  Here the full 6 x 6 grid Lambda(y1, y2), y in {1.55, 1.7, 1.85, 2.0, 2.3, 2.6} (B = 0, family torus at rho; one-
loop engine, K = 4 windings; frozen in data/w33_pass11118_grid_integrands_model57.json) is interpolated (bicubic spline)
and the gradient flow with the moduli-space metric (kinetic term dy^2 / (4 y^2) per torus),
    dy_i / dt = -4 y_i^2 dLambda / dy_i,
is integrated from nine starting points until the flow leaves the grid at y = 1.55 (just above the common critical radius
1.512).  The boundary crossed first tells which tachyon the moduli meet.
Scope: overdamped (gradient) flow -- the direction of steepest descent, not the cosmological trajectory; family modulus
fixed at rho; the twisted sectors add a moduli-independent constant.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import RectBivariateSpline

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
GRID = ROOT / "data" / "w33_pass11118_grid_integrands_model57.json"
OUT = ROOT / "data" / "w33_pass11118_wilson_line_flow.json"
YCRIT = 1.512
STARTS = [(2.6, 2.6), (2.3, 2.6), (2.6, 2.3), (2.0, 2.6), (2.6, 2.0), (2.3, 2.3), (2.0, 2.0), (1.85, 2.3), (2.3, 1.85)]


def grid():
    import w33_pass11106_one_loop_orbifold_point as P6
    d = json.loads(GRID.read_text())
    L = {}
    for r in d['runs']:
        s, c, t = P6.fitted(d['pts'], [v[0] for v in r['vals']])
        L[tuple(r['pair'])] = -(s + c + t) / 2
    ys = sorted({k[0] for k in L})
    return np.array(ys), np.array([[L[(a, b)] for b in ys] for a in ys])


def flow(dt=2e-6, nmax=400000, starts=None):
    ys, Z = grid()
    sp = RectBivariateSpline(ys, ys, Z, kx=3, ky=3)
    out = []
    for y0 in (starts or STARTS):
        y = np.array(y0, float)
        for _ in range(nmax):
            g = np.array([sp(y[0], y[1], dx=1)[0, 0], sp(y[0], y[1], dy=1)[0, 0]])
            y = y - dt * 4 * y ** 2 * g
            if y.min() <= ys[0]:
                break
        out.append(dict(start=list(y0), end=[float(y[0]), float(y[1])],
                        first_boundary='torus1 (neutral)' if y[0] <= y[1] else 'torus2 (charged)'))
    asym = [float(Z[i, j] - Z[j, i]) for i in range(len(ys)) for j in range(len(ys)) if i < j]
    return ys, Z, out, asym


def summarize():
    ys, Z, out, asym = flow()
    sep = {}
    for top in (2.6, 2.3, 2.0):
        deltas = (0.0, 0.005, 0.01, 0.015, 0.02, 0.04)
        o = flow(starts=[(top, top - d) for d in deltas])[2]
        sep[str(top)] = [(d, r['first_boundary']) for d, r in zip(deltas, o)]
    res = dict(pass_id=11118, separatrix=sep, ys=ys.tolist(), Lambda=Z.tolist(), flows=out,
               neutral_first=sum(1 for f in out if f['first_boundary'].startswith('torus1')), flows_total=len(out),
               antisymmetric_part_sign=dict(negative=sum(1 for a in asym if a < 0), positive=sum(1 for a in asym if a > 0)),
               minimum_inside=bool(np.unravel_index(np.argmin(Z), Z.shape) not in [(0, k) for k in range(len(ys))] + [(k, 0) for k in range(len(ys))]))
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    print(json.dumps({k: v for k, v in res.items() if k not in ('Lambda',)}, indent=1))
    return res


if __name__ == "__main__":
    summarize()
