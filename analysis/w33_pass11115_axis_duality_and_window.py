#!/usr/bin/env python3
"""Pass 11115: the T-duality turnaround of the Wilson-line potential never falls inside a tachyon-free window.

Idea: at large volume Lambda grows like the volume in ANY T-dual description, so along the B = 0 axis of the Wilson-line
tori V(y) must turn around somewhere at small radius; if the model were tachyon-free there, the Wilson-line moduli would
be stabilised.  Tests (explicit enumerator of Pass 11107 with windings/momenta to 7 and cached torus tables; one-loop
engine with K = 6 windings below Im T = 0.9, converged to 1e-9 at Im T = 0.5):
  A. axis tachyon levels, both Wilson-line tori at i y, y = 0.3 ... 1.0 (data/w33_pass11115_axis_tachyons.json):
       * model 10 shows an approximate Fricke duality y -> 1/(3y) (paired levels agree within 0.03 at the nearest
         samples), deepest at the self-dual radius 1/sqrt3 (Delta = -1/6): the duality-forced extremum lies INSIDE the
         tachyonic region; the other models show no such symmetry on the sampled points;
       * models 2 and 53: tachyonic on the whole segment, increasingly so;
       * model 77 (a full survivor): tachyon-free at y = 0.87, 0.9, 1.0, tachyonic at y = 0.85, 0.8 (Delta = -0.0093, charge
         +-1/3): onset at ~sqrt3/2.
  B. model 77's potential on its tachyon-free window (data/w33_pass11115_axis_integrands_model77.json):
       Lambda_beta = 993.7, 633.2, 350.8, 238.6, 180.9, 146.7, 116.0 at y = 2.5, 2.0, 1.5, 1.25, 1.1, 1.0, 0.9 --
       monotone down to the onset.
Reading: the turnaround that duality requires does not occur inside a tachyon-free window in any tested model -- model
10's duality extremum sits at its tachyonic self-dual radius, and model 77's potential falls monotonically to its onset.
The Wilson-line moduli of the survivors are not stabilised inside the tachyon-free region.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
AXIS = ROOT / "data" / "w33_pass11115_axis_tachyons.json"
POT77 = ROOT / "data" / "w33_pass11115_axis_integrands_model77.json"
OUT = ROOT / "data" / "w33_pass11115_axis_duality_and_window.json"


def fricke_check(levels, tol=0.03):
    """compare Delta(y) with Delta(1/(3y)) where both are sampled (nearest sample within 0.05)"""
    ys = sorted(levels)
    pairs = []
    for y in ys:
        yd = 1 / (3 * y)
        near = min(ys, key=lambda z: abs(z - yd))
        if abs(near - yd) < 0.05 and levels[y] is not None and levels[near] is not None:
            pairs.append((y, near, levels[y], levels[near]))
    return pairs, all(abs(a - b) < tol for _, _, a, b in pairs) if pairs else None


def summarize():
    import w33_pass11106_one_loop_orbifold_point as P6
    ax = json.loads(AXIS.read_text())
    per = {}
    for m, y, d in ax:
        per.setdefault(m, {})[y] = d
    fr = {m: fricke_check(v) for m, v in per.items()}
    d = json.loads(POT77.read_text())
    prof = []
    for r in d['runs']:
        s, c, t = P6.fitted(d['pts'], [v[0] for v in r['vals']])
        prof.append((complex(eval(r['Tw'])[0]).imag, -(s + c + t) / 2))
    prof.sort()
    res = dict(pass_id=11115, axis_levels=per,
               fricke={m: dict(pairs=v[0], symmetric=v[1]) for m, v in fr.items()},
               model77_window_start=min(y for y, dl in per['77'].items() if dl is None),
               model77_first_tachyonic=max(y for y, dl in per['77'].items() if dl is not None),
               model77_lambda_profile=prof,
               model77_monotone=all(a[1] < b[1] for a, b in zip(prof, prof[1:])))
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    print(json.dumps({k: v for k, v in res.items() if k != 'axis_levels'}, indent=1))
    return res


if __name__ == "__main__":
    summarize()
