#!/usr/bin/env python3
"""Pass 11122: the one-loop potential of the six neutral-exit survivors (Pass 11119) along the B = 0 diagonal of the
two Wilson-line tori, family torus at rho -- does the flow lead to the neutral boundary, and how steeply?

Engine: analysis/w33_so16_one_loop.py (Witten-twisted beta sectors, K = 3 Wilson-line cosets), integrated on the
fundamental-domain grid of Pass 11106 with the tail fit c0 + c1 exp(-a tau2).  The theta-twisted sectors do not depend
on the Wilson-line moduli (Pass 11106), so the moduli dependence of Lambda = -(1/2) M^4 I is carried entirely by
I_beta: Delta Lambda = -(1/2) Delta I_beta.  Frozen integrands: data/w33_pass11122_six_beta_integrands.json
(3 radii per model: y_c + 0.07 just above the neutral onset of Pass 11119, 2.0 and 2.5).

Results (Lambda_beta in units of M^4; y = Im T of both Wilson-line tori):
  * all six decrease monotonically from 2.5 through 2.0 to just above the neutral onset: the flow reaches the neutral
    boundary from the whole diagonal, with no barrier and no interior minimum;
  * the profile is the volume law Lambda_beta = a y^2 + b with a = 160.0-160.8 and b = -4.8 to -12.6 in every model
    (Lambda/y^2 = 156.7-159.2 and rising with y): the drive is the large-volume (n_B - n_F) term, common to the six,
    and the model-dependent part is a small NEGATIVE correction that grows near the onset -- the potential falls
    slightly faster than the volume there, never flatter;
  * dLambda/dy at the onset is 575-615 M^4 per unit Im T: a steep, O(1)-in-string-units force; the approach is not
    slow roll.
Scope: B = 0 diagonal, 3 radii per model, family torus at rho; K = 3 cosets; moduli-independent theta sectors omitted
(they shift Lambda by a constant, so the absolute sign of Lambda is not addressed here -- Pass 11106 has it positive).
"""
from __future__ import annotations

import json
import sys

import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11106_one_loop_orbifold_point as P6  # noqa: E402

DATA = ROOT / "data" / "w33_pass11122_six_beta_integrands.json"
OUT = ROOT / "data" / "w33_pass11122_six_potentials.json"


def summarize():
    d = json.loads(DATA.read_text())
    prof = {}
    for r in d['runs']:
        s, c, t = P6.fitted(d['pts'], [v[0] for v in r['vals']])
        prof.setdefault(r['model'], []).append((r['y'], -0.5 * (s + c + t)))
    models = {}
    for m, rows in prof.items():
        rows.sort()
        (y0, L0), (y1, L1), (y2, L2) = rows
        models[m] = dict(profile=rows, monotone_toward_onset=bool(L0 < L1 < L2),
                         slope_near_onset=(L1 - L0) / (y1 - y0), slope_far=(L2 - L1) / (y2 - y1),
                         volume_fit=[float(x) for x in np.polyfit([y * y for y, _ in rows], [L for _, L in rows], 1)],
                         lambda_over_y2=[L / (y * y) for y, L in rows])
    res = dict(pass_id=11122, models=models,
               all_monotone=all(v['monotone_toward_onset'] for v in models.values()),
               volume_law_all=all(155 < v['volume_fit'][0] < 165 and -15 < v['volume_fit'][1] < 0 for v in models.values()),
               subvolume_negative_all=all(a < b < c for a, b, c in (v['lambda_over_y2'] for v in models.values())))
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
