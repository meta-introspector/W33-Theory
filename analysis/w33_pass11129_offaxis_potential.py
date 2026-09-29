#!/usr/bin/env python3
"""Pass 11129: the one-loop potential off the B = 0 axis -- the B-field is a near-flat direction; the potential does not
select it during the fall to the neutral onset.

Scan: analysis/w33_pass11129_scan_offaxis.py (the Pass 11122 engine with both Wilson-line tori at x + i y); one sqrt3
model (36621) and one golden model (40521), y = y_c + 0.07 and 2.0, x = 0.25 and 0.5; x = 0 from Pass 11122 (same grid,
same K = 3, same tail fit, so the differences are free of the common quadrature error).
Frozen: data/w33_pass11129_offaxis_beta_integrands.json.

Results (Lambda_beta in M^4):
  * Lambda(x) - Lambda(0) = A (1 - cos 2 pi x)/2 within the sampling: the x = 0.5 difference is 1.8-2.0 x the x = 0.25 one
    in all four (model, y) pairs;
  * just above the onset A = +0.011 (36621) and +0.027 (40521): B = 0 is a MINIMUM there, so the B-field is pulled back to
    zero where it matters -- the neutral exit (a B = 0 feature, Pass 11125) is the one the flow reaches;
  * at y = 2.0 A = -0.004, -0.006: B = 0 is a slight maximum;
  * |A| / (dLambda/dy) <= 5e-5: B is effectively frozen during the radial fall (the B-dependence comes only from winding
    states, exponentially suppressed in Im T).
Scope: two models, two radii, x in {0, 0.25, 0.5} (reflection x -> -x); A is ~2e-5 of Lambda, so its sign is reliable
only because the grid and fit are shared; no dynamical (kinetic) analysis.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11106_one_loop_orbifold_point as P6  # noqa: E402

DATA = ROOT / "data" / "w33_pass11129_offaxis_beta_integrands.json"
BASE = ROOT / "data" / "w33_pass11122_six_potentials.json"
OUT = ROOT / "data" / "w33_pass11129_offaxis_potential.json"


def summarize():
    d = json.loads(DATA.read_text())
    base = json.loads(BASE.read_text())['models']
    rows = {}
    for r in d['runs']:
        s, c, t = P6.fitted(d['pts'], [v[0] for v in r['vals']])
        L0 = dict((round(y, 3), v) for y, v in base[r['model']]['profile'])[round(r['y'], 3)]
        rows.setdefault(f"{r['model']}|{r['y']}", {})[str(r['x'])] = -0.5 * (s + c + t) - L0
    fits = {}
    for k, v in rows.items():
        m = k.split('|')[0]
        slope = abs(base[m]['slope_near_onset'])
        fits[k] = dict(diff_0p25=v['0.25'], diff_0p5=v['0.5'], ratio=v['0.5'] / v['0.25'], A=v['0.5'],
                       A_over_radial_slope=abs(v['0.5']) / slope)
    near = {k: f for k, f in fits.items() if not k.endswith('|2.0')}
    far = {k: f for k, f in fits.items() if k.endswith('|2.0')}
    res = dict(pass_id=11129, fits=fits,
               cosine_shape=all(1.6 < f['ratio'] < 2.2 for f in fits.values()),
               B0_minimum_near_onset=all(f['A'] > 0 for f in near.values()),
               B0_maximum_at_2=all(f['A'] < 0 for f in far.values()),
               max_A_over_radial=max(f['A_over_radial_slope'] for f in fits.values()))
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
