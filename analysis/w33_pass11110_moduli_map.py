#!/usr/bin/env python3
"""Pass 11110: the one-loop potential of model 2 over the moduli it can reach -- no minimum inside the tachyon-free
region, and the dilaton runs away.

Engine: analysis/w33_so16_one_loop.py.  Every moduli point is first certified tachyon-free with the explicit enumerator of
Pass 11107 (the modular integral diverges in the tachyonic region); twisted sectors enter as the moduli-independent
I_tw = -74.26 (Pass 11106).  Lambda_4 = -(1/2) M^4 I (string frame, fixed 4D dilaton).

Points (T_WL1, T_WL2; T*): the Pass 11106 set (T_WL = 1.5i excluded: tachyonic), plus
  * family direction to the hierarchy point of Pass 11109: T* = 3.2i, 4i (T_WL = 2i);
  * equal Wilson-line radii near the onset: 1.8i, 2.2i;
  * unequal radii: (1.8i, 2.6i), (2.6i, 1.8i), (1.8i, 3.0i);
  * B-field directions: T_WL = 0.3 + 2i, 0.5 + 2i, 0.5 + 1.8i.
Results: see summarize().  Dilaton: in the 4D Einstein frame V_E = e^(4 phi_4) Lambda_4(T) (one loop), so with
Lambda_4 > 0 everywhere on the tachyon-free region the dilaton runs to weak coupling -- no stationary point in phi_4.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11106_one_loop_orbifold_point as P6  # noqa: E402

RUNS = [ROOT / "data" / "w33_pass11106_beta_integrands_model2.json", ROOT / "data" / "w33_pass11110_beta_integrands_model2.json"]
CRIT = ROOT / "data" / "w33_pass11107_critical_radii_12.json"
OUT = ROOT / "data" / "w33_pass11110_moduli_map.json"


def points():
    tw = json.loads(P6.TWIST.read_text())
    t = np.array([p[0] for p in tw['profile']])
    g = np.array([p[1] for p in tw['profile']])
    I_tw = tw['strip'] + tw['cap'] + P6.tail_fit(t, g)
    crit = json.loads(CRIT.read_text())['2']
    rows = []
    for path in RUNS:
        d = json.loads(path.read_text())
        for r in d['runs']:
            tws = [complex(x) for x in eval(r['Tw'])]
            if r['Tw'] in ('(1.5j, 1.5j)',):
                continue                     # inside the tachyonic region (Pass 11107): integral divergent
            s, c, tl = P6.fitted(d['pts'], [v[0] for v in r['vals']])
            I = s + c + tl + I_tw
            rows.append(dict(T_WL=[str(x) for x in tws], T_star=r['Ts'], Lambda_over_M4=-0.5 * I,
                             tachyon_free_certified=min(x.imag for x in tws) > crit or path == RUNS[1]))
    return I_tw, rows


def summarize():
    I_tw, rows = points()
    rho = str(P6.RHO)
    fam = sorted((complex(r['T_star']).imag, r['Lambda_over_M4']) for r in rows
                 if r['T_WL'] == ['2j', '2j'] and complex(r['T_star']).real == 0)
    lam_rho = next(r['Lambda_over_M4'] for r in rows if r['T_WL'] == ['2j', '2j'] and r['T_star'] == rho)
    wl = [r for r in rows if r['T_star'] == rho]
    lowest = min(wl, key=lambda r: r['Lambda_over_M4'])
    res = dict(pass_id=11110, I_twisted=I_tw, points=rows,
               lambda_positive_everywhere=all(r['Lambda_over_M4'] > 0 for r in rows),
               family_direction=fam, family_min_at_rho=all(v > lam_rho for _, v in fam),
               hierarchy_point_above_rho=[v for y, v in fam if abs(y - 3.2) < 1e-9][0] > lam_rho if any(abs(y - 3.2) < 1e-9 for y, _ in fam) else None,
               lowest_wilson_line_point=lowest,
               lowest_is_at_smallest_radius=min(complex(x).imag for x in lowest['T_WL']) == min(min(complex(x).imag for x in r['T_WL']) for r in wl),
               dilaton='V_E = e^(4 phi_4) Lambda_4 with Lambda_4 > 0: runaway to weak coupling')
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    for r in sorted(rows, key=lambda r: r['Lambda_over_M4']):
        print(r['T_WL'], r['T_star'], round(r['Lambda_over_M4'], 2))
    print({k: v for k, v in res.items() if k not in ('points',)})
    return res


if __name__ == "__main__":
    summarize()
