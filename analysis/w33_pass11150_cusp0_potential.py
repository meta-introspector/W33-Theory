#!/usr/bin/env python3
"""Pass 11150: the one-loop potential near the cusp at 0 (the tachyon-free region below the golden disk) tracks its Fricke
image and runs into the disk -- no minimum there either.

Scan: analysis/w33_pass11150_scan_cusp0.py -- model 40521 (golden), both Wilson-line tori at i y, family at rho;
y = 0.20 and 0.17 (tachyon-free: below the golden disk's lower edge 0.2205, Pass 11141) with K = 5 cosets (K = 5 and 7
agree to 1e-8 at y = 0.2, 0.15), and the Fricke images y' = 1/(3y) = 1.667, 1.961 with K = 3.
Frozen: data/w33_pass11150_cusp0_beta_integrands.json.
Results (Lambda_beta / M^4):
  * y = 0.20: 439.30  vs  Fricke image 1.667: 436.92  (0.5%)
  * y = 0.17: 610.58  vs  Fricke image 1.961: 608.81  (0.3%)
  * at fixed tau the partition function matches its Fricke image to 3-4e-4 (y = 0.2) and 6-7e-5 (y = 0.15);
  * Lambda FALLS as y rises toward the disk edge (dLambda/dy ~ -5.7e3 between 0.17 and 0.20): from the cusp-0 side the
    Wilson-line modulus rolls up into the tachyonic disk; toward y -> 0 Lambda grows (the image of large volume).
So the near-cusp-0 tachyon-free region of the six golden-type models is, to 0.3-0.5%, the Fricke mirror of the
large-volume runaway, and holds no critical point: Pass 11141's conclusion (no symmetric or other stabilisation of the
Wilson-line modulus) extends to it.  Scope: one model, B = 0, two radii.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11106_one_loop_orbifold_point as P6  # noqa: E402

DATA = ROOT / "data" / "w33_pass11150_cusp0_beta_integrands.json"
OUT = ROOT / "data" / "w33_pass11150_cusp0_potential.json"


def summarize():
    d = json.loads(DATA.read_text())
    lam = {}
    for r in d['runs']:
        s, c, t = P6.fitted(d['pts'], [v[0] for v in r['vals']])
        lam[round(r['T'][1], 4)] = -0.5 * (s + c + t)
    pairs = {'0.2': (lam[0.2], lam[1.6667]), '0.17': (lam[0.17], lam[1.9608])}
    res = dict(pass_id=11150, model=d['model'], lambda_beta={str(k): v for k, v in lam.items()},
               fricke_relative_difference={k: abs(a - b) / b for k, (a, b) in pairs.items()},
               slope_toward_disk=(lam[0.2] - lam[0.17]) / 0.03,
               rolls_into_disk=lam[0.2] < lam[0.17])
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
