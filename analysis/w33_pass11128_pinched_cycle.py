#!/usr/bin/env python3
"""Pass 11128: shrinking the neutral Wilson-line torus to zero -- the tachyon deepens linearly to the -1/2 ground state
of a T-dual tachyonic string; there is no tachyon-free island and no new orbifold vacuum at small radius.

Scan: analysis/w33_pass11128_scan_small_radius.py (explicit Witten-twisted enumerator of Pass 11107 with winding ranges
scaled as 2.5/y; the neutral torus at i y, y = 1.6 ... 0.2, the other Wilson-line torus at 3i, the family torus at rho).
Frozen: data/w33_pass11128_small_radius_six.json.

Results (all six; exact to 1e-10):
  * the lowest tachyon is SM-neutral at every radius down to y = 0.2 (no charged tachyon ever appears on this line);
  * its level is LINEAR in the radius: Delta = -1/2 + p_R^2/2 with
        p_R^2 = y1/sqrt3                          (36621, 17224),
        p_R^2 = y1/sqrt3 + 1/(3 sqrt3 y2)          (46043, 24165, 40521, 5904),
    i.e. a pure unit winding on the neutral torus (A2 lattice, B = 0), plus -- in the golden models -- a Wilson-line-
    shifted momentum 1/3 on the other torus; level matching then forces l^2 = 1;
  * on the diagonal y1 = y2 = y the golden onset is the root of y^2 - sqrt3 y + 1/3 = 0:
        y_c = (sqrt3 + sqrt(5/3))/2 = (3 + sqrt5)/(2 sqrt3) = phi^2/sqrt3  (Pass 11123's golden radius, now derived);
  * as y -> 0, Delta -> -1/2 and the number of tachyonic states grows (4 for y >= 0.5, 6 for 0.25-0.4, 8 at 0.2):
    the winding tower becomes the momentum tower of a decompactifying T-dual circle, carrying the Delta = -1/2 ground
    state of a tachyonic non-supersymmetric heterotic string (the SO(16)xSO(16) string on a small circle with Wilson
    lines is T-dual to the tachyonic ten-dimensional strings; Fraiman-Grana-Parra De Freitas-Sethi, arXiv:2307.13745).
Reading: the pinched-cycle endpoint is not a nearby tachyon-free orbifold. Along the whole line the instability only
deepens; no radius below y_c stops it.  Scope: B = 0, one torus shrunk at a time, the other at 3i; the gauge group of the
T-dual string is not identified here (its tachyon has l^2 = 1, the vector class of the SO(32) / SO(16)xE8 tachyonic
strings).
"""
from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "w33_pass11128_small_radius_six.json"
OUT = ROOT / "data" / "w33_pass11128_pinched_cycle.json"
S3 = math.sqrt(3)
PHI = (1 + math.sqrt(5)) / 2


def summarize():
    d = json.loads(DATA.read_text())
    fits, counts, neutral_only = {}, {}, True
    for m, v in d.items():
        shifts = set()
        for y, c in v.items():
            neutral_only &= set(c) == {'neutral'}
            shifts.add(round(c['neutral'][1] - (-0.5 + float(y) / (2 * S3)), 9))
            counts.setdefault(m, {})[y] = c['neutral'][0]
        fits[m] = sorted(shifts)
    golden = sorted(m for m, s in fits.items() if s == [round(1 / (18 * S3), 9)])
    root = (S3 + math.sqrt(5 / 3)) / 2
    res = dict(pass_id=11128, neutral_only_down_to_0p2=neutral_only, linear_law_shift=fits, golden_models=golden,
               sqrt3_models=sorted(m for m, s in fits.items() if s == [0.0]), tachyon_counts=counts,
               golden_root=root, golden_root_is_phi2_over_sqrt3=abs(root - PHI ** 2 / S3) < 1e-12,
               quadratic_residual=root ** 2 - S3 * root + 1 / 3,
               deepest=min(c['neutral'][1] for v in d.values() for c in v.values()))
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
