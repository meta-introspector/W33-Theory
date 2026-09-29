#!/usr/bin/env python3
"""Pass 11136: the tachyon domains in the Wilson-line modulus plane are exact hyperbolic disks, centred on the special
points of Gamma_0(3) -- and the neutral disk always lies above the charged one.

Scan: analysis/w33_pass11136_scan_domains.py (explicit enumerator; both Wilson-line tori at T = B + i y, family at rho;
22 points per model incl. B = 1.0, 1.5 for periodicity).  Frozen: data/w33_pass11136_domain_samples.json.
Closed forms (Delta = -1/2 + p_R^2 / 2; b = B - nearest integer, c = B - floor(B) - 1/2), verified at every sample to 1e-9:
  neutral, sqrt3 models:   p_R^2 = (b^2 + y^2)/(sqrt3 y)                      (unit winding, A2 torus)
  neutral, golden models:  p_R^2 = min[(b^2 + y^2)/(sqrt3 y) + 1/(3 sqrt3 y),  charged form]
  charged (both):          p_R^2 = (c^2 + y^2 + 1/12)/(sqrt3 y) + 1/3
Tachyonic <=> p_R^2 < 1, i.e. Euclidean disks, period 1 in B:
  * sqrt3 neutral:  b^2 + (y - sqrt3/2)^2 < 3/4 -- tangent to the real axis at the cusp 0: the HORODISK
    Im(-1/T) > 1/sqrt3, the S-dual of the large-volume region (the T-dual reading of Pass 11128, now exact in B);
  * golden neutral: b^2 + (y - sqrt3/2)^2 < 5/12 -- the hyperbolic disk centred at i/sqrt3, the FRICKE point of level 3
    (fixed by T -> -1/(3T)), with hyperbolic radius R = 2 ln(phi): e^R = phi^2 is where the golden ratio comes from;
  * charged:        c^2 + (y - 1/sqrt3)^2 < 1/4 -- the hyperbolic disk centred at (3 + i sqrt3)/6, the order-3 ELLIPTIC
    point of Gamma_0(3), with radius R = ln(2 + sqrt3);
    in the golden models the charged-disk level also contains SM-neutral states.
  * Consequence: at every B the neutral boundary lies above the charged one: the minimum gap is at B = 1/2,
    y_N = sqrt3/2 + sqrt(5/12 - 1/4) = 1.2743 (golden; 1.5 for sqrt3) vs y_C = 1/sqrt3 + 1/2 = 1.0774 -- a descending
    Wilson-line modulus meets the SM-neutral instability first on the whole B circle, exactly.
Scope: the two models sampled (the charged form coincides in both); single-winding states (higher windings matter only
near the real axis); equal T on both Wilson-line tori.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "w33_pass11136_domain_samples.json"
OUT = ROOT / "data" / "w33_pass11136_tachyon_domains.json"
S3 = math.sqrt(3)
PHI = (1 + math.sqrt(5)) / 2


def neutral(B, y, golden):
    b = B - round(B)
    v = (b * b + y * y) / (S3 * y)
    return min(v + 1 / (3 * S3 * y), charged(B, y)) if golden else v


def charged(B, y):
    c = B - math.floor(B) - 0.5
    return (c * c + y * y + 1 / 12) / (S3 * y) + 1 / 3


def summarize():
    d = json.loads(DATA.read_text())
    dev, absent_ok = {}, True
    for m, rows in d.items():
        g = '40521' in m
        dn = max(abs(2 * (r['neutral'] + 0.5) - neutral(r['B'], r['y'], g)) for r in rows)
        dc = max((abs(2 * (r['charged'] + 0.5) - charged(r['B'], r['y'])) for r in rows if 'charged' in r), default=0)
        absent_ok &= all(charged(r['B'], r['y']) >= 1 - 1e-9 for r in rows if 'charged' not in r)
        dev[m] = dict(neutral=dn, charged=dc)
    yN_half = S3 / 2 + math.sqrt(5 / 12 - 1 / 4)
    yC_top = 1 / S3 + 0.5
    res = dict(pass_id=11136, max_deviation=dev, charged_absent_consistent=absent_ok,
               golden_top=S3 / 2 + math.sqrt(5 / 12), golden_top_is_phi2_over_sqrt3=abs(S3 / 2 + math.sqrt(5 / 12) - PHI ** 2 / S3) < 1e-12,
               golden_hyperbolic_centre=math.sqrt(3 / 4 - 5 / 12), fricke_point=1 / S3,
               golden_eR=(S3 / 2 + math.sqrt(5 / 12)) * S3, phi_squared=PHI ** 2,
               charged_hyperbolic_centre=math.sqrt(1 / 3 - 1 / 4), elliptic_point_im=S3 / 6,
               charged_eR=yC_top / math.sqrt(1 / 3 - 1 / 4), two_plus_sqrt3=2 + S3,
               sqrt3_horodisk_tangent=abs(S3 / 2 - math.sqrt(3 / 4)) < 1e-15,
               min_gap_neutral_over_charged=yN_half - yC_top)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
