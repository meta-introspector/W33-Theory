#!/usr/bin/env python3
"""Pass 11141: every symmetry-forced critical point of the Wilson-line potential is tachyonic or a runaway -- the one-loop
potential cannot stabilise the Wilson-line modulus at a self-dual point.

Inputs: Pass 11139 (the duality group is Gamma_0(3): cusps infinity (width 1) and 0 (width 3), one elliptic point of
order 3, (3 + i sqrt3)/6; genus 0); Pass 11136 (exact tachyon disks); Pass 11122 (Lambda_beta = volume law, monotone in
Im T); Pass 11129 (B near-flat); data/w33_pass11146_no_holonomy_21.json (tachyons at small radius on the diagonal, all
21 models).
A Gamma_0(3)-invariant potential is automatically stationary at the elliptic point and extremal toward the cusps.
Results:
  * the elliptic point (3 + i sqrt3)/6 is EXACTLY the centre of the charged tachyon disk: p_R^2 = 2/3 there, Delta = -1/6
    (a charged tachyon) -- in every model;
  * the cusp at infinity is the large-volume runaway (Lambda grows as the volume, the moduli roll inward);
  * the cusp at 0, approached along the diagonal of both Wilson-line tori: tachyonic at Im T = 0.2 in 15/21 models
    (winding tachyons, count growing); in the other 6 (46043, 24165, 40521, 5904, 5285, 2080) the neutral tachyon
    becomes SHALLOWER below Im T ~ 0.6 and is gone at 0.2 -- the golden disk's lower edge (the other root of
    y^2 - sqrt3 y + 1/3 = 0, y = 0.2205): there the region near cusp 0 is tachyon-free, the approximate Fricke image
    (y <-> 1/(3y), 0.04% in 40521) of the large-volume region, and by that approximate symmetry a runaway into the disk
    from below.  (The single-torus zero-radius limit is the U(16) string, Pass 11140 -- a different direction.)
  * the Fricke point i/sqrt3 (centre of the golden neutral disk) is only an approximate symmetry point (Pass 11139) and is
    itself tachyonic;
  * the tachyon-free region of the Gamma_0(3) fundamental domain is a neighbourhood of the cusp at infinity above the
    disks, where Lambda_beta is monotone in Im T with a near-flat B direction: no critical point.
So symmetric (self-dual-point) stabilisation of the Wilson-line modulus is excluded in the neutral-exit vacua: every
point the duality group singles out is tachyonic or a runaway.  Scope: the beta-sector potential on the diagonal of the
two Wilson-line tori; the absolute Lambda includes moduli-independent sectors (positive, Pass 11106).
"""
from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "data"
OUT = D / "w33_pass11141_no_symmetric_stabilization.json"
S3 = math.sqrt(3)


def charged_pR2(B, y):
    c = B - math.floor(B) - 0.5
    return (c * c + y * y + 1 / 12) / (S3 * y) + 1 / 3


def summarize():
    grp = json.loads((D / "w33_pass11139_wilson_line_duality_group.json").read_text())
    dom = json.loads((D / "w33_pass11136_tachyon_domains.json").read_text())
    nh = json.loads((D / "w33_pass11146_no_holonomy_21.json").read_text())
    ell = (0.5, S3 / 6)
    pr2 = charged_pR2(*ell)
    free02 = sorted(m for m, v in nh.items() if v['0.2']['with_wl'] == 0)
    res = dict(pass_id=11141, group=grp['group'], elliptic_point=list(ell), elliptic_pR2=pr2, elliptic_Delta=-0.5 + pr2 / 2,
               elliptic_is_charged_disk_centre=abs(dom['charged_hyperbolic_centre'] - ell[1]) < 1e-12,
               cusp0_diagonal_tachyonic=21 - len(free02), cusp0_diagonal_tachyon_free_at_0p2=free02,
               golden_lower_root=(S3 - math.sqrt(5 / 3)) / 2,
               fricke_only_approximate='Fricke3' in grp['approximate'])
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
