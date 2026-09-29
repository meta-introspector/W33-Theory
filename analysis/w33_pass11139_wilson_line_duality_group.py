#!/usr/bin/env python3
"""Pass 11139: the duality group of the Wilson-line modulus is Gamma_0(3) -- Fricke and S are broken; the exact group
explains the Gamma_0(3) points of the tachyon disks (Pass 11136).

Test (analysis/w33_pass11139_scan_duality_group.py): the Witten-twisted (beta) partition function Z_beta(tau; T) of the
full engine (both Wilson-line tori at T, family torus at rho, K = 3 and K = 4 Wilson-line cosets) compared with
Z_beta(tau; gamma T) at 3 values of tau and 1-2 points T per map, for one sqrt3 model (36621) and one golden model (40521).
Frozen: data/w33_pass11139_duality_group_tests.json.
An exact symmetry shows a deviation that falls with the coset cutoff K; a broken one a K-independent deviation.
Results (max relative deviation, K = 3 -> K = 4):
  * T -> T + 1: 0 (exact);  T -> T + 3: 1.8e-7 at both K (floating-point, large Re T);
  * T -> T/(3T + 1) (generator of Gamma_0(3)): 5e-6 -> 7e-10 (36621), 3e-6 -> 1e-9 (40521): EXACT;
  * T -> (-T + 1)/(-3T + 2) (the order-3 elliptic element of Gamma_0(3)): 3e-4 -> 2.7e-7: EXACT;
  * S: T -> -1/T: 3% (broken);  T -> T/(T + 1) (Gamma^0(3)): 18-19% (broken);
  * Fricke T -> -1/(3T): 0.34% (36621), 0.043% (40521), K-independent: APPROXIMATE only (cf. Pass 11115's approximate
    Fricke duality in model 10).
So the duality group is Gamma_0(3) (the level-3 congruence subgroup expected for Z3 Wilson lines), not SL(2,Z) and
not its Fricke extension.  Its only elliptic point (order 3), (3 + i sqrt3)/6, is EXACTLY the centre of the charged
tachyon disk of Pass 11136; the golden neutral disk is centred at the Fricke point i/sqrt3, which is a symmetry point
only approximately.
Scope: two models, the beta sector at fixed tau (the moduli-independent twisted sectors are invariant trivially).
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "w33_pass11139_duality_group_tests.json"
OUT = ROOT / "data" / "w33_pass11139_wilson_line_duality_group.json"


def summarize():
    d = json.loads(DATA.read_text())
    table = {}
    for lab, maps in d.items():
        for name, e in maps.items():
            dev = {3: 0.0, 4: 0.0}
            for key, v in e.items():
                K = int(key.split('|')[2][1:])
                zT = complex(*v['T']); zg = complex(*v['gT'])
                dev[K] = max(dev[K], abs(zg / zT - 1))
            table.setdefault(name, {})[lab] = dev
    exact = sorted(n for n, v in table.items() if all(x[4] < 1e-6 for x in v.values()))
    broken = sorted(n for n, v in table.items() if all(x[4] > 1e-2 for x in v.values()))
    approx = sorted(n for n in table if n not in exact and n not in broken)
    res = dict(pass_id=11139, deviations={n: {m: [v[3], v[4]] for m, v in x.items()} for n, x in table.items()},
               exact=exact, broken=broken, approximate=approx, group='Gamma_0(3)')
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    r = summarize()
    print(r['exact'], r['broken'], r['approximate'])
