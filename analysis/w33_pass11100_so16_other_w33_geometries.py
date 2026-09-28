#!/usr/bin/env python3
"""Pass 11100: SO(16)xSO(16) completions of the W(3,3) shifts on the other geometries -- only T6/Z3 works.

Pass 11095 found 104 tachyon-free three-generation SO(16)xSO(16) models built on the W(3,3) A8 shift of T6/Z3.  Here the
same construction is run on the other geometries that carry W(3,3) shifts:
  * T6/(Z3xZ3): the two distinct (V1, V2) W(3,3) bases of the supersymmetric Z3xZ3 scan (V1 = the A8 Kac pair);
  * T6/Z6-I:    the 58 W(3,3) base shifts (Holotrade 5b3f3ad scan);
  * T6/Z12-I:   the 117 distinct base shifts that gave supersymmetric Standard Models in the Z12-I scan.
For each base, a Witten shift V0 = (a; b) with a, b in (1/2)(norm-4 E8 vectors) and V0.Vk in Z for every shift is chosen
(standard-like first; the first one the non-SUSY orbifolder loads); the W(3,3) shifts are held fixed and the Wilson
lines drawn at random (analysis/orbifolder_n0_drivers/nsoscan.cpp), with the non-SUSY orbifolder's tachyon and SM tests.
Every base found a loadable V0 (2 + 58 + 117).
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EVID = ROOT / "data" / "w33_pass11100_so16_other_geometries_evidence.json"
P95 = ROOT / "data" / "w33_pass11095_so16xso16_a8_standard_models.json"
OUT = ROOT / "data" / "w33_pass11100_so16_other_w33_geometries.json"


def main():
    ev = json.loads(EVID.read_text())
    s95 = json.loads(P95.read_text())["summary"]["scan"]
    z3 = dict(tries=sum(v["tries"] for v in s95.values()), tachyon_free=sum(v["tachyon_free"] for v in s95.values()),
              sm_draws=sum(v["sm"] for v in s95.values()), sm_inequivalent=sum(v["sm_inequivalent"] for v in s95.values()),
              sm_and_tachyon_free=sum(v["sm_and_tachyon_free"] for v in s95.values()))
    table = {
        "T6/Z3 (Pass 11095)": z3,
        "T6/(Z3xZ3)": dict(bases=ev["Z3xZ3"]["distinct"]["bases"], tries=ev["Z3xZ3"]["distinct"]["tries"],
                           tachyon_free=ev["Z3xZ3"]["distinct"]["tachyon_free"], sm_draws=ev["Z3xZ3"]["distinct"]["sm_draws"],
                           sm_and_tachyon_free=ev["Z3xZ3"]["distinct"]["sm_and_tachyon_free"]),
        "T6/Z6-I": dict(bases=ev["Z6-I"]["bases"], tries=ev["Z6-I"]["totals"]["tries"], tachyon_free=ev["Z6-I"]["totals"]["tachyon_free"],
                        sm_draws=ev["Z6-I"]["totals"]["sm"], sm_and_tachyon_free=ev["Z6-I"]["totals"]["sm_and_tachyon_free"]),
        "T6/Z12-I": dict(bases=ev["Z12-I"]["bases"], tries=ev["Z12-I"]["totals"]["tries"], tachyon_free=ev["Z12-I"]["totals"]["tachyon_free"],
                         sm_draws=ev["Z12-I"]["totals"]["sm"], sm_and_tachyon_free=ev["Z12-I"]["totals"]["sm_and_tachyon_free"]),
    }
    for k, v in table.items():
        v["tachyon_free_fraction"] = round(v["tachyon_free"] / v["tries"], 4)
    res = dict(pass_id=11100, table=table,
               only_z3_gives_tachyon_free_sm=all(v["sm_and_tachyon_free"] == 0 for k, v in table.items() if "Z3 (" not in k))
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
