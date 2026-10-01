#!/usr/bin/env python3
"""Pass 11247: what the frozen data can and cannot say about R-symmetries of the Z6-I models (the W33 vacuum's class).

Pass 11241 decides mu protection by the full vacuum symmetry that the frozen charges encode.  For Z6-II the charges of
Pass 10967 include the geometric R-symmetries Z6^R x Z3^R x Z2^R (W charge -1 each); for Z6-I they do not.  This pass
certifies that gap from the data and records exactly which inputs would close it:
  * per field, the shifted H-momentum q_sh (three components) and the oscillator numbers N^i, Nbar^i in each twisted
    sector -- the R-charge R^i = q_sh^i - N^i + Nbar^i (mod the order of the plane rotation);
  * for Z6-I, twist v = (1/6, 1/6, -1/3): rotations of orders (6, 6, 3) give candidate Z6^R x Z6^R x Z3^R;
these come from the orbifolder spectrum files (raw .sp dumps), which are not in this repository.
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DISC = ROOT / "data" / "w33_pass10967_space_group_discrete_charges.json.gz"
OUT = ROOT / "data" / "w33_pass11247_z6i_r_symmetry_data.json"


def run():
    d = json.load(gzip.open(DISC, "rt"))
    z6i = sorted(k for k in d if k.startswith("Z6I_"))
    z6ii = sorted(k for k in d if k.startswith("Z6II_"))
    return dict(pass_id=11247, z6i_models=len(z6i), z6i_with_R=sum(1 for k in z6i if d[k]["R"]),
                z6i_nonR_orders=sorted({tuple(d[k]["nonR_orders"]) for k in z6i}),
                z6ii_models=len(z6ii), z6ii_with_R=sum(1 for k in z6ii if d[k]["R"]),
                z6ii_R=sorted({tuple(d[k]["R"]) for k in z6ii}),
                required=["q_sh per field per twisted sector", "oscillator numbers N^i, Nbar^i",
                          "Z6-I plane-rotation orders (6, 6, 3)"],
                status="BLOCKED: Z6-I R charges are not in the repository; Pass 11241's Z6-I verdict covers U(1)s, "
                       "their discrete remnants and the point-group Z6 only")


def main():
    res = run()
    OUT.write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
