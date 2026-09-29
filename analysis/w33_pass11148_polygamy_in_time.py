#!/usr/bin/env python3
"""Pass 11148: entanglement across time is polygamous -- a qutrit is maximally entangled with its past AND its future,
which no arrangement of three qutrits in space can approach.

Scan: analysis/w33_pass11148_scan_polygamy.py.  Frozen: data/w33_pass11148_polygamy.json.
Temporal: a qutrit evolving by the identity (or any Clifford gate) from a maximally mixed start: the two-time
pseudo-density operator of EVERY pair of times is SWAP/3 (Pass 11143's temporal Bell line), temporal negativity 1.  So
for three times t1 < t2 < t3: N(t1,t2) + N(t2,t3) = 2 and the sum over all three pairs = 3.
Spatial: maximise the same sums of two-qutrit negativities over three-qutrit states (negativity is convex, so pure states
suffice; 16 Nelder-Mead + Powell restarts each):
  * N_AB + N_AC <= 1.0327955589886 (numerical optimum) = 4/sqrt15 to 1e-15 (an integer-relation search gives
    15x^3 - 15x^2 - 16x + 16 = (x - 1)(15x^2 - 16): an identification, not a proof) -- versus 2 in time;
  * N_AB + N_AC + N_BC <= 1.1394 -- versus 3 in time;
  * references: the totally antisymmetric state gives 1/3 per pair; GHZ gives 0 per pair.
Reading: in space a qutrit cannot share its maximal entanglement (monogamy); across time it can, with its whole history --
the temporal Bell lines Delta_12, Delta_23, Delta_13 hold simultaneously.  This is the precise sense in which
self-entanglement across time is not a resource that can be split between partners: it is persistence, not a shared
bond.  (Temporal non-monogamy of pseudo-density operators is known in general; the qutrit numbers and the W(3,3)
reading are specific.)
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "w33_pass11148_polygamy.json"
OUT = ROOT / "data" / "w33_pass11148_polygamy_in_time.json"


def summarize():
    d = json.loads(DATA.read_text())
    sp = d['spatial_max']
    res = dict(pass_id=11148, spatial_AB_AC=sp['((0, 1), (0, 2))'], spatial_all_pairs=sp['((0, 1), (0, 2), (1, 2))'],
               temporal_AB_AC=d['temporal_AB_AC'], temporal_all_pairs=d['temporal_all'], references=d['refs'],
               polygamy_gap_two_pairs=d['temporal_AB_AC'] - sp['((0, 1), (0, 2))'],
               polygamy_gap_all_pairs=d['temporal_all'] - sp['((0, 1), (0, 2), (1, 2))'],
               ab_ac_equals_4_over_sqrt15=abs(sp['((0, 1), (0, 2))'] - 4 / 15 ** 0.5) < 1e-12)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
