#!/usr/bin/env python3
"""Pass 11154: the structure of the 4/sqrt15 optimum -- the best a qutrit can share its entanglement between two partners
in space is carried by a state whose every invariant is a fifteenth.

Scan: analysis/w33_pass11154_scan_polygamy_optimum.py (30 Nelder-Mead + Powell restarts of max N_AB + N_AC over pure
three-qutrit states; the optimum is frozen, with its invariants, in data/w33_pass11154_polygamy_optimum_state.json).
Identified closed forms (all to ~1e-8, the optimiser's precision):
  * single-party spectra: A (7, 4, 4)/15; B and C (11, 2, 2)/15 (B <-> C symmetric);
  * N_AB = N_AC = 2/sqrt15 (so N_AB + N_AC = 4/sqrt15), N_BC = 1/15;
  * the partial transpose of rho_AB has spectrum (1 - sqrt15)/15 (twice), -2/15, 2/15 (three times), (1 + sqrt15)/15
    (twice), 7/15 -- the three negative eigenvalues sum to 2/sqrt15.
Compare time (Pass 11148): the same qutrit, persisting under the identity, has negativity 1 with its past AND its future.
Status: the optimum's invariants are identified exactly; that 4/sqrt15 is the global maximum is supported by 46 random
restarts (Passes 11148 and 11154) but not proved -- an analytic proof, e.g. from the rational spectra, is open.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "w33_pass11154_polygamy_optimum_state.json"
OUT = ROOT / "data" / "w33_pass11154_polygamy_optimum_structure.json"


def summarize():
    d = json.loads(DATA.read_text())
    devs = {k: v for k, v in d.items() if k.endswith('_dev')}
    res = dict(pass_id=11154, deviations=devs, all_closed_forms_within_1e_7=all(v < 1e-7 for v in devs.values()),
               spectra={'A': '(7,4,4)/15', 'B': '(11,2,2)/15', 'C': '(11,2,2)/15'},
               negativities={'AB': '2/sqrt15', 'AC': '2/sqrt15', 'BC': '1/15'})
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
