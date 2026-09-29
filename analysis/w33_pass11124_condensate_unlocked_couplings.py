#!/usr/bin/env python3
"""Pass 11124: what the SM-neutral tachyon condensate unlocks in the six neutral-exit survivors -- charged-lepton,
down-quark and neutrino-Dirac Yukawas forbidden at every order become allowed; up quarks and mu are untouched.

Scan: analysis/w33_pass11124_scan_dumps.py (orbifolder dumps of Pass 11108; FastOrders MILP of the selection-rule engine
with the hidden-singlet SM-neutral scalars as VEV fields, then again with the tachyon T and its conjugate added).
Frozen outputs:
  data/w33_pass11124_unlock_spacegroup.json    T's winding class in its torus's space-group class coordinate (correct);
  data/w33_pass11124_unlock_separate_winding.json  T's winding as a separate exactly-conserved Z3 (first encoding);
  data/w33_pass11124_target_classes.json       continuous-charge targets: in span(S), need q_T, or outside;
  data/w33_pass11124_rank_and_control.json     rank of the singlet charges with/without q_T; positive control.

Results:
  * q_T is NOT in the span of the singlet-VEV charges (rank 5 -> 6 or 4 -> 5 in all six): the condensate breaks a U(1)
    the fractional-charge-removing vacuum left unbroken (Pass 11117's one U(1));
  * every coupling target that lies in span(S, q_T) but not span(S) needs EXACTLY +-1 unit of q_T -- one condensate
    insertion (102 distinct targets over the six, coefficients all +-1);
  * with the space-group encoding they are all unlocked: charged-lepton Yukawas in 6/6 (27-108 entries), down-quark
    Yukawas in 4/6 (81 entries each; not in 36621, 46043), neutrino Dirac in 6/6 (54-324), at orders 2-6; up-type
    Yukawas and the mu term: 0 in all six;
  * 1854 of the 1863 unlocked entries contain a theta-twisted SM field, where the mod-(1 - theta) winding rule certainly
    applies; the other 9 (40521, nu Dirac) are conditional on a twisted singlet in the minimal solution;
  * CORRECTION within the pass: the first encoding (winding as a separate exact Z3) gave 0 unlocked couplings in all
    six -- an artefact: winding mod (1 - theta) Lambda is the same Z3 as the fixed-point class and is compensated by
    twisted fields.  It is kept as the strict-rule bound (exact for couplings with no twisted field);
  * control: a probe with charge 6 q_T is forbidden without T and allowed at order 6 with it, in all six.
Reading (CORRECTED by Passes 11127 and 11130): the unlocked entries are forbidden ENTRIES of down/lepton matrices that
already have allowed ones -- the light (tree-level) doublets H6-8 have lepton and down Yukawas at tree level, and their
hierarchy ([0,1,1]) is unchanged by the condensate.  The original reading ("a new source for exactly the Yukawas this
class lacks") was an over-read.  Where the unlocked entries matter is the conjugate of the top Higgs (Hd3-5): with one
light Higgs, m_b/m_t drops from eps^3 to eps^2 in 4 of the six (Pass 11130).  The selection-rule statements above stand
(Pass 11127 re-derived the counts with the hidden-gauge check: identical).
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "data"
OUT = D / "w33_pass11124_condensate_unlocked_couplings.json"
KINDS = ('up', 'down', 'lepton', 'nu_dirac', 'mu_HuHd')


def summarize():
    sg = json.loads((D / "w33_pass11124_unlock_spacegroup.json").read_text())
    sep = json.loads((D / "w33_pass11124_unlock_separate_winding.json").read_text())
    tc = json.loads((D / "w33_pass11124_target_classes.json").read_text())
    rk = json.loads((D / "w33_pass11124_rank_and_control.json").read_text())
    new = {m: {k: v['results'][k]['newly_allowed'] for k in KINDS} for m, v in sg.items()}
    orders = sorted({int(o) for v in sg.values() for k in KINDS for o in v['results'][k]['new_orders']})
    total = sum(sum(x.values()) for x in new.values())
    twisted = sum(v['results'][k]['new_with_theta_twisted_sm_field'] for v in sg.values() for k in KINDS)
    coeffs = [c for v in tc.values() for k in KINDS for c in v['classes'][k].get('T_coeffs', [])]
    res = dict(pass_id=11124, newly_allowed=new, orders=orders, total_new=total, new_with_twisted_sm_field=twisted,
               lowered=sum(v['results'][k]['order_lowered'] for v in sg.values() for k in KINDS),
               separate_winding_total=sum(v['results'][k]['newly_allowed'] for v in sep.values() for k in KINDS),
               qT_outside_singlet_span=all(v['rank_with_T'] == v['rank_singlets'] + 1 for v in rk.values()),
               control=all(v['control_order_without_T'] is None and v['control_order_with_T'] == 6 for v in rk.values())
               and all(v['control'][:2] == [None, 6] for v in sg.values()),
               needs_T_targets=len(coeffs), needs_T_coeffs=sorted(set(coeffs)),
               lepton_all_six=all(x['lepton'] > 0 for x in new.values()),
               down_models=sorted(m for m, x in new.items() if x['down'] > 0),
               up_and_mu_untouched=all(x['up'] == 0 and x['mu_HuHd'] == 0 for x in new.values()))
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
