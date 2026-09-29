#!/usr/bin/env python3
"""Pass 11131: the condensate's selection-rule effect across all 21 neutral-exit survivors (Pass 11119) -- never an
up-type Yukawa or a mu term; down/lepton entries in most; one model unlocks nothing.

Scan: analysis/w33_pass11131_scan_unlock_21.py -- as Pass 11124 (tachyon winding class in its torus's space-group class
coordinate), but with the neutral state selected independently (lowest-Delta SM-neutral tachyon on its own torus, from
the Pass 11119 torus-resolved classification) and the hidden-gauge check enforced.  Frozen:
data/w33_pass11131_unlock_scan_21.json.

Results:
  * up-type Yukawas and the mu term: 0 unlocked in 21/21;
  * lepton entries unlocked in 17/21, down-quark entries in 14/21, neutrino Dirac in 20/21 (A8SM_20260982_47048 unlocks
    nothing at all);
  * the hidden-gauge check removes no entry in any model; the positive control (a 6 q_T probe: forbidden without T,
    order 6 with it) passes in 21/21;
  * the six of Pass 11124 reproduce its counts exactly, with the neutral state chosen independently.
Reading: "up and mu untouched" is class-wide; the down-sector lowering that matters for a single light Higgs (Pass 11130)
is possible in the 14 models with unlocked down entries, verified only in the four of the six.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "w33_pass11131_unlock_scan_21.json"
PREV = ROOT / "data" / "w33_pass11124_unlock_spacegroup.json"
OUT = ROOT / "data" / "w33_pass11131_unlocked_across_21.json"
KINDS = ('up', 'down', 'lepton', 'nu_dirac', 'mu_HuHd')


def summarize():
    d = json.loads(DATA.read_text())
    prev = json.loads(PREV.read_text())
    cnt = {k: sum(1 for v in d.values() if v['results'][k]['newly_allowed'] > 0) for k in KINDS}
    res = dict(pass_id=11131, n_models=len(d), models_with_unlocked=cnt,
               counts={m: {k: v['results'][k]['newly_allowed_hidden_ok'] for k in KINDS} for m, v in d.items()},
               hidden_check_removes_none=all(v['results'][k]['newly_allowed'] == v['results'][k]['newly_allowed_hidden_ok']
                                             for v in d.values() for k in KINDS),
               control_all=all(list(v['control']) == [None, 6] for v in d.values()),
               nothing_unlocked=sorted(m for m, v in d.items() if all(v['results'][k]['newly_allowed'] == 0 for k in KINDS)),
               reproduces_11124=all(d[m]['results'][k]['newly_allowed'] == prev[m]['results'][k]['newly_allowed']
                                    for m in prev for k in KINDS))
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    r = summarize()
    print(json.dumps({k: v for k, v in r.items() if k != 'counts'}, indent=1))
