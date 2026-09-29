#!/usr/bin/env python3
"""Passes 11132 and 11133: the conjugate-Higgs (single light Higgs) down and lepton textures in all 21 neutral-exit
survivors -- the bottom is lowered to eps^2 in exactly the 14 with unlocked down entries, the tau can be put at the same
order in all 14, and no choice of the two VEV scales splits the light generations: the conjugate sector is anarchic.

Scan: analysis/w33_pass11132_scan_conjugate_textures.py.  For each tree-level up doublet Hu_k (the candidate light
Higgs h), its conjugate Hd (U(1) charges and space-group classes negated -- a bijection in 21/21), and the down (3 x 7)
and lepton (8 x 3) order matrices of Q d^c h*, L e^c h*: exact cubic rules, hidden-gauge check, then (a) singlets only,
(b) singlets + condensate with a weighted MILP, T insertion costing r = log eps_T / log eps_S (r = 0.5, 1, 2, 3).
Frozen: data/w33_pass11132_conjugate_textures_21.json.

11132 results:
  * the conjugation map (tree-level up doublets -> Hd) is a bijection onto three distinct doublets in 21/21;
  * without the condensate the bottom is at eps^3 in 21/21;
  * with it (r = 1) the bottom drops to eps^2 in EXACTLY the 14 models with unlocked down entries (Pass 11131);
  * in all 14 at least two of the three candidate Higgs (all three in 24165) put the tau at eps^2 too, so
    m_b / m_tau = O(1) as observed; otherwise the tau sits at eps (m_b / m_tau ~ eps).
11133 results (two scales):
  * the lowered entries use exactly one condensate insertion and one singlet (n_S, n_T) = (1, 1), or three singlets
    and no condensate; every entry of the 3 x 3 block has the SAME options, so for every r the three singular values
    share one exponent (min(1 + r, 3)): 0 cases with three distinct exponents over 21 models x 3 Higgs x 2 sectors x 4 r;
  * the conjugate sector is therefore ANARCHIC -- it can set the third-generation scale but not the light-generation
    hierarchy; with the light-doublet sector's [0,1,1] (Pass 11127), neither sector produces m_d << m_s << m_b.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "w33_pass11132_conjugate_textures_21.json"
UNL = ROOT / "data" / "w33_pass11131_unlock_scan_21.json"
OUT = ROOT / "data" / "w33_pass11132_bottom_tau_across_21.json"
RS = ('0.5', '1.0', '2.0', '3.0')


def summarize():
    d = json.loads(DATA.read_text())
    unl = json.loads(UNL.read_text())
    down_unlocked = sorted(m for m, v in unl.items() if v['results']['down']['newly_allowed'] > 0)
    lowered, tau_same, distinct, bijection, bottom3 = [], {}, 0, True, True
    for m, v in d.items():
        rows = v['rows']
        bijection &= len(rows) == 3 and all('conj_Hd' in r for r in rows.values()) \
            and len({r['conj_Hd'] for r in rows.values()}) == 3
        bottom3 &= all(r['down']['exp0'] == [3.0, 3.0, 3.0] for r in rows.values())
        if all(r['down']['expT']['1.0'] == [2.0, 2.0, 2.0] for r in rows.values()):
            lowered.append(m)
            tau_same[m] = sum(1 for r in rows.values() if r['lepton']['expT']['1.0'][0] == 2.0)
        for r in rows.values():
            for s in ('down', 'lepton'):
                for rr in RS:
                    e = r[s]['expT'][rr]
                    distinct += bool(e) and len(set(e)) == 3
    anarchic_down = all(len(set(r['down']['expT'][rr])) == 1 for v in d.values() for r in v['rows'].values() for rr in RS)
    res = dict(pass_id=[11132, 11133], n_models=len(d), conjugation_bijection=bijection, bottom_eps3_without_T=bottom3,
               lowered_models=sorted(lowered), lowered_equals_down_unlocked=sorted(lowered) == down_unlocked,
               tau_at_eps2_higgs_count=tau_same, three_distinct_exponent_cases=distinct,
               down_anarchic_every_r=anarchic_down)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
