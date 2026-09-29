#!/usr/bin/env python3
"""Pass 11130: with ONE light Higgs -- the doublet that gives the top -- the bottom and tau come from its conjugate, and
there the condensate lowers m_b / m_t from eps^3 to eps^2 (4 of the six); b and tau then sit at the same order.

Data:
  * data/w33_pass11130_higgs_conjugation_six.json (scan analysis/w33_pass11130_scan_higgs_conjugation.py): in all six
    the doublets with a tree-level up coupling are Hu6, Hu7, Hu8 (theta^2-twisted, one per fixed point of the family
    torus); their complex conjugates are EXACTLY Hd3, Hd5, Hd4 (U(1) charges and space-group classes negated) -- not the
    tree-level down doublets Hd6-8, which are the conjugates of Hu0-2;
  * data/w33_pass11127_textures_six.json: the down and lepton order matrices of every Hd with and without the condensate.
In a non-supersymmetric model with a single light doublet h = Hu_k, the bottom and tau Yukawas are Q d^c h* and
L e^c h*, i.e. the Hd = conj(Hu_k) textures.  Singular-value exponents (eps_VEV = eps_geo = eps, one condensate insertion
counted as one power of eps):

  model          | down (no T -> T)      | lepton (no T -> T), conj of Hu6 / Hu7 / Hu8
  36621, 46043   | [3,3,3] -> [3,3,3]    | [1,3,3] / [2,2,3] / [2,2,3]   (unchanged)
  24165          | [3,3,3] -> [2,2,2]    | [2,3,3] -> [2,2,2] (all three)
  40521,5904,17224| [3,3,3] -> [2,2,2]   | [1,3,3] -> [1,2,2] (Hu6); [2,3,3] -> [2,2,2] (Hu7, Hu8)

Estimate against the SM running masses at 2e16 GeV (Xing-Zhang-Zhou, PRD 77 (2008) 113016: m_b/m_t ~ 0.0135,
m_tau/m_t ~ 0.0228): with the condensate, n_b = n_tau = 2 needs eps_b ~ 0.116 and eps_tau ~ 0.151 -- the SAME order, so
m_b / m_tau = O(1) as observed (0.59), with a combined VEV insertion <T><S>/M^2 ~ 0.014-0.023; without the condensate
n_b = 3 would need eps ~ 0.24 with the tau at eps^2 (a factor ~4 mismatch between the sectors' natural eps).
What does NOT work: all three down (and two lepton) exponents are equal -- m_s ~ m_d ~ m_b at this order; the light
generations are not reproduced (the Pass 11102 degeneracy recurs in the conjugate sector).
Scope: exponents only (random O(1) coefficients; one common eps for S and T insertions); the size of <T> is not
computed (Pass 11123: beyond perturbation theory); the SM-at-2e16 values are an order-of-magnitude reference.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONJ = ROOT / "data" / "w33_pass11130_higgs_conjugation_six.json"
TEX = ROOT / "data" / "w33_pass11127_textures_six.json"
OUT = ROOT / "data" / "w33_pass11130_single_higgs_bottom_tau.json"
RATIO_B, RATIO_TAU = 0.0135, 0.0228            # m_b/m_t, m_tau/m_t in the SM at 2e16 GeV (Xing-Zhang-Zhou 2008)


def summarize():
    conj = json.loads(CONJ.read_text())
    tex = json.loads(TEX.read_text())
    res = dict(pass_id=11130, models={})
    for m, c in conj.items():
        hig = {r['k']: r for r in tex[m]['higgs']}
        rows = {}
        for hu, hd in c['conj_of_tree_neg'].items():
            assert len(hd) == 1
            r = hig[hd[0]]
            rows[hu] = dict(conj_Hd=hd[0], down_0=r['down_0']['exp'], down_T=r['down_T']['exp'],
                            lepton_0=r['lepton_0']['exp'], lepton_T=r['lepton_T']['exp'])
        res['models'][m] = dict(tree_up=c['tree_up_Hu'], rows=rows,
                                conj_not_tree_down=all(r['conj_Hd'] not in (6, 7, 8) for r in rows.values()))
    lowered = sorted(m for m, v in res['models'].items()
                     if all(r['down_0'] == [3.0, 3.0, 3.0] and r['down_T'] == [2.0, 2.0, 2.0] for r in v['rows'].values()))
    res.update(tree_up_is_678=all(v['tree_up'] == [6, 7, 8] for v in res['models'].values()),
               conjugates_are_heavy_doublets=all(v['conj_not_tree_down'] for v in res['models'].values()),
               bottom_lowered_3_to_2=lowered,
               eps_needed=dict(b_n2=RATIO_B ** 0.5, tau_n2=RATIO_TAU ** 0.5, b_n3=RATIO_B ** (1 / 3), tau_n1=RATIO_TAU),
               light_generations_degenerate=all(r['down_T'][1] == r['down_T'][2] for v in res['models'].values()
                                                for r in v['rows'].values()))
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    r = summarize()
    print(json.dumps({k: v for k, v in r.items() if k != 'models'}, indent=1))
    for m, v in r['models'].items():
        print(m, {k: (x['conj_Hd'], x['down_T'], x['lepton_0'], x['lepton_T']) for k, x in v['rows'].items()})
