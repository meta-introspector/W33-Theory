#!/usr/bin/env python3
"""Pass 11145: TESTED AND REFUTED -- the condensate does not make the neutrino Dirac texture anarchic.  Instead the
light-Higgs Dirac texture copies the up-quark one at tree level (an SO(10)-like m_D = m_u pattern), and the anarchic
charged-lepton sector that would give large lepton mixing gives large quark mixing too.

Hypothesis (new direction): the conjugate sector is anarchic (Passes 11133, 11138); anarchy is a known, viable
explanation of large leptonic mixing (Hall-Murayama-Weiner); if the condensate also made the neutrino Dirac matrix
anarchic, the class would predict large PMNS and (with hierarchical up quarks) small CKM.
Scan: analysis/w33_pass11145_scan_neutrino_dirac.py -- L N^c h for the three tree-level top doublets h, all hidden-
singlet SM-neutral fermions N (39-51 per model), the Pass 11102 engine with and without the condensate; the charged-lepton
texture of conj(h) from Pass 11142.  Frozen: data/w33_pass11145_neutrino_textures.json.
Results (36621, 24165, 40521; singular-value exponents on the twisted lepton doublets):
  * neutrino Dirac: [0, 1, 1] at TREE level for every light Higgs -- one same-fixed-point entry (a top-like Dirac mass)
    and a degenerate pair -- identical with and without the condensate: the unlocked nu-Dirac entries (Passes 11124,
    11131) sit below existing tree-level ones.  The pattern is the up-quark one (Pass 11102): m_D(nu) ~ m_u in texture,
    as in SO(10);
  * charged leptons (conjugate sector): [3,3,3] (36621) or [2,2,2] (24165, 40521) -- anarchic;
  * so the hypothesis fails: the neutrino side is not anarchic but degenerate like the up quarks; lepton mixing would be
    large through the anarchic charged-lepton rotation, but the down quarks come from the same conjugate sector, so quark
    mixing would be large as well -- in conflict with the small CKM angles.
Scope: exponents only; the Majorana (seesaw) sector is not computed; with a top-like Dirac mass, the light neutrino scale
needs M_R ~ m_t^2 / 0.05 eV ~ 6e14 GeV.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "w33_pass11145_neutrino_textures.json"
OUT = ROOT / "data" / "w33_pass11145_neutrino_anarchy.json"


def summarize():
    d = json.loads(DATA.read_text())
    tree_up_like, unchanged, lepton_anarchic = True, True, True
    for m, v in d.items():
        for k, r in v['rows'].items():
            tree_up_like &= r['nuT']['exp'][:3] == [0.0, 1.0, 1.0] and r['nuT']['min_order'] == 0
            unchanged &= r['nu0']['exp'] == r['nuT']['exp']
            e = r.get('charged_lepton_conj')
            lepton_anarchic &= e is not None and len(set(e)) == 1
    res = dict(pass_id=11145, models=sorted(d), neutrino_dirac_tree_up_like=tree_up_like,
               condensate_leaves_neutrino_texture=unchanged, charged_leptons_anarchic=lepton_anarchic,
               hypothesis_neutrino_anarchy='refuted', seesaw_scale_GeV=173.0 ** 2 / 0.05e-9)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
