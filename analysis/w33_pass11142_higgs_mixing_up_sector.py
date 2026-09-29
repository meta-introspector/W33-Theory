#!/usr/bin/env python3
"""Pass 11142: no mixture of Higgs doublets splits charm from up -- the up-sector degeneracy survives every
two-doublet light Higgs.

Data: data/w33_pass11142_all_doublet_textures.json (up order matrices Q u^c Hu_k for all nine Hu, down/lepton for all nine
Hd, with the condensate; Pass 11102 engine) and data/w33_pass11142_mixing_search.json (search).
Search: light Higgs h = Hu_k + eps^a Hu_j with Hu_k one of the three tree-level (top) doublets, Hu_j any other doublet,
a = 1, 2, 3; entrywise order min(n^k_ij, a + n^j_ij); singular-value exponents with random O(1) coefficients; the
geometric factor of cubic entries counted as eps^g with g = 1 (small) or g = 0 (O(1), as at T* = rho where it is 1/2).
Results (models 36621, 24165, 40521, 2233; 72 mixtures x 2 conventions each):
  * the up exponents are ALWAYS [0, 1, 1] (g = 1) or [0, 0, 0] (g = 0): 0 of 576 mixtures give three distinct
    exponents with the top at order 0;
  * reason: every admixture enters at order >= a >= 1 on top of an entry pattern whose charm/up pair is already set at
    tree level by the point reflection (Pass 11102); admixtures shift the pair together at leading order.
With Passes 11102 (every vacuum alignment), 11114 (every Higgs direction in the tree triplet), 11127 (the condensate) and
11133 (two VEV scales), this closes the up-sector hierarchy question at the level of selection rules: m_c = m_u at
leading order in every construction tried.  Scope: exponents only; subleading O(eps^a) splittings are not hierarchies.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "w33_pass11142_mixing_search.json"
OUT = ROOT / "data" / "w33_pass11142_higgs_mixing_up_sector.json"


def summarize():
    d = json.loads(DATA.read_text())
    total = sum(sum(v['patterns'].values()) for v in d.values())
    hits = sum(len(v['hierarchical_hits']) for v in d.values())
    pats = sorted({p.split('|')[1].replace('-0.0', '0.0') for v in d.values() for p in v['patterns']})
    res = dict(pass_id=11142, models=sorted(d), mixtures_tested=total, hierarchical_hits=hits, exponent_patterns=pats,
               tree_doublets={m: v['tree'] for m, v in d.items()})
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
