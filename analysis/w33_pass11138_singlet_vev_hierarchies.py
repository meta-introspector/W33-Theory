#!/usr/bin/env python3
"""Pass 11138: no pattern of singlet VEVs splits the light generations -- the conjugate-sector degeneracy is structural.

Scan: analysis/w33_pass11138_scan_singlet_vevs.py.  For three models with the bottom at eps^2 (24165, 40521, 2233 =
model 57) and each tree-level top doublet h, the down (Q d^c h*) and lepton (L e^c h*) textures of the single-Higgs
sector (Pass 11130), with the condensate, but now every singlet carries its OWN scale: entry cost = sum_s w_s n_s + n_T,
w_s drawn uniformly in [1, 3] independently for every distinct singlet charge vector (12 draws), minimised by MILP;
singular-value exponents as before.  Frozen: data/w33_pass11138_singlet_vev_draws.json.
Results (3 models x 3 Higgs x 12 draws):
  * the three down singular values stay degenerate in every draw: maximal exponent spread 0.12 (a factor eps^0.12,
    about 1.3 at eps = 0.1) -- a hierarchy needs a spread of order 1;
  * the two lighter charged leptons stay degenerate (spread <= 0.10); the tau separates (by >= 1.1) only for one of the
    three Higgs choices, as in Pass 11130;
  * so the degeneracy of the conjugate sector (Pass 11133's "anarchy") is not an artefact of equal VEVs: it survives every
    singlet VEV pattern sampled.  The generation label (the fixed-point class of the SM fields) is compensated by the
    same singlet monomials for all three generations.
Scope: 12 random weight vectors per Higgs (not exhaustive); exponents only.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "w33_pass11138_singlet_vev_draws.json"
OUT = ROOT / "data" / "w33_pass11138_singlet_vev_hierarchies.json"


def summarize():
    d = json.loads(DATA.read_text())
    down_spread, light_lepton_spread, draws = 0.0, 0.0, 0
    tau_split = {}
    for m, v in d.items():
        for k, r in v.items():
            ds = [max(e) - min(e) for e in r['down'] if e]
            ls = [max(e[1:]) - min(e[1:]) for e in r['lepton'] if e]
            draws += len(ds)
            down_spread = max(down_spread, max(ds)); light_lepton_spread = max(light_lepton_spread, max(ls))
            tau_split[f"{m}|{k}"] = min(e[1] - e[0] for e in r['lepton'] if e)
    res = dict(pass_id=11138, draws=draws, max_down_spread=down_spread, max_light_lepton_spread=light_lepton_spread,
               tau_split_min=tau_split, hierarchy_found=down_spread >= 0.5 or light_lepton_spread >= 0.5)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
