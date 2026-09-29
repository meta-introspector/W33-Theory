#!/usr/bin/env python3
"""Pass 11127: the condensate-unlocked Yukawas do not change the lepton (or down-quark) mass hierarchy of the light Higgs
-- m_mu = m_e survives, and CORRECTS the reading of Pass 11124.

Scan: analysis/w33_pass11127_scan_textures.py -- for every H_d doublet of the six, the full 8 x 3 lepton and 3 x 7 down
order matrices with the Pass 11102 engine (hidden-gauge invariance, EXACT cubic rules, H-momentum mod 3 beyond), once
with the singlet vacuum and once with the condensate T added (as in Pass 11124).  Singular-value exponents with
eps_VEV = eps_geo = eps (Pass 11102 convention, random O(1) coefficients).  Frozen: data/w33_pass11127_textures_six.json.

Results (all six):
  * the three LIGHT doublets (those with a tree-level down coupling, H6-H8 in every model) already have tree-level
    lepton Yukawas: one same-fixed-point entry (geometric factor 1) and two distinct-fixed-point entries on the twisted
    leptons -- exponents [0, 1, 1] for down AND lepton, with and without T;
  * the cubic lepton entries of a light doublet form a permutation of the three twisted generations, symmetric under the
    point reflection about the heavy generation's fixed point -- the Pass 11102 mechanism that forces m_c = m_u -- so
    m_mu = m_e at tree level, and the condensate (a Delta(54) singlet entering only at order >= 2) cannot split them at
    leading order;
  * the condensate changes singular-value exponents ONLY for the heavier doublets H3-H5 (down [3,3,3] -> [2,2,2] in four
    models; lepton [2,3,3] or [1,3,3] -> [2,2,2] or [1,2,2]); never for a light doublet;
  * with hidden-gauge invariance and exact cubic rules the entries that exist with T but not without are counted per
    model (field 'new_entries'): lepton 27/27/81/108/81/81 and down 0/0/81/81/81/81 -- IDENTICAL to Pass 11124, so the
    hidden-gauge check (omitted there) removes none: the 11124 counts stand.
CORRECTION to Pass 11124: its reading "a new source for exactly the Yukawas this class lacks (down/lepton)" was an
over-read -- the lepton and down Yukawas exist at tree level for the light doublets (Pass 11101 already gave a heavy
bottom); what the condensate unlocks are individual forbidden ENTRIES of matrices that already have allowed ones, and
they do not change the light-Higgs hierarchy.  The selection-rule statements of 11124 (one insertion; q_T outside the
singlet span; up and mu untouched) stand.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11102_higher_order_quark_hierarchy as H  # noqa: E402

DATA = ROOT / "data" / "w33_pass11127_textures_six.json"
OUT = ROOT / "data" / "w33_pass11127_lepton_texture_under_condensation.json"


def cubic_block(rec):
    """the lepton entries of order 0 (cubic), re-indexed to a 3x3 block over the rows that carry them"""
    o, g = rec['orders'], rec['geo']
    cub = {tuple(map(int, k.split(','))): v for k, v in o.items() if v == 0}
    rows = sorted({i for i, _ in cub})
    return {f"{rows.index(i)},{j}": 0 for (i, j) in cub}, {f"{rows.index(i)},{j}": g[f'{i},{j}'] for (i, j) in cub}, rows


def summarize():
    d = json.loads(DATA.read_text())
    res = dict(pass_id=11127, models={})
    all_light_unchanged, all_sym, heavy_only = True, True, True
    for m, v in d.items():
        light = [r for r in v['higgs'] if r['down_0']['min_order'] == 0]
        rec = dict(light_doublets=[r['k'] for r in light], new_entries={}, exps_light={}, reflection_symmetric={})
        for r in light:
            ex = {s: r[s]['exp'] for s in ('down_0', 'down_T', 'lepton_0', 'lepton_T')}
            rec['exps_light'][str(r['k'])] = ex
            all_light_unchanged &= ex['down_0'] == ex['down_T'] == ex['lepton_0'] == ex['lepton_T'] == [0.0, 1.0, 1.0]
            o, g, rows = cubic_block(r['lepton_T'])
            sym = H.reflection_symmetric(dict(orders=o, geo=g)) if len(rows) == 3 else None
            rec['reflection_symmetric'][str(r['k'])] = sym
            all_sym &= sym is True
        for s in ('down', 'lepton'):
            rec['new_entries'][s] = sum(len(set(r[s + '_T']['orders']) - set(r[s + '_0']['orders'])) for r in v['higgs'])
        changed = [(r['k'], s) for r in v['higgs'] for s in ('down', 'lepton') if r[s + '_0']['exp'] != r[s + '_T']['exp']]
        rec['exponent_changes'] = changed
        heavy_only &= all(k in (3, 4, 5) for k, _ in changed)
        res['models'][m] = rec
    res.update(light_exponents_unchanged_0_1_1=all_light_unchanged, light_lepton_reflection_symmetric=all_sym,
               exponent_changes_only_heavy_doublets=heavy_only,
               new_entries_total={s: sum(r['new_entries'][s] for r in res['models'].values()) for s in ('down', 'lepton')})
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    r = summarize()
    print(json.dumps({k: v for k, v in r.items() if k != 'models'}, indent=1))
    for m, v in r['models'].items():
        print(m, v['light_doublets'], v['new_entries'], v['reflection_symmetric'], v['exponent_changes'])
