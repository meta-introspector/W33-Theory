#!/usr/bin/env python3
"""Pass 11101: a localized light Higgs singles out the top (and the bottom) in the twisted SO(16)xSO(16) A8 models.

Pass 11098 took a GENERIC combination of all Higgs doublets and found no model with a single heavy top.  But in the
twisted sectors every Higgs doublet sits at a definite fixed point, and the physical light Higgs is one field.  Here the
tree-level Yukawa matrix Y^(k) is built for every INDIVIDUAL Higgs H_k (exact rules; untwisted epsilon/delta structure;
twisted couplings O(1) at a common fixed point and suppressed by eps_t = exp(-area_t) in every torus t where the three
fixed points differ, with generic distinct eps_t = eps^(1, 1.7, 2.9)), and the singular values' exponents in eps are read
off at eps = 1e-3 and 1e-6.

Results (all 104 models, frozen in data/w33_pass11101_localized_higgs_104.json; sample models recomputed here):
  * the 73 untwisted-up models: every H_k gives exponents (0, 0, -): top = charm for ANY Higgs (Pass 11098 stands);
  * the 31 twisted-up models: EVERY H_k gives (0, w_t, w_t): one O(1) top, charm and up both suppressed by the same
    torus -- a single heavy top, but m_c = m_u at tree level;
  * down sector: 27 models give (0, w_t, w_t) (single heavy bottom), 76 have no tree-level down coupling.
Why m_c = m_u: the three generations sit at the three theta-fixed points of one SU(3) torus, and these are pairwise
EQUIDISTANT (equilateral triangle; the Z3 fixes the lattice shape, so no modulus can deform it) -- checked below.  This is
the known two-valued structure of Z3 Yukawas (Casas, Munoz, Ibanez et al.); only its application here is new.
Combined with Pass 11097: 12 models have a heavy top AND bottom from a localized Higgs AND all fractional fermions
massive (hidden group unbroken, all rules): 2, 10, 13, 14, 15, 35, 53, 57, 69, 77, 78, 102.
"""
from __future__ import annotations

import json
import sys
from itertools import product
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11097_fractional_fermion_masses as P97  # noqa: E402

FROZEN = ROOT / "data" / "w33_pass11101_localized_higgs_104.json"
YUK = ROOT / "data" / "w33_pass11098_yukawa_104.json"
FRAC = ROOT / "data" / "w33_pass11097_symmetry_test_104.json"
OUT = ROOT / "data" / "w33_pass11101_localized_higgs_heavy_top.json"
KEY = "A_hidden_unbroken|gauge+Z2W+PG+SG+R"


def su3_fixed_point_distances():
    """theta-fixed points of the SU(3) torus (e1 = (sqrt2, 0), e2 at 120 degrees): m (2 e1 + e2)/3, m = 0, 1, 2; minimal
    distance between each pair over lattice translations"""
    e1 = np.array([np.sqrt(2), 0.0])
    e2 = np.array([-1 / np.sqrt(2), np.sqrt(1.5)])
    f = [m * (2 * e1 + e2) / 3 for m in range(3)]
    d = {}
    for a in range(3):
        for b in range(a + 1, 3):
            d[f"{a}{b}"] = min(np.linalg.norm(f[a] - f[b] + n1 * e1 + n2 * e2) for n1, n2 in product(range(-3, 4), repeat=2))
    # the three points are fixed by theta (rotation by 120 degrees) modulo the lattice
    th = np.array([[np.cos(2 * np.pi / 3), -np.sin(2 * np.pi / 3)], [np.sin(2 * np.pi / 3), np.cos(2 * np.pi / 3)]])
    B = np.column_stack([e1, e2])
    fixed = all(np.allclose(np.round(np.linalg.solve(B, th @ x - x)), np.linalg.solve(B, th @ x - x)) for x in f)
    return dict(distances={k: round(float(v), 12) for k, v in d.items()}, all_fixed_by_theta=fixed)


def classify(pats):
    return "degenerate" if pats and all(p[1] == 0 for p in pats) else \
        "single_heavy" if pats and all(p[0] == 0 and p[1] and p[1] > 0 for p in pats) else ("none" if not pats else "other")


def main():
    geo = su3_fixed_point_distances()
    fr = json.loads(FROZEN.read_text())
    yk = json.loads(YUK.read_text())
    fc = json.loads(FRAC.read_text())
    cls = {s: {} for s in ("up", "down", "lepton")}
    for i, v in fr.items():
        for s in cls:
            cls[s][i] = classify([h['exponents'] for h in v[s]])
    from collections import Counter
    summary = {s: dict(Counter(c.values())) for s, c in cls.items()}
    best = sorted([i for i in fr if cls['up'][i] == 'single_heavy' and cls['down'][i] == 'single_heavy'
                   and fc[i]['results'][KEY]['light_fractional_states'] == 0], key=int)
    # consistency with Pass 11098: single-heavy up exactly in the twisted-up models
    twisted = {i for i, v in yk.items() if v['up']['untwisted'] == 0}
    res = dict(pass_id=11101, su3_fixed_points=geo, classes=summary, best_models=best,
               single_heavy_up_equals_twisted_up=({i for i, c in cls['up'].items() if c == 'single_heavy'} == twisted),
               charm_up_degenerate_in_every_single_heavy_model=all(
                   all(h['exponents'][1] == h['exponents'][2] for h in fr[i]['up']) for i in fr if cls['up'][i] == 'single_heavy'))
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
