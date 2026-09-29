#!/usr/bin/env python3
"""Pass 11112: in models 53 and 57 the Standard-Model-neutral and the charged winding tachyons live on DIFFERENT
Wilson-line tori -- and in model 57 the one-loop potential shrinks the neutral torus first.

Pass 11107: in 10 of the 12 survivors every winding tachyon has fractional electric charge; only 53 and 57 also have
SM-neutral ones (Y = 0, SU(2) and colour singlets, hidden singlets), exactly degenerate with the charged ones when the two
Wilson-line tori carry the same modulus.  In a non-supersymmetric string no D-term fixes the condensation direction.
Findings (explicit enumerator of Pass 11107, cached torus tables):
  A. the degeneracy survives any common B-field: at T_WL1 = T_WL2 = x + iy (x = 0 ... 1/2, three radii each) the lowest
     neutral and charged levels agree to all digits (data/w33_pass11112_bfield_degeneracy.json);
  B. its origin: the neutral tachyons wind ONLY in Wilson-line torus 1 (N = (+-1, 0), (+-2, 0)), the charged ones ONLY in
     torus 2 (N = (0, +-1), (0, +-2)).  Shrinking torus 1 alone (T_WL2 = 3i) produces only neutral tachyons, shrinking
     torus 2 alone only charged ones, with identical levels (-0.0638 / -0.1215 / -0.1792 at Im T = 1.4 / 1.2 / 1.0 in
     57; -0.0959 / -0.1536 / -0.2113 in 53) (data/w33_pass11112_torus_resolved_tachyons.json);
  C. which torus reaches its critical radius first is decided by the potential (equal kinetic metrics 1/(4 y^2)).  Model
     57 (one-loop engine, data/w33_pass11112_pair_integrands_model57.json):
        Lambda_beta(1.6, 1.8) = 445.86 < Lambda_beta(1.8, 1.6) = 447.15,
        Lambda_beta(1.6, 2.0) = 496.90 < Lambda_beta(2.0, 1.6) = 499.17
     -- the asymmetry sits in the exact strip + cap parts of the modular integral (-1.21, -2.10), not in the fitted tail
     (-0.08, -0.17): the potential pulls the NEUTRAL torus down faster.  Model 53 likewise: Lambda_beta(1.8, 2.0) - (2.0, 1.8)
     = -0.94, (1.8, 2.2) - (2.2, 1.8) = -1.67 (strip + cap -0.85, -1.48).
Reading: in BOTH models with neutral tachyons (53, 57) the small-radius instability is reached first along an SM-neutral
direction.  Its condensation
breaks only extra U(1)s -- electromagnetism and the SM gauge group survive.  (Model 2, whose tachyons are all charged,
prefers its other torus, Pass 11110: the asymmetry is model-specific.)  The endpoint of the condensation is not computed.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
TORUS = ROOT / "data" / "w33_pass11112_torus_resolved_tachyons.json"
BFIELD = ROOT / "data" / "w33_pass11112_bfield_degeneracy.json"
PAIRS = {m: ROOT / "data" / f"w33_pass11112_pair_integrands_model{m}.json" for m in ('57', '53')}
OUT = ROOT / "data" / "w33_pass11112_neutral_tachyon_corridor.json"


def pair_lambdas(path):
    import w33_pass11106_one_loop_orbifold_point as P6
    d = json.loads(path.read_text())
    out = {}
    for r in d['runs']:
        s, c, t = P6.fitted(d['pts'], [v[0] for v in r['vals']])
        out[tuple(r['pair'])] = dict(strip_cap=-(s + c) / 2, tail=-t / 2, Lambda=-(s + c + t) / 2)
    return out


def summarize():
    torus = json.loads(TORUS.read_text())
    bf = json.loads(BFIELD.read_text())
    sep = []
    for m, T1, T2, lv in torus:
        small1 = complex(T1).imag < complex(T2).imag
        sep.append((lv['neutral'] is not None and lv['charged'] is None) if small1 else (lv['charged'] is not None and lv['neutral'] is None))
    degenerate = all(r[2]['neutral'] == r[2]['charged'] for r in bf)
    res = dict(pass_id=11112, bfield_degenerate=degenerate, bfield_points=len(bf),
               torus_separation=all(sep), torus_points=len(sep))
    for m, path in PAIRS.items():
        if not path.exists():
            continue
        L = pair_lambdas(path)
        comps = []
        keys = sorted(L)
        for a in keys:
            b = (a[1], a[0])
            if a[0] < a[1] and b in L:
                comps.append(dict(small_torus1=list(a), small_torus2=list(b),
                                  dLambda=L[a]['Lambda'] - L[b]['Lambda'],
                                  dLambda_strip_cap=L[a]['strip_cap'] - L[b]['strip_cap'],
                                  dLambda_tail=L[a]['tail'] - L[b]['tail']))
        res[f'model{m}'] = dict(Lambda={f"{k[0]}:{k[1]}": v['Lambda'] for k, v in L.items()}, comparisons=comps,
                                neutral_torus_shrinks_first=all(c['dLambda'] < 0 and c['dLambda_strip_cap'] < 0 for c in comps))
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    summarize()
