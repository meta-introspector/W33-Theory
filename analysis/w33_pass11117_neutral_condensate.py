#!/usr/bin/env python3
"""Pass 11117: what the Standard-Model-neutral tachyon condensate of models 53 and 57 breaks -- exactly one U(1), nothing
of the Standard Model, and nothing of the family symmetry.

The neutral winding tachyons of Pass 11112 (Wilson-line torus 1 shrunk to 1.4i, torus 2 at 3i, family torus at rho;
explicit enumerator of Pass 11107; charges from the orbifolder U(1) basis of the field dumps, frozen in
data/w33_pass11117_neutral_condensate_charges.json):
  * every neutral tachyon has Y = 0 and is a singlet of SU(3)_c, SU(2)_L and every hidden non-abelian factor;
  * their U(1) charge vectors are +-one vector: rank 1 -- the condensate breaks EXACTLY ONE U(1) combination; it has a
    component along the anomalous U(1) (charge +-4), already massive by the Green-Schwarz mechanism, and along
    non-anomalous U(1)s (model 57: 8 U(1)s; 53: 9);
  * the states wind in the Wilson-line torus with no fixed-point label: the physical state is the Z3-invariant
    combination of winding images, a singlet of the family Delta(54) -- the condensate cannot lift m_c = m_u;
  * self-couplings: the states carry charge +-l_A and winding classes N = +-1, +-2 (N = 1 and -2 share l_A, so do -1
    and 2); a term needs zero total gauge charge, so only |T|^2-type combinations (even in the charge) are allowed --
    no cubic term: the condensation starts as a second-order (Higgs-like) instability.
Scope: the endpoint needs the quartic coefficients (string four-point amplitudes of the winding states) and the
back-reaction on the moduli; neither is computed.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "w33_pass11117_neutral_condensate_charges.json"
OUT = ROOT / "data" / "w33_pass11117_neutral_condensate.json"


def summarize():
    d = json.loads(DATA.read_text())
    res = dict(pass_id=11117)
    for m, v in d.items():
        qs = {tuple(r['q']) for r in v['states']}
        charge_sums = {tuple(round(a + b, 6) for a, b in zip(x, y)) for x in qs for y in qs}
        cubic = any(all(abs(sum(t)) < 1e-9 for t in zip(a, b, c)) for a in qs for b in qs for c in qs)
        res[m] = dict(n_u1=v['n_u1'], broken_u1_rank=v['broken_u1_rank'], hypercharge_zero=v['hypercharge_zero'],
                      nonabelian_singlet=v['nonabelian_singlet'], charged_under_anomalous=v['charged_under_anomalous'],
                      charge_vectors=len(qs), opposite_pair=len(qs) == 2 and any(all(abs(x) < 1e-9 for x in s) for s in charge_sums),
                      cubic_allowed=cubic, levels=sorted({round(r['Delta'], 4) for r in v['states']}))
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    summarize()
