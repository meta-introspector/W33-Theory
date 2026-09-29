#!/usr/bin/env python3
"""Pass 11140: at zero radius of the neutral torus the massless gauge algebra is su(16) + u(1) -- the T-dual end of the
neutral exit is the non-supersymmetric U(16) heterotic string (2 tachyons, matching Pass 11134).

Scan: analysis/w33_pass11140_scan_zero_radius_roots.py.  As the neutral Wilson-line torus shrinks (Pass 11128) a state
stays massless iff its momentum on that torus can be cancelled exactly: l = pi + N W (N mod 3, W the torus's Wilson
line) with l^2 = 2, Wilson-line offset pi.W + N W^2/2 integral, and the Witten projection of the NS sector
(k = 0: phase +1; k = 1: phase -1, as for the tachyon).  Validation: N = 0 gives the 224 roots of so(16)+so(16), and the
Witten-twisted sector gives no vectors at N = 0 (its 256 norm-2 states are the (128,1)+(1,128) fermions).  Three
candidate phase rules for winding states agree.  Without the offset condition the vector set is NOT a root system (not
reflection-closed) -- the condition is required, not optional.  Frozen: data/w33_pass11140_zero_radius_roots_21.json.
Results (21/21 neutral-exit survivors, other Wilson-line torus decompactified):
  * 240 roots, rank 15, one irreducible component: A15 = su(16) (plus u(1) to rank 16), reflection-closed;
  * this is the gauge algebra of the charge lattice Upsilon^(1)_16 ~ su16 + u1 of Fraiman-Grana-Parra De Freitas-Sethi
    (arXiv:2307.13745, eq. 3.15), the U(16) non-supersymmetric heterotic string, which has 2^1 = 2 tachyons -- exactly
    the TWO norm-1 states (one complex species) found at zero radius in Pass 11134;
  * with the other Wilson-line torus kept finite, the algebra breaks to products of su(n): e.g. su(8)+su(4)+su(4),
    su(7)+su(5)+su(4), su(10)+su(4)+su(2) (model-dependent; the finite_other_torus_algebras field lists them).
CORRECTS Pass 11128's remark (and resolves Pass 11134's open identification): the T-dual string is not the SO(32) or
SO(16)xE8 tachyonic string but the U(16) one.
Scope: the ten-dimensional parent (the Z3 orbifold projection of the family directions is not imposed).
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "w33_pass11140_zero_radius_roots_21.json"
OUT = ROOT / "data" / "w33_pass11140_zero_radius_u16.json"


def summarize():
    d = json.loads(DATA.read_text())
    key = 'a: pi.V0 in Z|offset_int=False'
    su16 = all(v[key]['n'] == 240 and v[key]['not_closed'] == 0 and v[key]['comps'] == [[15, 240]] for v in d.values())
    rules_agree = all(len({(v[k]['n'], str(v[k]['comps'])) for k in v if k.endswith('offset_int=False')}) == 1 for v in d.values())
    finite = sorted({str(v['a: pi.V0 in Z|offset_int=True']['comps']) for v in d.values()})
    res = dict(pass_id=11140, n_models=len(d), su16_in_all=su16, phase_rules_agree=rules_agree,
               finite_other_torus_algebras=finite, identified_string='U(16) non-supersymmetric heterotic (Upsilon^(1)_16)',
               tachyon_count_match=2)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
