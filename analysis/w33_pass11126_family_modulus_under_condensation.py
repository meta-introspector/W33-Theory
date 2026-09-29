#!/usr/bin/env python3
"""Pass 11126: the neutral condensate cannot move the family modulus off rho -- so it cannot touch the family hierarchy.

Every tachyon of the six neutral-exit survivors (just below the diagonal onset, at Im T_WL = 1, and on each Wilson-line
torus alone; explicit enumerator of Pass 11107, frozen in data/w33_pass11126_family_torus_quantum_numbers.json) has ZERO
momentum and winding on the family (Wilson-line-free) torus; each model has exactly one complex pair (4 real states),
living on one Wilson-line torus.  Consequences:
  * the condensate does not couple to the family modulus T* at tree level (its vertex operators contain no family-torus
    lattice momentum), so T* enters the potential only through the Narain factor Gamma_{2,2}(T*, rho) of the Witten sectors
    (Pass 11106): the one-loop potential stays exactly SL(2,Z)-invariant in T*, with critical points at i and rho;
  * the condensate is a Delta(54) singlet (Pass 11117), so the family degeneracy m_c = m_u and the tree-level ratio
    m_{c,u}/m_t = |Y_dist/Y_same|(T*) (Pass 11109) are untouched; even a jump of T* between the two fixed points would give
    1/2 (rho) or 0.366 (i), never 0.0036 (needs Im T* ~ 3.2).
The condensate is a Wilson-line-sector phenomenon; the family problem is decided elsewhere.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "w33_pass11126_family_torus_quantum_numbers.json"
OUT = ROOT / "data" / "w33_pass11126_family_modulus_under_condensation.json"


def summarize():
    d = json.loads(DATA.read_text())
    zero = all(v['family_torus_states'] in ([], [[0, 0, 0, 0]]) for m in d.values() for v in m.values())
    counts = {m: {k: v['n'] for k, v in r.items()} for m, r in d.items()}
    single_torus = {m: [k for k in ('torus1_only', 'torus2_only') if r[k]['n']] for m, r in d.items()}
    res = dict(pass_id=11126, family_torus_quantum_numbers_zero=zero, tachyon_state_counts=counts,
               tachyon_torus=single_torus, one_complex_pair_each=all(v['below_crit'] == 4 for v in counts.values()),
               ratio_at_fixed_points={'rho': 0.5, 'i': (3 ** 0.5 - 1) / 2}, needed=0.0036)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    summarize()
