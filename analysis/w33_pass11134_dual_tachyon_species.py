#!/usr/bin/env python3
"""Pass 11134: the ground-state tachyon at zero radius is ONE complex species, intrinsically a winding (T-dual
momentum) state -- not a vector-class multiplet.

Scan: analysis/w33_pass11134_scan_dual_tachyon.py.  As the neutral torus shrinks (Pass 11128) every winding class
becomes light; the Delta -> -1/2 ground states are the vectors l = pi + V0 + N W (W the neutral torus's Wilson line,
N mod 3 since 3W lies in the lattice) with l^2 = 1 and the beta phase -1 of Pass 11107.
Frozen: data/w33_pass11134_dual_tachyon_multiplet_21.json.

Results (21/21 neutral-exit survivors):
  * exactly TWO such vectors: one with N = 1, its conjugate with N = 2; none with N = 0 -- the tachyon exists only
    through the Wilson-line winding, i.e. it is a momentum state of the T-dual circle, one complex species (its KK tower
    is the 4 -> 6 -> 8 growth of Pass 11128);
  * l^2 = 1 is split across the two E8s: l1^2 = 2/9 (13 models), 5/9 (5), 4/9 (2), 7/9 (1) -- ninths, the Z3 Wilson line;
  * so it is NOT the 32 of the SO(32) or the 16 of the SO(16)xE8 tachyonic string: whatever ten-dimensional multiplet it
    descends from, the Wilson lines keep one component.
CORRECTION to the Pass 11128 remark "l^2 = 1, the vector class of the SO(32) / SO(16)xE8 tachyonic strings": l^2 = 1
holds, but the species content does not identify a ten-dimensional string; the identification stays open.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "w33_pass11134_dual_tachyon_multiplet_21.json"
OUT = ROOT / "data" / "w33_pass11134_dual_tachyon_species.json"


def summarize():
    d = json.loads(DATA.read_text())
    split = Counter(next(iter(v['split'])).split('=')[1] for v in d.values())
    res = dict(pass_id=11134, n_models=len(d),
               one_complex_species=all(v['total'] == 2 and v['by_N'] == {'0': 0, '1': 1, '2': 1} for v in d.values()),
               three_W_in_lattice=all(v['three_W_in_lattice'] for v in d.values()),
               l1sq_distribution=dict(split))
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
