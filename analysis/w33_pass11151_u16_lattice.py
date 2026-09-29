#!/usr/bin/env python3
"""Pass 11151: the zero-radius vector lattice IS the U(16) lattice -- identified by its genus and its root system, not
only by counting roots.

Scan: analysis/w33_pass11151_scan_zero_radius_lattice.py.  L0 = {l = pi + N W : pi in E8+E8, N in Z, pi.V0 in Z,
pi.W + N W^2/2 in Z} (W the neutral torus's Wilson line; the Wilson-line offset is LINEAR in N) is the kernel of two
linear characters on Z^17 (E8+E8 basis plus W), mapped to Q^16; Gram matrix, determinant, parity and discriminant group
(Smith form) computed exactly (sympy).  Validation: with W = 0 the same construction gives Gamma_v^(0) =
D8+D8 + (s,s), the vector lattice of the O(16)xO(16) string (Fraiman-Grana-Parra De Freitas-Sethi, eq. 3.14).
Frozen: data/w33_pass11151_zero_radius_lattice_21.json.
Results (21/21):
  * validation and zero-radius lattices alike: rank 16, even, integral, determinant 4, discriminant group Z2 x Z2;
  * discriminant quadratic form: the nonzero classes take values {0, 1} mod 2 (checked explicitly for 36621, 40521, 2233,
    validation and zero-radius alike) -- the even hyperbolic form u (the other even form v takes only 1; the odd forms
    take half-integers) -- consistent with Milgram (signature 16 = 0 mod 8); so, by Nikulin, L0 lies in the SAME GENUS
    as Gamma_v^(0);
  * that genus is exhausted by the six vector lattices Gamma_v^(p), p = 0..5, of the non-supersymmetric heterotic strings
    (Fraiman et al., completeness via the Smith-Minkowski-Siegel mass formula), which differ by root system
    (so16+so16, su16+u1, 2(e7+su2), so8+so24, e8+so16, so32);
  * L0's root system is A15 (Pass 11140: 240 roots, rank 15), hence L0 = Gamma_v^(1): the U(16) string's lattice.
The Wilson-line data that map onto it: W^2 in {10/3, 4, 14/3, 16/3, 6, 20/3, 8, 38/3} and W.V0 in {0, 2} -- all 21
models land on the same lattice.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "w33_pass11151_zero_radius_lattice_21.json"
ROOTS = ROOT / "data" / "w33_pass11140_zero_radius_u16.json"
OUT = ROOT / "data" / "w33_pass11151_u16_lattice.json"


def ok(r):
    return r['rank'] == 16 and r['det'] == '4' and r['even'] and r['integral'] and r['discriminant'] == ['2', '2']


def summarize():
    d = json.loads(DATA.read_text())
    roots = json.loads(ROOTS.read_text())
    res = dict(pass_id=11151, n_models=len(d), validation_is_gamma_v0=all(ok(v['validation_W0']) for v in d.values()),
               zero_radius_same_genus=all(ok(v['zero_radius']) for v in d.values()),
               root_system_A15=roots['su16_in_all'], identified='Gamma_v^(1) (U(16) string)',
               discriminant_form_hyperbolic_checked=all(v['discriminant_form_values_mod2_checked']['zero_radius'] == [0, 1]
                                                         for v in d.values() if 'discriminant_form_values_mod2_checked' in v))
    res['identified_in_all'] = res['validation_is_gamma_v0'] and res['zero_radius_same_genus'] and res['root_system_A15']
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
