#!/usr/bin/env python3
"""Pass 11111: three generations cannot carry the two-qutrit Pauli group 3^(1+4) (whose commutation geometry is W(3,3))
-- neither as a family group nor glued to colour.  This explains the 104/104 single-qutrit count of Pass 11105.

A. Representation theory (computed).  The extraspecial group 3^(1+2k) (the k-qutrit Pauli group, built here as explicit
   3^k x 3^k clock/shift matrices) has 3^(2k) linear characters and exactly two irreducible representations on which the
   centre acts nontrivially, each of dimension 3^k (the defining representation is irreducible: <chi, chi> = 1; the
   degrees satisfy 3^(2k) + 2 (3^k)^2 = |G|).  In the orbifold the centre is the point-group twist, acting on a theta^l
   field as w^l != 1 on twisted matter.  So twisted generations carrying a faithful 3^(1+2k) come in multiples of 3^k:
   THREE generations force k = 1 -- one qutrit, one Wilson-line-free torus (Pass 11105: 104/104).  With two free tori
   (k = 2) twisted families come in 9s; untwisted families (l = 0) see the Heisenberg group only through its abelian
   quotient, so the symplectic (W(3,3)) structure is invisible on them.
B. The colour loophole (tested).  Quarks carry generation (x) colour: 3 x 3 = 9, the faithful dimension of 3^(1+4), and
   SU(3)_colour contains its own clock/shift H27 with centre Z(SU(3)).  The central product H27_family o H27_colour = 3^(1+4)
   would act on quarks if the point-group phase w^l agreed with a Standard-Model gauge-centre element on every quark:
   w^l = w^(a t) exp(2 pi i b Y) (-1)^(2 s T3)  (t = colour triality).  In all 12 survivors the SM quarks are
   Q (3, 1/6, l = 2), ubar (3bar, -2/3, l = 2), dbar (3bar, 1/3, l = 2): no (a, b, s) exists (doubling the Q condition
   contradicts the dbar condition mod 1).  The orbifold twist is not a gauge-centre element on quarks, so the family H27
   cannot be glued to colour: the 243 group is not realised this way either.
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction as F
from itertools import product
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
OUT = ROOT / "data" / "w33_pass11111_no_243_family_group.json"
W = np.exp(2j * np.pi / 3)
QUARKS = {"Q": (2, 1, F(1, 6), True), "ubar": (2, 2, F(-2, 3), False), "dbar": (2, 2, F(1, 3), False)}   # (l, triality, Y, doublet)


def pauli_group(k):
    X, Z = np.roll(np.eye(3), 1, axis=0), np.diag([1, W, W * W])
    gens = []
    for i in range(k):
        for P in (X, Z):
            M = np.array([[1]], dtype=complex)
            for j in range(k):
                M = np.kron(M, P if j == i else np.eye(3))
            gens.append(M)
    key = lambda M: tuple(np.round(M, 6).flatten().tolist())
    I = np.eye(3 ** k, dtype=complex)
    seen, fr = {key(I): I}, [I]
    while fr:
        nx = []
        for A in fr:
            for g in gens:
                B = g @ A
                if key(B) not in seen:
                    seen[key(B)] = B
                    nx.append(B)
        fr = nx
    return list(seen.values())


def part_a():
    out = {}
    for k in (1, 2):
        G = pauli_group(k)
        chi = np.array([np.trace(g) for g in G])
        centre = [g for g in G if np.allclose(g, g[0, 0] * np.eye(g.shape[0]))]
        out[f"k={k}"] = dict(order=len(G), centre_order=len(centre), defining_dim=3 ** k,
                             defining_irreducible=bool(abs((np.abs(chi) ** 2).sum() / len(G) - 1) < 1e-9),
                             degree_sum_check=3 ** (2 * k) + 2 * 9 ** k == len(G))
    out["three_generations_allow_k"] = [k for k in (1, 2, 3) if 3 % (3 ** k) == 0]
    return out


def gluing_solutions(quarks=QUARKS):
    sols = []
    for a, b, s in product(range(3), range(6), range(2)):
        if all((F(l, 3) - F(a * t, 3) - b * Y - (F(s, 2) if dbl else 0)) % 1 == 0 for l, t, Y, dbl in quarks.values()):
            sols.append((a, b, s))
    return sols


def part_b_spectra(ours, theirs):
    """record the (l, triality, Y) of the SM quarks in the 12 survivors"""
    import w33_pass11097_fractional_fermion_masses as P97
    import w33_pass11102_higher_order_quark_hierarchy as H
    import w33_so16_selection_rules as S
    US, TH = S.parse_ours(ours), S.parse_theirs(theirs)
    rec = {}
    for m in H.BEST:
        us, th = US[m], TH[m]
        P97.prepare(us, th)
        c, oc, ow, qcol = S.sm_data(th, us)
        qb = P97.conj(qcol, 'A2')
        ferm = [f for f in us['fields'] if f['m'] == 2]
        get = lambda col, w, Y: sorted({f['l'] for f in ferm if f['col'] == col and f['w'] == w and f['Y'] == Y})
        rec[str(m)] = dict(Q=get(qcol, 2, F(1, 6)), ubar=get(qb, 1, F(-2, 3)), dbar=get(qb, 1, F(1, 3)))
    return rec


def main(ours=None, theirs=None):
    res = dict(pass_id=11111, A=part_a(), B_gluing_solutions_for_l2_quarks=gluing_solutions(),
               B_control_trivial_twist=gluing_solutions({k: (0,) + v[1:] for k, v in QUARKS.items()}))
    if ours:
        res['B_quark_sectors'] = part_b_spectra(ours, theirs)
    elif OUT.exists() and 'B_quark_sectors' in json.loads(OUT.read_text()):
        res['B_quark_sectors'] = json.loads(OUT.read_text())['B_quark_sectors']
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    main(*sys.argv[1:3]) if len(sys.argv) > 2 else main()
