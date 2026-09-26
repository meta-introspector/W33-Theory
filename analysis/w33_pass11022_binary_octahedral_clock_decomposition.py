#!/usr/bin/env python3
"""Pass 11022: binary-octahedral decomposition of the exact signed clock carrier.

Pass 11021 gives the minimal 24D noncentral H27 carrier under the split signed
GL2(3) action.  This packet decomposes that exact representation over C.

The central element -I splits the carrier as 12+12.  The + sector descends to
PGL2(3)=S4.  The - sector is exactly three copies of the faithful 4D spinorial
irrep of GL2(3), the binary-octahedral double cover of S4.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "analysis") not in sys.path:
    sys.path.insert(0, str(ROOT / "analysis"))

import w33_pass10975_e6_cubic_reynolds_clock_weld as P75
import w33_pass11020_signed_clock_extension_split as P78
import w33_pass11021_minimal_signed_clock_carrier as P79
OUT = ROOT / "data/w33_pass11022_binary_octahedral_clock_decomposition.json"
P11021 = ROOT / "data/w33_pass11021_minimal_signed_clock_carrier.json"

CLASS_REPS = (
    (1, 0, 0, 1),
    (0, 2, 1, 1),
    (2, 0, 0, 2),
    (0, 2, 1, 2),
    (0, 2, 1, 0),
    (0, 1, 1, 2),
    (0, 1, 1, 1),
    (2, 0, 0, 1),
)
CLASS_ORDERS = (1, 6, 2, 3, 4, 8, 8, 2)
CLASS_SIZES = (1, 8, 1, 8, 6, 6, 6, 12)


def mm(a, b):
    x, y, z, w = a
    p, q, r, s = b
    return ((x*p+y*r) % 3, (x*q+y*s) % 3,
            (z*p+w*r) % 3, (z*q+w*s) % 3)
def minv(a):
    x, y, z, w = a
    det = (x*w-y*z) % 3
    t = pow(det, -1, 3)
    return ((t*w) % 3, (-t*y) % 3, (-t*z) % 3, (t*x) % 3)


def conjugacy_classes(group):
    out = []
    for x in CLASS_REPS:
        out.append(frozenset(mm(mm(g, x), minv(g)) for g in group))
    assert [len(c) for c in out] == list(CLASS_SIZES)
    assert len(frozenset().union(*out)) == 48
    return tuple(out)


def mon_trace(element, subset):
    p, signs = element
    return sum(
        -1 if signs[i] else 1
        for i in subset if p[i] == i
    )


def projective_perm(m):
    pts = ((0, 1), (1, 0), (1, 1), (1, 2))
    idx = {p: i for i, p in enumerate(pts)}
    out = []
    for a, b in pts:
        x = (m[0]*a + m[1]*b) % 3
        y = (m[2]*a + m[3]*b) % 3
        if x:
            t = pow(x, -1, 3)
            p = (1, (t*y) % 3)
        else:
            p = (0, 1)
        out.append(idx[p])
    return tuple(out)


def perm_trace(p):
    return sum(i == j for i, j in enumerate(p))


def exact_character_table():
    rt2 = sp.sqrt(2)
    rows = (
        ("trivial", (1, 1, 1, 1, 1, 1, 1, 1)),
        ("det", (1, 1, 1, 1, 1, -1, -1, -1)),
        ("quotient_2d", (2, -1, 2, -1, 2, 0, 0, 0)),
        ("spin_2a", (2, 1, -2, -1, 0, -rt2, rt2, 0)),
        ("spin_2b", (2, 1, -2, -1, 0, rt2, -rt2, 0)),
        ("standard3_twist", (3, 0, 3, 0, -1, 1, 1, -1)),
        ("standard3", (3, 0, 3, 0, -1, -1, -1, 1)),
        ("spin_4", (4, -1, -4, 1, 0, 0, 0, 0)),
    )
    # Verify this is an orthonormal complete character table.
    inner = []
    for _, a in rows:
        line = []
        for _, b in rows:
            ip = sp.simplify(sum(
                n*x*sp.conjugate(y)
                for n, x, y in zip(CLASS_SIZES, a, b)
            ) / 48)
            line.append(ip)
        inner.append(line)
    assert sp.Matrix(inner) == sp.eye(8)
    assert sum(int(r[1][0])**2 for r in rows) == 48
    return rows


def multiplicities(character, rows):
    out = []
    for name, row in rows:
        m = sp.simplify(sum(
            n*x*sp.conjugate(y)
            for n, x, y in zip(CLASS_SIZES, character, row)
        ) / 48)
        assert m.is_Integer and m >= 0
        out.append((name, int(m)))
    return tuple(out)


def rank(vectors):
    return sp.Matrix.hstack(*[sp.Matrix(v) for v in vectors]).rank()


def act_vec(element, v):
    return P79.act_vec(element, v)
def payload():
    parent = json.loads(P11021.read_text(encoding="utf-8"))
    assert parent["minimal_closure"]["augmentation_span_rank"] == 24

    e6_to_h, dirs, rows75, perms, section = P79.split_section()
    p_to_m = {r[3]: tuple(map(int, r[0])) for r in rows75}
    group = frozenset(p_to_m[p] for p in perms)
    assert len(group) == 48
    classes = conjugacy_classes(group)
    section_by_m = {p_to_m[p]: (p, s) for p, s in section}

    noncentral = tuple(sorted(
        i for i, h in e6_to_h.items() if h[:2] != (0, 0)
    ))
    center = tuple(sorted(
        i for i, h in e6_to_h.items() if h[:2] == (0, 0)
    ))

    chi24 = []
    chi3 = []
    for c in classes:
        vals24 = {mon_trace(section_by_m[m], noncentral) for m in c}
        vals3 = {mon_trace(section_by_m[m], center) for m in c}
        assert len(vals24) == len(vals3) == 1
        chi24.append(next(iter(vals24)))
        chi3.append(next(iter(vals3)))
    assert chi24 == [24, 0, 0, 6, 0, 0, 0, 2]
    assert chi3 == [3, 3, 3, 3, 3, 1, 1, 1]

    table = exact_character_table()
    mult24 = multiplicities(tuple(chi24), table)
    mult3 = multiplicities(tuple(chi3), table)
    assert mult24 == (
        ("trivial", 2),
        ("det", 1),
        ("quotient_2d", 0),
        ("spin_2a", 0),
        ("spin_2b", 0),
        ("standard3_twist", 1),
        ("standard3", 2),
        ("spin_4", 3),
    )
    assert mult3 == (
        ("trivial", 2), ("det", 1), ("quotient_2d", 0),
        ("spin_2a", 0), ("spin_2b", 0),
        ("standard3_twist", 0), ("standard3", 0), ("spin_4", 0),
    )

    # Identify the pulled-back S4 standard 3D character from the P1(F3) action.
    chi_p1 = [perm_trace(projective_perm(m)) for m in CLASS_REPS]
    chi_std = [x - 1 for x in chi_p1]
    assert chi_std == list(dict(table)["standard3"])
    det_char = [
        1 if (m[0]*m[3]-m[1]*m[2]) % 3 == 1 else -1
        for m in CLASS_REPS
    ]
    assert det_char == list(dict(table)["det"])

    zmat = (2, 0, 0, 2)
    z = section_by_m[zmat]
    assert P78.mon_order(z) == 2
    zp, zs = z

    def projector(v, eps):
        zv = act_vec(z, v)
        return tuple(sp.Rational(a + eps*b, 2) for a, b in zip(v, zv))

    fibres = []
    indicators = []
    for d in dirs:
        fibre = tuple(sorted(
            i for i, h in e6_to_h.items()
            if h[:2] != (0, 0) and P75.norm_dir(h[:2]) == d
        ))
        assert len(fibre) == 6
        fibres.append(fibre)
        v = [0] * 27
        for i in fibre:
            v[i] = 1
        indicators.append(tuple(v))

    fibre_splits = []
    for fibre in fibres:
        tr = mon_trace(z, fibre)
        assert tr == 0
        fibre_splits.append(((len(fibre)+tr)//2, (len(fibre)-tr)//2))
    assert fibre_splits == [(3, 3)] * 4

    augmentation = [
        tuple(a-b for a, b in zip(indicators[i], indicators[3]))
        for i in range(3)
    ]

    def projected_orbit_rank(seeds, eps):
        vecs = [
            projector(act_vec(g, v), eps)
            for g in section for v in seeds
        ]
        return rank(vecs)

    ranks = {
        "indicator_plus": projected_orbit_rank(indicators, +1),
        "indicator_minus": projected_orbit_rank(indicators, -1),
        "augmentation_plus": projected_orbit_rank(augmentation, +1),
        "augmentation_minus": projected_orbit_rank(augmentation, -1),
    }
    assert set(ranks.values()) == {12}

    # Central-character projectors split chi24 into its +/- pieces.
    class_index = {}
    for i, c in enumerate(classes):
        for m in c:
            class_index[m] = i
    chi_plus = []
    chi_minus = []
    for i, rep in enumerate(CLASS_REPS):
        j = class_index[mm(zmat, rep)]
        chi_plus.append((chi24[i] + chi24[j]) // 2)
        chi_minus.append((chi24[i] - chi24[j]) // 2)
    assert chi_plus == [12, 3, 12, 3, 0, 0, 0, 2]
    assert chi_minus == [12, -3, -12, 3, 0, 0, 0, 0]

    mult_plus = multiplicities(tuple(chi_plus), table)
    mult_minus = multiplicities(tuple(chi_minus), table)
    assert mult_minus == (
        ("trivial", 0), ("det", 0), ("quotient_2d", 0),
        ("spin_2a", 0), ("spin_2b", 0),
        ("standard3_twist", 0), ("standard3", 0), ("spin_4", 3),
    )
    assert mult_plus == (
        ("trivial", 2), ("det", 1), ("quotient_2d", 0),
        ("spin_2a", 0), ("spin_2b", 0),
        ("standard3_twist", 1), ("standard3", 2), ("spin_4", 0),
    )

    checks = {
        "GL23_order48": len(group) == 48,
        "eight_conjugacy_classes": len(classes) == 8,
        "character_table_orthonormal": True,
        "chi24_exact": chi24 == [24, 0, 0, 6, 0, 0, 0, 2],
        "central_split_12_12": chi_plus[0] == chi_minus[0] == 12,
        "four_fibres_split_3_3": fibre_splits == [(3, 3)] * 4,
        "clock_indicators_generate_both_12s":
            ranks["indicator_plus"] == ranks["indicator_minus"] == 12,
        "augmentation_generates_both_12s":
            ranks["augmentation_plus"] == ranks["augmentation_minus"] == 12,
        "minus_is_three_spin4": dict(mult_minus)["spin_4"] == 3,
        "no_spin2_in_24":
            dict(mult24)["spin_2a"] == dict(mult24)["spin_2b"] == 0,
        "center3_is_two_trivial_plus_det":
            dict(mult3)["trivial"] == 2 and dict(mult3)["det"] == 1,
        "P1_standard_character_identified": chi_std == list(dict(table)["standard3"]),
    }
    assert all(checks.values())

    return {
        "schema": "w33.pass11022.binary-octahedral-clock-decomposition.v1",
        "status": "PASS",
        "headline": (
            "The minimal 24D signed clock carrier decomposes under GL2(3) as "
            "2*1 + det + (3std x det) + 2*3std + 3*4spin. The central -I "
            "splits it canonically into a 12D S4-descending sector and a 12D "
            "sector consisting exactly of three copies of the faithful 4D "
            "spinorial irrep."
        ),
        "group": {
            "name": "GL2(3)",
            "order": 48,
            "projective_quotient": "PGL2(3) ~= S4",
            "class_orders": list(CLASS_ORDERS),
            "class_sizes": list(CLASS_SIZES),
            "character_table_source_check":
                "independently matched against GAP CharacterTable(GL(2,3))",
            "external_identification":
                "classically isomorphic to the binary octahedral double cover 2.S4",
        },
        "character_24": chi24,
        "decomposition_24": {k: v for k, v in mult24},
        "central_minus_I": {
            "plus_dimension": chi_plus[0],
            "minus_dimension": chi_minus[0],
            "plus_character": chi_plus,
            "minus_character": chi_minus,
            "plus_decomposition": {k: v for k, v in mult_plus},
            "minus_decomposition": {k: v for k, v in mult_minus},
            "per_clock_fibre_plus_minus": [list(x) for x in fibre_splits],
            "projected_orbit_ranks": ranks,
        },
        "center3": {
            "character": chi3,
            "decomposition": {k: v for k, v in mult3},
            "central_minus_I_eigenvalue": "+1 on all three dimensions",
        },
        "representation_reading": {
            "S4_descending_12":
                "2*trivial + det + standard3_twist + 2*standard3",
            "spinorial_12": "3*spin_4",
            "clock_seed":
                "the three coarse augmentation generators have nonzero projections "
                "whose exact group orbits span both 12D sectors",
            "fibre_law":
                "each six-coordinate projective clock direction resolves as 3+3 "
                "under the central double-cover involution",
        },
        "boundary": (
            "The words 'spinorial' and 'binary octahedral' are representation-"
            "theoretic: central -I acts as -1 on the 4D irrep and the group is "
            "the double cover of S4. This does not identify these 12 modes with "
            "physical fermions, spacetime spinors, or measured particle species. "
            "The ADE/McKay E7 association is classical external context, not a "
            "derived E7 gauge theory."
        ),
        "checks": checks,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--output", type=Path, default=OUT)
    a = ap.parse_args()
    p = payload()
    text = json.dumps(p, indent=2, sort_keys=True) + "\n"
    if a.check:
        if not a.output.exists() or a.output.read_text(encoding="utf-8") != text:
            raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True, exist_ok=True)
        a.output.write_text(text, encoding="utf-8")
    print(json.dumps({
        "status": p["status"],
        "plus": p["central_minus_I"]["plus_dimension"],
        "minus": p["central_minus_I"]["minus_dimension"],
        "minus_decomposition": p["central_minus_I"]["minus_decomposition"],
        "fibre_splits": p["central_minus_I"]["per_clock_fibre_plus_minus"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
