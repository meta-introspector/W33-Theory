#!/usr/bin/env python3
"""Pass 11025: Schur-cover correction plus native E6 cubic tensor-quiver firewall."""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import w33_pass10975_e6_cubic_reynolds_clock_weld as P75
import w33_pass11021_minimal_signed_clock_carrier as P21
import w33_pass11022_binary_octahedral_clock_decomposition as P22
import w33_pass11023_mckay_e7_clock_saturation as P23

OUT = ROOT / "data/w33_pass11025_native_cubic_mckay_firewall.json"
P11022 = ROOT / "data/w33_pass11022_binary_octahedral_clock_decomposition.json"


def payload():
    signed, e6_to_h, parent, affine, p72 = P75.load_inputs()
    e6, dirs, rows, perms, section = P21.split_section()
    p_to_m = {r[3]: tuple(map(int, r[0])) for r in rows}
    group = sorted({p_to_m[p] for p in perms})
    classes = P22.conjugacy_classes(group)
    ci = {m: i for i, C in enumerate(classes) for m in C}

    ident = (1, 0, 0, 1)
    involutions = sum(g != ident and P22.mm(g, g) == ident for g in group)
    assert involutions == 13

    table = dict(P22.exact_character_table())
    order8_values = list(table["spin_2a"][5:7])
    assert all(sp.simplify(x*x + 2) == 0 for x in order8_values)

    A = P23.tensor_adjacency("spin_2a")
    B = P23.tensor_adjacency("spin_2b")
    assert A != A.T and B == A.T
    assert len(P23.directed_edges(A)) == 14
    assert len(P23.underlying_edges(A)) == 11

    p22cert = json.loads(P11022.read_text(encoding="utf-8"))
    m = P23.vec(p22cert["decomposition_24"])
    nonlinear = sp.Matrix([0, 0, 1, 1, 1, 1, 1, 1])
    assert A*m == 3*nonlinear and B*m == 3*nonlinear

    chi24 = [24, 0, 0, 6, 0, 0, 0, 2]
    chip = [12, 3, 12, 3, 0, 0, 0, 2]
    chim = [12, -3, -12, 3, 0, 0, 0, 0]
    def sym2(chi, m):
        m2 = P22.mm(m, m)
        return (chi[ci[m]] ** 2 + chi[ci[m2]]) // 2

    def sym3(chi, m):
        m2 = P22.mm(m, m)
        m3 = P22.mm(m2, m)
        return (
            chi[ci[m]] ** 3
            + 3 * chi[ci[m]] * chi[ci[m2]]
            + 2 * chi[ci[m3]]
        ) // 6

    inv3 = sum(sym3(chi24, m) for m in group) // 48
    inv3p = sum(sym3(chip, m) for m in group) // 48
    inv3m = sum(sym3(chim, m) for m in group) // 48
    cross_pmm = sum(chip[ci[m]] * sym2(chim, m) for m in group) // 48
    cross_mpp = sum(chim[ci[m]] * sym2(chip, m) for m in group) // 48
    assert (inv3, inv3p, inv3m, cross_pmm, cross_mpp) == (71, 23, 0, 48, 0)

    noncentral = set(i for i, h in e6_to_h.items() if h[:2] != (0, 0))
    center = set(range(27)) - noncentral
    triad_center_hist = Counter(sum(i in center for i in tri) for tri in signed)
    assert triad_center_hist == Counter({0: 32, 1: 12, 3: 1})
    z = next(g for g in section if p_to_m[g[0]] == (2, 0, 0, 2))
    p, signs = z
    seen = set()
    pairs = []
    for i in sorted(noncentral):
        if i in seen:
            continue
        j = p[i]
        assert j != i and p[j] == i and signs[i] == signs[j]
        seen |= {i, j}
        pairs.append((i, j, -1 if signs[i] else 1))
    assert len(pairs) == 12

    # Expand the noncentral restriction in a central-eigenbasis.  A common
    # factor 1/(2 sqrt(2)) is stripped; remaining coefficients are integers.
    idx = {}
    for a, (i, j, eps) in enumerate(pairs):
        idx[i] = (a, eps, False)
        idx[j] = (a, eps, True)

    poly = defaultdict(int)
    for tri, coeff in signed.items():
        if not set(tri) <= noncentral:
            continue
        terms = [(1, ())]
        for i in tri:
            a, eps, is_j = idx[i]
            opts = [(1, (0, a)), ((-1 if is_j else 1), (1, a))]
            if is_j:
                opts = [(eps * c, v) for c, v in opts]
            terms = [(c*d, vs+(v,)) for c, vs in terms for d, v in opts]
        for c, vs in terms:
            poly[tuple(sorted(vs))] += coeff * c
    poly = {k: v for k, v in poly.items() if v}
    minus_count = Counter(sum(x[0] for x in k) for k in poly)
    coeff_hist = {
        str(m): {str(k): int(v) for k, v in sorted(Counter(
            c for mon, c in poly.items() if sum(x[0] for x in mon) == m
        ).items())}
        for m in sorted(minus_count)
    }
    assert minus_count == Counter({2: 48, 0: 16})
    assert all(sum(x[0] for x in mon) % 2 == 0 for mon in poly)
    assert set(abs(c) for c in poly.values()) == {2}

    checks = {
        "group_order_48": len(group) == 48,
        "gl23_has_13_nonidentity_involutions": involutions == 13,
        "order8_faithful_character_values_square_to_minus2":
            all(sp.simplify(x*x + 2) == 0 for x in order8_values),
        "faithful_tensor_quiver_is_directed": A != A.T,
        "conjugate_faithful_quiver_is_transpose": B == A.T,
        "directed_edges_14_underlying_edges_11":
            len(P23.directed_edges(A)) == 14 and len(P23.underlying_edges(A)) == 11,
        "clock_saturation_survives_character_correction":
            A*m == 3*nonlinear and B*m == 3*nonlinear,
        "symmetric_cubic_invariant_dimension_71": inv3 == 71,
        "plus_cubic_invariants_23": inv3p == 23,
        "plus_times_sym2minus_invariants_48": cross_pmm == 48,
        "odd_central_parity_invariants_zero": inv3m == 0 and cross_mpp == 0,
        "native_noncentral_triads_32": triad_center_hist[0] == 32,
        "native_eigenbasis_monomials_64": len(poly) == 64,
        "native_only_even_minus_parity": set(minus_count) == {0, 2},
        "native_16_plus3_48_plusminus2": minus_count == Counter({0: 16, 2: 48}),
    }
    assert all(checks.values())
    return {
        "schema": "w33.pass11025.schur-cover-cubic-firewall.v2",
        "status": "PASS",
        "headline": (
            "Route 2 produces a correction and a surviving theorem. GL2(3) is "
            "the plus Schur cover, not binary octahedral; the faithful order-8 "
            "characters are +/-i*sqrt(2), so the corrected tensor quiver is directed "
            "with 14 arrows and an 11-edge underlying graph, not affine E7. The "
            "clock saturation identity survives. Independently, Sym^3(V24)^GL2(3) "
            "has dimension 71, so the native E6 cubic is highly non-unique under "
            "the finite symmetry and requires its extra incidence/sign data."
        ),
        "cover_and_tensor_correction": {
            "group": "GL2(3) = 2^+S4",
            "nonidentity_involutions": involutions,
            "binary_octahedral_identification": False,
            "classical_SL2_McKay_applicable": False,
            "faithful_order8_character_values": ["-i*sqrt(2)", "+i*sqrt(2)"],
            "spin2a_directed_edges": len(P23.directed_edges(A)),
            "underlying_undirected_edges": len(P23.underlying_edges(A)),
            "conjugate_quiver_relation": "B = A^T",
            "clock_saturation_survives": True,
        },
        "invariant_cubic_space": {
            "total_dimension": inv3,
            "Sym3_Vplus": inv3p,
            "Vplus_tensor_Sym2_Vminus": cross_pmm,
            "Sym3_Vminus": inv3m,
            "Vminus_tensor_Sym2_Vplus": cross_mpp,
        },
        "native_E6_restriction": {
            "total_signed_triads": len(signed),
            "triads_by_number_of_central_coordinates": {
                str(k): int(v) for k, v in sorted(triad_center_hist.items())
            },
            "noncentral_triads": triad_center_hist[0],
            "central_eigenbasis_distinct_monomials": len(poly),
            "monomials_by_minus_count": {
                str(k): int(v) for k, v in sorted(minus_count.items())
            },
            "integer_coefficients_after_common_factor_removed": coeff_hist,
            "common_factor_removed": "1/(2*sqrt(2))",
        },
        "tensor_quiver_firewall": {
            "conclusion": (
                "The corrected GL2(3) tensor quiver and exact finite symmetry do "
                "not select a unique cubic interaction: there are 71 independent "
                "invariant symmetric cubics on V24. The native E6 incidence/sign "
                "tensor is additional structure not recoverable from the tensor "
                "quiver alone."
            ),
            "generalized_mckay_reading": (
                "GL2(3) is a finite subgroup of GL2(C), not SL2(C), so the relevant "
                "generalized GL2 McKay/reconstruction-algebra setting is not the "
                "classical ADE correspondence. No unique cubic potential follows "
                "from the finite tensor quiver by itself."
            ),
        },
        "boundary": (
            "This is an invariant-theory, Schur-cover and native-tensor audit. "
            "It corrects the prior affine-E7/binary-octahedral wording while "
            "preserving the exact 24D decomposition and saturation identities. "
            "It does not rule out a richer GL2 reconstruction algebra or a "
            "quiver-with-potential formulation using additional E6 data."
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
        "inv3": p["invariant_cubic_space"]["total_dimension"],
        "native_monomials": p["native_E6_restriction"]["central_eigenbasis_distinct_monomials"],
        "parity": p["native_E6_restriction"]["monomials_by_minus_count"],
    }, sort_keys=True))


if __name__ == "__main__":
    raise SystemExit(main())
