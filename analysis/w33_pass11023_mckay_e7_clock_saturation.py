#!/usr/bin/env python3
"""Pass 11023: exact McKay-E7 placement of the signed clock carrier.

Using the verified GL2(3) character table from Pass 11022, reconstruct the
McKay graph by tensoring irreducibles with either faithful 2D spinor.
The graph is affine E7.  The actual 24D clock multiplicity vector has the
strong saturation law S tensor V24 = 3 times the sum of all nonlinear irreps.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "analysis") not in sys.path:
    sys.path.insert(0, str(ROOT / "analysis"))

import w33_pass11022_binary_octahedral_clock_decomposition as P80

OUT = ROOT / "data/w33_pass11023_mckay_e7_clock_saturation.json"
P11022 = ROOT / "data/w33_pass11022_binary_octahedral_clock_decomposition.json"
NAMES = (
    "trivial", "det", "quotient_2d", "spin_2a", "spin_2b",
    "standard3_twist", "standard3", "spin_4",
)


def inner(a, b):
    return sp.simplify(sum(
        n*x*sp.conjugate(y)
        for n, x, y in zip(P80.CLASS_SIZES, a, b)
    ) / 48)


def tensor_adjacency(spin_name):
    rows = dict(P80.exact_character_table())
    spin = rows[spin_name]
    matrix = []
    for name in NAMES:
        product = tuple(sp.expand(x*y) for x, y in zip(spin, rows[name]))
        mult = [int(inner(product, rows[target])) for target in NAMES]
        assert all(x >= 0 for x in mult)
        matrix.append(mult)
    A = sp.Matrix(matrix)
    assert A == A.T
    return A


def edges(A):
    return tuple(
        (NAMES[i], NAMES[j], int(A[i, j]))
        for i in range(8) for j in range(i + 1, 8)
        if A[i, j]
    )
def graph_checks(A):
    E = edges(A)
    assert len(E) == 7
    assert all(m == 1 for _, _, m in E)
    deg = [sum(int(A[i, j]) for j in range(8)) for i in range(8)]
    assert sorted(deg) == [1, 1, 1, 2, 2, 2, 2, 3]

    # Explicit connectedness.
    seen = {0}
    front = [0]
    while front:
        i = front.pop()
        for j in range(8):
            if A[i, j] and j not in seen:
                seen.add(j)
                front.append(j)
    assert len(seen) == 8

    # Remove the affine/trivial leaf: finite E7 is a six-chain with one branch.
    finite = list(range(1, 8))
    finite_deg = {
        i: sum(int(A[i, j]) for j in finite) for i in finite
    }
    assert sorted(finite_deg.values()) == [1, 1, 1, 2, 2, 2, 3]
    assert finite_deg[NAMES.index("spin_4")] == 3
    return deg


def vec(mapping):
    return sp.Matrix([int(mapping.get(n, 0)) for n in NAMES])
def payload():
    parent = json.loads(P11022.read_text(encoding="utf-8"))
    rows = dict(P80.exact_character_table())

    A = tensor_adjacency("spin_2a")
    B = tensor_adjacency("spin_2b")
    degA = graph_checks(A)
    degB = graph_checks(B)

    dims = sp.Matrix([int(rows[n][0]) for n in NAMES])
    assert list(dims) == [1, 1, 2, 2, 2, 3, 3, 4]
    assert A*dims == 2*dims
    assert B*dims == 2*dims

    m = vec(parent["decomposition_24"])
    mplus = vec(parent["central_minus_I"]["plus_decomposition"])
    mminus = vec(parent["central_minus_I"]["minus_decomposition"])
    assert m == mplus + mminus
    assert list(m) == [2, 1, 0, 0, 0, 1, 2, 3]

    nonlinear = sp.Matrix([0, 0, 1, 1, 1, 1, 1, 1])
    odd = sp.Matrix([0, 0, 0, 1, 1, 0, 0, 1])
    even_nontrivial = sp.Matrix([0, 0, 1, 0, 0, 1, 1, 0])
    for X in (A, B):
        assert X*m == 3*nonlinear
        assert X*mplus == 3*odd
        assert X*mminus == 3*even_nontrivial
        assert X*nonlinear == dims - nonlinear
        assert X*(X + sp.eye(8))*m == 3*dims

    checks = {
        "spin2a_graph_tree_8_nodes_7_edges": len(edges(A)) == 7,
        "spin2b_graph_tree_8_nodes_7_edges": len(edges(B)) == 7,
        "finite_E7_branch_at_spin4":
            degA[NAMES.index("spin_4")] == 3,
        "dimension_vector_is_affine_null_mark":
            A*dims == 2*dims and B*dims == 2*dims,
        "clock_tensor_spin2a_saturates_nonlinear":
            A*m == 3*nonlinear,
        "clock_tensor_spin2b_saturates_nonlinear":
            B*m == 3*nonlinear,
        "plus_tensors_to_all_odd_nodes":
            A*mplus == 3*odd and B*mplus == 3*odd,
        "minus_tensors_to_all_even_nonlinear_nodes":
            A*mminus == 3*even_nontrivial and B*mminus == 3*even_nontrivial,
        "green_identity":
            A*(A + sp.eye(8))*m == 3*dims,
    }
    assert all(checks.values())
    return {
        "schema": "w33.pass11023.mckay-e7-clock-saturation.v1",
        "status": "PASS",
        "headline": (
            "The binary-octahedral McKay graph reconstructed from the exact "
            "GL2(3) character table is affine E7. For either faithful 2D "
            "spinor S, the actual 24D signed clock carrier obeys "
            "S tensor V24 = 3 times the direct sum of every nonlinear irrep."
        ),
        "nodes": [
            {
                "name": n,
                "dimension": int(dims[i]),
                "clock_multiplicity": int(m[i]),
                "plus_multiplicity": int(mplus[i]),
                "minus_multiplicity": int(mminus[i]),
            }
            for i, n in enumerate(NAMES)
        ],
        "spin2a_mckay_edges": [list(x) for x in edges(A)],
        "spin2b_mckay_edges": [list(x) for x in edges(B)],
        "affine_dimension_vector": list(map(int, dims)),
        "clock_multiplicity_vector": list(map(int, m)),
        "exact_identities": {
            "full":
                "S tensor V24 = 3*(quotient_2d + spin_2a + spin_2b + "
                "standard3_twist + standard3 + spin_4)",
            "plus":
                "S tensor Vplus = 3*(spin_2a + spin_2b + spin_4)",
            "minus":
                "S tensor Vminus = 3*(quotient_2d + standard3_twist + standard3)",
            "adjacency_vector":
                "A*m = 3*n_nonlin",
            "affine_mark":
                "A*d = 2*d",
            "second_order":
                "A*(A+I)*m = 3*d",
        },
        "interpretation": (
            "Tensoring by a defining binary-octahedral spinor flips central "
            "parity. The clock carrier is balanced so that each 12D parity half "
            "maps to three uniform copies of all three nontrivial nodes on the "
            "opposite McKay bipartition."
        ),
        "boundary": (
            "This is an identity in the complex representation ring of the "
            "finite group GL2(3). Affine E7 here is the classical McKay graph, "
            "not an E7 Lie-algebra field content, gauge symmetry, interaction "
            "Lagrangian, or continuum unification result."
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
        "m": p["clock_multiplicity_vector"],
        "edges": len(p["spin2a_mckay_edges"]),
        "full_identity": p["exact_identities"]["full"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
