#!/usr/bin/env python3
"""Pass 11023: exact E7-shaped tensor graph for the signed clock carrier.

Using the verified GL2(3) character table from Pass 11022, reconstruct the
representation tensor graph for either faithful 2D irrep. The unlabeled graph
has affine-E7 shape, and the actual 24D clock multiplicity vector obeys the
strong saturation law S tensor V24 = 3 times the sum of all nonlinear irreps.

Pass 11025 corrects the terminology: GL2(3)=2^+S4 is not the binary-octahedral
2^-S4 subgroup of SL2(C), so this is not the classical ADE McKay graph.
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
    return sp.Matrix(matrix)


def directed_edges(A):
    return tuple(
        (NAMES[i], NAMES[j], int(A[i, j]))
        for i in range(8) for j in range(8)
        if A[i, j]
    )


def underlying_edges(A):
    return tuple(sorted({
        tuple(sorted((NAMES[i], NAMES[j])))
        for i in range(8) for j in range(8) if A[i, j]
    }))


def quiver_checks(A):
    E = directed_edges(A)
    assert len(E) == 14
    assert all(m == 1 for _, _, m in E)
    outdeg = [sum(int(A[i, j]) for j in range(8)) for i in range(8)]
    indeg = [sum(int(A[i, j]) for i in range(8)) for j in range(8)]
    assert sorted(outdeg) == [1, 1, 1, 2, 2, 2, 2, 3]
    assert sorted(indeg) == [1, 1, 1, 2, 2, 2, 2, 3]
    assert len(underlying_edges(A)) == 11

    # Explicit strong connectedness.
    for start in range(8):
        seen = {start}
        front = [start]
        while front:
            i = front.pop()
            for j in range(8):
                if A[i, j] and j not in seen:
                    seen.add(j)
                    front.append(j)
        assert len(seen) == 8
    return outdeg, indeg


def vec(mapping):
    return sp.Matrix([int(mapping.get(n, 0)) for n in NAMES])
def payload():
    parent = json.loads(P11022.read_text(encoding="utf-8"))
    rows = dict(P80.exact_character_table())

    A = tensor_adjacency("spin_2a")
    B = tensor_adjacency("spin_2b")
    outA, inA = quiver_checks(A)
    outB, inB = quiver_checks(B)
    assert A != A.T
    assert B == A.T

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
        "spin2a_directed_quiver_14_edges": len(directed_edges(A)) == 14,
        "spin2b_directed_quiver_14_edges": len(directed_edges(B)) == 14,
        "two_faithful_quivers_are_transposes": B == A.T,
        "faithful_tensor_matrix_is_not_symmetric": A != A.T,
        "underlying_undirected_graph_has_11_edges": len(underlying_edges(A)) == 11,
        "dimension_vector_tensor_eigenvalue_2":
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
        "schema": "w33.pass11023.gl23-directed-tensor-clock-saturation.v3",
        "status": "PASS",
        "headline": (
            "For either faithful 2D irrep of GL2(3), the corrected cyclotomic "
            "character table gives a directed non-symmetric tensor quiver. The "
            "two faithful quivers are transposes, each has 14 directed edges and "
            "an 11-edge underlying undirected graph. The clock saturation identity "
            "S tensor V24 = 3 times every nonlinear irrep survives exactly, but "
            "the earlier affine-E7/McKay graph claim does not."
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
        "spin2a_directed_edges": [list(x) for x in directed_edges(A)],
        "spin2b_directed_edges": [list(x) for x in directed_edges(B)],
        "legacy_key_note": (
            "Previous certificate keys called these McKay edges. Pass 11025 "
            "corrects that terminology; these are directed GL2(3) tensor-quiver edges."
        ),
        "underlying_undirected_edges": [list(x) for x in underlying_edges(A)],
        "irrep_dimension_vector": list(map(int, dims)),
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
            "dimension_identity":
                "A*d = 2*d",
            "second_order":
                "A*(A+I)*m = 3*d",
        },
        "interpretation": (
            "Tensoring by either faithful two-dimensional GL2(3) irrep flips "
            "central parity. The clock carrier is balanced so that each 12D "
            "parity half maps to three uniform copies of the three nontrivial "
            "irreducible types of the opposite central parity. This saturation "
            "is a representation-ring identity and does not require an E7 graph."
        ),
        "terminology_correction": (
            "The legacy filename retains 'mckay_e7' for provenance. With the "
            "cyclotomic order-8 character values corrected to +/-i*sqrt(2), the "
            "faithful tensor matrix is directed and non-symmetric; its underlying "
            "undirected graph has 11 edges, not the 7-edge affine-E7 tree. "
            "GL2(3)=2^+S4 is also not the binary-octahedral 2^-S4 subgroup of SL2(C)."
        ),
        "boundary": (
            "This is an identity in the complex representation ring of GL2(3). "
            "The surviving saturation identities do not define an affine-E7 graph, "
            "a classical SU(2) McKay correspondence, or an E7 gauge theory."
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
        "edges": len(p["spin2a_directed_edges"]),
        "full_identity": p["exact_identities"]["full"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
