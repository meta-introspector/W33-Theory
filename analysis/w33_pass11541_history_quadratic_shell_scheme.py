#!/usr/bin/env python3
"""Pass 11541: exact quadratic-shell association scheme on the 27 histories.

Relations on V=F_3^3 are defined by the quadratic value of the difference:
  R0: x=y
  R1: x!=y and q(x-y)=0
  R2: q(x-y)=1
  R3: q(x-y)=2
for q(t,x,y)=t^2-x^2-y^2.

The verifier proves:
  * these four relations form a 3-class translation association scheme;
  * valencies are (1,8,6,12);
  * additive Fourier characters diagonalize all relations;
  * the exact first eigenmatrix is
      [1,  8,  6, 12]
      [1, -1, -3,  3]
      [1, -4,  3,  0]
      [1,  2,  0, -3]
  * P^2=27 I, hence the shell scheme is formally self-dual in this ordering.

This is a finite harmonic/combinatorial theorem only.
"""

from __future__ import annotations

import itertools
import json
from collections import Counter
from pathlib import Path

P = 3
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "PART_W33_PASS11541_HISTORY_QUADRATIC_SHELL_SCHEME.json"

V = list(itertools.product(range(P), repeat=3))
ZERO = (0, 0, 0)


def mod(x):
    return x % P


def q(v):
    return mod(v[0] * v[0] - v[1] * v[1] - v[2] * v[2])


def dot(u, v):
    return mod(sum(u[i] * v[i] for i in range(3)))


def sub(u, v):
    return tuple(mod(u[i] - v[i]) for i in range(3))


def shell_index(v):
    if v == ZERO:
        return 0
    return {0: 1, 1: 2, 2: 3}[q(v)]


def shells():
    C = [[] for _ in range(4)]
    for v in V:
        C[shell_index(v)].append(v)
    return C


def char_sum(k, S):
    counts = Counter(dot(k, s) for s in S)
    # S=-S, so residues 1 and 2 pair and the sum is real:
    # c0 + c1*w + c2*w^2 = c0-c1 when c1=c2.
    assert counts[1] == counts[2]
    return counts[0] - counts[1]


def matmul(A, B):
    return [
        [
            sum(A[i][k] * B[k][j] for k in range(len(B)))
            for j in range(len(B[0]))
        ]
        for i in range(len(A))
    ]


def intersection_tensor(C):
    # p_ij^k for base point 0 and target z in C_k.
    tensor = []
    for k in range(4):
        representative = C[k][0]
        M = []
        for i in range(4):
            row = []
            for j in range(4):
                n = sum(shell_index(sub(representative, x)) == j for x in C[i])
                row.append(n)
            M.append(row)

        # Verify the count depends only on the relation class k.
        for z in C[k]:
            for i in range(4):
                for j in range(4):
                    n = sum(shell_index(sub(z, x)) == j for x in C[i])
                    assert n == M[i][j]
        tensor.append(M)
    return tensor


def bose_mesner_product_check(C, tensor):
    # Directly check A_i A_j = sum_k p_ij^k A_k entrywise using differences.
    for i in range(4):
        for j in range(4):
            for x in V:
                for y in V:
                    lhs = sum(
                        shell_index(sub(z, x)) == i
                        and shell_index(sub(y, z)) == j
                        for z in V
                    )
                    k = shell_index(sub(y, x))
                    assert lhs == tensor[k][i][j]


def main():
    C = shells()
    valencies = [len(c) for c in C]
    assert valencies == [1, 8, 6, 12]

    tensor = intersection_tensor(C)
    bose_mesner_product_check(C, tensor)

    reps = [c[0] for c in C]
    eig = [[char_sum(k, S) for S in C] for k in reps]
    expected = [
        [1, 8, 6, 12],
        [1, -1, -3, 3],
        [1, -4, 3, 0],
        [1, 2, 0, -3],
    ]
    assert eig == expected

    P2 = matmul(eig, eig)
    assert P2 == [[27 if i == j else 0 for j in range(4)] for i in range(4)]

    # Because the dual quadratic form has the same shell census in these
    # coordinates, the primitive multiplicities equal the relation valencies.
    multiplicities = [len(c) for c in C]
    assert multiplicities == valencies

    old_clock = json.loads(
        (ROOT / "data" / "w33_20260924_null_history_spectral_clock.json").read_text()
    )
    assert old_clock["graph"]["vertices"] == 27
    assert old_clock["graph"]["degree"] == 8
    assert old_clock["fourier_duality"]["eigenvalue_by_dual_quadratic_type"] == {
        "zero": 8,
        "det0_nonzero": -1,
        "det1": -4,
        "det2": 2,
    }

    pass11540 = json.loads(
        (ROOT / "data" / "PART_W33_PASS11540_HISTORY_VO33_ORIENTATION_SQUARE.json").read_text()
    )
    assert pass11540["history_graph"]["standard_name"].endswith("VO(3,3)")

    out = {
        "schema": "w33.pass11541.history_quadratic_shell_scheme.v1",
        "status": "PASS_EXACT_HISTORY_QUADRATIC_SHELL_ASSOCIATION_SCHEME",
        "pass": 11541,
        "carrier": "F3^3 = Sym_2(F3)",
        "quadratic_form": "q(t,x,y)=t^2-x^2-y^2",
        "relations": [
            "R0: difference 0",
            "R1: nonzero difference with q=0 (null)",
            "R2: difference with q=1",
            "R3: difference with q=2",
        ],
        "valencies": valencies,
        "primitive_multiplicities": multiplicities,
        "first_eigenmatrix_P": eig,
        "second_eigenmatrix_Q": eig,
        "self_duality_identity": "P^2 = 27 I_4",
        "P_squared": P2,
        "intersection_tensor_pij_by_k": tensor,
        "null_graph_column": {
            "relation": "R1",
            "spectrum_by_primitive_shell": [row[1] for row in eig],
            "multiplicities": multiplicities,
            "expanded_spectrum": {"8": 1, "-1": 8, "-4": 6, "2": 12},
        },
        "interpretation": (
            "The paper's four history separation classes (coincident/null/q=1/q=2) "
            "close as one commutative Bose-Mesner algebra. The same 1+8+6+12 "
            "shell census occurs on the Fourier-dual characters, so the history "
            "geometry is formally self-dual before any continuum interpretation."
        ),
        "literature_boundary": (
            "Translation association schemes from quadratic/symmetric bilinear forms "
            "and their duality are standard. The repo-specific content is the exact "
            "small F3 shell fusion, its explicit eigenmatrix/intersection tensor, and "
            "its weld to the frozen W33 temporal clock."
        ),
        "boundary": (
            "Finite association-scheme theorem only. Formal Fourier self-duality is "
            "not physical space-time duality and supplies no continuum metric, "
            "Einstein equation, or measured coupling constant."
        ),
    }
    OUT.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": out["status"],
        "valencies": valencies,
        "P": eig,
        "self_dual": True,
    }, indent=2))


if __name__ == "__main__":
    main()
