#!/usr/bin/env python3
"""Pass 11540: identify the 27-history null graph with VO(3,3).

The verifier is deliberately standard-library only.  It rebuilds the ternary
quadratic history space from scratch and exhaustively enumerates its finite
orthogonal group, so the group identifications do not depend on GAP/Sage.

Exact finite statements verified here:
  * history adjacency is q(x-y)=0 on F_3^3, hence the parabolic affine-polar
    graph VO(3,3);
  * O(q), SO(q), and Omega(q) have orders 48, 24, and 12;
  * SO(q) acts faithfully as S4 on the four projective null directions;
  * its even-permutation kernel is the derived subgroup A4 = Omega(3,3);
  * the affine orders are 1296, 648, 324;
  * the two independent C2 characters (orthogonal determinant and null-frame
    parity / spinor-norm kernel) give a C2 x C2 orientation square.

The standard classical-group name "spinor norm" is external mathematical
classification: for odd q, Omega(n,q) is the kernel of the spinor norm on
SO(n,q).  This script verifies the exact kernel subgroup without importing a
spinor-norm implementation.

Physical firewall: no continuum Lorentz group, speed of light, gravity, or
space-time dynamics is inferred.
"""

from __future__ import annotations

import itertools
import json
from collections import Counter, deque
from pathlib import Path

P = 3
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "PART_W33_PASS11540_HISTORY_VO33_ORIENTATION_SQUARE.json"

V = list(itertools.product(range(P), repeat=3))
ZERO = (0, 0, 0)
I3 = (1, 0, 0, 0, 1, 0, 0, 0, 1)


def mod(x: int) -> int:
    return x % P


def add(u, v):
    return tuple(mod(u[i] + v[i]) for i in range(3))


def sub(u, v):
    return tuple(mod(u[i] - v[i]) for i in range(3))


def neg(u):
    return tuple(mod(-x) for x in u)


def q(u) -> int:
    # Diagonal form equivalent over F_3 to det([[t+x,y],[y,t-x]]).
    return mod(u[0] * u[0] - u[1] * u[1] - u[2] * u[2])


def dot(u, v) -> int:
    return mod(sum(u[i] * v[i] for i in range(3)))


def det3(M) -> int:
    return mod(
        M[0] * (M[4] * M[8] - M[5] * M[7])
        - M[1] * (M[3] * M[8] - M[5] * M[6])
        + M[2] * (M[3] * M[7] - M[4] * M[6])
    )


def mv(M, v):
    return (
        mod(M[0] * v[0] + M[1] * v[1] + M[2] * v[2]),
        mod(M[3] * v[0] + M[4] * v[1] + M[5] * v[2]),
        mod(M[6] * v[0] + M[7] * v[1] + M[8] * v[2]),
    )


def mm(A, B):
    C = [0] * 9
    for i in range(3):
        for j in range(3):
            C[3 * i + j] = mod(
                sum(A[3 * i + k] * B[3 * k + j] for k in range(3))
            )
    return tuple(C)


def matrix_order(M, cap=48):
    X = I3
    for n in range(1, cap + 1):
        X = mm(X, M)
        if X == I3:
            return n
    raise AssertionError("matrix order exceeded cap")


def permutation_parity(p):
    inv = sum(p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p)))
    return -1 if inv % 2 else 1


def enumerate_orthogonal_group():
    O = []
    SO = []
    for entries in itertools.product(range(P), repeat=9):
        d = det3(entries)
        if d == 0:
            continue
        if all(q(mv(entries, v)) == q(v) for v in V):
            O.append(tuple(entries))
            if d == 1:
                SO.append(tuple(entries))
    return O, SO


def inverse_in_group(M, G):
    for N in G:
        if mm(M, N) == I3 and mm(N, M) == I3:
            return N
    raise AssertionError("inverse not found")


def subgroup_generated(generators):
    H = {I3}
    queue = deque([I3])
    while queue:
        h = queue.popleft()
        for g in generators:
            x = mm(h, g)
            if x not in H:
                H.add(x)
                queue.append(x)
    return H


def projective_null_pairs(nulls):
    unseen = set(nulls)
    pairs = []
    while unseen:
        v = min(unseen)
        w = neg(v)
        pair = tuple(sorted((v, w)))
        pairs.append(pair)
        unseen.remove(v)
        unseen.remove(w)
    return tuple(sorted(pairs))


def projective_perm(M, pairs):
    def pair_index(v):
        for i, pair in enumerate(pairs):
            if v in pair:
                return i
        raise AssertionError("image left null cone")

    return tuple(pair_index(mv(M, pair[0])) for pair in pairs)


def adjacency():
    nulls = [v for v in V if v != ZERO and q(v) == 0]
    A = {v: set() for v in V}
    for x in V:
        for n in nulls:
            A[x].add(add(x, n))
    return A, nulls


def distance_profile(A, start=ZERO):
    d = {start: 0}
    queue = deque([start])
    while queue:
        x = queue.popleft()
        for y in A[x]:
            if y not in d:
                d[y] = d[x] + 1
                queue.append(y)
    return dict(sorted(Counter(d.values()).items()))


def fourier_spectrum(nulls):
    # For F_3 characters, the imaginary parts cancel in +/- pairs.  If
    # c_r(k) is the number of null steps with k.n=r, then
    # lambda(k)=c_0-c_1 because c_1=c_2 and omega+omega^2=-1.
    rows = []
    for k in V:
        c = Counter(dot(k, n) for n in nulls)
        assert c[1] == c[2]
        lam = c[0] - c[1]
        rows.append((k, q(k), lam))
    return rows


def affine_action_is_graph_automorphism(A, O):
    # Every x -> Mx+a preserves q(x-y), but check all 1296 actions explicitly.
    checked = 0
    for M in O:
        for a in V:
            for x in V:
                image_neighbors = {add(mv(M, y), a) for y in A[x]}
                xx = add(mv(M, x), a)
                assert image_neighbors == A[xx]
            checked += 1
    return checked


def main():
    A, nulls = adjacency()
    assert len(V) == 27
    assert len(nulls) == 8
    assert set(len(A[v]) for v in V) == {8}
    edges = sum(map(len, A.values())) // 2
    assert edges == 108
    assert distance_profile(A) == {0: 1, 1: 8, 2: 18}

    fourier = fourier_spectrum(nulls)
    eig = Counter(lam for _, _, lam in fourier)
    assert eig == Counter({2: 12, -1: 8, -4: 6, 8: 1})
    shell_to_eig = {}
    for qq in (0, 1, 2):
        vals = {lam for k, qk, lam in fourier if k != ZERO and qk == qq}
        assert len(vals) == 1
        shell_to_eig[str(qq)] = vals.pop()

    O, SO = enumerate_orthogonal_group()
    assert len(O) == 48
    assert len(SO) == 24

    pairs = projective_null_pairs(nulls)
    assert len(pairs) == 4

    perms_SO = {projective_perm(M, pairs) for M in SO}
    assert len(perms_SO) == 24  # faithful full S4 action

    parity_kernel = {
        M for M in SO if permutation_parity(projective_perm(M, pairs)) == 1
    }
    assert len(parity_kernel) == 12

    inv = {M: inverse_in_group(M, SO) for M in SO}
    commutators = {
        mm(mm(mm(inv[A], inv[B]), A), B)
        for A in SO
        for B in SO
    }
    derived = subgroup_generated(commutators)
    assert len(derived) == 12
    assert derived == parity_kernel

    so_orders = Counter(matrix_order(M) for M in SO)
    omega_orders = Counter(matrix_order(M) for M in derived)
    assert so_orders == Counter({2: 9, 3: 8, 4: 6, 1: 1})
    assert omega_orders == Counter({3: 8, 2: 3, 1: 1})

    char_pairs = Counter()
    for M in O:
        det_bit = "+" if det3(M) == 1 else "-"
        frame_bit = "+" if permutation_parity(projective_perm(M, pairs)) == 1 else "-"
        char_pairs[det_bit + frame_bit] += 1
    assert char_pairs == Counter({"++": 12, "+-": 12, "-+": 12, "--": 12})

    minus_I = (2, 0, 0, 0, 2, 0, 0, 0, 2)
    assert minus_I in O
    assert det3(minus_I) == 2
    assert projective_perm(minus_I, pairs) == tuple(range(4))
    # This proves the determinant and projective-null-frame parity characters
    # are independent: -I flips the former while fixing the latter.

    affine_checked = affine_action_is_graph_automorphism(A, O)
    assert affine_checked == 27 * 48

    old_bigcell = json.loads(
        (ROOT / "data" / "w33_20260924_history_bigcell_q43_compactification.json").read_text()
    )
    old_orientation = json.loads(
        (ROOT / "data" / "PART_W33_PASS9741_9748_ORIENTATION_CHARACTER_WELD.json").read_text()
    )
    old_clock = json.loads(
        (ROOT / "data" / "w33_20260924_null_history_spectral_clock.json").read_text()
    )
    assert old_bigcell["stabilizer_actions"]["PGSp_bell_line_action_order"] == 1296
    assert old_bigcell["stabilizer_actions"]["PSp_bell_line_action_order"] == 648
    assert old_orientation["W33_line_orientation"]["order_kernel"] == 324
    assert old_clock["graph"]["vertices"] == 27
    assert old_clock["graph"]["degree"] == 8

    out = {
        "schema": "w33.pass11540.history_vo33_orientation_square.v1",
        "status": "PASS_EXACT_VO33_ORTHOGONAL_ORIENTATION_SQUARE",
        "pass": 11540,
        "history_graph": {
            "standard_name": "parabolic affine orthogonal polar graph VO(3,3)",
            "carrier": "F3^3",
            "quadratic_form": "q(t,x,y)=t^2-x^2-y^2",
            "adjacency": "x~y iff x!=y and q(x-y)=0",
            "vertices": 27,
            "degree": 8,
            "edges": 108,
            "distance_profile_from_origin": {"0": 1, "1": 8, "2": 18},
            "adjacency_spectrum": {"8": 1, "2": 12, "-1": 8, "-4": 6},
            "fourier_shell_eigenvalue": {
                "zero_character": 8,
                "nonzero_q0": -1,
                "q1": -4,
                "q2": 2,
            },
        },
        "projective_null_boundary": {
            "nonzero_null_vectors": 8,
            "projective_null_directions": [
                [list(v) for v in pair] for pair in pairs
            ],
            "direction_count": 4,
            "identification": "Q(2,3) = P1(F3), the same four-set used by the Hesse/tetracode clock",
        },
        "orthogonal_group": {
            "O_3_3_order": 48,
            "SO_3_3_order": 24,
            "Omega_3_3_order": 12,
            "SO_order_census": dict(sorted(so_orders.items())),
            "Omega_order_census": dict(sorted(omega_orders.items())),
            "SO_on_four_null_directions": "faithful S4",
            "Omega_on_four_null_directions": "A4 = even permutations",
            "Omega_equals_derived_SO": True,
            "standard_isomorphisms": [
                "SO(3,3) ~= PGL(2,3) ~= S4",
                "Omega(3,3) ~= PSL(2,3) ~= A4",
            ],
        },
        "affine_group_ladder": {
            "AffO_order": 27 * 48,
            "AffSO_order": 27 * 24,
            "AffOmega_order": 27 * 12,
            "AffO": "3^3:O(3,3)",
            "AffSO": "3^3:SO(3,3) ~= 3^3:S4",
            "AffOmega": "3^3:Omega(3,3) ~= 3^3:A4",
            "all_1296_affine_isometries_explicitly_checked": True,
        },
        "orientation_square": {
            "character_pair_census": dict(sorted(char_pairs.items())),
            "outer_global_bit": {
                "quotient": "AffO/AffSO ~= C2",
                "character": "determinant of the 3D orthogonal linear part",
                "repo_match": "PGSp Bell stabilizer 1296 -> PSp Bell stabilizer 648; this is the existing global history-cycle/antiunitary/chirality reversal bit",
            },
            "inner_null_frame_bit": {
                "quotient": "AffSO/AffOmega ~= C2",
                "character": "parity of the S4 action on the four projective null directions",
                "standard_classical_name": "spinor norm quotient; Omega(3,3)=ker(spinor norm) inside SO(3,3)",
                "repo_match": "existing 3^3:S4 -> 3^3:A4 line-orientation quotient 648 -> 324",
            },
            "combined_quotient": "AffO/AffOmega ~= C2 x C2",
            "independence_witness": "-I has orthogonal determinant -1 but fixes all four projective null directions",
        },
        "repo_weld": {
            "existing_bigcell_orders_recovered": {"PGSp": 1296, "PSp": 648},
            "existing_line_orientation_kernel_recovered": 324,
            "existing_fourier_clock_spectrum_recovered": True,
            "new_content": "standard VO(3,3) name plus the orthogonal determinant/spinor-norm separation of the two C2 orientation bits",
        },
        "boundary": (
            "Exact finite quadratic-graph and group theorem. The identification of "
            "Omega with the spinor-norm kernel is standard finite orthogonal-group "
            "theory; this verifier independently reconstructs that subgroup as the "
            "derived/even-null-frame kernel. No continuum Lorentz symmetry, physical "
            "speed of light, Einstein dynamics, or gravity is inferred."
        ),
    }

    OUT.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": out["status"],
        "orders": [
            out["affine_group_ladder"]["AffO_order"],
            out["affine_group_ladder"]["AffSO_order"],
            out["affine_group_ladder"]["AffOmega_order"],
        ],
        "characters": out["orientation_square"]["character_pair_census"],
        "quotient": out["orientation_square"]["combined_quotient"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
