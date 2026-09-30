#!/usr/bin/env python3
"""Pass 11193: an optimal two-qutrit Clifford compiler from one W33 perfect gate.

Pass 11177 identifies the three double cosets for the action of Sp(4,3) on
the 45 W33-compatible tensor factorisations.  If H is the stabiliser of one
factorisation and p is the symplectic tetracode representative, then

    H                  local/swap gates,
    H p H              perfect gates,
    H p h_* p H        partially entangling gates

for one fixed h_* in H.  Thus every two-qutrit Clifford uses at most two
copies of the same perfect gate p, and the three classes have optimal perfect
gate depths 0, 1, and 2.  This file constructs and verifies every word.

The result is a Clifford compiler.  Approximate universality enters only after
adjoining the already-certified E6 cubic |0>-controlled-X primitive of Pass
10944 (Roy--van de Wetering--Yeh, arXiv:2307.10095).

Tetracode provenance is cited rather than reclaimed: Pass 10946 identifies the
clock evaluation code, Passes 10954--10955 its clock/spin carriers, Pass 10970
its projective code/lattice tower, and Pass 11091 separates its orientation
character.  BT927 and Passes 9173--9196 delimit the distinct E8 glue problem.
Pass 11035 supplies a separate cocycle-to-W33 descent, while Pass 11164 treats
the distinct three-qutrit F9 tetracode analogue.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import w33_pass11156_perfect_spacetime_gates as PG  # noqa: E402
import w33_pass11177_perfect_gates_tritangent as TG  # noqa: E402

OUT = ROOT / "data" / "w33_pass11193_optimal_perfect_gate_compiler.json"


def key(M: np.ndarray) -> tuple[int, ...]:
    return tuple(int(x) for x in (np.asarray(M, dtype=np.int64) % 3).ravel())


def mat(k: tuple[int, ...]) -> np.ndarray:
    return np.array(k, dtype=np.int64).reshape(4, 4)


def det2(M: np.ndarray) -> int:
    M = np.asarray(M, dtype=np.int64) % 3
    return int(M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0]) % 3


def relation(S: np.ndarray) -> str:
    """Rank-three relation relative to the standard two-qutrit split."""
    S = np.asarray(S, dtype=np.int64) % 3
    blocks = [[S[0:2, 0:2], S[0:2, 2:4]],
              [S[2:4, 0:2], S[2:4, 2:4]]]
    zero = [[not b.any() for b in row] for row in blocks]
    block_monomial = (
        (not zero[0][0] and zero[0][1] and zero[1][0] and not zero[1][1])
        or (zero[0][0] and not zero[0][1] and not zero[1][0] and zero[1][1])
    )
    if block_monomial:
        return "local"
    if all(PG.rank3(b) == 2 for row in blocks for b in row):
        return "perfect"
    return "partial"


def local_abi(S: np.ndarray) -> dict:
    """Write an element of H as local SL(2,3)^2 with an optional qutrit swap."""
    S = np.asarray(S, dtype=np.int64) % 3
    Z = np.zeros((2, 2), dtype=np.int64)
    W = np.block([[Z, np.eye(2, dtype=np.int64)],
                  [np.eye(2, dtype=np.int64), Z]])
    if not S[0:2, 2:4].any() and not S[2:4, 0:2].any():
        swapped = False
        A, B = S[0:2, 0:2], S[2:4, 2:4]
        rebuilt = np.block([[A, Z], [Z, B]])
        expression = "diag(A,B)"
    else:
        swapped = True
        A, B = S[0:2, 2:4], S[2:4, 0:2]
        rebuilt = np.block([[A, Z], [Z, B]]) @ W
        expression = "diag(A,B) SWAP"
    assert det2(A) == det2(B) == 1
    assert np.array_equal(rebuilt % 3, S)
    return {
        "expression": expression,
        "swap": swapped,
        "A": A.tolist(),
        "B": B.tolist(),
    }


def fixed_perfect_gate() -> np.ndarray:
    """Pass 11156's tetracode similitude followed by local time reversal."""
    tetracode = np.kron(np.array([[1, 1], [1, 2]], dtype=np.int64),
                          np.eye(2, dtype=np.int64)) % 3
    local_time_reversal = np.diag([1, 2, 1, 2]).astype(np.int64)
    p = (tetracode @ local_time_reversal) % 3
    assert np.array_equal((p.T @ TG.J @ p) % 3, TG.J % 3)
    assert relation(p) == "perfect"
    return p


def double_coset(H: list[np.ndarray], seed: np.ndarray) -> dict[tuple[int, ...], tuple[int, int]]:
    """Enumerate H seed H and retain one deterministic pair of local factors."""
    out: dict[tuple[int, ...], tuple[int, int]] = {}
    for i, L in enumerate(H):
        left = (L @ seed) % 3
        for j, R in enumerate(H):
            out.setdefault(key(left @ R), (i, j))
    return out


def graph_hecke_certificate() -> dict:
    """Verify the rank-three Hecke identity on all 45 W33 factorisations."""
    _, facs = TG.factorisations()
    octets = [frozenset().union(*f) for f in facs]
    n = len(facs)
    A = np.array([[int(i != j and not (octets[i] & octets[j]))
                   for j in range(n)] for i in range(n)], dtype=np.int64)
    I = np.eye(n, dtype=np.int64)
    J = np.ones((n, n), dtype=np.int64)
    B = J - I - A
    lhs = A @ A
    rhs = 12 * I + 3 * A + 3 * B
    assert np.array_equal(lhs, rhs)
    eig = np.rint(np.linalg.eigvalsh(A)).astype(int)
    spectrum = {int(x): int((eig == x).sum()) for x in sorted(set(eig.tolist()))}
    assert spectrum == {-3: 24, 3: 20, 12: 1}
    return {
        "factorisations": n,
        "relations": {"identity": 1, "perfect": 12, "partial": 32},
        "identity": "A^2 = 12 I + 3 A + 3 (J-I-A)",
        "spectrum": {str(k): v for k, v in spectrum.items()},
        "two_perfect_paths_per_target": {
            "identity_target": 12,
            "perfect_target": 3,
            "partial_target": 3,
        },
        "uniform_two_step_probabilities": {
            "local": str(Fraction(12, 12**2)),
            "perfect": str(Fraction(12 * 3, 12**2)),
            "partial": str(Fraction(32 * 3, 12**2)),
        },
        "diameter": 2,
    }


def yukawa_sector_walk_certificate() -> dict:
    """Quotient the perfect-gate walk by Pass 11185's E6 Yukawa sectors.

    Relative to one complete factorisation frame, the 45 tritangent planes
    split as five 1.10.10 and forty 16.16.10 cubic monomials.  This verifies
    that the perfect-relation walk is lumpable on that 5+40 partition and
    freezes its exact two-state dynamics.
    """
    _, facs = TG.factorisations()
    octets = [frozenset().union(*f) for f in facs]
    n = len(facs)
    A = np.array([[int(i != j and not (octets[i] & octets[j]))
                   for j in range(n)] for i in range(n)], dtype=np.int64)
    frames = [c for c in itertools.combinations(range(n), 5)
              if all(A[a, b] for a, b in itertools.combinations(c, 2))]
    assert len(frames) == 27
    phi = frozenset(frames[0])
    sector = ["1.10.10" if i in phi else "16.16.10" for i in range(n)]
    assert Counter(sector) == {"1.10.10": 5, "16.16.10": 40}

    labels = ["1.10.10", "16.16.10"]
    quotient_rows = []
    for src in labels:
        profiles = set()
        for i in range(n):
            if sector[i] != src:
                continue
            profiles.add(tuple(sum(int(A[i, j]) for j in range(n)
                                   if sector[j] == dst) for dst in labels))
        assert len(profiles) == 1
        quotient_rows.append(list(profiles.pop()))
    assert quotient_rows == [[4, 8], [1, 11]]

    # Dividing by the valency 12 gives an exact Markov chain.  Its stationary
    # distribution is simply the sector census, and the only transient mode
    # decays by 1/4 per perfect step.
    P = np.array(quotient_rows, dtype=np.float64) / 12.0
    eig = sorted(np.linalg.eigvals(P).real.tolist(), reverse=True)
    assert np.allclose(eig, [1.0, 0.25])
    return {
        "source": "Pass 11185 E6 Yukawa-sector dictionary",
        "sector_sizes": {"1.10.10": 5, "16.16.10": 40},
        "perfect_neighbour_quotient": quotient_rows,
        "transition_matrix": [["1/3", "2/3"], ["1/12", "11/12"]],
        "stationary_distribution": ["1/9", "8/9"],
        "eigenvalues": ["1", "1/4"],
        "closed_form": {
            "start_1.10.10__probability_1.10.10_after_k": "1/9 + (8/9)(1/4)^k",
            "start_16.16.10__probability_1.10.10_after_k": "1/9 - (1/9)(1/4)^k",
        },
        "full_45_state_walk_eigenvalues": {"1": 1, "1/4": 20, "-1/4": 24},
        "scope": (
            "This is an exact walk on cubic monomial labels induced by uniformly random local dressing of p. "
            "It supplies no Yukawa amplitudes, masses, Hamiltonian, or physical probability law."
        ),
    }


def compile_all() -> dict:
    group = sorted(PG.sp43(), key=key)
    H = sorted((S for S in group if relation(S) == "local"), key=key)
    p = fixed_perfect_gate()
    assert len(group) == 51840 and len(H) == 1152 and PG.order(p) == 4
    assert all(relation(h) == "local" for h in H)

    depth0 = {key(h): i for i, h in enumerate(H)}
    depth1 = double_coset(H, p)
    assert len(depth1) == 13824
    assert all(relation(mat(k)) == "perfect" for k in depth1)

    # A single fixed middle local factor reaches the third double coset.
    h_star_index = next(i for i, h in enumerate(H) if relation(p @ h @ p) == "partial")
    h_star = H[h_star_index]
    q = (p @ h_star @ p) % 3
    depth2 = double_coset(H, q)
    assert len(depth2) == 36864
    assert all(relation(mat(k)) == "partial" for k in depth2)

    K0, K1, K2 = set(depth0), set(depth1), set(depth2)
    assert not (K0 & K1 or K0 & K2 or K1 & K2)
    assert K0 | K1 | K2 == {key(S) for S in group}

    reconstructions = 0
    depth_bytes = bytearray()
    for S in group:
        k = key(S)
        if k in depth0:
            word = [H[depth0[k]]]
            depth = 0
        elif k in depth1:
            i, j = depth1[k]
            word = [H[i], p, H[j]]
            depth = 1
        else:
            i, j = depth2[k]
            word = [H[i], p, h_star, p, H[j]]
            depth = 2
        product = np.eye(4, dtype=np.int64)
        for g in word:
            product = (product @ g) % 3
        assert np.array_equal(product, S)
        reconstructions += 1
        depth_bytes.extend(bytes(k)); depth_bytes.append(depth)

    # The usual qutrit SUM symplectic is partial, so its depth-two word is optimal.
    SUM = np.array([[1, 0, 0, 0],
                    [0, 1, 1, 0],
                    [0, 0, 1, 0],
                    [1, 0, 0, 1]], dtype=np.int64) % 3
    assert relation(SUM) == "partial" and key(SUM) in depth2
    li, ri = depth2[key(SUM)]
    sum_word = [H[li], p, h_star, p, H[ri]]
    prod = np.eye(4, dtype=np.int64)
    for g in sum_word:
        prod = (prod @ g) % 3
    assert np.array_equal(prod, SUM)

    middle_fusion = Counter(relation((p @ h @ p) % 3) for h in H)
    assert middle_fusion == {"local": 96, "perfect": 288, "partial": 768}
    one_qutrit_blocks = set()
    swap_count = 0
    for h in H:
        abi = local_abi(h)
        one_qutrit_blocks.add(tuple(np.array(abi["A"]).ravel()))
        one_qutrit_blocks.add(tuple(np.array(abi["B"]).ravel()))
        swap_count += int(abi["swap"])
    assert len(one_qutrit_blocks) == 24 and swap_count == 576

    pass10944 = json.loads((ROOT / "data/w33_pass10944_five_front_computational_closure.json").read_text())
    universal = pass10944["front2_reversible_cubic_clock"]["universality"]
    assert universal["status"] == "APPROXIMATELY_UNIVERSAL_WITH_EXISTING_QUTRIT_FOURIER_GATE"

    return {
        "group": {
            "name": "Sp(4,3) two-qutrit Clifford symplectic quotient",
            "order": len(group),
            "projective_order": len(group) // 2,
            "local_factorisation_stabiliser_H_order": len(H),
            "H_structure": "(SL(2,3) x SL(2,3)) : C2",
            "one_qutrit_SL23_factors": len(one_qutrit_blocks),
            "H_elements_using_SWAP": swap_count,
        },
        "fixed_gate": {
            "name": "p = tetracode similitude times local time reversal",
            "matrix": p.tolist(),
            "order": PG.order(p),
            "relation": relation(p),
            "four_qutrit_Choi_type": "AME(4,3) perfect tensor",
        },
        "normal_form": {
            "theorem": "Sp(4,3) = H disjoint_union H p H disjoint_union H p h_* p H",
            "h_star": h_star.tolist(),
            "h_star_local_ABI": local_abi(h_star),
            "p_h_star_p_relation": relation(q),
            "optimal_depth_counts": {"0": len(K0), "1": len(K1), "2": len(K2)},
            "optimal_projective_depth_counts": {"0": len(K0) // 2, "1": len(K1) // 2, "2": len(K2) // 2},
            "maximum_perfect_gate_depth": 2,
            "mean_perfect_gate_depth": str(Fraction(len(K1) + 2 * len(K2), len(group))),
            "all_group_elements_reconstructed": reconstructions,
            "depth_table_sha256": hashlib.sha256(depth_bytes).hexdigest(),
            "minimality": (
                "Depth zero is exactly H; depth at most one is exactly H union HpH. "
                "The partial double coset is disjoint from both and therefore requires exactly two p gates."
            ),
        },
        "hecke_fusion": {
            "p_H_p_middle_factor_counts": dict(sorted(middle_fusion.items())),
            "p_H_p_middle_factor_fractions": {
                k: str(Fraction(v, len(H))) for k, v in sorted(middle_fusion.items())
            },
            "graph": graph_hecke_certificate(),
        },
        "E6_yukawa_sector_walk": yukawa_sector_walk_certificate(),
        "SUM_compilation": {
            "target_matrix": SUM.tolist(),
            "optimal_perfect_gate_depth": 2,
            "word": "L p h_* p R",
            "L": H[li].tolist(),
            "L_local_ABI": local_abi(H[li]),
            "h_star": h_star.tolist(),
            "h_star_local_ABI": local_abi(h_star),
            "R": H[ri].tolist(),
            "R_local_ABI": local_abi(H[ri]),
            "product_verified": True,
        },
        "universal_computation_boundary": {
            "Clifford_result": (
                "H and the single fixed perfect gate p generate every two-qutrit Clifford, with optimal p-depth <= 2."
            ),
            "non_Clifford_resource": "Pass 10944 E6_CUBIC_TICK / exact |0>-controlled qutrit X",
            "status_with_existing_E6_cubic_and_Fourier": universal["status"],
            "literature_theorem": universal["literature_theorem"],
            "boundary": (
                "The perfect gate alone remains Clifford and is not computationally universal. "
                "No physical Hamiltonian, error rate, or device implementation is inferred."
            ),
        },
        "repo_prior_art": [
            "analysis/w33_pass11156_perfect_spacetime_gates.py",
            "analysis/PASS11177_PERFECT_GATES_TRITANGENT.md",
            "analysis/w33_pass11180_mereology.py",
            "analysis/PASS11185_E6_YUKAWA_SECTORS.md",
            "analysis/PASS10944_RESERVATION.md",
            "analysis/PASS10946_CLOCK_CODE_CONE_OBJECTWISE.md",
            "analysis/PASS10954_REGULAR_C8_CLOCK_COMPLETION.md",
            "analysis/PASS10955_D4_HALFSPIN_CLOCK_BRIDGE.md",
            "analysis/PASS10970_PROJECTIVE_CLOCK_CODE_LATTICE_TOWER.md",
            "analysis/PASS11091_ONE_GLOBAL_ORIENTATION_BIT.md",
            "analysis/PASS11035_CUBIC_COCYCLE_WEYL_W33_DESCENT.md",
            "analysis/PASS11164_THREE_QUTRIT_TETRACODE.md",
            "analysis/BT927_e8_lift_artifact_reconciliation.md",
            "analysis/PASS9173_9196_RESERVATION.md",
        ],
        "external_checks": [
            "Ian Tan, arXiv:2601.19677: LU uniqueness and local symmetries of AME(4,3)",
            "Roy, van de Wetering, Yeh, arXiv:2307.10095: odd-prime |0>-controlled-X plus Fourier universality",
        ],
    }


def main() -> dict:
    out = {
        "schema": "w33.pass11193.optimal_perfect_gate_compiler.v1",
        "status": "PASS_OPTIMAL_W33_PERFECT_GATE_DEPTH_TWO_COMPILER",
        **compile_all(),
    }
    OUT.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": out["status"],
        "optimal_depth_counts": out["normal_form"]["optimal_depth_counts"],
        "maximum_depth": out["normal_form"]["maximum_perfect_gate_depth"],
        "SUM_depth": out["SUM_compilation"]["optimal_perfect_gate_depth"],
        "hecke_identity": out["hecke_fusion"]["graph"]["identity"],
    }, indent=2))
    return out


if __name__ == "__main__":
    main()
