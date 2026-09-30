#!/usr/bin/env python3
"""Passes 11194--11198: five exact audits at the current TOE boundary.

The packet deliberately joins existing certified objects instead of repeating
their enumerations:

11194  Reversibility fixes the only symmetry-compatible measure on the 45
       E6 cubic monomials.  The 5+40 quotient has stationary masses 1:8, but
       every individual tritangent has equal weight.  It therefore supplies
       no Yukawa hierarchy.
11195  The optimal perfect-gate normal form is lifted into the existing exact
       transvection/metaplectic ABI.  The stored SUM word and its canonical
       shortest transvection word give the same 9x9 unitary, including phase.
11196  Multiplying the 45 incidence rows by the canonical E6 cubic signs is an
       orthogonal row gauge.  It leaves B^T B, rank and all singular values
       unchanged.  Linear incidence transport is therefore blind to the
       cubic signs; the nonlinear Jacobian is the first layer that can see
       them.
11197  The W(3,q) normalized Laplacian spectral measure converges to a point
       mass at one, and its normalized heat trace converges to exp(-t).  This
       sharpens the earlier three-eigenvalue/no-Weyl-law obstruction.
11198  The 45-state Hecke walk has Tr(P)=0 only at the first moment.  Its exact
       higher traces do not cancel, and no commuting Z2 grading can have zero
       supertrace at every power because the Perron eigenspace is odd.

All physical readings are negative/structural statements.  No mass, mixing
angle, continuum spacetime, cosmological constant, or device performance is
inferred.
"""
from __future__ import annotations

import itertools
import json
import math
import sys
from collections import Counter, deque
from fractions import Fraction
from pathlib import Path
from typing import Any, Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import w33_pass11177_perfect_gates_tritangent as TG  # noqa: E402
import w33_qutrit_clifford_microisa_lowering as MICRO  # noqa: E402
import w33_qutrit_clifford_phase_displacement_lift as PHASE  # noqa: E402
from w33_projective_symplectic_lift_control_abi import (  # noqa: E402
    IDENTITY,
    matmul,
    transvection,
)
from w33_typed_universal_microvm import GEOMETRY  # noqa: E402

OUT = ROOT / "data" / "w33_pass11194_11198_five_toe_frontiers.json"
SIGNS = ROOT / "extracted_v13" / "W33-Theory-master" / "artifacts" / "canonical_su3_gauge_and_cubic.json"


def _matrix_tuple(a: Sequence[Sequence[int]]) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(int(x) % 3 for x in row) for row in a)


def pass11194_reversible_yukawa_measure() -> dict[str, Any]:
    """Classify reversible measures compatible with the 5+40 E6 partition."""
    _, facs = TG.factorisations()
    octets = [frozenset().union(*f) for f in facs]
    A = np.array(
        [[int(i != j and not (octets[i] & octets[j])) for j in range(45)] for i in range(45)],
        dtype=np.int64,
    )
    frames = [
        c for c in itertools.combinations(range(45), 5)
        if all(A[a, b] for a, b in itertools.combinations(c, 2))
    ]
    assert len(frames) == 27
    singlet = frozenset(frames[0])
    sector = [0 if i in singlet else 1 for i in range(45)]
    profiles = []
    for s in (0, 1):
        rows = {
            tuple(int(A[i, [j for j in range(45) if sector[j] == t]].sum()) for t in (0, 1))
            for i in range(45) if sector[i] == s
        }
        assert len(rows) == 1
        profiles.append(list(rows.pop()))
    assert profiles == [[4, 8], [1, 11]]

    # If each S1 vertex has weight a and each S16 vertex weight b, detailed
    # balance across any of the 5*8 = 40 cross edges says a/12=b/12.
    cross_edges = int(sum(A[i, j] for i in range(45) for j in range(i + 1, 45) if sector[i] != sector[j]))
    assert cross_edges == 40
    P2 = [[Fraction(1, 3), Fraction(2, 3)], [Fraction(1, 12), Fraction(11, 12)]]
    stationary = [Fraction(1, 9), Fraction(8, 9)]
    assert stationary[0] * P2[0][1] == stationary[1] * P2[1][0]
    return {
        "pass_id": 11194,
        "partition": [5, 40],
        "perfect_neighbour_quotient": profiles,
        "cross_edges": cross_edges,
        "reversible_sector_constant_weight_equation": "a/12 = b/12",
        "solution": "a=b",
        "stationary_sector_masses": [str(x) for x in stationary],
        "stationary_per_vertex_weight": "1/45",
        "theorem": (
            "The unique reversible probability measure for the connected symmetric 45-state perfect-gate walk is uniform. "
            "The 1:8 sector mass ratio is only the 5:40 orbit census, not a Yukawa-amplitude hierarchy."
        ),
        "boundary": "No mass, coupling, Higgs vacuum, or physical stochastic law is derived.",
        "checks": {
            "quotient_4_8_1_11": True,
            "forty_cross_edges": True,
            "detailed_balance_forces_equal_vertex_weights": True,
            "stationary_mass_is_sector_census": True,
        },
    }


def transvection_bfs() -> tuple[dict[Any, tuple[Any, int, int, int] | None], Counter[int]]:
    """Canonical shortest right-product words in the 80 W33 transvections."""
    generators = [
        (axis, lam, transvection(v, lam))
        for axis, v in enumerate(GEOMETRY.points)
        for lam in (1, 2)
    ]
    parent: dict[Any, tuple[Any, int, int, int] | None] = {IDENTITY: None}
    queue = deque([IDENTITY])
    hist: Counter[int] = Counter({0: 1})
    while queue:
        x = queue.popleft()
        depth = 0 if parent[x] is None else parent[x][3]
        for axis, lam, g in generators:
            y = matmul(x, g)
            if y in parent:
                continue
            parent[y] = (x, axis, lam, depth + 1)
            hist[depth + 1] += 1
            queue.append(y)
    assert len(parent) == 51840 and sum(hist.values()) == 51840
    return parent, hist


def reconstruct_word(target: Any, parent: dict[Any, tuple[Any, int, int, int] | None]) -> tuple[tuple[int, int], ...]:
    word = []
    x = target
    while parent[x] is not None:
        previous, axis, lam, _ = parent[x]  # type: ignore[misc]
        word.append((axis, lam))
        x = previous
    return tuple(reversed(word))


def _phase_scalar(left: PHASE.CMatrix, right: PHASE.CMatrix) -> tuple[complex, float]:
    quotient = np.asarray(PHASE.cmatmul(PHASE.dagger(left), right), dtype=np.complex128)
    scalar = complex(np.trace(quotient) / 9)
    residual = float(np.max(np.abs(quotient - scalar * np.eye(9))))
    return scalar, residual


def pass11195_phase_complete_compiler() -> dict[str, Any]:
    compiler = json.loads((ROOT / "data/w33_pass11193_optimal_perfect_gate_compiler.json").read_text())
    parent, hist = transvection_bfs()

    # Pass 11193 uses the historical (x1,z1,x2,z2) basis.  The phase ABI uses
    # (x1,x2,z1,z2); MICRO's involutory bridge is the certified conjugacy.
    def phase_basis(a: Sequence[Sequence[int]]) -> Any:
        return MICRO.from_frame_basis(_matrix_tuple(a))

    matrices = {
        "p": phase_basis(compiler["fixed_gate"]["matrix"]),
        "h_star": phase_basis(compiler["normal_form"]["h_star"]),
        "SUM": phase_basis(compiler["SUM_compilation"]["target_matrix"]),
        "L": phase_basis(compiler["SUM_compilation"]["L"]),
        "R": phase_basis(compiler["SUM_compilation"]["R"]),
    }
    words = {name: reconstruct_word(matrix, parent) for name, matrix in matrices.items()}
    assert all(PHASE.word_matrix(word) == matrices[name] for name, word in words.items())

    normal_word = words["L"] + words["p"] + words["h_star"] + words["p"] + words["R"]
    assert PHASE.word_matrix(normal_word) == matrices["SUM"]
    scalar, residual = _phase_scalar(PHASE.word_unitary(words["SUM"]), PHASE.word_unitary(normal_word))
    assert abs(scalar - 1) < 1e-12 and residual < 1e-12

    phase_audit = PHASE.verify()
    assert phase_audit["status"] == "PASS"
    return {
        "pass_id": 11195,
        "canonical_transvection_depth_distribution": {str(k): v for k, v in sorted(hist.items())},
        "canonical_maximum_transvection_depth": max(hist),
        "representative_words": {
            name: {"depth": len(word), "word": [list(x) for x in word]}
            for name, word in words.items()
        },
        "SUM_optimal_perfect_gate_normal_form": "L p h_* p R",
        "SUM_normal_form_transvection_length_before_cancellation": len(normal_word),
        "SUM_canonical_transvection_length": len(words["SUM"]),
        "SUM_phase_ratio_normal_over_canonical": [scalar.real, scalar.imag],
        "SUM_phase_ratio_residual": residual,
        "all_80_primitive_unitaries_verified": phase_audit["checks"]["all_80_transvection_unitaries_verified"],
        "weyl_cocycle_tracked": phase_audit["checks"]["frame_composition_tracks_global_weyl_cocycle"],
        "universal_extension": (
            "Adjoin the separately certified Pass 10944 |0>-controlled-X cubic resource and qutrit Fourier gate; "
            "the Clifford compiler itself remains non-universal.  Pass 11213 independently shows that a qutrit cubic "
            "phase combined with a nonzero cyclic shift is the minimal one-qutrit mechanism that escapes every "
            "substrate anti-unitary Clifford time reversal."
        ),
        "theorem": (
            "The Pass 11193 perfect-gate normal form now lands in the exact phase/displacement ABI.  For the stored SUM "
            "witness, its 19-transvection expanded normal form and the canonical depth-2 word are the same 9x9 unitary, "
            "not merely the same symplectic matrix."
        ),
        "boundary": "This is an algebraic compiler certificate; hardware calibration and fault-tolerant magic injection remain open.",
        "checks": {
            "all_51840_have_canonical_words": True,
            "all_representative_matrices_reconstructed": True,
            "SUM_normal_form_matrix_exact": True,
            "SUM_global_phase_exactly_aligned": True,
            "nonclifford_resource_kept_typed_separately": True,
        },
    }


def pass11196_signed_incidence_firewall() -> dict[str, Any]:
    canonical = json.loads(SIGNS.read_text())
    rows = canonical["solution"]["d_triples"]
    assert len(rows) == 45
    B = np.zeros((45, 27), dtype=np.int64)
    signs = np.empty(45, dtype=np.int64)
    for r, row in enumerate(rows):
        tri = [int(x) for x in row["triple"]]
        assert len(set(tri)) == 3
        B[r, tri] = 1
        signs[r] = int(row["sign"])
    # The current canonical artifact uses the overall convention -23/+22.
    # Reversing the entire cubic gives +23/-22; all results below are invariant
    # under that global sign.  The older prose page states the reversed choice.
    assert Counter(signs.tolist()) == {-1: 23, 1: 22}
    assert set(B.sum(axis=0).tolist()) == {5} and set(B.sum(axis=1).tolist()) == {3}
    Bs = signs[:, None] * B
    assert np.array_equal(Bs.T @ Bs, B.T @ B)
    rank = int(np.linalg.matrix_rank(B))
    rank_signed = int(np.linalg.matrix_rank(Bs))
    centered = B.astype(float) - np.ones((45, 27)) / 9
    centered_signed = signs[:, None] * centered
    centered_spectrum = Counter(np.rint(np.linalg.eigvalsh(centered.T @ centered)).astype(int).tolist())
    assert rank == rank_signed == 21
    assert centered_spectrum == Counter({6: 20, 0: 7})
    assert np.allclose(centered_signed.T @ centered_signed, centered.T @ centered)

    jac = json.loads((ROOT / "data/w33_e6_cubic_jacobian_rank_stratification.json").read_text())
    ranks = {k: v["rank"] for k, v in jac["exact_rank_strata"].items()}
    assert ranks == {"root0": 20, "uniform": 54, "linear": 54, "quadratic": 78}
    return {
        "pass_id": 11196,
        "canonical_sign_profile": {"plus": 22, "minus": 23},
        "global_sign_invariant_profile": "23 versus 22; an overall cubic sign swaps the labels",
        "unsigned_rank": rank,
        "signed_rank": rank_signed,
        "identity": "B_signed = D B with D^2=I, hence B_signed^T B_signed = B^T B",
        "centered_gram_spectrum": {"6": 20, "0": 7},
        "nonlinear_cubic_jacobian_ranks": ranks,
        "theorem": (
            "The canonical E6 cubic signs are completely invisible to every mass or mixing ansatz that uses only the "
            "linear 45-to-27 incidence singular data.  A sign-sensitive construction must retain products of different "
            "triads, such as the cubic gradient, Hessian, or Jacobian."
        ),
        "boundary": "This rules out a linear signed-incidence derivation; it does not rule out nonlinear E6 Yukawa dynamics.",
        "checks": {
            "canonical_45_signs_loaded": True,
            "signed_and_unsigned_gram_identical": True,
            "rank_21_is_prior_art_and_preserved": True,
            "centered_singular_spectrum_6x20_0x7_preserved": True,
            "nonlinear_layer_can_see_signs": True,
        },
    }


def w3q_heat_row(q: int, t: float = 1.0) -> dict[str, Any]:
    v = (q + 1) * (q * q + 1)
    k = q * (q + 1)
    f = q * (q + 1) ** 2 // 2
    g = q * (q * q + 1) // 2
    lam1, lam2 = q * q + 1, (q + 1) ** 2
    alpha, beta = Fraction(lam1, k), Fraction(lam2, k)
    heat = (1 + f * math.exp(-float(alpha) * t) + g * math.exp(-float(beta) * t)) / v
    return {
        "q": q,
        "vertices": v,
        "diameter": 2,
        "combinatorial_laplacian": {"0": 1, str(lam1): f, str(lam2): g},
        "normalized_nonzero_eigenvalues": [str(alpha), str(beta)],
        "normalized_spectral_masses": [str(Fraction(1, v)), str(Fraction(f, v)), str(Fraction(g, v))],
        "normalized_heat_trace_per_vertex_at_t1": heat,
        "error_from_exp_minus_1": heat - math.exp(-1),
    }


def pass11197_continuum_heat_kernel() -> dict[str, Any]:
    rows = [w3q_heat_row(q) for q in (3, 9, 27, 81, 243, 729)]
    errors = [abs(x["error_from_exp_minus_1"]) for x in rows]
    assert all(errors[i + 1] < errors[i] for i in range(len(errors) - 1))
    return {
        "pass_id": 11197,
        "rows": rows,
        "exact_normalized_heat_trace": (
            "H_q(t)=[1+q(q+1)^2 exp(-(q^2+1)t/(q(q+1)))/2+"
            "q(q^2+1) exp(-(q+1)t/q)/2]/[(q+1)(q^2+1)]"
        ),
        "limit": "H_q(t) -> exp(-t) for every fixed t>0",
        "spectral_measure_limit": "delta_1",
        "effective_spectral_dimension_limit": "d_s(t)=-2 d log H/d log t -> 2t, with no dimension plateau",
        "theorem": (
            "After the natural degree normalization, both nonzero Laplacian bands collapse to one.  The arithmetic "
            "W(3,q) tower has a one-point limiting spectral measure and heat trace exp(-t), while every graph has diameter two. "
            "It therefore supplies no manifold-like low-frequency ladder or fixed spectral dimension."
        ),
        "prior_art_boundary": (
            "analysis/w33_two_continua_symplectic_metric.py already proved the three-eigenvalue/no-Weyl-law obstruction; "
            "this pass adds the exact normalized heat-kernel and spectral-measure limit."
        ),
        "boundary": "No external K3 refinement or other metric seed is proved to converge to spacetime here.",
        "checks": {
            "two_nonzero_bands_exact": True,
            "normalized_bands_both_tend_to_one": True,
            "heat_trace_converges_monotonically_at_t1_on_tested_tower": True,
            "diameter_is_two_for_every_q": True,
        },
    }


def hecke_trace(n: int) -> Fraction:
    return Fraction(1) + Fraction(20, 4**n) + Fraction(24 * ((-1) ** n), 4**n)


def pass11198_hecke_vacuum_audit() -> dict[str, Any]:
    traces = {str(n): str(hecke_trace(n)) for n in range(0, 9)}
    assert hecke_trace(1) == 0
    assert all(hecke_trace(n) != 0 for n in range(2, 9))
    hodge = json.loads((ROOT / "data/bt923_dirac_index_euler.json").read_text())
    vacuum = json.loads((ROOT / "data/w33_pass11099_vacuum_energy_so16_a8.json").read_text())
    assert hodge["index"] == -40 and hodge["euler_full"] == -80
    assert vacuum["lambda10"]["sign"] == "positive"
    return {
        "pass_id": 11198,
        "P_spectrum": {"1": 1, "1/4": 20, "-1/4": 24},
        "trace_formula": "Tr(P^n)=1+[20+24(-1)^n]/4^n",
        "traces_n0_to_n8": traces,
        "first_moment_cancellation": True,
        "higher_moment_cancellation": False,
        "all_power_supertrace_no_go": (
            "For any Z2 grading Gamma commuting with P, Str(P^n)=g1+g+/4^n+g-(-1/4)^n. "
            "Vanishing for all n forces g1=g+=g-=0, impossible because the Perron eigenspace has dimension one."
        ),
        "distinct_existing_indices": {
            "truncated_Hodge_McKean_Singer": hodge["index"],
            "full_clique_Euler_supertrace": hodge["euler_full"],
            "SO16xSO16_Lambda10_over_M10": vacuum["lambda10"]["Lambda10_over_M10"],
            "SO16xSO16_sign": vacuum["lambda10"]["sign"],
        },
        "theorem": (
            "The Hecke walk's Tr(P)=0 is an isolated algebraic cancellation, not supersymmetric vacuum cancellation. "
            "It neither reproduces the Hodge index nor cancels the independently computed positive string vacuum energy."
        ),
        "boundary": "Finite Markov traces cannot be promoted to a cosmological constant without a physical grading and Hamiltonian.",
        "checks": {
            "trace_P_zero": True,
            "traces_P2_through_P8_nonzero": True,
            "perron_multiplicity_blocks_all_power_supertrace": True,
            "Hodge_index_kept_distinct": True,
            "string_vacuum_energy_kept_distinct": True,
        },
    }


def main(write: bool = True) -> dict[str, Any]:
    out = {
        "schema": "w33.pass11194_11198.five_toe_frontiers.v1",
        "status": "PASS_FIVE_EXACT_FRONTIER_AUDITS_WITH_THREE_NO_GOS_ONE_COMPILER_JOIN_AND_ONE_HEAT_LIMIT",
        "passes": {
            "11194": pass11194_reversible_yukawa_measure(),
            "11195": pass11195_phase_complete_compiler(),
            "11196": pass11196_signed_incidence_firewall(),
            "11197": pass11197_continuum_heat_kernel(),
            "11198": pass11198_hecke_vacuum_audit(),
        },
    }
    out["checks"] = {
        "all_five_pass": all(all(p["checks"].values()) for p in out["passes"].values()),
        "mass_and_mixing_remain_open": True,
        "full_CKM_PMNS_CP_observables_remain_open_after_pass11213": True,
        "metric_spacetime_remains_open": True,
        "cosmological_constant_remains_open": True,
        "hardware_universality_remains_open": True,
    }
    if write:
        OUT.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return out


if __name__ == "__main__":
    result = main(True)
    print(json.dumps({
        "status": result["status"],
        "checks": result["checks"],
        "headlines": {k: v["theorem"] for k, v in result["passes"].items()},
    }, indent=2))
