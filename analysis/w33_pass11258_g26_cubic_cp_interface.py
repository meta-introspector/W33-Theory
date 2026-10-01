#!/usr/bin/env python3
"""Pass 11258: a conditional cubic-phase CP interface for a G26 control.

The standard G26 target point is normalized as a qutrit control state.  After
the two-qutrit gate (T tensor T) SUM, postselecting the same control state
induces a target Kraus operator.  Its polar unitary has a nonzero, explicitly
rephasing-invariant Jarlskog quartet.  This is an interface theorem, not a CKM
fit: no map from standard G26 coordinates to the physical 81 has been built.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/w33_pass11258_g26_cubic_cp_interface.json"


def sum_gate():
    U = np.zeros((9, 9), complex)
    for a in range(3):
        for b in range(3):
            U[3*a + ((a+b) % 3), 3*a+b] = 1
    return U


def polar_unitary(K):
    left, singular, right = np.linalg.svd(K)
    return left @ right, singular


def jarlskog(W):
    return float(np.imag(W[0, 0] * W[1, 1] * np.conj(W[0, 1]) * np.conj(W[1, 0])))


def induced(point, left_phase=True, right_phase=True):
    zeta9 = np.exp(2j * np.pi / 9)
    T = np.diag([1, zeta9, np.conj(zeta9)])
    eye = np.eye(3, dtype=complex)
    U = np.kron(T if left_phase else eye, T if right_phase else eye) @ sum_gate()
    v = np.array(point, dtype=complex)
    v /= np.linalg.norm(v)
    K = np.einsum("a,atbs,b->ts", v.conj(), U.reshape(3, 3, 3, 3), v)
    W, singular = polar_unitary(K)
    return K, W, singular


def main(write=True):
    K, W, singular = induced((1, 2, 3), True, True)
    J = jarlskog(W)
    controls = {}
    for name, phases in {
        "SUM_only": (False, False),
        "control_T_only": (True, False),
        "target_T_only": (False, True),
        "both_T": (True, True),
    }.items():
        _, unitary, values = induced((1, 2, 3), *phases)
        controls[name] = {"jarlskog": jarlskog(unitary), "singular_values": [float(x) for x in values]}

    _, Wsymmetric, _ = induced((1, 1, 1), True, True)
    rng = np.random.default_rng(11258)
    errors = []
    for _ in range(32):
        rows = np.diag(np.exp(2j*np.pi*rng.random(3)))
        cols = np.diag(np.exp(2j*np.pi*rng.random(3)))
        errors.append(abs(jarlskog(rows @ W @ cols) - J))
    conjugate_error = abs(jarlskog(W.conj()) + J)

    out = {
        "schema": "w33.pass11258.g26-cubic-cp-interface.v1",
        "status": "PASS_CONDITIONAL_CUBIC_PHASE_CP_INTERFACE",
        "control": {
            "standard_g26_point": [1, 2, 3],
            "normalized_as_qutrit_state": True,
            "two_qutrit_gate": "(T tensor T) SUM, T=diag(1,zeta9,zeta9^-1)",
            "reduction": "K_v=<v|U|v>, followed by polar unitary W",
        },
        "result": {
            "kraus_singular_values": [float(x) for x in singular],
            "polar_unitary_absolute_squares": np.abs(W).astype(float).__pow__(2).tolist(),
            "jarlskog": J,
            "symmetric_point_jarlskog": jarlskog(Wsymmetric),
            "conjugate_jarlskog": jarlskog(W.conj()),
            "maximum_32_trial_row_column_rephasing_error": float(max(errors)),
            "conjugation_sign_flip_error": float(conjugate_error),
        },
        "gate_controls": controls,
        "interpretation": {
            "mechanism": "A real asymmetric control amplitude and a cubic phase produce an effective target CP phase after conditioning.",
            "why_control_T_can_survive": "Control postselection converts a phase that would be locally removable before conditioning into relative coefficients of the target Kraus operator.",
            "physical_ckm_claim": False,
            "missing_map": "No certified semisimple map from these standard G26 coordinates into the W33 3 tensor 27 matter carrier exists yet.",
        },
        "checks": {
            "full_rank_kraus": bool(np.min(singular) > 1e-10),
            "cp_nonzero": abs(J) > 1e-6,
            "symmetric_control_cp_zero": abs(jarlskog(Wsymmetric)) < 1e-12,
            "rephasing_invariant": max(errors) < 1e-12,
            "complex_conjugation_flips_sign": conjugate_error < 1e-12,
            "sum_only_cp_zero": abs(controls["SUM_only"]["jarlskog"]) < 1e-12,
            "target_only_cp_zero": abs(controls["target_T_only"]["jarlskog"]) < 1e-12,
        },
    }
    assert all(out["checks"].values())
    if write:
        OUT.write_text(json.dumps(out, indent=2) + "\n")
    return out


if __name__ == "__main__":
    print(json.dumps(main(True), indent=2))
