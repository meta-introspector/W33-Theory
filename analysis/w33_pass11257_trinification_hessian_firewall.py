#!/usr/bin/env python3
"""Pass 11257: pull the controlled Hessian through the old kernel, as an audit.

The 81-by-3 matrix below is the Pass 11218 abelian-centralizer basis.  Pass
11255 proves that it is not a semisimple Cartan embedding.  We nevertheless
compute the unique minimum-Frobenius-norm 81-dimensional quadratic form whose
compression is the Pass 11256 Hessian.  It is a useful counterfeit: its
dependence on the chosen Schlaefli reference line quantifies exactly why no
1+10+16 mass reading is licensed by the old data.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/w33_pass11257_trinification_hessian_firewall.json"


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def main(write=True):
    old = load(ROOT / "analysis/w33_pass11218_e8_cubic_cartan_kernel.py", "p11218_for_11257")
    common = load(ROOT / "analysis/w33_pass11255_11259_g26_common.py", "p11255_common_11257")
    geom = load(ROOT / "analysis/w33_pass4992_4999_common.py", "p4992_common_11257")

    parent = old.load_parent()
    records = parent.build_records()
    A = sp.Matrix(parent.jacobian(old.dense_from_sparse(old.BACKGROUND), records))
    kernel = [old.primitive_integer_vector(v) for v in A.nullspace()]
    K = np.array(sp.Matrix.hstack(*(sp.Matrix(v) for v in kernel)), dtype=float)
    _, _, _, Hq = common.target_potential_data((1, 2, 3))
    H = np.array(Hq, dtype=float)

    gram_inverse = np.linalg.inv(K.T @ K)
    pullback = K @ gram_inverse @ H @ gram_inverse @ K.T
    residual = float(np.max(np.abs(K.T @ pullback @ K - H)))

    G = geom.build_base()["G27"]
    profiles = []
    support_audits = []
    for ref in range(27):
        one = [ref]
        ten = sorted(G.neighbors(ref))
        sixteen = sorted(set(range(27)) - set(one) - set(ten))
        assert (len(one), len(ten), len(sixteen)) == (1, 10, 16)
        blocks27 = [one, ten, sixteen]
        blocks81 = [[3*i+a for i in block for a in range(3)] for block in blocks27]
        energies = []
        supports = []
        for j in range(3):
            energies.append([float(np.sum(K[idx, j] ** 2)) for idx in blocks81])
            supports.append([int(np.count_nonzero(K[idx, j])) for idx in blocks81])
        block_traces = [float(np.trace(pullback[np.ix_(idx, idx)])) for idx in blocks81]
        cross = [[float(np.linalg.norm(pullback[np.ix_(a, b)], "fro")) for b in blocks81]
                 for a in blocks81]
        profiles.append({
            "reference_line": ref,
            "block_traces_1_10_16": block_traces,
            "cross_block_frobenius": cross,
        })
        support_audits.append({
            "reference_line": ref,
            "kernel_squared_energy_by_basis_and_1_10_16": energies,
            "kernel_support_by_basis_and_1_10_16": supports,
        })

    keys = {
        tuple(round(x, 10) for x in p["block_traces_1_10_16"])
        for p in profiles
    }
    traces = np.array([p["block_traces_1_10_16"] for p in profiles])
    nonzero = np.linalg.eigvalsh(pullback)
    nonzero = nonzero[nonzero > 1e-9]
    out = {
        "schema": "w33.pass11257.trinification-hessian-firewall.v1",
        "status": "PASS_CONDITIONAL_PULLBACK_COMPUTED_BUT_PHYSICAL_TRINIFICATION_PULLBACK_BLOCKED",
        "construction": {
            "formula": "M=K(K^T K)^(-1) H (K^T K)^(-1)K^T",
            "property": "K^T M K=H",
            "maximum_compression_residual": residual,
            "rank": int(np.linalg.matrix_rank(pullback, tol=1e-9)),
            "nonzero_eigenvalues": [float(x) for x in nonzero],
            "meaning": "Unique symmetric lift in the span of K with minimum Frobenius norm.",
        },
        "schlaefli_1_10_16_audit": {
            "references_checked": 27,
            "distinct_block_trace_profiles": len(keys),
            "block_trace_ranges": [[float(traces[:, i].min()), float(traces[:, i].max())] for i in range(3)],
            "profiles": profiles,
            "kernel_support_audits": support_audits,
        },
        "scope": {
            "counterfeit_reason": "K contains a nonzero nilpotent vector and is not a Vinberg Cartan embedding.",
            "reference_dependence": "The 1+10+16 block traces vary with the selected cubic-surface line.",
            "required_for_physics": [
                "an explicit semisimple Cartan embedding C^3 -> 3 tensor 27",
                "a physical kinetic metric",
                "a chosen SO(10) reference line or a proof of reference independence",
            ],
            "mass_ratio_claim": False,
            "mixing_angle_claim": False,
        },
        "checks": {
            "compression_identity": residual < 1e-8,
            "lift_rank_three": len(nonzero) == 3,
            "all_27_reference_lines_checked": len(profiles) == 27,
            "reference_dependence_detected": len(keys) > 1,
        },
    }
    assert all(out["checks"].values())
    if write:
        OUT.write_text(json.dumps(out, indent=2) + "\n")
    return out


if __name__ == "__main__":
    print(json.dumps(main(True), indent=2))
