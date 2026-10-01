#!/usr/bin/env python3
"""Pass 11256: an explicit controlled potential on the standard G26 quotient.

This builds a reproducible orbit-distance potential from the three basic
invariants of degrees 6, 12, and 18.  The coordinates are the standard G26
reflection coordinates; Pass 11255 forbids identifying them with the old
Pass 11218 centralizer without a new semisimple embedding.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/w33_pass11256_g26_invariant_potential.json"


def load_common():
    path = ROOT / "analysis/w33_pass11255_11259_g26_common.py"
    spec = importlib.util.spec_from_file_location("p11255_common", path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def main(write=True):
    common = load_common()
    fs, vals, Jn, H = common.target_potential_data((1, 2, 3))
    J, determinant = common.invariant_jacobian()
    det_at_target = int(determinant.subs(dict(zip(common.VARS, (1, 2, 3)))))
    charpoly = sp.factor(H.charpoly().as_expr())
    eig = np.linalg.eigvalsh(np.array(H, dtype=float))
    factors = sp.factor_list(determinant)[1]
    factor_degrees = sorted((sp.Poly(f, common.VARS).total_degree(), int(e)) for f, e in factors)
    weighted_degree = sum(d * e for d, e in factor_degrees)
    out = {
        "schema": "w33.pass11256.g26-invariant-potential.v1",
        "status": "PASS_CONTROLLED_G26_ORBIT_POTENTIAL_CONSTRUCTED",
        "basic_invariants": {
            "degrees": [6, 12, 18],
            "formulas": [str(f) for f in fs],
            "target_point": [1, 2, 3],
            "target_values": [int(v) for v in vals],
        },
        "reflection_jacobian": {
            "degree": int(sp.Poly(determinant, common.VARS).total_degree()),
            "factorization": str(determinant),
            "factor_degrees_and_multiplicities": factor_degrees,
            "weighted_degree_check": weighted_degree,
            "complex_reflection_hyperplanes_counted_with_distinct_linear_factors": 21,
            "value_at_target": det_at_target,
        },
        "potential": {
            "formula": "V=sum_i (u_i/u_i(p)-1)^2",
            "scope": "A dimensionless invariant orbit-distance potential, not a fitted physical scalar potential.",
            "normalized_jacobian": common.rational_matrix(Jn),
            "hessian": common.rational_matrix(H),
            "hessian_characteristic_polynomial": str(charpoly),
            "hessian_trace": str(sp.factor(sp.trace(H))),
            "hessian_determinant": str(sp.factor(H.det())),
            "hessian_eigenvalues": [float(x) for x in eig],
        },
        "checks": {
            "degrees_are_6_12_18": [sp.Poly(f, common.VARS).total_degree() for f in fs] == [6, 12, 18],
            "jacobian_degree_is_33": weighted_degree == 33,
            "target_is_regular": det_at_target != 0,
            "hessian_is_positive_definite": bool(np.all(eig > 0)),
            "hessian_identity": H == sp.simplify(2 * Jn.T * Jn),
        },
        "firewall": {
            "pass11218_coordinates_used": False,
            "physical_mass_claim": False,
            "required_next_map": "An explicit semisimple Cartan embedding into 3 tensor 27 with a kinetic metric.",
        },
    }
    assert all(out["checks"].values())
    if write:
        OUT.write_text(json.dumps(out, indent=2) + "\n")
    return out


if __name__ == "__main__":
    print(json.dumps(main(True), indent=2))
