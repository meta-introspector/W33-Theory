#!/usr/bin/env python3
"""Pass 11259: discriminate G26 ramification, CP, and supertrace claims.

For a potential centred on an orbit, its Hessian at that orbit is 2 J^T J.
We use this exact identity to compare regular and reflection-discriminant
strata.  We then exhaust all boson/fermion signings of the three controlled
normalized modes.  This is a finite diagnostic; it is not a cosmological-
constant calculation because no physical spectrum or statistics assignment
has been derived.
"""
from __future__ import annotations

import importlib.util
import itertools
import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/w33_pass11259_g26_discriminant_supertrace.json"


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def stratum_record(common, cp, point):
    fs = common.invariants()
    J, detJ = common.invariant_jacobian()
    sub = dict(zip(common.VARS, point))
    values = [int(f.subs(sub)) for f in fs]
    Jp = J.subs(sub)
    H = 2 * Jp.T * Jp
    eig = np.linalg.eigvalsh(np.array(H, dtype=float))
    _, W, _ = cp.induced(point, True, True)
    return {
        "point": list(point),
        "invariants": values,
        "reflection_jacobian": int(detJ.subs(sub)),
        "jacobian_rank": int(Jp.rank()),
        "orbit_distance_hessian_rank": int(H.rank()),
        "orbit_distance_hessian_eigenvalues": [float(x) for x in eig],
        "orbit_distance_hessian_trace": str(sp.trace(H)),
        "conditional_cp_jarlskog": cp.jarlskog(W),
    }


def main(write=True):
    common = load(ROOT / "analysis/w33_pass11255_11259_g26_common.py", "p11255_common_11259")
    cp = load(ROOT / "analysis/w33_pass11258_g26_cubic_cp_interface.py", "p11258_for_11259")
    points = [(1, 2, 3), (1, 1, 0), (1, 0, 0), (1, 1, 1), (1, -1, 0), (1, 1, 2)]
    strata = [stratum_record(common, cp, p) for p in points]

    _, _, _, Hnormalized = common.target_potential_data((1, 2, 3))
    modes = np.linalg.eigvalsh(np.array(Hnormalized, dtype=float))
    signings = [s for s in itertools.product((-1, 1), repeat=3) if len(set(s)) > 1]
    moments = []
    for power in range(1, 7):
        candidates = [(abs(sum(s*x**power for s, x in zip(sign, modes))), sign) for sign in signings]
        value, sign = min(candidates)
        moments.append({"power": power, "minimum_absolute_supertrace": float(value), "signs": list(sign)})
    vacuum = min(
        (abs(0.5*sum(s*np.sqrt(x) for s, x in zip(sign, modes))), sign)
        for sign in signings
    )
    reflection_cp_counterexample = next(r for r in strata if r["point"] == [1, 1, 2])

    out = {
        "schema": "w33.pass11259.g26-discriminant-supertrace.v1",
        "status": "PASS_DISCRIMINANT_ZERO_MODES_PROVED_SUPERTRACE_CANCELLATION_ABSENT",
        "strata": strata,
        "theorems": {
            "ramification_zero_modes": "For V_p=sum_i(u_i-u_i(p))^2, Hess_p(V_p)=2 J(p)^T J(p); hence Hessian nullity equals quotient-map ramification nullity.",
            "cp_is_not_the_discriminant": "The reflection point (1,1,2) has det(J)=0 but a nonzero conditional CP quartet.",
            "supertrace_result": "No nontrivial +/- grading of the three controlled normalized modes cancels any tested moment 1 through 6.",
        },
        "normalized_regular_modes": [float(x) for x in modes],
        "exhaustive_nontrivial_signings": len(signings),
        "minimum_supertrace_moments": moments,
        "vacuum_proxy": {
            "formula": "abs((1/2) sum_i sign_i sqrt(lambda_i))",
            "minimum": float(vacuum[0]),
            "signs": list(vacuum[1]),
        },
        "scope": {
            "physical_statistics_assigned": False,
            "physical_complete_spectrum": False,
            "cosmological_constant_claim": False,
            "meaning": "The controlled three-mode G26 Hessian supplies no internal Hecke-like supertrace cancellation.",
        },
        "checks": {
            "generic_target_regular": strata[0]["jacobian_rank"] == 3,
            "listed_reflection_strata_singular": all(r["reflection_jacobian"] == 0 for r in strata[1:]),
            "singular_strata_have_zero_modes": all(r["orbit_distance_hessian_rank"] < 3 for r in strata[1:]),
            "reflection_cp_counterexample_nonzero": abs(reflection_cp_counterexample["conditional_cp_jarlskog"]) > 1e-6,
            "all_six_supertraces_nonzero": all(m["minimum_absolute_supertrace"] > 1e-9 for m in moments),
            "vacuum_proxy_nonzero": bool(vacuum[0] > 1e-9),
        },
    }
    assert all(out["checks"].values())
    if write:
        OUT.write_text(json.dumps(out, indent=2) + "\n")
    return out


if __name__ == "__main__":
    print(json.dumps(main(True), indent=2))
