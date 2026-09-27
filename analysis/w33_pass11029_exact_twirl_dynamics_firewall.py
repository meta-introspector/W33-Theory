#!/usr/bin/env python3
"""Pass 11029: exact twirl dynamics and the 3D-clock equivariance firewall."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import w33_pass10975_e6_cubic_reynolds_clock_weld as P75
import w33_pass11021_minimal_signed_clock_carrier as P21

OUT = ROOT / "data/w33_pass11029_exact_twirl_dynamics_firewall.json"


def vec_rank(vs):
    return sp.Matrix.hstack(*[sp.Matrix(v) for v in vs]).rank()


def payload():
    e6_to_h, dirs, rows, perms, section = P21.split_section()
    noncentral = tuple(sorted(i for i, h in e6_to_h.items() if h[:2] != (0, 0)))
    pos = {i: j for j, i in enumerate(noncentral)}
    fibres = []
    for d in dirs:
        fibre = sorted(
            i for i, h in e6_to_h.items()
            if h[:2] != (0, 0) and P75.norm_dir(h[:2]) == d
        )
        v = [0] * 27
        for i in fibre:
            v[i] = 1
        fibres.append(tuple(v))

    aug = [
        tuple(a - b for a, b in zip(fibres[i], fibres[3]))
        for i in range(3)
    ]
    assert vec_rank(aug) == 3

    p_to_m = {r[3]: tuple(map(int, r[0])) for r in rows}
    z = next(g for g in section if p_to_m[g[0]] == (2, 0, 0, 2))
    z_aug = [P21.act_vec(z, v) for v in aug]
    assert vec_rank(aug + z_aug) == 6

    # Exact signed-permutation matrices on the 24D noncentral carrier.
    mats = []
    traces = []
    for p, s in section:
        U = sp.zeros(24)
        for i in noncentral:
            U[pos[p[i]], pos[i]] = -1 if s[i] else 1
        assert U.T * U == sp.eye(24)
        mats.append(U)
        traces.append(int(sp.trace(U)))
    assert len(mats) == 48
    R = sum(mats, sp.zeros(24)) / 48
    assert R * R == R
    assert R.T == R
    vector_fixed_dim = int(R.rank())
    assert vector_fixed_dim == 2

    # Conjugation twirl on End(V): rank equals character inner product.
    commutant_dim = sum(t * t for t in traces) // 48
    assert commutant_dim == 19

    # The canonical finite-group twirl is CPTP/unital because its Kraus
    # operators are U_g/sqrt(48); L=T-I exponentiates to
    # exp(tL)=e^-t I +(1-e^-t)T, a convex combination of channels.
    closure_rank = P21.orbit_span_rank(section, aug)[0]
    assert closure_rank == 24

    checks = {
        "signed_group_order_48": len(mats) == 48,
        "all_representation_matrices_orthogonal": True,
        "vector_reynolds_idempotent": R * R == R,
        "vector_fixed_dimension_2": vector_fixed_dim == 2,
        "operator_twirl_commutant_dimension_19": commutant_dim == 19,
        "clock_augmentation_dimension_3": vec_rank(aug) == 3,
        "central_minusI_breaks_clock_subspace": vec_rank(aug + z_aug) == 6,
        "clock_signed_orbit_closes_to_24": closure_rank == 24,
    }
    assert all(checks.values())
    return {
        "schema": "w33.pass11029.exact-twirl-dynamics-firewall.v1",
        "status": "PASS",
        "headline": (
            "The canonical exact signed-GL2(3) twirl defines a CPTP unital "
            "conditional expectation on End(V24), and L=T-I gives a CPTP "
            "Markov semigroup. Its fixed operator algebra has dimension 19; "
            "the vector Reynolds projector has rank 2. No GL2(3)-equivariant "
            "idempotent can have the naive 3D clock augmentation as image, "
            "because that 3D subspace is not invariant and its signed orbit "
            "already spans all 24 dimensions."
        ),
        "carrier": {
            "dimension": 24,
            "group": "exact signed GL2(3)",
            "group_order": 48,
            "vector_fixed_dimension": vector_fixed_dim,
            "operator_commutant_dimension": commutant_dim,
        },
        "canonical_channel": {
            "twirl": "T(X)=(1/48) sum_g U_g X U_g^T",
            "kraus_family": "K_g=U_g/sqrt(48)",
            "properties": ["completely_positive", "trace_preserving", "unital", "idempotent"],
            "lindbladian": "L=T-I",
            "semigroup": "exp(tL)=exp(-t) I + (1-exp(-t)) T",
            "fixed_algebra_dimension": commutant_dim,
        },
        "clock_firewall": {
            "coarse_augmentation_dimension": 3,
            "augmentation_plus_central_image_rank": 6,
            "full_signed_orbit_span_rank": closure_rank,
            "theorem": (
                "The image of an equivariant idempotent is invariant. Since the "
                "naive three-dimensional clock augmentation is not invariant under "
                "the exact signed group, no exact signed-GL2(3)-equivariant "
                "conditional expectation or linear projector can select precisely "
                "that subspace."
            ),
            "consequence": (
                "A physical mechanism selecting the coarse 3D clock must either "
                "break/reduce the exact signed symmetry, use an explicit quotient "
                "or readout map, or introduce additional dynamics beyond the "
                "canonical group twirl."
            ),
        },
        "boundary": (
            "The CPTP statement concerns the finite 24-mode representation and "
            "its operator algebra. It is not a microscopic open-system derivation, "
            "a measured relaxation rate, or a proof that nature implements this "
            "twirl. The no-go applies to exact full signed-GL2(3) equivariance."
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
        "vector_fixed": p["carrier"]["vector_fixed_dimension"],
        "commutant": p["carrier"]["operator_commutant_dimension"],
        "clock_orbit": p["clock_firewall"]["full_signed_orbit_span_rank"],
    }, sort_keys=True))


if __name__ == "__main__":
    raise SystemExit(main())
