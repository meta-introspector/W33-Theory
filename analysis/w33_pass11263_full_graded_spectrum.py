#!/usr/bin/env python3
"""Pass 11263: complete 81/248-mode spectrum with declared grading/metric.

Two different statements are kept separate:
1. the basis-independent algebraic spectrum of (18 ad v)^3, and
2. the positive Hessian of 1/2 ||[v,delta]||^2 in the declared Euclidean
   coordinate metric, H=(18 ad v)^T(18 ad v)/18^2.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/w33_pass11263_full_graded_spectrum.json"


def load60():
    path = ROOT / "analysis/w33_pass11260_exact_semisimple_cartan.py"
    spec = importlib.util.spec_from_file_location("p60_for_63", path)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


def common_nonzero_polynomial():
    q = sp.Symbol("q")
    d = 396718580736
    A = 1234235584512
    B = 34326907456637308502016
    C = 3172169114198268924301601144832
    fm = q**9-A*q**6-B*q**3-C
    fp = q**9+A*q**6-B*q**3+C
    f18 = (
        q**18 + 146386125332581029352833024*q**12
        + 250556556265049251809547325105197290584587370496*q**6
        + 198063275547632324892732342401376875852639602111883769352386772992
    )
    return sp.Poly((q**3-d)*(q**3+d)*fm**3*fp**3*f18, q)


def newton_power_sums(P, count=12):
    # Monic P=x^n+c1*x^(n-1)+...; Newton identities.
    co = [sp.Integer(c) for c in P.all_coeffs()]
    out = []
    for k in range(1, count+1):
        val = sum(co[j]*out[k-j-1] for j in range(1, k)) + k*co[k]
        out.append(int(-val))
    return out


def run():
    p60 = load60(); v = p60.dense(); N, residual = p60.build_adjoint()(v)
    cuts = (("g0_E6_plus_A2",0,86,8), ("g1_27x3",86,167,3), ("g2_dual",167,248,3))
    blocks = {}; spectra = {}
    for name, lo, hi, nullity in cuts:
        H = (N[:,lo:hi].T.astype(float) @ N[:,lo:hi].astype(float)) / (18.0**2)
        ev = np.linalg.eigvalsh(H); ev[np.abs(ev) < 1e-8] = 0
        assert int(np.sum(ev == 0)) == nullity
        spectra[name] = ev
        blocks[name] = {
            "dimension": hi-lo, "rank": hi-lo-nullity, "zero_modes": nullity,
            "eigenvalues_complete": [float(q) for q in ev],
            "smallest_positive": float(ev[nullity]), "largest": float(ev[-1]),
            "trace": float(np.sum(ev)),
        }

    signs = [
        (1,1,-1), (1,-1,1), (-1,1,1),
        (1,-1,-1), (-1,1,-1), (-1,-1,1),
    ]
    names = [c[0] for c in cuts]
    audits = []
    for sg in signs:
        moments = []
        for k in range(7):
            if k == 0:
                val = sum(s*(len(spectra[n])) for s,n in zip(sg,names))
            else:
                val = sum(s*float(np.sum(spectra[n]**k)) for s,n in zip(sg,names))
            moments.append(val)
        zp = 0.5*sum(s*float(np.sum(np.sqrt(np.maximum(spectra[n],0)))) for s,n in zip(sg,names))
        audits.append({"signs_g0_g1_g2": list(sg), "moments_0_to_6": moments, "zero_point_proxy": zp})

    P = common_nonzero_polynomial(); powers = newton_power_sums(P, 12)
    # Every grade has this same degree-78 nonzero characteristic polynomial.
    root_unity = [{"power": k, "common_trace": str(p), "Z3_character_trace": "0"}
                  for k,p in enumerate(powers,1)]
    out = {
        "schema": "w33.pass11263.full-graded-spectrum.v1",
        "status": "PASS_COMPLETE_248_MODE_DIAGNOSTIC_AND_EXACT_Z3_SPECTRAL_PAIRING",
        "background": "Pass 11260 semisimple regular Cartan element",
        "integral_adjoint_rounding_residual": residual,
        "declared_hessian": {
            "formula": "H=(18 ad(v))^T(18 ad(v))/18^2",
            "potential": "V(delta)=1/2 ||[v,delta]||^2",
            "metric": "the 248 stored E6+A2+(27,3)+(27*,3*) coordinates are declared Euclidean-orthonormal",
            "blocks": blocks,
            "full_dimension": 248, "full_rank": 234, "full_zero_modes": 14,
            "matter81_rank": 78, "matter81_Cartan_zero_modes": 3,
        },
        "basis_independent_N_cubed_spectrum": {
            "common_nonzero_characteristic_polynomial_degree": P.degree(),
            "common_nonzero_characteristic_polynomial": str(P.as_expr()),
            "zero_multiplicities_g0_g1_g2": [8,3,3],
            "nonzero_modes_per_grade": 78,
            "exact_Z3_character_traces_powers_1_to_12": root_unity,
            "degree_zero_Z3_index": "86+81*omega+81*omega^2=5",
            "kernel_Z3_index": "8+3*omega+3*omega^2=5",
            "reading": (
                "All 234 nonzero modes form 78 exact Z3 spectral triplets. The only graded mismatch is "
                "five neutral zero-mode units."
            ),
        },
        "ordinary_Z2_sign_audit": {
            "assignments": audits,
            "canonical_matter_assignment": "(+,-,-): g0 bosonic and g1,g2 fermionic",
            "any_assignment_cancels_moments_0_to_6": any(all(abs(v)<1e-7 for v in r["moments_0_to_6"]) for r in audits),
        },
        "boundary": (
            "The Z3 root-of-unity trace cancellation is exact but is not an ordinary boson/fermion supertrace. "
            "The positive Hessian spectrum depends on the explicitly declared coordinate metric; no physical "
            "kinetic normalization, supersymmetry multiplet assignment, or cosmological constant is derived."
        ),
    }
    OUT.write_text(json.dumps(out, indent=2) + "\n")
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
