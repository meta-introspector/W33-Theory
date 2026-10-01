#!/usr/bin/env python3
"""Pass 11255: factor the restricted Pfaffian and audit the Cartan claim.

The Pass 11218 kernel is an exact three-dimensional abelian centralizer.  The
new calculation proves two further facts:

1. its degree-39 principal Pfaffian factors exactly;
2. the first displayed kernel vector is nonzero ad-nilpotent of index five.

The second fact is decisive: a Vinberg Cartan subspace consists of commuting
semisimple elements.  Therefore the old kernel is not a Cartan subspace, and
its Pfaffian is not the G26 reflection Jacobian in disguised coordinates.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/w33_pass11255_restricted_pfaffian_cartan_audit.json"


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def pfaffian_exact(M):
    """Exact skew Gaussian elimination over Q."""
    a = [[Fraction(int(x)) for x in row] for row in M]
    n = len(a)
    pf = Fraction(1)
    for k in range(0, n, 2):
        jj = next((j for j in range(k + 1, n) if a[k][j]), None)
        if jj is None:
            return 0
        if jj != k + 1:
            a[k + 1], a[jj] = a[jj], a[k + 1]
            for row in a:
                row[k + 1], row[jj] = row[jj], row[k + 1]
            pf = -pf
        q = a[k][k + 1]
        pf *= q
        for i in range(k + 2, n):
            for j in range(i + 1, n):
                value = a[i][j] - (a[k][i] * a[k + 1][j] - a[k][j] * a[k + 1][i]) / q
                a[i][j], a[j][i] = value, -value
    assert pf.denominator == 1
    return int(pf)


def build_slice():
    old = load(ROOT / "analysis/w33_pass11218_e8_cubic_cartan_kernel.py", "p11218_old")
    parent = old.load_parent()
    records = parent.build_records()
    A = sp.Matrix(parent.jacobian(old.dense_from_sparse(old.BACKGROUND), records))
    kernel = [old.primitive_integer_vector(v) for v in A.nullspace()]
    keep = [i for i in range(81) if i not in (1, 2, 16)]
    mats = []
    for v in kernel:
        M = parent.jacobian(v, records)
        mats.append([[int(M[i][j]) for j in keep] for i in keep])
    return old, parent, records, kernel, mats


def candidate_pfaffian(x, y, z, pf_abs):
    h6 = y**6 + 18*y**5*z - 4320*y**3*z**3 - 77760*y**2*z**4 - 653184*y*z**5 - 2239488*z**6
    q2 = y**2 + 24*y*z + 72*z**2
    linear = 96*x + 5*y - 120*z
    return -(pf_abs // 5) * h6**3 * q2**10 * linear


def exact_pfaffian_identity(mats, pf_abs):
    """Check an exact unisolvent triangular grid for degree <=39."""
    checked = 0
    max_digits = 0
    for x in range(40):
        for y in range(40 - x):
            M = [[x*mats[0][i][j] + y*mats[1][i][j] + mats[2][i][j]
                  for j in range(78)] for i in range(78)]
            got = pfaffian_exact(M)
            want = candidate_pfaffian(x, y, 1, pf_abs)
            assert got == want
            checked += 1
            max_digits = max(max_digits, len(str(abs(got))))
    assert checked == 820
    return checked, max_digits


def ad_nilpotence_certificate(kernel):
    toe = load(ROOT / "tools/toe_e8_z3graded_bracket_jacobi.py", "toe_e8_pass11255")
    B = np.load(ROOT / "artifacts/e6_27rep_basis_export/E6_basis_78.npy").astype(np.complex128)
    projector = toe.E6Projector(B)
    bracket = toe.E8Z3Bracket(
        e6_projector=projector,
        cubic_triads=toe._load_signed_cubic_triads(),
        scale_g1g1=1.0,
        scale_g2g2=-1.0 / 6.0,
        scale_e6=1.0,
        scale_sl3=1.0 / 6.0,
    )
    sl3 = []
    for i in range(3):
        for j in range(3):
            if i != j:
                m = np.zeros((3, 3), complex); m[i, j] = 1; sl3.append(m)
    sl3 += [np.diag([1, -1, 0]).astype(complex), np.diag([0, 1, -1]).astype(complex)]
    sl3 = np.array(sl3)
    sf = sl3.reshape(8, -1)
    sleft = np.linalg.inv(sf @ sf.conj().T) @ sf.conj()
    zero = toe.E8Z3.zero()

    basis = []
    for i in range(78):
        basis.append(toe.E8Z3(B[i], zero.sl3, zero.g1, zero.g2))
    for i in range(8):
        basis.append(toe.E8Z3(zero.e6, sl3[i], zero.g1, zero.g2))
    for grade in (1, 2):
        for i in range(81):
            m = np.zeros((27, 3), complex); m.flat[i] = 1
            basis.append(toe.E8Z3(zero.e6, zero.sl3, m if grade == 1 else zero.g1,
                                   m if grade == 2 else zero.g2))

    def coords(e):
        rhs = np.einsum("aij,ji->a", B, e.e6)
        ce6 = projector.gram_inv @ rhs
        csl3 = sleft @ e.sl3.reshape(-1)
        return np.r_[ce6, csl3, e.g1.ravel(), e.g2.ravel()]

    X = np.array(kernel[0], dtype=complex).reshape(27, 3)
    v = toe.E8Z3(zero.e6, zero.sl3, X, zero.g2)
    ad = np.column_stack([coords(bracket.bracket(v, b)) for b in basis])
    scaled = 18 * ad
    N = np.rint(scaled.real).astype(np.int64)
    assert np.max(np.abs(scaled - N)) < 1e-9
    powers = []
    Q = N.copy()
    for k in range(1, 6):
        powers.append({"power": k, "nonzero_entries": int(np.count_nonzero(Q)),
                       "max_abs_entry": int(np.max(np.abs(Q)))})
        if k < 5:
            Q = Q @ N
    assert powers[3]["nonzero_entries"] == 1
    assert powers[4]["nonzero_entries"] == 0
    return powers, int(np.count_nonzero(N)), int(np.linalg.matrix_rank(N.astype(float)))


def main(write=True, full_grid=True):
    old, parent, records, kernel, mats = build_slice()
    prior = old.main(False)
    pf_abs = int(prior["rank_witness"]["absolute_pfaffian"])
    assert pfaffian_exact(mats[1]) == -pf_abs
    checked, digits = exact_pfaffian_identity(mats, pf_abs) if full_grid else (0, 0)
    powers, ad_nnz, ad_rank = ad_nilpotence_certificate(kernel)
    out = {
        "schema": "w33.pass11255.restricted-pfaffian-cartan-audit.v1",
        "status": "PASS_PFAFFIAN_FACTORED_AND_PASS11218_CARTAN_LABEL_RETRACTED",
        "restricted_pfaffian": {
            "coordinates": "x,y,z multiply the three stored Pass11218 kernel vectors",
            "degree": 39,
            "formula": "-(5112848641047920640/5)*H6(y,z)^3*Q2(y,z)^10*(96x+5y-120z)",
            "H6": "y^6+18y^5z-4320y^3z^3-77760y^2z^4-653184yz^5-2239488z^6",
            "Q2": "y^2+24yz+72z^2",
            "factor_degrees_with_multiplicity": [[6, 3], [2, 10], [1, 1]],
            "exact_unisolvent_grid_points": checked,
            "largest_evaluated_integer_digits": digits,
        },
        "semisimplicity_audit": {
            "tested_vector": "first stored kernel basis vector",
            "scaled_adjoint": "N=18 ad(v), an exact integer 248x248 matrix",
            "scaled_adjoint_nonzero_entries": ad_nnz,
            "scaled_adjoint_rank": ad_rank,
            "powers": powers,
            "nilpotency_index": 5,
            "decisive_fact": "N^4 is nonzero and N^5 is zero, so v is nonzero nilpotent and not semisimple.",
        },
        "correction": {
            "retained": "The exact Jacobian rank 78, three-dimensional abelian centralizer, determinant witness, and generic-rank conclusion remain valid.",
            "withdrawn": "The displayed centralizer is a Vinberg Cartan subspace and therefore carries the little-Weyl G26 action.",
            "reason": "A Cartan subspace consists of commuting semisimple elements; this centralizer contains the certified nonzero nilpotent vector.",
            "consequence": "Its degree-39 Pfaffian cannot be used as the G26 reflection Jacobian or as physical vacuum coordinates.",
        },
        "checks": {
            "degree39_factorization_exact": checked == 820 if full_grid else True,
            "ad_fourth_power_nonzero": powers[3]["nonzero_entries"] > 0,
            "ad_fifth_power_zero": powers[4]["nonzero_entries"] == 0,
            "old_rank78_witness_retained": prior["explicit_regular_background"]["exact_rank"] == 78,
        },
    }
    if write:
        OUT.write_text(json.dumps(out, indent=2) + "\n")
    return out


if __name__ == "__main__":
    print(json.dumps(main(True, True), indent=2))
