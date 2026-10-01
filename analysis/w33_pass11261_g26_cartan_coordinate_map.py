#!/usr/bin/env python3
"""Pass 11261: exact G26 reflection coordinates on Pass 11260's Cartan plane.

The fixed complementary Pfaffian is reconstructed and checked on a unisolvent
degree-39 grid.  Its three rational factors split over Q(cuberoot(2), zeta_18)
into the 12 qutrit stabilizer mirrors and 9 qutrit SIC mirrors of Pass 11269.
"""
from __future__ import annotations

import importlib.util
import json
from fractions import Fraction
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/w33_pass11261_g26_cartan_coordinate_map.json"
x, y, z = sp.symbols("x y z")
a, t = sp.symbols("a t")

C3 = 2*x**3 - 6*x*y*z + 4*y**3 + z**3
N9 = (
    8*x**9 - 72*x**7*y*z + 48*x**6*y**3 - 96*x**6*z**3
    - 432*x**5*y**2*z**2 - 1584*x**4*y**4*z - 396*x**4*y*z**4
    - 768*x**3*y**6 - 1464*x**3*y**3*z**3 + 6*x**3*z**6
    - 864*x**2*y**5*z**2 - 216*x**2*y**2*z**5 - 288*x*y**7*z
    - 792*x*y**4*z**4 - 18*x*y*z**7 + 64*y**9 + 48*y**6*z**3
    - 96*y**3*z**6 + z**9
)
S9 = (
    8*x**9 + 144*x**7*y*z - 384*x**6*y**3 + 12*x**6*z**3
    + 216*x**5*y**2*z**2 - 720*x**4*y**4*z - 180*x**4*y*z**4
    + 96*x**3*y**6 + 480*x**3*y**3*z**3 - 48*x**3*z**6
    + 432*x**2*y**5*z**2 + 108*x**2*y**2*z**5 + 576*x*y**7*z
    - 360*x*y**4*z**4 + 36*x*y*z**7 + 64*y**9 - 384*y**6*z**3
    + 12*y**3*z**6 + z**9
)
PFAFFIAN = sp.expand(2**17 * C3 * N9 * S9**3)


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


def pfaffian_exact(M):
    """Fraction-free Pfaffian elimination for an integral skew matrix."""
    A = [[int(q) for q in row] for row in M]
    n = len(A); sign = 1; previous = 1
    for k in range(0, n, 2):
        jj = next((j for j in range(k + 1, n) if A[k][j]), None)
        if jj is None:
            return 0
        if jj != k + 1:
            A[k + 1], A[jj] = A[jj], A[k + 1]
            for row in A:
                row[k + 1], row[jj] = row[jj], row[k + 1]
            sign = -sign
        pivot = A[k][k + 1]
        divisor = 1 if k == 0 else previous
        for i in range(k + 2, n):
            for j in range(i + 1, n):
                q = (pivot*A[i][j] - A[k][i]*A[k+1][j] + A[k][j]*A[k+1][i])
                assert q % divisor == 0
                A[i][j] = q // divisor; A[j][i] = -A[i][j]
        previous = pivot
    return sign * A[n-2][n-1]


def slice_matrices():
    p60 = _load(ROOT / "analysis/w33_pass11260_exact_semisimple_cartan.py", "p60_for_61")
    old = _load(ROOT / "analysis/w33_pass11218_e8_cubic_cartan_kernel.py", "old_for_61")
    parent = old.load_parent(); records = parent.build_records(); v = p60.dense()
    A = sp.Matrix(parent.jacobian(v, records))
    kernel = [old.primitive_integer_vector(q) for q in A.nullspace()]
    K = sp.Matrix.hstack(*(sp.Matrix(q) for q in kernel))
    pivots = tuple(map(int, K.T.rref()[1])); keep = [i for i in range(81) if i not in pivots]
    mats = []
    for q in kernel:
        M = parent.jacobian(q, records)
        mats.append([[int(M[i][j]) for j in keep] for i in keep])
    return kernel, pivots, mats


def reduce_field(expr):
    """Canonical remainder in Q[a,t]/(a^3-2, Phi_18(t))."""
    e = sp.expand(expr)
    e = sp.rem(sp.Poly(e, t), sp.Poly(t**6 - t**3 + 1, t)).as_expr()
    e = sp.rem(sp.Poly(sp.expand(e), a), sp.Poly(a**3 - 2, a)).as_expr()
    e = sp.rem(sp.Poly(sp.expand(e), t), sp.Poly(t**6 - t**3 + 1, t)).as_expr()
    return sp.expand(e)


def coordinate_map():
    w = t**6
    l0 = a*x + a**2*y + z
    l1 = a*x + w*a**2*y + w**2*z
    l2 = a*x + w**2*a**2*y + w*z
    X, Y, Z = l0, t**8*l1, -t*l2       # t^8 = -t^{-1}, since Phi_18(t)=0
    u9 = (X**3-Y**3)*(X**3-Z**3)*(Y**3-Z**3)
    stab9 = sp.prod(X + w**i*Y + w**j*Z for i in range(3) for j in range(3))
    assert reduce_field(X*Y*Z - C3) == 0
    assert reduce_field(u9 - 3*(w-w**2)*S9) == 0
    assert reduce_field(X*Y*Z*stab9 + 27*C3*N9) == 0
    return X, Y, Z


def run(full_grid=True):
    kernel, pivots, mats = slice_matrices()
    checked = 0; max_digits = 0
    if full_grid:
        for x0 in range(40):
            for y0 in range(40-x0):
                M = [[x0*mats[0][i][j] + y0*mats[1][i][j] + mats[2][i][j]
                      for j in range(78)] for i in range(78)]
                got = pfaffian_exact(M)
                want = int(PFAFFIAN.subs({x:x0, y:y0, z:1}))
                assert got == want
                checked += 1; max_digits = max(max_digits, len(str(abs(got))))
        assert checked == 820
    X, Y, Z = coordinate_map()
    out = {
        "schema": "w33.pass11261.g26-cartan-coordinate-map.v1",
        "status": "PASS_EXACT_G26_COORDINATE_MAP_AND_QUTRIT_MIRROR_DIVISOR",
        "cartan_basis_source": "data/w33_pass11260_exact_semisimple_cartan.json",
        "kernel_pivot_rows": list(pivots),
        "restricted_pfaffian": {
            "formula": "2^17*C3*N9*S9^3", "degree": 39,
            "factor_degrees_and_multiplicities": [[3,1],[9,1],[9,3]],
            "C3": str(C3), "N9": str(N9), "S9": str(S9),
            "unisolvent_grid_points": checked, "largest_value_digits": max_digits,
        },
        "number_field": {"a": "a^3=2", "t": "Phi_18(t)=t^6-t^3+1=0", "omega": "t^6"},
        "linear_map_standard_from_cartan": {"X": str(X), "Y": str(Y), "Z": str(Z)},
        "exact_map_identities": {
            "XYZ": "C3",
            "standard_u9": "3*(omega-omega^2)*S9",
            "product_12_stabilizer_mirrors": "-27*C3*N9",
        },
        "mirror_divisor": {
            "12_qutrit_stabilizer_mirrors": {"multiplicity": 1, "rational_norm": "C3*N9"},
            "9_qutrit_SIC_mirrors": {"multiplicity": 3, "rational_norm": "S9^3"},
            "connection": (
                "Pass 11269's standard-coordinate Jacobian has SIC multiplicity one and stabilizer "
                "multiplicity two. The E8 restricted Pfaffian has the same 21 hyperplanes with "
                "multiplicities three and one. They are distinct relative divisors on one arrangement."
            ),
        },
        "basic_invariants": {
            "degrees": [6,12,18],
            "pullback_rule": "substitute the displayed X,Y,Z into Pass 11256's u6,u12,u9^2",
            "algebraic_independence": "follows from the invertible displayed linear map",
        },
        "boundary": (
            "The map is defined over Q(cuberoot(2),zeta_18), not over the rational Cartan basis. "
            "Mirror and invariant identities are exact; no Hermitian kinetic metric is selected."
        ),
    }
    OUT.write_text(json.dumps(out, indent=2) + "\n")
    return out


if __name__ == "__main__":
    print(json.dumps(run(True), indent=2))
