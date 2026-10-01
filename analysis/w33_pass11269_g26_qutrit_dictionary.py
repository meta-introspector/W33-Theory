"""Pass 11269: a G26 Cartan dictionary -- in Codex's own coordinates the 21 mirrors are the 9 SIC and the 12 stabilizer
states of a qutrit, and the reflection Jacobian is a product of SIC and MUB overlaps.

Background (cited, not re-derived): G26 = C2 x G25 with G25 = 3^{1+2}.SL(2,3) (data/w33_extended_clifford_g26_no_go.json,
literature: Magma handbook); its collineation group is the Hessian group of order 216 = the projective one-qutrit
Clifford group ASL(2,3) (2026-09-23_extended_clifford_hesse_null_cone.md); the Hesse configuration's 9 flexes are the
qutrit SIC (Hughston 2007; Bengtsson et al.).  Passes 11255-11259 (parallel session) use Maschke's standard basic
invariants u6, u12, u9^2 in coordinates (x, y, z) and need the reflection Jacobian for vacuum ramification.

Checks here (all in those coordinates, reading (x, y, z) as qutrit amplitudes):
  1. det(d(u6, u12, u9^2)/d(x, y, z)) = 19440 * prod_{9 SIC} <phi|v> * prod_{12 stabilizer} <s|v>^2 (bilinear overlaps
     with the conjugate kets), proved as a polynomial identity in Q[w]/(w^2+w+1);
  2. the 9 order-2 mirror normals form a SIC (|<phi_a|phi_b>|^2 = 1/4) and a Weyl orbit; the 12 order-3 mirror normals
     are the 4 MUBs (overlaps 0, 1/3) = all qutrit stabilizer states;
  3. the reflections (order 3 about stabilizer states, order 2 about SIC vectors) are qutrit Clifford unitaries and
     generate a group of order 1296 whose collineation image has order 216;
  4. EXACT: u6 = (1/6) sum_a (n_a.v)^6, u12 = (1/6) sum_a (n_a.v)^12, u9 = prod_a (n_a.v) over the 9 SIC normals
     n_a -- every basic invariant is a SIC moment (stabilizer power sums do not work: least-squares residual 0.54);
  5. Codex's control point (1, 2, 3) is regular: no SIC or stabilizer overlap vanishes.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11255_11259_g26_common as G  # noqa: E402

OUT = ROOT / "data" / "w33_pass11269_g26_qutrit_dictionary.json"
w = sp.Rational(-1, 2) + sp.sqrt(3) * sp.I / 2
WN = np.exp(2j * np.pi / 3)
SIC = [(1, -w ** k, 0) for k in range(3)] + [(1, 0, -w ** k) for k in range(3)] + [(0, 1, -w ** k) for k in range(3)]
STAB = [(1, 0, 0), (0, 1, 0), (0, 0, 1)] + [(1, w ** a, w ** b) for a in range(3) for b in range(3)]


def num(vs):
    return [np.array([complex(sp.N(c)) for c in v]) for v in vs]


def exact_jacobian_check():
    """Prove the degree-33 factorisation symbolically over Q[w]/(w^2+w+1)."""
    x, y, z = G.VARS
    J, _ = G.invariant_jacobian()
    detJ = sp.expand(J.det())
    W = sp.Symbol("W")
    modulus = sp.Poly(W ** 2 + W + 1, W)
    sic = [(1, -W ** k, 0) for k in range(3)] + [(1, 0, -W ** k) for k in range(3)] + [(0, 1, -W ** k) for k in range(3)]
    stab = [(1, 0, 0), (0, 1, 0), (0, 0, 1)] + [(1, W ** a, W ** b) for a in range(3) for b in range(3)]

    def reduce_w(expr):
        return sp.expand(sp.rem(sp.Poly(sp.expand(expr), W), modulus).as_expr())

    rhs = sp.Integer(19440)
    for n in sic:
        rhs = reduce_w(rhs * (n[0] * x + n[1] * y + n[2] * z))
    for n in stab:
        form = n[0] * x + n[1] * y + n[2] * z
        rhs = reduce_w(rhs * form * form)
    difference = sp.expand(rhs - detJ)
    return difference == 0, len(sp.Poly(rhs, x, y, z).terms())


def reflection_group():
    S = [v / np.linalg.norm(v) for v in num(SIC)]
    T = [v / np.linalg.norm(v) for v in num(STAB)]
    gens = [np.eye(3) - 2 * np.outer(v.conj(), v) for v in S]              # order-2 reflections, mirror <phi|v> = 0
    gens += [np.eye(3) + (WN - 1) * np.outer(v.conj(), v) for v in T]      # order-3 reflections
    X = np.roll(np.eye(3), 1, axis=0)
    Z = np.diag([1, WN, WN ** 2])

    def is_clifford(M):
        for P in (X, Z):
            Q = M @ P @ M.conj().T
            ok = False
            for a in range(3):
                for b in range(3):
                    W_ = np.linalg.matrix_power(X, a) @ np.linalg.matrix_power(Z, b)
                    r = np.trace(W_.conj().T @ Q) / 3
                    if abs(abs(r) - 1) < 1e-9 and np.allclose(Q, r * W_):
                        ok = True
            if not ok:
                return False
        return True

    def key(M):
        return tuple((np.round(M.ravel(), 6) + 0.0).view(float) + 0.0)

    def pkey(M):
        k = np.argmax(np.abs(M.ravel()) > 1e-9)
        Mn = M / (M.ravel()[k] / abs(M.ravel()[k]))
        return tuple((np.round(Mn.ravel(), 6) + 0.0).view(float) + 0.0)

    els = [np.eye(3, dtype=complex)]
    arr = np.array(els)
    frontier = [np.eye(3, dtype=complex)]
    while frontier and len(els) < 5000:
        new = []
        for M in frontier:
            for g in gens:
                N = g @ M
                if np.abs(arr - N[None]).reshape(len(arr), -1).max(axis=1).min() > 1e-6:
                    els.append(N)
                    arr = np.array(els)
                    new.append(N)
        frontier = new
    # collineations: identify up to unimodular scalars
    reps = []
    for M in els:
        k = int(np.argmax(np.abs(M.ravel()) > 1e-9))
        Mn = M / (M.ravel()[k] / abs(M.ravel()[k]))
        if not reps or np.abs(np.array(reps) - Mn[None]).reshape(len(reps), -1).max(axis=1).min() > 1e-6:
            reps.append(Mn)
    proj = reps
    return dict(order=len(els), collineation_order=len(proj),
                all_reflections_clifford=all(is_clifford(g) for g in gens),
                all_elements_clifford=all(is_clifford(M) for M in els),
                sic_overlaps=sorted({round(abs(np.vdot(a, b)) ** 2, 9) for a in S for b in S}),
                stab_overlaps=sorted({round(abs(np.vdot(a, b)) ** 2, 9) for a in T for b in T}))


def invariant_dictionary():
    """exact polynomial identities: u6 = (1/6) sum_a (n_a.v)^6, u12 = (1/6) sum_a (n_a.v)^12, u9 = prod_a (n_a.v), over
    the 9 SIC normals n_a (unnormalised; |n_a|^2 = 2).  Checked symbolically after reduction w^2 = -1 - w."""
    x, y, z = G.VARS
    u6, u12, _ = G.invariants()
    u9 = sp.expand((x ** 3 - y ** 3) * (x ** 3 - z ** 3) * (y ** 3 - z ** 3))
    W = sp.Symbol("W")
    sicW = [(1, -W ** k, 0) for k in range(3)] + [(1, 0, -W ** k) for k in range(3)] + [(0, 1, -W ** k) for k in range(3)]
    forms = [n[0] * x + n[1] * y + n[2] * z for n in sicW]

    def reduce_w(e):
        return sp.expand(sp.rem(sp.Poly(sp.expand(e), W), sp.Poly(W ** 2 + W + 1, W)).as_expr())

    s6 = reduce_w(sum(f ** 6 for f in forms) / 6)
    s12 = reduce_w(sum(f ** 12 for f in forms) / 6)
    pr = sp.Integer(1)
    for f in forms:
        pr = sp.expand(pr * f)
    pr = reduce_w(pr)
    return dict(u6_is_sixth_SIC_power_sum=bool(sp.expand(s6 - u6) == 0),
                u12_is_twelfth_SIC_power_sum=bool(sp.expand(s12 - u12) == 0),
                u9_vs_SIC_product=str(sp.simplify(pr / u9)) if sp.expand(pr) != 0 else "0")


def control_point(p=(1, 2, 3)):
    v = np.array(p, dtype=complex)
    so = [abs(n @ v) for n in num(SIC)]
    to = [abs(n @ v) for n in num(STAB)]
    return dict(point=list(p), min_SIC_overlap=float(min(so)), min_stab_overlap=float(min(to)),
                regular=bool(min(so) > 1e-9 and min(to) > 1e-9))


def run():
    res = dict(pass_id=11269)
    ok, terms = exact_jacobian_check()
    res["jacobian_equals_19440_SIC_x_stab2_exact"] = ok
    res["jacobian_symbolic_identity"] = "Q[w]/(w^2+w+1)[x,y,z]"
    res["jacobian_nonzero_monomials"] = terms
    print("jacobian", ok, flush=True)
    res["reflection_group"] = reflection_group()
    print(res["reflection_group"], flush=True)
    res["invariants"] = invariant_dictionary()
    print(res["invariants"], flush=True)
    res["codex_control_point"] = control_point()
    res["degree_check"] = dict(SIC_mirrors=9, stab_mirrors=12, jacobian_degree=9 * 1 + 12 * 2, sum_degrees_minus_3=5 + 11 + 17)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
