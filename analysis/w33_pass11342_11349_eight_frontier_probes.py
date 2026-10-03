"""Eight scoped TOE probes, with exact identities and independently replayed numerics.

This packet does not infer observed parameters or a complete UV theory.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
import sympy as sp
from numpy.polynomial.legendre import Legendre, leggauss
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parents[1]


def flavor_barriers():
    previous = json.loads((ROOT / "data/w33_pass11335_larger_complex_flavor_search.json").read_text())
    rows = []
    for row in previous["rows"]:
        if row["S_portal_barrier"] > 0 or row["K_portal_barrier"] > 0:
            continue
        singular = np.array(row["A_singular_values"])
        gamma = float(np.sum(singular**2))
        box = float(np.sum(singular**4))
        # The physical r_S=1/5 background and positive five-coupling form of 11340.
        linear = 344 / 3 - 4 * gamma
        constant = 30 - 4 * box / 5
        minimum = constant - linear**2 / (4 * 124)
        assert row["residual"] < 1e-8 and minimum > 2
        rows.append({"n": row["n"], "start": row["start"],
                     "beta_residual": row["residual"], "barrier_minimum": minimum})
    assert len(rows) == 2 and all(r["n"] == 4 for r in rows)
    return {"status": "PASS", "scope": "Both sampled old-bound escape witnesses, not every complex texture", "rows": rows,
            "certificate": "For these witnesses the same physical S background gives dz/dtau >=124z^2-(344/3-4gamma)z+30-(4/5)box; its minimum is positive."}


def nongaussian_mediator():
    # X=x.sigma, Y=y.sigma, A0=a n.sigma, B0=a m.sigma.
    # TrX²=2|x|² and ||[X,Y]||_F²=8|x cross y|².
    a = sp.Rational(1, 10)
    t = sp.nroots(4 * sp.Symbol("t")**3 + sp.Symbol("t") - a)
    t = float(next(sp.re(z) for z in t if abs(float(sp.im(z))) < 1e-12 and float(sp.re(z)) > 0))
    min_energy = 2 * t*t + 4 * t**4 - 4 * float(a) * t
    # Coercivity: 3r^4-8|x cross y|² >= r^4, r²=|x|²+|y|².
    rng = np.random.default_rng(11343)
    margins = []
    for _ in range(100):
        x, y = rng.normal(size=(2, 3))
        r2 = x@x + y@y
        margins.append(3*r2*r2 - 8*np.linalg.norm(np.cross(x,y))**2 - r2*r2)
    assert min(margins) > -1e-12
    return {"status": "PASS", "scope": "Positive-mass two-adjoint renormalizable mediator sector at fixed two-family Gram spectra; no three-family CP or UV asymptotic freedom",
            "potential": "V=|x|²+|y|²+3(|x|²+|y|²)²-8|x cross y|²-2a(x.n+y.m), a=1/10",
            "global_argument": "For fixed |x|,|y|, source reward is largest when x||n,y||m and negative commutator term is smallest when x perpendicular y. Both are attained together by n perpendicular m. Equal radii then minimize by symmetry; t+4t³=a.",
            "selected_Gram_Bloch_angle_degrees": 90, "two_family_mixing_angle_degrees": 45,
            "positive_quartic_margin_min_sample": float(min(margins)), "mediator_radius": t,
            "global_minimum_energy": min_energy,
            "boundary": "The mediator quartic is bounded below exactly. Fixed Gram spectra are supplied; observed hierarchy, CKM phase, full scalar sector, and RG control remain open."}


def ward_normal_freedom():
    d = json.loads((ROOT / "data/w33_pass11337_quotient_hard_ward_columns.json").read_text())
    T = np.array(d["Goldstone_tangents"])
    K = np.array(d["hard_Ward_columns"])
    W = np.array(d["symmetric_minimal_completion"])
    assert T.shape == (22, 8) and np.linalg.matrix_rank(T) == 8
    U, _, _ = np.linalg.svd(T, full_matrices=True)
    N = U[:, 8:]
    perturb = N @ np.diag(np.linspace(-2, 3, 14)) @ N.T
    assert np.max(abs((W+perturb)@T-K)) < 1e-8
    assert np.linalg.norm(N.T@perturb@N) > 1
    return {"status": "PASS", "scope": "Linear-algebraic information limit of the actual eight Ward columns, not a computed 1PI normal block",
            "quotient_real_dimension": 22, "Ward_rank": 8, "normal_dimension": 14,
            "unfixed_symmetric_normal_entries": 105,
            "perturbed_Ward_residual": float(np.max(abs((W+perturb)@T-K))),
            "normal_change_norm": float(np.linalg.norm(N.T@perturb@N))}


def transverse_rank():
    import sys
    sys.path.insert(0, str(ROOT / "analysis"))
    from w33_pass11289_cycle_gram_gluing_flux import cycle_basis
    B, C = cycle_basis()
    L = B @ B.T
    evals = np.linalg.eigvalsh(L)
    assert B.shape == (80, 160) and C.shape[1] == 81
    assert abs(evals[0]) < 1e-10 and evals[1] > 0
    # At aligned identity spatial coframes, three relative boost components
    # have linearized shift equations (B B.T) tensor I_3.
    return {"status": "PASS", "scope": "Aligned-frame transverse boost linearization on the actual graph; not nonlinear Dirac closure",
            "shift_boost_Jacobian": "(B B.T) tensor I3",
            "rank": 3 * int(np.linalg.matrix_rank(B)), "common_boost_nullity": 3,
            "smallest_positive_graph_eigenvalue": float(evals[1]),
            "native_sites": B.shape[0], "native_edges": B.shape[1], "cycle_rank": C.shape[1]}


def wall_data(h):
    b = sp.symbols("b", positive=True)
    h = sp.Rational(h)
    bw = sp.Rational(5, 4)
    q = (1 + sp.sqrt(1 + 4*h*bw)) / 2
    c = q / 5
    rho = q*q/2 - h - 2*c
    C = rho/3 + q*q/2 + h
    F = C/b - rho*b*b/3 - q*q/(2*b*b) - h
    return b, bw, h, q, c, rho, C, sp.factor(F)


def exact_wall_zero(h):
    b, bw, h, q, c, rho, C, F = wall_data(h)
    chi = 1/b - 1/bw
    v = -F/(q*c*b*b)
    p = b**4*c*c*sp.diff(v, b)/2 + q*c*chi
    rchi = sp.factor(-sp.diff(F*sp.diff(chi, b), b) + q*c*sp.diff(v, b))
    rv = sp.factor(-sp.diff(p, b) + h*b*b*c*c*v/F)
    assert sp.simplify(rchi) == sp.simplify(rv) == 0
    assert sp.simplify(v.subs(b, 1)) == 0
    assert sp.simplify(chi.subs(b, bw)) == 0
    assert sp.simplify(sp.diff(v, b).subs(b, bw)) == 0
    assert sp.simplify(q*q-q-h*bw) == 0
    return {"h": str(h), "q": str(q), "rho": str(sp.simplify(rho)),
            "F": str(F), "chi": str(chi), "V_phi": str(sp.factor(v)),
            "chi_equation_residual": str(rchi), "metric_connection_equation_residual": str(rv)}


def odd_wall_ritz(h, n, quadrature=320):
    b0, _, _, q0, c0, _, _, F0 = wall_data(h)
    q, c = float(q0), float(c0)
    z, w = leggauss(quadrature)
    b = 1 + (z+1)/8
    w = w/8
    F = sp.lambdify(b0, F0, "numpy")(b)
    U, DU, V, DV = [], [], [], []
    for k in range(n):
        p = Legendre.basis(k)
        value, deriv = p(z), 8*p.deriv()(z)
        U.append((1.25-b)*value)
        DU.append(-value+(1.25-b)*deriv)
        V.append((b-1)*value)
        DV.append(value+(b-1)*deriv)
    U, DU, V, DV = (np.array(x).T for x in (U, DU, V, DV))
    A = DU.T@((w*F)[:, None]*DU)
    B = q*c*U.T@(w[:, None]*DV)
    D = DV.T@((w*b**4*c*c/2)[:, None]*DV)
    D += V.T@((w*float(h)*b*b*c*c/F)[:, None]*V)
    M = U.T@(w[:, None]*U)
    return eigh(A-B@np.linalg.solve(D, B.T), M, eigvals_only=True)[:4]


def winding_vectors():
    # Independent coordinate-metric control of the KK/Stueckelberg argument.
    b0, f0, c0, v0 = sp.symbols("b f c v", positive=True)
    metric = sp.Matrix([[1/f0, 0, 0],
                        [0, f0/c0**2+b0**2*v0**2, b0**2*v0],
                        [0, b0**2*v0, b0**2]])
    assert sp.simplify(metric.inv()[2, 2] - (1/b0**2+c0**2*v0**2/f0)) == 0
    assert sp.simplify(metric.det()-b0**2/c0**2) == 0
    exact = exact_wall_zero(sp.Rational(1, 100))
    lo, hi = odd_wall_ritz(sp.Rational(1, 100), 12), odd_wall_ritz(sp.Rational(1, 100), 18)
    assert abs(hi[0]) < 1e-9 and hi[1] > 8 and max(abs(lo-hi)) < 1e-4
    return {"status": "PASS", "scope": "Exact odd mixed vector zero on the backreacted h=1/100 wall; two torus fibers, not total wall stability",
            "exact_profile": exact, "odd_Ritz_eigenvalues": hi.tolist(),
            "direct_metric_inverse_control": "g^yy=b^-2+c^2 Vphi^2/F, det g=b^2/c^2 on the (b,phi,y) block",
            "basis12_18_error": float(max(abs(lo-hi))),
            "interpretation": "Each torus fiber retains an odd mixed zero. Constant Wilson scalars remain zero by gauge holonomy and global connection flux constraint. The two shape modes were lifted in 11341; remaining sectors and global spectral positivity need separate checks."}


def cp_word_degree():
    def cyclic(w):
        return {w[i:]+w[:i] for i in range(len(w))}
    first = None
    for degree in range(1, 7):
        failures = []
        for bits in itertools.product("AB", repeat=degree):
            word = "".join(bits)
            if word[::-1] not in cyclic(word):
                failures.append(word)
        if failures:
            first = degree, failures[0]
            break
    assert first is not None and first[0] == 6
    coefficients = {}
    for terms in itertools.product((("AB", 1), ("BA", -1)), repeat=3):
        word = "".join(term[0] for term in terms)
        key = min(cyclic(word))
        coefficients[key] = coefficients.get(key, 0) + np.prod([term[1] for term in terms])
    assert coefficients["AABBAB"] == 3 and coefficients["AABABB"] == -3
    assert coefficients["ABABAB"] == 0
    # Cross-check against the earlier three-ray Bargmann owner, Pass9949--9956.
    u = sp.Matrix([1, 0, 0])
    v = sp.Matrix([1, 1, 0])/sp.sqrt(2)
    w = sp.Matrix([1, sp.I, 1])/sp.sqrt(3)
    P, Q, R = (z*z.conjugate().T for z in (u, v, w))
    comm = P*Q-Q*P
    assert sp.simplify(sp.trace(comm**3)) == 0
    triangle = sp.simplify(sp.trace(P*Q*R))
    assert triangle == (1+sp.I)/6
    return {"status": "PASS", "scope": "Two Hermitian three-family Gram spurions and single-trace invariants; standard CP-invariant structure, not a novel theorem of W33",
            "first_non_reversal_equivalent_trace_word_degree_in_Grams": first[0],
            "example_word": first[1],
            "commutator_identity": "Im Tr([A,B]^3) = -6 Im Tr(AABABB), by cyclic expansion and Hermitian word reversal",
            "Yukawa_field_degree": 12,
            "two_rank_one_projector_commutator_cubic": "0 for any two rank-one projectors, since their commutator has eigenvalues +i lambda,-i lambda,0",
            "three_ray_Bargmann_control": str(triangle),
            "bridge_boundary": "The earlier three-ray Bargmann phase carries orientation, but a W33 magic-axis point is not itself a qutrit state projector. An explicit map to full flavor Grams is required before any CP claim.",
            "implication": "A renormalizable quartic potential built only from these two Gram spurions cannot itself contain the Jarlskog CP-odd invariant. A CP-even square of that invariant occurs at still higher order; additional fields or dynamics are needed for controlled chiral selection."}


def scale_obstruction():
    return {"status": "PASS", "scope": "Dependency audit, not an exclusion of quantum dimensional transmutation",
            "dimensionless_inputs": ["W33 incidence", "b_wall/b_pole=5/4", "Dirac integer=1", "axion winding integers", "q/c=5"],
            "free_input": "A physical length or reference coupling scale is absent from these algebraic constraints.",
            "conclusion": "The displayed dimensionless spectra and wall relations alone cannot determine an absolute mass, Newton constant, or observed cosmological constant. A running boundary condition, measured input, or independent scale-setting dynamics is required."}


def zero_bare_lambda_wall():
    b, bw, h, q, c, rho, C, F = wall_data(sp.Rational(16, 45))
    assert sp.simplify(q-sp.Rational(4, 3)) == 0
    assert sp.simplify(rho) == 0
    assert sp.simplify(F-8*(b-1)*(5-2*b)/(45*b*b)) == 0
    assert sp.simplify(F.subs(b, bw)-sp.Rational(16, 225)) == 0
    assert sp.simplify(q/c*(1-1/bw)-1) == 0
    assert sp.simplify(4*sp.sqrt(F.subs(b, bw))/bw-sp.Rational(64, 75)) == 0
    assert sp.simplify(F.subs(b, 1)) == 0
    assert sp.simplify(sp.diff(F,b).subs(b,1)-2*c) == 0
    assert sp.simplify(sp.diff(F,b).subs(b,bw)-2*F.subs(b,bw)/bw) == 0
    assert sp.simplify(-sp.diff(F,b,2)/2-sp.diff(F,b)/b-rho-q*q/(2*b**4)) == 0
    assert sp.simplify(-sp.diff(F,b)/b-F/b**2-rho+q*q/(2*b**4)-h/b**2) == 0
    vector = exact_wall_zero(h)
    return {"status": "PASS", "scope": "Exact zero-bare-rho point in a supplied Einstein-Maxwell-axion-wall model; not a selected cosmological constant",
            "h": str(h), "q": str(q), "c": str(c), "rho": str(rho), "F": str(F),
            "F_wall": "16/225", "wall_tension": "64/75", "Dirac_integer": 1,
            "shape_gap_Rayleigh_bounds": ["512/1125", "32/45"],
            "mixed_vector_zero": vector,
            "boundary": "Setting bare rho to zero requires the supplied h=16/45; wall tension, matter coefficient, and overall physical scale are still inputs."}


def payload():
    sections = {
        "11342_complex_flavor": flavor_barriers(),
        "11343_nongaussian_mediator": nongaussian_mediator(),
        "11344_Ward_normal_block": ward_normal_freedom(),
        "11345_transverse_frame_rank": transverse_rank(),
        "11346_winding_vector_zero": winding_vectors(),
        "11347_CP_word_degree": cp_word_degree(),
        "11348_scale_obstruction": scale_obstruction(),
        "11349_zero_bare_Lambda_wall": zero_bare_lambda_wall(),
    }
    assert all(x["status"] == "PASS" for x in sections.values())
    return {"status": "PASS", "scope": "Eight finite or exact probes with their own explicit boundaries; no TOE completion", "sections": sections,
            "prior_owners": ["analysis/w33_pass11335_larger_complex_flavor_search.py", "analysis/w33_pass11336_gaussian_mediator_alignment.py", "analysis/w33_pass11337_quotient_hard_ward_columns.py", "analysis/w33_pass11338_native_boost_constraint_slice.py", "analysis/w33_pass11339_coupled_maxwell_metric_modes.py", "analysis/w33_pass11341_winding_lifted_wall.py"],
            "primary_sources": ["https://arxiv.org/abs/hep-th/9803116", "https://arxiv.org/abs/1410.7774", "https://cds.cern.ch/record/161795"]}


if __name__ == "__main__":
    result = payload()
    path = ROOT / "data/w33_pass11342_11349_eight_frontier_probes.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(result["status"], list(result["sections"]))
