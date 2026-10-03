"""Eight scoped physical-frontier checks on the current W33 certificates.

The CP formula and metric-cycle rank statements are exact.  The wall Ritz
numbers are trial-space values; none of these checks is a TOE prediction.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def rays_to_flavor():
    # Exactly the three normalized qutrit rays of the prior Pass 9949--9956.
    rays = [sp.Matrix([1, 0, 0]),
            sp.Matrix([1, 1, 0]) / sp.sqrt(2),
            sp.Matrix([1, sp.I, 1]) / sp.sqrt(3)]
    X = sp.Matrix.hstack(*rays)
    P = [u * u.conjugate().T for u in rays]
    G = X.conjugate().T * X
    bargmann = sp.simplify(sp.trace(P[0] * P[1] * P[2]))
    det_gram = sp.simplify(G.det())
    a = (1, 2, 4)
    b = (2, 3, 5)
    A = sum((a[i] * P[i] for i in range(3)), sp.zeros(3))
    B = sum((b[i] * P[i] for i in range(3)), sp.zeros(3))
    minors = [(a[i] * b[j] - a[j] * b[i]) for i, j in ((0, 1), (0, 2), (1, 2))]
    J = sp.simplify(sp.trace((A * B - B * A) ** 3) / sp.I)
    rhs = sp.simplify(-6 * det_gram * sp.im(bargmann) * sp.prod(minors))
    assert bargmann == (1 + sp.I) / 6 and det_gram == sp.Rational(1, 6)
    assert minors == [-1, -3, -2] and J == rhs == 1
    assert A.det() > 0 and B.det() > 0
    assert sp.discriminant(A.charpoly().as_expr(), sp.Symbol("lambda")) != 0
    assert sp.discriminant(B.charpoly().as_expr(), sp.Symbol("lambda")) != 0
    # Real positive weights give full-rank positive Yukawa Grams.
    An = np.array(A.evalf(), dtype=complex)
    Bn = np.array(B.evalf(), dtype=complex)
    ea, Ua = np.linalg.eigh(An)
    eb, Ub = np.linalg.eigh(Bn)
    V = Ua.conj().T @ Ub
    jmix = float(np.imag(V[0, 0] * V[1, 1] * V[0, 1].conjugate() * V[1, 0].conjugate()))
    assert min(ea) > 0 and min(eb) > 0 and abs(jmix) > 1e-6
    return {
        "status": "PASS", "scope": "Exact algebraic ray-to-full-rank Gram construction, not a dynamical or observed CKM prediction",
        "identity": "Tr([A,B]^3)/i=-6 det(G) ImTr(P1P2P3) c12 c13 c23; c_ij=a_i b_j-a_j b_i",
        "proof": "A=X diag(a) Xdag, B=X diag(b) Xdag. Set K=diag(a)Gdiag(b)-diag(b)Gdiag(a). Then [A,B]=X K Xdag, det[A,B]=det(G)det(K), det(K)=-2i c12 c13 c23 Im(G12 G23 G31), and Tr([A,B]^3)=3det[A,B].",
        "nonzero_criterion": "The cubic is nonzero exactly when the three rays span C3, their Bargmann loop has nonzero imaginary part, and all three pairwise weight minors are nonzero.",
        "Bargmann": str(bargmann), "det_Gram": str(det_gram), "weight_minors": minors,
        "commutator_cubic_over_i": str(J), "A_eigenvalues": ea.tolist(), "B_eigenvalues": eb.tolist(),
        "mixing_J": jmix,
        "boundary": "The rays are prior-owned; the positive weights are supplied. Their conjugate orientation reverses J at equal CP-even energy. No W33 dynamics selects these weights, orientation, masses or observed mixing."
    }


def heavy_radial_loop():
    # Compute the radial invariant and family tangents in ONE live chart.
    import w33_pass11287_quotient_self_energy_matrix as Q
    p, N, _, *_ = Q.geometry()
    z = N.conj().T @ p
    v = np.sqrt(2) * np.r_[z.real, z.imag]  # d(phi^dag phi)/dx at the vacuum
    F = Q.Q.L.generators()[78:]
    T = np.array([np.sqrt(2) * np.r_[(N.conj().T @ (1j * f @ p)).real,
                                      (N.conj().T @ (1j * f @ p)).imag] for f in F]).T
    assert T.shape == (22, 8) and np.linalg.matrix_rank(T, tol=1e-9) == 8
    assert abs(v @ v - 2) < 1e-10 and np.max(abs(T.T @ v)) < 1e-10
    # A supplied real heavy scalar has mass^2=M^2+g(I-I0), I=phi^dag phi.
    # At mu^2=M^2/e its CW tadpole vanishes and f''(M^2)=2.
    m2, mu2, g = sp.symbols("m2 mu2 g", positive=True)
    x = sp.symbols("x", positive=True)
    f = x*x * (sp.log(x / mu2) - sp.Rational(3, 2))
    assert sp.simplify(sp.diff(f, x).subs({x: m2, mu2: m2 / sp.E})) == 0
    assert sp.simplify(sp.diff(f, x, 2).subs({x: m2, mu2: m2 / sp.E})) == 2
    H = np.outer(v, v) / (32 * np.pi**2)  # g=1, one real scalar
    assert np.max(abs(H @ T)) < 1e-12
    assert abs(np.trace(H) - 1 / (16 * np.pi**2)) < 1e-12
    return {
        "status": "PASS", "scope": "One specified invariant heavy-scalar 1PI normal contribution, not the 105-entry native E6 hard block",
        "invariant": "I=phi^dag phi; mS^2=M^2+g(I-I0), one real spectator; mu^2=M^2/e",
        "loop": "V1=mS^4[log(mS^2/mu^2)-3/2]/(64pi^2)",
        "normal_Hessian": "Hloop=g^2 v v^T/(32pi^2), v=dI/dx, ||v||^2=2; one nonzero eigenvalue g^2/(16pi^2)",
        "radial_norm_squared": float(v @ v), "Ward_column_residual": float(np.max(abs(H @ T))),
        "normal_eigenvalue_g1": float(np.trace(H)),
        "boundary": "The heavy field and g are supplied. This computes one scheme-dependent normal direction and proves Ward columns do not fix it; it is not native UV matching or a mass prediction."
    }


def graph_data():
    from w33_pass11289_cycle_gram_gluing_flux import cycle_basis
    B, C = cycle_basis()
    assert B.shape == (80, 160) and C.shape == (160, 81)
    assert np.all(B @ C == 0)
    return B.astype(np.int64), C.astype(np.int64)


def rank_mod(A, prime=1000003):
    A = np.array(A, dtype=np.int64) % prime
    row = 0
    pivots = []
    for col in range(A.shape[1]):
        k = next((i for i in range(row, A.shape[0]) if A[i, col]), None)
        if k is None:
            continue
        A[[row, k]] = A[[k, row]]
        A[row] = A[row] * pow(int(A[row, col]), -1, prime) % prime
        for i in range(row + 1, A.shape[0]):
            c = int(A[i, col])
            if c:
                A[i] = (A[i] - c * A[row]) % prime
        pivots.append(col)
        row += 1
        if row == A.shape[0]:
            break
    return row, pivots


def metric_cycle_rank():
    B, C = graph_data()
    U = abs(B)
    rng = np.random.default_rng(11363)
    eta = rng.permutation(80).astype(np.int64)
    N = rng.integers(100, 201, size=80, dtype=np.int64)
    r = B.T @ eta
    L = U.T @ N
    prime = 1000003
    assert np.max(abs(r / L)) < 1
    # At this *actual stationary point*, j/sqrt(1+j^2)=r/L and
    # C^T L j/sqrt(1+j^2)=C^T B^T eta=0.
    d = (r % prime) * np.array([pow(int(t), -1, prime) for t in L], dtype=np.int64) % prime
    G = (U @ (d[:, None] * C)) % prime
    rank, pivots = rank_mod(G, prime)
    checker = np.r_[np.ones(40, dtype=np.int64), -np.ones(40, dtype=np.int64)]
    assert rank == 78
    assert np.all((G.T @ N) % prime == 0) and np.all((G.T @ checker) % prime == 0)
    # Equal lapses have an additional exact null vector eta: on each edge
    # (eta_i+eta_j)(eta_j-eta_i)=eta_j^2-eta_i^2.
    Geq = U @ (r[:, None] * C)
    equal_rank, _ = rank_mod(Geq, prime)
    assert np.all(Geq.T @ eta == 0) and equal_rank == 77
    assert np.linalg.matrix_rank(np.column_stack([np.ones(80), checker, eta])) == 3
    return {
        "status": "PASS", "scope": "Exact rank theorem for the named native pair-metric cycle model, not the preferred multivielbein Dirac algebra",
        "graph": "80 native sites, 160 edges, 81 integral cycles",
        "generic_rank_mod_1000003": rank, "generic_pivot_columns": pivots,
        "equal_lapse_rank_mod_1000003": equal_rank,
        "hidden_equal_lapse_null": "eta; because (U^T eta)_e(B^T eta)_e=(B^T eta^2)_e and C^T B^T=0",
        "generic_nulls": "N and point/line checkerboard; both exact",
        "eta": eta.tolist(), "N": N.tolist(),
        "stationary_construction": "L=U^T N>0, r=B^T eta, j=(r/L)/sqrt(1-(r/L)^2), p=Bj. The cycle stationarity is exact and A=C^T diag(L/(1+j^2)^(3/2)) C>0. Thus H_lapse=-G A^-1 G^T has rank78 at this point; a nonzero rational minor makes this rank generic in the algebraic family.",
        "boundary": "Rank78 is a ghost-constraint obstruction for the named pair-metric cycle action. It does not prove a full Dirac count, nor refute the distinct vertex-frame multivielbein action."
    }


def wound_wall_spectrum():
    from w33_pass11342_11349_eight_frontier_probes import exact_wall_zero, odd_wall_ritz
    h = sp.Rational(16, 45)
    exact = exact_wall_zero(h)
    low = odd_wall_ritz(h, 12)
    high = odd_wall_ritz(h, 18)
    assert abs(high[0]) < 1e-9 and high[1] > 10
    assert np.max(abs(low - high)) < 1e-4
    return {
        "status": "PASS", "scope": "Exact odd mixed zero and converged Ritz trial-space values at zero bare rho; no full wall stability claim",
        "h": str(h), "exact_profile": exact,
        "odd_schur_Ritz_basis18": high.tolist(), "basis12_to18_error": float(np.max(abs(low - high))),
        "boundary": "Ritz eigenvalues are upper bounds within a restricted odd m=0 Schur block and cannot exclude negative modes in the full breathing, wall-bending, axion or tensor sectors."
    }


def wall_family_and_scale():
    b, q = sp.symbols("b q", positive=True)
    h = sp.Rational(4, 5) * (q*q - q)
    c = q / 5
    rho = sp.factor(q*q / 2 - h - 2*c)
    C = rho / 3 + q*q / 2 + h
    F = sp.factor(C/b - rho*b*b/3 - q*q/(2*b*b) - h)
    expected = q*(b-1)*(q*(3*b**3+3*b*b-21*b+15)-4*b**3-4*b*b+20*b)/(30*b*b)
    assert sp.simplify(F-expected) == 0
    assert sp.simplify(rho-q*(4-3*q)/10) == 0
    H1 = sp.factor((expected*30*b*b/(q*(b-1))).subs(q, 1))
    H4 = sp.factor((expected*30*b*b/(q*(b-1))).subs(q, sp.Rational(4, 3)))
    assert sp.simplify(H1-(15-b-b*b-b**3)) == 0
    assert sp.simplify(H4-(20-8*b)) == 0
    R = sp.factor(4*rho+2*h/b**2)
    assert sp.simplify(R.subs(q, sp.Rational(4, 3))-32/(45*b*b)) == 0
    return {
        "status": "PASS", "scope": "Exact one-parameter supplied winding-wall family and bare-rho uniqueness, not observed vacuum-energy selection",
        "h_of_q": str(h), "rho_of_q": str(rho), "F_of_b_q": str(F),
        "positive_F_proof": "For 1<=q<=4/3 and 1<b<=5/4, F=q(b-1)H/(30b^2). H is affine in q; H(q=1)=15-b-b^2-b^3>0 and H(q=4/3)=20-8b>0.",
        "unique_positive_q_zero_bare_rho": "q=4/3, h=16/45",
        "Ricci_scalar_at_zero_bare_rho": str(sp.factor(R.subs(q, sp.Rational(4, 3)))),
        "scale_boundary": "The dimensionless solution fixes neither the physical length L nor mass units: curvature scales as L^-2 and vacuum-energy density as L^-4. Bare rho=0 still has strictly positive Ricci scalar from axion stress."
    }


def transverse_cycle_holonomy():
    B, C = graph_data()
    edges = np.flatnonzero(C[:, 0])
    assert len(edges) == 8 and np.all(abs(C[edges, 0]) == 1)
    degrees = abs(B[:, edges]).sum(axis=1)
    assert np.count_nonzero(degrees == 2) == 8 and np.count_nonzero(degrees) == 8
    sx = sp.Matrix([[0, 1], [1, 0]])
    sy = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    sz = sp.diag(1, -1)
    Kx, Ky = sx/2, sy/2
    assert Kx*Ky-Ky*Kx == sp.I*sz/2
    t = sp.symbols("t", real=True)
    I = sp.eye(2)
    # Keep the second-order term of each exponential; omitting it would
    # introduce spurious -Kx^2-Ky^2 terms into the loop coefficient.
    def E(K):
        return I + t*K + t*t*K*K/2
    product = E(Kx)*E(Ky)*E(-Kx)*E(-Ky)
    second = product.applyfunc(lambda z: sp.expand(z).coeff(t, 2))
    assert second == Kx*Ky-Ky*Kx
    return {
        "status": "PASS", "scope": "Exact noncommuting transverse-boost holonomy on a native octagon; necessary flatness condition, not full Dirac closure",
        "native_cycle_edges": edges.tolist(),
        "holonomy_order_two": "exp(tKx)exp(tKy)exp(-tKx)exp(-tKy)=I+t^2 i sigma_z/2+O(t^3)",
        "vertex_frame_flatness": "For edge transports Lambda_i^-1 Lambda_j, ordered product on every closed graph cycle is exactly I. Independent transverse edge boosts with the displayed octagon product therefore cannot come from vertex frames without compensating rotations.",
        "boundary": "This rules out independent-edge transverse boost assignments in a vertex-frame construction. It neither counts secondary constraints nor proves generic multivielbein ghost freedom."
    }


def physical_claim_scope():
    site = (ROOT / "docs/index.html").read_text(encoding="utf-8")
    fallback = (ROOT / "docs/404.html").read_text(encoding="utf-8")
    paper = (ROOT / "papers/forty_points/main.tex").read_text(encoding="utf-8")
    for page in (site, fallback):
        assert "Legacy SRG → Physics Formula Catalogue" in page
        assert "40×3 + 4 − 2 = <strong>122</strong>, not 127" in page
        assert "246 (dimensionless)" in page and "125 (dimensionless)" in page
        assert "−127" not in page
        assert "A graph calculation alone does not supply a field-theory map" in page
        assert "physical labels are not automatically verified by arithmetic checks" in page
    assert "passes11361-11368-ray-cp-gravity-wall-scope" in site
    assert "sec09_physics" in paper
    return {
        "status": "PASS", "scope": "Source-level claim audit, not a physical parameter fit",
        "conflict": "The legacy site labelled dimensionless integer coincidences in GeV and log10(CC/MPl^4)=-127 as exact predictions; its own arithmetic gives 40*3+4-2=122. The current paper and Pass 11348 leave mass/CC scale dynamics open.",
        "corrected_integer": 122,
        "resolution": "The legacy catalogue now calls these arithmetic correspondences, corrects 127 to 122, and states that GeV and vacuum-energy identifications require missing dynamics.",
        "boundary": "No observed mass, Higgs VEV, or vacuum energy is derived by the corrected wording."
    }


def payload():
    sections = {
        "11361_ray_CP_bridge": rays_to_flavor(),
        "11362_invariant_hard_radial_loop": heavy_radial_loop(),
        "11363_11364_exact_metric_cycle_rank": metric_cycle_rank(),
        "11365_wound_wall_spectrum": wound_wall_spectrum(),
        "11366_wall_family_scale": wall_family_and_scale(),
        "11367_transverse_octagon_holonomy": transverse_cycle_holonomy(),
        "11368_legacy_claim_scope": physical_claim_scope(),
    }
    assert all(v["status"] == "PASS" for v in sections.values())
    return {"status": "PASS", "scope": "Eight reserved passes represented by seven certificate sections; five targets and three additional probes remain scoped, not a complete TOE",
            "sections": sections,
            "prior_owners": ["analysis/w33_pass9949_9956_bargmann_orientation_readout.py",
                             "analysis/w33_pass11337_quotient_hard_ward_columns.py",
                             "analysis/w33_pass11323_native_cycle_lapse_hessian.py",
                             "analysis/w33_pass11342_11349_eight_frontier_probes.py"],
            "primary_sources": ["https://cds.cern.ch/record/161795", "https://arxiv.org/abs/1410.7774",
                                "https://arxiv.org/abs/2004.01868"]}


if __name__ == "__main__":
    result = payload()
    target = ROOT / "data/w33_pass11361_11368_physical_frontier.json"
    target.write_text(json.dumps(result, indent=2) + "\n")
    print(result["status"], list(result["sections"]))
