"""Independent controls for the exact physical-frontier claims."""
import sys
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11361_11368_physical_frontier as F


def test_three_ray_cp_orientation_and_weight_boundary():
    u = sp.Matrix([1, 0, 0])
    v = sp.Matrix([1, 1, 0]) / sp.sqrt(2)
    w = sp.Matrix([1, sp.I, 2]) / sp.sqrt(6)  # different third ray from the producer
    P = [z*z.conjugate().T for z in (u, v, w)]
    X = sp.Matrix.hstack(u, v, w)
    A = P[0] + 2*P[1] + 4*P[2]
    B = 2*P[0] + 3*P[1] + 5*P[2]
    J = sp.simplify(sp.trace((A*B-B*A)**3)/sp.I)
    barg = sp.im(sp.trace(P[0]*P[1]*P[2]))
    assert J == -6*sp.simplify((X.conjugate().T*X).det())*barg*(-1)*(-3)*(-2)
    assert J == 1
    conjugate_J = sp.simplify(sp.trace((A.conjugate()*B.conjugate()-B.conjugate()*A.conjugate())**3)/sp.I)
    assert conjugate_J == -J
    parallel_weights = 2*A
    assert sp.trace((A*parallel_weights-parallel_weights*A)**3) == 0


def test_equal_lapse_hidden_null_and_generic_second_prime():
    B, C = F.graph_data()
    U = abs(B)
    rng = np.random.default_rng(11363)
    eta = rng.permutation(80).astype(np.int64)
    N = rng.integers(100, 201, size=80, dtype=np.int64)
    r = B.T@eta
    Geq = U@(r[:, None]*C)
    assert np.all(Geq.T@eta == 0)
    assert F.rank_mod(Geq, 1000033)[0] == 77
    L = U.T@N
    p = 1000033
    d = (r%p)*np.array([pow(int(x), -1, p) for x in L], dtype=np.int64)%p
    G = (U@(d[:, None]*C))%p
    assert F.rank_mod(G, p)[0] == 78


def test_zero_bare_rho_does_not_make_the_wall_flat():
    d = F.wall_family_and_scale()
    b = sp.symbols("b", positive=True)
    R = sp.sympify(d["Ricci_scalar_at_zero_bare_rho"], locals={"b": b})
    assert sp.simplify(R.subs(b, sp.Rational(5, 4))-sp.Rational(512, 1125)) == 0
    rows = F.wound_wall_spectrum()
    assert abs(rows["odd_schur_Ritz_basis18"][0]) < 1e-9
    assert rows["odd_schur_Ritz_basis18"][1] > 10


def test_supplied_heavy_loop_is_normal_to_live_family_orbit():
    d = F.heavy_radial_loop()
    assert d["Ward_column_residual"] < 1e-10
    assert abs(d["normal_eigenvalue_g1"]-1/(16*np.pi**2)) < 1e-12


def test_native_octagon_and_corrected_public_scope():
    h = F.transverse_cycle_holonomy()
    assert len(h["native_cycle_edges"]) == 8
    d = F.physical_claim_scope()
    assert d["corrected_integer"] == 122
    assert "not a physical parameter fit" in d["scope"]
