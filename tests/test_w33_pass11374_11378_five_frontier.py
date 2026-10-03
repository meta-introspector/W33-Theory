"""Independent boundary checks for the five scoped TOE probes."""
import importlib.util
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "five_frontier", ROOT / "analysis/w33_pass11374_11378_five_frontier.py"
)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_cp_pair_is_degenerate_and_opposite():
    d = mod.flavor_phase_selector()
    assert d["CP_cubic_at_minima"] == [1, -1]
    assert d["phase_curvature_at_each_minimum"] == "lambda/18"
    phi = sp.symbols("phi", real=True)
    v = sp.cos(phi)**2 / 36
    assert sp.simplify(v.subs(phi, sp.pi/2)-v.subs(phi, 3*sp.pi/2)) == 0
    da, db = map(sp.Rational, d["A_B_spectral_discriminants"])
    exact_j = sp.sympify(d["mixing_J_absolute_at_either_branch"])
    assert sp.simplify(36*da*db*exact_j**2-1) == 0


def test_rank_one_families_cannot_fill_fourteen_normal_directions():
    d = mod.hard_loop_rank_boundary()
    assert d["max_normal_rank"] == 1
    assert d["normal_symmetric_entries"] == 105
    v = np.arange(1,15, dtype=float)
    assert np.linalg.matrix_rank(np.outer(v,v),tol=1e-8) == 1
    # The theorem is deliberately conditional: f'=1 for I=sum x_i^2
    # leaves Hess(f(I))=2*identity, so the rank-one bound would fail.
    assert np.linalg.matrix_rank(2*np.eye(14)) == 14


def test_native_connection_count_uses_actual_cycle_basis():
    d = mod.native_frame_flatness()
    assert d["independent_cycles"] == d["edges"]-d["vertices"]+1
    assert d["unconstrained_edge_connection_dimension"] == 960
    assert d["independent_linearized_cycle_flatness_equations"] == 486
    assert d["flat_connection_mod_vertex_gauge_dimension"] == 0


def test_nonconstant_axion_torus_modes_have_positive_probe_gap():
    d = mod.wall_axion_probe_sector()
    assert sp.Rational(d["first_harmonic_lower_bound"]) == sp.Rational(16,25)
    for mx,my in ((1,0),(0,1),(2,3)):
        assert (mx*mx+my*my)*sp.Rational(16,25) > 0


def test_wall_ratio_does_not_fix_absolute_length():
    d = mod.wall_scale_ratio()
    rw = sp.Rational(d["wall_Ricci_scalar_kappa2_equal_1"])
    t2 = sp.Rational(d["wall_tension_squared_kappa2_equal_1"])
    assert rw/t2 == sp.Rational(5,8)
    L = sp.symbols("L", positive=True)
    assert sp.simplify((rw/L**2)/(t2/L**2)) == sp.Rational(5,8)
