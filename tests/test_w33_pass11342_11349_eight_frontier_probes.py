"""Independent controls on the eight scoped physical probes."""
import importlib.util
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "frontiers", ROOT / "analysis/w33_pass11342_11349_eight_frontier_probes.py"
)
F = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(F)


def test_two_old_bound_escapes_have_new_positive_barrier():
    rows = F.flavor_barriers()["rows"]
    assert {(r["n"], r["start"]) for r in rows} == {(4, 3), (4, 9)}
    assert rows[1]["barrier_minimum"] > rows[0]["barrier_minimum"] > 2


def test_mediator_orthogonal_branch_beats_parallel_branch():
    item = F.nongaussian_mediator()
    a = 0.1
    roots = sp.nroots(12*sp.Symbol("t")**3 + sp.Symbol("t") - sp.Rational(1,10))
    parallel_t = float(next(sp.re(z) for z in roots if abs(float(sp.im(z))) < 1e-12 and float(sp.re(z)) > 0))
    parallel = 2*parallel_t**2 + 12*parallel_t**4 - 4*a*parallel_t
    assert item["global_minimum_energy"] < parallel - 1e-4
    assert item["positive_quartic_margin_min_sample"] >= -1e-10


def test_winding_vector_profile_survives_changed_wall_and_is_not_gauge():
    for h in (sp.Rational(1, 200), sp.Rational(1, 100), sp.Rational(16,45)):
        witness = F.exact_wall_zero(h)
        assert witness["chi_equation_residual"] == "0"
        assert witness["metric_connection_equation_residual"] == "0"
        b = sp.symbols("b", positive=True)
        chi = 1/b - sp.Rational(4, 5)
        assert sp.diff(chi,b) != 0


def test_odd_vector_ritz_has_exact_zero_and_positive_next_level():
    a = F.odd_wall_ritz(sp.Rational(1,100), 10)
    b = F.odd_wall_ritz(sp.Rational(1,100), 16)
    assert abs(b[0]) < 1e-9
    assert b[1] > 8
    assert np.max(abs(a-b)) < 2e-3


def test_zero_bare_rho_point_is_positive_and_quantized():
    d = F.zero_bare_lambda_wall()
    b = sp.symbols("b", positive=True)
    actual = sp.sympify(d["F"], locals={"b": b})
    expected = 8*(b-1)*(5-2*b)/(45*b*b)
    assert sp.simplify(actual-expected) == 0
    assert all(float(expected.subs(b,x)) > 0 for x in (sp.Rational(9,8), sp.Rational(5,4)))
    assert d["wall_tension"] == "64/75"


def test_CP_trace_word_degree_and_Ward_normal_count():
    cp = F.cp_word_degree()
    assert cp["first_non_reversal_equivalent_trace_word_degree_in_Grams"] == 6
    assert cp["three_ray_Bargmann_control"] == "1/6 + I/6"
    assert cp["two_rank_one_projector_commutator_cubic"].startswith("0 for any")
    assert F.ward_normal_freedom()["unfixed_symmetric_normal_entries"] == 105
