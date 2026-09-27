"""Regression for Pass 11024: Z3-graded superpotential on the Z3xZ3 parity vacua, explicit SUSY roots, the
hidden-representation index (massless charge-1/3 states), hidden condensates, and the exact-lattice certificate of
the Pass 10978-10980 integer programs."""
import importlib.util
import json
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
C = json.loads((ROOT / "data" / "w33_pass11024_z3xz3_susy_vacua_fractional_charge_index.json").read_text())
SPEC = importlib.util.spec_from_file_location("exo", ROOT / "analysis" / "w33_exact_monomial_orders.py")
EXO = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(EXO)
C3 = "z3z3|Z3Z3_0001_c3__SM_20260934_739"
C4 = "z3z3|Z3Z3_0001_c4__SM_20260935_3729"


def test_z3_graded_superpotential_all_orders():
    A = C["A"]["summary"]
    assert A["vacua"] == 28 and A["z3_graded_all_orders"] == 28
    assert A["k=1"] == 9 and A["k=2"] == 19
    assert A["prediction_matches_10979"] == 28
    for vs in C["A"]["vacua"].values():
        for v in vs:
            assert v["lattice_index"] == 1 and v["periodic"] and v["z3_graded"]


def test_susy_roots_exist_everywhere():
    B = C["B"]["summary"]
    assert B["vacua"] == 28 and B["with_verified_susy_roots"] == 28
    lo, hi = B["max_abs_median_range"]
    assert 0.6 < lo <= hi < 1.3


def test_fractional_charge_index():
    for name in (C3, C4):
        ix = C["C"]["models"][name]["index"]
        assert ix["fractional_states"] == 48 and ix["unconfinable_fractional"] == 6
        beta = C["C"]["models"][name]["beta"]
        assert all(b["infrared_free"] for b in beta.values() if b.get("N"))
    # the index is not automatic: it vanishes in many other W(3,3) models (control)
    assert C["C"]["controls"]["Z6-I/II"]["zero_unconfinable"] > 0
    assert C["C"]["models"]["z3z3|Z3Z3_0001_c2__SM_20260933_4158"]["index"]["fractional_states"] == 0


def test_su4_meson_rescue_never_pairs_everything():
    for name in (C3, C4):
        r = C["D"]["su4_rescue"][name]
        assert r["excess"] == 42 and r["same_class"] == 36 and r["need_su4"] == 6
        assert not any(k.startswith("all_paired") for k in r["tally"])
        assert all(v["massless_unconfined_fractional"] >= 6 for v in r["vacua"])


def test_ilp_certificates():
    E = C["E"]
    assert E["p10979"]["differences"] == 0 and E["p10979"]["checked"] == 72
    assert E["p10978"]["differences"] == 0 and E["p10978"]["c1063_first_linear"] == [4, "n_36"]
    assert E["p10980"]["identical_summary"]
    assert E["fast_parity"]["tests"] == E["fast_parity"]["agree"] == 217


def test_exact_solver_known_positive_and_brute_force():
    # toy lattice: charges (u1 | Z3), W charge (0 | 1/3); field a=(1|0), b=(-1|1/3), c=(0|2/3)
    vecs = [[Fraction(1), Fraction(0)], [Fraction(-1), Fraction(1, 3)], [Fraction(0), Fraction(2, 3)]]
    E = EXO.ExactOrders(vecs, 1)
    for tgt in ([0, Fraction(1, 3)], [2, 0], [0, 0], [-3, Fraction(2, 3)]):
        assert E.order(tgt) == EXO.brute_force_order(vecs, 1, tgt, 12)
    # known positive in the coupling convention: X.Y.S^e allowed iff sum e v = w - v(X) - v(Y)
    assert E.coupling([[Fraction(-1), Fraction(0)]], [0, Fraction(1, 3)]) == E.order([1, Fraction(1, 3)])


def test_hidden_condensate_census_and_the_near_miss():
    D = C["D"]["census"]
    assert D["z3z3|Z3Z3_0001_c2__SM_20260933_4158"]["census"]["verdict"].startswith("closed: Farkas")
    c1 = D["z3z3|Z3Z3_0001_c1__SM_20260932_2822"]
    assert c1["curing_slots"] == [1] and c1["census"]["realizable_curing_rays"] == 4
    assert c1["composite_vacua"]["summary"] == {"viable": 1, "none": 3}
    (cand,) = c1["candidates"]
    assert cand["parity_composites_even"] == "viable" and cand["parity_aligned"] == "none"
    assert cand["light_fractional_multiplets"] == 158 and cand["W_on_vacuum"] is None


def test_su4_condensed_parity_vacua_keep_massless_fractional_charges():
    D = C["D"]["census"]
    for name, n_rays in (("z3z3|Z3Z3_0001_c3__SM_20260934_739", 39), ("z3z3|Z3Z3_0001_c4__SM_20260935_3729", 15)):
        r = D[name]
        assert r["curing_slots"] == [0]                      # only breaking hidden SU(4) can cure the index
        assert r["index_with_one_factor_broken"]["0"]["fractional"] == 0
        assert r["index_with_one_factor_broken"]["1"]["fractional"] == 48
        assert r["census"]["realizable_curing_rays"] == n_rays
        assert r["composite_vacua"]["summary"] == {"viable": n_rays}
        for c in r["candidates"]:
            assert c["parity_aligned"] == "viable" and c["W_on_vacuum"] is None
            assert c["light_fractional_multiplets"] >= 120
            assert c["light_fractional_protected_by_unbroken_symmetry"] >= 106
            assert c["unbroken_u1s"] == 1
            order = 1
            for d in c["unbroken_discrete_invariant_factors"]:
                order *= d
            assert order in (3 ** 4, 3 ** 6)   # a 3-group: Z3^2xZ9^2 or Z3^4xZ9 (3^6), Z3^2xZ9 (3^4)
