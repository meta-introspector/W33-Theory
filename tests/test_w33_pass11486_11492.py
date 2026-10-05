"""Regression for Passes 11486-11492 (magic-axis law for z1 in Im(M - I); exact depth 8; six qutrits; cut audit; time-odd
Molien series; time-odd degree by dimension)."""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def load(name):
    return json.load(open(ROOT / "data" / name))


def test_11486_uniform_proof_n2():
    d = load("w33_pass11486_blind_classes_phases.json")
    b = d["n2"]
    assert b["blind classes"] == 168
    for k in ("e = 3Q + s^3", "l_of_v_zero", "l_prop_Bv", "B|R nondegenerate", "F' true", "theorem sound",
              "F' fully proved", "Lemma 2 criterion = exhaustive verdicts", "every admissible L fixes v"):
        assert b[k] == 168, k
    p = d["n2_partial_hyperplane"]
    assert p["partial hyperplane classes"] == p["theorem sound"] == p["F' fully proved"] == 48
    t = d["n2_transvection"]
    assert t["Pass 11457 partial classes"] == t["transvection, Mz1 = z1, Wall q1 != 0"] == t["F' true (theorem sound)"] == 24
    assert 352 + 704 + 48 + 168 + 24 == 1296


def test_11486_n3_n4():
    d = load("w33_pass11486_blind_classes_phases.json")
    c = d["n3"]["counts"]
    assert c["blind orbits"] == 33 and c["theorem sound on decided"] == c["decided"]
    n4 = d["n4_unreached"]
    assert n4["count"] == 317 and sum(n4["classes"].values()) == 317
    assert n4["classes"]["blind (q = 0): fully proved by Pass 11486"] == 68
    assert n4["classes"]["partial hyperplane (q != 0): fully proved by Pass 11486"] == 90
    assert n4["mass_fraction"]["outside the phase theorems (z1 not in R_M, Rad != 0)"] < 0.1
    a3 = d["n3_all_theorems"]
    assert a3["soundness"]["sound"] == a3["soundness"]["decided"] == 102 and abs(a3["fully_proved_mass"] - 0.7968) < 1e-3


def test_11486_transvection_fast():
    import w33_pass11350_linear_decider as L
    import w33_pass11486_blind_classes_phases as B
    D = L.Decider(1)
    M = np.array([[1, 0], [1, 1]], dtype=np.int64)           # a transvection of the one-qutrit plane fixing z1 = (0, 1)
    h = B.transvection_case(D, M)
    assert h is None or isinstance(h["holds"], bool)


def test_11487_extension():
    d = load("w33_pass11487_magic_axis_extension.json")
    n2 = d["n2"]
    assert n2["sound"] == n2["hyperplane classes"] == n2["F' true on class (all omega(a,v)!=0 frames violate)"] == 752
    assert n2["fully certified"] == 704 and n2["no level hyperplane (q = 0: magnitude-blind)"] == 168
    s3 = d["n3_soundness"]
    assert s3["sound"] == s3["checked"] > 0
    assert d["n4_unreached"]["count"] == 317


def test_11487_level_hyperplane_fast():
    import w33_pass11350_linear_decider as L
    import w33_pass11487_magic_axis_extension as X
    import w33_pass11421_four_qutrits as F4
    lw = F4.LightWeyl(1)
    G = np.eye(3, dtype=complex)
    m = X.weyl_moduli(lw, G)                       # identity: only W(0) has a nonzero coefficient
    assert abs(m[0] - 1) < 1e-12 and np.allclose(m[1:], 0)
    assert L is not None


def test_11488_exact_depth8():
    from fractions import Fraction
    d = load("w33_pass11488_exact_depth8.json")
    lv = d["levels"]
    for k, v in {"1": "11/12", "2": "19/32", "3": "347/768", "4": "623/2048", "5": "11003/49152", "6": "20327/131072",
                 "7": "353927/3145728"}.items():
        assert lv[k]["P"] == v
    assert [lv[str(k)]["distinct_operators"] // 216 for k in range(1, 7)] == [1, 5, 26, 136, 772, 4440]
    P8 = Fraction(lv["8"]["P"])
    assert 0 < P8 < Fraction(353927, 3145728)
    assert d["overlap_gap"]["min_reversible"] > 3 - 1e-9 and d["overlap_gap"]["max_violator"] < 3 - 1e-4


def test_11489_six_qutrits():
    d = load("w33_pass11489_six_qutrits.json")
    e = d["n6"]["estimate"]
    assert 0.05 < e["P_estimate"] < 0.25 and e["stderr"] < 0.03
    assert min(d["completed_by_kind"].values()) >= 50


def test_11490_cut_audit():
    d = load("w33_pass11490_cut_audit.json")
    for k, r in d["deep_words"].items():
        assert r["deciders_disagree"] == 0
        assert r["overlap_gap"]["min_deficit_violator"] > 1e-5
        if r["overlap_gap"]["max_deficit_reversible"] is not None:
            assert r["overlap_gap"]["max_deficit_reversible"] < 1e-12
        assert r["j6_gap"]["violators_below_1e9"] == 0
    assert min(r["j6_gap"]["min_j6_violator"] for r in d["deep_words"].values()) < 1e-8     # the thin margin
    ls = d["level_statistics"]
    assert ls["between_1e12_and_1e6"] == 0 and ls["min_above_1e9"] > 1e-3


def test_11491_molien_closed_forms():
    d = load("w33_pass11491_time_odd_molien.json")
    assert all(d["proved"].values())
    assert d["periods"]["values_checked"] >= d["periods"]["values_needed"] == 60
    assert d["sequences"]["odd"][:16] == [0] * 6 + [1, 2, 3, 5, 8, 11, 16, 22, 29, 38]
    assert d["lowest_odd_degree"] == 6 and d["odd_not_multiples_of_h6"]["first_new_odd_degree"] == 7
    assert (d["design_check"]["D2"], d["design_check"]["D3"]) == (1, 2)
    eg = d["explicit_generators"]
    assert [eg[k]["rank"] for k in ("6", "7", "8")] == [1, 2, 3]
    assert all(eg[k]["first_dropped_sv"] < 1e-12 < eg[k]["last_kept_sv"] for k in ("6", "7", "8"))
    lo, hi = eg["g7_p4p2p_residual_vs_h6_all_24_oriented_triples"]
    assert lo > 0.02 and hi - lo < 1e-9


def test_11491_series_identity_symbolic():
    import sympy as sp
    import w33_pass11491_time_odd_molien as M
    t = M.t
    F = M.FORMS
    assert sp.simplify(F["total"] - F["even"] - F["odd"]) == 0
    assert sp.simplify(F["twisted"] - F["even"] + F["odd"]) == 0
    for f in F.values():                            # the shared reciprocity H(1/t) = -t^3 H(t)
        assert sp.simplify(f.subs(t, 1 / t) + t ** 3 * f) == 0


def test_11492_by_dimension():
    d = load("w33_pass11492_time_odd_by_dimension.json")
    assert d["lowest_odd_degree_by_d"] == {"2": 9, "3": 6, "5": 6, "7": 5, "11": 5}
    assert d["count_at_lowest_odd_degree"] == {"2": 1, "3": 1, "5": 4, "7": 2, "11": 16}
    assert d["control_p5_two_constructions_agree"] is True
    assert all(v["proved"] for v in d["d2"]["closed_forms"].values())
    assert all(v["proved"] for v in d["d3"]["closed_forms"].values())
    o, e = d["d2"]["odd"], d["d2"]["even"]
    assert all(o[k] == (e[k - 9] if k >= 9 else 0) for k in range(len(o)))   # qubit: odd = pseudoscalar x even
    assert d["d3"]["odd"][:41] == load("w33_pass11491_time_odd_molien.json")["sequences"]["odd"]
    tq = d["two_qutrits"]
    assert tq["group_order"] == 81 * 51840 and tq["odd"][:9] == [0, 0, 0, 0, 0, 0, 1, 6, 25]
    assert d["two_qutrit_restriction"]["max_dev_from_1_over_108"] < 1e-8
