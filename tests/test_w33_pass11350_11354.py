"""Regression for Passes 11350-11354 (frozen certificates plus fast exact recomputations)."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def load(name):
    return json.loads((ROOT / "data" / name).read_text())


def test_11350_linear_decider_certificate():
    d = load("w33_pass11350_linear_decider.json")
    assert d["g_is_e_plus_t_k"] and d["two_qutrit_classes"] == 51840
    assert d["agree_with_pass_11330_counts"] == 51838 and d["undecided_cap"] == 2 and d["mismatches"] == []
    assert d["frame_spot_checks_agree"] == 300


def test_11350_decider_fast():
    import w33_pass11330_orbit_census as O
    import w33_pass11350_linear_decider as L
    D = L.Decider(2)
    Ms, keys = O.all_symplectic(D.wl)
    counts = np.load(ROOT / "data" / "w33_pass11330_class_counts.npy")
    rng = np.random.default_rng(7)
    for i in list(rng.choice(len(Ms), 25, replace=False)) + list(np.flatnonzero(counts > 0)[:15]):
        g = D.good_frames(Ms[i])
        if g is not None:
            assert int((~g).sum()) == counts[i]


def test_11351_one_ninth():
    d = load("w33_pass11351_one_ninth.json")
    assert d["noncollinear_classes"] == 34992 and d["classes_with_k_nonzero_solutions"] == 0
    t = d["symplectic_solution_count_vs_verdict"]
    assert t["3"] == {"bad": 3888, "good": 3888} and t["1"] == {"good": 23328}
    assert all("bad" not in t[k] for k in ("1", "4", "6", "24"))
    assert d["three_solutions_differ_by_transvection"] == 7776
    assert d["three_solution_symmetries"] == {"-M complementary": 7776, "M^-1 same": 7776, "JMJ same": 7776}
    p = d["minus_identity_pairing_by_cell"]
    assert p["fix"] == {"10": 648} and p["rev"] == {"01": 648}


def test_11352_three_qutrits():
    d = load("w33_pass11352_three_qutrit_geometry.json")
    c = d["uniform_combined"]
    assert c["decided"] == 31994 and c["bad"] == 4265 and c["z_vs_one_eighth"] > 4
    assert d["weyl_cross_check"]["agree"] == d["weyl_cross_check"]["frame_checks"] == 660
    assert d["targeted_fixed"]["Mz1 = z1"]["bad"] == 150
    assert d["targeted_reversed"]["Mz1 = -z1"]["good"] == 149


def test_11353_level_statistics():
    d = load("w33_pass11353_level_statistics.json")
    e = d["ensembles"]
    assert abs(e["n=5:violating"]["r_mean"] - 0.5996) < 0.01
    assert abs(e["n=5:reversible"]["r_mean"] - 0.5307) < 0.01
    assert e["n=5:clifford"]["degenerate_spacing_fraction"] > 0.3
    assert d["theta_squared_values"] == ["(1+0j)", "(1-0j)"] or set(d["theta_squared_values"]) <= {"(1+0j)", "(1-0j)"}


def test_11354_three_eighths():
    d = load("w33_pass11354_decay_rate.json")
    r = d["reps"]
    assert abs(r["(3,3)"]["rho"] - 0.375) < 1e-9 and abs(r["(5,2)"]["rho"] - 0.375) < 1e-9
    assert r["(1,1)"]["clifford_invariants"] == r["(1,1)"]["su3_invariants"] == 1
    assert all(r[k]["rho"] == 0.0 for k in ("(3,0)", "(2,2)", "(4,1)", "(6,0)"))


def test_11352_affine_failure():
    d = load("w33_pass11352_affine_failure.json")
    rows = d["non_affine_classes"]
    assert d["scanned_seeds"] == 32000 and len(rows) == 2
    assert all(r["good_frames"] == 567 and r["weyl_agree"] == r["weyl_checks"] == 12 for r in rows)
    assert all(r["good_set_is_union_of_those_hyperplanes"] for r in rows)
