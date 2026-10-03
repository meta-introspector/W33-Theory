"""Regression for Passes 11330-11334 (frozen certificates plus fast exact recomputations)."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def load(name):
    return json.loads((ROOT / "data" / name).read_text())


def test_11330_unconditional_census():
    d = load("w33_pass11330_orbit_census.json")
    assert all(d["generators_commute_with_T1"]) and d["orbits"] == 4110
    assert d["validation_direct_vs_orbit_agree"] == d["validation_samples"] == 1000
    assert d["violating_cliffords_exact"] == 385344 and d["probability_exact"] == "223/2430"
    assert d["class_frame_counts"] == {"0": 45360, "54": 5160, "72": 24, "81": 1296}
    assert d["classes_violating_the_affine_law"] == 0 and d["five_frame_screen_exact"]
    counts = np.load(ROOT / "data" / "w33_pass11330_class_counts.npy")
    assert counts.shape == (51840,) and int(counts.sum()) == 385344


def test_11330_conjugation_formula_fast():
    import w33_pass11252_exact_reversibility as R
    import w33_pass11330_orbit_census as O
    wl = R.Weyl(2)
    rng = np.random.default_rng(5)
    gens = O.k_generators(wl)
    assert all(O.check_generators(wl, gens, rng))
    N, d = gens[0]
    M = np.array([[1, 1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    a = np.array([2, 1, 0, 1])
    C = wl.W[wl.index(a)] @ R.weil(wl, M, rng)
    D = wl.W[wl.index(d)] @ R.weil(wl, N, rng)
    Mp = (N @ M @ R._inv_mod3(N)) % 3
    C2 = wl.W[wl.index((N @ a + d - Mp @ d) % 3)] @ R.weil(wl, Mp, rng)
    r = np.trace(C2.conj().T @ (D @ C @ D.conj().T)) / 9
    assert abs(abs(r) - 1) < 1e-9


def test_11331_independent_method():
    d = load("w33_pass11331_twisted_fixed_point.json")
    assert d["N_T_equals_stab_z1"] and d["stab_z1"] == 648
    assert d["violating_cliffords_by_twisted_method"] == 385344 and d["agrees_with_pass_11330_every_class"]
    s = d["symplectic_level"]
    assert s["distinct_symplectic_pairs"] == 1944 and s["all_frames_violate_iff_symplectically_unsolvable"]
    assert s["split_solvable_by_frame_count"]["unsolvable:81"] == 1296


def test_11332_limit():
    d = load("w33_pass11332_reversibility_limit.json")
    assert d["T9_is_identity"] and d["T9m_reversible"]
    assert d["haar_reversible"] == 0 and d["lie_algebra_dimension"] == 8


def test_11332_t9_fast():
    import w33_pass11213_cubic_t_violation as P1
    assert np.allclose(np.linalg.matrix_power(P1.T1, 9), np.eye(3))


def test_11331_w33_geometry():
    d = load("w33_pass11331_twisted_fixed_point.json")["w33_geometry_of_bad_set"]
    assert d["Mz1 = z1"] == {"bad": 648, "all81": 0} and d["Mz1 = -z1"]["good"] == 648
    assert d["Mz1, M^-1 z1 collinear with z1 on different lines"] == {"good": 11664, "all81": 0}
    assert d["same line, other"] == {"good": 1296, "all81": 1296, "bad": 1296}
    assert d["Mz1 not collinear with z1"]["bad"] == 3888 and d["Mz1 not collinear with z1"]["good"] == 31104
    assert sum(v.get("bad", 0) for v in d.values()) == 6480


def test_11333_three_qutrits():
    d = load("w33_pass11333_three_qutrit_one_gate.json")
    assert d["flag_classes"] == 1100 and d["flagged_bad"] == 122
    assert abs(d["one_eighth_z_score"]) < 2
    assert all(c["nonviolating_affine"] and c["violating_frames"] == 486 for c in d["full_frame_checks_of_bad_classes"])


def test_11334_f9_law():
    d = load("w33_pass11334_f9_three_leg_law.json")
    assert d["f9_perfect_classes"] == 64 and d["three_leg_law_exhaustive"]
    assert d["magic_legs"]["3"] == {"violating": 165888, "cases": 165888, "fraction": 1.0}
    assert d["magic_legs"]["1"]["violating"] * 12 == d["magic_legs"]["1"]["cases"]
    assert d["covers_unreduced_cases"] == 13436928
