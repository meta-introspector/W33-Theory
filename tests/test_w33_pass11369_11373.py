"""Regression for Passes 11369-11373 (J6 incomplete on PU(3); twisted-indicator degree law; flavour vs arrow;
irrep contraction rates; exact three-qutrit one-gate fraction)."""

import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def load(name):
    return json.load(open(ROOT / "data" / name))


def test_11369_j6_incomplete_j8_clean():
    d = load("w33_pass11369_j6_completeness.json")
    g6, g8 = d["global_d3_J6"], d["global_d3_J8"]
    certified6 = g6["spurious_zeros"] - g6["spurious_not_certified_by_next_witness"]
    assert certified6 > 100                                   # J6 has many certified non-reversible zeros
    assert g8["spurious_zeros"] == g8["spurious_not_certified_by_next_witness"]   # none certified for J8
    assert g8["max_distance_at_zeros"] < 0.01
    c = d["control_d5_J4"]
    assert c["spurious_zeros"] - c["spurious_not_certified_by_next_witness"] > 0   # the search is not blind
    gen = [s for s in d["local_d3_J6"] if s.get("stratum_dim") == 5]
    assert len(gen) == 1 and gen[0]["rank_equals_codim"]
    deep = d["deep_words"]
    assert sum(v["words"] for v in deep.values()) > 1_390_000
    assert all(v["j6_zero_but_violating"] == 0 and v["j6_positive_but_reversible"] == 0 for v in deep.values())
    r = d["refined_spurious_zero"]
    assert r["residual_norms"][-1] < 1e-12 and r["dF_rank"] == 4 and r["J8"] > 0.1 and r["distance_to_reversible"] > 0.1


def test_11369_spurious_point_fast():
    import w33_pass11357_jarlskog_degree_by_dimension as P7
    import w33_pass11369_j6_completeness as M
    Cl = P7.clifford_group(3)
    U, J = M.descend(M.haar(3, np.random.default_rng(2)), 3, Cl, M.hermitian_basis(3))
    assert abs(J) < 1e-11 and P7.J(U, 4, Cl) > 0.1 and M.rev_distance(U, Cl) > 0.1


def test_11370_degree_law():
    d = load("w33_pass11370_twisted_indicator_degree_law.json")
    law = d["degree_law"]
    assert {k: v["witness_degree"] for k, v in law.items()} == {"2": 10, "3": 6, "5": 4, "7": 4}
    assert d["matches_pass_11357"] and d["d7_prediction_confirmed"]
    assert law["3"]["cumulative"]["3"]["odd"] == 7 and law["2"]["cumulative"]["5"]["odd"] == 1
    assert d["octahedral_spin_restrictions"]["5"]["mult"] == {"E": 1, "T1": 2, "T2": 1}
    assert all(d["octahedral_spin_restrictions"][str(j)]["odd"] == 0 for j in range(5))


def test_11370_counts_fast():
    import w33_pass11357_jarlskog_degree_by_dimension as P7
    import w33_pass11370_twisted_indicator_degree_law as T
    rows = T.counts(P7.clifford_group(3), 3)
    assert sum(r["odd"] for r in rows) == 7 and sum(r["odd"] for r in rows if r["k"] < 3) == 0


def test_11371_flavour_vs_arrow():
    d = load("w33_pass11371_flavour_vs_arrow.json")
    assert d["det_commutator_identity_max_dev"] < 1e-12 and d["theorem1_max_dev"] < 1e-12
    assert d["borel_cliffords"] == 54 and d["nonmonomial_cliffords_all_maximal_J"]
    w = d["words"]
    assert w["1"]["arrow_by_flavour_cp"] == {"0": {"reversible": 36, "violating": 18}, "max": {"reversible": 162},
                                             "other": {}}
    assert w["3"]["arrow_by_flavour_cp"]["max"] == {"reversible": 22680, "violating": 648}
    assert d["only_subgroup_types"] and all("K" not in t or t == "TKI" for v in w.values() for t in v["types_seen"])


def test_11371_altland_zirnbauer():
    az = load("w33_pass11371_flavour_vs_arrow.json")["altland_zirnbauer"]
    assert az["1"]["az_classes"] == {"BDI": 18, "AI": 180, "AIII": 18}
    assert az["2"]["az_classes"] == {"BDI": 270, "AI": 2808, "AIII": 54, "A": 2052}
    for v in az.values():
        assert "D" not in v["az_classes"] and v["chiral_pairing_failures"] == 0
        assert not any(k.endswith("= -1") for k in v["squares"])         # odd dimension: T^2 = C^2 = +1
        assert v["reversible_words"] == v["reversible_words_with_a_reversal_squaring_to_plus_one"]


def test_11371_two_qutrit_violators_not_all_chiral():
    c = load("w33_pass11371_flavour_vs_arrow.json")["two_qutrit_one_gate_chirality"]
    assert c["violating, unpaired"] > 100 and c["violating, paired"] > 100
    assert c["reversible, paired"] / (c["reversible, paired"] + c["reversible, unpaired"]) < 0.1


def test_11371_transpose_even_fast():
    import w33_pass11371_flavour_vs_arrow as F
    rng = np.random.default_rng(3)
    U, _ = np.linalg.qr(rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)))
    assert abs(F.jcp(U.T) - F.jcp(U)) < 1e-14 and abs(F.jcp(U.conj()) + F.jcp(U)) < 1e-14


def test_11372_rates():
    d = load("w33_pass11372_irrep_contraction_rates.json")
    ir = d["irreps"]
    assert d["lie_algebra_dim_of_closure_of_<Cl,T Cl T^-1>"] == 8
    assert abs(ir["(3,3)"]["rho"] - 3 / 8) < 1e-9 and abs(ir["(5,5)"]["rho"] - 11 / 24) < 1e-9
    assert abs(ir["(9,0)"]["rho"] - 1 / 2) < 1e-9 and abs(ir["(12,0)"]["rho"] - 7 / 12) < 1e-9
    assert abs(ir["(7,7)"]["rho"] - (1 / 8 + np.sqrt(7 / 32))) < 1e-9
    assert all(v["rho"] < 1 for v in ir.values())


def test_11372_moments_family_designs():
    d = load("w33_pass11372_irrep_contraction_rates.json")
    m = d["trace_moment_law"]
    assert m["max_dev"] < 1e-9 and m["t=3"]["nonzero_spectrum_multiplicities"] == {"1.0000000000": 6, "0.3750000000": 1}
    assert m["t=4"]["nonzero_spectrum_multiplicities"] == {"1.0000000000": 23, "0.3750000000": 16}
    cp = d["exact_characteristic_polynomials"]
    assert cp["(18,0)"]["integer_polynomial"] == [864, -984, 230, 9]
    assert cp["(8,8)"]["cubic"] == [2592, -1332, -585, 140]
    assert len(cp["(7,7)"]["eigenvalues_solving_it"]) == 2
    r = {int(k): v["rate"] for k, v in d["design_rates"].items()}
    assert r[1] == r[2] == 0 and abs(r[3] - 3 / 8) < 1e-9 and abs(r[5] - 11 / 24) < 1e-9 and abs(r[6] - 7 / 12) < 1e-9
    assert abs(r[8] - np.roots([2592, -1332, -585, 140]).real.max()) < 1e-9
    assert abs(d["p0_family"]["(18,0)"]["moduli"][0] - np.roots([864, -984, 230, 9]).real.max()) < 1e-7
    x = d["independent_coherent_state_crosscheck"]
    assert abs(x["(3,3)"][1] - 3 / 8) < 1e-7 and abs(x["(5,5)"][1] - 11 / 24) < 1e-7


def test_11372_moment_law_fast():
    import w33_pass11213_cubic_t_violation as P1
    CF = np.array(P1.clifford1())
    U1 = np.einsum('cij,jk->cik', CF, P1.T1)
    assert abs(np.mean(np.abs(np.einsum('cii->c', U1)) ** 6) - (6 + 3 / 8)) < 1e-12


def test_11373_exact_three_qutrit_fraction():
    d = load("w33_pass11373_three_qutrit_exact_fraction.json")
    assert d["n2_matches_223_2430"] and d["n2_bad_classes_is_one_eighth"]
    r = d["n=3"]
    assert r["orbits"] == 2308 and r["sum_of_orbit_sizes"] == 9170703360
    assert r["violating_fraction"] == "240857/2388204" and r["bad_class_fraction"] == "437/3276"
    assert r["invariance_controls"] == 2 * r["slow_classes_by_frame_orbits"] == 218
    c = r["cells"]
    assert c["Mz1 = z1"]["bad_class_fraction"] == "1" and c["Mz1 = -z1"]["bad_class_fraction"] == "0"
    assert c["same line, M^2 z1 = z1"]["bad_class_fraction"] == "1"
    assert c["same line, M^2 z1 = -z1"]["bad_class_fraction"] == "0"
    assert c["collinear, different lines"]["bad_class_fraction"] == "1/9"
    assert c["non-collinear"]["bad_class_fraction"] == "10/81"
    assert c["same line, other"]["bad_class_fraction"] == "77/342"
    total = sum(v["classes"] for v in c.values())
    assert Fraction(sum(v["bad_frames"] for v in c.values()), total * 729) == Fraction(240857, 2388204)


def test_11373_orbit_files_fast():
    lines = (ROOT / "data" / "w33_pass11373_orbits_n2.txt").read_text().splitlines()
    assert len(lines) == 200 and sum(int(x.rsplit(";", 1)[1]) for x in lines) == 51840


def test_11373_reversal_route():
    r = load("w33_pass11373_three_qutrit_exact_fraction.json")["reversal_route_n2"]
    assert r["classes"] == {"sigma exists: True, best single-sigma frames: 27": 27,
                            "sigma exists: True, best single-sigma frames: 81": 621}
