import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"


def load(name):
    return json.loads((DATA / name).read_text())


def test_11260_semisimple_cartan():
    d = load("w33_pass11260_exact_semisimple_cartan.json")
    assert d["cartan_checks"]["dimension"] == 3
    assert d["cartan_checks"]["pairwise_brackets_zero"]
    assert d["background"]["cubic_jacobian_rank"] == 78
    assert d["scaled_adjoint"]["rank_N"] == d["scaled_adjoint"]["rank_N_cubed"] == 234
    assert [d["N_cubed_blocks"][g]["rank"] for g in ("g0", "g1", "g2")] == [78, 78, 78]


def test_11261_mirror_multiplicities():
    d = load("w33_pass11261_g26_cartan_coordinate_map.json")
    assert d["restricted_pfaffian"]["formula"] == "2^17*C3*N9*S9^3"
    assert d["restricted_pfaffian"]["unisolvent_grid_points"] == 820
    m = d["mirror_divisor"]
    assert m["9_qutrit_SIC_mirrors"]["multiplicity"] == 3
    assert m["12_qutrit_stabilizer_mirrors"]["multiplicity"] == 1
def test_11262_antiunitary_completion():
    d = load("w33_pass11262_g26_mub_s4_yukawa_bridge.json")
    m = d["mub_action"]
    assert sorted(map(len, m["components"])) == [3, 3, 3, 3]
    assert (m["complex_linear_image"], m["complex_linear_image_order"]) == ("A4", 12)
    assert m["all_linear_permutations_even"]
    assert m["conjugation_is_odd"]
    assert (m["completed_image"], m["completed_image_order"]) == ("S4", 24)


def test_11263_z3_spectral_triplets():
    d = load("w33_pass11263_full_graded_spectrum.json")
    h = d["declared_hessian"]
    assert (h["full_dimension"], h["full_rank"], h["full_zero_modes"]) == (248, 234, 14)
    assert (h["matter81_rank"], h["matter81_Cartan_zero_modes"]) == (78, 3)
    s = d["basis_independent_N_cubed_spectrum"]
    assert s["nonzero_modes_per_grade"] == 78
    assert s["zero_multiplicities_g0_g1_g2"] == [8, 3, 3]
    assert [r["Z3_character_trace"] for r in s["exact_Z3_character_traces_powers_1_to_12"]] == ["0"] * 12
    assert not d["ordinary_Z2_sign_audit"]["any_assignment_cancels_moments_0_to_6"]
