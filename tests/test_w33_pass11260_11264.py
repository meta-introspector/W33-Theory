"""Regression locks for the exact Cartan-to-mu packet, Passes 11260--11264."""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"


def load(name):
    return json.loads((DATA / name).read_text())


def test_pass11260_is_an_exact_semisimple_cartan_plane():
    d = load("w33_pass11260_exact_semisimple_cartan.json")
    assert d["cartan_checks"] == {
        "dimension": 3,
        "pairwise_brackets_zero": True,
        "background_equals_minus_basis_2": True,
        "pivot_rows": [3, 4, 6],
        "coordinate_minor": 4,
    }
    assert d["background"]["cubic_jacobian_rank"] == 78
    assert d["scaled_adjoint"]["rank_N"] == 234
    assert d["scaled_adjoint"]["rank_N_cubed"] == 234
    assert [d["N_cubed_blocks"][g]["rank"] for g in ("g0", "g1", "g2")] == [78] * 3
    assert all(
        set(d["repeated_degree9_factor_nullities"][g].values()) == {27}
        for g in ("g0", "g1", "g2")
    )


def test_pass11261_map_and_relative_mirror_multiplicities():
    d = load("w33_pass11261_g26_cartan_coordinate_map.json")
    assert d["restricted_pfaffian"]["formula"] == "2^17*C3*N9*S9^3"
    assert d["restricted_pfaffian"]["factor_degrees_and_multiplicities"] == [
        [3, 1], [9, 1], [9, 3]
    ]
    assert d["restricted_pfaffian"]["unisolvent_grid_points"] == 820
    assert d["exact_map_identities"] == {
        "XYZ": "C3",
        "standard_u9": "3*(omega-omega^2)*S9",
        "product_12_stabilizer_mirrors": "-27*C3*N9",
    }
    mirrors = d["mirror_divisor"]
    assert mirrors["12_qutrit_stabilizer_mirrors"]["multiplicity"] == 1
    assert mirrors["9_qutrit_SIC_mirrors"]["multiplicity"] == 3


def test_pass11262_s4_is_the_antiunitary_completion():
    d = load("w33_pass11262_g26_mub_s4_yukawa_bridge.json")
    m = d["mub_action"]
    assert sorted(map(len, m["components"])) == [3, 3, 3, 3]
    assert (m["complex_linear_image"], m["complex_linear_image_order"]) == ("A4", 12)
    assert m["all_linear_permutations_even"]
    assert m["conjugation_is_odd"]
    assert (m["completed_image"], m["completed_image_order"]) == ("S4", 24)
    assert len(d["flavon_dictionary"]["twelve_oriented_chords"]) == 12
    assert d["yukawa_audit"]["epsilon_does_not_determine_angles"]


def test_pass11263_complete_spectrum_and_z3_triplets():
    d = load("w33_pass11263_full_graded_spectrum.json")
    h = d["declared_hessian"]
    assert (h["full_dimension"], h["full_rank"], h["full_zero_modes"]) == (248, 234, 14)
    assert (h["matter81_rank"], h["matter81_Cartan_zero_modes"]) == (78, 3)
    blocks = h["blocks"]
    assert [(v["dimension"], v["rank"], v["zero_modes"]) for v in blocks.values()] == [
        (86, 78, 8), (81, 78, 3), (81, 78, 3)
    ]
    s = d["basis_independent_N_cubed_spectrum"]
    assert s["nonzero_modes_per_grade"] == 78
    assert s["zero_multiplicities_g0_g1_g2"] == [8, 3, 3]
    traces = s["exact_Z3_character_traces_powers_1_to_12"]
    assert [r["power"] for r in traces] == list(range(1, 13))
    assert [r["Z3_character_trace"] for r in traces] == ["0"] * 12
    assert not d["ordinary_Z2_sign_audit"]["any_assignment_cancels_moments_0_to_6"]


def test_pass11264_raw_rcharge_archive_and_exhaustive_audit():
    path = DATA / "w33_pass11264_z6i_r_charges.json.gz"
    compressed = path.read_bytes()
    plain = gzip.decompress(compressed)
    raw = json.loads(plain)
    assert len(raw) == 87
    assert sum(map(len, raw.values())) == 30980
    assert hashlib.sha256(plain).hexdigest() == (
        "2ca5326ef217282a8038900f296153efabe99c480abe0769149c42c661001f61"
    )
    d = load("w33_20261001_z6i_rcharge_lattice_quotient.json")
    assert d["extraction"]["uncompressed_sha256"] == hashlib.sha256(plain).hexdigest()
    assert d["audit"]["models_total"] == 87
    assert d["audit"]["models_dflat"] == 23
    assert d["audit"]["vacua"] == 6695116
    assert d["audit"]["mu_protected_vacua"] == 6695116
    assert d["audit"]["clean_vacua"] == 0
    assert d["audit"]["unique_support_lattices"] == 1464
    assert d["audit"]["exhaustive"]
