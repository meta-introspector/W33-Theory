import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "p11022",
    ROOT / "analysis" / "w33_pass11022_binary_octahedral_clock_decomposition.py",
)
P = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(P)

CERT_PATH = ROOT / "data" / "w33_pass11022_binary_octahedral_clock_decomposition.json"
CERT = json.loads(CERT_PATH.read_text(encoding="utf-8"))


def test_certificate_replays_exactly():
    assert P.payload() == CERT


def test_exact_24d_character_and_decomposition():
    assert CERT["character_24"] == [24, 0, 0, 6, 0, 0, 0, 2]
    d = CERT["decomposition_24"]
    assert d["trivial"] == 2
    assert d["det"] == 1
    assert d["standard3_twist"] == 1
    assert d["standard3"] == 2
    assert d["spin_4"] == 3
def test_central_minus_i_splits_12_plus_12():
    z = CERT["central_minus_I"]
    assert z["plus_dimension"] == 12
    assert z["minus_dimension"] == 12
    assert z["per_clock_fibre_plus_minus"] == [[3, 3]] * 4
    assert set(z["projected_orbit_ranks"].values()) == {12}


def test_odd_half_is_exactly_three_spin4():
    z = CERT["central_minus_I"]
    minus = z["minus_decomposition"]
    assert minus["spin_4"] == 3
    assert sum(v for k, v in minus.items() if k != "spin_4") == 0


def test_center_three_is_two_trivial_plus_det():
    d = CERT["center3"]["decomposition"]
    assert d["trivial"] == 2
    assert d["det"] == 1
    assert sum(v for k, v in d.items() if k not in {"trivial", "det"}) == 0


def test_projective_standard_is_identified_internally():
    assert CERT["checks"]["P1_standard_character_identified"] is True
    assert CERT["checks"]["augmentation_generates_both_12s"] is True
