import importlib.util, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("p11029", ROOT/"analysis"/"w33_pass11029_exact_twirl_dynamics_firewall.py")
P = importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(P)
CERT = json.loads((ROOT/"data"/"w33_pass11029_exact_twirl_dynamics_firewall.json").read_text())

def test_certificate_replays():
    assert P.payload() == CERT

def test_canonical_twirl_fixed_dimensions():
    assert CERT["carrier"]["vector_fixed_dimension"] == 2
    assert CERT["carrier"]["operator_commutant_dimension"] == 19

def test_clock_is_not_exact_invariant_image():
    f = CERT["clock_firewall"]
    assert f["coarse_augmentation_dimension"] == 3
    assert f["augmentation_plus_central_image_rank"] == 6
    assert f["full_signed_orbit_span_rank"] == 24

def test_channel_semigroup_is_frozen():
    c = CERT["canonical_channel"]
    assert "completely_positive" in c["properties"]
    assert c["lindbladian"] == "L=T-I"
