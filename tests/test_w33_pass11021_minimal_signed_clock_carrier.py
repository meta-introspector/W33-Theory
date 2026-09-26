import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "p11021",
    ROOT / "analysis" / "w33_pass11021_minimal_signed_clock_carrier.py",
)
P = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(P)

CERT_PATH = ROOT / "data" / "w33_pass11021_minimal_signed_clock_carrier.json"
CERT = json.loads(CERT_PATH.read_text(encoding="utf-8"))


def test_certificate_replays_exactly():
    assert P.payload() == CERT


def test_minimal_signed_clock_carrier_is_24d():
    m = CERT["minimal_closure"]
    assert m["four_indicator_span_rank"] == 24
    assert m["augmentation_span_rank"] == 24
    assert m["center_span_rank"] == 3
def test_coordinate_split_and_fibre_refinement():
    d = CERT["coordinate_decomposition"]
    f = CERT["clock_fibres"]
    assert d["support_orbit_sizes"] == [1, 2, 8, 16]
    assert d["noncentral_dimension"] == 24
    assert d["central_dimension"] == 3
    assert f["intersection_with_8_and_16"] == [[2, 4]] * 4


def test_quadratic_phase_graph_is_unique():
    q = CERT["quadratic_phase_orbit"]
    assert q["orbit_size"] == 8
    assert q["unique_coefficients_mod3"] == [0, 1, 0, 1, 1, 0]
    assert q["equation"] == "c = a*b + a + b (mod 3)"
    assert q["complement_noncentral_orbit_size"] == 16


def test_character_diagnostics_are_exact():
    c = CERT["character_diagnostics"]
    assert c["noncentral24"]["invariant_dimension"] == 2
    assert c["noncentral24"]["commutant_dimension"] == 19
    assert c["center3"]["invariant_dimension"] == 2
    assert c["full27"]["invariant_dimension"] == 4
