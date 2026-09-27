import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "p11023",
    ROOT / "analysis" / "w33_pass11023_mckay_e7_clock_saturation.py",
)
P = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(P)

CERT_PATH = ROOT / "data" / "w33_pass11023_mckay_e7_clock_saturation.json"
CERT = json.loads(CERT_PATH.read_text(encoding="utf-8"))


def test_certificate_replays_exactly():
    assert P.payload() == CERT


def test_corrected_tensor_quiver_is_directed_not_e7_tree():
    assert len(CERT["spin2a_directed_edges"]) == 14
    assert len(CERT["spin2b_directed_edges"]) == 14
    assert len(CERT["underlying_undirected_edges"]) == 11
    assert CERT["checks"]["two_faithful_quivers_are_transposes"] is True
    assert CERT["checks"]["faithful_tensor_matrix_is_not_symmetric"] is True


def test_clock_tensor_saturation_survives_correction():
    ident = CERT["exact_identities"]["full"]
    assert ident.startswith("S tensor V24 = 3*(")
    assert CERT["checks"]["clock_tensor_spin2a_saturates_nonlinear"] is True
    assert CERT["checks"]["clock_tensor_spin2b_saturates_nonlinear"] is True


def test_central_halves_saturate_opposite_parity_types():
    assert CERT["exact_identities"]["plus"] == (
        "S tensor Vplus = 3*(spin_2a + spin_2b + spin_4)"
    )
    assert CERT["exact_identities"]["minus"] == (
        "S tensor Vminus = 3*(quotient_2d + standard3_twist + standard3)"
    )


def test_tensor_dimension_and_second_order_identities():
    assert CERT["checks"]["dimension_vector_tensor_eigenvalue_2"] is True
    assert CERT["checks"]["green_identity"] is True
    assert CERT["irrep_dimension_vector"] == [1, 1, 2, 2, 2, 3, 3, 4]
    assert CERT["clock_multiplicity_vector"] == [2, 1, 0, 0, 0, 1, 2, 3]
