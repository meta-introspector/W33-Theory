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


def test_mckay_graph_is_affine_e7_shape():
    assert len(CERT["spin2a_mckay_edges"]) == 7
    assert len(CERT["spin2b_mckay_edges"]) == 7
    assert CERT["affine_dimension_vector"] == [1, 1, 2, 2, 2, 3, 3, 4]
def test_clock_spinor_tensor_saturates_all_nonlinear_nodes():
    ident = CERT["exact_identities"]["full"]
    assert ident.startswith("S tensor V24 = 3*(")
    assert CERT["checks"]["clock_tensor_spin2a_saturates_nonlinear"] is True
    assert CERT["checks"]["clock_tensor_spin2b_saturates_nonlinear"] is True


def test_central_halves_saturate_opposite_bipartitions():
    assert CERT["exact_identities"]["plus"] == (
        "S tensor Vplus = 3*(spin_2a + spin_2b + spin_4)"
    )
    assert CERT["exact_identities"]["minus"] == (
        "S tensor Vminus = 3*(quotient_2d + standard3_twist + standard3)"
    )


def test_graph_vector_identities():
    assert CERT["checks"]["dimension_vector_is_affine_null_mark"] is True
    assert CERT["checks"]["green_identity"] is True
    assert CERT["clock_multiplicity_vector"] == [2, 1, 0, 0, 0, 1, 2, 3]
