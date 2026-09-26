import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "p11020",
    ROOT / "analysis" / "w33_pass11020_signed_clock_extension_split.py",
)
P = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(P)

CERT_PATH = ROOT / "data" / "w33_pass11020_signed_clock_extension_split.json"
CERT = json.loads(CERT_PATH.read_text(encoding="utf-8"))

def test_certificate_replays_exactly():
    assert P.payload() == CERT

def test_signed_extension_splits():
    e = CERT["extension"]
    assert e["support_order"] == 48
    assert e["sign_kernel_order"] == 64
    assert e["full_monomial_extension_order"] == 3072
    assert e["splits"] is True
    assert e["split_subgroup_order"] == 48
def test_explicit_split_preserves_generator_orders():
    w = CERT["explicit_split_witness"]
    assert w["support_generator_orders"] == w["lift_generator_orders"]
    assert w["support_generator_orders"] == [3, 2]
    assert len(w["lift_generator_sign_vectors"]) == 2
    assert all(len(v) == 27 for v in w["lift_generator_sign_vectors"])

def test_central_minus_i_blocks_naive_clock_subspace():
    f = CERT["clock_reduction_firewall"]
    assert f["signed_repairs_checked"] == 64
    assert f["fibrewise_constant_signed_repairs"] == 0
    assert f["split_section_minusI_order"] == 2

def test_only_identity_in_split_section_is_fibrewise_constant():
    f = CERT["clock_reduction_firewall"]
    assert f["split_section_elements_preserving_fibre_constant_module"] == 1
    assert CERT["checks"]["split_projection_bijective"] is True
    assert CERT["checks"]["split_section_preserves_signed_cubic"] is True
