import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "analysis/w33_pass11218_e8_cubic_cartan_kernel.py"


def load_module():
    spec = importlib.util.spec_from_file_location("pass11218", SOURCE)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_pass11218_exact_abelian_centralizer_and_generic_rank():
    out = load_module().main(write=False)
    assert out["status"].startswith("PASS_GENERIC_CUBIC_JACOBIAN_RANK78")
    assert out["explicit_regular_background"]["exact_rank"] == 78
    assert out["abelian_centralizer"]["dimension"] == 3
    assert out["abelian_centralizer"]["basis_pairwise_brackets_zero"] is True
    assert out["abelian_centralizer"]["is_cartan_subspace"] is False
    assert out["pass11255_correction"]["reason"] == "the first stored kernel vector has ad-nilpotency index five"
    assert out["rank_witness"]["principal_minor_size"] == 78
    assert out["rank_witness"]["determinant_factorization"] == {
        "2": 76,
        "3": 24,
        "5": 2,
        "7": 2,
    }
    assert out["vinberg_classification"]["little_weyl_group"] == "Shephard-Todd G26"
    assert out["vinberg_classification"]["invariant_degrees"] == [6, 12, 18]
    assert out["g26_firewall"]["compatible_with_prior_no_go"] is True
