import importlib.util
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load(stem):
    path = ROOT / "analysis" / f"{stem}.py"
    spec = importlib.util.spec_from_file_location(f"test_{stem}", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_pass11255_exact_nilpotence_and_stored_polynomial_identity():
    out = load("w33_pass11255_restricted_pfaffian_cartan_audit").main(False, False)
    assert out["semisimplicity_audit"]["nilpotency_index"] == 5
    assert out["checks"]["ad_fourth_power_nonzero"] is True
    assert out["checks"]["ad_fifth_power_zero"] is True
    stored = json.loads((ROOT / "data/w33_pass11255_restricted_pfaffian_cartan_audit.json").read_text())
    assert stored["restricted_pfaffian"]["exact_unisolvent_grid_points"] == 820
    assert stored["checks"]["degree39_factorization_exact"] is True


def test_pass11256_standard_g26_invariants_and_hessian():
    out = load("w33_pass11256_g26_invariant_potential").main(False)
    assert out["basic_invariants"]["degrees"] == [6, 12, 18]
    assert out["reflection_jacobian"]["degree"] == 33
    assert out["reflection_jacobian"]["complex_reflection_hyperplanes_counted_with_distinct_linear_factors"] == 21
    assert out["checks"]["hessian_is_positive_definite"] is True
    assert out["firewall"]["pass11218_coordinates_used"] is False


def test_pass11257_all_schlaefli_references_are_audited():
    out = load("w33_pass11257_trinification_hessian_firewall").main(False)
    audit = out["schlaefli_1_10_16_audit"]
    assert audit["references_checked"] == 27
    assert audit["distinct_block_trace_profiles"] == 27
    assert out["scope"]["mass_ratio_claim"] is False


def test_pass11258_rephasing_invariant_conditional_cp():
    out = load("w33_pass11258_g26_cubic_cp_interface").main(False)
    assert abs(out["result"]["jarlskog"]) > 0.009
    assert out["checks"]["rephasing_invariant"] is True
    assert out["checks"]["complex_conjugation_flips_sign"] is True
    assert out["interpretation"]["physical_ckm_claim"] is False


def test_pass11259_discriminant_and_supertrace_boundaries():
    out = load("w33_pass11259_g26_discriminant_supertrace").main(False)
    assert out["checks"]["singular_strata_have_zero_modes"] is True
    assert out["checks"]["reflection_cp_counterexample_nonzero"] is True
    assert out["checks"]["all_six_supertraces_nonzero"] is True
    assert out["scope"]["cosmological_constant_claim"] is False
