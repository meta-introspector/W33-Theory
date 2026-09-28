"""Regression for Pass 11097: fractional fermions of the SO(16)xSO(16) A8 models -- massive in 33/104 (all rules)."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11097_fractional_fermion_masses as P  # noqa: E402

C = json.loads(P.OUT.read_text())


def test_symmetry_test_counts():
    s = C["symmetry_test"]
    assert s[P.KEY]["models_all_fractional_massive"] == 33 and s[P.KEY]["models_no_light_hidden_singlet"] == 42
    assert s["B_hidden_broken|gauge+Z2W+PG+SG"]["models_all_fractional_massive"] == 66


def test_bottleneck_orders():
    assert C["bottleneck_order_T_star"] == {"1": 6, "3": 15, "4": 6, "5": 4, "7": 2}


def test_sample_models_recomputed():
    US, TH = P.sample_models()
    frozen = json.loads(P.FRAC.read_text())
    for i in US:
        st = P.symmetry_test(US[i], TH[i])
        assert all(st["results"][k]["light_fractional_states"] == frozen[str(i)]["results"][k]["light_fractional_states"]
                   for k in st["results"])
    assert P.cubic_T1(US[12], TH[12]) == (0, 0)          # model 12: every fractional fermion massive at cubic order
