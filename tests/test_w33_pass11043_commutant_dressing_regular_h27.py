import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11043_commutant_dressing_regular_h27.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11043_commutant_dressing_regular_h27.json").read_text())
def test_replay(): assert P.payload()==C
def test_regularized():
    assert C["before_after"]["dressed_E6_27"]=={"V_omega":3,"V_omega2":3,"one_dimensional_total":9}
def test_character():
    assert all(r["chi_V_tensor_A9"]==r["chi_regular"] for r in C["character_table_rows"])
