import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11070_dual_parabolic_flag_diamond.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11070_dual_parabolic_flag_diamond.json").read_text())
def test_replay(): assert P.payload()==C
def test_diamond(): assert C["radical_weld"]["generated_order"]==81 and C["radical_weld"]["intersection_order"]==9
def test_flag(): assert C["flag_borel"]["order"]==162 and C["flag_borel"]["sylow3_order"]==81
