import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11040_exceptional_cubic_temporal_echo.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11040_exceptional_cubic_temporal_echo.json").read_text())
def test_replay(): assert P.payload()==C
def test_basis_cases():
    assert C["exceptional_level"]["basis_cases"]==360
    assert C["exceptional_level"]["nonzero_cases"]==180
def test_orientation(): assert C["gate_level"]["orientation_from_difference_order"] is False
