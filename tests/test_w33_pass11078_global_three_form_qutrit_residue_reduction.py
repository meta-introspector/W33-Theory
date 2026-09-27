import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11078_global_three_form_qutrit_residue_reduction.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11078_global_three_form_qutrit_residue_reduction.json").read_text())
def test_replay(): assert P.payload()==C
def test_count(): assert C["census"]["exhaustive_contraction_identities"]==787320
def test_residue(): assert C["checks"]["all40_residue_forms_rank2"] and C["census"]["projective_directions_per_residue"]==4
