import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11044_minimal_latent_dimension_and_54_budget.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11044_minimal_latent_dimension_and_54_budget.json").read_text())
def test_replay(): assert P.payload()==C
def test_minimal():
    assert C["classification"]["latent_irrep_content"]["dimension"]==9
def test_budget():
    assert C["compiler_budget"]["maximal_common_K_dimension"]==27
    assert C["compiler_budget"]["forced_retyped_dimension"]==54
