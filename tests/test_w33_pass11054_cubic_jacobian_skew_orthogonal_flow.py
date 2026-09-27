import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11054_cubic_jacobian_skew_orthogonal_flow.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11054_cubic_jacobian_skew_orthogonal_flow.json").read_text())
def test_replay(): assert P.payload()==C
def test_basis(): assert C["basis_jacobians"]["independent_span_dimension"]==81 and C["basis_jacobians"]["exact_rank_each"]==20
def test_flow(): assert C["checks"]["orthogonal_exponential_follows_exactly"]
