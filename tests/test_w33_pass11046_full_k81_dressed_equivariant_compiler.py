import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11046_full_k81_dressed_equivariant_compiler.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11046_full_k81_dressed_equivariant_compiler.json").read_text())
def test_replay(): assert P.payload()==C
def test_full_rank():
    assert all(x["T81_rank"]==81 and x["T81_det_mod_p"] for x in C["split_prime_certificates"])
def test_54_structure():
    assert C["why_54_retypings"]=={
      "equivariant":27,"retyped_central":27,"retyped_abelian":27,"total_retyped":54,
      "explanation":C["why_54_retypings"]["explanation"]}
