import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11036_triqutrit_cubic_correlation_descent.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11036_triqutrit_cubic_correlation_descent.json").read_text())
def test_replay(): assert P.payload()==C
def test_ladder(): assert C["cubic_echo"]["successive_exponents"]==["abc","bc","c","1"]
def test_rank_firewall():
    assert C["slant_rank_audit"]["individual_pair_block_ranks"]==[4,4,4]
    assert C["slant_rank_audit"]["joint_fixed_three_slants_rank"]==7
