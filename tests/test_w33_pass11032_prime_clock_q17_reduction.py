import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11032_prime_clock_q17_reduction.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11032_prime_clock_q17_reduction.json").read_text())
def test_replay(): assert P.payload()==C
def test_q17_exact():
    assert C["q17"]["minimum_glue_norm"]==16
    assert C["q17"]["normalized_exact_census"]["words_checked"]==24137568
def test_boundary(): assert C["tower_status"]["all_prime_proof"] is False
