import importlib.util,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
P=ROOT/"analysis/w33_20261001_u81_cocycle_blind_54_transducer.py"
S=importlib.util.spec_from_file_location("u81_54",P)
M=importlib.util.module_from_spec(S);S.loader.exec_module(M)

def test_certificate_replay():
    got=M.payload()
    frozen=json.loads((ROOT/"data/w33_20261001_u81_cocycle_blind_54_transducer.json").read_text())
    assert got==frozen
    assert got["orientations"]["plus"]["cocycle_blind_factors"]==[6,7,8,9]
    assert got["orientations"]["minus"]["cocycle_blind_factors"]==[6,7,8,9]
    assert [x["factor"] for x in got["orientations"]["plus"]["perfect_survivors"]]==[7,8]
    assert [x["factor"] for x in got["orientations"]["minus"]["perfect_survivors"]]==[6,9]
    assert all(got["checks"].values())
