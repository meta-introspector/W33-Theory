import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11084_pgsp_u81_four_sheet_cover.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11084_pgsp_u81_four_sheet_cover.json").read_text())
def test_replay(): assert P.payload()==C
def test_cover(): assert C["orders"]["PGSp_over_U81"]==640 and C["deck_generators"]["quotient"]=="V4"
