import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11048_induced_h27_hesse_latent_modules.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11048_induced_h27_hesse_latent_modules.json").read_text())
def test_replay(): assert P.payload()==C
def test_twelve(): assert len(C["twelve_sectors"])==12
def test_regular_sums(): assert all(x["direct_sum"]=="Reg(H27)" for x in C["four_directions"].values())
