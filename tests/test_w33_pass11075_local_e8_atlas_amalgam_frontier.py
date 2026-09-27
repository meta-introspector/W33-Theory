import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11075_local_e8_atlas_amalgam_frontier.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11075_local_e8_atlas_amalgam_frontier.json").read_text())
def test_replay(): assert P.payload()==C
def test_local_e8(): assert C["local_chart"]["E8_root_count"]==240 and C["local_chart"]["reflection_closed"]
def test_frontier(): assert C["atlas"]["overlap_H1"]=="Z^81" and C["checks"]["global_amalgam_not_silently_asserted"]
