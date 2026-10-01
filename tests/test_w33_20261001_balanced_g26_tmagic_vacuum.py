import importlib.util,json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
P=ROOT/"analysis/w33_20261001_balanced_g26_tmagic_vacuum.py"
S=importlib.util.spec_from_file_location("tmagic_vac",P)
M=importlib.util.module_from_spec(S)
assert S.loader is not None
S.loader.exec_module(M)

def test_balanced_tmagic_certificate_replay():
    got=M.payload()
    frozen=json.loads((ROOT/"data/w33_20261001_balanced_g26_tmagic_vacuum.json").read_text())
    assert json.loads(json.dumps(got,default=str))==frozen
    assert got["exact_local"]["hessian_eigenvalues"]==[30,30,54,54]
    assert got["finite_orbit"]["T_state_orbit_size"]==72
    assert got["finite_orbit"]["T_state_stabilizer_order"]==3
    assert got["finite_orbit"]["MUB_axis_counts"]==[18,18,18,18]
    assert got["exact_local"]["unit_basic_invariants"]==["0","0","-1/729"]
