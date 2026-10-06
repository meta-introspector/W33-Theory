import json
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"analysis"/"w33_pass11541_history_quadratic_shell_scheme.py"
DATA=ROOT/"data"/"PART_W33_PASS11541_HISTORY_QUADRATIC_SHELL_SCHEME.json"

def test_pass11541_history_quadratic_shell_scheme():
    subprocess.run([sys.executable,str(SCRIPT)],cwd=ROOT,check=True,timeout=60)
    c=json.loads(DATA.read_text())
    assert c["status"]=="PASS_EXACT_HISTORY_QUADRATIC_SHELL_ASSOCIATION_SCHEME"
    assert c["valencies"]==[1,8,6,12]
    assert c["primitive_multiplicities"]==[1,8,6,12]
    assert c["first_eigenmatrix_P"]==[
        [1,8,6,12],
        [1,-1,-3,3],
        [1,-4,3,0],
        [1,2,0,-3],
    ]
    assert c["P_squared"]==[
        [27,0,0,0],
        [0,27,0,0],
        [0,0,27,0],
        [0,0,0,27],
    ]
