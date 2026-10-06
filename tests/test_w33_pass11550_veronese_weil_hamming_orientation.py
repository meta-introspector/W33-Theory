import json
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"analysis"/"w33_pass11550_veronese_weil_hamming_orientation.py"
DATA=ROOT/"data"/"PART_W33_PASS11550_VERONESE_WEIL_HAMMING_ORIENTATION.json"


def test_pass11550_veronese_weil_hamming_orientation():
    subprocess.run([sys.executable,str(SCRIPT)],cwd=ROOT,check=True,timeout=60)
    c=json.loads(DATA.read_text())
    assert c["status"]=="PASS_VERONESE_WEIL_HAMMING_ORIENTATION_WELD"
    assert c["veronese_sheet"]["chi_value"]==-1
    assert c["veronese_sheet"]["old_weil_phase"]=="-i"
    assert c["negative_veronese_sheet"]["chi_value"]==1
    assert c["negative_veronese_sheet"]["old_weil_phase"]=="+i"
    assert c["closed_formula"]["Weil_phase_equals"]=="i * chi"
    assert c["closed_formula"]["objectwise_matches"]==8
    assert c["cubic_invariant"]["PSp_law"]=="p(Mz)=p(z)"
    assert c["cubic_invariant"]["full_O_law"]=="p(Mz)=sigma(M)*p(z)"
    assert c["tetracode_orientation_firewall"]["veronese_kills_projective_sign"] is True
