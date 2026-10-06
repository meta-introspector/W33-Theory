import json
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"analysis"/"w33_pass11547_history_hamming_canonical_clock.py"
DATA=ROOT/"data"/"PART_W33_PASS11547_HISTORY_HAMMING_CANONICAL_CLOCK.json"


def test_pass11547_history_hamming_canonical_clock():
    subprocess.run([sys.executable,str(SCRIPT)],cwd=ROOT,check=True,timeout=60)
    c=json.loads(DATA.read_text())

    assert c["status"]=="PASS_EXACT_FINITE_MINKOWSKI_HAMMING_CONJUGACY"
    assert c["hamming_scheme"]["standard_name"]=="ternary Hamming association scheme H(3,3)"
    assert c["exact_relation_dictionary"]["q=+1 / paper timelike shell"]=="Hamming distance 1"
    assert c["exact_relation_dictionary"]["q=-1=2 / paper spacelike shell"]=="Hamming distance 2"
    assert c["exact_relation_dictionary"]["q=0 nonzero / paper null shell"]=="Hamming distance 3"

    assert c["automorphism_group"]["order"]==1296
    assert c["automorphism_group"]["linear_stabilizer_order"]==48
    assert c["automorphism_group"]["pass11540_group_recovered"] is True

    assert c["canonical_hamming_frame"]["projective_orthogonal_frame_count"]==4
    assert c["canonical_hamming_frame"]["unique_all_q1_frame"] is True

    assert c["clock_factorization"]["phase_law"]=="on Fourier Hamming weight w, U has phase omega^{-w}"
    assert c["clock_factorization"]["factorization"]=="in the canonical Hamming mereology, U=C tensor C tensor C"

    assert c["qutrit_dimension_uniqueness"]["positive_distance_residue_injective_by_dimension"]["3"] is True
    assert c["qutrit_dimension_uniqueness"]["positive_distance_residue_injective_by_dimension"]["4"] is False
