"""Regression for Pass 11224: AME(10,3) = Glynn's code; |Aut_LC| = 2880; one orbit."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def test_one_state_recomputed():
    import w33_pass11175_scan_five_qutrit as F
    import w33_pass11190_ame10_sign_structure as P
    import w33_pass11224_ame10_f9_structure as A
    d = json.loads(P.TABU.read_text())
    G = F.to_mat(np.array(d["graphs"][0]))
    st = A.f9_structures(G)
    assert len(st) == 4
    tr1 = [s for s in st if all(int(np.trace(A.ORD8[k]) % 3) == 1 for k in s)]
    rows, _ = A.to_f9_code(G, tr1[0])
    assert A.schur_square_dim(rows) == (5, 10)


def test_frozen():
    d = json.loads((ROOT / "data" / "w33_pass11224_ame10_f9_structure.json").read_text())
    assert d["all_f9_linear"] and not d["all_grs"] and len(d["states"]) == 71
    assert all(s["schur_square_dim"] == 10 and s["f9_dimension"] == 5 for s in d["states"])
    g = d["glynn"]
    assert g["grs_hermitian_weightings"] == 0 and g["glynn_hermitian_weightings"] == 2
    assert g["aut_lc_order"] == 2880 and g["permutation_group_order"] == 720 and g["sharply_3_transitive"]
    assert g["single_orbit"] and g["orbit_size"] == 79888260016373760
