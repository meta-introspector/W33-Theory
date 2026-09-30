"""Regression for Pass 11190: sign structure of AME(10,3) stabilizer states (one graph recomputed; the CP-SAT proof and
automorphism group from the frozen certificate)."""
import itertools
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def test_one_graph_recomputed():
    import w33_pass11175_scan_five_qutrit as F
    import w33_pass11190_ame10_sign_structure as P
    w = json.loads(P.TABU.read_text())["graphs"][0]
    r = P.analyse_code(F.to_mat(np.array(w)))
    assert r["sign_sums_ok"] and r["local_graphs_symmetric"] and r["local_graph_types"] == {"prism+K1": 120}
    assert r["uniform_six_sets"] == 30 and r["uniform_four_sets_are_S3_4_10"] and r["duality_violations"] == 0
    assert r["signs_predict_gate_dets"]
    assert sum(v for k, v in r["patterns"].items() if k.startswith("cross+perm")) == 180


def test_frozen_proof():
    d = json.loads((ROOT / "data" / "w33_pass11190_ame10_sign_structure.json").read_text())
    assert d["local_graph_types"] == {"prism+K1": 8520} and d["all_S3_4_10"] and d["signs_predict_all_gate_dets"]
    assert d["degree_lemma"]["graphs_with_degrees_0_3_6"] == 1052
    q = d["quadric_types"]
    assert q["empty"]["span_type"] == "O-" and q["K33+K1"]["span_type"] == "O-"
    assert d["csp_some_K4"] == "INFEASIBLE" and d["csp_K4_at_012_one_worker"] == "INFEASIBLE"
    assert d["csp_prism_only_control"] in ("OPTIMAL", "FEASIBLE")
    p = d["csp_patterns"]
    assert p["C4+C6"] == p["double-cross"] == p["all"] == "INFEASIBLE"
    assert p["C10"] in ("OPTIMAL", "FEASIBLE") and p["cross+perm"] in ("OPTIMAL", "FEASIBLE")
    assert d["structures_with_fixed_steiner_system"][1] == 2


def test_block_containment_counts():
    blocks_max = 1
    # two 4-subsets of a 5-set share 3 points, so a Steiner S(3,4,10) puts at most one block in any 5-set
    for a, b in itertools.combinations(itertools.combinations(range(5), 4), 2):
        assert len(set(a) & set(b)) == 3
    assert 30 * 6 == 180 and 252 - 180 == 72 and blocks_max == 1
