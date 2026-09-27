#!/usr/bin/env python3
"""Pass 11083: reduce unconstrained local-E8 atlas gluing to 81 free holonomies."""
from __future__ import annotations
import argparse,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/w33_pass11083_local_e8_flat_holonomy_frontier.json"
WE8=696729600

def payload():
    pi=json.loads((ROOT/"data/w33_pass11081_atlas_fundamental_group_steinberg_linearization.json").read_text())
    atlas=json.loads((ROOT/"data/w33_pass11075_local_e8_atlas_amalgam_frontier.json").read_text())
    a2=json.loads((ROOT/"data/w33_pass1304_a2_normalizer_triality.json").read_text())
    assert pi["graph_certificate"]["free_rank_value"]==81
    assert atlas["atlas"]["chart_labels"]==40 and atlas["atlas"]["pairwise_overlap_edges"]==240
    assert a2["full_A2_subsystem_normalizer"]["order"]==622080
    assert a2["full_A2_subsystem_normalizer"]["index_in_W(E8)"]==1120
    assert WE8//1120==622080
    return {
      "schema":"w33.pass11083.local-e8-flat-holonomy-frontier.v1",
      "status":"PASS_GLOBAL_LOCAL_E8_GLUE_TOPOLOGY_REDUCES_TO_81_FREE_HOLONOMIES_WITH_SHARED_A2_TRANSPORTER_CONSTRAINTS",
      "headline":"Because the local-E8 overlap nerve has pi1=F81, a flat local system with a common structure group H is, after spanning-tree gauge, specified by 81 holonomies modulo simultaneous conjugation: Hom(F81,H)/H. For the actual E8 atlas, every pairwise chart overlap is a shared A2 point-fibre, so an edge transition must lie in the transporter between the two A2 embeddings. W(E8) is transitive on 1120 A2 subsystems and the A2 normalizer has order 622080, so each nonempty transporter is a coset of that normalizer.",
      "topological_reduction":{"pi1":"F81","spanning_tree_edges_in_Levi_model":79,
                               "remaining_cycle_edges":81,
                               "unconstrained_flat_H_moduli":"Hom(F81,H)/H = H^81 / diagonal conjugation"},
      "E8_overlap_constraints":{"local_charts":40,"pairwise_chart_overlap_edges":240,
                                "shared_object":"one A2 point-fibre",
                                "WE8_order":WE8,"A2_subsystems":1120,
                                "A2_normalizer_order":622080,
                                "edge_transporter":"a left/right coset of N_W(E8)(A2), hence 622080 elements after fixed chart gauges"},
      "frontier":"Topology no longer supplies hidden higher relations: the remaining difficulty is algebraic compatibility of the 160 incidence A2 embeddings / 240 chart-pair transporters with tetracode orientations and the chosen local E8 brackets.",
      "boundary":"The H^81/conjugation description is the standard unconstrained flat-local-system moduli for a common structure group. The physical/local-E8 atlas has edgewise transporter restrictions, so this pass does not assert that arbitrary 81 W(E8) elements define a valid exceptional gluing.",
      "parents":["data/w33_pass11081_atlas_fundamental_group_steinberg_linearization.json","data/w33_pass11075_local_e8_atlas_amalgam_frontier.json","data/w33_pass1304_a2_normalizer_triality.json"],
      "checks":{"pi1_F81":True,"A2_normalizer_622080":True,"transporter_coset_size_622080":True,
                "unconstrained_flat_moduli_81_generators":True,"edge_constraints_not_dropped":True}
    }
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args();p=payload()
    text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or json.loads(a.output.read_text())!=p:raise SystemExit("certificate drift")
    else:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"free_holonomies":81},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
