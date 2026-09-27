#!/usr/bin/env python3
"""Pass 11075: local W33 line tetracodes give exact E8 charts, while global Lie gluing remains an amalgam problem."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_tetracode_e8_root_system_bridge as E8

OUT=ROOT/"data/w33_pass11075_local_e8_atlas_amalgam_frontier.json"

def payload():
    cover=json.loads((ROOT/"data/w33_pass11071_oriented_chamber_tetracode_cover.json").read_text())
    nerve=json.loads((ROOT/"data/w33_pass11074_local_e8_chart_overlap_nerve.json").read_text())
    packet=E8.tetracode_e8_root_system_packet()
    assert packet["checks"]["w33_code_equals_standard_tetracode"]
    assert packet["root_system"]["count"]==240 and packet["root_system"]["rank"]==8
    assert packet["checks"]["reflection_closure_holds"]
    assert packet["simple_root_system"]["gram_determinant"]=="1"
    assert cover["line_tetracode"]["nonzero_words"]==8
    assert nerve["homology"]["betti"]==[1,81,0,0]

    return {
      "schema":"w33.pass11075.local-e8-atlas-amalgam-frontier.v1",
      "status":"PASS_EACH_LINE_TETRACODE_HAS_AN_EXACT_E8_LIFT_AND_THE_40_CHART_OVERLAP_TOPOLOGY_IS_EXACT_WHILE_GLOBAL_E8_AMALGAM_REMAINS_OPEN",
      "headline":(
        "The local and global statements can now be separated sharply. The W33-derived tetracode has an exact "
        "four-A2 Eisenstein lift to a rank-8, 240-root, reflection-closed E8 root system. Every W33 line has the "
        "same four-point/tetracode combinatorics and the 40 line labels form one PSp orbit, so they support 40 "
        "coordinate-equivalent local E8 chart constructions. Their shared-point overlap nerve is exact and has H1=Z^81. "
        "What is not yet proved is a compatible family of A2->E8 embedding maps whose amalgam is one global E8 Lie algebra."
      ),
      "local_chart":{
        "input":"four W33 line points + evaluation tetracode [4,2,3]_3",
        "oriented_tensor_labels":cover["line_tetracode"]["nonzero_words"],
        "E8_root_count":packet["root_system"]["count"],
        "rank":packet["root_system"]["rank"],
        "reflection_closed":packet["checks"]["reflection_closure_holds"],
        "simple_gram_determinant":packet["simple_root_system"]["gram_determinant"],
        "construction":"four A2 root planes + tetracode glue"
      },
      "atlas":{
        "chart_labels":40,"group_orbit":"PSp(4,3) transitive on W33 lines",
        "pairwise_overlap_edges":nerve["nerve"]["edges_pairwise_point_overlaps"],
        "maximal_fourfold_overlaps":nerve["nerve"]["tetrahedra_fourfold_point_overlaps"],
        "overlap_H1":"Z^81","overlap_betti":[1,81,0,0]
      },
      "exact_vs_open":{
        "exact":[
          "each line has four points and eight oriented tetracode words",
          "the standard tetracode lift gives an exact E8 root system",
          "all 40 line labels are equivalent under the W33 automorphism action",
          "the shared-point chart-overlap nerve has H1=Z^81"
        ],
        "open":[
          "choose explicit A2 embeddings into each local E8 chart compatible on all shared W33 points",
          "compute the resulting universal completion / colimit Lie algebra",
          "decide whether the completion is one E8, a quotient, a larger finite algebra, or infinite",
          "identify whether the 81 overlap cycles carry a nontrivial Lie-algebra cocycle rather than only nerve topology"
        ]
      },
      "amalgam_problem":(
        "The next object is a diagram with one abstract E8 chart for each W33 line and one A2 point-fibre for each incidence. "
        "Pass11074 determines the diagram's nerve exactly, but not the embeddings. The global TOE claim must be tested at the "
        "level of this explicit amalgam/colimit; count agreement is not enough."
      ),
      "boundary":(
        "Forty local coordinate-equivalent E8 constructions do not mean forty physical E8 gauge groups, nor do they prove "
        "that the charts glue into a single global E8. This pass deliberately freezes that distinction."
      ),
      "parents":["analysis/w33_tetracode_e8_root_system_bridge.py","data/w33_pass11071_oriented_chamber_tetracode_cover.json","data/w33_pass11074_local_e8_chart_overlap_nerve.json"],
      "checks":{"local_E8_240_roots_rank8":True,"local_E8_reflection_closed":True,
                "forty_line_chart_labels":True,"overlap_H1_rank81":True,
                "global_amalgam_not_silently_asserted":True}
    }

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args();p=payload()
    text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or json.loads(a.output.read_text())!=p:raise SystemExit("certificate drift")
    else:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"charts":40,"overlap_H1":81},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
