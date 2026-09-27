#!/usr/bin/env python3
"""Pass 11070: the two order-648 W33 parabolics meet in the U81 chamber weld."""
from __future__ import annotations
import argparse,json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/w33_pass11070_dual_parabolic_flag_diamond.json"

def payload():
    p1047=json.loads((ROOT/"data/w33_pass1047_two_648_stabilizers.json").read_text())
    p4519=json.loads((ROOT/"data/PART_W33_PASS4519_FLAG_BOREL_SYLOW3_NORMALIZER.json").read_text())
    p5105=json.loads((ROOT/"data/PART_W33_PASS5105_U81_DUAL_TORSOR_CONTROLLER.json").read_text())

    assert p1047["group_order"]==25920
    assert p1047["point_stabilizer"]["order"]==p1047["line_stabilizer"]["order"]==648
    assert p4519["flag"]["order"]==p4519["normalizer"]["order"]==162
    assert p4519["normalizer"]["equals_flag_stabilizer"]
    assert p4519["sylow3"]["order"]==81 and p4519["sylow3"]["index_in_flag"]==2
    assert p5105["state_H27"]["order"]==27 and p5105["state_H27"]["center"]==3
    assert p5105["program_F3_3"]["order"]==27 and p5105["program_F3_3"]["center"]==27
    assert p5105["weld"]["intersection_order"]==9
    assert p5105["weld"]["intersection_equals_U_derived"]
    assert p5105["weld"]["generated_product_order"]==81
    assert p5105["U"]=={"order":81,"center":3,"derived":9,"order_census":{"1":1,"3":44,"9":36}}

    flags=40*4
    assert flags==25920//162==160
    product_order=27*27//9
    assert product_order==81

    return {
      "schema":"w33.pass11070.dual-parabolic-flag-diamond.v1",
      "status":"PASS_POINT648_AND_LINE648_MEET_IN_THE_FLAG_BOREL_AND_THEIR_ORDER27_RADICALS_WELD_TO_U81",
      "headline":(
        "The dual order-648 stabilizers are not competing descriptions of one local group. "
        "For an incident point p on line L, P=G_p and Q=G_L intersect in the flag stabilizer "
        "B=G_{p<L} of order 162. Its unique Sylow-3 subgroup U has order 81. The canonical "
        "normal order-27 radical on the point side is extraspecial H27; the line-side radical "
        "is flat F3^3. Inside U they intersect in U' of order 9 and generate all of U."
      ),
      "ambient":{"group":"PSp(4,3)","order":25920,"points":40,"lines":40,"flags":flags},
      "parabolics":{
        "point":{"order":648,"structure":"3^(1+2):2A4","radical":"H27","radical_order":27,
                 "center_order":p1047["point_stabilizer"]["center_order"]},
        "line":{"order":648,"structure":"3^3:S4","radical":"F3^3","radical_order":27,
                "center_order":p1047["line_stabilizer"]["center_order"]}
      },
      "flag_borel":{
        "identity":"B=G_p intersect G_L",
        "order":162,"index":160,
        "sylow3":"U81","sylow3_order":81,"quotient_order":2,
        "normalizer_identity":"B=N_G(U81)"
      },
      "radical_weld":{
        "H_state_order":27,"H_program_order":27,
        "intersection_order":9,"intersection":"U81 derived subgroup",
        "generated_order":product_order,"generated_group":"U81",
        "center_U81_order":3,
        "equation":"|H27 F3^3|=27*27/9=81"
      },
      "structural_reading":(
        "The point radical supplies the noncommutative qutrit-state torsor and the line radical "
        "supplies the commuting history/program torsor. A chamber does not choose one or the other: "
        "its Sylow-3 update group is their nontrivial product, glued along the common order-9 derived layer."
      ),
      "boundary":(
        "This is a synthesis of previously certified finite-group computations. The words state, "
        "program, history and update are architectural readings; no continuum spacetime or hardware timing is inferred."
      ),
      "parents":[
        "data/w33_pass1047_two_648_stabilizers.json",
        "data/PART_W33_PASS4519_FLAG_BOREL_SYLOW3_NORMALIZER.json",
        "data/PART_W33_PASS5105_U81_DUAL_TORSOR_CONTROLLER.json"
      ],
      "checks":{
        "dual_parabolics_order648":True,"flag_intersection_order162":True,
        "flag_is_sylow3_normalizer":True,"U81_order81":True,
        "point_radical_H27_order27":True,"line_radical_F3_3_order27":True,
        "radical_intersection_order9":True,"radicals_generate_U81":True
      }
    }

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args();p=payload()
    text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or json.loads(a.output.read_text())!=p:raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"diamond":[648,162,81,27,27,9]},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
