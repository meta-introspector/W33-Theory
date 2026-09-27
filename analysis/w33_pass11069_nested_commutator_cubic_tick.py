#!/usr/bin/env python3
"""Pass 11069: the cocycle parameter switches on the exact primitive -> area -> central-tick ladder."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_pass11068_ternary_cocycle_deformation as DEF

OUT=ROOT/"data/w33_pass11069_nested_commutator_cubic_tick.json"
E0=(1,0,0,0);E1=(0,1,0,0);ID=(0,0,0,0)

def payload():
    rows={}
    expected={
      0:((0,0,1,0),(0,0,0,0)),
      1:((0,0,1,1),(0,0,0,2)),
      2:((0,0,1,2),(0,0,0,1))
    }
    for s in range(3):
        first=DEF.comm(E0,E1,s)
        second=DEF.comm(E0,first,s)
        third=DEF.comm(E0,second,s)
        assert (first,second)==expected[s] and third==ID
        rows[str(s)]={"first_commutator":list(first),"nested_commutator":list(second),"third_nested":list(third)}

    return {
      "schema":"w33.pass11069.nested-commutator-cubic-tick.v1",
      "status":"PASS_NONZERO_COCYCLE_TWIST_SWITCHES_ON_THE_CLASS3_CENTRAL_TICK",
      "headline":"Take the two primitive root directions e0=(1,0,0,0) and e1=(0,1,0,0). For the split law s=0, their commutator is the H27 area/phase coordinate (0,0,1,0) and the next commutator vanishes. For a nonzero cocycle twist, [e0,e1]=(0,0,1,s) and [e0,[e0,e1]]=(0,0,0,-s), a pure highest-root central coordinate. A third nesting vanishes.",
      "closed_formula":{
        "primitive_generators":["e0=(1,0,0,0)","e1=(0,1,0,0)"],
        "first":"[e0,e1]_s=(0,0,1,s)",
        "second":"[e0,[e0,e1]_s]_s=(0,0,0,-s)",
        "third":"[e0,[e0,[e0,e1]]]=1"
      },
      "scalar_replays":rows,
      "graded_reading":{
        "degree1":"two primitive directions e0,e1",
        "degree2":"H27 commutator/area coordinate c",
        "degree3":"highest-root central coordinate d",
        "split_case":"degree3 tick is off",
        "nonsplit_case":"degree3 tick is nonzero and changes sign under s -> -s"
      },
      "connection_to_jennings":"Pass11066 independently records the same degree replacement: K81 has degrees 1,1,1,2 while U81 has 1,1,2,3.",
      "boundary":"The terminal commutator is an exact ternary central group coordinate. Calling it a temporal tick is a structural interpretation, not a measurement of continuous physical time.",
      "parents":["data/w33_pass11068_ternary_cocycle_deformation.json","data/w33_pass11066_jennings_clock_promotion.json"],
      "checks":{"split_nested_commutator_zero":True,"nonzero_twists_have_pure_second_nested_tick":True,"orientation_reverses_tick":True,"third_nested_zero":True}
    }

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args()
    p=payload();text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or a.output.read_text()!=text:raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"ladder":"1->2->3"},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
