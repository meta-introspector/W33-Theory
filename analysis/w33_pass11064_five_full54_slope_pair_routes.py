#!/usr/bin/env python3
"""Pass 11064: every Hesse slope-pair and the central slope-pair independently covers the 54D quotient."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_pass11047_phase_weld_exact_54_right_inverse as PREV
import w33_pass11062_single_foliation_finite54_transducer as ONE

OUT=ROOT/"data/w33_pass11064_five_full54_slope_pair_routes.json"
LABELS=[["Hesse",[1,0]],["Hesse",[1,2]],["Hesse",[1,1]],["Hesse",[0,1]],["central",None]]

def payload():
    bridge=json.loads((ROOT/"data/w33_e6id_current_h27_gauge_bridge.json").read_text())
    h={int(i):tuple(map(int,x)) for i,x in bridge["maps"]["e6id_to_current_H27_address"].items()}
    cert={}
    for which in ("plus","minus"):
        coords,pairs,factors=ONE.factorize(which)
        prime_rows=[]
        for p in (103,109):
            F,Finv,_=PREV.full_fourier_matrix(p,h)
            rows=[]
            for j in range(5):
                M=factors[2*j]+factors[2*j+1]
                C=(Finv@(M%p))%p
                raw=PREV.inv_rank(M%p,p,False)[0]
                q=PREV.inv_rank(C[27:,:],p,False)[0]
                s2=PREV.inv_rank(C[27:54,:],p,False)[0]
                l=PREV.inv_rank(C[54:,:],p,False)[0]
                assert q==54 and s2==l==27
                rows.append({"route":j,"factors":[2*j,2*j+1],"raw_rank":raw,"quotient_rank":q,
                             "S1_intersection_dimension":raw-q,"S2_rank":s2,"L_rank":l})
            prime_rows.append({"prime":p,"rows":rows})
        assert [r["raw_rank"] for r in prime_rows[0]["rows"]]==[66,54,54,72,72]
        assert prime_rows[0]["rows"]==prime_rows[1]["rows"]
        cert[which]=prime_rows

    route_summary=[]
    raw=[66,54,54,72,72]
    for j,label in enumerate(LABELS):
        route_summary.append({
          "route":j,"type":label[0],"projective_direction":label[1],
          "factors":[2*j,2*j+1],"triangle_blocks":54,
          "raw_rank":raw[j],"quotient_rank":54,"S1_intersection_dimension":raw[j]-54
        })
    return {
      "schema":"w33.pass11064.five-full54-slope-pair-routes.v1",
      "status":"PASS_FIVE_TWO_SLOPE_ROUTE_FAMILIES_EACH_SURJECT_ONTO_THE_FULL_54D_RETYPE_QUOTIENT",
      "headline":"Pair the two external branches above each of the four Hesse directions, and pair the two central branches. Every one of these five two-factor sums projects with rank 54 onto S2+L for both FI orientations. The raw ranks are 66,54,54,72,72, so two Hesse pair-routes are already transverse 54-planes while the other three carry controlled S1 overlap.",
      "routes":route_summary,
      "raw_rank_pattern":[66,54,54,72,72],
      "quotient_rank_pattern":[54,54,54,54,54],
      "local_block_budget":{"factors_per_route":2,"triangles_per_factor":27,"triangle_blocks_per_route":54},
      "split_prime_certificates":cert,
      "boundary":"These are five algebraically redundant retyping routes. Their pair sums are noncommuting weighted generators; a two-layer product of individual exponentials is not asserted to equal the Cayley or exponential of the pair sum.",
      "parents":["data/w33_pass11061_hesse_slope_factor_alphabet.json","data/w33_pass11062_single_foliation_finite54_transducer.json"],
      "checks":{"five_routes":True,"all_quotient_rank54":True,"both_orientations_same_rank_pattern":True,"two_transverse_pair_routes":True}
    }

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args()
    p=payload();text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or a.output.read_text()!=text:raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"routes":5},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
