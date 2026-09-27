#!/usr/bin/env python3
"""Pass 11061: the ten diagonal-weld factors are four Hesse directions x two slopes plus two central slopes."""
from __future__ import annotations
import argparse,json,sys
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_pass11056_diagonal_weld_ten_triangle_factors as FCT

OUT=ROOT/"data/w33_pass11061_hesse_slope_factor_alphabet.json"

def pdir(a,b):
    a%=3;b%=3
    if a:
        z=pow(a,-1,3);return (1,b*z%3)
    if b:return (0,1)
    return None

def payload():
    coords,v,D,pairs,factors=FCT.geometry()
    rows=[]
    for i,(s,si) in enumerate(pairs):
        (a,b,c),ext=s
        d=pdir(a,b)
        rows.append({
          "factor":i,
          "generator":[list(s[0]),s[1]],
          "inverse":[list(si[0]),si[1]],
          "projective_hesse_direction":list(d) if d is not None else None,
          "external_branch":1 if ext==1 else -1,
          "central":d is None
        })
    outer=[r for r in rows if not r["central"]]
    central=[r for r in rows if r["central"]]
    assert len(outer)==8 and len(central)==2

    dirs=sorted({tuple(r["projective_hesse_direction"]) for r in outer})
    assert dirs==[(0,1),(1,0),(1,1),(1,2)]
    for d in dirs:
        rr=[r for r in outer if tuple(r["projective_hesse_direction"])==d]
        assert len(rr)==2 and {r["external_branch"] for r in rr}=={-1,1}

    old=json.loads((ROOT/"data/w33_pass11048_induced_h27_hesse_latent_modules.json").read_text())
    olddirs={pdir(*x["direction"]) for x in old["four_directions"].values()}
    assert set(dirs)==olddirs

    def symp(r,s):
        a,b,_=r["generator"][0];d,e,_=s["generator"][0]
        return (a*e-b*d)%3
    edges=[]
    deg=[0]*10
    for i in range(10):
        for j in range(i+1,10):
            if symp(rows[i],rows[j])==0:
                edges.append([i,j]);deg[i]+=1;deg[j]+=1
    assert len(edges)==21
    assert Counter(deg)==Counter({3:8,9:2})

    return {
      "schema":"w33.pass11061.hesse-slope-factor-alphabet.v1",
      "status":"PASS_TEN_WELD_FACTORS_EQUAL_FOUR_HESSE_DIRECTIONS_TIMES_TWO_SLOPES_PLUS_TWO_CENTRAL_SLOPES",
      "headline":"The ten order-three triangle factors of the diagonal cubic weld are not an unstructured ten-set. Eight are exactly the four projective qutrit/Hesse directions in P1(F3), each lifted with the two nonzero external C3 slopes. The remaining two are the two pure-center/external slopes. Thus 10=4x2+2, and the four noncentral direction labels are exactly those used by the induced-H27 Hesse compiler.",
      "decomposition":{
        "formula":"10 = 4*2 + 2",
        "indexing_set":"P1(F3) x {+,-} disjoint_union {central+ , central-}",
        "hesse_directions":[[0,1],[1,0],[1,1],[1,2]],
        "outer_factors":8,"central_factors":2
      },
      "factor_rows":rows,
      "commutation_graph":{
        "vertices":10,"symplectic_zero_edges":21,"symplectic_nonzero_pairs":24,
        "degree_histogram":{"3":8,"9":2},
        "structure":"K2 join (4 disjoint K2): the central pair commutes with all, and each Hesse slope pair commutes internally"
      },
      "boundary":"The +/- labels are the frozen external C3 branches in the current generator convention. Their physical interpretation as time orientation or chirality requires additional dynamics.",
      "parents":["data/w33_pass11056_diagonal_weld_ten_triangle_factors.json","data/w33_pass11048_induced_h27_hesse_latent_modules.json"],
      "checks":{"eight_outer_plus_two_central":True,"four_hesse_directions_exact":True,"two_branches_per_hesse_direction":True,"commutation_graph_21_24":True}
    }

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args()
    p=payload();text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or a.output.read_text()!=text:raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"outer":8,"central":2},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
