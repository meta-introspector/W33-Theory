#!/usr/bin/env python3
"""Pass 11082: every pair of Bell-line charts has the same 18/9/9/4 transition census."""
from __future__ import annotations
import argparse,itertools,json
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/w33_pass11082_bell_chart_transition_census.json"

def canon(v):
    v=tuple(int(x)%3 for x in v)
    for x in v:
        if x:
            z=pow(x,-1,3);return tuple((z*y)%3 for y in v)
    raise ValueError
def symp(x,y):return (x[0]*y[2]-x[2]*y[0]+x[1]*y[3]-x[3]*y[1])%3
def payload():
    pts=sorted({canon(v) for v in itertools.product(range(3),repeat=4) if any(v)})
    adj=[[False]*40 for _ in range(40)]
    for i,j in itertools.combinations(range(40),2):
        if symp(pts[i],pts[j])==0:adj[i][j]=adj[j][i]=True
    lines=[frozenset(c) for c in itertools.combinations(range(40),4)
           if all(adj[i][j] for i,j in itertools.combinations(c,2))]
    assert len(lines)==40
    meet=[[False]*40 for _ in range(40)]
    for i,j in itertools.combinations(range(40),2):
        if lines[i]&lines[j]:meet[i][j]=meet[j][i]=True
    charts=[]
    for b in range(40):
        boundary={b}|{j for j in range(40) if meet[b][j]}
        finite=set(range(40))-boundary
        assert (len(finite),len(boundary))==(27,13)
        charts.append((finite,boundary))
    hist=Counter()
    for b,c in itertools.combinations(range(40),2):
        hb,bb=charts[b];hc,bc=charts[c]
        key=(len(hb&hc),len(hb&bc),len(bb&hc),len(bb&bc))
        hist[key]+=1
    assert hist==Counter({(18,9,9,4):780})
    return {
      "schema":"w33.pass11082.bell-chart-transition-census.v1",
      "status":"PASS_EVERY_DISTINCT_BELL_CHART_PAIR_HAS_THE_SAME_18_9_9_4_FINITE_BOUNDARY_TRANSITION_CENSUS",
      "headline":"For any two distinct Bell-line centers B and C, the 40 global line contexts split uniformly into 18 finite in both charts, 9 finite in B but boundary in C, 9 boundary in B but finite in C, and 4 boundary in both. All 780 unordered chart pairs have exactly the same census.",
      "transition_partition":{"finite_B_finite_C":18,"finite_B_boundary_C":9,
                              "boundary_B_finite_C":9,"boundary_B_boundary_C":4,
                              "sum":40},
      "pair_census":{"unordered_chart_pairs":780,"distinct_transition_profiles":1,
                     "profile":"18/9/9/4","pairs_with_profile":780},
      "design_derivation":"The boundary blocks form 2-(40,13,4). Thus two boundaries meet in 4; each has 9 exclusive contexts; their complement has 40-(13+13-4)=18 contexts.",
      "interpretation":"Changing Bell polarization has a fixed finite cost: exactly nine currently finite histories cross the old/new chart boundary in each direction, while eighteen remain finite and four remain boundary.",
      "boundary":"This is a finite chart-incidence theorem. The words finite and boundary refer to the Lagrangian big-cell atlas, not a continuum causal horizon.",
      "parents":["data/w33_pass11079_bell_chart_history_boundary_design.json"],
      "checks":{"all780_pairs_checked":True,"unique_profile":True,"profile_18_9_9_4":True}
    }
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args();p=payload()
    text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or json.loads(a.output.read_text())!=p:raise SystemExit("certificate drift")
    else:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"profile":[18,9,9,4]},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
