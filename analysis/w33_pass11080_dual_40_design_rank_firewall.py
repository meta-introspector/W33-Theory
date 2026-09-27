#!/usr/bin/env python3
"""Pass 11080: point-side and line-side 2-(40,13,4) designs have different ternary ranks."""
from __future__ import annotations
import argparse,itertools,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/w33_pass11080_dual_40_design_rank_firewall.json"
def canon(v):
    v=tuple(int(x)%3 for x in v)
    for x in v:
        if x:
            z=pow(x,-1,3);return tuple((z*y)%3 for y in v)
    raise ValueError
def symp(x,y):return (x[0]*y[2]-x[2]*y[0]+x[1]*y[3]-x[3]*y[1])%3
def rank_mod(A,p=3):
    M=np.array(A,dtype=np.int64)%p;n,m=M.shape;r=0
    for c in range(m):
        z=next((i for i in range(r,n) if M[i,c]%p),None)
        if z is None:continue
        M[[r,z]]=M[[z,r]]
        iv=pow(int(M[r,c]),-1,p);M[r]=(M[r]*iv)%p
        for i in range(n):
            if i!=r and M[i,c]%p:
                q=int(M[i,c]);M[i]=(M[i]-q*M[r])%p
        r+=1
        if r==n:break
    return r
def payload():
    pts=sorted({canon(v) for v in itertools.product(range(3),repeat=4) if any(v)})
    Ap=np.zeros((40,40),dtype=np.int64)
    for i,j in itertools.combinations(range(40),2):
        if symp(pts[i],pts[j])==0:Ap[i,j]=Ap[j,i]=1
    lines=[frozenset(c) for c in itertools.combinations(range(40),4)
           if all(Ap[i,j] for i,j in itertools.combinations(c,2))]
    assert len(lines)==40
    Al=np.zeros((40,40),dtype=np.int64)
    for i,j in itertools.combinations(range(40),2):
        if len(lines[i]&lines[j])==1:Al[i,j]=Al[j,i]=1
    I=np.eye(40,dtype=np.int64);J=np.ones((40,40),dtype=np.int64)
    Np=I+Ap;Nl=I+Al
    assert np.array_equal(Np@Np.T,9*I+4*J)
    assert np.array_equal(Nl@Nl.T,9*I+4*J)
    rp,rl=rank_mod(Np),rank_mod(Nl)
    assert (rp,rl)==(11,15)
    leech=json.loads((ROOT/"data/PART_W33_PASS7717_7724_LEECH_DUAL40_TERNARY_CODE.json").read_text())
    pol=json.loads((ROOT/"data/PART_W33_PASS7709_7716_LEECH_DUAL40_PROJECTIVE_POLARITY.json").read_text())
    par=json.loads((ROOT/"data/w33_pass11070_dual_parabolic_flag_diamond.json").read_text())
    assert leech["incidence_rank_F3"]==11
    assert pol["cross_design"]["parameters"]=="symmetric 2-(40,13,4)"
    assert par["parabolics"]["point"]["radical"]=="H27" and par["parabolics"]["line"]["radical"]=="F3^3"
    return {
      "schema":"w33.pass11080.dual-40-design-rank-firewall.v1",
      "status":"PASS_POINT_AND_LINE_TEMPORAL_BOUNDARY_DESIGNS_SHARE_PARAMETERS_BUT_ARE_NOT_THE_SAME_TERNARY_DESIGN",
      "headline":"The point closed-neighborhood design and the line closed-neighborhood/temporal-boundary design are both symmetric 2-(40,13,4) designs and have the same rational Gram identity 9I+4J, but their F3 incidence ranks are 11 and 15 respectively. Ternary rank is an isomorphism invariant, so the two designs are not permutation-equivalent. This is the design-level shadow of the non-self-dual W33 parabolics.",
      "point_side":{"objects":"W33 points","design":"2-(40,13,4)","rank_F3":11,
                    "code":"[40,11,13]_3 from the existing Leech/PG(3,3) certificate",
                    "local_radical":"H27"},
      "line_side":{"objects":"W33 lines / Bell-chart centers","design":"2-(40,13,4)","rank_F3":15,
                   "meaning":"13-context temporal boundary design","local_radical":"F3^3"},
      "shared_integer_geometry":{"gram":"N N^T=9I+4J","rational_singular_values":"13 once and 3 with multiplicity 39","same_BIBD_parameters":True},
      "firewall":{"rank_difference":"11 != 15","designs_isomorphic":False,
                  "consequence":"The Leech dual-40 point-hyperplane design cannot be silently identified with the Bell-line temporal-boundary design merely from 40,13,4 counts. An explicit point-line duality datum would be required, and the rank obstruction rules out a permutation design isomorphism."},
      "structural_reading":"The earlier dual-648 theorem has a coding/design avatar: point memory and line history have identical coarse design parameters but different characteristic-three linearizations.",
      "boundary":"This is an exact p-rank obstruction. It does not rule out weaker correspondences, correspondences after changing coefficients, or a larger categorical relation between the two sides.",
      "parents":["data/w33_pass11079_bell_chart_history_boundary_design.json","data/PART_W33_PASS7717_7724_LEECH_DUAL40_TERNARY_CODE.json","data/w33_pass11070_dual_parabolic_flag_diamond.json"],
      "checks":{"both_symmetric_2_40_13_4":True,"same_rational_gram":True,"point_rank11":True,"line_rank15":True,"nonisomorphic_over_F3":True}
    }
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args();p=payload()
    text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or json.loads(a.output.read_text())!=p:raise SystemExit("certificate drift")
    else:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"ranks":[11,15]},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
