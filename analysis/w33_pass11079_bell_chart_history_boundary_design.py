#!/usr/bin/env python3
"""Pass 11079: Bell-line temporal charts form dual symmetric designs and a regular 39-simplex."""
from __future__ import annotations
import argparse,itertools,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/w33_pass11079_bell_chart_history_boundary_design.json"
def canon(v):
    v=tuple(int(x)%3 for x in v)
    for x in v:
        if x:
            z=pow(x,-1,3);return tuple((z*y)%3 for y in v)
    raise ValueError
def symp(x,y):return (x[0]*y[2]-x[2]*y[0]+x[1]*y[3]-x[3]*y[1])%3
def rank_mod(A,p):
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
    adj=[[False]*40 for _ in range(40)]
    for i,j in itertools.combinations(range(40),2):
        if symp(pts[i],pts[j])==0:adj[i][j]=adj[j][i]=True
    lines=[frozenset(c) for c in itertools.combinations(range(40),4)
           if all(adj[i][j] for i,j in itertools.combinations(c,2))]
    assert len(lines)==40
    A=np.zeros((40,40),dtype=np.int64)
    for i,j in itertools.combinations(range(40),2):
        if len(lines[i]&lines[j])==1:A[i,j]=A[j,i]=1
    N=np.eye(40,dtype=np.int64)+A
    H=np.ones((40,40),dtype=np.int64)-N
    I=np.eye(40,dtype=np.int64);J=np.ones((40,40),dtype=np.int64)
    assert set(map(int,N.sum(1)))=={13} and set(map(int,H.sum(1)))=={27}
    assert np.array_equal(N@N.T,9*I+4*J)
    assert np.array_equal(H@H.T,9*I+18*J)
    assert np.array_equal(N+H,J)
    W=40*N-13*J;WH=40*H-27*J
    assert np.array_equal(WH,-W)
    assert np.array_equal(W@W.T,360*(40*I-J))
    assert rank_mod(W,101)==39
    assert rank_mod(N,3)==15 and rank_mod(H,3)==14
    return {
      "schema":"w33.pass11079.bell-chart-history-boundary-design.v1",
      "status":"PASS_FORTY_BELL_LINE_CHARTS_FORM_ANTIPODAL_REGULAR_SIMPLEX_DESIGNS",
      "headline":"For each W33 line B, the 27 lines disjoint from B are its finite history chart and the complementary 13 lines (B plus the 12 meeting B) are its boundary. Across all 40 centers these are symmetric 2-(40,27,18) and 2-(40,13,4) designs. After centering, the 40 boundary indicators are the vertices of a regular 39-simplex; the history indicators are exactly their antipodes.",
      "history_design":{"parameters":"2-(40,27,18) symmetric","block_size":27,"pair_intersection":18,
                        "gram":"H H^T = 9 I + 18 J","rank_F3":14},
      "boundary_design":{"parameters":"2-(40,13,4) symmetric","block_size":13,"pair_intersection":4,
                         "gram":"N N^T = 9 I + 4 J","rank_F3":15},
      "simplex":{"boundary_centered":"W=40N-13J","history_centered":"WH=40H-27J=-W",
                 "gram":"W W^T = 360(40I-J)","dimension":39,
                 "normalized_off_diagonal_inner_product":"-1/39",
                 "reading":"no Bell-line chart is metrically preferred inside the centered chart-frame geometry"},
      "chart_partition":{"for_each_center":"40 = 27 finite histories + 13 boundary contexts","boundary":"center plus its 12 intersecting W33 lines"},
      "boundary":"The regular-simplex statement is finite design geometry. It is a precise covariance property of the chart atlas, not a derivation of continuum general covariance.",
      "checks":{"history_design_exact":True,"boundary_design_exact":True,"centered_antipodes":True,"regular_39_simplex":True}
    }
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args();p=payload()
    text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or json.loads(a.output.read_text())!=p:raise SystemExit("certificate drift")
    else:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"designs":["40,27,18","40,13,4"],"simplex":39},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
