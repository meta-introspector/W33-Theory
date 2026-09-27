#!/usr/bin/env python3
"""Pass 11056: diagonal cubic weld support is ten order-3 Cayley triangle factors."""
from __future__ import annotations
import argparse,itertools,json,sys
from collections import Counter
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_e6_cubic_fourier54_alignment as ALIGN

OUT=ROOT/"data/w33_pass11056_diagonal_weld_ten_triangle_factors.json"

def hmul(x,y):
    a,b,c=x; d,e,f=y
    return ((a+d)%3,(b+e)%3,(c+f-d*b)%3)

def hinv(x):
    a,b,c=x
    return ((-a)%3,(-b)%3,(-c-a*b)%3)

def kmul(x,y): return (hmul(x[0],y[0]),(x[1]+y[1])%3)
def kinv(x): return (hinv(x[0]),(-x[1])%3)
IDENT=((0,0,0),0)

def rank_mod(A,p):
    M=np.array(A,dtype=np.int64)%p; n,m=M.shape; r=0
    for c in range(m):
        z=next((i for i in range(r,n) if M[i,c]%p),None)
        if z is None: continue
        M[[r,z]]=M[[z,r]]
        iv=pow(int(M[r,c]%p),-1,p); M[r]=(M[r]*iv)%p
        for i in range(r+1,n):
            q=int(M[i,c]%p)
            if q: M[i]=(M[i]-q*M[r])%p
        r+=1
        if r==n: break
    return r

def base():
    bridge=json.loads((ROOT/"data/w33_e6id_current_h27_gauge_bridge.json").read_text())
    h={int(i):tuple(map(int,x)) for i,x in bridge["maps"]["e6id_to_current_H27_address"].items()}
    coords=[(h[eid],phase) for eid in range(27) for phase in range(3)]
    return coords,ALIGN.ordered_records()

def diagonal(which="plus"):
    coords,records=base()
    if which=="plus": v=[1+((h[2]+phase)%3) for h,phase in coords]
    elif which=="minus": v=[1+((h[2]-phase)%3) for h,phase in coords]
    else: raise ValueError(which)
    D=np.zeros((81,81),dtype=np.int64)
    for u,x,o,c in records: D[o,x]+=c*v[u]
    return coords,v,D

def components(F):
    adj=[set(np.where(F[i]!=0)[0].tolist()) for i in range(81)]
    seen=set(); out=[]
    for i in range(81):
        if i in seen: continue
        stack=[i];seen.add(i);cc=[]
        while stack:
            q=stack.pop();cc.append(q)
            for z in adj[q]:
                if z not in seen:seen.add(z);stack.append(z)
        out.append(sorted(cc))
    return out

def geometry():
    coords,v,D=diagonal("plus"); idx={x:i for i,x in enumerate(coords)}
    ie=idx[IDENT]
    connection=[coords[j] for j in range(81) if D[ie,j]!=0]
    S=set(connection); assert len(S)==20
    for i,x in enumerate(coords):
        xi=kinv(x)
        for j,y in enumerate(coords):
            assert (D[i,j]!=0)==(kmul(xi,y) in S)

    seen=set(); pairs=[]
    for s in connection:
        if s in seen: continue
        si=kinv(s); assert si in S and si!=s
        seen.add(s);seen.add(si);pairs.append((s,si))
    assert len(pairs)==10
    for s,_ in pairs:
        assert kmul(kmul(s,s),s)==IDENT

    factors=[]
    for s,si in pairs:
        F=np.zeros_like(D)
        pair={s,si}
        for i,x in enumerate(coords):
            xi=kinv(x)
            for j,y in enumerate(coords):
                if kmul(xi,y) in pair: F[i,j]=D[i,j]
        assert np.array_equal(F.T,-F)
        cc=components(F); assert len(cc)==27 and all(len(x)==3 for x in cc)
        assert all(rank_mod(F[np.ix_(x,x)],103)==2 and rank_mod(F[np.ix_(x,x)],109)==2 for x in cc)
        factors.append(F)
    assert np.array_equal(sum(factors),D)
    return coords,v,D,pairs,factors

def payload():
    coords,v,D,pairs,factors=geometry()
    assert np.array_equal(D.T,-D)
    assert Counter(map(int,np.count_nonzero(D,axis=1)))==Counter({20:81})
    weights=Counter(abs(int(D[i,j])) for i,j in np.argwhere(np.triu(D!=0,1)))
    assert weights==Counter({1:270,2:270,3:270})
    rows=[]
    for n,((s,si),F) in enumerate(zip(pairs,factors)):
        w=Counter(abs(int(F[i,j])) for i,j in np.argwhere(np.triu(F!=0,1)))
        assert w==Counter({1:27,2:27,3:27})
        rows.append({
          "factor":n,
          "generator":[list(s[0]),s[1]],
          "inverse":[list(si[0]),si[1]],
          "order":3,"triangles":27,"vertices_per_triangle":3,
          "exact_rank":54,
          "absolute_edge_weight_histogram":{"1":27,"2":27,"3":27}
        })
    return {
      "schema":"w33.pass11056.diagonal-weld-ten-triangle-factors.v1",
      "status":"PASS_DIAGONAL_WELD_SUPPORT_IS_A_20_REGULAR_K81_CAYLEY_GRAPH_SPLIT_INTO_TEN_TRIANGLE_FACTORS",
      "headline":"The center-plus-external cubic generator has support equal to a 20-regular Cayley graph on K=H27 x C3. Its 20 connection elements form ten inverse pairs of order-three elements. Each pair gives a rank-54 skew factor made of 27 disjoint triangles, and the ten factors sum exactly to the rank-66 diagonal weld generator.",
      "support_graph":{
        "group":"K=H27 x C3_external","vertices":81,"degree":20,"undirected_edges":810,
        "connection_elements":20,"inverse_pairs":10,"all_connection_orders":3,
        "absolute_edge_weight_histogram":{"1":270,"2":270,"3":270}
      },
      "triangle_factors":rows,
      "factorization":{
        "formula":"D_plus=sum_{j=0}^9 F_j",
        "factors":10,"triangles_per_factor":27,"rank_each_factor":54,
        "parallel_three_mode_blocks_per_factor":27
      },
      "boundary":"This is an exact support and generator decomposition. The weighted factors do not generally commute, so exp(sum F_j) is not asserted to equal the product of exp(F_j).",
      "parents":["data/w33_pass11054_cubic_jacobian_skew_orthogonal_flow.json","data/w33_e6_cubic_diagonal_phase_weld.json"],
      "checks":{
        "support_is_K81_Cayley":True,"twenty_connection_elements":True,
        "ten_inverse_pairs":True,"all_order_three":True,
        "each_factor_27_disjoint_triangles":True,"each_factor_rank54":True,
        "ten_factors_sum_exactly_to_D_plus":True
      }
    }

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args()
    p=payload();text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or a.output.read_text()!=text:raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"factors":10,"triangles":270},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
