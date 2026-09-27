#!/usr/bin/env python3
"""Pass 11076: explicit chain map from the 40-chart nerve to the W33 Levi building."""
from __future__ import annotations
import argparse,itertools,json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/w33_pass11076_atlas_nerve_steinberg_chain_map.json"

def canon(v):
    v=tuple(int(x)%3 for x in v)
    for x in v:
        if x:
            z=pow(x,-1,3);return tuple((z*y)%3 for y in v)
    raise ValueError("zero")
def symp(x,y):
    return (x[0]*y[2]-x[2]*y[0]+x[1]*y[3]-x[3]*y[1])%3
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
def nullspace(A,p):
    M=np.array(A,dtype=np.int64)%p;n,m=M.shape;r=0;piv=[]
    for c in range(m):
        z=next((i for i in range(r,n) if M[i,c]%p),None)
        if z is None:continue
        M[[r,z]]=M[[z,r]]
        iv=pow(int(M[r,c]),-1,p);M[r]=(M[r]*iv)%p
        for i in range(n):
            if i!=r and M[i,c]%p:
                q=int(M[i,c]);M[i]=(M[i]-q*M[r])%p
        piv.append(c);r+=1
        if r==n:break
    free=[c for c in range(m) if c not in piv]
    B=np.zeros((m,len(free)),dtype=np.int64)
    for j,f in enumerate(free):
        B[f,j]=1
        for i,c in enumerate(piv):B[c,j]=(-M[i,f])%p
    return B
def geometry():
    pts=sorted({canon(v) for v in itertools.product(range(3),repeat=4) if any(v)})
    adj=[[False]*40 for _ in range(40)]
    for i,j in itertools.combinations(range(40),2):
        if symp(pts[i],pts[j])==0:adj[i][j]=adj[j][i]=True
    lines=[frozenset(c) for c in itertools.combinations(range(40),4)
           if all(adj[i][j] for i,j in itertools.combinations(c,2))]
    assert len(pts)==len(lines)==40
    p2l={p:[] for p in range(40)}
    for li,L in enumerate(lines):
        for p in L:p2l[p].append(li)
    assert set(map(len,p2l.values()))=={4}
    return pts,lines,p2l
def payload():
    pts,lines,p2l=geometry()
    edges=set();tris=set();edge_point={}
    for p,ls0 in p2l.items():
        ls=tuple(sorted(ls0))
        for a,b in itertools.combinations(ls,2):
            e=(a,b);edges.add(e)
            assert e not in edge_point or edge_point[e]==p
            edge_point[e]=p
        tris.update(itertools.combinations(ls,3))
    edges=sorted(edges);tris=sorted(tris)
    assert (len(edges),len(tris))==(240,160)
    ei={e:i for i,e in enumerate(edges)}
    D1=np.zeros((40,240),dtype=np.int64)
    for j,(a,b) in enumerate(edges):D1[a,j]=-1;D1[b,j]=1
    D2=np.zeros((240,160),dtype=np.int64)
    for j,(a,b,c) in enumerate(tris):
        for e,s in [((b,c),1),((a,c),-1),((a,b),1)]:
            D2[ei[tuple(sorted(e))],j]+=s
    flags=[(p,l) for p in range(40) for l in p2l[p]]
    fi={f:i for i,f in enumerate(flags)}
    DL=np.zeros((80,160),dtype=np.int64)
    for j,(p,l) in enumerate(flags):DL[p,j]=-1;DL[40+l,j]=1
    T=np.zeros((160,240),dtype=np.int64)
    for j,(a,b) in enumerate(edges):
        p=edge_point[(a,b)]
        T[fi[(p,a)],j]=-1;T[fi[(p,b)],j]=1
    S=np.zeros((80,40),dtype=np.int64)
    for l in range(40):S[40+l,l]=1
    assert np.array_equal(DL@T,S@D1)
    assert not np.any(T@D2)
    replay={}
    for p in (2,3,5,7,101):
        r1,r2,rL=rank_mod(D1,p),rank_mod(D2,p),rank_mod(DL,p)
        Z=nullspace(D1,p);ri=rank_mod((T@Z)%p,p)
        assert (r1,r2,rL,ri)==(39,120,79,81)
        assert Z.shape==(240,201)
        replay[str(p)]={"nerve_D1_rank":39,"nerve_D2_rank":120,"levi_boundary_rank":79,
                        "nerve_cycle_dimension":201,"induced_H1_rank":81,"levi_H1_dimension":81}
    return {
      "schema":"w33.pass11076.atlas-nerve-steinberg-chain-map.v1",
      "status":"PASS_CANONICAL_CHART_NERVE_TO_LEVI_CHAIN_MAP_IS_AN_H1_ISOMORPHISM",
      "headline":"Each nerve edge joins two line charts meeting at a unique W33 point p. Map the oriented chart edge La->Lb to the two-step Levi path La->p->Lb, i.e. -e_(p,La)+e_(p,Lb). Triangle boundaries collapse to zero. The induced map on H1 has rank 81 over every tested field and is an isomorphism onto Levi H1.",
      "chain_map":{
        "nerve_C0":40,"nerve_C1":240,"nerve_C2":160,
        "levi_C0":80,"levi_C1_flags":160,
        "edge_formula":"[La,Lb] -> -[p,La] + [p,Lb], where p=La intersect Lb",
        "boundary_identity":"D_Levi T = S D_nerve",
        "triangle_identity":"T D2 = 0",
        "matrix_nonzeros":480
      },
      "homology":{"nerve_betti1":81,"levi_betti1":81,"induced_map_rank":81,
                  "kernel_on_nerve_cycles_dimension":120,"nerve_boundary_dimension":120,
                  "isomorphism":True},
      "equivariance":"The map is canonical from incidence: every PSp(4,3) automorphism preserves the unique shared point p, so T commutes with the group action.",
      "field_replays":replay,
      "boundary":"This identifies the chart-overlap cycle module with building H1. It does not yet provide Lie-algebra transition maps between local E8 charts.",
      "parents":["data/w33_pass11074_local_e8_chart_overlap_nerve.json","data/PART_W33_20260901_DOUBLE_BUILDING_STEINBERG_HOMOLOGY.json"],
      "checks":{"chain_boundary_identity":True,"triangle_boundaries_killed":True,"H1_rank81_all_fields":True,"canonical_PSp_equivariance":True}
    }
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args();p=payload()
    text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or json.loads(a.output.read_text())!=p:raise SystemExit("certificate drift")
    else:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"H1":81},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
