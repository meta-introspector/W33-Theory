#!/usr/bin/env python3
"""Pass 11074: the 40 W33 line-chart overlap nerve has H1=81 and the homotopy type of the Levi graph."""
from __future__ import annotations
import argparse,itertools,json
from collections import Counter
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/w33_pass11074_local_e8_chart_overlap_nerve.json"

def canon(v):
    v=tuple(int(x)%3 for x in v)
    for x in v:
        if x:
            z=pow(x,-1,3)
            return tuple((z*y)%3 for y in v)
    raise ValueError("zero")

def symp(x,y):
    return (x[0]*y[2]-x[2]*y[0]+x[1]*y[3]-x[3]*y[1])%3

def rank_mod(A,p):
    M=np.array(A,dtype=np.int64)%p;n,m=M.shape;r=0
    for c in range(m):
        z=next((i for i in range(r,n) if M[i,c]%p),None)
        if z is None:continue
        M[[r,z]]=M[[z,r]]
        iv=pow(int(M[r,c]%p),-1,p);M[r]=(M[r]*iv)%p
        for i in range(n):
            if i!=r and M[i,c]%p:
                q=int(M[i,c]);M[i]=(M[i]-q*M[r])%p
        r+=1
        if r==n:break
    return r

def geometry():
    pts=sorted({canon(v) for v in itertools.product(range(3),repeat=4) if any(v)})
    assert len(pts)==40
    adj=[[False]*40 for _ in range(40)]
    for i,j in itertools.combinations(range(40),2):
        if symp(pts[i],pts[j])==0:adj[i][j]=adj[j][i]=True
    lines=[frozenset(c) for c in itertools.combinations(range(40),4)
           if all(adj[i][j] for i,j in itertools.combinations(c,2))]
    assert len(lines)==40
    p2l={p:[] for p in range(40)}
    for li,L in enumerate(lines):
        for p in L:p2l[p].append(li)
    assert set(map(len,p2l.values()))=={4}
    return pts,lines,p2l

def payload():
    pts,lines,p2l=geometry()
    edges=set();tris=set();tets=set()
    for p,ls0 in p2l.items():
        ls=tuple(sorted(ls0));tets.add(ls)
        edges.update(tuple(sorted(x)) for x in itertools.combinations(ls,2))
        tris.update(tuple(sorted(x)) for x in itertools.combinations(ls,3))
    edges=sorted(edges);tris=sorted(tris);tets=sorted(tets)
    assert (len(lines),len(edges),len(tris),len(tets))==(40,240,160,40)

    hist=Counter()
    for a,b in itertools.combinations(tets,2):hist[len(set(a)&set(b))]+=1
    assert hist==Counter({0:540,1:240})

    ei={e:i for i,e in enumerate(edges)};ti={t:i for i,t in enumerate(tris)}
    D1=np.zeros((40,240),dtype=np.int64)
    for j,(a,b) in enumerate(edges):D1[a,j]=-1;D1[b,j]=1
    D2=np.zeros((240,160),dtype=np.int64)
    for j,(a,b,c) in enumerate(tris):
        for e,s in [((b,c),1),((a,c),-1),((a,b),1)]:D2[ei[tuple(sorted(e))],j]+=s
    D3=np.zeros((160,40),dtype=np.int64)
    for j,t in enumerate(tets):
        for i in range(4):
            face=tuple(x for k,x in enumerate(t) if k!=i)
            D3[ti[tuple(sorted(face))],j]+=(-1)**i
    assert not np.any(D1@D2) and not np.any(D2@D3)

    replay={}
    for p in (2,3,5,7,101):
        r1,r2,r3=rank_mod(D1,p),rank_mod(D2,p),rank_mod(D3,p)
        betti=[40-r1,240-r1-r2,160-r2-r3,40-r3]
        assert (r1,r2,r3)==(39,120,40) and betti==[1,81,0,0]
        replay[str(p)]={"boundary_ranks":[r1,r2,r3],"betti":[1,81,0,0]}

    degree=Counter()
    for a,b in edges:degree[a]+=1;degree[b]+=1
    assert set(degree.values())=={12}
    tetra_per_chart=Counter(x for T in tets for x in T)
    assert set(tetra_per_chart.values())=={4}

    return {
      "schema":"w33.pass11074.local-e8-chart-overlap-nerve.v1",
      "status":"PASS_40_LINE_CHART_OVERLAP_NERVE_HAS_F_VECTOR_40_240_160_40_AND_EXACT_H1_RANK81",
      "headline":(
        "Treat each W33 line as one local tetracode/E8 chart label, and declare a family of charts "
        "to overlap when their W33 lines share one point/A2 label. The maximal overlaps are the "
        "four line charts through a point, hence 40 tetrahedra. The resulting simplicial nerve has "
        "f-vector (40,240,160,40), boundary ranks (39,120,40), and homology (1,81,0,0)."
      ),
      "nerve":{
        "vertices_line_charts":40,"edges_pairwise_point_overlaps":240,
        "triangles_triple_point_overlaps":160,"tetrahedra_fourfold_point_overlaps":40,
        "f_vector":[40,240,160,40],"euler_characteristic":-80,
        "chart_overlap_degree":12,"point_overlaps_per_chart":4,
        "tetrahedron_pair_intersection_histogram":{"0":540,"1":240}
      },
      "homology":{
        "boundary_ranks":[39,120,40],"betti":[1,81,0,0],
        "H1":"Z^81","higher_reduced_homology":"0",
        "split_prime_and_control_replays":replay
      },
      "homotopy_reading":(
        "The 40 maximal tetrahedra intersect one another only in vertices. After barycentric subdivision, "
        "each tetrahedron collapses relative to its four chart vertices onto the star from its point-barycenter. "
        "Doing this independently for all 40 tetrahedra yields exactly the 80-vertex, 160-edge W33 point-line Levi graph. "
        "Thus the chart-overlap nerve has the Levi graph homotopy type and its cycle rank is 160-80+1=81."
      ),
      "local_e8_boundary":(
        "The nerve records where local line charts would share their point/A2 labels. It does not yet supply "
        "Lie-algebra embedding maps between the E8 charts on those overlaps."
      ),
      "parents":["data/w33_pass11071_oriented_chamber_tetracode_cover.json","analysis/w33_tetracode_e8_root_system_bridge.py"],
      "checks":{"W33_40_points_40_lines":True,"forty_maximal_tetrahedra":True,
                "f_vector_exact":True,"boundary_ranks_39_120_40":True,
                "H1_rank81":True,"higher_homology_zero":True,
                "maximal_simplices_intersect_at_most_in_one_vertex":True}
    }

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args();p=payload()
    text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or json.loads(a.output.read_text())!=p:raise SystemExit("certificate drift")
    else:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"f":[40,240,160,40],"b1":81},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
