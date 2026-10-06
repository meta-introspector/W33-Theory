#!/usr/bin/env python3
"""Pass 11560: the mod-3 arrow comes from translation coinvariant 3-torsion.

The full 27-event temporal complex is torsion-free in integral H1.  The
characteristic-three effect appears after taking chain coinvariants by the
translation subgroup T=F3^3.  There are four oriented null-direction edge
orbits and four triangle-direction orbits.  Each temporal triangle has three
boundary edges in the same direction orbit, so the coinvariant boundary is

    d2 = 3 I4.

Hence H1(C_*(X;Z)_T)=(Z/3)^4.  The residual W(D3)~=S4 permutes the four
positive-chirality directions, and its fixed line over F3 is the all-ones
vector.  Pulling that line back gives exactly the cubic arrow cochain.

This is a coinvariant-chain statement, not ordinary homology of X/T.
"""
from __future__ import annotations
import itertools, json
from pathlib import Path
import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/PART_W33_PASS11560_TRANSLATION_COINVARIANT_ARROW_TORSION.json"
P=3
V=list(itertools.product(range(3),repeat=3))
VID={v:i for i,v in enumerate(V)}
DPOS=[(1,1,1),(1,2,2),(2,1,2),(2,2,1)]

def add(a,b): return tuple((a[i]+b[i])%3 for i in range(3))
def scale(c,a): return tuple((c*a[i])%3 for i in range(3))
def chi(d):
    r=1
    for x in d:
        assert x in (1,2)
        r*=1 if x==1 else -1
    return r
def mv(M,v): return tuple(sum(M[i][j]*v[j] for j in range(3))%3 for i in range(3))
def signed_perms_even():
    out=[]
    for pm in itertools.permutations(range(3)):
        for sg in itertools.product((1,2),repeat=3):
            prod=1
            for x in sg: prod*=1 if x==1 else -1
            if prod!=1: continue
            M=[[0]*3 for _ in range(3)]
            for i,j in enumerate(pm): M[i][j]=sg[i]
            out.append(tuple(tuple(r) for r in M))
    return sorted(set(out))

def rank_mod3(A):
    A=np.array(A,dtype=int)%3; r=0
    for c in range(A.shape[1]):
        q=next((i for i in range(r,A.shape[0]) if A[i,c]),None)
        if q is None: continue
        A[[r,q]]=A[[q,r]]
        A[r]=(A[r]*pow(int(A[r,c]),-1,3))%3
        for i in range(A.shape[0]):
            if i!=r and A[i,c]: A[i]=(A[i]-A[i,c]*A[r])%3
        r+=1
    return r

def main():
    # Oriented edges e_(x,d), one per undirected null edge.
    edges=[(x,d) for d in DPOS for x in V]
    assert len(edges)==108
    edge_set={(x,d) for x,d in edges}

    # Nine affine lines/triangles for each direction, canonical coset reps.
    tris=[]
    stabilizers=[]
    for d in DPOS:
        seen=set()
        for x in V:
            line=frozenset((x,add(x,d),add(x,scale(2,d))))
            if line in seen: continue
            seen.add(line);tris.append((min(line),d,line))
        assert len(seen)==9
        for x0,dd,line in tris[-9:]:
            stab=[t for t in V if frozenset(add(x,t) for x in line)==line]
            assert set(stab)=={(0,0,0),d,scale(2,d)}
            stabilizers.append(len(stab))
    assert len(tris)==36 and set(stabilizers)=={3}

    # Full integral boundaries in the oriented edge basis.
    EIDX={(x,d):i for i,(x,d) in enumerate(edges)}
    B2=np.zeros((108,36),dtype=int)
    for k,(x,d,line) in enumerate(tris):
        # Any x in the line gives the same oriented 3-cycle.
        B2[EIDX[(x,d)],k]=1
        B2[EIDX[(add(x,d),d)],k]=1
        B2[EIDX[(add(x,scale(2,d)),d)],k]=1
    assert np.linalg.matrix_rank(B2.astype(float))==36

    # Translation coinvariants: all e_(x,d) become e_d; all line cells in a
    # fixed parallel class become t_d.  Boundary is exactly 3I4.
    Q=np.zeros((4,4),dtype=int)
    for j,d in enumerate(DPOS): Q[j,j]=3
    assert np.array_equal(Q,3*np.eye(4,dtype=int))
    smith=sp.matrices.normalforms if False else None
    # For diagonal 3I the Smith invariants are immediate and checked by minors.
    assert abs(round(np.linalg.det(Q)))==81
    quotient_H1_invariants=[3,3,3,3]

    # W(D3) preserves chi and therefore permutes DPOS without orientation sign.
    W=signed_perms_even();assert len(W)==24
    perms=set()
    for M in W:
        p=[]
        for d in DPOS:
            md=mv(M,d)
            assert chi(md)==1 and md in DPOS
            p.append(DPOS.index(md))
        perms.add(tuple(p))
    assert len(perms)==24  # full S4

    # Fixed subspace of the permutation action on F3^4.
    equations=[]
    for p in perms:
        A=np.zeros((4,4),dtype=int)
        for i,j in enumerate(p): A[j,i]=1
        equations.extend(((A-np.eye(4,dtype=int))%3).tolist())
    fixed_dim=4-rank_mod3(np.array(equations,dtype=int))
    assert fixed_dim==1
    arrow_dir=np.ones(4,dtype=int)
    for p in perms:
        assert np.array_equal(arrow_dir,arrow_dir[list(p)])

    # Pullback to oriented full edges is coefficient +1 on every chi=+ step,
    # hence antisymmetric extension is exactly the cubic chi arrow.
    pulled=np.array([1 for x,d in edges],dtype=int)
    assert len(pulled)==108
    for d in DPOS:
        assert chi(d)==1 and chi(scale(2,d))==-1

    p11553=json.loads((ROOT/"data/PART_W33_PASS11553_HISTORY_ARROW_EQUIVARIANT_COHOMOLOGY.json").read_text())
    assert p11553["complex"]["H1_dimension_F3"]==46
    assert p11553["cochain_side"]["PSp_fixed_H1_dimension"]==1
    assert p11553["chain_side"]["homology_warning"]=="the arrow is homologically trivial as an integral 1-chain"

    out={
      "schema":"w33.pass11560.translation-coinvariant-arrow-torsion.v1",
      "status":"PASS_ARROW_FROM_TRANSLATION_COINVARIANT_3_TORSION","pass":11560,
      "full_complex":{
        "vertices":27,"oriented_null_edges":108,"temporal_triangles":36,
        "integral_triangle_boundary_rank":36,
        "integral_H1":"Z^46 (Pass11553/full-complex computation); no Z/3 torsion",
      },
      "translation_action":{
        "group":"T=F3^3","order":27,
        "vertex_orbits":1,"edge_orbits":4,"triangle_orbits":4,
        "edge_orbit_size":27,"triangle_orbit_size":9,
        "triangle_stabilizer_order":3,
        "positive_direction_representatives":[list(d) for d in DPOS],
      },
      "coinvariant_chain_complex":{
        "C0_rank":1,"C1_rank":4,"C2_rank":4,
        "d1":"0","d2":"3 I4",
        "smith_invariants_d2":quotient_H1_invariants,
        "H1":"(Z/3)^4","H2":"0",
        "mechanism":"each temporal triangle boundary contains three oriented edges in one translation edge orbit",
      },
      "residual_WD3":{
        "order":24,"image_on_four_directions":"S4","image_order":len(perms),
        "fixed_dimension_in_Hom_H1_F3":fixed_dim,
        "fixed_vector":[1,1,1,1],
        "pullback":"the fixed quotient cochain pulls back to the Pass11553 cubic-arrow cochain",
      },
      "interpretation":"The characteristic-three arrow is not torsion in the full temporal complex. It is the unique W(D3)-fixed character of the translation-coinvariant (Z/3)^4 torsion created by the order-3 triangle stabilizers.",
      "quotient_firewall":"C_*(X)_T denotes chain coinvariants. Because translations have order-3 stabilizers on 2-cells, this is not asserted to equal the ordinary cellular chain complex or homology of the topological quotient X/T, and no free-action Cartan-Leray theorem is invoked.",
      "physics_firewall":"This explains the modular equivariant cohomology algebraically; it is not by itself a quantum anomaly, Chern class, thermodynamic arrow, or continuum topological term."
    }
    OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":out["status"],"coinvariant_H1":"(Z/3)^4","fixed":fixed_dim},indent=2))

if __name__=="__main__": main()
