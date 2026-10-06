#!/usr/bin/env python3
"""Pass 11553: the history arrow is the unique PSp-fixed mod-3 H^1 class."""
from __future__ import annotations
import itertools, json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/PART_W33_PASS11553_HISTORY_ARROW_EQUIVARIANT_COHOMOLOGY.json"
V=list(itertools.product(range(3),repeat=3)); VID={v:i for i,v in enumerate(V)}
def sub(a,b): return tuple((a[i]-b[i])%3 for i in range(3))
def chi(d):
    r=1
    for x in d: r*=1 if x==1 else -1
    return r

E=[(i,j) for i in range(27) for j in range(i+1,27) if all(sub(V[j],V[i]))]
EI={e:k for k,e in enumerate(E)}
T=[t for t in itertools.combinations(range(27),3)
   if all(tuple(sorted(e)) in EI for e in itertools.combinations(t,2))]
B1=np.zeros((27,len(E)),dtype=int)
B2=np.zeros((len(E),len(T)),dtype=int)
for k,(i,j) in enumerate(E): B1[i,k]=-1; B1[j,k]=1
for k,(a,b,c) in enumerate(T):
    for x,y,sgn in [(b,c,1),(a,c,-1),(a,b,1)]: B2[EI[(x,y)],k]=sgn

def rref(A,p=3):
    A=np.array(A,dtype=int)%p; m,n=A.shape; r=0; piv=[]
    for c in range(n):
        q=next((i for i in range(r,m) if A[i,c]%p),None)
        if q is None: continue
        A[[r,q]]=A[[q,r]]
        A[r]=(A[r]*pow(int(A[r,c]),-1,p))%p
        for i in range(m):
            if i!=r and A[i,c]%p: A[i]=(A[i]-A[i,c]*A[r])%p
        piv.append(c); r+=1
    return A,piv

def rank(A): return len(rref(A)[1])
def nullbasis(A):
    R,piv=rref(A); free=[j for j in range(A.shape[1]) if j not in piv]; out=[]
    for f in free:
        x=np.zeros(A.shape[1],dtype=int); x[f]=1
        for i,p in enumerate(piv): x[p]=(-R[i,f])%3
        out.append(x)
    return np.array(out,dtype=int).T if out else np.zeros((A.shape[1],0),dtype=int)
def colbasis(A):
    _,piv=rref(A)
    return np.array(A,dtype=int)[:,piv]
def invmod(A):
    A=np.array(A,dtype=int)%3; n=A.shape[0]
    X=np.c_[A,np.eye(n,dtype=int)]; r=0
    for c in range(n):
        q=next(i for i in range(r,n) if X[i,c]%3)
        X[[r,q]]=X[[q,r]]
        X[r]=(X[r]*pow(int(X[r,c]),-1,3))%3
        for i in range(n):
            if i!=r and X[i,c]%3: X[i]=(X[i]-X[i,c]*X[r])%3
        r+=1
    return X[:,n:]%3
def edge_action(perm,x):
    out=np.zeros(len(E),dtype=int)
    for k,(a,b) in enumerate(E):
        u,v=perm[a],perm[b]
        if u<v: out[EI[(u,v)]]+=x[k]
        else: out[EI[(v,u)]]-=x[k]
    return out%3
def affine_perm(M,t):
    def mv(v): return tuple(sum(M[i][j]*v[j] for j in range(3))%3 for i in range(3))
    return tuple(VID[tuple((a+b)%3 for a,b in zip(mv(v),t))] for v in V)

def main():
    assert len(E)==108 and len(T)==36 and np.array_equal(B1@B2,np.zeros((27,36),dtype=int))
    arrow=np.array([chi(sub(V[j],V[i])) for i,j in E],dtype=int)
    assert np.all(B1@arrow==0)
    tri_coeff=np.zeros(36,dtype=int)
    for t in range(36):
        es=np.flatnonzero(B2[:,t])
        vals={int(B2[e,t]*arrow[e]) for e in es}
        assert len(vals)==1
        tri_coeff[t]=vals.pop()
    assert np.array_equal(B2@tri_coeff,arrow)
    tri_integrals=B2.T@arrow
    assert set(map(int,tri_integrals))=={-3,3}
    assert np.all(tri_integrals%3==0)

    Z=nullbasis(B2.T); B=colbasis(B1.T%3)
    assert Z.shape==(108,72) and B.shape==(108,26)
    cols=[B[:,i] for i in range(B.shape[1])]; H=[]; rr=26
    for j in range(Z.shape[1]):
        cand=np.column_stack(cols+[Z[:,j]]); nr=rank(cand)
        if nr>rr: cols.append(Z[:,j]); H.append(Z[:,j]); rr=nr
        if rr==72: break
    H=np.column_stack(H); Q=np.column_stack([B,H]); assert H.shape==(108,46)
    _,rows=rref(Q.T); rows=rows[:72]
    Minv=invmod(Q[rows,:])
    def coeff(v): return (Minv@(v[rows]%3))%3

    I=((1,0,0),(0,1,0),(0,0,1))
    swap12=((0,1,0),(1,0,0),(0,0,1))
    cycle=((0,1,0),(0,0,1),(1,0,0))
    flip12=((2,0,0),(0,2,0),(0,0,1))
    gens=[affine_perm(I,e) for e in [(1,0,0),(0,1,0),(0,0,1)]]
    gens += [affine_perm(M,(0,0,0)) for M in [swap12,cycle,flip12]]
    acts=[]
    for g in gens:
        A=np.zeros((46,46),dtype=int)
        for j in range(46): A[:,j]=coeff(edge_action(g,H[:,j]))[26:]
        acts.append(A%3)
    fixed=nullbasis(np.vstack([(A-np.eye(46,dtype=int))%3 for A in acts]))
    ac=coeff(arrow%3)[26:]%3
    assert np.any(ac)
    assert all(np.all(((A-np.eye(46,dtype=int))@ac)%3==0) for A in acts)
    assert fixed.shape[1]==1

    out={"schema":"w33.pass11553.history_arrow_equivariant_cohomology.v1",
      "status":"PASS_UNIQUE_PSP_FIXED_MOD3_H1_ARROW_CLASS","pass":11553,
      "complex":{"vertices":27,"edges":108,"triangles":36,"rank_d0_F3":26,"rank_d1_F3":36,
        "H1_dimension_F3":46,"H2_dimension_F3":0},
      "chain_side":{"integral_boundary_zero":True,"is_boundary_of_all_temporal_triangles":True,
        "triangle_chain_coefficient_histogram":{"+1":int(np.sum(tri_coeff==1)),"-1":int(np.sum(tri_coeff==-1))},
        "formula":"arrow is the boundary of a +/-1 sum of all 36 temporal triangles; the coefficient pattern depends on vertex-order gauge",
        "homology_warning":"the arrow is homologically trivial as an integral 1-chain"},
      "cochain_side":{"triangle_integrals_over_Z":sorted(set(map(int,tri_integrals))),
        "cocycle_mod3":True,"exact_mod3":False,"PSp_fixed_H1_dimension":1,"arrow_spans_fixed_line":True,
        "PSp_generators_checked":"3 translations + S3 coordinate generators + one even sign flip"},
      "interpretation":"The same edge field is an integral boundary on the chain side but the unique PSp-invariant nonzero mod-3 H^1 class on the cochain side. Its +/-3 triangle circulation is exactly what disappears modulo 3.",
      "topology_firewall":"H^2 over F3 vanishes for this 2-complex, so this result does not by itself produce a Chern class. A Maslov/transgression interpretation would require additional structure.",
      "boundary":"Exact finite simplicial/cohomological statement. It is not thermodynamic irreversibility, continuum time orientation, or a characteristic-class derivation of physical CPT."}
    OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":out["status"],"H1":46,"fixed":1,"triangle_values":sorted(set(map(int,tri_integrals)))},indent=2))
if __name__=="__main__": main()
