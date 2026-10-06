#!/usr/bin/env python3
"""Pass 11556: the actual PSp history action is a directed fission of H(3,3).

In Pass 11547 the 27 histories became F3^3 and the four Hamming relations.
Pass 11548 identified the actual PSp linear group with W(D3), the even signed
permutations.  Its orbits on nonzero differences have sizes 6,12,4,4:
Hamming weights 1,2 and the two chi tetrahedra inside weight 3.

This producer proves that those five orbitals form a commutative non-symmetric
Schurian translation association scheme.  Its symmetrisation is H(3,3), while
the antisymmetric difference A_+ - A_- is exactly the Pass11555 arrow operator
D.  The dual four-dimensional pair is the old +/-i Weil pair.

Abstract non-symmetric fissions of the 27-point Hamming parent are prior art:
Chia--Kok, Discrete Mathematics 306 (2006), Table 2 / case (3,3).
The repository-specific content is the objectwise PSp/arrow/Weil realization.
"""
from __future__ import annotations

import itertools
import json
from collections import Counter, deque
from pathlib import Path

import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data"/"PART_W33_PASS11556_PSP_HAMMING_DIRECTED_FISSION_SCHEME.json"
P=3
V=list(itertools.product(range(3),repeat=3))
VID={v:i for i,v in enumerate(V)}
N=len(V)


def mod(x): return int(x)%3

def sub(a,b): return tuple(mod(a[i]-b[i]) for i in range(3))

def wt(d): return sum(x!=0 for x in d)

def chi(d):
    assert all(x in (1,2) for x in d)
    z=1
    for x in d: z*=1 if x==1 else -1
    return z

def relation(d):
    w=wt(d)
    if w==0: return 0
    if w==1: return 1
    if w==2: return 2
    return 3 if chi(d)==1 else 4

def mv(M,v):
    return tuple(mod(sum(M[i][j]*v[j] for j in range(3))) for i in range(3))

def signed_permutation_group():
    out=[]
    for pm in itertools.permutations(range(3)):
        for sg in itertools.product((1,2),repeat=3):
            M=[[0]*3 for _ in range(3)]
            for i,j in enumerate(pm): M[i][j]=sg[i]
            out.append(tuple(tuple(r) for r in M))
    return sorted(set(out))

def sigma(M):
    z=1
    for row in M:
        x=next(x for x in row if x)
        z*=1 if x==1 else -1
    return z

def affine_perm(M,t):
    return tuple(VID[tuple(mod(a+b) for a,b in zip(mv(M,v),t))] for v in V)

def relation_action(g, mats):
    # Which relation matrix does conjugation by permutation g send each A_i to?
    out=[]
    for i,A in enumerate(mats):
        B=np.zeros_like(A)
        for x in range(N):
            for y in range(N):
                if A[x,y]: B[g[x],g[y]]=1
        matches=[j for j,C in enumerate(mats) if np.array_equal(B,C)]
        assert len(matches)==1
        out.append(matches[0])
    return tuple(out)

def connected_directed(A):
    # Strong connectivity via forward and reverse BFS.
    def reach(M):
        seen={0};q=deque([0])
        while q:
            x=q.popleft()
            for y in np.flatnonzero(M[x]):
                y=int(y)
                if y not in seen:seen.add(y);q.append(y)
        return len(seen)==N
    return reach(A) and reach(A.T)

def main():
    mats=[np.zeros((N,N),dtype=np.int64) for _ in range(5)]
    for i,x in enumerate(V):
        for j,y in enumerate(V):
            mats[relation(sub(y,x))][i,j]=1

    assert np.array_equal(sum(mats),np.ones((N,N),dtype=np.int64))
    assert np.array_equal(mats[0],np.eye(N,dtype=np.int64))
    assert np.array_equal(mats[1].T,mats[1])
    assert np.array_equal(mats[2].T,mats[2])
    assert np.array_equal(mats[3].T,mats[4])
    valencies=[int(A[0].sum()) for A in mats]
    assert valencies==[1,6,12,4,4]

    # Full exact association-scheme closure.
    pint=np.zeros((5,5,5),dtype=np.int64)
    for i in range(5):
        for j in range(5):
            M=mats[i]@mats[j]
            for k in range(5):
                vals=np.unique(M[mats[k].astype(bool)])
                assert len(vals)==1
                pint[i,j,k]=int(vals[0])
            assert np.array_equal(M,sum(int(pint[i,j,k])*mats[k] for k in range(5)))
    assert all(np.array_equal(mats[i]@mats[j],mats[j]@mats[i])
               for i in range(5) for j in range(5))

    # The actual W(D3) base stabilizer has precisely these difference orbits.
    O=signed_permutation_group()
    W=[M for M in O if sigma(M)==1]
    assert len(O)==48 and len(W)==24
    classes=[{d for d in V if relation(d)==r} for r in range(5)]
    base_orbits=[]
    unseen=set(V)
    while unseen:
        d=min(unseen);orb={mv(M,d) for M in W}
        base_orbits.append(orb);unseen-=orb
    assert sorted(base_orbits,key=lambda s:(len(s),min(s))) != []  # nonvacuity
    assert {frozenset(x) for x in base_orbits}=={frozenset(x) for x in classes}

    # Affine PSp orbitals are exactly the five relations.
    psp={affine_perm(M,t) for M in W for t in V}
    pgsp={affine_perm(M,t) for M in O for t in V}
    assert len(psp)==648 and len(pgsp)==1296
    assert all(relation_action(g,mats)==(0,1,2,3,4) for g in psp)
    outer_actions=Counter(relation_action(g,mats) for g in pgsp-psp)
    assert outer_actions==Counter({(0,1,2,4,3):648})

    # Exact first eigenmatrix from Fourier character sums.
    r3=sp.sqrt(-3)
    E=sp.Matrix([
      [1, 6,12,4,4],
      [1, 3, 0,-2,-2],
      [1, 0,-3,1,1],
      [1,-3, 3,(-1-3*r3)/2,(-1+3*r3)/2],
      [1,-3, 3,(-1+3*r3)/2,(-1-3*r3)/2],
    ])
    C=sp.eye(5); C[3,3]=0;C[4,4]=0;C[3,4]=1;C[4,3]=1
    assert sp.simplify(E*E-27*C)==sp.zeros(5)
    Q=sp.simplify(27*E.inv())
    assert sp.simplify(Q-E*C)==sp.zeros(5)

    # Verify all 27 character sums against E row classes.
    omega=sp.exp(2*sp.pi*sp.I/3)
    dual_class_counts=Counter()
    for k in V:
        row=[]
        for r in range(5):
            z=0
            for d in classes[r]:
                phase=sum(k[i]*d[i] for i in range(3))%3
                z+=omega**phase
            z=sp.simplify(sp.expand_complex(z))
            row.append(z)
        matches=[a for a in range(5)
                 if all(sp.simplify(row[b]-E[a,b])==0 for b in range(5))]
        assert len(matches)==1
        dual_class_counts[matches[0]]+=1
    assert dual_class_counts==Counter({0:1,1:6,2:12,3:4,4:4})

    # Hamming fusion and arrow/Weil fission.
    A3=mats[3]+mats[4]
    D=mats[3]-mats[4]
    assert np.array_equal(A3.T,A3) and np.all(A3.sum(axis=1)==8)
    p11555=json.loads((ROOT/"data/PART_W33_PASS11555_TEMPORAL_LOCAL_DYNAMICS_DIRAC.json").read_text())
    assert p11555["operator_D"]["exact_square_identity"]=="3 D^2 = L(L-6I)(L-12I)"
    Dold=np.zeros((N,N),dtype=np.int64)
    for i,x in enumerate(V):
        for j,y in enumerate(V):
            if i!=j and wt(sub(y,x))==3: Dold[i,j]=chi(sub(y,x))
    assert np.array_equal(D,Dold)

    # Full Hamming automorphism group is the affine signed-permutation group,
    # already certified in Pass11547.  Therefore the color-preserving subgroup
    # of the fission has exactly 648 elements: the other 648 swap R+ and R-.
    p11547=json.loads((ROOT/"data/PART_W33_PASS11547_HISTORY_HAMMING_CANONICAL_CLOCK.json").read_text())
    p11548=json.loads((ROOT/"data/PART_W33_PASS11548_HISTORY_ORIENTATION_SQUARE_CORRECTION.json").read_text())
    assert p11547["automorphism_group"]["order"]==1296
    assert p11548["corrected_subgroups"]["repo_PSp_linear"]["order"]==24
    assert len(pgsp)==p11547["automorphism_group"]["order"]

    # Each nontrivial relation is connected/strongly connected: primitive scheme.
    connectivity=[connected_directed(A) for A in mats[1:]]
    assert connectivity==[True,True,True,True]

    out={
      "schema":"w33.pass11556.psp_hamming_directed_fission_scheme.v1",
      "status":"PASS_PSP_ORBITAL_FISSION_OF_H33","pass":11556,
      "scheme":{
        "carrier":"F3^3","order":27,"classes":4,"rank":5,
        "type":"primitive commutative non-symmetric Schurian translation association scheme",
        "relations":{
          "R0":"difference Hamming weight 0",
          "R1":"difference Hamming weight 1",
          "R2":"difference Hamming weight 2",
          "Rplus":"weight 3 and chi(d)=+1",
          "Rminus":"weight 3 and chi(d)=-1",
        },
        "valencies":valencies,
        "transpose_map":[0,1,2,4,3],
        "all_nontrivial_relations_connected":connectivity,
        "intersection_numbers_pijk":pint.tolist(),
      },
      "schurian_orbitals":{
        "group":"3^3:W(D3) = affine repo PSp history stabilizer",
        "order":648,
        "base_stabilizer_order":24,
        "base_difference_orbit_sizes":sorted(map(len,base_orbits)),
        "orbitals_equal_five_relations":True,
        "color_preserving_automorphism_group_order":648,
      },
      "hamming_fusion":{
        "fusion":"R3 = Rplus union Rminus",
        "result":"ternary Hamming scheme H(3,3)",
        "full_Hamming_automorphism_order":1296,
        "outer_coset_size":648,
        "outer_relation_action":"fix R0,R1,R2 and swap Rplus <-> Rminus",
      },
      "eigenstructure":{
        "first_eigenmatrix":[[str(sp.simplify(E[i,j])) for j in range(5)] for i in range(5)],
        "multiplicities":[1,6,12,4,4],
        "second_eigenmatrix":[[str(sp.simplify(Q[i,j])) for j in range(5)] for i in range(5)],
        "formal_self_duality":"P^2=27*C where C swaps the two chiral classes; Q=P*C",
        "splitting_field":"Q(sqrt(-3))",
      },
      "arrow_weil_weld":{
        "null_adjacency":"A3=Aplus+Aminus",
        "arrow_operator":"D=Aplus-Aminus",
        "D_is_exact_Pass11555_operator":True,
        "dual_chiral_multiplicities":[4,4],
        "reading":"the two four-dimensional primitive character classes are the chi=+/- weight-3 Fourier sheets; Pass11550 identifies them with the conjugate +/-i Weil sheets",
      },
      "prior_art":{
        "abstract_fission":"Chia and Kok, Discrete Mathematics 306 (2006) 3189-3222, Table 2 includes the 27-point Hamming spectra and a feasible class-4 non-symmetric case (3,3)",
        "classification":"Hanaki-Miyamoto classify 502 association schemes of order 27; no abstract novelty claim is made here",
        "repo_specific_new_content":"objectwise identification of the actual PSp history orbitals with the H(3,3) directed fission and exact weld D=Aplus-Aminus to the temporal arrow/Weil operator",
      },
      "boundary":"Exact finite orbital/association-scheme theorem. It does not imply continuum causality, thermodynamic irreversibility, a physical Dirac equation, or gravitational dynamics.",
    }
    OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":out["status"],"valencies":valencies,
                      "aut":648,"outer_swap":648,"dual":dict(dual_class_counts)},indent=2))

if __name__=="__main__":
    main()
