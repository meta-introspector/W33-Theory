#!/usr/bin/env python3
"""Pass 11540 (corrected by Pass 11548): VO(3,3) history theorem.

Pass 11540's graph/classical-group theorem survives unchanged:
  * the 27-event null graph is VO(3,3);
  * |O(3,3)|=48, |SO(3,3)|=24, |Omega(3,3)|=12;
  * affine orthogonal orders are 1296, 648, 324;
  * O/Omega has a C2 x C2 character quotient.

Pass 11548 corrected one order-only inference from the original version:
the repository's PSp Bell stabilizer is NOT the affine SO subgroup.  In the
Pass-11547 Hamming coordinates, O(3,3)=C2^3:S3 and

  repo PSp linear = ker(sigma) = W(D3) ~= S4,
  SO(3,3)          = ker(det)   ~= S4,

where sigma is the product of coordinate signs and det=sigma*pi with pi the
coordinate-permutation parity.  Their intersection is A4=Omega(3,3).

This producer keeps the original VO(3,3) verification executable while
publishing the corrected orientation dictionary from Pass 11548.
"""
from __future__ import annotations

import itertools
import json
from collections import Counter, deque
from pathlib import Path

P=3
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data"/"PART_W33_PASS11540_HISTORY_VO33_ORIENTATION_SQUARE.json"
V=list(itertools.product(range(3),repeat=3))
ZERO=(0,0,0)
I3=(1,0,0,0,1,0,0,0,1)


def mod(x): return x%P
def add(u,v): return tuple(mod(u[i]+v[i]) for i in range(3))
def q(u): return mod(u[0]*u[0]-u[1]*u[1]-u[2]*u[2])
def dot(u,v): return mod(sum(u[i]*v[i] for i in range(3)))
def neg(u): return tuple(mod(-x) for x in u)


def det3(M):
    return mod(
        M[0]*(M[4]*M[8]-M[5]*M[7])
       -M[1]*(M[3]*M[8]-M[5]*M[6])
       +M[2]*(M[3]*M[7]-M[4]*M[6])
    )


def mv(M,v):
    return (
        mod(M[0]*v[0]+M[1]*v[1]+M[2]*v[2]),
        mod(M[3]*v[0]+M[4]*v[1]+M[5]*v[2]),
        mod(M[6]*v[0]+M[7]*v[1]+M[8]*v[2]),
    )


def mm(A,B):
    C=[0]*9
    for i in range(3):
        for j in range(3):
            C[3*i+j]=mod(sum(A[3*i+k]*B[3*k+j] for k in range(3)))
    return tuple(C)


def enumerate_orthogonal_group():
    O=[];SO=[]
    for M in itertools.product(range(3),repeat=9):
        if det3(M)==0: continue
        if all(q(mv(M,v))==q(v) for v in V):
            M=tuple(M);O.append(M)
            if det3(M)==1: SO.append(M)
    return O,SO


def inverse_in_group(M,G):
    for N in G:
        if mm(M,N)==I3 and mm(N,M)==I3:return N
    raise AssertionError


def subgroup_generated(gens):
    H={I3};Q=deque([I3])
    while Q:
        h=Q.popleft()
        for g in gens:
            x=mm(h,g)
            if x not in H:
                H.add(x);Q.append(x)
    return H


def projective_null_pairs(nulls):
    unseen=set(nulls);pairs=[]
    while unseen:
        v=min(unseen);w=neg(v)
        pairs.append(tuple(sorted((v,w))))
        unseen.remove(v);unseen.remove(w)
    return tuple(sorted(pairs))


def pair_perm(M,pairs):
    def idx(v):
        for i,p in enumerate(pairs):
            if v in p:return i
        raise AssertionError
    return tuple(idx(mv(M,p[0])) for p in pairs)


def parity(p):
    inv=sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
    return -1 if inv%2 else 1


def adjacency():
    nulls=[v for v in V if v!=ZERO and q(v)==0]
    A={v:set() for v in V}
    for x in V:
        for n in nulls:A[x].add(add(x,n))
    return A,nulls


def main():
    A,nulls=adjacency()
    assert len(nulls)==8
    assert {len(A[v]) for v in V}=={8}
    assert sum(map(len,A.values()))//2==108

    eig=Counter()
    for k in V:
        c=Counter(dot(k,n) for n in nulls)
        assert c[1]==c[2]
        eig[c[0]-c[1]]+=1
    assert eig==Counter({2:12,-1:8,-4:6,8:1})

    O,SO=enumerate_orthogonal_group()
    assert len(O)==48 and len(SO)==24
    pairs=projective_null_pairs(nulls)
    assert len({pair_perm(M,pairs) for M in SO})==24

    even={M for M in SO if parity(pair_perm(M,pairs))==1}
    inv={M:inverse_in_group(M,SO) for M in SO}
    comm={mm(mm(mm(inv[A],inv[B]),A),B) for A in SO for B in SO}
    Omega=subgroup_generated(comm)
    assert len(Omega)==12 and Omega==even

    correction=json.loads(
        (ROOT/"data"/"PART_W33_PASS11548_HISTORY_ORIENTATION_SQUARE_CORRECTION.json").read_text()
    )
    assert correction["status"]=="PASS_OBJECTWISE_CORRECTION_OF_PASS11540_CHARACTER_ASSIGNMENT"
    assert correction["corrected_subgroups"]["repo_PSp_linear"]["order"]==24
    assert correction["corrected_subgroups"]["repo_PSp_linear"]["not_equal_SO"] is True

    old_big=json.loads(
        (ROOT/"data"/"w33_20260924_history_bigcell_q43_compactification.json").read_text()
    )
    assert old_big["stabilizer_actions"]["PGSp_bell_line_action_order"]==1296
    assert old_big["stabilizer_actions"]["PSp_bell_line_action_order"]==648

    out={
      "schema":"w33.pass11540.history_vo33_orientation_square.v2",
      "status":"PASS_EXACT_VO33_WITH_CHARACTER_ASSIGNMENT_CORRECTED_BY_11548",
      "pass":11540,
      "correction_notice":{
        "corrected_by_pass":11548,
        "withdrawn":"repo PSp = affine SO and global history bit = orthogonal determinant",
        "replacement":"repo PSp = 3^3:ker(sigma)=3^3:W(D3); global history bit is sigma; det=sigma*pi",
      },
      "history_graph":{
        "standard_name":"parabolic affine orthogonal polar graph VO(3,3)",
        "carrier":"F3^3",
        "quadratic_form":"q(t,x,y)=t^2-x^2-y^2",
        "adjacency":"x~y iff x!=y and q(x-y)=0",
        "vertices":27,"degree":8,"edges":108,
        "adjacency_spectrum":{"8":1,"2":12,"-1":8,"-4":6},
      },
      "orthogonal_group":{
        "O_3_3_order":48,
        "SO_3_3_order":24,
        "Omega_3_3_order":12,
        "SO_on_four_null_directions":"faithful S4",
        "Omega_on_four_null_directions":"A4 = even permutations",
        "Omega_equals_derived_SO":True,
        "standard_isomorphisms":[
          "SO(3,3) ~= PGL(2,3) ~= S4",
          "Omega(3,3) ~= PSL(2,3) ~= A4",
        ],
      },
      "affine_orthogonal_ladder":{
        "AffO_order":1296,
        "AffSO_order":648,
        "AffOmega_order":324,
        "warning":"AffSO has order 648 but is not the repository PSp Bell stabilizer",
      },
      "repo_stabilizer_ladder":{
        "PGSp_order":1296,
        "PGSp_linear":"O(3,3)=C2^3:S3",
        "PSp_order":648,
        "PSp_linear":"ker(sigma)=W(D3)~=S4, distinct from SO(3,3)",
        "oriented_kernel_order":324,
        "oriented_linear_kernel":"A4=Omega(3,3)=PSp_linear intersect SO(3,3)",
      },
      "orientation_square_corrected":{
        "sigma":"product of Hamming-coordinate signs; kernel is repo PSp linear; global history/outer bit",
        "pi":"coordinate-permutation parity = projective-null-line permutation parity; repo 648->324 bit",
        "determinant":"sigma*pi; kernel is SO(3,3), not repo PSp",
        "quotient":"O(3,3)/A4 ~= C2 x C2",
        "source":"Pass 11548",
      },
      "projective_null_boundary":{
        "direction_count":4,
        "projective_null_directions":[[list(v) for v in p] for p in pairs],
      },
      "repo_weld":{
        "existing_bigcell_orders_recovered":{"PGSp":1296,"PSp":648},
        "pass11548_objectwise_character_correction_loaded":True,
      },
      "boundary":"Exact finite VO(3,3) theorem with corrected subgroup dictionary. No continuum Lorentz, physical CPT, speed-of-light, Einstein-dynamics or gravity claim.",
    }
    OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
      "status":out["status"],
      "orders":[48,24,12,1296,648,324],
      "repo_PSp":"ker(sigma), not SO",
      "correction":11548,
    },indent=2,sort_keys=True))


if __name__=="__main__":
    main()
