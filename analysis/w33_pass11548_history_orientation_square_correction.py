#!/usr/bin/env python3
"""Pass 11548: objectwise correction of Pass 11540's history orientation square.

Pass 11540 got the finite graph and group orders right but identified the
repo's 648-element PSp Bell stabilizer with the wrong index-two subgroup of
O(3,3).  The order-24 linear PSp action must be reconstructed from the
repository's actual GL(2,3) congruence action, not inferred from its order.

In the Pass-11547 Hamming coordinates, O(3,3) is the 48-element signed
permutation group C2^3:S3.  It has two independent sign characters:
  sigma(M) = product of the three coordinate signs,
  pi(M)    = parity of the coordinate permutation.
Their product is the ordinary 3x3 determinant:
  det(M) = sigma(M) pi(M).

Exact corrected dictionary:
  * repo PSp Bell stabilizer linear part = ker(sigma) = W(D3) ~= S4;
  * repo PGSp Bell stabilizer linear part = O(3,3), order 48;
  * SO(3,3) = ker(det), another order-24 S4, NOT the repo PSp subgroup;
  * PSp_linear intersect SO(3,3) = A4 = Omega(3,3), order 12;
  * the 648->324 repo line-orientation bit is pi restricted to ker(sigma);
  * the 1296->648 global history/antiunitary bit is sigma;
  * det = sigma*pi is the third nontrivial character of O/A4 ~= C2 x C2.

The eight oriented null vectors in Hamming coordinates are the cube corners
(+/-1)^3.  ker(sigma) preserves their product, splitting them into two
tetrahedra of four; -I flips the tetrahedra.  This reproduces the repo's
previously frozen 4+4 PSp / fused-8 PGSp result objectwise.

This pass supersedes only the character assignment in Pass 11540.  The
VO(3,3) identification, orders 48/24/12 and affine orders 1296/648/324 remain.
"""
from __future__ import annotations

import itertools
import json
from collections import Counter, deque
from pathlib import Path

P=3
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data"/"PART_W33_PASS11548_HISTORY_ORIENTATION_SQUARE_CORRECTION.json"
I3=((1,0,0),(0,1,0),(0,0,1))
V=list(itertools.product(range(3),repeat=3))
ZERO=(0,0,0)


def mod(x): return x%P


def det2(A):
    return mod(A[0][0]*A[1][1]-A[0][1]*A[1][0])


def det3(M):
    return mod(
        M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
       -M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
       +M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0])
    )


def mm(A,B):
    return tuple(tuple(
        mod(sum(A[i][k]*B[k][j] for k in range(3)))
        for j in range(3)
    ) for i in range(3))


def mv(M,v):
    return tuple(mod(sum(M[i][j]*v[j] for j in range(3))) for i in range(3))


def qh(v):
    return mod(sum(x*x for x in v))


def Hsym(s):
    """Pass-11547 map for q(a,b,c)=ac-b^2."""
    a,b,c=s
    return (mod(2*a+b+c),mod(2*a+2*b+c),mod(2*a+2*c))


HINV={Hsym(v):v for v in V}
assert len(HINV)==27


def sym2_action(A,s):
    a,b,c=s
    S=((a,b),(b,c))
    # A S A^T, expanded with tiny integer loops.
    AS=tuple(tuple(mod(sum(A[i][k]*S[k][j] for k in range(2)))
                   for j in range(2)) for i in range(2))
    R=tuple(tuple(mod(sum(AS[i][k]*A[j][k] for k in range(2)))
                  for j in range(2)) for i in range(2))
    assert R[0][1]==R[1][0]
    return (R[0][0],R[0][1],R[1][1])


def gl2_projective_reps():
    reps={}
    for z in itertools.product(range(3),repeat=4):
        A=(z[:2],z[2:])
        if det2(A)==0: continue
        neg=tuple(mod(-x) for x in z)
        key=min(tuple(z),neg)
        reps[key]=(key[:2],key[2:])
    assert len(reps)==24
    return list(reps.values())


def linear_matrix(f):
    e=((1,0,0),(0,1,0),(0,0,1))
    cols=[f(x) for x in e]
    return tuple(tuple(cols[j][i] for j in range(3)) for i in range(3))


def repo_psp_linear():
    mats=set()
    for A in gl2_projective_reps():
        def f(z,A=A):
            return Hsym(sym2_action(A,HINV[z]))
        M=linear_matrix(f)
        assert all(f(v)==mv(M,v) for v in V)
        mats.add(M)
    assert len(mats)==24
    return mats


def signed_permutation_group():
    out=set()
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1,2),repeat=3):
            M=[[0]*3 for _ in range(3)]
            for i,j in enumerate(perm):
                M[i][j]=signs[i]
            out.add(tuple(tuple(r) for r in M))
    assert len(out)==48
    return out


def parity(p):
    inv=sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
    return -1 if inv%2 else 1


def chars(M):
    perm=[]; signs=[]
    for row in M:
        nz=[j for j,x in enumerate(row) if x]
        assert len(nz)==1
        j=nz[0]
        perm.append(j)
        signs.append(row[j])
    assert sorted(perm)==[0,1,2]
    pi=parity(perm)
    sigma=1
    for s in signs:
        sigma*=1 if s==1 else -1
    det=1 if det3(M)==1 else -1
    assert det==sigma*pi
    return dict(
        coordinate_permutation=tuple(perm),
        coordinate_signs=tuple(signs),
        coordinate_permutation_parity=pi,
        sign_product=sigma,
        determinant=det,
    )


def inverse_in_group(M,G):
    for N in G:
        if mm(M,N)==I3 and mm(N,M)==I3:
            return N
    raise AssertionError("inverse absent")


def subgroup_generated(gens):
    H={I3}; q=deque([I3])
    while q:
        h=q.popleft()
        for g in gens:
            x=mm(h,g)
            if x not in H:
                H.add(x);q.append(x)
    return H


def derived(G):
    inv={A:inverse_in_group(A,G) for A in G}
    comm={
        mm(mm(mm(inv[A],inv[B]),A),B)
        for A in G for B in G
    }
    return subgroup_generated(comm)


def projective_null_pairs():
    full=[v for v in V if all(v)]
    unseen=set(full); pairs=[]
    while unseen:
        v=min(unseen)
        w=tuple(mod(-x) for x in v)
        pair=tuple(sorted((v,w)))
        pairs.append(pair)
        unseen.remove(v);unseen.remove(w)
    return tuple(sorted(pairs))


NULL_PAIRS=projective_null_pairs()


def null_line_perm(M):
    def idx(v):
        for i,pair in enumerate(NULL_PAIRS):
            if v in pair:return i
        raise AssertionError
    return tuple(idx(mv(M,p[0])) for p in NULL_PAIRS)


def oriented_null_sign(v):
    assert all(x in (1,2) for x in v)
    out=1
    for x in v:
        out*=1 if x==1 else -1
    return out


def orbit(G,x):
    return {mv(M,x) for M in G}


def order(M,cap=24):
    X=I3
    for n in range(1,cap+1):
        X=mm(X,M)
        if X==I3:return n
    raise AssertionError("order cap")


def main():
    O=signed_permutation_group()
    PSp=repo_psp_linear()

    # Every reconstructed repo action is an orthogonal signed permutation.
    assert PSp <= O
    assert all(qh(mv(M,v))==qh(v) for M in O for v in V)

    ker_sigma={M for M in O if chars(M)["sign_product"]==1}
    SO={M for M in O if chars(M)["determinant"]==1}
    ker_pi={M for M in O if chars(M)["coordinate_permutation_parity"]==1}
    assert len(ker_sigma)==len(SO)==len(ker_pi)==24

    # The core correction.
    assert PSp==ker_sigma
    assert PSp!=SO

    A4=PSp & SO
    assert len(A4)==12
    assert A4=={M for M in O if chars(M)["sign_product"]==1
                                  and chars(M)["coordinate_permutation_parity"]==1}
    assert derived(PSp)==A4
    assert derived(SO)==A4
    assert Counter(order(M) for M in A4)==Counter({3:8,2:3,1:1})

    # Projective four-null-direction parity is exactly coordinate permutation parity.
    assert all(parity(null_line_perm(M))==chars(M)["coordinate_permutation_parity"]
               for M in O)

    # Character square: all four (sigma,pi) sectors have size 12.
    square=Counter()
    for M in O:
        c=chars(M)
        square[(c["sign_product"],c["coordinate_permutation_parity"])]+=1
    assert set(square.values())=={12}
    assert len(square)==4

    # The oriented null cube splits 4+4 under PSp by product sign; full O fuses.
    plus=orbit(PSp,(1,1,1))
    minus=orbit(PSp,(1,1,2))
    assert len(plus)==len(minus)==4 and plus.isdisjoint(minus)
    assert {oriented_null_sign(v) for v in plus}=={1}
    assert {oriented_null_sign(v) for v in minus}=={-1}
    assert len(orbit(O,(1,1,1)))==8

    minusI=((2,0,0),(0,2,0),(0,0,2))
    assert minusI in O and minusI not in PSp
    assert chars(minusI)=={
        "coordinate_permutation":(0,1,2),
        "coordinate_signs":(2,2,2),
        "coordinate_permutation_parity":1,
        "sign_product":-1,
        "determinant":-1,
    }
    assert {oriented_null_sign(mv(minusI,v)) for v in plus}=={-1}

    # Match frozen repo observations.
    old_big=json.loads((ROOT/"data"/"w33_20260924_history_bigcell_q43_compactification.json").read_text())
    old_or=json.loads((ROOT/"data"/"PART_W33_PASS9741_9748_ORIENTATION_CHARACTER_WELD.json").read_text())
    old_44=json.loads((ROOT/"data"/"w33_20260924_fixed9_chiral_null_fourplusfour.json").read_text())
    assert old_big["stabilizer_actions"]["PSp_bell_line_action_order"]==27*len(PSp)==648
    assert old_big["stabilizer_actions"]["PGSp_bell_line_action_order"]==27*len(O)==1296
    assert old_or["W33_line_orientation"]["order_kernel"]==27*len(A4)==324
    assert old_44["status"]=="PASS_FIXED9_IS_VACUUM_PLUS_TWO_ORIENTED_NULL_FOURS"

    out={
      "schema":"w33.pass11548.history_orientation_square_correction.v1",
      "status":"PASS_OBJECTWISE_CORRECTION_OF_PASS11540_CHARACTER_ASSIGNMENT",
      "pass":11548,
      "supersedes":{
        "pass11540_claim":"repo PSp affine subgroup = 3^3:SO(3,3) and global bit = orthogonal determinant",
        "verdict":"WITHDRAWN",
        "survives":[
          "history graph = VO(3,3)",
          "|O(3,3)|=48, |SO(3,3)|=24, |Omega(3,3)|=12",
          "affine orders 1296,648,324",
          "O/Omega has a C2 x C2 character square",
        ],
      },
      "hamming_orthogonal_group":{
        "O_order":48,
        "model":"signed permutation group C2^3:S3",
        "characters":{
          "sigma":"product of three coordinate signs",
          "pi":"parity of coordinate permutation = parity on four projective null directions",
          "determinant":"sigma*pi",
        },
        "character_sector_counts":{
          f"sigma={s:+d},pi={p:+d}":n for (s,p),n in sorted(square.items())
        },
      },
      "corrected_subgroups":{
        "repo_PSp_linear":{
          "definition":"ker(sigma): even number of coordinate sign flips",
          "order":24,
          "standard_name":"W(D3), the even signed permutation group",
          "isomorphism":"S4",
          "not_equal_SO":True,
        },
        "SO_3_3":{
          "definition":"ker(determinant)=ker(sigma*pi)",
          "order":24,
          "isomorphism":"S4",
          "distinct_from_repo_PSp":True,
        },
        "third_index_two_kernel":{
          "definition":"ker(pi)",
          "order":24,
          "structure":"C2 x A4 (equivalently C2^3:C3)",
        },
        "common_A4_Omega":{
          "definition":"ker(sigma) intersect ker(pi) = repo_PSp intersect SO",
          "order":12,
          "isomorphism":"A4 = Omega(3,3)",
          "derived_of_repo_PSp":True,
          "derived_of_SO":True,
        },
      },
      "corrected_orientation_dictionary":{
        "global_history_outer_bit":{
          "affine_quotient":"(3^3:O)/(3^3:ker(sigma)) = 1296/648 = C2",
          "character":"sigma = product of Hamming-coordinate sign flips",
          "repo_meaning":"existing history-cycle orientation / unitary-vs-antiunitary / chirality outer bit",
          "outer_witness":"-I has sigma=-1 and swaps the two oriented-null tetrahedra",
        },
        "inner_line_orientation_bit":{
          "affine_quotient":"(3^3:ker(sigma))/(3^3:A4) = 648/324 = C2",
          "character":"pi = parity of coordinate permutation = S4/A4 sign",
          "repo_meaning":"existing W33 line-orientation kernel 3^3:A4",
        },
        "orthogonal_determinant":{
          "character":"det=sigma*pi",
          "meaning":"third nontrivial character on O/A4; its kernel SO(3,3) is NOT the repo PSp subgroup",
        },
        "quotient":"O(3,3)/A4 = C2 x C2",
      },
      "oriented_null_tetrahedra":{
        "plus_product":[list(v) for v in sorted(plus)],
        "minus_product":[list(v) for v in sorted(minus)],
        "PSp_orbits":[4,4],
        "PGSp_O_orbit":8,
        "closed_form_label":"chi(d)=product_i sign(d_i) for d_i in {+1,-1}",
      },
      "spinor_norm_firewall":{
        "valid_standard_statement":"Omega(3,3) is the spinor-norm kernel inside SO(3,3)",
        "invalid_old_shortcut":"the repo 648->324 quotient itself is SO->Omega",
        "correct_statement":"repo 648->324 is ker(sigma)->A4 via pi; A4 also equals Omega because it is the intersection with SO",
      },
      "boundary":"Exact finite subgroup correction. No continuum Lorentz, physical CPT theorem, weak-chirality selection or gravity claim is added.",
    }
    OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
      "status":out["status"],
      "PSp":"ker(sigma)",
      "SO":"ker(sigma*pi)",
      "A4":len(A4),
      "orbits":[len(plus),len(minus)],
      "sector_counts":out["hamming_orthogonal_group"]["character_sector_counts"],
    },indent=2,sort_keys=True))


if __name__=="__main__":
    main()
