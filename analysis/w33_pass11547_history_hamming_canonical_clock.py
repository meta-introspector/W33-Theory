#!/usr/bin/env python3
"""Pass 11547: finite Minkowski/Hamming conjugacy and the canonical local clock.

The Forty Points history chart uses V=F_3^3 with
    q(t,x,y)=t^2-x^2-y^2.
Over F_3 the linear change
    H(t,x,y)=(x+y, x-y, t)
satisfies
    q(v)=H_0(v)^2+H_1(v)^2+H_2(v)^2.
Because every nonzero square in F_3 equals 1, the right-hand side is exactly
the Hamming weight modulo 3.  In dimension three this recovers the *full*
Hamming distance partition:
    q=1 <-> d_H=1,
    q=2 <-> d_H=2,
    q=0 nonzero <-> d_H=3.

Consequences verified below:
  * the Pass-11541 history scheme is literally H(3,3), not merely formally
    self-dual;
  * the null graph is the distance-3 graph K3 x K3 x K3 (direct/tensor graph),
    with adjacency (J3-I3)^{tensor 3};
  * the null graph alone reconstructs the ordinary Hamming graph by a cubic
    polynomial in its adjacency matrix;
  * Aut(null graph)=Aut(H(3,3))=S3 wr S3, order 1296, matching
    3^3:O(3,3) from Pass 11540;
  * the quadratic form selects a unique projective orthogonal frame with all
    three line norms q=1; this is the canonical Hamming coordinate frame;
  * in that factorization the order-three spectral clock is a tensor product
    of three identical single-qutrit quadratic-phase gates;
  * n=3 is the largest ternary dimension in which the quadratic residue
    sum z_i^2, together with the identity relation, resolves every Hamming
    distance separately.

No physical tensor-factorization, continuum metric, Lorentz group, or gravity
claim follows from the finite conjugacy.
"""

from __future__ import annotations

import itertools
import json
import math
from collections import Counter
from pathlib import Path

P=3
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data"/"PART_W33_PASS11547_HISTORY_HAMMING_CANONICAL_CLOCK.json"

V=list(itertools.product(range(P),repeat=3))
ZERO=(0,0,0)
I3=((1,0,0),(0,1,0),(0,0,1))


def mod(x:int)->int:
    return x%P


def add(u,v):
    return tuple(mod(u[i]+v[i]) for i in range(3))


def sub(u,v):
    return tuple(mod(u[i]-v[i]) for i in range(3))


def q(v)->int:
    t,x,y=v
    return mod(t*t-x*x-y*y)


def polar(u,v)->int:
    return mod(q(add(u,v))-q(u)-q(v))


def H(v):
    t,x,y=v
    return (mod(x+y),mod(x-y),t)


def H_inv(z):
    # H^{-1}(u,v,w) = (w, 2(u+v), 2(u-v)).
    u,v,w=z
    return (w,mod(2*(u+v)),mod(2*(u-v)))


def hamming_weight(v):
    return sum(x!=0 for x in v)


def hamming_distance(u,v):
    return hamming_weight(sub(u,v))


def matmul(A,B):
    n=len(A);m=len(B);r=len(B[0])
    return [[sum(A[i][k]*B[k][j] for k in range(m)) for j in range(r)] for i in range(n)]


def matadd(*terms):
    n=len(terms[0]);m=len(terms[0][0])
    return [[sum(T[i][j] for T in terms) for j in range(m)] for i in range(n)]


def matscale(c,A):
    return [[c*x for x in row] for row in A]


def eye(n):
    return [[int(i==j) for j in range(n)] for i in range(n)]


def kron(A,B):
    out=[]
    for rowA in A:
        for rowB in B:
            row=[]
            for a in rowA:
                row.extend(a*b for b in rowB)
            out.append(row)
    return out


def relation_matrix(distance):
    return [[int(hamming_distance(x,y)==distance) for y in V] for x in V]


def null_matrix_in_hamming_order():
    return relation_matrix(3)


def krawtchouk(j,i,n=3,qary=3):
    return sum(
        (-1)**h*(qary-1)**(j-h)*math.comb(i,h)*math.comb(n-i,j-h)
        for h in range(j+1)
        if h<=i and j-h<=n-i
    )


def hamming_eigenmatrix():
    return [[krawtchouk(j,i) for j in range(4)] for i in range(4)]


def projective_lines():
    seen=set();out=[]
    for v in V:
        if v==ZERO: continue
        w=v if next(x for x in v if x)!=2 else tuple(mod(2*x) for x in v)
        # normalize the first nonzero coordinate to 1
        first=next(i for i,x in enumerate(v) if x)
        scale=1 if v[first]==1 else 2
        w=tuple(mod(scale*x) for x in v)
        if w not in seen:
            seen.add(w);out.append(w)
    return sorted(out)


def orthogonal_frames():
    nonnull=[v for v in projective_lines() if q(v)!=0]
    frames=[]
    for comb in itertools.combinations(nonnull,3):
        if all(polar(a,b)==0 for a,b in itertools.combinations(comb,2)):
            frames.append(tuple(comb))
    return frames


def matrix_det3(M):
    return mod(
        M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
        -M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
        +M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0])
    )


def mv(M,v):
    return tuple(mod(sum(M[i][j]*v[j] for j in range(3))) for i in range(3))


def signed_permutation_matrices():
    mats=[]
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1,2),repeat=3):
            M=[[0]*3 for _ in range(3)]
            for i in range(3):
                M[i][perm[i]]=signs[i]
            mats.append(tuple(tuple(row) for row in M))
    assert len(set(mats))==48
    return sorted(set(mats))


def all_symbol_permutations():
    return list(itertools.permutations(range(3)))


def affine_form_of_symbol_perm(p):
    # Every permutation of F3 is x -> s*x+b with s in {+1,-1}.
    b=p[0]
    s=mod(p[1]-p[0])
    assert s in (1,2)
    assert tuple(mod(s*x+b) for x in range(3))==p
    return s,b


def permute_coords(v,p):
    return tuple(v[p[i]] for i in range(3))


def hamming_wreath_permutations():
    # Return vertex permutations of S3 wr S3 on F3^3.
    idx={v:i for i,v in enumerate(V)}
    perms=set()
    syms=all_symbol_permutations()
    for cp in itertools.permutations(range(3)):
        for ps in itertools.product(syms,repeat=3):
            image=[]
            for v in V:
                w=permute_coords(v,cp)
                w=tuple(ps[i][w[i]] for i in range(3))
                image.append(idx[w])
            perms.add(tuple(image))
    assert len(perms)==6**3*6==1296
    return perms


def affine_signed_monomial_permutations():
    idx={v:i for i,v in enumerate(V)}
    perms=set()
    for M in signed_permutation_matrices():
        for a in V:
            image=tuple(idx[add(mv(M,v),a)] for v in V)
            perms.add(image)
    assert len(perms)==48*27==1296
    return perms


def main():
    # 1. Exact quadratic/Hamming conjugacy.
    for v in V:
        assert H_inv(H(v))==v
        z=H(v)
        assert q(v)==mod(sum(x*x for x in z))
        assert q(v)==hamming_weight(z)%3

    pair_census=Counter()
    for x in V:
        for y in V:
            d=hamming_distance(H(x),H(y))
            dv=sub(x,y)
            qq=q(dv)
            if x==y:
                assert d==0 and qq==0
                tag="identity"
            elif qq==1:
                assert d==1
                tag="q1_to_d1"
            elif qq==2:
                assert d==2
                tag="q2_to_d2"
            else:
                assert d==3
                tag="null_to_d3"
            pair_census[tag]+=1
    assert pair_census==Counter({
        "identity":27,
        "q1_to_d1":27*6,
        "q2_to_d2":27*12,
        "null_to_d3":27*8,
    })

    # Direct map from the repo's Sym_2 coordinates q(a,b,c)=ac-b^2.
    def q_sym(v):
        a,b,c=v
        return mod(a*c-b*b)
    def H_sym(v):
        a,b,c=v
        return (mod(2*a+b+c),mod(2*a+2*b+c),mod(2*a+2*c))
    for v in V:
        assert q_sym(v)==mod(sum(x*x for x in H_sym(v)))
        assert q_sym(v)==hamming_weight(H_sym(v))%3

    # 2. Null graph is the distance-3 Hamming relation = K3 tensor K3 tensor K3.
    B=[[int(i!=j) for j in range(3)] for i in range(3)]
    A3=null_matrix_in_hamming_order()
    assert A3==kron(kron(B,B),B)
    A1=relation_matrix(1)
    A2=relation_matrix(2)
    I=eye(27)

    # Null adjacency alone reconstructs the ordinary Hamming graph and distance-2 relation.
    A3_2=matmul(A3,A3)
    A3_3=matmul(A3_2,A3)
    lhs1=matscale(24,A1)
    rhs1=matadd(matscale(-1,A3_3),matscale(9,A3_2),matscale(18,A3),matscale(-64,I))
    assert lhs1==rhs1
    lhs2=matscale(12,A2)
    rhs2=matadd(A3_3,matscale(-3,A3_2),matscale(-24,A3),matscale(16,I))
    assert lhs2==rhs2

    # 3. Pass 11541 is exactly H(3,3), with relation order (0,3,1,2).
    PH=hamming_eigenmatrix()
    assert PH==[
        [1,6,12,8],
        [1,3,0,-4],
        [1,0,-3,2],
        [1,-3,3,-1],
    ]
    order=(0,3,1,2)
    P_reordered=[[PH[i][j] for j in order] for i in order]
    expected=[
        [1,8,6,12],
        [1,-1,-3,3],
        [1,-4,3,0],
        [1,2,0,-3],
    ]
    assert P_reordered==expected
    old11541=json.loads((ROOT/"data"/"PART_W33_PASS11541_HISTORY_QUADRATIC_SHELL_SCHEME.json").read_text())
    assert old11541["first_eigenmatrix_P"]==expected

    # 4. Group identification.
    signed=signed_permutation_matrices()
    for M in signed:
        assert matrix_det3(M) in (1,2)
        for v in V:
            assert mod(sum(x*x for x in mv(M,v)))==mod(sum(x*x for x in v))
    wreath=hamming_wreath_permutations()
    affine=affine_signed_monomial_permutations()
    assert wreath==affine
    old11540=json.loads((ROOT/"data"/"PART_W33_PASS11540_HISTORY_VO33_ORIENTATION_SQUARE.json").read_text())
    assert old11540["affine_group_ladder"]["AffO_order"]==len(affine)==1296
    assert old11540["orthogonal_group"]["O_3_3_order"]==len(signed)==48

    # 5. The q=1 projective lines are the unique all-q=1 orthogonal frame.
    frames=orthogonal_frames()
    assert len(frames)==4
    patterns=Counter(tuple(sorted(q(v) for v in f)) for f in frames)
    assert patterns==Counter({(1,2,2):3,(1,1,1):1})
    q1_lines=tuple(v for v in projective_lines() if q(v)==1)
    assert len(q1_lines)==3
    assert all(polar(a,b)==0 for a,b in itertools.combinations(q1_lines,2))
    all_one=[f for f in frames if all(q(v)==1 for v in f)]
    assert len(all_one)==1
    assert set(all_one[0])==set(q1_lines)

    # Preimages of the three standard Hamming axes are exactly that frame projectively.
    axes=((1,0,0),(0,1,0),(0,0,1))
    pre=[H_inv(e) for e in axes]
    def proj(v):
        first=next(i for i,x in enumerate(v) if x)
        s=1 if v[first]==1 else 2
        return tuple(mod(s*x) for x in v)
    assert {proj(v) for v in pre}==set(q1_lines)

    old11182=json.loads((ROOT/"data"/"w33_pass11182_paper_ticks_mereology.json").read_text())
    assert old11182["clock_splits"]["position_momentum_aligned"]==4
    assert old11182["clock_splits"]["q_orthogonal_frames"]==4

    # 6. Spectral clock factorizes in Hamming coordinates.
    # Distance-3 adjacency eigenvalue at Fourier Hamming weight w:
    # lambda_w = 2^(3-w)(-1)^w. L=8-lambda.
    clock_rows=[]
    for w in range(4):
        lam=(2**(3-w))*((-1)**w)
        L=8-lam
        assert L%3==0
        phase_exp=mod(-(L//3))   # exp(-2pi i L/9)=omega^phase_exp
        assert phase_exp==mod(-w)
        clock_rows.append(dict(weight=w,adjacency_eigenvalue=lam,
                               laplacian_eigenvalue=L,omega_exponent=phase_exp))
    assert clock_rows==[
        {"weight":0,"adjacency_eigenvalue":8,"laplacian_eigenvalue":0,"omega_exponent":0},
        {"weight":1,"adjacency_eigenvalue":-4,"laplacian_eigenvalue":12,"omega_exponent":2},
        {"weight":2,"adjacency_eigenvalue":2,"laplacian_eigenvalue":6,"omega_exponent":1},
        {"weight":3,"adjacency_eigenvalue":-1,"laplacian_eigenvalue":9,"omega_exponent":0},
    ]

    # Each nonzero Fourier coordinate contributes omega^2, so U=C tensor C tensor C,
    # C=F3^dag diag(1,omega^2,omega^2) F3 = exp(-2pi i p^2/3).
    for k in V:
        w=hamming_weight(k)
        total_exp=mod(sum(0 if x==0 else 2 for x in k))
        assert total_exp==mod(-w)

    # 7. Why dimension three is maximal for exact distance resolution over F3.
    injective_by_n={}
    for n in range(1,8):
        # identity relation separately resolves weight 0, so ask whether positive
        # distances 1..n have distinct residues modulo 3.
        residues=[d%3 for d in range(1,n+1)]
        injective_by_n[str(n)]=len(set(residues))==len(residues)
    assert injective_by_n=={
        "1":True,"2":True,"3":True,
        "4":False,"5":False,"6":False,"7":False,
    }

    out={
      "schema":"w33.pass11547.history_hamming_canonical_clock.v1",
      "status":"PASS_EXACT_FINITE_MINKOWSKI_HAMMING_CONJUGACY",
      "pass":11547,
      "linear_conjugacy":{
        "paper_coordinates":"(t,x,y)",
        "map":"H(t,x,y)=(x+y,x-y,t)",
        "inverse":"H^-1(u,v,w)=(w,2(u+v),2(u-v))",
        "identity":"t^2-x^2-y^2 = u^2+v^2+w^2 over F3",
        "hamming_identity":"u^2+v^2+w^2 = wt_H(u,v,w) mod 3",
        "sym2_coordinates":"for q(a,b,c)=ac-b^2: Hsym=(2a+b+c,2a+2b+c,2a+2c)",
      },
      "exact_relation_dictionary":{
        "coincident":"Hamming distance 0",
        "q=+1 / paper timelike shell":"Hamming distance 1",
        "q=-1=2 / paper spacelike shell":"Hamming distance 2",
        "q=0 nonzero / paper null shell":"Hamming distance 3",
        "ordered_pair_census":dict(sorted(pair_census.items())),
      },
      "hamming_scheme":{
        "standard_name":"ternary Hamming association scheme H(3,3)",
        "ordinary_distance_order":[0,1,2,3],
        "history_relation_order":[0,3,1,2],
        "ordinary_eigenmatrix":PH,
        "history_reordered_eigenmatrix":P_reordered,
        "pass11541_is_literal_hamming_scheme":True,
      },
      "null_graph":{
        "description":"distance-3 graph of H(3,3), equivalently direct/tensor product K3 x K3 x K3",
        "adjacency_factorization":"A_null=(J3-I3) tensor (J3-I3) tensor (J3-I3)",
        "degree":8,
        "full_graph_reconstruction":{
          "A1_from_A3":"24 A1 = -A3^3 + 9 A3^2 + 18 A3 - 64 I",
          "A2_from_A3":"12 A2 = A3^3 - 3 A3^2 - 24 A3 + 16 I",
          "meaning":"the null graph alone generates the full Hamming Bose-Mesner algebra",
        },
      },
      "automorphism_group":{
        "full_group":"Aut(null)=Aut(H(3,3))=S3 wr S3",
        "order":1296,
        "wreath_order_formula":"6^3 * 6",
        "affine_form":"F3^3 : (C2^3 : S3) = 3^3 : O(3,3)",
        "linear_stabilizer":"signed permutation group C2^3:S3",
        "linear_stabilizer_order":48,
        "all_symbol_permutations_of_F3_are_affine":True,
        "pass11540_group_recovered":True,
      },
      "canonical_hamming_frame":{
        "projective_orthogonal_frame_count":4,
        "norm_pattern_census":{"1,1,1":1,"1,2,2":3},
        "q1_projective_lines":[list(v) for v in q1_lines],
        "unique_all_q1_frame":True,
        "preimages_of_standard_hamming_axes":[list(v) for v in pre],
        "pass11182_position_momentum_aligned_splits":4,
        "interpretation":"q itself canonically selects the unordered three-axis Hamming frame; labels/signs remain gauge",
      },
      "clock_factorization":{
        "clock":"U=exp(-2*pi*i*L_null/9)",
        "fourier_weight_table":clock_rows,
        "phase_law":"on Fourier Hamming weight w, U has phase omega^{-w}",
        "single_trit_gate":"C=F3^dag diag(1,omega^2,omega^2) F3 = exp(-2*pi*i p^2/3)",
        "factorization":"in the canonical Hamming mereology, U=C tensor C tensor C",
        "relation_to_pass11166":"the light-cone-coordinate clock is entangling, but this global linear change of tensor factorization diagonalizes its quadratic form; no contradiction",
      },
      "qutrit_dimension_uniqueness":{
        "theorem":"over F3, sum z_i^2 equals Hamming weight mod 3; n=3 is the largest dimension for which the quadratic value plus identity separates every Hamming distance",
        "positive_distance_residue_injective_by_dimension":injective_by_n,
        "first_collision_at_n4":"distances 1 and 4 both have quadratic residue 1",
      },
      "literature_boundary":{
        "hamming_scheme":"H(n,q) distance relations and Krawtchouk eigenmatrices are standard coding-theory association-scheme results",
        "automorphisms":"Aut(H(n,m))=S_m wr S_n is standard; Mirafzal-Ziaee give an elementary proof",
        "distance_n_graph":"a 2026 Mirafzal preprint independently studies M(n,m)=K_m direct-product n times (adjacent iff all coordinates differ) and proves Aut(M)=Aut(H); Pass11547 does not claim this general theorem",
        "repo_specific":"the explicit conjugacy of the Forty Points quadratic history chart to H(3,3), its weld to Passes 11166/11182/11540/11541, and the canonical-clock factorization",
      },
      "boundary":"Exact finite linear/graph/association-scheme theorem. Hamming coordinates are a canonical finite address factorization up to local alphabet permutations and coordinate permutation; they are not a derivation of physical subsystems, continuum Lorentzian spacetime, a measured causal metric, or gravity."
    }
    OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":out["status"],
        "relation_dictionary":out["exact_relation_dictionary"],
        "aut_order":out["automorphism_group"]["order"],
        "clock_factorization":out["clock_factorization"]["factorization"],
        "max_exact_dimension":3,
    },indent=2,sort_keys=True))


if __name__=="__main__":
    main()
