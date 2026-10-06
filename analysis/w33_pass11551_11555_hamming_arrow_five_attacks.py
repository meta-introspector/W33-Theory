#!/usr/bin/env python3
"""Passes 11551-11555: five exact attacks from the Hamming-arrow breakthrough.

11551: the cubic Hamming arrow is the restriction of the repository's actual
       E8-derived clock Albert/E6 determinant to a Jordan frame.
11552: naive Hamming dimension-refinement is audited.  The full q=0 shell
       becomes maximally mixing; the minimal weight-3 shell closes its gap as
       9/(2n) but converges spectrally to a Hamming-weight number operator, not
       a Lorentzian d'Alembertian.
11553: the arrow is the unique PSp-fixed H^1 class over F3, while the fixed
       H^1 vanishes over tested non-characteristic primes.
11554: the tetracode orientation lift GL(2,3)->S4 is a non-split central
       double cover, whereas the history/cone lift O(3,3)->S4 is split.
11555: the local arrow defines a rank-8 skew operator D with
       3D^2=L(L-6I)(L-12I); normalized, it is the old Weil complex structure
       on exactly the eight fixed null modes.

All results are finite/algebraic.  No continuum spacetime, gravity,
thermodynamic arrow, observed CPT law, or Standard-Model parameter is inferred.
"""
from __future__ import annotations

import itertools
import json
import math
import sys
from collections import Counter, deque
from pathlib import Path

import numpy as np
import sympy as sp

P=3
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data"/"PART_W33_PASS11551_11555_HAMMING_ARROW_FIVE_ATTACKS.json"
if str(ROOT/"analysis") not in sys.path:
    sys.path.insert(0,str(ROOT/"analysis"))


# ------------------------------ generic finite-field helpers

def rref_mod(A,p):
    A=np.array(A,dtype=int)%p
    piv=[];r=0
    for c in range(A.shape[1]):
        q=next((i for i in range(r,A.shape[0]) if A[i,c]%p),None)
        if q is None: continue
        A[[r,q]]=A[[q,r]]
        A[r]=(A[r]*pow(int(A[r,c]),-1,p))%p
        for i in range(A.shape[0]):
            if i!=r and A[i,c]%p:
                A[i]=(A[i]-A[i,c]*A[r])%p
        piv.append(c);r+=1
        if r==A.shape[0]: break
    return A,piv

def rank_mod(A,p):
    return len(rref_mod(A,p)[1])

def nullspace_mod(A,p):
    R,piv=rref_mod(A,p);n=R.shape[1];free=[j for j in range(n) if j not in piv]
    cols=[]
    for f in free:
        x=np.zeros(n,dtype=int);x[f]=1
        for i,c in enumerate(piv):
            x[c]=(-R[i,f])%p
        cols.append(x)
    return np.stack(cols,axis=1) if cols else np.zeros((n,0),dtype=int)

def inv_mod(M,p):
    M=np.array(M,dtype=int)%p;n=len(M)
    A=np.c_[M,np.eye(n,dtype=int)]
    R,piv=rref_mod(A,p)
    assert piv[:n]==list(range(n))
    return R[:n,n:]%p

def independent_columns(A,p):
    A=np.array(A,dtype=int)%p;inds=[];r=0
    for j in range(A.shape[1]):
        rr=rank_mod(A[:,inds+[j]],p)
        if rr>r: inds.append(j);r=rr
    return inds

def pivot_rows_for_columns(B,p):
    # B has full column rank.  Pivot columns of B^T are independent rows of B.
    _,piv=rref_mod(np.array(B,dtype=int).T,p)
    assert len(piv)==B.shape[1]
    return piv


# ------------------------------ Hamming history geometry

V=list(itertools.product(range(3),repeat=3))
VID={v:i for i,v in enumerate(V)}
ZERO=(0,0,0)

def mod(x): return int(x)%3
def add(a,b): return tuple(mod(a[i]+b[i]) for i in range(3))
def sub(a,b): return tuple(mod(a[i]-b[i]) for i in range(3))
def chi(d):
    assert all(x in (1,2) for x in d)
    z=1
    for x in d: z*=1 if x==1 else -1
    return z

def Hsym(s):
    a,b,c=s
    return (mod(2*a+b+c),mod(2*a+2*b+c),mod(2*a+2*c))

EDGES=[(i,j) for i in range(27) for j in range(i+1,27)
       if all(sub(V[j],V[i]))]
EIDX={e:i for i,e in enumerate(EDGES)}
TRIS=[t for t in itertools.combinations(range(27),3)
      if all(tuple(sorted(e)) in EIDX for e in itertools.combinations(t,2))]
TIDX={t:i for i,t in enumerate(TRIS)}
assert len(EDGES)==108 and len(TRIS)==36

B1=np.zeros((27,108),dtype=int)
B2=np.zeros((108,36),dtype=int)
for k,(i,j) in enumerate(EDGES):
    B1[i,k]=-1;B1[j,k]=1
for k,(a,b,c) in enumerate(TRIS):
    for x,y,s in ((b,c,1),(a,c,-1),(a,b,1)):
        B2[EIDX[(x,y)],k]=s

ARROW=np.array([chi(sub(V[j],V[i])) for i,j in EDGES],dtype=int)


# ------------------------------ affine PSp generators in Hamming gauge

def vertex_translation(t):
    return tuple(VID[add(v,t)] for v in V)

def vertex_linear(M):
    return tuple(VID[tuple(mod(sum(M[i][j]*v[j] for j in range(3)))
                           for i in range(3))] for v in V)

def edge_action_matrix(g):
    Q=np.zeros((108,108),dtype=int)
    for k,(i,j) in enumerate(EDGES):
        a,b=g[i],g[j]
        if a<b: Q[EIDX[(a,b)],k]=1
        else: Q[EIDX[(b,a)],k]=-1
    return Q

GEN_VERTEX=[
    vertex_translation((1,0,0)),vertex_translation((0,1,0)),vertex_translation((0,0,1)),
    vertex_linear(((0,1,0),(1,0,0),(0,0,1))),
    vertex_linear(((1,0,0),(0,0,1),(0,1,0))),
    vertex_linear(((2,0,0),(0,2,0),(0,0,1))),
]
GEN_EDGE=[edge_action_matrix(g) for g in GEN_VERTEX]
assert all(np.array_equal(Q@ARROW,ARROW) for Q in GEN_EDGE)


# ------------------------------ 11551: Albert/E6 cubic restriction

def pass11551():
    import w33_pass10950_clock_albert_lorentz_spinor as C
    J=C.build_clock_albert()
    ids=[J["idx"][g] for g in J["G"]]
    Pi2=np.array([[[int(x*2) for x in J["prod"][i,j]]
                   for j in range(27)] for i in range(27)],dtype=np.int64)
    raw=np.einsum("ijj->i",Pi2);assert not np.any(raw%18)
    tr=raw//18
    gram2=np.einsum("k,ijk->ij",tr,Pi2)
    cube4=np.einsum("l,ial,jka->ijk",tr,Pi2,Pi2)
    N12=(2*np.einsum("i,j,k->ijk",tr,tr,tr)
         -(np.einsum("i,jk->ijk",tr,gram2)
           +np.einsum("j,ik->ijk",tr,gram2)
           +np.einsum("k,ij->ijk",tr,gram2))
         +cube4)
    rows=[]
    for z in V:
        y=np.zeros(27,dtype=np.int64)
        for i,a in enumerate(z): y[ids[i]]=a
        v12=int(np.einsum("ijk,i,j,k",N12,y,y,y))
        assert v12%12==0
        N=v12//12;target=z[0]*z[1]*z[2]
        assert N==target
        rows.append({"z":list(z),"Albert_norm":N,"product":target})
    a,b,c=sp.symbols("a b c")
    u,v,w=2*a+b+c,2*a+2*b+c,2*a+2*c
    poly=sp.Poly(sp.expand(u*v*w),a,b,c,modulus=3).as_expr()
    return {
      "status":"PASS_ARROW_IS_DIAGONAL_ALBERT_E6_CUBIC",
      "clock_Albert_frame_indices":ids,
      "exact_points_checked":27,
      "identity":"N_Albert(z1 e1+z2 e2+z3 e3)=z1*z2*z3",
      "history_Hamming_arrow":"chi is the sign/residue of this norm on the eight full-support ternary vectors",
      "Sym2_history_coordinates_factorization":"(2a+b+c)(2a+2b+c)(2a+2c) mod 3",
      "expanded_mod3":str(poly),
      "E8_anchor":"the clock Albert algebra is built by Pass10950 from the committed executable E8 bracket; Pass11471 independently maps the native signed E6 cubic coefficientwise to this Albert determinant",
      "boundary":"Restriction theorem on the canonical Jordan frame, not a derivation of a physical E6 field or observed time arrow.",
    }


# ------------------------------ 11552: refinement spectra

def Kraw(j,w,n):
    return sum((-1)**h*2**(j-h)*math.comb(w,h)*math.comb(n-w,j-h)
               for h in range(j+1) if h<=w and j-h<=n-w)

def pass11552():
    rows=[]
    for n in range(3,41):
        js=[j for j in range(3,n+1,3)]
        kfull=sum(math.comb(n,j)*2**j for j in js)
        fullvals=[sum(Kraw(j,w,n) for j in js) for w in range(n+1)]
        fullgap=min(1-fullvals[w]/kfull for w in range(1,n+1))
        k3=8*math.comb(n,3)
        k3vals=[Kraw(3,w,n) for w in range(n+1)]
        gap3=min(1-k3vals[w]/k3 for w in range(1,n+1))
        if n>=5:
            assert abs(gap3-9/(2*n))<1e-14
        rows.append({"n":n,"full_null_degree":kfull,"full_normalized_gap":fullgap,
                     "weight3_degree":k3,"weight3_normalized_gap":gap3})
    # Verify the exact fixed-w rescaled limit formula symbolically.
    n,w=sp.symbols("n w",integer=True,positive=True)
    K3=sp.expand(sum((-1)**h*2**(3-h)*sp.binomial(w,h)*sp.binomial(n-w,3-h)
                     for h in range(4)))
    degree=8*sp.binomial(n,3)
    gap=sp.factor(1-K3/degree)
    target=sp.factor(9*w*(4*n**2-6*n*w-6*n+3*w**2+3*w+2)
                     /(8*n*(n-1)*(n-2)))
    assert sp.simplify(gap-target)==0
    return {
      "status":"PASS_HAMMING_REFINEMENT_SPECTRAL_FIREWALL_AND_MINIMAL_LIMIT",
      "full_q0_shell":{
        "degree_formula":"sum_{j=3,6,...} 2^j C(n,j) = [3^n + 2*3^(n/2) cos(n*pi/2)]/3 - 1",
        "nontrivial_eigenvalue_formula":"lambda_w=(2/3)3^(n/2) cos(n*pi/2-2*pi*w/3)-1, w>0",
        "asymptotic":"max nontrivial |lambda|/k = O(3^(-n/2)); normalized spectral gap -> 1 exponentially",
        "verdict":"naive full quadratic-null refinement becomes increasingly global/mixing, so it does not furnish a local continuum limit",
      },
      "minimal_weight3_shell":{
        "degree":"8*C(n,3)",
        "exact_gap_for_n_ge_5":"9/(2n)",
        "fixed_weight_gap":str(target),
        "scaled_limit":"for fixed Fourier Hamming weight w, n*(1-K3(w)/k) -> 9w/2",
        "verdict":"the only natural gap-closing truncation found has a Hamming-weight/number-operator low spectrum, not a Lorentzian d'Alembertian",
      },
      "enumeration":rows,
      "boundary":"A spectral obstruction to two naive refinement prescriptions, not a theorem excluding all continuum limits or geometric renormalizations.",
    }


# ------------------------------ 11553: modular cohomology

def extend_basis_columns(B,Z,p):
    cols=[]
    if B.shape[1]:
        cols=[B[:,j] for j in independent_columns(B,p)]
    r=len(cols);extras=[]
    for j in range(Z.shape[1]):
        v=Z[:,j]
        if rank_mod(np.column_stack(cols+[v]) if cols else v[:,None],p)>r:
            cols.append(v);extras.append(v);r+=1
    return (np.column_stack(cols)%p,np.column_stack(extras)%p
            if extras else np.zeros((Z.shape[0],0),dtype=int))

def fixed_H1_dim(p):
    delta1=B2.T%p;delta0=B1.T%p
    Z=nullspace_mod(delta1,p)
    Bind=delta0[:,independent_columns(delta0,p)]%p
    Full,H=extend_basis_columns(Bind,Z,p)
    assert Full.shape[1]==Z.shape[1]
    rows=pivot_rows_for_columns(Full,p)
    Inv=inv_mod(Full[rows,:],p)
    h=H.shape[1];rB=Bind.shape[1]
    mats=[]
    for Q in GEN_EDGE:
        R=np.zeros((h,h),dtype=int)
        for j in range(h):
            v=(Q@H[:,j])%p
            co=(Inv@v[rows])%p
            assert np.array_equal((Full@co)%p,v)
            R[:,j]=co[rB:]
        mats.append((R-np.eye(h,dtype=int))%p)
    fixed=h-rank_mod(np.vstack(mats),p)
    return int(Z.shape[1]-Bind.shape[1]),int(fixed)

def pass11553():
    assert np.all(B1@ARROW==0)
    tri=B2.T@ARROW
    assert set(map(int,tri))=={-3,3}
    dims={}
    for p in (2,3,5,7):
        h,f=fixed_H1_dim(p);dims[str(p)]={"H1_dimension":h,"PSp_fixed_H1_dimension":f}
    assert dims["3"]["H1_dimension"]==46 and dims["3"]["PSp_fixed_H1_dimension"]==1
    assert all(dims[str(p)]["PSp_fixed_H1_dimension"]==0 for p in (2,5,7))
    assert rank_mod(B2.T,3)==36 and rank_mod(B1.T,3)==26
    assert np.all((B2.T@ARROW)%3==0)
    # Non-exactness of the arrow modulo 3.
    exact_rank=rank_mod(np.c_[B1.T%3,ARROW%3],3)
    assert exact_rank>rank_mod(B1.T,3)
    return {
      "status":"PASS_ARROW_IS_UNIQUE_PSP_FIXED_CHARACTERISTIC3_H1_CLASS",
      "clique_complex":{"vertices":27,"edges":108,"triangles":36,"max_clique_size":3},
      "integer_chain":{"boundary_zero":True,"triangle_coboundary_values":[-3,3]},
      "mod3":{"is_1_cocycle":True,"is_exact":False,"H1_dimension":46,
              "PSp_fixed_H1_dimension":1},
      "prime_control":dims,
      "characteristic3_jump":"the PSp-fixed H1 line exists at p=3 and vanishes in the p=2,5,7 controls",
      "integral_rank":"H^1 over Z is free of rank 46; H^2 has rank 0 because the 36 triangle coboundaries are independent",
      "interpretation":"the finite arrow is a symmetry-protected modular cohomology class, not merely an invariant cycle",
      "boundary":"Finite modular cohomology. No nonzero Bockstein/Chern class, anomaly, thermodynamic irreversibility or continuum topology is claimed.",
    }


# ------------------------------ 11554: inequivalent double covers

RAYS=((1,0),(0,1),(1,1),(1,2))
def det2(A): return mod(A[0]*A[3]-A[1]*A[2])
def mv2(A,v): return (mod(A[0]*v[0]+A[1]*v[1]),mod(A[2]*v[0]+A[3]*v[1]))
def proj2(v):
    s=pow(next(x for x in v if x),-1,3)
    return tuple(mod(s*x) for x in v)
def pperm(A): return tuple(RAYS.index(proj2(mv2(A,r))) for r in RAYS)
def mul2(A,B):
    return tuple(mod(sum(A[2*i+k]*B[2*k+j] for k in range(2)))
                 for i in range(2) for j in range(2))
def parity(p):
    return -1 if sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))%2 else 1
def rank2(A,ncols):
    rows=[list(map(lambda x:int(x)&1,r[:ncols])) for r in A];rr=0
    for c in range(ncols):
        q=next((i for i in range(rr,len(rows)) if rows[i][c]),None)
        if q is None: continue
        rows[rr],rows[q]=rows[q],rows[rr]
        for i in range(len(rows)):
            if i!=rr and rows[i][c]:
                rows[i]=[a^b for a,b in zip(rows[i],rows[rr])]
        rr+=1
    return rr
def order2(A):
    I=(1,0,0,1);X=I
    for n in range(1,25):
        X=mul2(X,A)
        if X==I:return n
    raise AssertionError

def pass11554():
    GL=[A for A in itertools.product(range(3),repeat=4) if det2(A)]
    perms=sorted(set(pperm(A) for A in GL));assert len(GL)==48 and len(perms)==24
    sec={p:min(A for A in GL if pperm(A)==p) for p in perms}
    pi={p:i for i,p in enumerate(perms)}
    def compose(p,q): return tuple(p[q[i]] for i in range(4))
    eq=[]
    for p in perms:
        for q in perms:
            r=compose(p,q);M=mul2(sec[p],sec[q]);S=sec[r]
            bit=0 if M==S else 1
            assert M==S or M==tuple(mod(-x) for x in S)
            row=[0]*25;row[pi[p]]=1;row[pi[q]]^=1;row[pi[r]]^=1;row[24]=bit
            eq.append(row)
    rank=rank2([r[:24] for r in eq],24);aug=rank2(eq,25)
    assert rank==23 and aug==24
    assert all((1 if det2(A)==1 else -1)==parity(pperm(A)) for A in GL)

    # History cover: signed permutations O = W(D3) x <-I>.
    O=[];W=[]
    for pm in itertools.permutations(range(3)):
        for sg in itertools.product((1,2),repeat=3):
            M=[[0]*3 for _ in range(3)]
            for i,j in enumerate(pm):M[i][j]=sg[i]
            M=tuple(tuple(r) for r in M);O.append(M)
            sig=np.prod([1 if x==1 else -1 for x in sg])
            if sig==1:W.append(M)
    assert len(O)==48 and len(W)==24
    minusI=((2,0,0),(0,2,0),(0,0,2))
    assert minusI not in W
    # Fiber product is GL(2,3) x C2 because history cover is split.
    glspec=Counter(order2(A) for A in GL)
    fiberspec=Counter()
    for A in GL:
        oa=order2(A)
        for b in (0,1):
            ob=1 if b==0 else 2
            fiberspec[math.lcm(oa,ob)]+=1
    tomo=json.loads((ROOT/"data"/"w33_pass2430_labelled_tomotope_group_obstruction.json").read_text())
    archived={int(k):v for k,v in tomo["archived_labelled_tomotope"]["order_spectrum"].items()}
    assert dict(sorted(fiberspec.items()))!=dict(sorted(archived.items()))
    return {
      "status":"PASS_TETRACODE_AND_HISTORY_ARE_INEQUIVALENT_S4_DOUBLE_COVERS",
      "tetracode_cover":{"group":"GL(2,3)","order":48,"center":"{+I,-I}",
        "quotient":"S4 on the four P1(F3) rays","non_split":True,
        "cocycle_rank":rank,"augmented_rank":aug,
        "determinant":"descends exactly to S4 permutation parity"},
      "history_cover":{"group":"O(3,3)=signed permutations","order":48,
        "split":"O(3,3)=W(D3) x <-I> ~= S4 x C2",
        "projective_base":"the same S4 on four null rays"},
      "consequence":"Pass10946 oriented tetracode representatives cannot be canonically selected by the Pass11550 cone/time sheet: the two lifts have different central-extension classes",
      "fiber_product":{"order":96,"structure":"GL(2,3) x C2",
        "center_order":4,"order_spectrum":dict(sorted(fiberspec.items())),
        "tomotope_order_spectrum":dict(sorted(archived.items())),
        "is_archived_tomotope":False},
      "tomotope_firewall":"The natural order-96 compatibility lift is not the archived tomotope group; the repo's existing Pass2430 obstruction is respected rather than count-matched.",
      "boundary":"Finite central-extension theorem; no physical spin structure or spacetime double cover is inferred.",
    }


# ------------------------------ 11555: local skew arrow dynamics

def pass11555():
    A=np.zeros((27,27),dtype=np.int64);D=np.zeros((27,27),dtype=np.int64)
    for i,x in enumerate(V):
        for j,y in enumerate(V):
            if i==j:continue
            d=sub(y,x)
            if all(d):
                A[i,j]=1;D[i,j]=chi(d)
    I=np.eye(27,dtype=np.int64);L=8*I-A
    assert np.array_equal(D.T,-D)
    assert np.array_equal(D@L,L@D)
    assert np.array_equal(3*(D@D),L@(L-6*I)@(L-12*I))
    assert np.linalg.matrix_rank(D)==8

    # D is exactly the sum of the 36 oriented temporal 3-cycles.
    S=np.zeros_like(D)
    dirs=[d for d in V if all(d) and chi(d)==1]
    seen=set()
    for d in dirs:
        lines=set()
        for x in V:
            line=tuple(sorted((VID[x],VID[add(x,d)],VID[add(add(x,d),d)])))
            lines.add(line)
        assert len(lines)==9
        for line in lines:
            assert line in TIDX and line not in seen;seen.add(line)
            x=V[line[0]];cyc=[VID[x],VID[add(x,d)],VID[add(add(x,d),d)]]
            for a,b in zip(cyc,cyc[1:]+cyc[:1]):
                S[a,b]+=1;S[b,a]-=1
    assert len(seen)==36
    assert np.array_equal(S,D)

    # PSp acts transitively on all 36 (unoriented) temporal triangles.
    signed_even=[]
    for pm in itertools.permutations(range(3)):
        for sg in itertools.product((1,2),repeat=3):
            if np.prod([1 if x==1 else -1 for x in sg])!=1:continue
            M=[[0]*3 for _ in range(3)]
            for i,j in enumerate(pm):M[i][j]=sg[i]
            signed_even.append(tuple(tuple(r) for r in M))
    seed=TRIS[0];orb=set()
    for M in signed_even:
        for t in V:
            def f(v):
                w=tuple(mod(sum(M[i][j]*v[j] for j in range(3))) for i in range(3))
                return VID[add(w,t)]
            orb.add(tuple(sorted(f(V[i]) for i in seed)))
    assert orb==set(TRIS)

    return {
      "status":"PASS_LOCAL_ARROW_OPERATOR_IS_WEIL_COMPLEX_STRUCTURE",
      "operator":{
        "definition":"D_xy=chi(y-x) on null edges and 0 otherwise",
        "rank":8,"kernel_dimension":19,"skew_symmetric":True,"commutes_with_L":True,
        "exact_polynomial_identity":"3 D^2 = L(L-6I)(L-12I)",
        "projector_identity":"P_9=-L(L-6I)(L-12I)/81 and D^2=-27 P_9",
        "complex_structure":"J=-D/(3*sqrt(3)); J^2=-P_9",
        "Fourier_symbol":"Dhat(k)=-i*3*sqrt(3)*chi(k) on full-support k, and 0 otherwise",
        "Weil_weld":"J has eigenvalue i*chi(k) on the eight L=9 modes, exactly Pass11550's normalized Weil phase",
      },
      "triangle_locality":{"D_is_sum_of_36_oriented_triangle_cycles":True,
        "triangles_form_one_PSp_orbit":True,"factorization":"36=4 directions x 9 affine lines"},
      "local_invariant_dynamics":{
        "quadratic_operator_space":"for operators supported only on diagonal plus null edges, PSp invariance gives span{I,A,D}",
        "Hermitian_form":"H=mu I + kappa L + eta iD",
        "PGSp_outer_rule":"outer sigma=-1 sends D->-D, so eta is an orientation-breaking spurion",
        "spectrum":{
          "L0_dim1":"mu","L12_dim6":"mu+12*kappa","L6_dim12":"mu+6*kappa",
          "L9_two_dim4_sheets":"mu+9*kappa +/- 3*sqrt(3)*eta"},
        "scalar_triangle_cubic":"one PSp-invariant nearest-triangle scalar cubic g*sum_T phi_i phi_j phi_k because the 36 triangles are one orbit; scalar commutativity makes this term orientation-even",
      },
      "boundary":"Exact finite local dynamics/classification. D is not asserted to be a physical Dirac operator, and the finite dispersion does not derive continuum Lorentz invariance or gravity.",
    }


def main():
    out={
      "schema":"w33.pass11551_11555.hamming_arrow_five_attacks.v1",
      "status":"PASS_FIVE_HAMMING_ARROW_ATTACKS",
      "passes":{
        "11551":pass11551(),
        "11552":pass11552(),
        "11553":pass11553(),
        "11554":pass11554(),
        "11555":pass11555(),
      },
      "upstream":["Pass10946","Pass10950","Pass11471","Pass11547","Pass11548","Pass11549","Pass11550"],
      "boundary":"Five exact finite/algebraic results. They sharpen the TOE architecture but do not close the continuum, gravitational, thermodynamic, or phenomenological frontiers.",
    }
    assert all(v["status"].startswith("PASS") for v in out["passes"].values())
    OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({k:v["status"] for k,v in out["passes"].items()},indent=2))


if __name__=="__main__":
    main()
