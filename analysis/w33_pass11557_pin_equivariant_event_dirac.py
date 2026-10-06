#!/usr/bin/env python3
"""Pass 11557: Pin-equivariant first-order Dirac operator on the 27-event history chart.

Pass 11555 found the scalar cubic-arrow operator D.  Here the eight null
steps are resolved into three vector-weighted difference operators B_i.
They satisfy D=B_1 B_2 B_3 exactly.  Contracting B_i with a graded Cl(3)
carrier gives a genuinely first-order finite operator Q whose square is
4L+D^2.  The three Weyl reflections of W(D3) lift to a 48-element Pin+
double cover with the same order census as the committed GL(2,3)=2^+S4
clock cover, and Q is invariant under the combined event/spin action.
"""
from __future__ import annotations
import itertools, json
from collections import Counter, deque
from pathlib import Path
import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/PART_W33_PASS11557_PIN_EQUIVARIANT_EVENT_DIRAC.json"

def small_clifford():
    I=sp.I
    I2=sp.eye(2)
    X=sp.Matrix([[0,1],[1,0]])
    Y=sp.Matrix([[0,-I],[I,0]])
    Z=sp.diag(1,-1)
    gam=[sp.kronecker_product(X,X),
         sp.kronecker_product(X,Y),
         sp.kronecker_product(X,Z)]
    grade=sp.kronecker_product(Z,I2)
    I4=sp.eye(4)
    for a in range(3):
        for b in range(3):
            rhs=(2 if a==b else 0)*I4
            assert sp.simplify(gam[a]*gam[b]+gam[b]*gam[a]-rhs)==sp.zeros(4)
        assert sp.simplify(grade*gam[a]+gam[a]*grade)==sp.zeros(4)
    assert grade*grade==I4
    return gam,grade

def roots_and_reflections(gam,grade):
    roots=[(1,-1,0),(0,1,-1),(0,1,1)]
    refs=[]; lifts=[]
    for a in roots:
        av=sp.Matrix(a)
        R=sp.eye(3)-av*av.T  # roots have norm^2=2
        p=sum((sp.Integer(a[j])*gam[j] for j in range(3)),sp.zeros(4))/sp.sqrt(2)
        S=sp.simplify(sp.I*grade*p)
        assert sp.simplify(S*S-sp.eye(4))==sp.zeros(4)
        for k in range(3):
            lhs=sp.simplify(S*gam[k]*S)
            rhs=sum((R[j,k]*gam[j] for j in range(3)),sp.zeros(4))
            assert sp.simplify(lhs-rhs)==sp.zeros(4)
        refs.append(np.array(R.tolist(),dtype=int))
        lifts.append(S)
    return roots,refs,lifts

def _nkey(M):
    A=np.array(M,dtype=complex)
    return tuple(np.round(A.real,12).ravel())+tuple(np.round(A.imag,12).ravel())

def enumerate_pin_cover(refs,lifts):
    gens=[np.array(S.evalf(50),dtype=complex) for S in lifts]
    I4=np.eye(4,dtype=complex); I3=np.eye(3,dtype=int)
    group={_nkey(I4):(I3,I4)}
    q=deque([(I3,I4)])
    while q:
        R,S=q.popleft()
        for rg,sg in zip(refs,gens):
            R2=R@rg;S2=S@sg;k=_nkey(S2)
            if k not in group:
                group[k]=(R2,S2);q.append((R2,S2))
    assert len(group)==48
    fibers=Counter(tuple(R.ravel()) for R,S in group.values())
    assert len(fibers)==24 and set(fibers.values())=={2}
    assert sum(np.array_equal(R,I3) for R,S in group.values())==2
    def order(S):
        P=I4.copy()
        for n in range(1,25):
            P=P@S
            if np.linalg.norm(P-I4)<1e-9:return n
        raise AssertionError("order")
    census=Counter(order(S) for R,S in group.values())
    assert census==Counter({2:13,8:12,3:8,6:8,4:6,1:1})
    return group,census

def event_operators():
    V=list(itertools.product(range(3),repeat=3));vid={v:i for i,v in enumerate(V)}
    Bs=[np.zeros((27,27),dtype=np.int64) for _ in range(3)]
    A=np.zeros((27,27),dtype=np.int64);D=np.zeros((27,27),dtype=np.int64)
    for i,x in enumerate(V):
        for j,y in enumerate(V):
            if i==j:continue
            d=tuple((y[k]-x[k])%3 for k in range(3))
            if all(d):
                e=[1 if z==1 else -1 for z in d]
                A[i,j]=1
                for k in range(3):Bs[k][i,j]=e[k]
                D[i,j]=e[0]*e[1]*e[2]
    L=8*np.eye(27,dtype=np.int64)-A
    assert all(np.array_equal(B.T,-B) for B in Bs)
    assert np.array_equal(Bs[0]@Bs[1]@Bs[2],D)
    assert np.array_equal(-sum(B@B for B in Bs),4*L+D@D)
    return V,vid,Bs,A,D,L

def main():
    gam,grade=small_clifford()
    roots,refs,lifts=roots_and_reflections(gam,grade)
    group,census=enumerate_pin_cover(refs,lifts)
    V,vid,Bs,A,D,L=event_operators()

    # Event covariance under all 24 orthogonal images.
    seen=set(); max_event=0
    for R,S in group.values():
        rk=tuple(R.ravel())
        if rk in seen: continue
        seen.add(rk)
        g=[]
        for x in V:
            y=tuple(int(sum(R[a,b]*x[b] for b in range(3))%3) for a in range(3))
            g.append(vid[y])
        P=np.zeros((27,27),dtype=np.int64)
        for x,y in enumerate(g):P[y,x]=1
        for k in range(3):
            lhs=P@Bs[k]@P.T
            rhs=sum(int(R[j,k])*Bs[j] for j in range(3))
            max_event=max(max_event,int(np.max(np.abs(lhs-rhs))))
    assert len(seen)==24 and max_event==0

    # Full combined covariance checked numerically on all 48 Pin lifts.
    gn=[np.array(g.evalf(40),dtype=complex) for g in gam]
    Q=1j*sum(np.kron(gn[i],Bs[i]) for i in range(3))
    max_q=0.0
    for R,S in group.values():
        g=[]
        for x in V:
            y=tuple(int(sum(R[a,b]*x[b] for b in range(3))%3) for a in range(3))
            g.append(vid[y])
        P=np.zeros((27,27),dtype=float)
        for x,y in enumerate(g):P[y,x]=1
        U=np.kron(S,P)
        err=np.linalg.norm(U@Q@U.conj().T-Q)
        max_q=max(max_q,float(err))
    assert max_q<1e-9

    p51=json.loads((ROOT/"data/w33_pass10951_clock_pin_spin_central_sign_bridge.json").read_text())
    old=Counter()
    for key,n in p51["lift_order_census"].items():
        order=int(key.split("order")[1]);old[order]+=int(n)
    assert old==census

    out={
      "schema":"w33.pass11557.pin-equivariant-event-dirac.v1",
      "status":"PASS_PIN_EQUIVARIANT_FIRST_ORDER_EVENT_DIRAC","pass":11557,
      "vector_differences":{
        "operators":["B1","B2","B3"],
        "definition":"B_i[x,y]=epsilon_i for y-x=epsilon in {+/-1}^3, else 0",
        "skew_symmetric":True,
        "cubic_factorization":"D = B1 B2 B3",
        "D_rank":int(np.linalg.matrix_rank(D.astype(float))),
      },
      "dirac_operator":{
        "formula":"Q = i sum_i Gamma_i tensor B_i",
        "spinor_dimension":4,
        "event_dimension":27,
        "total_dimension":108,
        "Hermitian":bool(np.linalg.norm(Q-Q.conj().T)<1e-10),
        "square_identity":"Q^2 = I4 tensor (4 L + D^2)",
        "first_order_reason":"each B_i has linear small-momentum symbol; the scalar D is their pseudoscalar product",
      },
      "pin_cover":{
        "projective_group":"W(D3) ~= S4",
        "projective_order":24,
        "lift_order":48,
        "fiber_size":2,
        "center_order":2,
        "simple_roots":[list(x) for x in roots],
        "reflection_lifts":"S_alpha = i Gamma_0 Gamma(alpha), with Gamma_0 anticommuting with Gamma_1..Gamma_3",
        "reflection_squares":"+I4",
        "order_census":{str(k):v for k,v in sorted(census.items())},
        "matches_Pass10951_GL23_order_census":True,
        "cover_reading":"same 2^+S4 lift type as the committed GL(2,3) clock cover",
      },
      "covariance":{
        "all_24_event_orthogonal_images_checked":True,
        "all_48_combined_pin_lifts_checked":True,
        "maximum_event_integer_residual":max_event,
        "maximum_Q_covariance_residual":max_q,
      },
      "two_component_firewall":"Ordinary conjugation on a single 2x2 complex Pauli algebra acts through SO(3,C), so it cannot realize the determinant-minus reflections present in this W(D3) realization. The graded 4-component carrier implements the required full O(3)-type action by ordinary conjugation.",
      "theorem":"The 27-event cubic arrow factors exactly into three vector finite differences. Their graded Clifford contraction is a Hermitian first-order operator Q with Q^2=4L+D^2. The actual W(D3) event symmetry lifts to a 48-element Pin+ cover with the same order census as repo GL(2,3)=2^+S4, and Q is invariant under all combined lifts.",
      "boundary":"Finite Pin/Clifford covariance only. Q is not yet a physical Dirac Hamiltonian, its coefficients and time interpretation are not selected dynamically, and no continuum Lorentzian metric or gravity is derived."
    }
    OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":out["status"],"cover":48,"rankD":out["vector_differences"]["D_rank"],"Q_residual":max_q},indent=2))

if __name__=="__main__": main()
