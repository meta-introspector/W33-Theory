#!/usr/bin/env python3
"""Pass 11558: embed the new event Pin carrier into the committed Albert Peirce-16.

Pass10961 gives exact rational Cl(9) gamma matrices on the Peirce 16.  The
first three, together with a grading built from five orthogonal remaining
gammas, generate a Cl(4,C)=M4(C) subalgebra.  Hence the 16-complex-dimensional
Peirce carrier is four copies of the 4D graded event-spinor module from
Pass11557.  The same Pin+ reflection lifts act on this actual Albert carrier.

This does NOT contradict Pass10959: the representation constructed here has a
different order-eight clock spectrum from the old Albert half-turn phase lift
whose 16D extension was proved impossible.
"""
from __future__ import annotations
import itertools, json
from collections import Counter, deque
from pathlib import Path
import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/PART_W33_PASS11558_ALBERT_PIN_EVENT_SPINOR_EMBEDDING.json"

def parse_matrix(rows):
    return sp.Matrix([[sp.sympify(x) for x in row] for row in rows])

def mrank(mats):
    cols=[sp.Matrix(M).reshape(M.rows*M.cols,1) for M in mats]
    return sp.Matrix.hstack(*cols).rank()

def nkey(M):
    A=np.array(M,dtype=complex)
    return tuple(np.round(A.real,10).ravel())+tuple(np.round(A.imag,10).ravel())

def main():
    cert=json.loads((ROOT/"data/w33_pass10961_albert_clifford9_gammas.json").read_text())
    g9=[parse_matrix(x) for x in cert["gamma9"]]
    I16=sp.eye(16)
    for i in range(9):
        for j in range(9):
            rhs=(2 if i==j else 0)*I16
            assert g9[i]*g9[j]+g9[j]*g9[i]==rhs

    gam=g9[:3]
    # Five orthogonal external gammas give an involutory grading that
    # anticommutes with the selected Cl(3).
    grade=sp.eye(16)
    for g in g9[3:8]: grade=grade*g
    assert grade*grade==I16
    for g in gam: assert grade*g+g*grade==sp.zeros(16)

    # The four Clifford generators span the full 16-dimensional Cl4 algebra.
    four=gam+[grade]
    basis=[]
    for mask in range(16):
        M=I16
        for k in range(4):
            if (mask>>k)&1: M=M*four[k]
        basis.append(M)
    assert mrank(basis)==16

    roots=[(1,-1,0),(0,1,-1),(0,1,1)]
    refs=[]; lifts=[]
    for a in roots:
        av=sp.Matrix(a); R=sp.eye(3)-av*av.T
        p=sum((sp.Integer(a[j])*gam[j] for j in range(3)),sp.zeros(16))/sp.sqrt(2)
        S=sp.simplify(sp.I*grade*p)
        assert sp.simplify(S*S-I16)==sp.zeros(16)
        for k in range(3):
            lhs=sp.simplify(S*gam[k]*S)
            rhs=sum((R[j,k]*gam[j] for j in range(3)),sp.zeros(16))
            assert sp.simplify(lhs-rhs)==sp.zeros(16)
        refs.append(np.array(R.tolist(),dtype=int));lifts.append(S)

    # Enumerate the finite cover numerically while all defining Clifford
    # identities above remain exact.
    gens=[np.array(S.evalf(35),dtype=complex) for S in lifts]
    I16n=np.eye(16,dtype=complex); I3=np.eye(3,dtype=int)
    group={nkey(I16n):(I3,I16n)}
    q=deque([(I3,I16n)])
    while q:
        R,S=q.popleft()
        for rg,sg in zip(refs,gens):
            R2=R@rg;S2=S@sg;k=nkey(S2)
            if k not in group:
                group[k]=(R2,S2);q.append((R2,S2))
    assert len(group)==48
    fibers=Counter(tuple(R.ravel()) for R,S in group.values())
    assert len(fibers)==24 and set(fibers.values())=={2}

    def order(S):
        P=I16n.copy()
        for n in range(1,25):
            P=P@S
            if np.linalg.norm(P-I16n)<1e-8:return n
        raise AssertionError
    census=Counter(order(S) for R,S in group.values())
    assert census==Counter({2:13,8:12,3:8,6:8,4:6,1:1})

    # Find one order-8 lift and classify its eighth-root spectrum.
    S8=next(S for R,S in group.values() if order(S)==8)
    vals=np.linalg.eigvals(S8)
    zeta=np.exp(1j*np.pi/4)
    exps=[]
    for v in vals:
        k=min(range(8),key=lambda j:abs(v-zeta**j))
        assert abs(v-zeta**k)<1e-7
        exps.append(k)
    ec=Counter(exps)
    assert ec==Counter({1:4,3:4,5:4,7:4})

    p59=json.loads((ROOT/"data/w33_pass10959_albert_clock_extension_obstruction.json").read_text())
    assert p59["dimension16_no_go"]["target_15_extension_exists"] is False
    assert p59["dimension16_no_go"]["target_37_extension_exists"] is False
    p61=json.loads((ROOT/"data/w33_pass10961_albert_clifford10_doubled_clock.json").read_text())
    assert p61["albert_cl9"]["gamma_count"]==9

    # Event-vector covariance uses the same reflection matrices as Pass11557.
    import sys
    sys.path.insert(0,str(ROOT/"analysis"))
    import w33_pass11557_pin_equivariant_event_dirac as P57
    V,vid,Bs,A,D,L=P57.event_operators()
    max_event=0
    for R in refs:
        perm=[]
        for x in V:
            y=tuple(int(sum(R[a,b]*x[b] for b in range(3))%3) for a in range(3))
            perm.append(vid[y])
        P=np.zeros((27,27),dtype=np.int64)
        for x,y in enumerate(perm):P[y,x]=1
        for k in range(3):
            lhs=P@Bs[k]@P.T
            rhs=sum(int(R[j,k])*Bs[j] for j in range(3))
            max_event=max(max_event,int(np.max(np.abs(lhs-rhs))))
    assert max_event==0

    out={
      "schema":"w33.pass11558.albert-pin-event-spinor-embedding.v1",
      "status":"PASS_EVENT_PIN_CARRIER_EMBEDS_FOURFOLD_IN_ALBERT16","pass":11558,
      "albert_clifford":{
        "source":"Pass10961 exact rational Cl(9) on the Peirce 16",
        "selected_gamma_count":3,
        "grading":"Gamma4*Gamma5*Gamma6*Gamma7*Gamma8 (five unused Cl9 vectors)",
        "grading_square":"+I16",
        "generated_Cl4_algebra_dimension":16,
        "module_reading":"over C, Cl4=M4(C); the 16D Peirce carrier is four copies of the irreducible 4D graded event-spinor module",
      },
      "pin_action":{
        "lift_order":48,"projective_order":24,"fiber_size":2,
        "order_census":{str(k):v for k,v in sorted(census.items())},
        "order8_eighth_root_multiplicities":{str(k):v for k,v in sorted(ec.items())},
        "order8_reading":"four copies of the faithful 4D central-odd GL(2,3) spectrum {1,3,5,7}",
        "same_reflection_covariance_as_Pass11557":True,
      },
      "clock_no_go_reconciliation":{
        "Pass10959_single16_old_phase_lift_no_go":True,
        "reason_no_contradiction":"Pass10959 forbids the two old Albert half-turn phase spectra {1,5}^8 or {3,7}^8. The new Clifford-reflection GL(2,3) action has {1,3,5,7} each multiplicity 4 and is a different representation.",
        "Pass10961_old_doubled_clock_grade_obstruction_still_valid":True,
      },
      "event_weld":{
        "same_three_event_Bi":True,
        "generator_event_covariance_residual":max_event,
        "Albert_event_operator":"Q_A = i sum_{i=1}^3 gamma_i^(Albert) tensor B_i",
        "equivariance":"proved on the three Pin/W(D3) generators, hence on the generated 48-element cover",
      },
      "theorem":"The exact Albert Peirce-16 already contains the four-component Pin event spinor fourfold: three rational Cl9 generators plus an internal grading form Cl4(C), and their reflection lifts generate the same 48-element 2^+S4 cover acting covariantly on the event differences. This supplies a native Albert carrier for the new first-order Q without identifying it with the old frozen order-eight clock operator.",
      "boundary":"Finite complex Clifford-module embedding. The choice of three of nine Albert spatial gamma directions is not dynamically selected; no physical family assignment, Lorentzian spacetime, fermion mass, or gravity follows."
    }
    OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":out["status"],"cover":48,"Cl4_dim":16,"order8":dict(ec)},indent=2))

if __name__=="__main__": main()
