#!/usr/bin/env python3
"""Pass 11047: exact 54D phase-weld tangent compiler onto the dressed S2+L target."""
from __future__ import annotations
import argparse, hashlib, itertools, json, sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_e6_cubic_fourier54_alignment as ALIGN

OUT=ROOT/"data/w33_pass11047_phase_weld_exact_54_right_inverse.json"


def mm(A,B,p): return (A@B)%p


def mpow(A,n,p):
    R=np.eye(A.shape[0],dtype=np.int64); B=A.copy()
    while n:
        if n&1: R=mm(R,B,p)
        B=mm(B,B,p); n//=2
    return R


def inv_rank(A,p,want_inverse=False):
    A=A.copy()%p; n,m=A.shape
    aug=np.concatenate([A,np.eye(n,dtype=np.int64)],axis=1) if want_inverse else A
    r=0; piv=[]
    for c in range(m):
        z=next((i for i in range(r,n) if aug[i,c]%p),None)
        if z is None: continue
        aug[[r,z]]=aug[[z,r]]
        iv=pow(int(aug[r,c]%p),-1,p); aug[r]=(aug[r]*iv)%p
        for i in range(n):
            if i!=r and aug[i,c]%p:
                q=int(aug[i,c]%p); aug[i]=(aug[i]-q*aug[r])%p
        piv.append(c); r+=1
        if r==n: break
    inv=None
    if want_inverse:
        assert n==m and r==n
        inv=aug[:,m:]
    return r,piv,inv


def det_mod(A,p):
    M=A.copy()%p; n,m=M.shape; assert n==m
    det=1
    for c in range(n):
        z=next((i for i in range(c,n) if M[i,c]%p),None)
        if z is None: return 0
        if z!=c:
            M[[c,z]]=M[[z,c]]; det=(-det)%p
        q=int(M[c,c]%p); det=det*q%p
        iv=pow(q,-1,p); M[c]=(M[c]*iv)%p
        for i in range(c+1,n):
            a=int(M[i,c]%p)
            if a: M[i]=(M[i]-a*M[c])%p
    return int(det%p)


def root_omega(p):
    roots=[x for x in range(2,p) if (x*x+x+1)%p==0]
    assert len(roots)==2
    return min(roots)
def full_fourier_matrix(p,e6_to_h):
    w=root_omega(p)
    X=np.array([[0,0,1],[1,0,0],[0,1,0]],dtype=np.int64)
    Z=np.diag([1,w,pow(w,2,p)]).astype(np.int64)

    def rho(h,power):
        a,b,c=h
        R=mm(mpow(Z,(power*a)%3,p),mpow(X,b,p),p)
        return (pow(w,(power*c)%3,p)*R)%p

    cols=[]; labels=[]
    # S1(t,r,i)
    for t,r,i in itertools.product(range(3),repeat=3):
        col=[]
        for eid in range(27):
            M=rho(e6_to_h[eid],1)
            for phase in range(3):
                col.append(int(M[i,r])*pow(w,t*phase,p)%p)
        cols.append(col); labels.append(("S1",t,r,i))
    # S2(t,r,i)
    for t,r,i in itertools.product(range(3),repeat=3):
        col=[]
        for eid in range(27):
            M=rho(e6_to_h[eid],2)
            for phase in range(3):
                col.append(int(M[i,r])*pow(w,t*phase,p)%p)
        cols.append(col); labels.append(("S2",t,r,i))
    # L(u,v,t), grouped by external character t as in the compiler.
    for t in range(3):
      for u in range(3):
       for v in range(3):
        col=[]
        for eid in range(27):
            a,b,c=e6_to_h[eid]
            for phase in range(3):
                col.append(pow(w,u*a+v*b+t*phase,p))
        cols.append(col); labels.append(("L",u,v,t))
    F=np.array(cols,dtype=np.int64).T%p
    assert F.shape==(81,81)
    r,_,Finv=inv_rank(F,p,True)
    assert r==81
    return F,Finv,labels
def jacobian(records,v,p):
    D=np.zeros((81,81),dtype=np.int64)
    for u,x,o,c in records:
        D[o,x]=(D[o,x]+c*v[u])%p
    return D


def right_inverse_certificate(p,e6_to_h,records,v):
    F,Finv,labels=full_fourier_matrix(p,e6_to_h)
    D=jacobian(records,v,p)
    C=mm(Finv,D,p)
    Q=C[27:81,:]  # canonical dressed retyping target S2(27)+L(27)
    rank,piv,_=inv_rank(Q,p,False)
    assert rank==54 and len(piv)==54
    B=Q[:,piv]
    det=det_mod(B,p); assert det!=0
    _,_,Binv=inv_rank(B,p,True)
    assert np.array_equal(mm(B,Binv,p),np.eye(54,dtype=np.int64)%p)

    R=np.zeros((81,54),dtype=np.int64)
    for j,col in enumerate(piv): R[col,:]=Binv[j,:]
    assert np.array_equal(mm(Q,R,p),np.eye(54,dtype=np.int64)%p)

    # Split the target: both dressed non-scalar sectors are individually complete.
    assert inv_rank(Q[:27,:],p,False)[0]==27
    assert inv_rank(Q[27:,:],p,False)[0]==27

    digest=hashlib.sha256(Binv.astype(np.int64).tobytes()).hexdigest()
    return {
      "prime":p,"omega":root_omega(p),"quotient_rank":rank,
      "pivot_columns":piv,"pivot_minor_det_mod_p":det,
      "right_inverse_sha256":digest,
      "S2_projection_rank":27,"L_projection_rank":27,
      "right_inverse_verified":True,
      "target_labels_first":list(labels[27]),
      "target_labels_last":list(labels[-1]),
    }
def payload():
    bridge=json.loads((ROOT/"data/w33_e6id_current_h27_gauge_bridge.json").read_text())
    e6_to_h={int(i):tuple(map(int,h))
              for i,h in bridge["maps"]["e6id_to_current_H27_address"].items()}
    records=ALIGN.ordered_records()
    coords=[(e6_to_h[eid],phase) for eid in range(27) for phase in range(3)]
    backgrounds={
      "center_plus_external":[1+((h[2]+phase)%3) for h,phase in coords],
      "center_minus_external":[1+((h[2]-phase)%3) for h,phase in coords],
    }
    cert={}
    for name,v in backgrounds.items():
        cert[name]=[right_inverse_certificate(p,e6_to_h,records,v) for p in (103,109)]

    # The deterministic pivot set should be stable across both split primes.
    pivot_stability={}
    for name,rows in cert.items():
        same=rows[0]["pivot_columns"]==rows[1]["pivot_columns"]
        pivot_stability[name]=same
        assert same

    dressed=json.loads((ROOT/"data/w33_pass11046_full_k81_dressed_equivariant_compiler.json").read_text())
    prior=json.loads((ROOT/"data/w33_e6_cubic_diagonal_phase_weld.json").read_text())
    assert dressed["why_54_retypings"]["total_retyped"]==54
    assert prior["weld_consequence"]["required_retyped_dimension"]==54

    checks={
      "both_diagonal_weld_orientations":set(cert)=={"center_plus_external","center_minus_external"},
      "rank54_all_four_replays":all(r["quotient_rank"]==54 for rows in cert.values() for r in rows),
      "same_54_pivot_columns_at_103_and_109":all(pivot_stability.values()),
      "nonzero_54x54_minor_all_replays":all(r["pivot_minor_det_mod_p"] for rows in cert.values() for r in rows),
      "explicit_right_inverse_all_replays":all(r["right_inverse_verified"] for rows in cert.values() for r in rows),
      "S2_and_L_each_rank27":all(r["S2_projection_rank"]==r["L_projection_rank"]==27 for rows in cert.values() for r in rows),
      "matches_dressed_54_target":dressed["why_54_retypings"]["total_retyped"]==54,
    }
    assert all(checks.values())
    return {
      "schema":"w33.pass11047.phase-weld-exact-54-right-inverse.v1",
      "status":"PASS_DIAGONAL_CUBIC_WELD_HAS_EXPLICIT_54D_TANGENT_RIGHT_INVERSE_ON_THE_DRESSED_COMPILER_TARGET",
      "headline":(
        "The diagonal E6 cubic weld now lands on the exact target identified by the "
        "minimal commutant dressing. In the full K-Fourier basis S1+S2+L, the "
        "Jacobian quotient onto S2(27)+L(27) has rank 54. A deterministic set of "
        "54 root-coordinate input columns gives a nonzero 54x54 minor at both split "
        "Eisenstein primes 103 and 109, and an explicit modular right inverse is "
        "verified for both diagonal orientations."
      ),
      "target":{
        "compatible_sector":"S1, dimension 27",
        "dressed_retyping_sector":"S2 + L, dimensions 27 + 27 = 54",
        "Pass11046_identification":{
          "S2":"latent V_omega tensor internal V_omega -> 3 V_omega2",
          "L":"latent V_omega2 tensor internal V_omega -> sum_9 one-dimensional characters",
        },
      },
      "certificates":cert,
      "pivot_stability":pivot_stability,
      "compiler_consequence":(
        "At the finite tangent level the same cubic phase weld that previously "
        "surjected onto an abstract 54D quotient now surjects onto the concrete "
        "non-scalar latent sectors required to regularize H27. The stable 54-column "
        "minor supplies a fixed set of matter perturbation directions from which a "
        "right inverse can be constructed. Thus the representation target and the "
        "cubic tangent mechanism are aligned, rather than merely dimension-matched."
      ),
      "exactness_argument":(
        "All Fourier and Jacobian entries lie in Z[omega]. The same 54-column minor "
        "is nonzero after reduction at the two split primes 103 and 109, so that "
        "minor is nonzero over Q(omega). Therefore the quotient map has a 54D "
        "characteristic-zero right inverse on those selected columns."
      ),
      "parents":[
        "data/w33_pass11046_full_k81_dressed_equivariant_compiler.json",
        "data/w33_e6_cubic_diagonal_phase_weld.json",
        "data/w33_e6_cubic_fourier54_alignment.json",
      ],
      "boundary":(
        "This is a tangent-level algebraic compiler. A right inverse of the cubic "
        "Jacobian is not yet a finite-time unitary pulse, a stable vacuum, or a "
        "fault-tolerant gate. Exponentiating the selected 54 directions while "
        "preserving the FI orientation is the next physical synthesis problem."
      ),
      "checks":checks,
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT)
    a=ap.parse_args(); p=payload()
    text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or a.output.read_text()!=text: raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(text)
    print(json.dumps({
      "status":p["status"],
      "plus_pivots":p["certificates"]["center_plus_external"][0]["pivot_columns"][:8],
      "minus_pivots":p["certificates"]["center_minus_external"][0]["pivot_columns"][:8],
      "stable":p["pivot_stability"]},sort_keys=True))


if __name__=="__main__":
    raise SystemExit(main())
