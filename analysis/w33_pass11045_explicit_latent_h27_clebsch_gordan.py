#!/usr/bin/env python3
"""Pass 11045: explicit 9x9 latent H27 action and 27x27 Clebsch-Gordan compiler."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/w33_pass11045_explicit_latent_h27_clebsch_gordan.json"


def mm(A,B,p): return (A@B)%p


def mpow(A,n,p):
    R=np.eye(A.shape[0],dtype=np.int64); B=A.copy()
    while n:
        if n&1: R=mm(R,B,p)
        B=mm(B,B,p); n//=2
    return R


def kron(A,B,p): return np.kron(A,B)%p


def blockdiag(blocks,p):
    n=sum(B.shape[0] for B in blocks)
    R=np.zeros((n,n),dtype=np.int64); k=0
    for B in blocks:
        d=B.shape[0]; R[k:k+d,k:k+d]=B%p; k+=d
    return R


def rank_det(A,p):
    M=A.copy()%p; n,m=M.shape; r=0; det=1
    for c in range(m):
        z=next((i for i in range(r,n) if M[i,c]%p),None)
        if z is None: continue
        if z!=r:
            M[[r,z]]=M[[z,r]]; det=(-det)%p
        piv=int(M[r,c]%p); det=det*piv%p
        iv=pow(piv,-1,p); M[r]=(M[r]*iv)%p
        for i in range(r+1,n):
            q=int(M[i,c]%p)
            if q: M[i]=(M[i]-q*M[r])%p
        r+=1
        if r==n: break
    return r,(det%p if n==m and r==n else 0)


def root_omega(p):
    roots=[x for x in range(2,p) if (x*x+x+1)%p==0]
    assert len(roots)==2
    return min(roots)


def qutrit_generators(p,w,power=1):
    X=np.array([[0,0,1],[1,0,0],[0,1,0]],dtype=np.int64)
    Z=np.diag([1,w,pow(w,2,p)]).astype(np.int64)
    if power==1: return X,Z
    return X,mpow(Z,power,p)


def latent_generators(p,w):
    X1,Z1=qutrit_generators(p,w,1)
    X2,Z2=qutrit_generators(p,w,2)
    # character line chi_00, chi_10, chi_20:
    AX=blockdiag([
      np.array([[1]]),np.array([[1]]),np.array([[1]]),X1,X2],p)
    AZ=blockdiag([
      np.array([[1]]),np.array([[w]]),np.array([[pow(w,2,p)]]),Z1,Z2],p)
    return AX,AZ
def target_generators(p,w):
    X1,Z1=qutrit_generators(p,w,1)
    X2,Z2=qutrit_generators(p,w,2)
    TX=blockdiag([X1,X1,X1,X2,X2,X2],p)
    TZ=blockdiag([Z1,Z1,Z1,Z2,Z2,Z2],p)
    charsX=[]; charsZ=[]
    for u in range(3):
        for v in range(3):
            charsX.append(pow(w,v,p)); charsZ.append(pow(w,u,p))
    TX=blockdiag([TX,np.diag(charsX).astype(np.int64)],p)
    TZ=blockdiag([TZ,np.diag(charsZ).astype(np.int64)],p)
    return TX,TZ


def cg_matrix(p,w):
    # Rows are target Fourier coordinates; columns are product basis (m,y).
    T=np.zeros((27,27),dtype=np.int64)
    row=0
    # 3 scalar-character blocks tensor V -> 3 standard V copies.
    for j in range(3):
        for i in range(3):
            y=(i-j)%3
            T[row,j*3+y]=1; row+=1
    # V tensor V -> 3 Vbar copies, labelled by difference d=y-x.
    for d in range(3):
        for i in range(3):
            x=(i-2*d)%3; y=(x+d)%3
            T[row,(3+x)*3+y]=1; row+=1
    # Vbar tensor V -> all nine one-dimensional characters chi_(u,v).
    for u in range(3):
        for v in range(3):
            for x in range(3):
                y=(x+u)%3
                T[row,(6+x)*3+y]=pow(w,v*x,p)
            row+=1
    assert row==27
    return T


def prime_certificate(p):
    w=root_omega(p)
    AX,AZ=latent_generators(p,w)
    X1,Z1=qutrit_generators(p,w,1)
    I9=np.eye(9,dtype=np.int64); I3=np.eye(3,dtype=np.int64)
    AC=blockdiag([
      np.array([[1]]),np.array([[1]]),np.array([[1]]),
      (w*np.eye(3,dtype=np.int64))%p,
      (pow(w,2,p)*np.eye(3,dtype=np.int64))%p],p)
    LX=kron(AX,I3,p); LZ=kron(AZ,I3,p); LC=kron(AC,I3,p)
    EX=kron(I9,X1,p); EZ=kron(I9,Z1,p)
    EC=(w*np.eye(27,dtype=np.int64))%p
    assert np.array_equal(mm(LX,EX,p),mm(EX,LX,p))
    assert np.array_equal(mm(LZ,EZ,p),mm(EZ,LZ,p))
    UX=mm(LX,EX,p); UZ=mm(LZ,EZ,p); UC=mm(LC,EC,p)

    TX,TZ=target_generators(p,w)
    TC=blockdiag([
      (w*np.eye(9,dtype=np.int64))%p,
      (pow(w,2,p)*np.eye(9,dtype=np.int64))%p,
      np.eye(9,dtype=np.int64)],p)
    T=cg_matrix(p,w)
    r,det=rank_det(T,p)
    assert r==27 and det!=0
    assert np.array_equal(mm(T,UX,p),mm(TX,T,p))
    assert np.array_equal(mm(T,UZ,p),mm(TZ,T,p))
    assert np.array_equal(mm(T,UC,p),mm(TC,T,p))

    I27=np.eye(27,dtype=np.int64)
    # H27 relations: X^3=Z^3=C^3=1 and Z X = C X Z.
    assert np.array_equal(mpow(UX,3,p),I27)
    assert np.array_equal(mpow(UZ,3,p),I27)
    assert np.array_equal(mpow(UC,3,p),I27)
    assert np.array_equal(mm(UZ,UX,p),mm(UC,mm(UX,UZ,p),p))
    # The target character is the regular one on all 27 H27 elements.
    char_rows=[]
    for a in range(3):
      for b in range(3):
       for c in range(3):
        U=mm(mpow(UC,c,p),mm(mpow(UZ,a,p),mpow(UX,b,p),p),p)
        tr=int(np.trace(U)%p)
        expected=27%p if (a,b,c)==(0,0,0) else 0
        assert tr==expected
        char_rows.append([a,b,c,tr])
    return {
      "prime":p,"omega":w,"CG_rank":r,"CG_det_mod_p":int(det),
      "latent_execution_commute":True,
      "dressed_H27_relations":True,
      "CG_intertwines_X_and_Z":True,
      "regular_character_all_27_elements":True,
      "character_rows":char_rows,
    }


def payload():
    certs=[prime_certificate(p) for p in (103,109)]
    prior=json.loads((ROOT/"data/w33_pass11043_commutant_dressing_regular_h27.json").read_text())
    assert prior["status"].startswith("PASS_9D_LATENT_DRESSING")
    checks={
      "two_split_prime_replays":len(certs)==2,
      "all_CG_ranks_27":all(c["CG_rank"]==27 for c in certs),
      "all_CG_determinants_nonzero":all(c["CG_det_mod_p"] for c in certs),
      "all_generator_intertwiners_exact_mod_p":all(c["CG_intertwines_X_and_Z"] for c in certs),
      "latent_action_commutes_with_old_execution":all(c["latent_execution_commute"] for c in certs),
      "regular_character_verified":all(c["regular_character_all_27_elements"] for c in certs),
    }
    assert all(checks.values())
    return {
      "schema":"w33.pass11045.explicit-latent-h27-clebsch-gordan.v1",
      "status":"PASS_EXPLICIT_M9_COMMUTANT_ACTION_AND_27D_CLEBSCH_GORDAN_REGULARIZATION",
      "headline":(
        "The latent representation is now explicit. In the 9D multiplicity factor, "
        "A9 is block diagonal chi00 + chi10 + chi20 + V + Vbar. It commutes with "
        "the old internal-qutrit execution because it acts in the M9 commutant. "
        "A closed-form 27x27 Clebsch-Gordan transform then sends A9 tensor V "
        "to 3V + 3Vbar + sum_9 chi, i.e. the regular H27 Fourier basis."
      ),
      "latent_generators":{
        "basis_blocks":["chi00","chi10","chi20","V_omega","V_omega2"],
        "X":"diag(1,1,1,X,X)",
        "Z":"diag(1,omega,omega^2,Z,Z^2)",
        "location":"M9 multiplicity commutant",
      },
      "CG_formula":{
        "scalar_times_V":"i=j+y; permutation to three V copies",
        "V_times_V":"d=y-x, i=x+2d; permutation to three Vbar copies",
        "Vbar_times_V":"u=y-x and Fourier transform in common coordinate x gives chi_(u,v)",
        "normalization":"the last 9x9 block becomes unitary over C after the usual 1/sqrt(3) Fourier normalization",
      },
      "split_prime_certificates":certs,
      "exactness_reading":(
        "The formulas are identities over Z[omega]; the split-prime replays certify "
        "that the displayed transform is invertible and intertwines both H27 generators. "
        "The nonzero determinants at 103 and 109 also witness a nonzero determinant "
        "before reduction in Q(omega)."
      ),
      "parents":[
        "data/w33_pass11043_commutant_dressing_regular_h27.json",
        "data/w33_pass11044_minimal_latent_dimension_and_54_budget.json",
      ],
      "boundary":(
        "The Clebsch-Gordan transform is a finite unitary basis change after Fourier "
        "normalization. The certificate does not yet supply a pulse schedule or "
        "Hamiltonian for synthesizing the latent A9 generators in hardware."
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
        if not a.output.exists() or a.output.read_text()!=text:
            raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(text)
    print(json.dumps({
      "status":p["status"],
      "certs":[[c["prime"],c["CG_det_mod_p"]] for c in p["split_prime_certificates"]]
    },sort_keys=True))


if __name__=="__main__":
    raise SystemExit(main())
