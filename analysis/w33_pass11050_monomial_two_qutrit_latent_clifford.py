#!/usr/bin/env python3
"""Pass 11050: induced A9 is a two-qutrit monomial Clifford representation."""
from __future__ import annotations
import argparse,itertools,json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/w33_pass11050_monomial_two_qutrit_latent_clifford.json"

def root_omega(p):
    r=[x for x in range(2,p) if (x*x+x+1)%p==0]
    assert len(r)==2
    return min(r)

def perm9(fn):
    M=np.zeros((9,9),dtype=np.int64)
    for b,c in itertools.product(range(3),repeat=2):
        bp,cp=fn(b,c);M[3*bp+cp,3*b+c]=1
    return M

def mpow(A,n,p):
    R=np.eye(A.shape[0],dtype=np.int64);B=A.copy()%p
    while n:
        if n&1:R=(R@B)%p
        B=(B@B)%p;n//=2
    return R
def prime_packet(p):
    w=root_omega(p);I=np.eye(9,dtype=np.int64)
    AX=perm9(lambda b,c:((b+1)%3,c))%p
    AC=perm9(lambda b,c:(b,(c+1)%3))%p
    S=perm9(lambda b,c:(b,(c+b)%3))%p
    Zb=np.diag([pow(w,b,p) for b in range(3) for c in range(3)]).astype(np.int64)
    Zc=np.diag([pow(w,c,p) for b in range(3) for c in range(3)]).astype(np.int64)
    AXc=AC
    assert np.array_equal(mpow(S,3,p),I)
    assert np.array_equal((S@AX)%p,(AC@AX@S)%p)
    assert np.array_equal((S@AC)%p,(AC@S)%p)

    Sinv=mpow(S,2,p)
    assert np.array_equal((S@AX@Sinv)%p,(AX@AC)%p)
    assert np.array_equal((S@AC@Sinv)%p,AC)
    assert np.array_equal((S@Zb@Sinv)%p,Zb)
    target=(mpow(Zb,2,p)@Zc)%p
    assert np.array_equal((S@Zc@Sinv)%p,target)
    theta_rows=[]
    for theta in range(3):
        AZ=(pow(w,theta,p)*S)%p
        assert np.array_equal(mpow(AZ,3,p),I)
        assert np.array_equal((AZ@AX)%p,(AC@AX@AZ)%p)
        assert np.array_equal((AZ@AC)%p,(AC@AZ)%p)

        max_mismatch=0
        traces=[]
        for a,b,c in itertools.product(range(3),repeat=3):
            R=(mpow(AZ,a,p)@mpow(AX,b,p)@mpow(AC,c,p))%p
            tr=int(np.trace(R)%p)
            line=sum(pow(w,theta*a+s*b,p) for s in range(3))%p
            vv=0
            if a==0 and b==0:
                vv=(3*pow(w,c,p)+3*pow(w,2*c,p))%p
            theo=(line+vv)%p
            max_mismatch=max(max_mismatch,(tr-theo)%p)
            assert tr==theo
            traces.append([a,b,c,tr])
        theta_rows.append({
          "theta":theta,
          "dual_affine_line":[[theta,s] for s in range(3)],
          "character_matches_induced_A9":True,
          "trace_rows":traces,
        })
    return {
      "prime":p,"omega":w,
      "theta_rows":theta_rows,
      "generator_relations":{
        "AZ_AX_equals_AC_AX_AZ":True,
        "AC_central":True,
        "orders":[3,3,3],
      },
      "clifford_conjugation":{
        "S_Xb_Sdag":"Xb Xc",
        "S_Xc_Sdag":"Xc",
        "S_Zb_Sdag":"Zb",
        "S_Zc_Sdag":"Zb^-1 Zc",
        "verified":True,
      },
    }

def payload():
    cert=[prime_packet(p) for p in (103,109)]
    checks={
      "all_three_theta_induced_characters":all(
          all(r["character_matches_induced_A9"] for r in c["theta_rows"]) for c in cert),
      "H27_relations_exact":all(c["generator_relations"]["AZ_AX_equals_AC_AX_AZ"] for c in cert),
      "shear_normalizes_two_qutrit_Paulis":all(c["clifford_conjugation"]["verified"] for c in cert),
      "only_permutation_plus_mu3_phase_needed":True,
    }
    assert all(checks.values())
    return {
      "schema":"w33.pass11050.monomial-two-qutrit-latent-clifford.v1",
      "status":"PASS_MINIMAL_LATENT_A9_HAS_A_TWO_QUTRIT_MONOMIAL_CLIFFORD_REALIZATION",
      "headline":(
        "Choose the induced direction L=<Z>. On the 9 cosets use coordinates "
        "|b,c> in F3^2. The latent H27 generators can be represented as "
        "X_A: b->b+1, C_A: c->c+1, and Z_A: c->c+b times a theta-dependent "
        "global mu3 phase. These are translations and a qutrit SUM shear, so the "
        "minimal A9 dressing is realizable inside the two-qutrit Clifford group."
      ),
      "coset_coordinates":{
        "basis":"|b,c>, (b,c) in F3^2",
        "X_A":"|b,c> -> |b+1,c>",
        "C_A":"|b,c> -> |b,c+1>",
        "Z_A_theta":"omega^theta |b,c> -> omega^theta |b,c+b>",
      },
      "gate_compiler":{
        "X_A":"single-qutrit X on latent coordinate b",
        "C_A":"single-qutrit X on latent coordinate c",
        "Z_A":"SUM(b -> c) plus a global mu3 phase",
        "generic_U9_required":False,
        "entangling_primitive_count_per_Z_generator":1,
      },
      "Clifford_normalizer":{
        "SUM_conjugation":{
          "X_b":"X_b X_c","X_c":"X_c",
          "Z_b":"Z_b","Z_c":"Z_b^-1 Z_c",
        },
        "consequence":"the latent H27 action normalizes the two-qutrit Weyl group",
      },
      "split_prime_certificates":cert,
      "parents":[
        "data/w33_pass11048_induced_h27_hesse_latent_modules.json",
        "data/w33_pass11045_explicit_latent_h27_clebsch_gordan.json",
      ],
      "boundary":(
        "The abstract 9D multiplicity subsystem is given an explicit C3 tensor C3 "
        "coordinate in this realization. That factorization is algebraically valid "
        "and Clifford-synthesizable, but a particular photonic device must still "
        "identify its physical modes with these latent coordinates."
      ),
      "checks":checks,
    }

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args();p=payload()
    text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or a.output.read_text()!=text:raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"generic_U9":False},sort_keys=True))

if __name__=="__main__":raise SystemExit(main())
