#!/usr/bin/env python3
"""Pass 11054: E6 cubic Jacobians are exact skew generators with orthogonal finite flow."""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_e6_cubic_fourier54_alignment as ALIGN

OUT=ROOT/"data/w33_pass11054_cubic_jacobian_skew_orthogonal_flow.json"

def rank_mod(A,p):
    M=np.array(A,dtype=np.int64)%p
    n,m=M.shape; r=0
    for c in range(m):
        z=next((i for i in range(r,n) if M[i,c]%p),None)
        if z is None: continue
        M[[r,z]]=M[[z,r]]
        iv=pow(int(M[r,c]%p),-1,p); M[r]=(M[r]*iv)%p
        for i in range(r+1,n):
            q=int(M[i,c]%p)
            if q: M[i]=(M[i]-q*M[r])%p
        r+=1
        if r==n: break
    return r

def basis_jacobians(records):
    out=[]
    for u in range(81):
        D=np.zeros((81,81),dtype=np.int64)
        for uu,x,o,c in records:
            if uu==u: D[o,x]+=c
        out.append(D)
    return out

def diagonal_backgrounds():
    bridge=json.loads((ROOT/"data/w33_e6id_current_h27_gauge_bridge.json").read_text())
    h={int(i):tuple(map(int,x)) for i,x in bridge["maps"]["e6id_to_current_H27_address"].items()}
    coords=[(h[eid],phase) for eid in range(27) for phase in range(3)]
    return {
      "center_plus_external":[1+((x[2]+phase)%3) for x,phase in coords],
      "center_minus_external":[1+((x[2]-phase)%3) for x,phase in coords],
    }

def combine(J,v):
    D=np.zeros((81,81),dtype=np.int64)
    for a,A in zip(v,J): D+=int(a)*A
    return D

def payload():
    records=ALIGN.ordered_records(); assert len(records)==1620
    J=basis_jacobians(records)
    assert all(np.array_equal(A.T,-A) for A in J)
    assert all(np.count_nonzero(np.triu(A,1))==10 for A in J)
    assert all(sorted(np.count_nonzero(A,axis=1).tolist())==[0]*61+[1]*20 for A in J)

    upper=np.triu_indices(81,1)
    span=np.array([A[upper] for A in J],dtype=np.int64).T
    assert rank_mod(span,103)==rank_mod(span,109)==81

    for A in J:
        rows=np.where(np.count_nonzero(A,axis=1)>0)[0]
        B=A[np.ix_(rows,rows)]
        assert rank_mod(B,103)==20 and rank_mod(B,109)==20

    diag={}
    for name,v in diagonal_backgrounds().items():
        D=combine(J,v)
        assert np.array_equal(D.T,-D)
        ranks=[rank_mod(D,p) for p in (103,109)]
        assert ranks==[66,66]
        diag[name]={"rank_mod_103":66,"rank_mod_109":66,"kernel_dimension":15}

    out={
      "schema":"w33.pass11054.cubic-jacobian-skew-orthogonal-flow.v1",
      "status":"PASS_E6_CUBIC_JACOBIANS_FORM_AN_81D_FAMILY_OF_EXACT_SKEW_GENERATORS",
      "headline":"Every E6 cubic Jacobian D_v is a real skew-symmetric 81x81 matrix. The 81 coordinate generators are linearly independent; each is exactly ten disjoint 2x2 rotation blocks, hence rank 20. Therefore exp(t D_v) lies in SO(81) for every real t, and the Cayley transform is orthogonal whenever defined.",
      "basis_jacobians":{
        "count":81,"independent_span_dimension":81,
        "edges_per_generator":10,"active_rows_per_generator":20,
        "zero_rows_per_generator":61,"exact_rank_each":20,
        "split_prime_span_ranks":{"103":81,"109":81}
      },
      "diagonal_weld_generators":diag,
      "finite_flow":{
        "exponential":"U_v(t)=exp(t D_v)",
        "orthogonality":"U_v(t)^T U_v(t)=I",
        "determinant":"+1",
        "reason":"D_v^T=-D_v, so exp(tD_v)^T=exp(-tD_v)",
        "all_ranks_even":True
      },
      "boundary":"This proves a finite orthogonal/unitary flow for the algebraic cubic Jacobian. It does not establish that a laboratory Hamiltonian implements D_v with the required coupling strength or isolation.",
      "parents":["analysis/w33_e6_cubic_fourier54_alignment.py","data/w33_e6_cubic_diagonal_phase_weld.json"],
      "checks":{
        "all_81_basis_generators_skew":True,
        "all_81_basis_generators_are_ten_matchings":True,
        "all_81_basis_generators_rank20":True,
        "basis_generator_span_rank81":True,
        "both_diagonal_welds_rank66":True,
        "orthogonal_exponential_follows_exactly":True
      }
    }
    return out

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT); a=ap.parse_args()
    p=payload(); text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or a.output.read_text()!=text: raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(text)
    print(json.dumps({"status":p["status"],"span":81,"basis_rank":20},sort_keys=True))
if __name__=="__main__": raise SystemExit(main())
