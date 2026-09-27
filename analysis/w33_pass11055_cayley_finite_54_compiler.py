#!/usr/bin/env python3
"""Pass 11055: Cayley-transform the diagonal cubic weld into a finite orthogonal 54D compiler."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_e6_cubic_fourier54_alignment as ALIGN
import w33_pass11047_phase_weld_exact_54_right_inverse as PREV

OUT=ROOT/"data/w33_pass11055_cayley_finite_54_compiler.json"

def backgrounds():
    bridge=json.loads((ROOT/"data/w33_e6id_current_h27_gauge_bridge.json").read_text())
    h={int(i):tuple(map(int,x)) for i,x in bridge["maps"]["e6id_to_current_H27_address"].items()}
    coords=[(h[eid],phase) for eid in range(27) for phase in range(3)]
    return h,{
      "center_plus_external":[1+((x[2]+phase)%3) for x,phase in coords],
      "center_minus_external":[1+((x[2]-phase)%3) for x,phase in coords],
    }

def payload():
    h,bg=backgrounds(); records=ALIGN.ordered_records()
    cert={}
    for name,v in bg.items():
        rows=[]
        for p in (103,109):
            F,Finv,_=PREV.full_fourier_matrix(p,h)
            D=PREV.jacobian(records,v,p)
            I=np.eye(81,dtype=np.int64)
            tested=[]
            for t in range(1,11):
                M=(I-t*D)%p
                r,_,Minv=PREV.inv_rank(M,p,True)
                assert r==81
                C=((I+t*D)%p@Minv)%p
                assert np.array_equal((C.T@C)%p,I%p)
                assert PREV.inv_rank((C-I)%p,p,False)[0]==66
                Q=(Finv@((C-I)%p))%p
                assert PREV.inv_rank(Q[27:,:],p,False)[0]==54
                tested.append(t)
            rows.append({"prime":p,"tested_nonzero_t":tested,"denominator_rank":81,"cayley_minus_identity_rank":66,"retyped_quotient_rank":54,"orthogonal":True})
        cert[name]=rows

    out={
      "schema":"w33.pass11055.cayley-finite-54-compiler.v1",
      "status":"PASS_CAYLEY_TRANSFORM_UPGRADES_THE_TANGENT_WELD_TO_A_FINITE_ORTHOGONAL_54D_COMPILER",
      "headline":"For either diagonal cubic weld, D is real skew-symmetric of rank 66 and projects with rank 54 onto S2+L. The Cayley gate C_t=(I+tD)(I-tD)^-1 is defined for every real t, is exactly orthogonal, and for every nonzero real t has Im(C_t-I)=Im(D). Hence the full 54D dressed quotient is preserved at finite amplitude, not only to first order.",
      "theorem":{
        "gate":"C_t=(I+tD)(I-tD)^-1",
        "real_denominator_never_singular":"det(I-tD)=product_j(1+t^2 lambda_j^2)>0",
        "orthogonality":"C_t^T C_t=I",
        "image_identity":"C_t-I=2t D (I-tD)^-1",
        "kernel_identity":"ker(C_t-I)=ker(D) for t!=0",
        "rank_for_diagonal_weld":66,
        "retyped_quotient_rank_for_diagonal_weld":54
      },
      "split_prime_replays":cert,
      "compiler_consequence":{
        "finite_gate_exists_algebraically":True,
        "finite_gate_covers_all_54_retyped_directions":True,
        "tangent_only_obstruction_remaining":False,
        "gate_is_real_orthogonal_and_complex_unitary":True
      },
      "boundary":"The Cayley transform is an exact finite unitary compiler in the algebraic carrier. It is not a proof that the physical device realizes the rational Cayley map as a single constant-Hamiltonian pulse; pulse synthesis, coupling calibration, loss and fault tolerance remain open.",
      "parents":["data/w33_pass11047_phase_weld_exact_54_right_inverse.json","data/w33_pass11054_cubic_jacobian_skew_orthogonal_flow.json"],
      "checks":{
        "both_orientations":True,
        "all_real_nonzero_t_preserve_image_by_identity":True,
        "orthogonal_by_skew_symmetry":True,
        "mod103_t1_to10_rank66_projection54":True,
        "mod109_t1_to10_rank66_projection54":True
      }
    }
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args()
    p=payload();text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or a.output.read_text()!=text:raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"finite_rank":66,"quotient_rank":54},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
