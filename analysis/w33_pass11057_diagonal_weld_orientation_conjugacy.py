#!/usr/bin/env python3
"""Pass 11057: the two diagonal FI welds are inverse-conjugate finite gates."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_pass11056_diagonal_weld_ten_triangle_factors as FCT
import w33_pass11047_phase_weld_exact_54_right_inverse as PREV

OUT=ROOT/"data/w33_pass11057_diagonal_weld_orientation_conjugacy.json"

def phase_reversal(coords):
    idx={x:i for i,x in enumerate(coords)}
    P=np.zeros((81,81),dtype=np.int64)
    for j,(h,phase) in enumerate(coords): P[idx[(h,(-phase)%3)],j]=1
    return P

def payload():
    coords,vp,Dp=FCT.diagonal("plus")
    coords2,vm,Dm=FCT.diagonal("minus"); assert coords2==coords
    P=phase_reversal(coords); I=np.eye(81,dtype=np.int64)
    assert np.array_equal(P@P,I)
    assert np.array_equal(Dm,-P@Dp@P.T)
    assert np.array_equal(Dm,P@Dp.T@P.T)

    replay=[]
    for p in (103,109):
        Ip=I%p; Pp=P%p; Ap=Dp%p; Am=Dm%p
        _,_,Ainv=PREV.inv_rank((Ip-Ap)%p,p,True)
        _,_,Minv=PREV.inv_rank((Ip-Am)%p,p,True)
        Cp=((Ip+Ap)%p@Ainv)%p
        Cm=((Ip+Am)%p@Minv)%p
        assert np.array_equal(Cm,(Pp@Cp.T@Pp.T)%p)
        replay.append({"prime":p,"t":1,"inverse_conjugacy_verified":True})

    return {
      "schema":"w33.pass11057.diagonal-weld-orientation-conjugacy.v1",
      "status":"PASS_EXTERNAL_PHASE_REVERSAL_EXCHANGES_THE_TWO_DIAGONAL_WELDS_AS_INVERSE_ORTHOGONAL_GATES",
      "headline":"Let R reverse only the external qutrit phase p->-p. The two diagonal cubic generators obey D_minus=-R D_plus R^T exactly. Consequently their finite Cayley gates satisfy C_minus(t)=R C_plus(-t) R^T=R C_plus(t)^T R^T. The two FI orientations are therefore exact inverse-conjugate gates, not merely equal-rank alternatives.",
      "generator_identity":{
        "phase_reversal":"R:(h,p)->(h,-p)","R_squared":"I",
        "relation":"D_minus=-R D_plus R^T","transpose_form":"D_minus=R D_plus^T R^T"
      },
      "finite_gate_identity":{
        "relation":"C_minus(t)=R C_plus(-t) R^T",
        "orthogonal_form":"C_minus(t)=R C_plus(t)^T R^T",
        "interpretation":"orientation reversal is external-phase conjugation plus time/generator reversal"
      },
      "split_prime_replays":replay,
      "boundary":"This identifies the exact algebraic relation between the two FI orientations. It does not select one orientation dynamically or infer spontaneous chirality.",
      "parents":["data/w33_pass11055_cayley_finite_54_compiler.json","data/w33_e6_cubic_diagonal_phase_weld.json"],
      "checks":{"phase_reversal_involution":True,"generator_inverse_conjugacy":True,"finite_cayley_inverse_conjugacy":True}
    }

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args()
    p=payload();text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or a.output.read_text()!=text:raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"relation":"inverse-conjugate"},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
