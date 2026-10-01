#!/usr/bin/env python3
"""Exact survival of the sparse 54D transducer under the K81->U81 cocycle twist."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_pass11062_single_foliation_finite54_transducer as F62
import w33_pass11068_ternary_cocycle_deformation as D68
import w33_pass11067_explicit_h27_central_cocycle_weld as C67
import w33_pass11047_phase_weld_exact_54_right_inverse as P47
OUT=ROOT/"data/w33_20261001_u81_cocycle_blind_54_transducer.json"

def k4(k):
    return tuple(k[0])+(k[1],)

def payload():
    bridge=json.loads((ROOT/"data/w33_e6id_current_h27_gauge_bridge.json").read_text())
    h={int(i):tuple(map(int,x)) for i,x in bridge["maps"]["e6id_to_current_H27_address"].items()}
    rows={}
    expected={"plus":[7,8],"minus":[6,9]}
    for which in ("plus","minus"):
        coords,pairs,factors=F62.factorize(which)
        four=[k4(x) for x in coords]
        blind=[]; certified=[]
        for i,(g0,_) in enumerate(pairs):
            g=k4(g0)
            is_blind=(g[0]==0)
            if is_blind:
                assert all(C67.kappa(x[:3],g[:3])==0 for x in four)
                assert all(D68.mul(x,g,s)==D68.mul(x,g,0) for x in four for s in (1,2))
                blind.append(i)
            if i in expected[which]:
                assert is_blind and [D68.order(g,s) for s in (0,1,2)]==[3,3,3]
                ranks=[]
                for p in (103,109):
                    _,Finv,_=P47.full_fourier_matrix(p,h)
                    M=factors[i]%p;C=(Finv@M)%p
                    ranks.append({"prime":p,"raw":P47.inv_rank(M,p,False)[0],
                                  "quotient":P47.inv_rank(C[27:,:],p,False)[0],
                                  "S2":P47.inv_rank(C[27:54,:],p,False)[0],
                                  "L":P47.inv_rank(C[54:,:],p,False)[0]})
                assert all(r["raw"]==r["quotient"]==54 and r["S2"]==r["L"]==27 for r in ranks)
                certified.append({"factor":i,"generator":list(g),"ranks":ranks})
        assert blind==[6,7,8,9]
        rows[which]={"cocycle_blind_factors":blind,"perfect_survivors":certified}
    out={
      "schema":"w33.20261001.u81-cocycle-blind-54-transducer.v1",
      "status":"PASS_PERFECT_54D_TRIANGLE_TRANSDUCERS_SURVIVE_THE_NONSPLIT_U81_COCYCLE_UNCHANGED",
      "cocycle":"kappa((a,b,c),(A,B,C))=A^2*b-2*A*c mod 3",
      "criterion":"right generators with A=0 are cocycle-blind for every carrier label",
      "orientations":rows,
      "theorem":(
        "All four Pass11062 perfect single-foliation generators have first H27 coordinate A=0. "
        "Hence kappa(h,g)=0 identically, so right multiplication by g is literally the same for "
        "s=0,+1,-1. Their weighted 27-triangle matrices therefore embed unchanged as U81 Cayley "
        "factors and retain raw rank 54 and quotient ranks S2=27, L=27 over both split primes."
      ),
      "finite_gate_consequence":(
        "The same one-parameter Cayley gate C_t(F) from Pass11062 remains a sparse finite retyping "
        "gate in the nonsplit chamber law for every nonzero real t; no cocycle-dependent retuning is needed."
      ),
      "nontriviality":"Only factors 6-9 are cocycle-blind; the other six weld directions genuinely feel the class-three twist.",
      "boundary":(
        "This closes the algebraic embedding of the certified sparse 54D factors into U81. It does not "
        "derive a laboratory Hamiltonian, pulse duration, fault-tolerance threshold, or a dynamical choice of cocycle orientation."
      ),
      "parents":[
        "data/w33_pass11062_single_foliation_finite54_transducer.json",
        "data/w33_pass11067_explicit_h27_central_cocycle_weld.json",
        "data/w33_pass11068_ternary_cocycle_deformation.json"
      ],
      "checks":{"blind_factors_6_to_9":True,"all_four_perfect_factors_blind":True,
                "twisted_right_actions_identical":True,"rank54_all_replays":True}
    }
    assert all(out["checks"].values())
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    a=ap.parse_args();p=payload();txt=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not OUT.exists() or OUT.read_text()!=txt:raise SystemExit("certificate drift")
    else:OUT.write_text(txt)
    print(json.dumps({"status":p["status"],"blind":[6,7,8,9]}))

if __name__=="__main__":
    raise SystemExit(main())
