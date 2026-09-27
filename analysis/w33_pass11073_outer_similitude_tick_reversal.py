#!/usr/bin/env python3
"""Pass 11073: an outer symplectic similitude reverses the cubic cocycle/tick, unlike the internal Borel deck involution."""
from __future__ import annotations
import argparse,itertools,json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/w33_pass11073_outer_similitude_tick_reversal.json"

def E(i,j):
    M=np.zeros((4,4),dtype=np.int64);M[i,j]=1;return M

I=np.eye(4,dtype=np.int64)%3
J=np.array([[0,0,1,0],[0,0,0,1],[-1,0,0,0],[0,-1,0,0]],dtype=np.int64)%3
X=[(E(0,1)-E(3,2))%3,E(1,3)%3,(E(0,3)+E(1,2))%3,E(0,2)%3]
S=np.diag([1,2,2,1]).astype(np.int64)%3

def hmul(x,y):
    a,b,c=x;A,B,C=y
    return ((a+A)%3,(b+B)%3,(c+C-A*b)%3)

def kappa(x,y):
    a,b,c=x;A,B,C=y
    return (A*A*b-2*A*c)%3

def psi(h):
    a,b,c=h;return ((-a)%3,(-b)%3,c)

def payload():
    assert np.array_equal((S.T@J@S)%3,(2*J)%3)
    assert np.array_equal((S@S)%3,I)
    signs=[]
    for A in X:
        C=(S@A@S)%3
        s=next(s for s in (1,2) if np.array_equal(C,(s*A)%3))
        signs.append(s)
    assert signs==[2,2,1,2]
    H=list(itertools.product(range(3),repeat=3))
    assert all(hmul(psi(x),psi(y))==psi(hmul(x,y)) for x in H for y in H)
    assert all(kappa(psi(x),psi(y))==(-kappa(x,y))%3 for x in H for y in H)

    ell=np.array([0,1],dtype=np.int64)
    SL=np.diag([1,2]).astype(np.int64)%3
    assert np.array_equal((ell@SL)%3,np.array([0,2]))

    internal=json.loads((ROOT/"data/w33_pass11072_internal_borel_deck_involution.json").read_text())
    assert internal["checks"]["tetracode_sign_swapped"] and internal["checks"]["highest_root_fixed"]

    return {
      "schema":"w33.pass11073.outer-similitude-tick-reversal.v1",
      "status":"PASS_OUTER_PGSP_INVOLUTION_REVERSES_THE_COCYCLE_AND_HIGHEST_ROOT_TICK_WHILE_THE_INTERNAL_BOREL_DECK_DOES_NOT",
      "headline":(
        "The diagonal similitude S=diag(1,-1,-1,1) satisfies S^T J S=-J, so its projective "
        "class preserves W33 isotropic incidence but lies in the PGSp outer coset rather than PSp. "
        "It fixes the same standard flag and also swaps the local tetracode +/- pair. On root coordinates "
        "it acts (a,b,c,d)->(-a,-b,c,-d); the induced H27 automorphism sends kappa to -kappa and the "
        "highest-root central tick changes sign."
      ),
      "matrix":{"S":"diag(1,-1,-1,1) over F3","similitude_multiplier":"-1","S_squared":"I",
                "group_location":"PGSp(4,3) minus PSp(4,3)","fixed_flag":"<e0> inside <e0,e1>"},
      "root_coordinate_action":{"a_short":"-a","b_long":"-b","c_alpha_plus_beta":"c","d_highest_root":"-d",
                                "sign_vector":["-","-","+","-"]},
      "cocycle_action":{"H27_automorphism":"psi(a,b,c)=(-a,-b,c)",
                        "identity":"kappa(psi(x),psi(y))=-kappa(x,y)","checked_pairs":729,
                        "highest_root_tick_reversed":True},
      "tetracode_sheet_action":{"action":"swap","same_local_pair_as_internal_deck":True},
      "orientation_firewall":{
        "internal_Borel_C2":{"tetracode_sign":"swapped","kappa":"fixed","highest_root_tick":"fixed"},
        "outer_PGSp_C2":{"tetracode_sign":"swapped","kappa":"negated","highest_root_tick":"negated"},
        "consequence":"tetracode/chamber sheet sign and cubic-tick chirality are distinct Z2 data"
      },
      "boundary":(
        "This is an exact finite symplectic/similitude distinction. Calling the outer involution time reversal "
        "or physical chirality requires a dynamical interpretation not supplied by the group calculation."
      ),
      "parents":["data/w33_pass11072_internal_borel_deck_involution.json","data/w33_pass11067_explicit_h27_central_cocycle_weld.json"],
      "checks":{"S_is_minus_symplectic_similitude":True,"S_is_involution":True,"root_signs_exact":True,
                "kappa_anti_invariant":True,"highest_root_reversed":True,"internal_outer_Z2_distinguished":True}
    }

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args();p=payload()
    text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or json.loads(a.output.read_text())!=p:raise SystemExit("certificate drift")
    else:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"firewall":"sheet_Z2 != tick_Z2"},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
