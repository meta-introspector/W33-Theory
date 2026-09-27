#!/usr/bin/env python3
"""Pass 11072: the internal chamber-sheet involution swaps tetracode sign but fixes the class-three central tick."""
from __future__ import annotations
import argparse,itertools,json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/w33_pass11072_internal_borel_deck_involution.json"

def E(i,j):
    M=np.zeros((4,4),dtype=np.int64);M[i,j]=1;return M

I=np.eye(4,dtype=np.int64)%3
J=np.array([[0,0,1,0],[0,0,0,1],[-1,0,0,0],[0,-1,0,0]],dtype=np.int64)%3
X=[(E(0,1)-E(3,2))%3,E(1,3)%3,(E(0,3)+E(1,2))%3,E(0,2)%3]
T=np.diag([1,2,1,2]).astype(np.int64)%3

def hmul(x,y):
    a,b,c=x;A,B,C=y
    return ((a+A)%3,(b+B)%3,(c+C-A*b)%3)

def kappa(x,y):
    a,b,c=x;A,B,C=y
    return (A*A*b-2*A*c)%3

def psi(h):
    a,b,c=h;return ((-a)%3,b,(-c)%3)

def payload():
    assert np.array_equal((T.T@J@T)%3,J)
    assert np.array_equal((T@T)%3,I)
    Tin=T
    signs=[]
    for A in X:
        C=(T@A@Tin)%3
        s=next(s for s in (1,2) if np.array_equal(C,(s*A)%3))
        signs.append(s)
    assert signs==[2,1,2,1]

    H=list(itertools.product(range(3),repeat=3))
    assert all(hmul(psi(x),psi(y))==psi(hmul(x,y)) for x in H for y in H)
    assert all(kappa(psi(x),psi(y))==kappa(x,y) for x in H for y in H)

    # Standard chamber: p=<e0> inside L=<e0,e1>.  T|L=diag(1,-1).
    ell=np.array([0,1],dtype=np.int64)
    TL=np.diag([1,2]).astype(np.int64)%3
    ell2=(ell@TL)%3
    assert np.array_equal(ell2,np.array([0,2]))

    return {
      "schema":"w33.pass11072.internal-borel-deck-involution.v1",
      "status":"PASS_INTERNAL_BOREL_C2_SWAPS_THE_TETRACODE_SHEET_BUT_FIXES_THE_HIGHEST_ROOT_TICK",
      "headline":(
        "The nontrivial C2 in the projective flag Borel can be represented by the symplectic "
        "diagonal matrix T=diag(1,-1,1,-1). It fixes the standard point-line flag and swaps "
        "the two nonzero covectors vanishing at the omitted point, hence swaps the local +/- "
        "tetracode words. But on C2 root coordinates it acts (a,b,c,d)->(-a,b,-c,d): "
        "the highest-root central coordinate d is fixed, and the H27 central cocycle kappa is invariant."
      ),
      "matrix":{"T":"diag(1,-1,1,-1) over F3","symplectic_multiplier":1,"T_squared":"I",
                "fixed_flag":"<e0> inside <e0,e1>"},
      "root_coordinate_action":{"a_short":"-a","b_long":"b","c_alpha_plus_beta":"-c","d_highest_root":"d",
                                "sign_vector":["-","+","-","+"]},
      "tetracode_sheet_action":{"omitted_point":"<e0>","vanishing_covectors":["(0,1)","(0,-1)"],
                                "action":"swap","deck_character":"nontrivial C2"},
      "cocycle_action":{"H27_automorphism":"psi(a,b,c)=(-a,b,-c)",
                        "identity":"kappa(psi(x),psi(y))=kappa(x,y)","checked_pairs":729,
                        "highest_root_tick_reversed":False},
      "firewall":(
        "The canonical two-sheet chamber/tetracode orientation is therefore not, by itself, "
        "the sign of the class-three central tick. The internal Borel deck map changes the former and fixes the latter."
      ),
      "parents":["data/w33_pass11071_oriented_chamber_tetracode_cover.json","data/w33_pass11067_explicit_h27_central_cocycle_weld.json"],
      "checks":{"T_is_symplectic":True,"T_is_involution":True,"root_signs_exact":True,
                "tetracode_sign_swapped":True,"kappa_invariant":True,"highest_root_fixed":True}
    }

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args();p=payload()
    text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or json.loads(a.output.read_text())!=p:raise SystemExit("certificate drift")
    else:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"root_signs":["-","+","-","+"]},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
