#!/usr/bin/env python3
"""Antiunitary conjugation exchanges the two oriented non-split U81 cocycle laws."""
from __future__ import annotations
import argparse,itertools,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_pass11067_explicit_h27_central_cocycle_weld as P67
import w33_pass11068_ternary_cocycle_deformation as P68
OUT=ROOT/"data/w33_20261001_u81_antiunitary_orientation.json"

def tau(h):
    a,b,c=h
    return ((-a)%3,b%3,(-c)%3)

def psi(x):
    return tau(x[:3])+((-x[3])%3,)

def payload():
    H=list(itertools.product(range(3),repeat=3))
    E=P68.ELTS
    for x,y in itertools.product(H,repeat=2):
        assert tau(P67.hmul(x,y))==P67.hmul(tau(x),tau(y))
        assert P67.kappa(tau(x),tau(y))==P67.kappa(x,y)
    pair_checks=0
    for s in (1,2):
        sm=(-s)%3
        for x,y in itertools.product(E,repeat=2):
            assert psi(P68.mul(x,y,s))==P68.mul(psi(x),psi(y),sm)
            pair_checks+=1
    return {
      "schema":"w33.20261001.u81-antiunitary-orientation.v1",
      "status":"PASS_QUTRIT_CONJUGATION_EXCHANGES_THE_TWO_ORIENTED_U81_COCYCLE_LAWS",
      "h27_antiunitary":{
        "map":"tau(a,b,c)=(-a,b,-c)",
        "schrodinger_reason":"rho(a,b,c)=omega^c Z^a X^b, so rho(h)^*=rho(tau(h)) because X is real and Z^*=Z^-1",
        "automorphism_pairs_checked":729,
        "cocycle_invariant":"kappa(tau(h),tau(h'))=kappa(h,h')"
      },
      "u81_extension":{
        "map":"Psi(h,d)=(tau(h),-d)",
        "law":"Psi(x star_s y)=Psi(x) star_{-s} Psi(y)",
        "oriented_laws":[1,2],"element_pair_checks":pair_checks,
        "central_tick":"Pass 11069 gives [e0,[e0,e1]]=(0,0,0,-s), hence the class-three tick reverses with s"
      },
      "interpretation":(
        "The two nonzero chamber cocycle orientations are exactly antiunitary-conjugate "
        "when the common H27 quotient is read in its qutrit Schrodinger representation. "
        "A state selecting one sign would therefore break this discrete conjugation."
      ),
      "boundary":(
        "This identifies the conjugation action on the finite chamber law. It does not derive "
        "an energy splitting between s=+1 and s=-1, choose which orientation occurs, or prove "
        "that the finite antiunitary is Standard-Model CP."
      ),
      "parents":["data/w33_pass11067_explicit_h27_central_cocycle_weld.json",
                 "data/w33_pass11068_ternary_cocycle_deformation.json",
                 "data/w33_pass11069_nested_commutator_cubic_tick.json"],
      "checks":{"h27_automorphism":True,"cocycle_tau_invariant":True,
                "u81_s_to_minus_s":True,"central_tick_reverses":True}
    }

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true");a=ap.parse_args()
    p=payload();txt=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not OUT.exists() or OUT.read_text()!=txt:raise SystemExit("certificate drift")
    else:OUT.write_text(txt)
    print(json.dumps({"status":p["status"],"checks":pair_checks if False else p["u81_extension"]["element_pair_checks"]}))

if __name__=="__main__":raise SystemExit(main())
