#!/usr/bin/env python3
"""Pass 11068: scalar multiples of the cocycle give split s=0 and two oriented non-split s=+/-1 laws."""
from __future__ import annotations
import argparse,itertools,json,sys
from collections import Counter,deque
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_pass11067_explicit_h27_central_cocycle_weld as COC

OUT=ROOT/"data/w33_pass11068_ternary_cocycle_deformation.json"
ELTS=list(itertools.product(range(3),repeat=4));ID=(0,0,0,0)

def mul(x,y,s):
    h=x[:3];g=y[:3]
    z=COC.hmul(h,g)
    return (z[0],z[1],z[2],(x[3]+y[3]+s*COC.kappa(h,g))%3)

def inv(x,s):
    for y in ELTS:
        if mul(x,y,s)==ID and mul(y,x,s)==ID:return y
    raise RuntimeError
def comm(x,y,s):return mul(mul(mul(x,y,s),inv(x,s),s),inv(y,s),s)

def closure(gens,s):
    S={ID};Q=deque([ID])
    while Q:
        a=Q.popleft()
        for g in gens:
            z=mul(a,g,s)
            if z not in S:S.add(z);Q.append(z)
    return S

def order(x,s):
    z=ID
    for n in range(1,30):
        z=mul(z,x,s)
        if z==ID:return n
    raise RuntimeError

def invariants(s):
    center={x for x in ELTS if all(mul(x,y,s)==mul(y,x,s) for y in ELTS)}
    der=closure({comm(x,y,s) for x in ELTS for y in ELTS},s)
    g3=closure({comm(x,y,s) for x in ELTS for y in der},s)
    g4=closure({comm(x,y,s) for x in ELTS for y in g3},s)
    return {
      "center_order":len(center),"derived_order":len(der),
      "lower_central_sizes":[81,len(der),len(g3),len(g4)],
      "order_census":{str(k):v for k,v in sorted(Counter(order(x,s) for x in ELTS).items())}
    }

def payload():
    rows={str(s):invariants(s) for s in range(3)}
    assert rows["0"]=={"center_order":9,"derived_order":3,"lower_central_sizes":[81,3,1,1],"order_census":{"1":1,"3":80}}
    for s in (1,2):
        assert rows[str(s)]=={"center_order":3,"derived_order":9,"lower_central_sizes":[81,9,3,1],"order_census":{"1":1,"3":44,"9":36}}

    phi=lambda x:(x[0],x[1],x[2],2*x[3]%3)
    for x,y in itertools.product(ELTS,repeat=2):
        assert phi(mul(x,y,1))==mul(phi(x),phi(y),2)

    return {
      "schema":"w33.pass11068.ternary-cocycle-deformation.v1",
      "status":"PASS_ONE_TERNARY_COCYCLE_PARAMETER_INTERPOLATES_SPLIT_K81_AND_THE_TWO_ORIENTED_NONSPLIT_U81_LAWS",
      "headline":"Scale the nontrivial cocycle by s in F3: (h,d) star_s (h',D)=(hh',d+D+s kappa(h,h')). At s=0 this is exactly the split K81 law. At s=1 and s=2=-1 it is the class-three U81 law with opposite central orientation. The two nonzero twists are isomorphic by central reversal d->-d, while neither is isomorphic to the split law.",
      "family_formula":"(h,d) star_s (h',D)=(hh',d+D+s kappa(h,h')), s in F3",
      "invariants":rows,
      "orientation_isomorphism":{
        "map":"phi(h,d)=(h,-d)",
        "source_s":1,"target_s":2,
        "verified_pairs":6561
      },
      "phase_structure":{
        "s0":"split scheduler K81, class 2, exponent 3",
        "s1":"non-split chamber law, class 3, exponent 9",
        "s2":"central-orientation reverse of s1"
      },
      "boundary":"The parameter s is a discrete F3 cocycle coefficient. Treating s=+/-1 as physical chirality is an interpretation requiring a dynamical selection mechanism.",
      "parents":["data/w33_pass11067_explicit_h27_central_cocycle_weld.json","data/w33_pass11065_split_nonsplit_order81_bridge.json"],
      "checks":{"s0_is_split_K81":True,"s1_s2_have_U81_invariants":True,"central_reversal_maps_s1_to_s2":True,"split_nonsplit_class_jump":True}
    }

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args()
    p=payload();text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or a.output.read_text()!=text:raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"s":[0,1,2]},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
