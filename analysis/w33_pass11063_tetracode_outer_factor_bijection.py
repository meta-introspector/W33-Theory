#!/usr/bin/env python3
"""Pass 11063: the eight noncentral weld factors are indexed by the eight nonzero tetracode words."""
from __future__ import annotations
import argparse,json,sys,itertools
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_pass11061_hesse_slope_factor_alphabet as ALPH

OUT=ROOT/"data/w33_pass11063_tetracode_outer_factor_bijection.json"
POINTS=[(1,0),(0,1),(1,1),(1,2)]

def norm(v):
    a,b=v[0]%3,v[1]%3
    if a:z=pow(a,-1,3)
    elif b:z=pow(b,-1,3)
    else:raise ValueError
    return (a*z%3,b*z%3)

def covector(d):
    x,y=d
    return norm((y,-x))

def evalword(c):
    a,b=c
    return tuple((a*x+b*y)%3 for x,y in POINTS)

def payload():
    a=ALPH.payload()
    outer=[r for r in a["factor_rows"] if not r["central"]]
    code=sorted({tuple((aa*x+bb*y)%3 for x,y in POINTS) for aa,bb in itertools.product(range(3),repeat=2)})
    nonzero=[w for w in code if any(w)]
    assert len(code)==9 and len(nonzero)==8
    assert all(sum(z==0 for z in w)==1 for w in nonzero)

    maps=[]
    used=set()
    for i,d in enumerate(POINTS):
        cp=covector(d);wp=evalword(cp);wm=tuple((-z)%3 for z in wp)
        assert wp in nonzero and wm in nonzero
        rr=[r for r in outer if tuple(r["projective_hesse_direction"])==d]
        by={r["external_branch"]:r for r in rr};assert set(by)=={-1,1}
        maps.append({"clock_index":i,"projective_direction":list(d),"canonical_covector":list(cp),
                     "plus_factor":by[1]["factor"],"plus_word":list(wp),
                     "minus_factor":by[-1]["factor"],"minus_word":list(wm)})
        used|={wp,wm}
    assert used==set(nonzero)

    central=[r["factor"] for r in a["factor_rows"] if r["central"]]
    assert central==[8,9]
    return {
      "schema":"w33.pass11063.tetracode-outer-factor-bijection.v1",
      "status":"PASS_EIGHT_NONCENTRAL_WELD_FACTORS_BIJECT_WITH_THE_EIGHT_NONZERO_TETRACODE_WORDS_AFTER_BELL_LINE_COORDINATE_CHOICE",
      "headline":"Choose the four Hesse directions in the same projective order used to evaluate the ternary tetracode: (1,0),(0,1),(1,1),(1,-1). A nonzero linear covector vanishes at exactly one of these four points, and its two signs give the two nonzero tetracode words with that omitted coordinate. Mapping the two external weld branches over that Hesse direction to those two signs gives an exact 8-to-8 bijection.",
      "tetracode":{
        "formula":"C={(a,b,a+b,a-b):a,b in F3}",
        "size":9,"nonzero_words":8,"weight_of_every_nonzero_word":3
      },
      "factor_word_bijection":maps,
      "central_factors_outside_tetracode8":central,
      "weld_alphabet":"8 tetracode-indexed noncentral factors + 2 central FI factors",
      "boundary":"The bijection is exact after the displayed coordinate/sign convention. It is an indexing bridge; it does not identify a triangle-factor matrix with an E8 tensor sector or prove physical chirality.",
      "parents":["data/w33_pass11061_hesse_slope_factor_alphabet.json","data/w33_tetracode_e8_root_system_bridge.json"],
      "checks":{"tetracode_9_words":True,"eight_nonzero_weight3":True,"outer_factor_bijection":True,"central_pair_separate":True}
    }

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args()
    p=payload();text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or a.output.read_text()!=text:raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"bijection":8},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
