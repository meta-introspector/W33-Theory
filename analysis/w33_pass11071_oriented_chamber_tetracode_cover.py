#!/usr/bin/env python3
"""Pass 11071: the 320 oriented chambers are the G/U81 double cover of flags and locally the eight tetracode words per line."""
from __future__ import annotations
import argparse,itertools,json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/w33_pass11071_oriented_chamber_tetracode_cover.json"
POINTS=[(1,0),(0,1),(1,1),(1,2)]

def word(c):
    a,b=c
    return tuple((a*x+b*y)%3 for x,y in POINTS)

def payload():
    d=json.loads((ROOT/"data/w33_pass11070_dual_parabolic_flag_diamond.json").read_text())
    old=json.loads((ROOT/"data/w33_pass11063_tetracode_outer_factor_bijection.json").read_text())
    G=d["ambient"]["order"];B=d["flag_borel"]["order"];U=d["flag_borel"]["sylow3_order"]
    assert (G//B,G//U,B//U)==(160,320,2)

    covs=[c for c in itertools.product(range(3),repeat=2) if c!=(0,0)]
    fibres=[];used=set()
    for i,p in enumerate(POINTS):
        cs=[c for c in covs if (c[0]*p[0]+c[1]*p[1])%3==0]
        assert len(cs)==2 and cs[1]==tuple((-x)%3 for x in cs[0])
        ws=[word(c) for c in cs]
        assert all(sum(x!=0 for x in w)==3 for w in ws)
        used.update(ws)
        fibres.append({"omitted_point_index":i,"projective_point":list(p),
                       "covectors":[list(c) for c in cs],"tetracode_words":[list(w) for w in ws]})
    assert len(used)==8
    oldwords=set()
    for row in old["factor_word_bijection"]:
        oldwords.add(tuple(row["plus_word"]));oldwords.add(tuple(row["minus_word"]))
    assert used==oldwords

    return {
      "schema":"w33.pass11071.oriented-chamber-tetracode-cover.v1",
      "status":"PASS_G_OVER_U81_IS_THE_320_SHEET_ORIENTED_CHAMBER_COVER_AND_EACH_LINE_FIBRE_IS_THE_EIGHT_NONZERO_TETRACODE_WORDS",
      "headline":(
        "Because the flag Borel B has order 162 and its unique Sylow-3 subgroup U81 has index two, "
        "the homogeneous-space map G/U81 -> G/B is a canonical two-sheeted cover of the 160 W33 flags, "
        "with 320 oriented chambers. On a fixed four-point line, a flag chooses the omitted point; exactly "
        "two nonzero covectors vanish there, and their evaluations are the opposite pair of weight-three "
        "tetracode words. Thus the eight oriented chambers over one line are exactly its eight nonzero tetracode words after coordinate choice."
      ),
      "homogeneous_cover":{
        "group":"PSp(4,3)","group_order":G,
        "unoriented_stabilizer":"B","unoriented_stabilizer_order":B,"flags":G//B,
        "oriented_stabilizer":"U81","oriented_stabilizer_order":U,"oriented_chambers":G//U,
        "sheets":B//U,"map":"G/U81 -> G/B"
      },
      "line_tetracode":{
        "projective_points":[list(p) for p in POINTS],
        "evaluation_formula":"(a,b) -> (a,b,a+b,a-b)",
        "nonzero_words":8,"weight":3,
        "oriented_chambers_per_line":8,
        "fibres":fibres
      },
      "global_count":{"lines":40,"oriented_per_line":8,"total":320},
      "coordinate_boundary":(
        "The identification with signed tetracode words uses a coordinate/sign convention on each line. "
        "The two-sheet homogeneous cover itself is canonical; a globally glued tetracode/E8 chart bundle requires additional transition data."
      ),
      "parents":["data/w33_pass11070_dual_parabolic_flag_diamond.json","data/w33_pass11063_tetracode_outer_factor_bijection.json"],
      "checks":{"flag_count160":True,"oriented_count320":True,"double_cover":True,
                "two_words_per_omitted_point":True,"eight_words_per_line":True,
                "matches_existing_tetracode_factor_words":True}
    }

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args();p=payload()
    text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or json.loads(a.output.read_text())!=p:raise SystemExit("certificate drift")
    else:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"flags":160,"oriented":320},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
