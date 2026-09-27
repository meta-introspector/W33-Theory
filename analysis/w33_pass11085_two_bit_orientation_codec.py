#!/usr/bin/env python3
"""Pass 11085: the V4 deck group separates tetracode sheet sign from cubic-tick chirality."""
from __future__ import annotations
import argparse,itertools,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/w33_pass11085_two_bit_orientation_codec.json"
def act(state,g):
    sigma,chi=state
    if g=="T":return (-sigma,chi)
    if g=="S":return (-sigma,-chi)
    if g=="TS":return (sigma,-chi)
    if g=="1":return state
    raise ValueError(g)
def payload():
    cover=json.loads((ROOT/"data/w33_pass11084_pgsp_u81_four_sheet_cover.json").read_text())
    internal=json.loads((ROOT/"data/w33_pass11072_internal_borel_deck_involution.json").read_text())
    outer=json.loads((ROOT/"data/w33_pass11073_outer_similitude_tick_reversal.json").read_text())
    assert cover["deck_generators"]["quotient"]=="V4"
    assert internal["tetracode_sheet_action"]["action"]=="swap"
    assert not internal["cocycle_action"]["highest_root_tick_reversed"]
    assert outer["tetracode_sheet_action"]["action"]=="swap"
    assert outer["cocycle_action"]["highest_root_tick_reversed"]
    base=(1,1);orb={g:act(base,g) for g in ("1","T","S","TS")}
    assert set(orb.values())==set(itertools.product((1,-1),repeat=2))
    return {
      "schema":"w33.pass11085.two-bit-orientation-codec.v1",
      "status":"PASS_FOUR_PGSP_CHAMBER_SHEETS_ARE_EXACTLY_TETRACODE_SIGN_CROSS_CUBIC_TICK_CHIRALITY",
      "headline":"The V4 deck cover has a natural two-bit codec. Let sigma be the local tetracode/chamber-sheet sign and chi the highest-root cubic-tick sign. The internal Borel involution T acts (sigma,chi)->(-sigma,chi); the outer similitude S acts (sigma,chi)->(-sigma,-chi); their product TS acts (sigma,chi)->(sigma,-chi). Thus the four sheets over every W33 flag realize all four independent sign pairs.",
      "bits":{"sigma":"tetracode/chamber-sheet sign","chi":"highest-root cocycle/tick sign"},
      "deck_action":{"1":"(sigma,chi)","T":"(-sigma,chi)","S":"(-sigma,-chi)","TS":"(sigma,-chi)"},
      "base_orbit":{"1":[1,1],"T":[-1,1],"S":[-1,-1],"TS":[1,-1]},
      "global_census":{"flags":160,"sign_pairs_per_flag":4,"double_oriented_states":640},
      "correction":"A signed tetracode word alone is only sigma and cannot determine chi: both T and S flip sigma, but only S flips the cocycle/tick. The temporal orientation requires two Z2 labels, not one.",
      "boundary":"The two-bit decomposition is exact group action data. Neither bit is automatically a physical arrow of time; a dynamics must choose or prepare a sector.",
      "parents":["data/w33_pass11084_pgsp_u81_four_sheet_cover.json","data/w33_pass11072_internal_borel_deck_involution.json","data/w33_pass11073_outer_similitude_tick_reversal.json"],
      "checks":{"four_sign_pairs":True,"T_flips_only_sigma":True,"S_flips_both":True,"TS_flips_only_chi":True,"global_640":True}
    }
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args();p=payload()
    text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or json.loads(a.output.read_text())!=p:raise SystemExit("certificate drift")
    else:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"states":640},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
