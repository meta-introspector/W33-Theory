#!/usr/bin/env python3
"""Pass 11077: the local-E8 atlas overlap H1 is the W33 Steinberg module and regular on U81."""
from __future__ import annotations
import argparse,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/w33_pass11077_atlas_h1_is_steinberg_regular_u81.json"
def payload():
    a=json.loads((ROOT/"data/w33_pass11076_atlas_nerve_steinberg_chain_map.json").read_text())
    st=json.loads((ROOT/"data/PART_W33_20260901_DOUBLE_BUILDING_STEINBERG_HOMOLOGY.json").read_text())
    u=json.loads((ROOT/"data/PART_W33_PASS5105_U81_DUAL_TORSOR_CONTROLLER.json").read_text())
    ch=json.loads((ROOT/"data/w33_pass11071_oriented_chamber_tetracode_cover.json").read_text())
    assert a["homology"]["isomorphism"] and a["homology"]["induced_map_rank"]==81
    assert st["buildings"]["W33_GQ33"]["cycleRank"]==81
    assert st["characterCertificate"]["H1_W33_degree"]==81 and st["characterCertificate"]["H1_W33_norm"]==1
    assert "Steinberg" in st["solomonTitsInterpretation"]
    assert u["U"]["order"]==81
    assert u["Steinberg"]["complex_restriction"]=="Reg(U81)"
    assert "free rank 1" in u["Steinberg"]["native_restriction"]
    assert ch["homogeneous_cover"]["oriented_stabilizer_order"]==81
    return {
      "schema":"w33.pass11077.atlas-h1-is-steinberg-regular-u81.v1",
      "status":"PASS_LOCAL_E8_ATLAS_OVERLAP_CYCLES_ARE_THE_STEINBERG81_AND_RESTRICT_REGULARLY_TO_THE_CHAMBER_U81",
      "headline":"Pass11076 canonically identifies the 81-dimensional H1 of the 40-chart overlap nerve with W33 Levi H1. The older Solomon-Tits certificate identifies that Levi H1 as the irreducible degree-81 Steinberg module. The certified chamber Sylow-3 subgroup U81 has order 81, and Steinberg restricted to U81 is exactly its regular module, over both C and the native field F3.",
      "module_chain":[
        "H1(local-E8 chart nerve) ~= H1(W33 Levi building)",
        "H1(W33 Levi building) = Steinberg_81",
        "Steinberg_81|U81 = Reg(U81) over C",
        "H1(F3)|U81 ~= F3[U81] free rank 1"
      ],
      "consequence":{
        "atlas_cycle_dimension":81,"chamber_update_group_order":81,
        "native_U81_module":"free rank one group algebra F3[U81]",
        "complex_U81_module":"regular representation",
        "dimension_match_is_not_accidental":True,
        "reading":"The 81 independent overlap/gluing cycles of the local-E8 atlas are the same protected building-homology carrier already used by the chamber update controller."
      },
      "boundary":"This is an exact module identification. It does not prove that a nontrivial Lie-algebra cocycle is present on all 81 cycles; explicit A2->E8 transition maps remain the missing global-amalgam datum.",
      "parents":["data/w33_pass11076_atlas_nerve_steinberg_chain_map.json","data/PART_W33_20260901_DOUBLE_BUILDING_STEINBERG_HOMOLOGY.json","data/PART_W33_PASS5105_U81_DUAL_TORSOR_CONTROLLER.json"],
      "checks":{"atlas_to_Levi_H1_iso":True,"Levi_H1_is_Steinberg81":True,"U81_order81":True,"complex_regular_restriction":True,"native_F3_free_rank1":True}
    }
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args();p=payload()
    text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or json.loads(a.output.read_text())!=p:raise SystemExit("certificate drift")
    else:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"module":"St81","U":"regular"},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
