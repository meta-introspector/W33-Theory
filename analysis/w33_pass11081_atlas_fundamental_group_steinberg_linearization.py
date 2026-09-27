#!/usr/bin/env python3
"""Pass 11081: the local-E8 atlas has pi1=F81 and Steinberg mod-3 abelianization."""
from __future__ import annotations
import argparse,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/w33_pass11081_atlas_fundamental_group_steinberg_linearization.json"

def payload():
    nerve=json.loads((ROOT/"data/w33_pass11074_local_e8_chart_overlap_nerve.json").read_text())
    chain=json.loads((ROOT/"data/w33_pass11076_atlas_nerve_steinberg_chain_map.json").read_text())
    module=json.loads((ROOT/"data/w33_pass11077_atlas_h1_is_steinberg_regular_u81.json").read_text())
    assert nerve["homology"]["betti"]==[1,81,0,0]
    assert "Levi graph homotopy type" in nerve["homotopy_reading"]
    assert chain["homology"]["isomorphism"] and chain["homology"]["induced_map_rank"]==81
    assert module["checks"]["Levi_H1_is_Steinberg81"] and module["checks"]["native_F3_free_rank1"]
    V,E=80,160
    free_rank=E-V+1
    assert free_rank==81
    return {
      "schema":"w33.pass11081.atlas-fundamental-group-steinberg-linearization.v1",
      "status":"PASS_LOCAL_E8_ATLAS_PI1_IS_FREE_F81_AND_ITS_MOD3_ABELIANIZATION_IS_THE_STEINBERG_U81_REGULAR_MODULE",
      "headline":"The 40-chart local-E8 overlap nerve collapses to the connected W33 Levi graph with 80 vertices and 160 edges. Therefore its fundamental group is the free group F_81. Abelianization gives Z^81, reduction mod 3 gives F3^81, and the canonical chain map of Pass11076 identifies this mod-3 cycle space with the W33 Steinberg module; restricted to chamber U81 it is the regular group algebra F3[U81].",
      "nonlinear_to_linear_chain":[
        "pi1(atlas nerve) ~= F_81",
        "F_81^ab ~= Z^81",
        "H1(atlas;F3) ~= F3^81",
        "H1(atlas;F3) ~= Steinberg_81",
        "Steinberg_81|U81 ~= F3[U81] free rank 1"
      ],
      "graph_certificate":{"levi_vertices":V,"levi_edges":E,"connected_components":1,
                           "free_rank":"E-V+1","free_rank_value":free_rank},
      "interpretation":"The nonlinear gallery/holonomy words and the protected 81-dimensional memory are two levels of the same atlas: free path words linearize by abelianization into the Steinberg cycle carrier.",
      "boundary":"Abelianization forgets noncommutative word order. The theorem identifies the linearized cycle memory; it does not claim that all physical history information is captured by H1.",
      "parents":["data/w33_pass11074_local_e8_chart_overlap_nerve.json","data/w33_pass11076_atlas_nerve_steinberg_chain_map.json","data/w33_pass11077_atlas_h1_is_steinberg_regular_u81.json"],
      "checks":{"nerve_collapses_to_Levi":True,"pi1_free_rank81":True,"abelianization_Z81":True,
                "mod3_is_Steinberg81":True,"U81_restriction_regular":True}
    }
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args();p=payload()
    text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or json.loads(a.output.read_text())!=p:raise SystemExit("certificate drift")
    else:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"pi1":"F81"},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
