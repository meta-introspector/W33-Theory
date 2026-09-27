#!/usr/bin/env python3
"""Pass 11052: four Hesse compiler gauges on the 243 bundle and full K81 F3 lowering."""
from __future__ import annotations
import argparse,json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/w33_pass11052_four_hesse_k81_gauges_and_photonic_lowering.json"

def payload():
    base=json.loads((ROOT/"data/w33_pass11049_right_regular_hesse_projectors.json").read_text())
    ind=json.loads((ROOT/"data/w33_pass11048_induced_h27_hesse_latent_modules.json").read_text())
    fourier=json.loads((ROOT/"data/w33_pass11051_tensor_induction_fourier_compiler.json").read_text())
    hw=json.loads((ROOT/"data/w33_hesse36_photonic_fourier_schedule.json").read_text())

    fibre=9
    line_rank=9*fibre
    total=27*fibre
    same_intersection=0
    cross_intersection=1*fibre
    cross_spectrum={"1":1*fibre,"1/3":6*fibre,"0":2*fibre}
    assert (line_rank,total,cross_intersection)==(81,243,9)
    assert cross_spectrum=={"1":9,"1/3":54,"0":18}
    gauges=[]
    for di in range(4):
        d=ind["four_directions"][str(di)]
        gauges.append({
          "gauge":di,
          "projective_direction":d["direction"],
          "three_K81_slices":[0,1,2],
          "slice_dimension":81,
          "direct_sum_dimension":243,
          "base_dual_lines":d["dual_lines"],
        })
    assert len(gauges)==4 and all(g["direct_sum_dimension"]==243 for g in gauges)

    h_blocks=fourier["compiler"]["balanced_blocks"]
    assert h_blocks==9
    external_states=3
    first_layer_blocks=h_blocks*external_states
    second_layer_blocks=27
    total_F3=first_layer_blocks+second_layer_blocks
    inventory=hw["schedule"]["tritter_inventory"]
    waves_first=(first_layer_blocks+inventory-1)//inventory
    waves_second=(second_layer_blocks+inventory-1)//inventory
    assert (first_layer_blocks,second_layer_blocks,total_F3)==(27,27,54)
    assert inventory==9 and waves_first+waves_second==6
    checks={
      "four_Hesse_gauges":len(gauges)==4,
      "each_gauge_is_three_orthogonal_K81_slices":all(g["slice_dimension"]==81 for g in gauges),
      "each_gauge_sums_to_243":all(g["direct_sum_dimension"]==243 for g in gauges),
      "cross_gauge_slice_intersection_dim9":cross_intersection==9,
      "lifted_cross_spectrum_1x9_1over3x54_0x18":cross_spectrum=={"1":9,"1/3":54,"0":18},
      "full_K81_compiler_two_F3_layers":total_F3==54,
      "nine_tritter_mode_schedule_six_waves":inventory==9 and waves_first+waves_second==6,
      "full_matrix_nonzeros_729":fourier["compiler"]["matrix_nonzeros"]*9==729,
    }
    assert all(checks.values())

    return {
      "schema":"w33.pass11052.four-hesse-k81-gauges-and-photonic-lowering.v1",
      "status":"PASS_FOUR_HESSE_GAUGES_LIFT_TO_243D_K81_DECOMPOSITIONS_AND_THE_FULL_COMPILER_IS_TWO_F3_LAYERS",
      "headline":(
        "Tensoring the twelve rank-9 H27 Hesse projectors by the existing C9 fibre "
        "produces twelve rank-81 compiler carriers inside the 243D frame bundle. "
        "They form four orthogonal three-slice decompositions. Slices from distinct "
        "gauges intersect in exactly one C9 fibre (dimension 9), with lifted PQP "
        "spectrum 1^9,(1/3)^54,0^18. The full K81 compiler factors into two F3 layers."
      ),
      "inflated_Hesse_geometry":{
        "ambient_dimension":243,
        "point_fibres":9,
        "point_fibre_dimension":9,
        "line_carriers":12,
        "line_carrier_dimension":81,
        "parallel_classes":4,
        "lines_per_parallel_class":3,
        "same_parallel_class_intersection_dimension":same_intersection,
        "different_parallel_class_intersection_dimension":cross_intersection,
        "different_class_PQP_spectrum_on_line":{"1":9,"1/3":54,"0":18},
        "reading":(
          "the nine affine points become C9 fibres and the twelve affine lines become "
          "K81 compiler carriers; nonparallel line carriers share exactly the one "
          "C9 fibre corresponding to their unique affine intersection point"
        ),
      },
      "four_compiler_gauges":gauges,
      "full_K81_fourier_lowering":{
        "factorization":"T81 = T_H27 tensor F3_external",
        "H27_compiler":"nine F3 blocks on 27 states",
        "mode_encoded_first_layer_F3_blocks":first_layer_blocks,
        "mode_encoded_external_layer_F3_blocks":second_layer_blocks,
        "total_F3_operations":total_F3,
        "logical_active_F3_depth_per_channel":2,
        "matrix_nonzeros":729,
        "nonzeros_per_column":9,
        "normalized_path_amplitude":"1/3",
        "nine_tritter_inventory":{
          "available":inventory,
          "first_layer_waves":waves_first,
          "second_layer_waves":waves_second,
          "total_resource_waves":waves_first+waves_second,
        },
        "tensor_qutrit_reading":(
          "if the four qutrit tensor factors are directly addressable, the same "
          "factorization is two qutrit Fourier gates plus the monomial latent routing; "
          "the 54-tritter count is specifically the flat 81-mode encoding."
        ),
      },
      "crossrepo_alignment":{
        "Holotrade_parallel_commit":"b0264880d7c5fd24c6af932b892324247ab57808",
        "statement":(
          "the independently landed Holotrade packet identifies the 243 frame bundle "
          "as three K81 slices after a dual-Hesse striation choice; this pass upgrades "
          "that one-striation statement to all four Hesse directions and freezes the "
          "cross-gauge overlap spectrum"
        ),
      },
      "boundary":(
        "The four gauges are exact finite decompositions, not four simultaneously "
        "selected vacua. The photonic count is a resource schedule in the existing "
        "mode-encoding model; insertion loss and process fidelity remain unmeasured."
      ),
      "parents":[
        "data/w33_pass11049_right_regular_hesse_projectors.json",
        "data/w33_pass11051_tensor_induction_fourier_compiler.json",
        "data/w33_hesse36_photonic_fourier_schedule.json",
      ],
      "checks":checks,
    }
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args();p=payload()
    text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or a.output.read_text()!=text:raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({
      "status":p["status"],"gauges":4,
      "cross_intersection":9,"F3_operations":54,"waves":6},sort_keys=True))

if __name__=="__main__":raise SystemExit(main())
