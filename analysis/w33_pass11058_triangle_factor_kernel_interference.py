#!/usr/bin/env python3
"""Pass 11058: triangle-factor kernels, noncommutation, and local synthesis structure."""
from __future__ import annotations
import argparse,itertools,json,sys
from collections import Counter
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_pass11056_diagonal_weld_ten_triangle_factors as FCT

OUT=ROOT/"data/w33_pass11058_triangle_factor_kernel_interference.json"

def rank_mod(A,p): return FCT.rank_mod(A,p)

def symp(s,t):
    a,b,_=s[0]; d,e,_=t[0]
    return (a*e-b*d)%3

def payload():
    coords,v,D,pairs,F=FCT.geometry()
    pair_hist=Counter(); zero_comm=0
    for i in range(10):
        for j in range(i+1,10):
            nd103=81-rank_mod(np.vstack([F[i],F[j]]),103)
            nd109=81-rank_mod(np.vstack([F[i],F[j]]),109)
            assert nd103==nd109
            sp=symp(pairs[i][0],pairs[j][0])
            expected=9 if sp==0 else 3
            assert nd103==expected
            pair_hist[("symplectic_zero" if sp==0 else "symplectic_nonzero",expected)]+=1
            C=F[i]@F[j]-F[j]@F[i]
            if rank_mod(C,103)==0: zero_comm+=1
    assert zero_comm==0
    assert pair_hist[("symplectic_zero",9)]==21
    assert pair_hist[("symplectic_nonzero",3)]==24

    stack=np.vstack(F)
    assert rank_mod(stack,103)==rank_mod(stack,109)==80
    vv=np.array(v,dtype=np.int64)
    assert all(np.count_nonzero(A@vv)==0 for A in F)
    assert rank_mod(D,103)==rank_mod(D,109)==66
    assert np.count_nonzero(D@vv)==0

    freq=Counter()
    for A in F:
        for cc in FCT.components(A):
            vals=[abs(int(A[a,b])) for a,b in itertools.combinations(cc,2) if A[a,b]!=0]
            assert len(vals)==3
            freq[sum(x*x for x in vals)]+=1
    assert sum(freq.values())==270

    return {
      "schema":"w33.pass11058.triangle-factor-kernel-interference.v1",
      "status":"PASS_TEN_LOCAL_TRIANGLE_FACTORS_SHARE_ONLY_THE_BACKGROUND_RAY_WHILE_THEIR_SUM_CREATES_14_EXTRA_DARK_MODES",
      "headline":"The ten rank-54 triangle factors are strongly noncommuting, but their kernel geometry is rigid. Every pair whose H27 connection directions commute has a 9D common kernel; every noncommuting pair has a 3D common kernel. All ten factors together share only one ray, exactly the diagonal background v_plus. Yet their sum D_plus has a 15D kernel. Thus fourteen dark modes are created by destructive interference among active factors rather than being inactive in every local layer.",
      "pair_kernel_law":{
        "symplectic_zero_pairs":21,"common_kernel_dimension":9,
        "symplectic_nonzero_pairs":24,"common_kernel_dimension_nonzero_case":3,
        "pairwise_commuting_factor_pairs":0
      },
      "global_kernel":{
        "common_kernel_all_ten_dimension":1,
        "common_kernel_generator":"v_plus",
        "rank_of_full_D_plus":66,"kernel_dimension_of_full_D_plus":15,
        "extra_dark_dimensions_created_by_interference":14
      },
      "local_synthesis":{
        "factors":10,"disjoint_triangle_blocks_per_factor":27,
        "triangle_block_exponentials_per_factor":27,
        "parallel_width_within_one_factor":27,
        "generator_layer_depth_for_one_factor_sweep":10,
        "total_triangle_blocks_per_sweep":270,
        "exact_product_equals_exp_of_sum":False,
        "reason":"no two weighted factors commute; product formulas require noncommuting synthesis/Trotter control"
      },
      "triangle_rotation_frequency_squared_histogram":{str(k):int(freq[k]) for k in sorted(freq)},
      "boundary":"This supplies a local three-mode decomposition of the generator and exact kernel/interference structure. It does not claim that a ten-layer product of factor exponentials equals the desired global exponential or Cayley gate.",
      "parents":["data/w33_pass11056_diagonal_weld_ten_triangle_factors.json","data/w33_pass11055_cayley_finite_54_compiler.json"],
      "checks":{
        "pair_kernel_law_verified_at_103_and_109":True,
        "no_two_factors_commute":True,
        "all_ten_common_kernel_is_exactly_background_ray":True,
        "full_sum_has_14_additional_interference_dark_modes":True,
        "each_factor_exponential_splits_into_27_parallel_three_mode_blocks":True
      }
    }

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args()
    p=payload();text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or a.output.read_text()!=text:raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"common_kernel":1,"sum_kernel":15,"interference_dark":14},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
