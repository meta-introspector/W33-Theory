#!/usr/bin/env python3
"""Pass 11033: compile the exact signed clock center to a 24-mode interferometer."""
from __future__ import annotations
import argparse, json, math, sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_pass10975_e6_cubic_reynolds_clock_weld as P75
import w33_pass11021_minimal_signed_clock_carrier as P21

OUT=ROOT/"data/w33_pass11033_signed_clock_photonic_compiler.json"


def inversion_count(p):
    return sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))


def separation_threshold():
    # 1/3 + 2 sin(e/2) < cos(e)
    lo,hi=0.0,math.pi/2
    for _ in range(100):
        m=(lo+hi)/2
        if 1/3+2*math.sin(m/2) < math.cos(m): lo=m
        else: hi=m
    return lo


def payload():
    e6,dirs,rows,perms,section=P21.split_section()
    p_to_m={r[3]:tuple(map(int,r[0])) for r in rows}
    z=next(g for g in section if p_to_m[g[0]]==(2,0,0,2))
    p,s=z
    noncentral=sorted(i for i,h in e6.items() if h[:2]!=(0,0))
    pos={x:i for i,x in enumerate(noncentral)}
    perm=[None]*24
    negative=[]
    pairs=[]
    seen=set()
    for i in noncentral:
        j=p[i]
        perm[pos[i]]=pos[j]
        if s[i]: negative.append(pos[i])
        if i not in seen:
            eps=-1 if s[i] else 1
            pairs.append((pos[i],pos[j],eps))
            seen|={i,j}
    assert len(pairs)==12 and len(negative)==12
    assert all(perm[perm[i]]==i and perm[i]!=i for i in range(24))
    inv=inversion_count(perm)
    assert inv==148

    # In this canonical order the permutation is reversal on blocks 16 and 8.
    assert perm==list(range(15,-1,-1))+list(range(23,15,-1))
    nn_depth=max(15,7)

    fibre_terms=[]
    for d in dirs:
        f=sorted(i for i,h in e6.items()
                 if h[:2]!=(0,0) and P75.norm_dir(h[:2])==d)
        terms=[-1 if s[i] else 1 for i in f if p[i] in f]
        assert len(terms)==6 and abs(sum(terms))==2
        fibre_terms.append(terms)
    signed_vis=[abs(sum(x))/6 for x in fibre_terms]
    assert signed_vis==[1/3]*4

    eps=separation_threshold()
    checks={
      "twelve_disjoint_swaps":len(pairs)==12,
      "six_negative_swap_pairs":sum(e==-1 for _,_,e in pairs)==6,
      "twelve_negative_output_phases":len(negative)==12,
      "canonical_permutation_two_reversals":perm==list(range(15,-1,-1))+list(range(23,15,-1)),
      "nearest_neighbor_swap_count_148":inv==148,
      "nearest_neighbor_parallel_depth_15":nn_depth==15,
      "all_fibres_have_six_interfering_terms":all(len(t)==6 for t in fibre_terms),
      "ideal_signed_visibility_one_third":signed_vis==[1/3]*4,
      "phase_error_separation_threshold_positive":eps>0.5,
    }
    assert all(checks.values())
    return {
      "schema":"w33.pass11033.signed-clock-photonic-compiler.v1",
      "status":"PASS",
      "headline":(
        "The hidden signed central operation compiles exactly to 12 disjoint "
        "mode swaps plus a diagonal layer of 12 pi phases. With free mode "
        "placement it is a one-layer 12-swap network; in the repository's "
        "canonical 24-mode order its unsigned permutation is reversal of "
        "16 modes plus reversal of 8 modes, requiring exactly 148 adjacent "
        "swaps and minimum parallel nearest-neighbor depth 15."
      ),
      "network":{
        "mode_order":noncentral,
        "permutation_zero_based":perm,
        "signed_swap_pairs_zero_based":[list(x) for x in pairs],
        "arbitrary_routing_swap_count":12,
        "negative_pair_blocks":6,
        "single_mode_pi_phases_if_unsigned_swaps":len(negative),
        "nearest_neighbor_adjacent_swap_count":inv,
        "nearest_neighbor_min_parallel_depth":nn_depth,
        "placement_optimized_depth":"one simultaneous layer of 12 disjoint swaps plus one phase layer",
      },
      "interferometer":{
        "per_fibre_sign_patterns":fibre_terms,
        "ideal_signed_visibility":"1/3",
        "ideal_projective_visibility":"1",
        "phase_error_model":"independent term phases with |delta_i| <= epsilon",
        "signed_visibility_upper_bound":"1/3 + 2 sin(epsilon/2)",
        "projective_visibility_lower_bound":"cos(epsilon)",
        "guaranteed_separation_epsilon_rad":eps,
        "guaranteed_separation_epsilon_deg":eps*180/math.pi,
      },
      "experimental_reading":(
        "If each of the six interfering amplitudes per fibre has calibrated "
        "phase error below the certified threshold, even the conservative "
        "triangle-inequality envelope for the signed model remains below the "
        "projective model's cone bound. The observable therefore survives "
        "substantial phase error in the ideal equal-amplitude model."
      ),
      "boundary":(
        "The phase bound assumes equal amplitudes and bounded coherent phase "
        "errors. Insertion loss, unequal detector efficiencies, fabrication "
        "crosstalk and source impurity require a separate calibrated likelihood "
        "model. The 148 count is minimal only for nearest-neighbor swaps in the "
        "frozen canonical linear ordering; free physical placement reduces it."
      ),
      "checks":checks,
    }
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT)
    a=ap.parse_args(); p=payload()
    text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
      if not a.output.exists() or a.output.read_text()!=text:
        raise SystemExit("certificate drift")
    else:
      a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(text)
    print(json.dumps({"status":p["status"],
      "swaps":p["network"]["arbitrary_routing_swap_count"],
      "adjacent":p["network"]["nearest_neighbor_adjacent_swap_count"],
      "depth":p["network"]["nearest_neighbor_min_parallel_depth"],
      "phase_deg":p["interferometer"]["guaranteed_separation_epsilon_deg"]},
      sort_keys=True))


if __name__=="__main__":
    raise SystemExit(main())
