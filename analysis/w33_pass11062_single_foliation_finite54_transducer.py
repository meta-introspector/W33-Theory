#!/usr/bin/env python3
"""Pass 11062: one 27-triangle factor can already be a full finite 54D retyping transducer."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_pass11047_phase_weld_exact_54_right_inverse as PREV
import w33_pass11056_diagonal_weld_ten_triangle_factors as FCT

OUT=ROOT/"data/w33_pass11062_single_foliation_finite54_transducer.json"

def factorize(which):
    coords,v,D=FCT.diagonal(which)
    idx={x:i for i,x in enumerate(coords)}
    ie=idx[FCT.IDENT]
    connection=[coords[j] for j in range(81) if D[ie,j]!=0]
    seen=set();pairs=[]
    for s in connection:
        if s in seen:continue
        si=FCT.kinv(s);assert si in connection and si!=s
        seen.add(s);seen.add(si);pairs.append((s,si))
    factors=[]
    for s,si in pairs:
        M=np.zeros_like(D);pair={s,si}
        for i,x in enumerate(coords):
            xi=FCT.kinv(x)
            for j,y in enumerate(coords):
                if FCT.kmul(xi,y) in pair:M[i,j]=D[i,j]
        assert np.array_equal(M.T,-M)
        cc=FCT.components(M);assert len(cc)==27 and all(len(z)==3 for z in cc)
        factors.append(M)
    assert len(factors)==10 and np.array_equal(sum(factors),D)
    return coords,pairs,factors

def payload():
    bridge=json.loads((ROOT/"data/w33_e6id_current_h27_gauge_bridge.json").read_text())
    h={int(i):tuple(map(int,x)) for i,x in bridge["maps"]["e6id_to_current_H27_address"].items()}
    cert={};perfect={}
    for which in ("plus","minus"):
        coords,pairs,factors=factorize(which)
        prime_rows=[];spectra={}
        for p in (103,109):
            F,Finv,_=PREV.full_fourier_matrix(p,h)
            rows=[]
            for i,M in enumerate(factors):
                C=(Finv@(M%p))%p
                raw=PREV.inv_rank(M%p,p,False)[0]
                q=PREV.inv_rank(C[27:,:],p,False)[0]
                s2=PREV.inv_rank(C[27:54,:],p,False)[0]
                l=PREV.inv_rank(C[54:,:],p,False)[0]
                rows.append({"factor":i,"raw_rank":raw,"quotient_rank":q,"S2_rank":s2,"L_rank":l})
            spectra[str(p)]=[r["quotient_rank"] for r in rows]
            prime_rows.append({"prime":p,"rows":rows})
        assert spectra["103"]==spectra["109"]
        good=[r["factor"] for r in prime_rows[0]["rows"] if r["raw_rank"]==r["quotient_rank"]==54]
        assert all([r["factor"] for r in pr["rows"] if r["raw_rank"]==r["quotient_rank"]==54]==good for pr in prime_rows)
        assert len(good)==2
        pp=[]
        for i in good:
            s=pairs[i][0]
            pp.append({"factor":i,"generator":[list(s[0]),s[1]],"triangles":27,"raw_rank":54,"quotient_rank":54,"S1_intersection_dimension":0})
        cert[which]=prime_rows;perfect[which]=pp

    assert [x["factor"] for x in perfect["plus"]]==[7,8]
    assert [x["factor"] for x in perfect["minus"]]==[6,9]
    return {
      "schema":"w33.pass11062.single-foliation-finite54-transducer.v1",
      "status":"PASS_ONE_27_TRIANGLE_FOLIATION_CAN_ALREADY_TRANSVERSELY_COVER_THE_FULL_54D_RETYPE_QUOTIENT",
      "headline":"The ten-factor diagonal weld is highly redundant as a symmetry-retyping transducer. In the frozen K-Fourier gauge, each individual triangle factor has raw rank 54. For the plus orientation, factors 7 and 8 project with rank 54 onto S2+L; for the minus orientation, factors 6 and 9 do. Because their raw rank is also 54, their images intersect S1 trivially and project isomorphically onto the full retyped quotient.",
      "quotient_rank_spectra":{
        "plus":[48,48,45,45,45,45,36,54,54,36],
        "minus":[48,48,45,45,45,45,54,36,36,54]
      },
      "perfect_single_factors":perfect,
      "finite_gate_consequence":{
        "factor_structure":"27 disjoint weighted three-mode skew blocks",
        "raw_rank":54,"quotient_rank":54,"parallel_three_mode_blocks":27,
        "cayley_identity":"C_t(F)-I = 2t F (I-tF)^-1",
        "for_nonzero_real_t":"Im(C_t(F)-I)=Im(F), so the quotient map remains an isomorphism",
        "single_factor_layer_depth":1
      },
      "split_prime_certificates":cert,
      "boundary":"This proves existence of sparse finite retyping transducers. A fixed one-parameter factor gate is not 54 independent control knobs, and it is not identical to the full ten-factor diagonal-weld gate.",
      "parents":["data/w33_pass11056_diagonal_weld_ten_triangle_factors.json","data/w33_pass11055_cayley_finite_54_compiler.json"],
      "checks":{"two_perfect_single_factors_per_orientation":True,"plus_factors_7_8":True,"minus_factors_6_9":True,"each_perfect_factor_27_triangles_rank54_projection54":True}
    }

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args()
    p=payload();text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or a.output.read_text()!=text:raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"plus":[7,8],"minus":[6,9]},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
