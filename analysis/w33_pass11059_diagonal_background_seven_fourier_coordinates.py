#!/usr/bin/env python3
"""Pass 11059: both diagonal weld backgrounds are exactly seven-sparse in K-Fourier coordinates."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_pass11047_phase_weld_exact_54_right_inverse as PREV
import w33_pass11056_diagonal_weld_ten_triangle_factors as FCT

OUT=ROOT/"data/w33_pass11059_diagonal_background_seven_fourier_coordinates.json"

def payload():
    bridge=json.loads((ROOT/"data/w33_e6id_current_h27_gauge_bridge.json").read_text())
    h={int(i):tuple(map(int,x)) for i,x in bridge["maps"]["e6id_to_current_H27_address"].items()}
    rows={}
    for name in ("plus","minus"):
        coords,v,D=FCT.diagonal(name)
        cert=[]
        for p in (103,109):
            F,Finv,labels=PREV.full_fourier_matrix(p,h)
            y=(Finv@(np.array(v,dtype=np.int64)%p))%p
            w=PREV.root_omega(p); inv3=pow(3,-1,p)
            t1=1 if name=="plus" else 2
            t2=2 if name=="plus" else 1
            expected={}
            a=((w-1)*inv3)%p
            b=((pow(w,2,p)-1)*inv3)%p
            for r in range(3): expected[("S1",t1,r,0)]=a
            for r in range(3): expected[("S2",t2,r,0)]=b
            expected[("L",0,0,0)]=2%p
            got={labels[i]:int(y[i]) for i in range(81) if y[i]%p}
            assert got==expected
            cert.append({"prime":p,"omega":w,"nonzero_coordinates":7,"support":[list(x) for x in expected]})
        rows[name]=cert
    return {
      "schema":"w33.pass11059.diagonal-background-seven-fourier-coordinates.v1",
      "status":"PASS_DIAGONAL_WELD_BACKGROUND_IS_EXACTLY_SEVEN_SPARSE_IN_THE_K81_FOURIER_BASIS",
      "headline":"The three-level diagonal background itself is exceptionally sparse in the same K-Fourier basis used by the compiler. For tau=c+p it has only seven nonzero coordinates: three S1 modes with external character t=1, three S2 modes with t=2, and the trivial L mode. The c-p orientation swaps t=1 and t=2. The cyclotomic coefficients are exact.",
      "exact_formula":{
        "center_plus_external":"v_plus = 2 L(0,0,0) + ((omega-1)/3) sum_r S1(1,r,0) + ((omega^2-1)/3) sum_r S2(2,r,0)",
        "center_minus_external":"v_minus = 2 L(0,0,0) + ((omega-1)/3) sum_r S1(2,r,0) + ((omega^2-1)/3) sum_r S2(1,r,0)",
        "nonzero_coordinates":7
      },
      "split_prime_certificates":rows,
      "structural_reading":"The background that unlocks all 54 retyped directions is not Fourier-generic. It is a coherent seven-mode seed coupling the trivial abelian mode to conjugate Schrodinger clock sectors.",
      "boundary":"This is an exact coordinate identity in the frozen K-Fourier gauge; it is not a dynamical vacuum-selection statement.",
      "parents":["data/w33_pass11047_phase_weld_exact_54_right_inverse.json","data/w33_pass11058_triangle_factor_kernel_interference.json"],
      "checks":{"plus_support7":True,"minus_support7":True,"cyclotomic_formula_verified_at_103_and_109":True}
    }

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args()
    p=payload();text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or a.output.read_text()!=text:raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"support":7},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
