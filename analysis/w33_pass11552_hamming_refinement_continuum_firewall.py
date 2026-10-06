#!/usr/bin/env python3
"""Pass 11552: Hamming refinement tower and continuum firewall."""
from __future__ import annotations
import json, math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/PART_W33_PASS11552_HAMMING_REFINEMENT_CONTINUUM_FIREWALL.json"

def K(n,j,w):
    return sum(((-1)**h)*(2**(j-h))*math.comb(w,h)*math.comb(n-w,j-h)
               for h in range(j+1) if h<=w and j-h<=n-w)

def full_shell(n):
    js=list(range(3,n+1,3))
    degree=sum(math.comb(n,j)*2**j for j in js)
    vals=[sum(K(n,j,w) for j in js) for w in range(n+1)]
    osc=0 if n%2 else 2*((-1)**(n//2))*3**(n//2)
    closed=(3**n+osc)//3-1
    assert degree==closed
    ratio=max(abs(v) for v in vals[1:])/degree
    return degree,vals,ratio

def main():
    rows=[]
    for n in range(3,41):
        deg,vals,ratio=full_shell(n)
        k3=8*math.comb(n,3); v1=K(n,3,1)
        branch_gap=1-v1/k3
        assert abs(branch_gap-9/(2*n))<1e-15
        ordinary=max(K(n,3,w) for w in range(1,n+1))
        ordinary_gap=1-ordinary/k3
        if n>=5: assert abs(ordinary_gap-9/(2*n))<1e-15
        bound=(1+2*(3**(n/2))/3)/deg
        assert ratio <= bound*(1+1e-12)
        rows.append({"n":n,"full_null_degree":deg,"full_null_max_abs_ratio":ratio,
          "full_null_absolute_gap":1-ratio,"root_filter_ratio_bound":bound,
          "distance3_degree":k3,"distance3_weight1_gap":branch_gap,
          "distance3_ordinary_gap":ordinary_gap})
    out={"schema":"w33.pass11552.hamming_refinement_continuum_firewall.v1",
      "status":"PASS_NAIVE_FULL_NULL_TOWER_IS_EXPANDER_LIKE","pass":11552,
      "exact_full_null_relation":{"carrier":"F3^n","connection":"nonzero vectors of Hamming weight divisible by 3",
        "degree":"(3^n + 2 Re((i sqrt(3))^n))/3 - 1",
        "character_eigenvalue":"lambda_w=(P_w(1)+P_w(omega)+P_w(omega^2))/3 - 1, P_w(t)=(1+2t)^(n-w)(1-t)^w",
        "nontrivial_bound":"|lambda_w| <= 1 + 2*3^(n/2-1)",
        "consequence":"max nontrivial |lambda|/degree -> 0 exponentially; normalized absolute spectral gap -> 1"},
      "nearest_null_distance3":{"degree":"8*C(n,3)","weight1_normalized_gap":"9/(2n)",
        "verified_as_ordinary_spectral_gap_for_n_at_least_5_through":40,
        "consequence":"this branch gap closes as 1/n, unlike the all-null relation"},
      "rows":rows,
      "continuum_verdict":"The naive all-q=0 Hamming tower is not a diffusive continuum refinement. A distance-3-only scaling can close a spectral gap, but increasing n changes ambient dimension and does not by itself produce fixed 3+1 Lorentzian spacetime.",
      "boundary":"Spectral theorem/firewall for a specific Hamming refinement family. It neither proves impossibility of every continuum construction nor supplies Lorentz symmetry, a metric, Einstein dynamics, or a physical lattice spacing."}
    OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":out["status"],"n40_abs_gap":rows[-1]["full_null_absolute_gap"],"distance3_gap_n40":rows[-1]["distance3_weight1_gap"]},indent=2))
if __name__=="__main__": main()
