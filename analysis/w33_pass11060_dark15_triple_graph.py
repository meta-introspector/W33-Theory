#!/usr/bin/env python3
"""Pass 11060: the diagonal-weld dark 15 is transverse to all three K-Fourier sectors."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_pass11047_phase_weld_exact_54_right_inverse as PREV
import w33_pass11056_diagonal_weld_ten_triangle_factors as FCT

OUT=ROOT/"data/w33_pass11060_dark15_triple_graph.json"

def nullspace_basis(A,p):
    M=np.array(A,dtype=np.int64)%p; n,m=M.shape; r=0;piv=[]
    for c in range(m):
        z=next((i for i in range(r,n) if M[i,c]%p),None)
        if z is None: continue
        M[[r,z]]=M[[z,r]]
        iv=pow(int(M[r,c]%p),-1,p);M[r]=(M[r]*iv)%p
        for i in range(n):
            if i!=r and M[i,c]%p:
                q=int(M[i,c]);M[i]=(M[i]-q*M[r])%p
        piv.append(c);r+=1
        if r==n:break
    free=[c for c in range(m) if c not in piv]
    B=[]
    for f in free:
        x=np.zeros(m,dtype=np.int64);x[f]=1
        for i,c in enumerate(piv):x[c]=(-M[i,f])%p
        B.append(x)
    return np.array(B,dtype=np.int64).T

def payload():
    bridge=json.loads((ROOT/"data/w33_e6id_current_h27_gauge_bridge.json").read_text())
    h={int(i):tuple(map(int,x)) for i,x in bridge["maps"]["e6id_to_current_H27_address"].items()}
    cert={}
    for name in ("plus","minus"):
        coords,v,D=FCT.diagonal(name); rows=[]
        for p in (103,109):
            F,Finv,_=PREV.full_fourier_matrix(p,h)
            A=(Finv@(D%p)@F)%p
            assert PREV.inv_rank(A,p,False)[0]==66
            N=nullspace_basis(A,p);assert N.shape==(81,15)
            proj=[FCT.rank_mod(N[a:b,:],p) for a,b in ((0,27),(27,54),(54,81))]
            pure=[27-PREV.inv_rank(A[:,a:b],p,False)[0] for a,b in ((0,27),(27,54),(54,81))]
            assert proj==[15,15,15] and pure==[0,0,0]
            rows.append({"prime":p,"kernel_dimension":15,"projection_ranks":{"S1":15,"S2":15,"L":15},"pure_sector_kernel_dimensions":{"S1":0,"S2":0,"L":0}})
        cert[name]=rows
    return {
      "schema":"w33.pass11060.dark15-triple-graph.v1",
      "status":"PASS_THE_DIAGONAL_WELD_DARK15_IS_A_TRIPLE_GRAPH_TRANSVERSE_TO_S1_S2_AND_L",
      "headline":"The 15D kernel of either diagonal rank-66 generator contains no nonzero vector lying purely in S1, S2, or L. Its projection onto each 27D Fourier sector has full rank 15 at both split primes. Thus the dark space is a genuinely mixed triple-graph subspace: every one of its 15 degrees of freedom is simultaneously visible in all three compiler sectors.",
      "sector_geometry":{
        "kernel_dimension":15,
        "pure_S1_kernel_dimension":0,"pure_S2_kernel_dimension":0,"pure_L_kernel_dimension":0,
        "projection_rank_to_S1":15,"projection_rank_to_S2":15,"projection_rank_to_L":15,
        "interpretation":"the kernel is the graph of injective linear maps between its 15D projections in any chosen pair of sectors"
      },
      "split_prime_certificates":cert,
      "connection_to_background_ray":"Pass 11058 isolates one common local-kernel ray v; Pass 11059 shows that ray already mixes S1, S2 and L. The remaining 14 dark dimensions are likewise not confined to any single Fourier sector.",
      "boundary":"Transversality is a finite linear-algebra statement in the frozen Fourier gauge. It does not by itself identify a protected physical code or decoherence-free subsystem.",
      "parents":["data/w33_pass11058_triangle_factor_kernel_interference.json","data/w33_pass11059_diagonal_background_seven_fourier_coordinates.json"],
      "checks":{"both_orientations":True,"kernel15":True,"no_pure_sector_dark_vectors":True,"projection_rank15_to_every_sector":True}
    }

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args()
    p=payload();text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or a.output.read_text()!=text:raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"kernel":15,"projections":[15,15,15]},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
