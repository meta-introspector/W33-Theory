#!/usr/bin/env python3
"""Cocycle-transport the sparse 54D foliation transducers from K81 to U81."""
from __future__ import annotations
import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_pass11047_phase_weld_exact_54_right_inverse as P47
import w33_pass11062_single_foliation_finite54_transducer as P62
import w33_pass11067_explicit_h27_central_cocycle_weld as P67
import w33_pass11068_ternary_cocycle_deformation as P68

OUT=ROOT/"data/w33_20261001_u81_54d_cocycle_transducer.json"

def hpow(h,k):
    z=(0,0,0)
    for _ in range(k): z=P67.hmul(z,h)
    return z

def quotient_coset_reps(hg):
    H=[(a,b,c) for a in range(3) for b in range(3) for c in range(3)]
    unseen=set(H); reps=[]
    while unseen:
        r=min(unseen)
        cos={P67.hmul(r,hpow(hg,k)) for k in range(3)}
        reps.append(r); unseen-=cos
    assert len(reps)==9
    return reps
def canonical_transport(g,s,coords):
    """Map split right-g cycles to star_s cycles with the same H27 coset and central seed."""
    idx={tuple(h)+(d,):i for i,(h,d) in enumerate(coords)}
    perm=np.full(81,-1,dtype=int); used=set()
    for rh in quotient_coset_reps(g[:3]):
        for d0 in range(3):
            xs=tuple(rh)+(d0,); xt=xs
            for _ in range(3):
                i=idx[xs]; j=idx[xt]
                assert perm[i] in (-1,j)
                perm[i]=j; used.add(j)
                xs=P68.mul(xs,g,0); xt=P68.mul(xt,g,s)
    assert (perm>=0).all() and len(used)==81
    return perm

def conjugate_factor(F0,perm):
    F=np.zeros_like(F0)
    ii,jj=np.nonzero(F0)
    for a,b in zip(ii,jj): F[perm[a],perm[b]]=F0[a,b]
    assert np.array_equal(F.T,-F)
    return F

def payload():
    bridge=json.loads((ROOT/"data/w33_e6id_current_h27_gauge_bridge.json").read_text())
    hmap={int(i):tuple(map(int,x)) for i,x in bridge["maps"]["e6id_to_current_H27_address"].items()}
    selected={"plus":[7,8],"minus":[6,9]}
    rows=[]
    for which,good in selected.items():
        coords,pairs,factors=P62.factorize(which)
        for factor in good:
            k=pairs[factor][0]; g=tuple(k[0])+(k[1],)
            orders=[P68.order(g,s) for s in (0,1,2)]
            assert orders==[3,3,3]
            for s in (1,2):
                perm=canonical_transport(g,s,coords)
                F=conjugate_factor(factors[factor],perm)
                cc=P62.FCT.components(F)
                assert len(cc)==27 and all(len(z)==3 for z in cc)
                primes=[]; pivot0=None
                for p in (103,109):
                    FF,Finv,_=P47.full_fourier_matrix(p,hmap)
                    C=(Finv@(F%p))%p; Q=C[27:,:]
                    raw=P47.inv_rank(F%p,p,False)[0]
                    q,piv,_=P47.inv_rank(Q,p,False)
                    s2=P47.inv_rank(Q[:27,:],p,False)[0]
                    lat=P47.inv_rank(Q[27:,:],p,False)[0]
                    assert (raw,q,s2,lat)==(54,54,27,27)
                    B=Q[:,piv]; det=P47.det_mod(B,p); assert det
                    if pivot0 is None: pivot0=piv
                    assert piv==pivot0
                    primes.append({"prime":p,"raw_rank":raw,"quotient_rank":q,
                                   "S2_rank":s2,"L_rank":lat,
                                   "pivot_minor_det_mod_p":det})
                rows.append({
                    "orientation":which,"factor":factor,"generator":list(g),
                    "cocycle_s":s,"generator_orders_s0_s1_s2":orders,
                    "triangles":27,"vertices_per_triangle":3,
                    "transport_permutation_sha256":hashlib.sha256(perm.astype(np.int64).tobytes()).hexdigest(),
                    "stable_pivot_columns":list(map(int,pivot0)),
                    "prime_replays":primes,
                })
    assert len(rows)==8
    return {
      "schema":"w33.20261001.u81-54d-cocycle-transducer.v1",
      "status":"PASS_SPARSE_54D_FOLIATION_TRANSDUCERS_SURVIVE_THE_K81_TO_U81_COCYCLE_WELD",
      "theorem":(
        "The four split K81 single-foliation transducers that are transverse to S2+L remain "
        "order-three 27-triangle factors for both nonzero cocycle twists. Canonical transport "
        "along the common H27 quotient conjugates each weighted split factor to the corresponding "
        "U81 Cayley foliation. Every transported factor has raw rank 54 and S2+L projection rank "
        "54, with S2 and L ranks 27+27, at both split primes 103 and 109."
      ),
      "transport":"same H27 right-coset representative, same central seed, same position in the g-cycle",
      "rows":rows,
      "exactness":(
        "Each transported factor is permutation-conjugate to the exact integer split factor, so "
        "its raw rank and 27 weighted triangle blocks are exact over characteristic zero. The same "
        "54 pivot columns give nonzero minors modulo 103 and 109, proving characteristic-zero "
        "transversality over Q(omega) for the frozen S2+L target."
      ),
      "finite_gate_consequence":(
        "For nonzero real t, the Cayley gate C_t(F_s) has Im(C_t(F_s)-I)=Im(F_s); therefore "
        "the finite one-layer 27-block transducer retains the full 54D retyping image after the "
        "class-three cocycle weld."
      ),
      "boundary":(
        "The transport is a canonical carrier identification over the common H27 quotient, not a "
        "group isomorphism K81~=U81. It proves a finite linear transducer compatible with the "
        "nonsplit chamber labels; it does not by itself specify an optical Hamiltonian."
      ),
      "parents":[
        "data/w33_pass11062_single_foliation_finite54_transducer.json",
        "data/w33_pass11067_explicit_h27_central_cocycle_weld.json",
        "data/w33_pass11068_ternary_cocycle_deformation.json"
      ],
      "checks":{"eight_nonsplit_replays":True,"all_order3":True,"all_27_triangles":True,
                "all_rank54":True,"all_S2_L_27_27":True,"stable_pivots_two_primes":True}
    }
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--check",action="store_true")
    a=ap.parse_args(); p=payload(); txt=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not OUT.exists() or OUT.read_text()!=txt: raise SystemExit("certificate drift")
    else:
        OUT.write_text(txt)
    print(json.dumps({"status":p["status"],"replays":len(p["rows"])},sort_keys=True))

if __name__=="__main__": raise SystemExit(main())
