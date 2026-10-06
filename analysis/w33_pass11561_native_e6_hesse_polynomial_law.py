#!/usr/bin/env python3
"""Pass 11561: tensor-level native E6 / Albert / Hesse polynomial-law weld.

Pass11551 identified xyz with the Albert norm on the clock Jordan frame.
Pass11471 supplies the exact Gaussian-integer carrier map T from that clock
Albert basis to the repository native signed E6 cubic tensor.  Pulling the
native tensor back to the three idempotents gives coefficient 1 on exactly
the six permutations of (0,1,2), hence native cubic contraction 6*x*y*z.

The Pass11114 Hesse flavor determinant then has, in characteristic three,
the coefficient-level polynomial-law decomposition

 H = (s-d)^3 N - s d^2 Ttrace^3,

with N=xyz and Ttrace=x+y+z.  The Frobenius trace-cube has zero mixed
third cross-effect in characteristic 3, so it must not be conflated with an
ordinary polarized symmetric trilinear tensor.  The two-plane <N,Ttrace^3>
breaks the W(D3) frame normalizer from 24 to the pure S3 permutations.
"""
from __future__ import annotations
import itertools, json, sys
from pathlib import Path
import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/PART_W33_PASS11561_NATIVE_E6_HESSE_POLYNOMIAL_LAW.json"
sys.path.insert(0,str(ROOT/"analysis"))
from w33_pass11384_11388_native_dynamics import native_tensors

def signed_permutations():
    out=[]
    for pm in itertools.permutations(range(3)):
        for sg in itertools.product((-1,1),repeat=3):
            M=np.zeros((3,3),dtype=int)
            for i,j in enumerate(pm): M[i,j]=sg[i]
            out.append(M)
    return out

def rank3(cols):
    A=np.array(cols,dtype=int).T%3
    r=0
    for c in range(A.shape[1]):
        q=next((i for i in range(r,A.shape[0]) if A[i,c]),None)
        if q is None: continue
        A[[r,q]]=A[[q,r]]
        A[r]=(A[r]*pow(int(A[r,c]),-1,3))%3
        for i in range(A.shape[0]):
            if i!=r and A[i,c]: A[i]=(A[i]-A[i,c]*A[r])%3
        r+=1
    return r

def degree3_vector(poly,x,y,z):
    P=sp.Poly(sp.expand(poly),x,y,z,modulus=3)
    mons=[(3,0,0),(2,1,0),(2,0,1),(1,2,0),(1,1,1),
          (1,0,2),(0,3,0),(0,2,1),(0,1,2),(0,0,3)]
    return np.array([int(P.coeff_monomial(x**a*y**b*z**c))%3 for a,b,c in mons],dtype=int)

def cross3_mod3(poly,x,y,z):
    # Third finite cross-effect on e1,e2,e3 at the origin.
    def val(a,b,c):
        return int(sp.expand(poly.subs({x:a,y:b,z:c})))%3
    ans=0
    for mask in range(8):
        a=(mask>>0)&1;b=(mask>>1)&1;c=(mask>>2)&1
        bits=a+b+c
        ans += ((-1)**(3-bits))*val(a,b,c)
    return ans%3

def main():
    bridge=json.loads((ROOT/"data/w33_pass11471_11475_frames_currents_matching.json").read_text())
    ab=bridge["frame"]["albert_bridge"]
    Traw=ab["native_from_clock_basis"]
    Tr=np.array(Traw["real"],dtype=float); Ti=np.array(Traw["imag"],dtype=float)
    assert np.array_equal(Tr,np.rint(Tr)) and np.array_equal(Ti,np.rint(Ti))
    T=np.rint(Tr).astype(np.int64)+1j*np.rint(Ti).astype(np.int64)
    ids=list(map(int,ab["clock_idempotent_indices"]))
    assert ids==[26,0,9]
    assert ab["cubic_equality"]=="12*N_native(T*y)=12*N_Albert(y)"

    _,tensors=native_tensors()
    native=tensors[1]
    V=T[:,ids]
    C=np.einsum("abc,ai,bj,ck->ijk",native,V,V,V,optimize=True)
    assert np.max(np.abs(C.imag))==0
    Cr=np.rint(C.real).astype(np.int64)
    expected=np.zeros((3,3,3),dtype=np.int64)
    for p in itertools.permutations(range(3)): expected[p]=1
    assert np.array_equal(Cr,expected)

    # Exact point controls: contraction is 6 xyz over the integral carrier.
    for zc in itertools.product(range(3),repeat=3):
        v=V@np.array(zc,dtype=np.int64)
        val=np.einsum("abc,a,b,c",native,v,v,v)
        assert val.imag==0 and int(round(val.real))==6*zc[0]*zc[1]*zc[2]

    x,y,z,s,d=sp.symbols("x y z s d")
    N=x*y*z
    Fermat=x**3+y**3+z**3
    Ttrace=x+y+z
    H=(s**3+2*d**3)*N-s*d**2*Fermat
    target=(s-d)**3*N-s*d**2*Ttrace**3
    diff=sp.Poly(sp.expand(H-target),x,y,z,s,d,modulus=3)
    assert diff.is_zero

    # Polarization/cross-effect firewall in characteristic 3.
    n_cross=cross3_mod3(N,x,y,z)
    trace_cross=cross3_mod3(Ttrace**3,x,y,z)
    assert n_cross==1 and trace_cross==0

    Nvec=degree3_vector(N,x,y,z)
    Tvec=degree3_vector(Ttrace**3,x,y,z)
    assert rank3([Nvec,Tvec])==2

    # Stabilizer of the polynomial-law two-plane under signed permutations.
    O=signed_permutations(); assert len(O)==48
    stabilizer_O=[]; stabilizer_W=[]
    for M in O:
        subs={
          x:int(M[0,0])*x+int(M[0,1])*y+int(M[0,2])*z,
          y:int(M[1,0])*x+int(M[1,1])*y+int(M[1,2])*z,
          z:int(M[2,0])*x+int(M[2,1])*y+int(M[2,2])*z,
        }
        Nv=degree3_vector(N.subs(subs,simultaneous=True),x,y,z)
        Tv=degree3_vector((Ttrace**3).subs(subs,simultaneous=True),x,y,z)
        stable=(rank3([Nvec,Tvec,Nv,Tv])==2)
        if stable: stabilizer_O.append(M)
        signprod=int(round(np.prod([next(v for v in row if v) for row in M])))
        if signprod==1 and stable: stabilizer_W.append(M)
    assert len(stabilizer_O)==12
    assert len(stabilizer_W)==6
    # The W-stabilizer is exactly the six unsigned coordinate permutations.
    assert all(set(M.ravel())<={0,1} for M in stabilizer_W)

    p51=json.loads((ROOT/"data/PART_W33_PASS11551_HAMMING_ARROW_ALBERT_DIAGONAL_CUBIC.json").read_text())
    assert p51["clock_albert_frame"]["exact_norm_law"]=="N_Albert(z1*e1+z2*e2+z3*e3)=z1*z2*z3"

    out={
      "schema":"w33.pass11561.native-e6-hesse-polynomial-law.v1",
      "status":"PASS_NATIVE_E6_HESSE_POLYNOMIAL_LAW_AND_STABILIZER","pass":11561,
      "native_carrier":{
        "Pass11471_map":"exact Gaussian-integer 27x27 native_from_clock_basis",
        "clock_idempotent_indices":ids,
        "map_entries_integral":True,
        "restricted_native_tensor":Cr.tolist(),
        "restricted_tensor_rule":"C_ijk=1 iff i,j,k are all distinct, else 0",
        "native_cubic_contraction":"d_native(T_frame z,T_frame z,T_frame z)=6 z1 z2 z3",
        "Albert_weld":"therefore (1/6) native contraction equals the Pass11551 Albert diagonal norm xyz over characteristic zero",
      },
      "char3_Hesse":{
        "Pass11114_integral_formula":"H=(s^3+2d^3)xyz-s d^2(x^3+y^3+z^3)",
        "polynomial_law_identity_mod3":"H=(s-d)^3 N - s d^2 T^3, N=xyz, T=x+y+z",
        "coefficientwise_mod3_verified":True,
        "N_third_cross_effect_on_e1_e2_e3":n_cross,
        "T_cubed_third_cross_effect_on_e1_e2_e3":trace_cross,
        "category_firewall":"in characteristic 3, T^3 is Frobenius/additive and its mixed third polarization vanishes; the correct statement is equality of cubic polynomial laws, not equality of ordinary symmetric trilinear tensors obtained by division by 3!.",
      },
      "frame_normalizer":{
        "full_signed_permutation_order":48,
        "WD3_order":24,
        "stabilizer_of_span_N_T3_in_full_signed_group":len(stabilizer_O),
        "stabilizer_of_span_N_T3_in_WD3":len(stabilizer_W),
        "WD3_stabilizer":"pure coordinate permutations S3",
        "symmetry_breaking":"W(D3) (24) -> S3 (6) when the Frobenius trace-cube direction is retained alongside the Albert norm",
      },
      "theorem":"The committed native E6 cubic, transported through the exact Pass11471 clock-to-native map, restricts coefficientwise to 6xyz on the same three Albert idempotents that carry the history arrow. The Pass11114 Hesse determinant reduces in characteristic three to an exact polynomial-law two-plane generated by that Albert norm and the Frobenius trace cube; this extra direction reduces the W(D3) frame normalizer from 24 to S3 of order 6.",
      "boundary":"This is a tensor/polynomial-law carrier theorem. It does not make the full Hesse flavor determinant E6-invariant, does not derive observed Yukawas or masses, and does not identify characteristic-three reduction with a physical low-energy limit."
    }
    OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":out["status"],"native_factor":6,"stab_W":len(stabilizer_W),"stab_O":len(stabilizer_O)},indent=2))

if __name__=="__main__": main()
