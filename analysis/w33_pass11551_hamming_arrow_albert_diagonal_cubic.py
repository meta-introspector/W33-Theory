#!/usr/bin/env python3
"""Pass 11551: the Hamming history-arrow cubic is an Albert diagonal norm."""
from __future__ import annotations
import itertools, json, sys
from pathlib import Path
import numpy as np
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
OUT=ROOT/"data/PART_W33_PASS11551_HAMMING_ARROW_ALBERT_DIAGONAL_CUBIC.json"
import w33_pass10950_clock_albert_lorentz_spinor as C

def main():
    J=C.build_clock_albert()
    ids=[J["idx"][g] for g in J["G"]]
    Pi2=np.array([[[int(x*2) for x in J["prod"][i,j]] for j in range(27)] for i in range(27)],dtype=np.int64)
    raw=np.einsum("ijj->i",Pi2); assert not np.any(raw%18)
    trace=raw//18
    gram2=np.einsum("k,ijk->ij",trace,Pi2)
    cube4=np.einsum("l,ial,jka->ijk",trace,Pi2,Pi2)
    norm12=(2*np.einsum("i,j,k->ijk",trace,trace,trace)
            -(np.einsum("i,jk->ijk",trace,gram2)+np.einsum("j,ik->ijk",trace,gram2)+np.einsum("k,ij->ijk",trace,gram2))
            +cube4)
    rows=[]
    for z in itertools.product(range(3),repeat=3):
        y=np.zeros(27,dtype=np.int64)
        for i,a in enumerate(z): y[ids[i]]=a
        v12=int(np.einsum("ijk,i,j,k",norm12,y,y,y)); assert v12%12==0
        N=v12//12; target=z[0]*z[1]*z[2]; assert N==target
        rows.append({"z":list(z),"N_Albert":N,"product":target,"mod3":N%3})
    a,b,c=sp.symbols("a b c")
    H=(2*a+b+c,2*a+2*b+c,2*a+2*c)
    history_poly=sp.Poly(sp.expand(H[0]*H[1]*H[2]),a,b,c,modulus=3).as_expr()
    p11550=json.loads((ROOT/"data/PART_W33_PASS11550_VERONESE_WEIL_HAMMING_ORIENTATION.json").read_text())
    assert p11550["closed_formula"]["Weil_phase_equals"]=="i * chi"
    out={
      "schema":"w33.pass11551.hamming_arrow_albert_diagonal_cubic.v1",
      "status":"PASS_HAMMING_ARROW_IS_ALBERT_DIAGONAL_CUBIC","pass":11551,
      "clock_albert_frame":{"idempotent_indices":ids,"frame_size":3,
        "exact_norm_law":"N_Albert(z1*e1+z2*e2+z3*e3)=z1*z2*z3",
        "all_ternary_vectors_checked":27,"integer_before_mod3":True},
      "history_weld":{"Pass11547_Hsym":["2a+b+c","2a+2b+c","2a+2c"],
        "pulled_back_cubic_mod3":str(history_poly),"null_cube":"all three Hamming coordinates nonzero",
        "arrow_reading":"on the null cube, N mod 3 = 1 is chi=+1 and N mod 3 = 2=-1 is chi=-1",
        "Weil_phase":"i*chi (Pass11550)"},
      "rows":rows,
      "theorem":"The finite Hamming history-arrow polynomial z1*z2*z3 is exactly the restriction of the repository clock-Albert cubic norm to its canonical diagonal three-idempotent frame.",
      "boundary":"Exact diagonal restriction of the committed Albert norm. It does not show that the full 27-event history carrier is the full 27-dimensional Albert algebra, nor derive E6 spacetime dynamics, gravity, or observed couplings."}
    OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":out["status"],"checked":27,"frame":ids},indent=2))
if __name__=="__main__": main()
