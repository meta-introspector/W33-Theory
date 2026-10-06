#!/usr/bin/env python3
"""Pass 11559: fixed-dimensional gauge/frame refinement of the event Dirac operator.

Use the tower (Z/3^r Z)^3 with the same eight body-diagonal steps.  The scalar
null Laplacian has a diffusive N^-2 gap at fixed dimension.  The vector
operators B_i have linear covariant symbols, and after parallel transport their
commutators expose gauge curvature at leading order.  A constant frame E turns
Q into a Dirac principal symbol for metric g=E^T E.

Variable frames reveal the remaining gravity gap explicitly: derivatives of E
appear, but no Levi-Civita/spin connection is selected automatically.
"""
from __future__ import annotations
import itertools, json, math, sys
from pathlib import Path
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/PART_W33_PASS11559_FIXED_DIMENSION_GAUGE_FRAME_REFINEMENT.json"
sys.path.insert(0,str(ROOT/"analysis"))
import w33_pass11557_pin_equivariant_event_dirac as P57

EPS=list(itertools.product((-1,1),repeat=3))

def moment(i,word):
    return sum(e[i]*math.prod(e[j] for j in word) for e in EPS)

def main():
    m0=[moment(i,()) for i in range(3)]
    m1=[[moment(i,(j,)) for j in range(3)] for i in range(3)]
    m2=[[[moment(i,(j,k)) for k in range(3)] for j in range(3)] for i in range(3)]
    assert m0==[0,0,0]
    assert m1==[[8 if i==j else 0 for j in range(3)] for i in range(3)]
    assert all(x==0 for A in m2 for B in A for x in B)

    cubic={}
    for i in range(3):
        rows={}
        for word in itertools.product(range(3),repeat=3):
            v=moment(i,word)
            if v: rows["".join(map(str,word))]=v
        cubic[str(i)]=rows

    # Fixed-dimensional normalized spectral gap for N=3^r.
    gaps=[]
    for r in range(1,7):
        N=3**r
        gap=1-math.cos(2*math.pi/N)
        gaps.append({
          "r":r,"N":N,"normalized_gap":gap,
          "N_squared_gap":N*N*gap,
          "target_2pi2":2*math.pi*math.pi,
        })
    assert abs(gaps[-1]["N_squared_gap"]-2*math.pi*math.pi)<1e-3

    # Exact constant-frame principal symbol.
    gam,_=P57.small_clifford()
    p=sp.symbols("p0:3",commutative=True)
    E=sp.Matrix([[1,1,0],[0,1,1],[1,0,1]])
    G=E.T*E
    symbol=sp.zeros(4)
    for a in range(3):
        q=sum(E[a,i]*p[i] for i in range(3))
        symbol += gam[a]*q
    scalar=(sp.Matrix(p).T*G*sp.Matrix(p))[0]
    assert sp.simplify(symbol*symbol-scalar*sp.eye(4))==sp.zeros(4)

    out={
      "schema":"w33.pass11559.fixed-dimension-gauge-frame-refinement.v1",
      "status":"PASS_FIXED_DIMENSION_COVARIANT_DIRAC_PRINCIPAL_SYMBOL","pass":11559,
      "tower":{
        "sites":"(Z/3^r Z)^3",
        "dimension":3,
        "local_steps":"epsilon in {+/-1}^3",
        "normalized_null_adjacency_symbol":"prod_i cos(theta_i)",
        "normalized_gap_exact":"1-cos(2*pi/N), N=3^r",
        "gap_asymptotic":"2*pi^2/N^2 + O(N^-4)",
        "rows":gaps,
      },
      "vector_difference_expansion":{
        "definition":"B_i(a)=sum_epsilon epsilon_i exp(a epsilon^j D_j)",
        "zeroth_moments":m0,
        "first_moments":m1,
        "second_moments":m2,
        "third_moments_nonzero":cubic,
        "principal_law":"B_i(a)=8 a D_i + O(a^3)",
        "normalized_law":"B_i/(8a)=D_i+O(a^2)",
      },
      "gauge_covariance":{
        "link_reading":"replace exp(a epsilon^j partial_j) by parallel-transported exp(a epsilon^j D_j)",
        "commutator_law":"[B_i,B_j]=64 a^2 [D_i,D_j]+O(a^4)",
        "U1_specialization":"[D_i,D_j]=i F_ij",
        "curvature_reading":"the first nontrivial commutator of vector event differences measures gauge curvature",
      },
      "frame_principal_symbol":{
        "operator":"Q_E/(8a)=i gamma^a E_a^i D_i + O(a^2)",
        "sample_integer_frame":[[int(x) for x in row] for row in E.tolist()],
        "metric_g_equals_EtE":[[int(x) for x in row] for row in G.tolist()],
        "exact_symbol_square_check":True,
        "constant_frame_square":"principal scalar part = -g^{ij}D_iD_j; antisymmetric Clifford part = -(1/4)[gamma^a,gamma^b]E_a^i E_b^j [D_i,D_j]",
      },
      "variable_frame_firewall":{
        "extra_term":"-gamma^a gamma^b E_a^i (partial_i E_b^j) D_j appears already at principal squaring order",
        "conclusion":"a varying frame alone does not select the Levi-Civita/spin connection; torsion/connection dynamics remain additional structure",
      },
      "theorem":"The fixed-dimensional 3-adic-resolution tower avoids the H(n,3) dimension-growth firewall: its scalar gap closes as N^-2, while the Pin event vector differences converge to covariant first derivatives. Their Clifford contraction has the Dirac principal symbol for g=E^T E and exposes gauge curvature through commutators.",
      "boundary":"Controlled local/asymptotic operator statement only. The tower is not yet a derived physical continuum, the connection and frame dynamics are supplied rather than selected, and no Einstein-Hilbert action, Newton scale, or Lorentzian signature is derived."
    }
    OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":out["status"],"Nmax":gaps[-1]["N"],"N2gap":gaps[-1]["N_squared_gap"]},indent=2))

if __name__=="__main__": main()
