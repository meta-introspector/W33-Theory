#!/usr/bin/env python3
"""Pass 11555: exact local temporal dynamics algebra on the 4x9 foliation.

The eight directed null steps split under the actual PSp linear group W(D3)
into the two chi tetrahedra.  Therefore every translation- and PSp-invariant
nearest-null kernel is two-dimensional.  Its reversal-even basis is the null
adjacency A; its reversal-odd basis is the cubic-arrow operator D.

D is rank eight, supported spectrally only on the old null/Weil sector, and
satisfies an exact polynomial square identity with the graph Laplacian.
"""
from __future__ import annotations
import itertools, json, math
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/PART_W33_PASS11555_TEMPORAL_LOCAL_DYNAMICS_DIRAC.json"
P=3
V=list(itertools.product(range(3),repeat=3)); N=len(V)

def sub(a,b): return tuple((a[i]-b[i])%3 for i in range(3))
def chi(d):
    r=1
    for x in d: r*=1 if x==1 else -1
    return r
def mv(M,v): return tuple(sum(M[i][j]*v[j] for j in range(3))%3 for i in range(3))
def sign_product(M):
    r=1
    for row in M:
        x=next(x for x in row if x)
        r*=1 if x==1 else -1
    return r

def signed_perms():
    out=[]
    for p in itertools.permutations(range(3)):
        for s in itertools.product((1,2),repeat=3):
            M=[[0]*3 for _ in range(3)]
            for i,j in enumerate(p): M[i][j]=s[i]
            out.append(tuple(tuple(r) for r in M))
    return sorted(set(out))

def main():
    S=[d for d in V if all(d)]
    G=[M for M in signed_perms() if sign_product(M)==1]
    assert len(S)==8 and len(G)==24
    unseen=set(S); orbits=[]
    while unseen:
        d=min(unseen); orb={mv(M,d) for M in G}; orbits.append(sorted(orb)); unseen-=orb
    assert sorted(map(len,orbits))==[4,4]
    assert sorted({chi(d) for d in o}.pop() for o in orbits)==[-1,1]

    A=np.zeros((N,N),dtype=np.int64); D=np.zeros((N,N),dtype=np.int64)
    for i,x in enumerate(V):
        for j,y in enumerate(V):
            if i==j: continue
            d=sub(y,x)
            if all(d):
                A[i,j]=1; D[i,j]=chi(d)
    I=np.eye(N,dtype=np.int64); L=8*I-A
    assert np.array_equal(D.T,-D)
    assert np.array_equal(D@L,L@D)
    assert np.array_equal(3*(D@D),L@(L-6*I)@(L-12*I))
    rankD=int(np.linalg.matrix_rank(D.astype(float))); assert rankD==8

    # Fourier eigenvalues, grouped by Hamming weight and cubic chirality.
    omega=np.exp(2j*np.pi/3)
    rows=[]
    for k in V:
        w=sum(x!=0 for x in k)
        ahat=sum(omega**(sum(k[i]*d[i] for i in range(3))%3) for d in S)
        dhat=sum(chi(d)*omega**(sum(k[i]*d[i] for i in range(3))%3) for d in S)
        expectedA=(2**(3-w))*((-1)**w)
        assert abs(ahat-expectedA)<1e-10
        if w<3: assert abs(dhat)<1e-10
        else:
            target=-1j*3*math.sqrt(3)*chi(k)
            assert abs(dhat-target)<1e-10
        rows.append({"k":list(k),"weight":w,"A_eigen":int(expectedA),
          "L_eigen":int(8-expectedA),"chi":chi(k) if w==3 else None,
          "D_eigen_imag":float(round(dhat.imag,12))})

    p11549=json.loads((ROOT/"data/PART_W33_PASS11549_HISTORY_ARROW_CUBIC_CHARACTER.json").read_text())
    p11550=json.loads((ROOT/"data/PART_W33_PASS11550_VERONESE_WEIL_HAMMING_ORIENTATION.json").read_text())
    assert p11549["triangle_foliation"]["factorization"]=="36 = 4 * 9"
    assert p11550["closed_formula"]["Weil_phase_equals"]=="i * chi"

    out={
      "schema":"w33.pass11555.temporal_local_dynamics_dirac.v1",
      "status":"PASS_EXACT_TEMPORAL_LOCAL_OPERATOR_ALGEBRA","pass":11555,
      "local_kernel_classification":{"directed_null_steps":8,"PSp_orbits":[4,4],
        "orbit_labels":"chi=+1 and chi=-1 tetrahedra",
        "translation_and_PSp_invariant_nearest_null_kernel_dimension":2,
        "reversal_even_dimension":1,"reversal_even_basis":"A, the distance-3/null adjacency",
        "reversal_odd_dimension":1,"reversal_odd_basis":"D[x,y]=chi(y-x)"},
      "operator_D":{"skew_symmetric":True,"rank":8,"kernel_dimension":19,
        "commutes_with_L":True,
        "exact_square_identity":"3 D^2 = L(L-6I)(L-12I)",
        "weight3_projector":"P_3 = -D^2/27 = -L(L-6I)(L-12I)/81",
        "partial_complex_structure":"J=-D/(3 sqrt(3)); J^2=-P_3",
        "Weil_weld":"on the eight weight-3 modes, J eigenvalue = i*chi = the Pass11550 normalized Weil phase"},
      "fourier_sectors":{"L_spectrum_by_Hamming_weight":{"0":0,"1":12,"2":6,"3":9},
        "multiplicities":{"0":1,"1":6,"2":12,"3":8},
        "D_law":"Dhat(k)=0 for weight<3; Dhat(k)=-i*3*sqrt(3)*chi(k) for weight=3",
        "weight3_chiral_split":"4 + 4"},
      "general_local_Hermitian_operator":{
        "formula":"K = m^2 I + alpha L + beta iD",
        "eigenvalues":["m^2 (w=0)","m^2+12 alpha (w=1, mult 6)",
          "m^2+6 alpha (w=2, mult 12)",
          "m^2+9 alpha + 3 sqrt(3) beta chi (w=3, 4+4)"],
        "reading":"beta is the unique nearest-null PSp-invariant chiral splitting coefficient; it acts only on the old Weil/null sector"},
      "triangle_foliation":"Pass11549 supplies 36=4*9 edge-disjoint oriented temporal triangles, so the odd local operator D is assembled from that same orientation field.",
      "boundary":"Exact finite local-operator classification. A kinetic term, alpha, beta, mass, continuum scaling and physical interpretation are not derived; D is a finite chiral first-order operator, not asserted to be the observed Dirac operator."
    }
    OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":out["status"],"rankD":rankD,"orbits":[len(o) for o in orbits]},indent=2))
if __name__=="__main__": main()
