#!/usr/bin/env python3
"""The literal W33 clique history does not carry a topological4-form flux."""
from pathlib import Path
import sys,json,itertools
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
from w33_pass11278_symplectic_ring_geometry import points,J
from w33_pass11275_polynomial_sm_higgs import rank_mod

def complex_data():
    P=points(3);adj=(P@J@P.T)%3==0;np.fill_diagonal(adj,False)
    edges=[(i,j) for i in range(40) for j in range(i+1,40) if adj[i,j]]
    tris=[(i,j,k) for i,j in edges for k in range(j+1,40) if adj[i,k] and adj[j,k]]
    tets=[t+(l,) for t in tris for l in range(t[-1]+1,40) if all(adj[i,l] for i in t)]
    assert len(edges)==240 and len(tris)==160 and len(tets)==40
    ti={t:i for i,t in enumerate(tris)};D=np.zeros((160,40),int)
    for j,t in enumerate(tets):
        for k in range(4):D[ti[t[:k]+t[k+1:]],j]=(-1)**k
    assert np.all(np.sum(abs(D),axis=1)==1)
    return edges,tris,tets,D

def payload():
    edges,tris,tets,D=complex_data();assert rank_mod(D)==40
    # Product with3-edge periodic history circle. D4 upper block=D3 tensor I3.
    d1=np.array([[-1,0,1],[1,-1,0],[0,1,-1]],int)
    D4=np.vstack([np.kron(D,np.eye(3,dtype=int)),-np.kron(np.eye(40,dtype=int),d1)])
    assert D4.shape==(600,120) and rank_mod(D4)==120
    volumes=np.ones(120,dtype=int);assert np.any(D4@volumes)
    return {'status':'PASS','result_scope':'PASS_EXACT_OBSTRUCTION_TO_NAIVE_W33_HISTORY_FLUX_SEQUESTERING',
      'explicit_complex':'K=clique complex of W(3,3), history H=K times a3-edge oriented circle.',
      'W33_cells':[40,240,160,40],'history_4_cells':120,'history_boundary_4_shape':list(D4.shape),'history_boundary_4_rank_mod101':120,
      'exact_top_boundary_injective':'Each triangular face belongs to exactly one tetrahedron, giving disjoint column supports and injective D3 over Z. The D3 tensor I3 block makes D4 injective over Z.',
      'ownership':'The clique f-vector, injective top boundary and bouquet-of81-circles homotopy are existing repo results, not new. This applies them to the specifically named history/4-form action.',
      'topological_obstruction':'H4(H;Z)=0 and H^4(H;R)=0. Every real4-cochain is exact; there is no intrinsic nonzero closed fundamental4-cycle or topological4-flux on this literal history complex.',
      'named_discrete_action':'S_top=<sigma(Lambda), F>, F=D4^T A3+F_background. Free variation of every3-cochain A3 requires D4 sigma(Lambda)=0.',
      'variational_obstruction':'Since D4 is injective, sigma(Lambda)=0 cellwise; the gauge equation does not produce a nonzero rigid global multiplier. The all-ones top-volume chain has nonzero boundary. Imposing the usual positive volume constraint therefore needs boundary data or a different closed geometry.',
      'completion_options':'Relative cohomology with fixed boundary3-form, a closed oriented4-complex, or an added background flux sector. Each is new physical input, not supplied by finite W33 incidence.',
      'conditional_covariant_continuum':'The known local sequestering action adds sigma(Lambda/mu4) F4 and hatsigma(kappa²/M²) Fhat4 to integral sqrt(-g)[kappa²R/2-Lambda-Lm]. On a suitable continuum with flux sectors it yields kappa²G=T-(1/4)g<TrT>-DeltaLambda*g; DeltaLambda is a flux-dependent integration constant. W33 has not selected it.',
      'prior_owners':['analysis/w33_pass11274_scale_and_sequestering.py','analysis/w33_pass1448_1454_hodge_maxwell_and_the_missing_star.md','analysis/w33_pass1944_1948_the_flux_reading_fails_its_own_test.md','analysis/W33_LEDGER_CONTINUATION_2026_09_21.md'],
      'primary_sources':['https://arxiv.org/html/1505.01492v2'],
      'boundary':['This rejects a specific literal clique-complex history construction, not every possible W33-derived spacetime.','No cosmological constant value, Newton scale, covariant spacetime dynamics or universal no-go theorem is claimed.']}
if __name__=='__main__':
    out=payload();(ROOT/'data/w33_pass11279_history_flux_obstruction.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
