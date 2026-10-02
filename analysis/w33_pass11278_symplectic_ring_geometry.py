#!/usr/bin/env python3
"""Literal mixed-characteristic projective lift; no spacetime metric inference."""
from pathlib import Path
import json,itertools
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
J=np.array([[0,0,1,0],[0,0,0,1],[-1,0,0,0],[0,-1,0,0]],int)

def points(n):
    pts=set()
    for v in itertools.product(range(n),repeat=4):
        ix=next((i for i,x in enumerate(v) if x%3),None)
        if ix is None:continue
        u=pow(v[ix],-1,n);pts.add(tuple(u*x%n for x in v))
    return np.array(sorted(pts),int)

def payload():
    p=points(9);assert len(p)==1080
    small=points(3);assert len(small)==40
    def norm(v):
        i=next(i for i,x in enumerate(v) if x%3);return tuple(int(x*pow(int(v[i]),-1,3)%3) for x in v)
    fibers={norm(v) for v in p};assert len(fibers)==40
    count={v:sum(norm(w)==v for w in p) for v in fibers};assert set(count.values())=={27}
    ortho=(p@J@p.T)%9==0;np.fill_diagonal(ortho,False);assert set(ortho.sum(1))=={116}
    # Congruence kernel I+3X modulo9: ten free entries because JX is symmetric.
    basis=[]
    for i in range(4):
        for j in range(i,4):
            S=np.zeros((4,4),int);S[i,j]=S[j,i]=1;basis.append(-J@S)
    assert len(basis)==10
    for coeff in itertools.product(range(3),repeat=10):
        X=sum((c*b for c,b in zip(coeff,basis)),np.zeros((4,4),int));g=np.eye(4,dtype=int)+3*X
        assert np.all((g.T@J@g-J)%9==0)
    # Integral symplectic shear is unipotent and cannot preserve positive metric:
    # T e2=e2+e0. Invariance on these vectors forces G00=0.
    T=np.eye(4,dtype=int);T[0,2]=1;assert np.array_equal(T.T@J@T,J)
    return {'status':'PASS','result_scope':'PASS_EXACT_SYMPLECTIC_RING_LIFT_AND_METRIC_OBSTRUCTION',
      'projective_count_formula':'|Pprimitive((Z/3^k)^4)|=40*3^(3(k-1)); reduction fibers have27 points per level.',
      'level2_points':len(p),'base_points':len(small),'level2_fiber_size':27,'level2_orthogonality_graph_degree':116,
      'level2_kernel_exhausted':59049,'symplectic_order_formula':'|Sp4(Z/3^k)|=51840*3^(10(k-1)); smooth symplectic group lifts, with sp4(F3) additive kernel at each successive reduction.',
      'splitting_boundary':'No claim that the mixed-characteristic extension splits. Pass4937 proves the split dual-number extension; equal order does not identify it with Z/9.',
      'limit':'The inverse limit is a3-adic symplectic module Z3^4. An Archimedean four-torus requires a separate metric, embedding and scaling choice.',
      'metric_obstruction':'An integral symplectic shear sends e2 to e2+e0 and fixes e0. From T^TGT=G the (0,2) entry forces G00=0. Thus no positive definite real quadratic metric is invariant under all integral symplectic lifts. Symplectic data alone cannot select the previous Euclidean Wilson operator.',
      'prior_owners':['analysis/w33_symplectic_basis_regular_lift.py','analysis/PASS4937_ADJOINT_DUAL_NUMBER_CONTROLLER.md','analysis/PASS444_HJELMSLEV_CONDUCTOR_GEOMETRY.md','analysis/w33_pass11273_wilson_uniform_tail.py'],
      'boundary':['Classical smooth-group and projective-ring counts are not new theorems. This implements their precise W33 bridge and metric obstruction.','No Lorentz signature, continuum Einstein action or physical spacetime dimension follows from this lift.']}
if __name__=='__main__':
    out=payload();(ROOT/'data/w33_pass11278_symplectic_ring_geometry.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
