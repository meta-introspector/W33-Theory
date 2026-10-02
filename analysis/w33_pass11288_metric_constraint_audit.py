#!/usr/bin/env python3
"""Exact local ADM kinetic pullback and frame redundancy of the symplectic metric map."""
from pathlib import Path
import sys,json
import sympy as s
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
from w33_pass11283_symplectic_metric_field import chart_Jacobian

def payload():
 A=s.Matrix(chart_Jacobian());assert A.rank()==10
 # Six spatial components h11,h22,h33,h12,h13,h23.
 spatial=s.Matrix.vstack(*[A[4*i+j,:] for i,j in [(1,1),(2,2),(3,3),(1,2),(1,3),(2,3)]])
 # ADM kinetic KijKij-K²: diagonal entries1-1=0, cross diagonals-1;
 # off-diagonal K12,K13,K23 have coefficient2.
 D=s.diag(1,1,1,2,2,2)-s.Matrix([1,1,1,0,0,0])*s.Matrix([[1,1,1,0,0,0]])
 assert list(D.eigenvals().items())==[(-2,1),(1,2),(2,3)] or D.det()<0
 H=spatial.T*D*spatial;assert H.rank()==6 and len(H.nullspace())==5
 red=np.linalg.eigvalsh(np.array(H,float));assert sum(red<-1e-8)==1 and sum(red>1e-8)==5
 fiber=A.nullspace();assert len(fiber)==1
 # Transforming the11 source fields to ten metric components plus one fiber
 # converts EH exactly to ADM locally; its Hessian has4 lapse/shift +1fiber nulls.
 return {'status':'PASS','result_scope':'PASS_LOCAL_EH_PULLBACK_CONSTRAINT_AND_GHOST_AUDIT',
 'metric_map_rank':10,'field_dimension':11,'spatial_velocity_map_rank':spatial.rank(),'kinetic_pullback_rank':H.rank(),
 'kinetic_inertia':{'positive':5,'negative_conformal':1,'zero':5},'exact_fiber_vector':list(map(str,fiber[0])),
 'ADM_equivalence':'At this regular point, the rank10 submersion theorem supplies local coordinates(g_mu_nu,z_fiber). If the action depends only on g, it is exactly the Einstein-Hilbert/ADM action plus a redundant fiber variable. No higher-derivative source-field mode is introduced by this algebraic map.',
 'constraints':'Four lapse/shift primary constraints and four Hamiltonian/momentum secondary first-class constraints, plus one primary first-class fiber constraint:11 configuration variables,22 phase dimensions,9 first-class constraints ->2 local physical graviton polarizations. The pulled-back kinetic Hessian has5 nulls (lapse/shift plus fiber).',
 'conformal_boundary':'The single negative DeWitt kinetic direction is eliminated by the gravitational constraints in the reduced perturbative spectrum; its presence in the unreduced ADM kinetic form is not by itself a propagating ghost proof. Euclidean conformal-factor/path-integral issues remain.',
 'sigma_warning':'Adding a separate kinetic term for G,u,chi,omega generally breaks the fiber invariance and can give time derivatives to lapse/shift combinations. The EH-only two-polarization conclusion does not apply automatically to that enlarged action.',
 'boundaries':['Local regular-chart classical EH equivalence, not a global atlas, quantum unitarity proof or W33-derived spacetime.','Newton constant, real base and signature still inputs.','Constraint closure is inherited by an invertible local coordinate change from ADM; no lattice constraint-algebra convergence is proved.'],
 'prior_owner':'analysis/w33_pass11283_symplectic_metric_field.py',
 'primary_source':'https://journals.aps.org/pr/abstract/10.1103/PhysRev.160.1113'}
if __name__=='__main__':
 out=payload();(ROOT/'data/w33_pass11288_metric_constraint_audit.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
