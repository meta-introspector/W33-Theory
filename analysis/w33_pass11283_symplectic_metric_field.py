#!/usr/bin/env python3
"""An explicit local Lorentz-metric chart over the W33 symplectic form.

The Archimedean base, field action, orientation and coefficients are additional
inputs. This does not identify a finite-field symmetry with a real Lie group.
"""
from pathlib import Path
import sys,json
import sympy as s
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
from w33_pass11278_symplectic_ring_geometry import J
from w33_pass11275_polynomial_sm_higgs import rank_mod

def metric(G,u,chi=0.,omega=0.):
    assert abs(u@G@u-1)<1e-9
    v=-J@G@u;wu=G@u;wv=G@v
    L=np.exp(-chi)*(-np.outer(wu,wu)+np.outer(wv,wv))+np.exp(chi)*(G-np.outer(wu,wu)-np.outer(wv,wv))
    return np.exp(2*omega)*L

def tangents():
    out=[]
    for typ in range(2):
        for i,j in [(0,0),(0,1),(1,1)]:
            A=np.zeros((2,2),int);A[i,j]=A[j,i]=1
            out.append(np.block([[A,np.zeros((2,2),int)],[np.zeros((2,2),int),-A]]) if typ==0 else np.block([[np.zeros((2,2),int),A],[A,np.zeros((2,2),int)]]))
    return out

def chart_Jacobian():
    u=np.eye(4,dtype=int)[:,0];cols=[]
    for D in tangents():
        du=-s.Rational(int(D[0,0]),2)*s.Matrix(u);dv=-s.Matrix(J)*(s.Matrix(D)*s.Matrix(u)+du)
        dw=s.Matrix(D)*s.Matrix(u)+du
        # chi=0: L=G-2(Gu)(Gu)^T, independent of the companion v.
        L=s.Matrix(D)-2*(dw*s.Matrix(u).T+s.Matrix(u)*dw.T);cols.append(list(L))
    for i in (1,2,3):
        du=s.eye(4)[:,i];L=-2*(du*s.Matrix(u).T+s.Matrix(u)*du.T);cols.append(list(L))
    cols.append(list(s.diag(1,1,-1,1))) # chi derivative: -P_uv+(I-P_uv).
    cols.append(list(2*s.diag(-1,1,1,1))) # conformal omega derivative.
    a=np.array(s.Matrix.hstack(*[s.Matrix(c) for c in cols]),int);return a

def payload():
    a=chart_Jacobian();assert rank_mod(a[:,:9])==8 and rank_mod(a[:,:10])==9 and rank_mod(a)==10
    Ds=tangents();gram=np.array([[np.trace(x@y) for y in Ds] for x in Ds],int);assert min(np.linalg.eigvalsh(gram))>0
    return {'status':'PASS','result_scope':'PASS_EXPLICIT_SYMPLECTIC_TO_LORENTZ_METRIC_FIELD_CHART',
      'positive_metric_field':'G symmetric positive, G J G=J; it ranges over Sp4(R)/U2 with6 local real coordinates. This introduces an Archimedean metric field, rather than selecting a fixed metric invariant under all symplectic shears.',
      'orientation_field':'u^T G u=1, v=-J G u. The two vectors are G-orthonormal. Let w_u=Gu,w_v=Gv, P=w_u w_u^T+w_v w_v^T.',
      'explicit_metric':'g=e^(2omega){e^(-chi)[-w_u w_u^T+w_v w_v^T]+e^chi[G-P]}. For every realchi,omega it has signature(-,+,+,+) and detg=-e^(8omega).',
      'field_variables':'6 metric coordinates +3 unit-vector coordinates +chi +omega =11 fields, with a1-dimensional local redundancy at the reference. The exact map has rank9 withoutomega and rank10 includingomega.',
      'Jacobian_at_reference':a.tolist(),'exact_full_rank':10,'exact_unimodular_rank':9,
      'why_chi_is_needed':'At chi=omega=0 the reflected compatible metric G-2(Gu)(Gu)^T has only rank8 as a map to determinant-1 Lorentz metrics. chi supplies the missing relative symplectic-plane scale; omega supplies the volume degree of freedom.',
      'sigma_model_metric_Gram':gram.tolist(),'declared_gravity_action':'S=integral d4x sqrt(-detg)[M² R[g]/2-Lambda0+Lm(g,Phi)]. g is the explicit field map above. Its rank10 means metric variations are locally unrestricted; with appropriate boundaries its Euler-Lagrange equations are the usual Einstein equations in this declared continuum theory.',
      'declared_positive_target_kinetic':'On a supplied Euclidean base, f²/8 Tr(G^-1 partialG G^-1 partialG) gives a positive6-dimensional target kinetic metric. This does not by itself prove a ghost-free Lorentzian gravity-plus-sigma model; the mapped Einstein action and its gauge constraints must be treated separately.',
      'covariance':'For symplectic S, G->S^-T G S^-1 and u->S u give g->S^-T g S^-1. A choice of compatibleG and u breaks the fixed symplectic symmetry; no invariant positive metric was assumed.',
      'boundaries':['This is an explicit local metric architecture on an additionally supplied real4-manifold. No W33-derived base, signature selection, Newton scale, cutoff or CC is proved.','The finite Sp4(F3) is not embedded into Sp4(R) by reading its entries as real numbers; reduction is a quotient operation, not a characteristic-changing embedding.','Siegel-space and Lie-algebra identifications are classical prior art. This witness names and checks the extra fields needed for the actual metric map.'],
      'prior_owners':['analysis/w33_pass11278_symplectic_ring_geometry.py','analysis/w33_BREAKTHROUGH_AdS4_Siegel.py','analysis/w33_BREAKTHROUGH_366_spacetime_emergence_Minkowski.py','analysis/w33_einstein_field_equations_from_spectral_action.py'],
      'primary_sources':['https://arxiv.org/abs/1504.03963']}
if __name__=='__main__':
    out=payload();(ROOT/'data/w33_pass11283_symplectic_metric_field.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
