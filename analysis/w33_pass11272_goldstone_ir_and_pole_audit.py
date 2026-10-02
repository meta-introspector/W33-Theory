#!/usr/bin/env python3
"""Actual soft scalar slopes plus resummation kernel and pole non-identifiability."""
from pathlib import Path
import sys,json
import numpy as np
from scipy.linalg import eigh
import sympy as s
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_20261001_completed_eft_one_loop as L

def scalar_matrix(q,step=1e-5):
    p=L.I.C.plane()@q/np.sqrt(3);a,b,g6,g12,H6,H12=L.field_data(p,step)
    pars=L.PARAM;acts=np.einsum('aij,j->ai',L.generators(),p)
    gs=np.r_[np.full(78,pars['g_trace']),np.full(8,pars['g_trace']/np.sqrt(6))]
    j=np.sqrt(2)*np.concatenate([acts.real,acts.imag],axis=1)*gs[:,None]
    M=j.T@j;x=np.sqrt(2)*np.r_[p.real,p.imag];norm=np.vdot(p,p).real
    M+=2*pars['lambda']*(np.outer(x,x)+(norm-1)*np.eye(162))
    for val,g,H,alpha in [(a,g6,H6,pars['alpha6']),(b,g12,H12,pars['alpha12'])]:
        J=np.r_[g,1j*g]/np.sqrt(2);hr=np.block([[H,1j*H],[1j*H,-H]])/2
        M+=2*alpha*np.real(np.outer(J.conj(),J)+val.conjugate()*hr)
    return M

def payload():
    q,dirs=L.tangent();M=scalar_matrix(q);ev,U=eigh(M);soft=U[:,abs(ev)<1e-9];assert soft.shape[1]==79
    h=.0005;qp=q+h*dirs[:,0];qm=q-h*dirs[:,0];qp/=np.linalg.norm(qp);qm/=np.linalg.norm(qm)
    dM=(scalar_matrix(qp)-scalar_matrix(qm))/(2*h);slopes=eigh(soft.T@dM@soft,eigvals_only=True)
    ircoef=float(np.sum(slopes**2)/(32*np.pi**2));assert ircoef>0
    x,mu,delta=s.symbols('G mu Delta',positive=True)
    f=x*x*(s.log(x/mu**2)-s.Rational(3,2))/(64*s.pi**2)
    assert s.simplify(s.diff(f,x,2)-s.log(x/mu**2)/(32*s.pi**2))==0
    # f(G+Delta) is a resummed soft kernel. Its expansion is explicitly recorded;
    # double-counting subtractions depend on the loop order of the hard Delta.
    second=s.simplify(s.diff(f.subs(x,x+delta),x,2).subs(x,0))
    rows=[{'hard_shift_input':d,'soft_curvature_at_G0':float(np.sum(slopes**2)*np.log(d)/(32*np.pi**2))} for d in (1e-2,1e-3,1e-4)]
    # Two inverse propagators with identical static curvature but different poles.
    # Gamma(p²)=Z p²-M²: the missing kinetic self-energy is observable at p²!=0.
    return {'status':'PASS','result_scope':'PASS_ACTUAL_SOFT_IR_COEFFICIENT_AND_POLE_INPUT_AUDIT','soft_scalar_dimension':79,'soft_direction':0,'projected_linear_soft_slopes':slopes.tolist(),'log_abs_angular_displacement_curvature_coefficient':ircoef,
      'IR_derivation':'For G_j(x)=a_j*x+O(x²), sum f(G_j) has d²/dx² = (sum a_j²)/(32pi²)*log|x|+finite. This positive coefficient times log|x| tends negative infinity; finite-step negative Hessians are not physical tachyon masses.',
      'resummed_kernel':'f(G+Delta_hard); expand and subtract terms already present at the declared loop order to avoid double counting. Delta_hard must be computed from hard1PI self-energy, not picked to make curvature positive.',
      'resummed_second_G_derivative_at_zero':str(second),'input_sensitivity_rows':rows,
      'pole_equation':'det[p² Z-M_tree²-Pi(p²)]=0 at the shifted stationary vacuum; Pi(0) alone is insufficient.',
      'same_static_curvature_counterexample':{'M_squared':1,'Z1':1,'pole1':1,'Z2':2,'pole2':'1/2'},
      'physical_poles_computed':False,'missing':['hardGoldstone1PI self-energy for the named paired EFT','BRST/ghost/gauge-fixing and counterterm completion','momentum-dependent physical scalar/vector/fermion self-energies and shifted stationarity'],
      'boundary':'The kernel and actual IR slopes are computed. This is not a full resummed effective action or a physical pole-mass prediction; unknown self-energies are explicit rather than silently replaced by a cutoff.',
      'prior_owners':['analysis/w33_20261001_completed_eft_one_loop.py'],'primary_source':'https://arxiv.org/abs/1406.2355'}
if __name__=='__main__':
    out=payload();(ROOT/'data/w33_pass11272_goldstone_ir_and_pole_audit.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
