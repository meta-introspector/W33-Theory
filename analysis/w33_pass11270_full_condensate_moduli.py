#!/usr/bin/env python3
"""Numerical physical E6 quotient Hessian; canonical and compensator Kähler control.

The canonical spectrum is prior art, Goh et al. arXiv:2505.07931 table4.
The executable W33 signed-tensor quotient, not a new mass law, is the deliverable.
"""
from pathlib import Path
import sys,json
import numpy as np
from scipy.linalg import null_space,eigh
from scipy.optimize import minimize_scalar
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_20261001_degree18_phase_completion as D
import w33_20261001_completed_eft_one_loop as L

def reference():
    z=np.exp(2j*np.pi/9);p=D.I.C.plane()@np.array([1,z,z.conjugate()])/3
    acts=np.einsum('aij,j->ia',L.generators()[:78],p)
    N=null_space(acts.conj().T,rcond=1e-10)
    assert N.shape==(81,11)
    assert np.max(abs(acts.conj().T@N))<1e-12
    return p,N

def superpotential(p,scale=3/14):
    i,g=D.evaluate(p);w=scale*np.exp(-np.log(-729*i)/3)
    return w,-w*g/(3*i)

def potential(p,kappa=0.,scale=3/14):
    """K=s+kappa*s², compensator m=1, Lambda^9=scale. All81 fields."""
    s=float(np.vdot(p,p).real);a=1+2*kappa*s;d=1+4*kappa*s
    if min(a,d)<=0:return np.inf
    w,g=superpotential(p,scale);eg=p@g
    vf=float(np.vdot(g,g).real)/a-2*kappa*abs(eg)**2/(a*d)
    soft=2*a/d*eg.real-6*w.real+(s*a*a/d-s-kappa*s*s)
    return float(vf+soft)

def hessian(p,N,kappa,h):
    dirs=np.column_stack([N,1j*N])/np.sqrt(2);dim=dirs.shape[1]
    v0=potential(p,kappa);out=np.zeros((dim,dim));grad=np.zeros(dim)
    for i in range(dim):
        e=h*dirs[:,i];vp=potential(p+e,kappa);vm=potential(p-e,kappa)
        grad[i]=(vp-vm)/(2*h);out[i,i]=(vp+vm-2*v0)/h**2
        for j in range(i):
            f=h*dirs[:,j]
            out[i,j]=out[j,i]=(potential(p+e+f,kappa)+potential(p-e-f,kappa)-potential(p+e-f,kappa)-potential(p-e+f,kappa))/(4*h*h)
    s=float(np.vdot(p,p).real);K=(1+2*kappa*s)*np.eye(81)+2*kappa*np.outer(p,p.conj())
    metric=2*np.real(dirs.conj().T@K@dirs)
    return grad,out,eigh(out,metric,eigvals_only=True)

def payload():
    p,N=reference();w,g=superpotential(p)
    assert np.max(abs(g+6*w*p.conj()))<1e-12
    Hf=np.column_stack([(superpotential(p+1e-5*N[:,j])[1]-superpotential(p-1e-5*N[:,j])[1])/(2e-5) for j in range(11)])
    ferm=np.linalg.svd(N.T@Hf,compute_uv=False)**2
    expected_f=np.array([81/49]*8+[324/49]*2+[81.]);assert np.max(abs(np.sort(ferm)-expected_f))<1e-5
    rows=[]
    for k in (0.,.01,-.01):
        radial=minimize_scalar(lambda r:potential(r*p,k),bounds=(.8,1.2),method='bounded',options={'xatol':1e-12})
        r=float(radial.x);hs=[]
        for step in (.001,.0005):
            gr,H,eig=hessian(r*p,N,k,step);hs.append({'step':step,'gradient_norm':float(np.linalg.norm(gr)),'eigenvalues':eig.tolist(),'Hessian':H.tolist()})
        rich=(4*np.array(hs[1]['Hessian'])-np.array(hs[0]['Hessian']))/3
        dirs=np.column_stack([N,1j*N])/np.sqrt(2);q=r*p;s=float(np.vdot(q,q).real)
        K=(1+2*k*s)*np.eye(81)+2*k*np.outer(q,q.conj());metric=2*np.real(dirs.conj().T@K@dirs)
        eig=eigh(rich,metric,eigvals_only=True)
        assert np.max(abs(eig[:8]))<.002 and min(eig[8:])>2
        if k==0:
            target=np.array([0.]*8+[162/49]*10+[486/49]*2+[486/7,648/7])
            assert np.max(abs(eig-target))<.003
        rows.append({'kappa':k,'stationary_radius':r,'potential':potential(q,k),'steps':hs,'extrapolated_eigenvalues':eig.tolist(),'positive_modes':14,'Goldstone_modes':8,'step_Hessian_difference_norm':float(np.linalg.norm(np.array(hs[0]['Hessian'])-np.array(hs[1]['Hessian']),2))})
        print('completed Kähler',k,flush=True)
    return {'status':'PASS','result_scope':'PASS_NUMERICAL_ELEVEN_COMPLEX_MODULUS_CONTROL','E6_complex_orbit_rank':70,'complex_physical_moduli':11,'canonical_scalar_masses_squared':{'0':8,'162/49':10,'486/49':2,'486/7':1,'648/7':1},'canonical_fermion_masses_squared':{'81/49':8,'324/49':2,'81':1},'computed_fermion_masses_squared':np.sort(ferm).tolist(),'rows':rows,
      'Kahler_action':'K=s+kappa*s^2; V=(W_i+m K_i) K^{i jbar} (Wbar_j+m Kbar_j)-m^2 K-3m(W+Wbar), m=1. Includes the compensator-induced correction, not only inverse-metric substitution.',
      'quotient_method':'N spans the complex orthogonal complement of the actual70-dimensional E6 orbit. Linear tangents are first-order D-flat; at ambient stationarity their second-order correction does not affect the Hessian. Canonical mass coordinates include the factor sqrt2; noncanonical masses solve the generalized metric eigenproblem.',
      'boundary':['Floating finite-step local test, not an interval stability proof or global minimum.','Only the stated radial Kähler family and two small coefficients are tested; general flavor-dependent Kähler corrections remain open.','The8 massless flavor Goldstone directions are expected for global SU3, and are not eight instabilities.','This E6-gauged/global-family near-SUSY theory is distinct from the paired EFT and still not an SM vacuum.'],
      'prior_owners':['analysis/w33_20261001_condensate_cartan_potential.py','analysis/w33_20261001_degree18_phase_completion.py'],'primary_source':'https://arxiv.org/html/2505.07931v1'}
if __name__=='__main__':
    out=payload();(ROOT/'data/w33_pass11270_full_condensate_moduli.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
