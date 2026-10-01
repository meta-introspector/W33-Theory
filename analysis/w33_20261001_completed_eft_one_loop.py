#!/usr/bin/env python3
"""Full declared constant-background one-loop determinant with the three EFT lifts.

162-real-field scalar Hessian,86 vectors, two original81-Weyl species and
conjugate pair. Numerical Landau/MSbar calculation on the exact D-flat plane.
This is a specified cutoff-EFT test, not a UV-complete quantum stability theorem.
"""
from pathlib import Path
from functools import lru_cache
import sys,json,itertools
import numpy as np
from scipy.linalg import eigh
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_20261001_global_e6_cartan_covariants as I
import w33_20261001_complete_gauge_scalar_fermion_loop as L

@lru_cache(None)
def generators():
    H,F,_,_=L.rational_hermitian_basis()
    def orth(b):
        b=np.array([np.array(x,complex) for x in b]);gram=np.einsum('aij,bji->ab',b,b).real
        w,v=np.linalg.eigh(gram);return np.einsum('ab,bij->aij',(v/np.sqrt(w)).T,b)
    he,hf=orth(H),orth(F)
    return np.array([np.kron(h,np.eye(3)) for h in he]+[np.kron(np.eye(27),f.T) for f in hf])

PARAM={'g_trace':.2,'y':.15,'lambda':.1,'alpha6':.02,'alpha12':.0005,'c5':.001,'c6':.002/16,'c12':.003/100,'cutoff':1.,'mu':1.}

def field_data(phi,step=1e-5):
    a,b,g6,g12=I.evaluate(phi);n=len(phi);hs=[]
    for j in range(n):
        e=np.zeros(n);e[j]=step;plus=I.evaluate(phi+e);minus=I.evaluate(phi-e)
        hs.append([(plus[k]-minus[k])/(2*step) for k in (2,3)])
    H6,H12=np.array(hs).transpose(1,2,0)
    H6=(H6+H6.T)/2;H12=(H12+H12.T)/2
    return a,b,g6,g12,H6,H12

def spectra(q,step=1e-5):
    E=I.C.plane()/np.sqrt(3);phi=E@q;a,b,g6,g12,H6,H12=field_data(phi,step)
    pars=PARAM;acts=np.einsum('aij,j->ai',generators(),phi)
    gs=np.r_[np.full(78,pars['g_trace']),np.full(8,pars['g_trace']/np.sqrt(6))]
    gauge=2*np.real((gs[:,None]*acts).conj()@(gs[:,None]*acts).T)
    vec=eigh(gauge,eigvals_only=True)
    j=np.sqrt(2)*np.concatenate([acts.real,acts.imag],axis=1)*gs[:,None]
    scalar=j.T@j;x=np.sqrt(2)*np.r_[phi.real,phi.imag];norm=float(np.vdot(phi,phi).real)
    scalar+=2*pars['lambda']*(np.outer(x,x)+(norm-1)*np.eye(162))
    for val,g,H,alpha in [(a,g6,H6,pars['alpha6']),(b,g12,H12,pars['alpha12'])]:
        J=np.r_[g,1j*g]/np.sqrt(2);hr=np.block([[H,1j*H],[1j*H,-H]])/2
        scalar+=2*alpha*np.real(np.outer(J.conj(),J)+val.conjugate()*hr)
    scal=eigh(scalar,eigvals_only=True)
    _,op=I.C.operators();mass=pars['y']*op(phi)+pars['c5']*np.outer(phi.conj(),phi.conj())+pars['c6']*np.outer(g6,g6)+pars['c12']*np.outer(g12,g12)
    ferm=np.linalg.svd(mass,compute_uv=False)**2
    # Scalar Hessian includes all81-field directions; finite-difference control below.
    return vec,scal,ferm

def determinant(q,step=1e-5,diagnostic=None):
    v,s,f=spectra(q,step)
    if diagnostic is not None:diagnostic.append({'negative_scalar_modes':int(sum(s < -1e-10)),'minimum_scalar_mass_squared':float(min(s))})
    def term(m,c):
        z=np.abs(m);nz=z>1e-13
        return float(np.sum(m[nz]**2*(np.log(z[nz]/PARAM['mu']**2)-c)))
    return (3*term(v,5/6)+term(s,1.5)-8*term(f,1.5))/(64*np.pi**2)

def tangent():
    z=np.exp(2j*np.pi/9);q=np.array([1,z,z.conjugate()])/np.sqrt(3)
    gs,_=I.C.gradient_data();basis=np.array([g.conj()/np.linalg.norm(g) for g in gs]).T
    assert np.max(abs(basis.conj().T@basis-np.eye(2)))<1e-12
    return q,np.column_stack([basis[:,0],1j*basis[:,0],basis[:,1],1j*basis[:,1]])/np.sqrt(2)

def local(h):
    q,D=tangent();value=determinant(q);cache={};diagnostics=[]
    def v(a):
        key=tuple(a)
        if key not in cache:
            p=q+D@np.array(a);p/=np.linalg.norm(p);cache[key]=determinant(p,diagnostic=diagnostics)
        return cache[key]
    H=np.zeros((4,4));grad=[]
    for i in range(4):
        a=np.zeros(4);a[i]=h
        H[i,i]=(v(a)+v(-a)-2*value)/h**2;grad.append((v(a)-v(-a))/(2*h))
        for j in range(i):
            b=np.zeros(4);b[j]=h
            H[i,j]=H[j,i]=(v(a+b)+v(-a-b)-v(a-b)-v(-a+b))/(4*h*h)
    return value,np.array(grad),H,diagnostics

def payload():
    q,_=tangent();v,sc,f=spectra(q);v2,sc2,f2=spectra(q,step=5e-6)
    serr=float(np.max(abs(sc-sc2)));assert serr<2e-7
    assert sum(v>1e-10)==78 and sum(sc>1e-9)==83 and sum(f>1e-10)==81
    assert np.max(abs(np.sqrt(np.sort(f)[:3])-np.array([.001,.002,.003])))<1e-10
    # 78 gauge-orbit Goldstones+one common phase; 78D+5 radial/angular scalar modes.
    rows=[]
    for h in (.001,.0005):
        val,grad,H,diag=local(h);rows.append({'step':h,'V1_at_T':val,'angular_gradient':grad.tolist(),'angular_Hessian':H.tolist(),'off_shell_scalar_diagnostics':diag});print('completed angular step',h,flush=True)
    H1,H2=[np.array(r['angular_Hessian']) for r in rows];rich=(4*H2-H1)/3
    tree=np.diag([16*PARAM['alpha6']]*2+[100*PARAM['alpha12']]*2)
    eig=np.linalg.eigvalsh(tree+rich)
    return {'status':'PASS_DECLARED_COMPLETED_EFT_FULL_ONE_LOOP_LOCAL_CONTROL','parameters':PARAM,
      'field_determinant':'3*Tr mV^4(log(mV^2/mu^2)-5/6)+Tr mS^4(log|mS^2/mu^2|-3/2)-8*Tr mF^4(log(mF^2/mu^2)-3/2), divided64pi^2',
      'counts_at_T':{'massive_vectors':78,'positive_real_scalars':83,'massive_singular_fermion_pairs':81,'remaining_scalar_zeros':79},
      'scalar_Hessian_step_replay_error':serr,'local_rows':rows,'Richardson_loop_Hessian':rich.tolist(),'tree_angular_Hessian':tree.tolist(),
      'finite_step_tree_plus_one_loop_angular_eigenvalues':eig.tolist(),'finite_step_positive_curvature':bool(min(eig)>0),'Hessian_step_difference_norm':float(np.linalg.norm(H2-H1,2)),
      'infrared_boundary':'Goldstone masses vanish at T and become negative linearly along nearby angular directions. The real unresummed determinant is infrared sensitive; the finite-step positive curvature is not a twice-differentiable vacuum-stability certificate. A Goldstone-resummed or physical quotient calculation is needed.',
      'scalar_method':'All162 canonical real fields, inverse-trace-Gram gauge generators; analytic potential Hessian chain rule using central differences of analytic holomorphic gradients. Backgrounds remain on the normalized exact D-flat plane.',
      'boundary':['Specified coefficients/cutoff/renormalization condition; numerical local curvature, not directed proof or global quantum minimum.','A positive curvature at T is not stationarity: the computed loop gradient must be included when solving for the shifted vacuum.','Off-shell negative scalar eigenvalues use log absolute value, i.e. real part of the determinant; imaginary instability contributions are not discarded as a stability proof.','The common phase is flat here; the new I18 selector is not added to this determinant.','All higher-dimensional counterterms consistent with the symmetry are available. This fixed EFT loop test does not determine their finite matching coefficients or establish UV completion.'],
      'prior_owners':['analysis/w33_20261001_global_e6_cartan_covariants.py','analysis/w33_20261001_complete_gauge_scalar_fermion_loop.py','analysis/w33_20261001_exact_cartan_mirror_completion.py'],
      'primary_sources':['https://arxiv.org/abs/hep-ph/0111209','https://arxiv.org/abs/1406.2355']}
if __name__=='__main__':
    out=payload();(ROOT/'data/w33_20261001_completed_eft_one_loop.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
