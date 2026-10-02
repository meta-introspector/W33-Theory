#!/usr/bin/env python3
"""All22 scalar/11 Weyl radial cuts and subtracted momentum dispersion.

Uses the canonical fixed-slice mass profiles, not a full nonlinear quotient action.
Other external self-energy matrix entries and heavy gauge thresholds remain open.
"""
from pathlib import Path
import sys,json
import numpy as np
from scipy.integrate import quad
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11272_resummed_condensate_radial_control as R

def channels(m=.01):
    F,A=R.matrices();evals,U=np.linalg.eigh(m*m*(F+A));g=U.T@(m*m*(-16*F-8*A)/np.sqrt(2))@U
    scalar=[]
    for i in range(6):
        for j in range(i,6):
            if abs(g[i,j])>1e-10:scalar.append({'kind':'scalar','m1_squared':float(evals[i]),'m2_squared':float(evals[j]),'weight':float(g[i,j]**2/(32*np.pi)*(1 if i==j else 2)),'multiplicity':1})
    partner=162/49*m*m;gp=(-16*729/49+8*81/7)*m*m/np.sqrt(2)
    scalar.append({'kind':'scalar','m1_squared':partner,'m2_squared':partner,'weight':float(8*gp*gp/(32*np.pi)),'multiplicity':8})
    gg=648/7*m*m/np.sqrt(2)
    scalar.append({'kind':'scalar','m1_squared':0.,'m2_squared':0.,'weight':float(8*gg*gg/(32*np.pi)),'multiplicity':8})
    ferm=[]
    for mass2,n in [(81/49*m*m,8),(324/49*m*m,2),(81*m*m,1)]:
        ferm.append({'kind':'Majorana','m1_squared':mass2,'m2_squared':mass2,'weight':float(n*32*mass2/(16*np.pi)),'multiplicity':n})
    return scalar+ferm

def threshold(c):return (np.sqrt(c['m1_squared'])+np.sqrt(c['m2_squared']))**2

def rho(c,t):
    if t<=threshold(c):return 0.
    if c['kind']=='Majorana':return c['weight']*t*(1-4*c['m1_squared']/t)**1.5
    b=1-2*(c['m1_squared']+c['m2_squared'])/t+(c['m1_squared']-c['m2_squared'])**2/t**2
    return c['weight']*np.sqrt(max(0,b))

def dispersion(c,s,s0,upper_factor=100.):
    th=threshold(c);upper=max(upper_factor*s,2*th+abs(s0))
    f=lambda t:rho(c,t)/(t-s0)**2
    if th<s<upper:
        main=quad(f,th,upper,weight='cauchy',wvar=s,epsabs=1e-11,limit=300)[0]
    else:main=quad(lambda t:f(t)/(t-s),th,upper,epsabs=1e-11,limit=300)[0]
    tail=quad(lambda t:f(t)/(t-s),upper,np.inf,epsabs=1e-11,limit=300)[0]
    return (s-s0)**2*(main+tail)/np.pi

def payload(m=.01):
    cs=channels(m);M2=648/7*m*m;M=np.sqrt(M2);s0=-m*m
    real=sum(dispersion(c,M2,s0) for c in cs);again=sum(dispersion(c,M2,s0,200.) for c in cs)
    assert abs(real-again)<2e-9
    scalar_im=sum(rho(c,M2) for c in cs if c['kind']=='scalar');ferm_im=sum(rho(c,M2) for c in cs if c['kind']=='Majorana');im=scalar_im+ferm_im
    gold=next(c for c in cs if c['kind']=='scalar' and c['m1_squared']==c['m2_squared']==0)
    assert abs(rho(gold,M2)/M-M**3/(8*np.pi))<1e-14
    assert im>0
    pole=complex(M2-real,-im)
    return {'status':'PASS','result_scope':'PASS_ALL_INTERNAL_MODULUS_RADIAL_CUTS_AND_SUBTRACTED_SELF_ENERGY',
      'internal_modes':{'real_scalars':22,'Weyl_Majorana_fermions':11,'massless_global_Goldstones':8},
      'external_sector':'The radial self-energy only. This includes every internal scalar and Weyl mode from the declared canonical fixed-slice condensate profiles; it does not calculate the full22x22 external self-energy matrix.',
      'vertices':'g_hab=(partial_r Mscalar²)_ab/sqrt2 in the fixed mass basis; y_hii=-8 Mfermion_i/sqrt2. Eight Goldstone vertices are M_radial²/sqrt2, independently fixed by the radial Ward identity.',
      'channels':cs,'m':m,'tree_radial_mass_squared':M2,'subtraction_point':s0,
      'scheme':'Inverse propagator s-Mtree²+Sigma_sub(s), with Sigma_sub(s0)=Sigmaprime_sub(s0)=0 at spacelike s0=-m². The mass/kinetic matching coefficients are explicit input; this is not a zero-momentum curvature subtraction.',
      'dispersion':'Sigma_sub(s)=(s-s0)²/pi PV integral_threshold^infinity rho(t)/[(t-s0)²(t-s)]dt + i rho(s). Two subtractions converge for scalar and fermion cuts; spacelike subtraction avoids massless-Goldstone zero-momentum logarithms.',
      'real_self_energy_at_tree_pole':float(real),'absorptive_scalar':float(scalar_im),'absorptive_fermion':float(ferm_im),'absorptive_total':float(im),
      'strict_one_loop_pole_squared':{'real':pole.real,'imag':pole.imag},'total_width_leading':float(im/M),'width_over_mass':float(im/M2),
      'Goldstone_width_exact_formula':'Gamma(h->8GG)=M_radial³/(8pi r0²), here r0=1; nonzero and independent of local subtraction counterterms.',
      'dispersion_upper_split_change':float(abs(real-again)),
      'boundaries':['The unstable radial excitation has a complex pole; the previous real one-scalar pole excluded these decay channels.','The fixed-slice derivative vertices define the declared reduced EFT. Nonlinear quotient Kähler geometry, heavy vector multiplets and UV threshold matching are not inherited or solved.','Only strict1loop evaluation at the tree pole is reported; no Dyson iteration, interval error enclosure or complete physical spectrum is claimed.'],
      'prior_owners':['analysis/w33_pass11270_full_condensate_moduli.py','analysis/w33_pass11272_resummed_condensate_radial_control.py','analysis/w33_pass11277_radial_momentum_pole.py'],
      'primary_sources':['https://arxiv.org/html/2505.07931v1','https://arxiv.org/abs/hep-ph/0502168']}
if __name__=='__main__':
    out=payload();(ROOT/'data/w33_pass11282_all_modulus_radial_self_energy.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'],out['strict_one_loop_pole_squared'])
