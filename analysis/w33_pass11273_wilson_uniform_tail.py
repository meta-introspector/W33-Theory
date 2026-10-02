#!/usr/bin/env python3
"""Uniform domination proof for the existing fixed-t Wilson curvature response.

Supplied4-torus only. No spontaneous spacetime dimension or Newton constant.
"""
from pathlib import Path
import sys,json
from fractions import Fraction
import mpmath as mp
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_20261001_wilson_gravity_refinement as W

def uniform_tail(K=64):
    mp.iv.dps=40;t=mp.iv.mpf(1)/5;a=t/5
    def moments(u):
        return [1+mp.iv.sqrt(mp.iv.pi/u),mp.iv.sqrt(mp.iv.pi)/(2*u**mp.iv.mpf('1.5'))+2/(mp.iv.exp(1)*u),3*mp.iv.sqrt(mp.iv.pi)/(4*u**mp.iv.mpf('2.5'))+8/(mp.iv.exp(2)*u**2)]
    S0,S2,S4=moments(a);T0,T2,T4=[mp.iv.exp(-a*K*K/2)*x for x in moments(a/2)]
    # E<=11|p|²; |response summand| <= (A r4+B r2+C) exp(-a r²).
    A=11*121*t*t*mp.iv.exp(2*t/5)
    B=mp.iv.mpf('13.5')*11*t+14*11*t*t*mp.iv.exp(2*t/5)
    C=3*t+4*t*t*mp.iv.exp(2*t/5)
    bound=4*(4*A*(T4*S0**3+3*T0*S4*S0**2)+B*(T2*S0**3+3*T0*S2*S0**2)+C*T0*S0**3)
    sign,man,exp,bc=bound._mpi_[1];rat=Fraction((-1 if sign else 1)*man)*Fraction(2)**exp
    return {'K':K,'upper_exact':str(rat),'upper_display':mp.libmp.to_str(bound._mpi_[1],35),'Gaussian_exponent':'a=t/5','precision':40}

def payload():
    # Broad independent controls of inequalities, not replacements for the proof.
    rng=np.random.default_rng(11273);worst=0.
    for n in (16,32,64):
        h=2*np.pi/n;ps=rng.integers(-n//2,n//2,size=(1000,4));x=h*ps
        E=np.sum(np.sin(x)**2,axis=1)/h**2+np.sum(1-np.cos(x),axis=1)**2/h**2;r2=np.sum(ps**2,axis=1)
        assert np.all(E>=.4*r2-1e-10) and np.all(E<=11*r2+1e-10)
        nextp=ps.copy();nextp[:,0]+=1;xp=h*nextp
        dv=np.column_stack([(np.sin(xp[:,0])-np.sin(x[:,0]))/h,(np.cos(x[:,0])-np.cos(xp[:,0]))/h]);worst=max(worst,float(np.max(np.linalg.norm(dv,axis=1))))
    assert worst<=1+1e-12
    tail=uniform_tail();assert Fraction(tail['upper_exact'])<Fraction(1,10**15)
    fiber=W.finite_fiber(.2)
    return {'status':'PASS','result_scope':'PASS_UNIFORM_WILSON_CURVATURE_DOMINATION_AND_FIXED_T_LIMIT','uniform_tail':tail,'max_neighbor_vector_difference_control':worst,
      'energy_proof':'E_h=sum sin(hp)^2/h²+(sum(1-cos(hp))/h)^2 >=sum4sin(hp/2)^2/h² >=(4/pi²)|p|² >=(2/5)|p|² on the principal Brillouin zone. Also E_h<11|p|² since sin x<=|x|,1-cos x<=x²/2 and h²|p|²<=4pi².',
      'neighbor_proof':'The5-vector changes by norm2|sin(h/2)|/h<=1 under a unit shift in p1. Therefore E_neighbor<=2E+2.',
      'response_bound':'H2<=27E/8+3/4; H1norm<=11E²+14E+4; divided difference <=t exp(-t min(E,E_neighbor)). min>= (2/5)(|p|-1)_+² >= |p|²/5-2/5. Bound by [t(13.5*11*r²+3)+t²(11*121*r4+14*11*r²+4)exp(2t/5)]exp(-t*r²/5).',
      'convergence_proof':'For every fixed integer momentum, sin(hp)/h->p and Wilson mass->0. The divided difference is a continuous integral, including equal energies. Both neighbors converge. The stated summable Gaussian-polynomial dominates every principal-zone summand uniformly; extend by zero outside the zone and use dominated convergence on Z4. Therefore the existing Wilson response converges to its continuum curvature response at fixed t>0.',
      'species_control':'The lower bound fails for the naive wilson0 stencil at corners; it is the Wilson term, not a presumed absence of doublers, that establishes uniform domination.',
      'W33_internal_heat_factor':fiber['heat_factor'],'product_limit':'Multiply the proved external response limit by the fixed finite W33 mass-fiber heat trace; finite internal dimension permits exact heat factorization.',
      'dimension_nonselection_proof':'For any fixed finite162-dimensional mass fiber, Theta_F(t)=162+O(t), so -2t*dlogTheta_F/dt->0. Multiplication by an external d-dimensional continuum heat trace therefore retains UV spectral dimensiond. The same W33 fiber can be attached to any suppliedd; this product ansatz cannot select4 dimensions.',
      'finite_fiber_scope':fiber['scope'],
      'boundary':['Analytic fixed-t limit and directed uniform high-momentum tail, not a numerical finite-N discretization-error enclosure or a uniform t->0 theorem.','Four dimensions, torus, conformal metric and spacing are supplied. W33 is the finite mass fiber; its finite spectrum alone does not derive these geometric inputs.','No Lorentzian gravity, Einstein dynamics, Newton scale or cosmological constant is derived.'],
      'prior_owners':['analysis/w33_20261001_wilson_gravity_refinement.py','analysis/w33_20261001_metric_dynamics_and_heat_bound.py']}
if __name__=='__main__':
    out=payload();(ROOT/'data/w33_pass11273_wilson_uniform_tail.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
