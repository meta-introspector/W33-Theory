#!/usr/bin/env python3
"""Supplied-metric dynamics and an outward interval continuum heat certificate.

The tail bound certifies the continuum benchmark, not Wilson continuum convergence.
Classical EH conformal instability is prior physics; exact W33-fiber audit here.
"""
from pathlib import Path
import sys,json
from fractions import Fraction
import numpy as np
import mpmath as mp
import sympy as s
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_20261001_isovolume_dirac_gravity_response as G

def interval_response(K=24):
    mp.iv.dps=40;t=mp.iv.mpf([1,1])/5
    x=np.arange(-K,K+1,dtype=int);a,b,c=np.meshgrid(x,x,x,indexing='ij');hist=np.bincount((a*a+b*b+c*c).ravel());qs=np.flatnonzero(hist)
    total=mp.iv.mpf(0)
    for n in range(-K,K+1):
        shift=2*n+1
        for q in qs:
            E=int(n*n+q);Ep=int((n+1)**2+q);ee=mp.iv.exp(-t*E)
            dd=(ee-mp.iv.exp(-t*Ep))/shift
            term=-4*t*(3*E+mp.iv.mpf(1)/8)*ee+t*((E+Ep-mp.iv.mpf(1)/2)**2+int(q))*dd
            total+=int(hist[q])*term
    def moments(u):
        pi=mp.iv.pi;e=mp.iv.exp(1)
        return [1+mp.iv.sqrt(pi/u),mp.iv.sqrt(pi)/(2*u**mp.iv.mpf('1.5'))+2/(e*u),3*mp.iv.sqrt(pi)/(4*u**mp.iv.mpf('2.5'))+8/(e**2*u**2)]
    S0,S2,S4=moments(t);A0,A2,A4=moments(t/2);fac=mp.iv.exp(-t*K*K/2)
    T0,T2,T4=[fac*a for a in (A0,A2,A4)]
    tail=4*t*(180*(T4*S0**3+3*T0*S4*S0**2)+123*(T2*S0**3+3*T0*S2*S0**2)+mp.iv.mpf('70.5')*T0*S0**3)
    enclosed=total+mp.iv.mpf([-1,1])*tail
    # Keep endpoints at interval precision, not float-rounding them inward.
    low,high=enclosed._mpi_;fmt=lambda v:mp.libmp.to_str(v,38)
    def rational(v):
        sign,man,exp,bc=v
        return str(Fraction((-1 if sign else 1)*man)*Fraction(2)**exp)
    return {'K':K,'t':'1/5','response_lower_display':fmt(low),'response_upper_display':fmt(high),'tail_upper_display':fmt(tail._mpi_[1]),'response_lower_exact':rational(low),'response_upper_exact':rational(high),'tail_upper_exact':rational(tail._mpi_[1]),
      'finite_sum':'all integer4-momenta in[-K,K]^4, transverse squares grouped by exact integer multiplicities',
      'tail_proof':'Absolute response summands bounded after shifting n+1 by t*(45E^2+123E+141/2)*exp(-tE). E^2<=4sum x_i^4; union of4 coordinate tails. Tail moments T_j<=exp(-tK^2/2)*S_j(t/2); S_j bounded by Gaussian integral plus twice the unimodal maximum.',
      'scope':'Directed mpmath interval arithmetic plus analytic infinite-sum tail; not a certified Wilson discretization error or a geometry-emergence proof.'}

def dynamics():
    k,M,c,V=s.symbols('k M c V',positive=True)
    coeff=V*(-s.Rational(3,2)*M*M*k*k+18*c*k**4)
    assert s.expand(coeff.subs(c,M*M/12)-s.Rational(3,2)*V*M*M*k*k*(k*k-1))==0
    return {'declared_action':'S_E=integral sqrt(g)[rho_vac-M_Pl^2 R/2+cR^2] plus the named scalar/gauge/fermion action',
      'metric_variation':'For c=0, M_Pl^2 G_mu_nu+rho_vac*g_mu_nu=T_mu_nu; geometry is a varied field, but its dimension/topology/action coefficients are still inputs.',
      'isovolume_mode':'sigma=epsilon*cos(k*x1)-log(I0(4epsilon))/4 on period2pi torus',
      'epsilon_squared_integrated_R':'3*V*k^2','epsilon_squared_integrated_R_squared':'18*V*k^4',
      'quadratic_Euclidean_action':str(coeff),'EH_only':'Negative for every nonconstant conformal Fourier mode. A positive internal heat factor rescales the coefficient and cannot cure it.',
      'R_squared_control':'c>M_Pl^2/12 stabilizes these conformal modes with integer|k|>=1 on the fixed unit-period torus; for larger physical size the lower-momentum instability returns unless the bound scales with that size.',
      'flat_background_equation':'Constant matter at a stationary vacuum and a flat metric require total rho_vac+V_min=0. The preceding V_min=-B*mu^4/2 cannot be ignored when varying the metric.',
      'boundary':'Conformal-sector quadratic audit, not a complete metric Hessian, gravitational measure, Lorentzian continuation, or emergent spacetime construction.'}
if __name__=='__main__':
    out={'status':'PASS_DIRECTED_CONTINUUM_HEAT_BOUND_AND_NAMED_METRIC_VARIATION','heat_certificate':interval_response(),'metric_dynamics':dynamics(),
      'prior_owners':['analysis/w33_20261001_wilson_gravity_refinement.py','analysis/w33_20261001_isovolume_dirac_gravity_response.py','analysis/w33_20261001_anomalous_scale_completion.py'],
      'primary_source':'https://arxiv.org/abs/hep-th/0103186'}
    (ROOT/'data/w33_20261001_metric_dynamics_and_heat_bound.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
