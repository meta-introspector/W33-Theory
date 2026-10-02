#!/usr/bin/env python3
"""Actual1-loop bubble and pole in the explicitly declared radial-only scalar EFT.

This is not the full11-modulus supersymmetric EFT or a Standard Model mass.
"""
from pathlib import Path
import json
import sympy as s
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[1]

def bubble(z):
    assert 0<=z<4
    return quad(lambda x:-np.log1p(-z*x*(1-x)),0,1,epsabs=1e-13)[0]

def payload(m=.01):
    r=s.symbols('r',positive=True);V=s.Rational(81,49)*r**-14-s.Rational(27,7)*r**-6
    M2=float(s.diff(V,r,2).subs(r,1)/2)*m*m
    g3=float(s.diff(V,r,3).subs(r,1)/(2*s.sqrt(2)))*m*m
    c=g3*g3/(32*np.pi**2)
    # Gamma_E(pE²)=pE²+Mcurv²-g3²/2 [I_E(pE²)-I_E(0)].
    # Wick rotation gives inverse s-Mcurv²+c B(s/Mloop²), with B positive below threshold.
    inv=lambda x:x-M2+c*bubble(x/M2)
    pole=brentq(inv,.5*M2,M2,xtol=1e-15)
    deriv=quad(lambda x:x*(1-x)/(1-pole/M2*x*(1-x)),0,1,epsabs=1e-13)[0]/M2
    residue=1/(1+c*deriv)
    assert 0<pole<M2 and 0<residue<1 and abs(inv(pole))<1e-13
    at_tree=M2-c*bubble(1)
    return {'status':'PASS','result_scope':'PASS_ONE_LOOP_MOMENTUM_POLE_IN_DECLARED_RADIAL_ONLY_SCALAR_EFT',
      'action':'L=1/2(partial h)²-V(1+h/sqrt2), V(r)=m²[81/(49r14)-27/(7r6)]. This is a separately declared single-real-scalar EFT extracted from the canonical radial potential.',
      'renormalization':'Tadpole fixed at r=1; curvature mass fixed at Vhh(0); local kinetic coefficient set to1 in this zero-momentum potential subtraction scheme. Momentum bubble kept explicitly; quartic tadpole is momentum independent and subtracted.',
      'm':m,'curvature_mass_squared':M2,'cubic_vertex':g3,'bubble_prefactor':c,
      'momentum_function':'B(z)=-int_0^1 log(1-z*x*(1-x)) dx; finite, positive below two-particle threshold z=4.',
      'inverse_propagator':'s-Mcurv²+(g3²/(32pi²)) B(s/Mcurv²)',
      'one_loop_pole_squared_strict':at_tree,'one_loop_Dyson_root_squared':pole,'residue_at_Dyson_root':residue,'pole_equation_residual':inv(pole),
      'relative_pole_shift':(pole-M2)/M2,'difference_Dyson_vs_strict':pole-at_tree,
      'scope':'A nonzero physical propagator pole in this declared scalar EFT, demonstrating why static curvature and pole differ. Dyson iteration includes selected higher-order terms; only the strict1loop expansion has1loop perturbative accuracy.',
      'excluded':['Ten other complex-modulus directions, Goldstone couplings, Weyl loops, gauge multiplets, noncanonical field-space metric and threshold matching.','This does not close the full11-modulus nonzero-pole frontier; no observed particle assignment or mass scale is predicted.'],
      'prior_owners':['analysis/w33_pass11272_resummed_condensate_radial_control.py','analysis/w33_pass11272_goldstone_ir_and_pole_audit.py'],
      'primary_sources':['https://arxiv.org/abs/hep-ph/0502168']}
if __name__=='__main__':
    out=payload();(ROOT/'data/w33_pass11277_radial_momentum_pole.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'],out['one_loop_Dyson_root_squared'])
