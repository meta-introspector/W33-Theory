#!/usr/bin/env python3
"""One-loop gauge running of the actual declared field inventory; scale boundaries."""
from pathlib import Path
import json
import sympy as s
ROOT=Path(__file__).resolve().parents[1]

def payload():
    # Long roots squared2, tr3(t_a t_b)=delta_ab/2, T_E6(27)=3.
    C6,T6,C3,T3=s.Integer(12),s.Integer(3),s.Integer(3),s.Rational(1,2)
    f6=4*3*T6;sc6=3*T6;f3=4*27*T3;sc3=27*T3
    b6=s.Rational(11,3)*C6-s.Rational(2,3)*f6-s.Rational(1,3)*sc6
    b3=s.Rational(11,3)*C3-s.Rational(2,3)*f3-s.Rational(1,3)*sc3
    assert b6==17 and b3==-s.Rational(59,2)
    # Trace-normalized prior couplings: g6old=sqrt3*g6std, g3old=g3std/sqrt2.
    assert s.sqrt(3)/s.sqrt(6)==1/s.sqrt(2)
    g,t=s.symbols('g t',positive=True)
    inv6=1/g**2+b6*t/(8*s.pi**2);inv3=1/g**2+b3*t/(8*s.pi**2)
    assert s.simplify(inv3-inv6)==-93*t/(16*s.pi**2)
    pole=-8*s.pi**2/(b3*g*g);strong=-8*s.pi**2/(b6*g*g)
    # General inventory scan: n Weyl81 multiplets, one complex scalar81.
    n=s.symbols('n',integer=True,nonnegative=True)
    b6n=41-6*n;b3n=s.Rational(13,2)-9*n
    assert b6n.subs(n,4)==b6 and b3n.subs(n,4)==b3
    d,mu,B,L=s.symbols('d mu B L',positive=True)
    V=L*d**4+B*d**4*(s.log(d*d/(mu*mu))-s.Rational(1,2))
    assert s.diff(V,mu)*mu==-2*B*d**4
    rho=s.symbols('rho')
    return {'status':'PASS_EXACT_DECLARED_UV_BETA_AND_SCALE_BOUNDARY',
      'conventions':'beta(g)=-b0*g^3/(16pi^2); long roots squared2; T27(E6)=3,T3(SU3)=1/2',
      'field_indices':{'E6_Weyl_sum':str(f6),'E6_complex_scalar_sum':str(sc6),'SU3_Weyl_sum':str(f3),'SU3_complex_scalar_sum':str(sc3)},
      'one_loop_b0':{'E6':str(b6),'SU3_family':str(b3)},
      'coupling_bridge':'g_E6_trace=sqrt3*g_E6_standard; g_SU3_trace=g_SU3_standard/sqrt2. The prior mirror-matched ratio is equality of the standard couplings at one scale.',
      'running_from_equal_g0':{'inverse_gE6_squared':str(inv6),'inverse_gSU3_squared':str(inv3),'inverse_difference':str(s.factor(inv3-inv6)),
       'SU3_Landau_log_mu_over_mu0':str(pole),'E6_strong_log_mu_over_mu0':str(strong)},
      'matching_boundary':'Equal standard couplings do not remain equal: the matched mirror gauge spectrum is a boundary condition at one scale, not an RG invariant relation. Yukawa beta functions and broken-phase thresholds are not computed here.',
      'inventory_obstruction':{'n_Weyl81':'bE6=41-6n; bSU3=13/2-9n','simultaneous_asymptotic_freedom':'requires n=0; even one Weyl81 makes the gauged family SU3 infrared free. Singlet spectators cannot change either one-loop gauge coefficient.'},
      'anomalous_potential_RG':{'explicit_scale_derivative':'mu*dV/dmu=-2B*d^4','required_quartic_running':'dL/dlogmu=2B at this order; a beta coefficient does not fix its integration constant.',
       'vacuum_counterterm':'Adding rho_vac leaves nongravitational stationary equations unchanged but changes the gravitational vacuum energy. Its finite renormalized value is independent input.',
       'scale_boundary':'The E6 strong scale mu0*exp[-8pi^2/(17g0^2)] is a candidate threshold, not a computed condensate or spectator/dilaton coupling. Absolute scale still requires a boundary datum.'},
      'prior_owners':['analysis/w33_20261001_complete_gauge_scalar_fermion_loop.py','analysis/w33_20261001_anomalous_scale_completion.py'],
      'primary_sources':['https://arxiv.org/abs/1809.06797','https://arxiv.org/abs/1012.5797','https://doi.org/10.1016/j.nuclphysb.2023.116266']}
if __name__=='__main__':
    out=payload();(ROOT/'data/w33_20261001_uv_running_scale_audit.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
