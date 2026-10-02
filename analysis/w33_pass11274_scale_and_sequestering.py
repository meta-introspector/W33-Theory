#!/usr/bin/env python3
"""Near-SUSY scale hierarchy, spectral CC obstruction, conditional sequestering."""
from pathlib import Path
import json
import sympy as s
ROOT=Path(__file__).resolve().parents[1]

def payload():
    # Canonical E6/global-family condensate spectrum: Goh et al. table4;
    # producer11270 independently checks its W33 realization.
    bos=[s.Rational(162,49)]*10+[s.Rational(486,49)]*2+[s.Rational(648,7),s.Rational(486,7)]+[s.Integer(0)]*8
    ferm=[s.Rational(81,49)]*8+[s.Rational(324,49)]*2+[s.Integer(81)]
    st2=s.factor(sum(bos)-2*sum(ferm));st4=s.factor(sum(x*x for x in bos)-2*sum(x*x for x in ferm));assert st2==0 and st4>0
    mu,L,m=s.symbols('mu Lambda m',positive=True);g=s.symbols('g',positive=True)
    strong=mu*s.exp(-8*s.pi**2/(27*g*g)) # SUSY E6 with3 chiral27:36-9=27.
    radial=(s.Rational(14,3)*L**9/m)**s.Rational(1,8)
    energy=s.factor(-s.Rational(72,7)*m*L**9/radial**6)
    assert s.simplify(radial/L-(s.Rational(14,3)*L/m)**s.Rational(1,8))==0
    # Four-volume average two epochs; common constant shift cancels exactly.
    c,u,v=s.symbols('c u v');weights=[s.Rational(1,4),s.Rational(3,4)];rho=[u,v]
    avg=sum(w*x for w,x in zip(weights,rho));residual=[s.expand(x-avg) for x in rho]
    shiftedavg=sum(w*(x+c) for w,x in zip(weights,rho));assert all(s.expand((x+c)-shiftedavg-r)==0 for x,r in zip(rho,residual))
    return {'status':'PASS','result_scope':'PASS_EXPLICIT_SCALE_INPUT_AND_CONDITIONAL_VACUUM_SHIFT_TEST',
      'near_SUSY_beta':'b_E6=3C_A-3T27=36-9=27, including gaugino and3 chiral supermultiplets. Distinct from b17 paired nonsupersymmetric EFT and b35 chiral nonsupersymmetric subset.',
      'Lambda_boundary_relation':str(strong),'radius_over_Lambda':'(14Lambda/(3m))^(1/8)','weak_coupling_control':'m/Lambda->0 gives r/Lambda->infinity; canonical UV-field Kähler treatment improves parametrically, but neither m nor the boundary coupling is selected by W33.',
      'stationary_vacuum_energy':str(energy),'vacuum_scaling':'Vmin proportional -(m^7 Lambda^9)^(1/4)',
      'canonical_supertrace_m_squared':str(st2),'canonical_supertrace_m_fourth':str(st4),'loop_vacuum_scale_dependence':'dV1/dlogmu= -STr(M^4)/(32pi²). Vanishing STr(M²) does not cancel the vacuum energy or its running.',
      'sequestering_declared_action':'S=int sqrt(-g)[Mpl²R/2-Lambda_g-lambda^4 L_m(lambda^-2 g^{mu nu},Phi)]+sigma(Lambda_g/(lambda^4 mu_s^4)). Lambda_g and lambda are spacetime-global variables.',
      'conditional_gravity_equation':'Mpl²G_mu_nu=T_mu_nu-(1/4)g_mu_nu<TrT>. A spacetime-constant shift T_mu_nu->T_mu_nu-Cg_mu_nu cancels exactly; residual history/flux and the global sigma constraints remain inputs.',
      'two_epoch_test':{'weights':list(map(str,weights)),'residual_density':list(map(str,residual)),'constant_shift_cancels':True},
      'W33_connection':'The condensate constant Vmin and constant vacuum loop shifts are canceled only if this additional global constraint is adopted. W33 has not produced lambda, Lambda_g, sigma or the spacetime averaging prescription.',
      'boundary':['A tested external mechanism, not a W33-derived CC solution.','History-dependent residual is not forced to the observed value. Nonconstant condensate stress and phase transitions remain gravitational sources.','Absolute mass and Newton scales still require a UV boundary coupling, SUSY-breaking sector and gravitational dynamics.'],
      'prior_owners':['analysis/w33_20261001_uv_running_scale_audit.py','analysis/w33_20261001_condensate_cartan_potential.py'],'primary_sources':['https://arxiv.org/html/2505.07931v1','https://arxiv.org/abs/1309.6562']}
if __name__=='__main__':
    out=payload();(ROOT/'data/w33_pass11274_scale_and_sequestering.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
