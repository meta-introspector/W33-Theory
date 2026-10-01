#!/usr/bin/env python3
"""A specified anomalous two-field scale completion, with exact stability checks.

The transmutation scale, renormalization conditions and EFT coefficients are
inputs. This is an explicit reduced mechanism, not a derived UV completion.
"""
from pathlib import Path
import json,argparse
import sympy as s
OUT=Path(__file__).resolve().parents[1]/'data/w33_20261001_anomalous_scale_completion.json'

def payload():
    rho,d,v,lam,B,mu=s.symbols('rho d v lambda B mu',positive=True)
    r=rho/s.sqrt(2)
    potential=lam*(r*r-v*v*d*d)**2+B*d**4*(s.log(d*d/(mu*mu))-s.Rational(1,2))
    at={rho:s.sqrt(2)*v*mu,d:mu}
    gradient=s.Matrix([s.diff(potential,x) for x in (rho,d)]).subs(at).applyfunc(s.simplify)
    H=s.hessian(potential,(rho,d)).subs(at).applyfunc(s.simplify)
    assert gradient==s.zeros(2,1)
    assert H[0,0]==4*lam*v*v*mu*mu
    assert s.factor(H.det())==32*B*lam*mu**4*v**2
    assert s.simplify(potential.subs(at))==-B*mu**4/2
    # Exactly declared control values; no fit to measured physical numbers.
    control=H.subs({v:s.Rational(1,5),lam:1,B:s.Rational(1,100),mu:1})
    vals=sorted(float(x) for x in control.eigenvals());assert vals[0]>0
    return {'schema':'w33.20261001.anomalous-scale-completion.v1',
      'status':'PASS_EXACT_REDUCED_RADIAL_DILATON_STABILITY_WITH_INPUT_TRANSMUTATION_SCALE',
      'named_potential':str(potential),'canonical_coordinates':'rho=sqrt(2)||Phi||, d=real dilaton with kinetic(1/2)(dd)^2',
      'stationary_point':'||Phi||=v mu, d=mu',
      'exact_Hessian':[[str(x) for x in row] for row in H.tolist()],
      'exact_leading_minor':str(H[0,0]),'exact_determinant':str(s.factor(H.det())),
      'stability':'Both radial/dilaton eigenvalues strictly positive for lambda,B,v,mu>0',
      'dimensionless_control_eigenvalues':vals,'vacuum_energy':'-B mu^4/2 before an independent constant',
      'possible_log_source':'n real singlet spectators X_a with m_X^2=h^2 d^2 give B_X=n h^4/(64pi^2). A d^4 counterterm and renormalization condition set the displayed logarithmic constant; the full matter contribution must also be included to determine B_total.',
      'angular_extension':'alpha6 |I6(Phi)|^2/d^8 + alpha12 |I12(Phi)|^2/d^20 + moment-square potential; I_d are the named theta-invariant extensions, alpha_d>0',
      'canonical_tree_angular_masses_squared':['16 alpha6 v^10 mu^2 twice','100 alpha12 v^22 mu^2 twice'],
      'local_loop_robustness':'At T a sufficient condition is min(16alpha6 v^10,100alpha12 v^22) > |4g^4/9-8|y|^4| v^2 ||H_F||/(128pi^2), with the fixed-radius mirror Hessian H_F. This is a local condition; global quantum selection is not proved.',
      'cutoff_links':'Lambda_heat=alpha d, Lambda_EFT=beta d; all three named scales share d, but alpha,beta,v remain free dimensionless parameters',
      'predicted_form_not_numbers':'m_light/m_heavy carries powers(v/beta)^1,(v/beta)^9,(v/beta)^21; M_Pl/mu follows the specified positive heat factor and alpha,kappa, not W33 arithmetic alone',
      'boundary':['An explicit reduced anomalous potential with exact local stability, not a UV-derived beta function or the complete quantum action.','The renormalization scale/condition is supplied; measured absolute scales and coupling ratios are not predicted.','The nonzero vacuum energy is not a small cosmological-constant solution.','The common phase remains flat unless separately gauged or lifted.','The global invariant circuits are built in the companion covariant packet; this reduced model does not compute their complete matter-loop corrections.'],
      'prior_owners':['analysis/w33_20261001_complete_gauge_scalar_fermion_loop.py','analysis/w33_20261001_exact_cartan_mirror_completion.py','analysis/w33_20261001_g26_global_polynomial_vacuum.py'],
      'external_sources':['https://arxiv.org/abs/hep-th/0512169','https://doi.org/10.1103/PhysRevD.7.1888']}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();out=payload()
    if a.check:assert out==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(out['status'])
