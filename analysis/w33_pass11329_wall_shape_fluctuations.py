"""Coupled Einstein-Maxwell-wall torus-shape sector: zero modes and spectral control."""
from pathlib import Path
import sys,json
import numpy as np
import sympy as s
from scipy.linalg import eigh
from numpy.polynomial.legendre import Legendre,leggauss
ROOT=Path(__file__).resolve().parents[1]

def spectrum(m,parity,n,quadrature=350):
 L=.25;c=.2;z,weights=leggauss(quadrature);x=(z+1)*L/2;b=1+x;weights*=L/2;F=8/(15*b)-b*b/30-1/(2*b*b);a=m/2;cols=[];ders=[]
 for k in range(n):
  P=Legendre.basis(k);v=P(z);dv=P.deriv()(z)*2/L;u=x**a*v;du=x**a*dv+(a*x**(a-1)*v if a else 0)
  if parity=='odd':du=du*(L-x)-u;u=u*(L-x)
  cols.append(u);ders.append(du)
 U=np.array(cols).T;D=np.array(ders).T;M=U.T@((weights*b*b)[:,None]*U);K=D.T@((weights*b*b*F)[:,None]*D)
 if m:K+=U.T@((weights*m*m*c*c*b*b/F)[:,None]*U)
 ev=eigh(K,M,eigvals_only=True);assert min(ev)>-1e-9
 return ev

def scalar_curvature_control():
 # Exact4D Ricci scalar for diagonal torus deformation with fixed determinant.
 x=s.symbols('x');r=s.Function('r')(x);b=s.Function('b')(x);u=s.Function('u')(x)
 # Flat fiber dimensions1+1+1: R=-2 sum ai''/ai -2 sum_i<j ai' aj'/(ai aj).
 f=[r,b*s.exp(u),b*s.exp(-u)];R=-2*sum(s.diff(v,x,2)/v for v in f)-2*sum(s.diff(f[i],x)*s.diff(f[j],x)/(f[i]*f[j]) for i in range(3) for j in range(i))
 f0=[r,b,b];R0=-2*sum(s.diff(v,x,2)/v for v in f0)-2*sum(s.diff(f0[i],x)*s.diff(f0[j],x)/(f0[i]*f0[j]) for i in range(3) for j in range(i));diff=s.simplify(R-R0);assert s.simplify(diff+2*s.diff(u,x)**2)==0
 return str(diff)

def payload():
 rows=[]
 for m in range(4):
  for parity in ['even','odd']:
   first=spectrum(m,parity,16)[:6];second=spectrum(m,parity,24)[:6];third=spectrum(m,parity,32)[:6];error=float(max(abs(second-third)));assert error<1e-5
   rows.append({'azimuthal_integer':m,'wall_parity':parity,'eigenvalues_basis32':third.tolist(),'basis16_to24_error':float(max(abs(first-second))),'basis24_to32_error':error,'multiplicity_each_radial_level':2*(1 if m==0 else 2)})
 assert abs(rows[0]['eigenvalues_basis32'][0])<1e-9
 return {'status':'PASS','result_scope':'PASS_FULL_COUPLED_TORUS_SHAPE_BLOCK_NONNEGATIVE_WITH_TWO_EXACT_ZERO_MODES_NOT_TOTAL_WALL_STABILITY','deformation':'Replace the torus metric b²I2 by b²Q, Q symmetric positive2x2 with detQ=1. Its two independent linear shape deformations are tracefree. Keep fields torus-translation invariant and allow xi/phi dependence.','coupled_sector':'The background Einstein-Maxwell-fourform-wall action has no shape potential: bulk/wall volume and magnetic2-flux depend on detQ only. Torus traceless metric equations decouple at linear order from the breathing mode, magnetic Maxwell fluctuation, base metric and wall displacement. These are genuine metric fluctuations of the full supplied action, not an added spectator scalar.','exact_Ricci_difference_diagonal_shape':scalar_curvature_control(),'quadratic_action':'For diagQ=(exp(2u),exp(-2u)), R=Rbackground-2|grad_base u|². Euclidean -kappa² R/2 therefore yields a positive quadratic gradient form. The second torus-shape polarization has the same linear operator. Pure-tension wall produces no shape mass and continuity of normal flux.','radial_operator':'b in[1,5/4], F=8/(15b)-b²/30-1/(2b²), c=1/5. -(b²F u_b)_b + m² c² b²/F u = lambda b² u. Regular pole u~(b-1)^(|m|/2). Reflection-even wall has u_b=0, odd wall u=0.','spectral_rows':rows,'positivity_proof':'b²F=(16b-b^4-15)/30 vanishes atb1 and has positive derivative(16-4b³)/30 on[1,5/4]. HenceF>0 inside. Multiply the radial equation by ubar and integrate. With the regular pole and even/odd wall conditions, lambda int b²|u|²=int b²F|uprime|²+m²c²b²/F|u|² >=0. This proves absence of negative modes throughout this entire shape sector for every integerm, beyond the sampled spectrum.','exact_zero_modes':{'count':2,'sector':'m0, reflection-even, constant torus shape','interpretation':'Two torus complex-structure moduli survive. Fixed coordinate periods and discrete SL2Z identifications do not remove their local continuous fluctuations. Canonical matter, wall tension and global volume/fourflux constraints do not lift them.'},'scope':['All coupled torus-shape block modes are nonnegative, with two exact flat directions. This is not full coupled fluctuation stability.','The breathing/base-metric/Maxwell/wall-bending block, four-form jump variations and conformal-factor treatment remain open.','Numerical Rayleigh-Ritz levels have basis controls, not interval-certified error bounds. They are Euclidean fluctuation eigenvalues, not observed Lorentzian pole masses.','The earlier one-junction negative radial curvature is not promoted to a physical mode. No nucleation determinant or observedCC is derived.'],'prior_owners':['analysis/w33_pass11324_quantized_warped_history_wall.py','analysis/w33_pass11302_sequestered_membrane_junction.py']}
if __name__=='__main__':
 d=payload();(ROOT/'data/w33_pass11329_wall_shape_fluctuations.json').write_text(json.dumps(d,indent=2)+'\n');print(d['status'],[(r['azimuthal_integer'],r['wall_parity'],r['eigenvalues_basis32'][:2]) for r in d['spectral_rows']])
