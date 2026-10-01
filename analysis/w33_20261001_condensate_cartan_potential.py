#!/usr/bin/env python3
"""Conditional near-SUSY condensate action on the canonical three-field slice.

Uses the scoped I18 denominator bridge and canonical Kähler norm on this slice.
Does not assert full11-modulus stability, non-SUSY continuation or SM charges.
"""
from pathlib import Path
import sys,json
import sympy as s
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11255_11259_g26_common as G

def payload():
    r,m,L=s.symbols('r m L',positive=True);V=36*L**18/r**14-18*m*L**9/r**6
    critical={L**9:s.Rational(3,14)*m*r**8}
    # Substitute by treating Lambda^9 as an independent scale datum.
    a=s.symbols('a',positive=True);Va=36*a*a/r**14-18*m*a/r**6
    assert s.factor(s.diff(Va,r).subs(a,s.Rational(3,14)*m*r**8))==0
    h=s.factor((s.diff(Va,r,2)/2).subs(a,s.Rational(3,14)*m*r**8));assert h==s.Rational(648,7)*m*m
    z=s.Symbol('zeta');mod=s.Poly(s.cyclotomic_poly(9,z),z,domain='EX')
    def red(x):return s.rem(s.Poly(s.expand(x),z,domain='EX'),mod).as_expr()
    conj=lambda x:red(s.conjugate(x).subs(s.conjugate(z),z**8))
    fs=G.invariants();vs=G.VARS;at=dict(zip(vs,[1/s.sqrt(3),z/s.sqrt(3),z**8/s.sqrt(3)]));I0=-s.Rational(1,729);W0=s.Rational(3,14)
    g=[red(s.diff(fs[2],x).subs(at)) for x in vs]
    H=[[red(s.diff(fs[2],x,y).subs(at)) for y in vs] for x in vs]
    J=[[[red(s.diff(fs[2],x,y,t).subs(at)) for t in vs] for y in vs] for x in vs]
    w=[red(-W0*g[i]/(3*I0)) for i in range(3)]
    wh=[[red(W0*(s.Rational(4,9)*g[i]*g[j]/I0**2-H[i][j]/(3*I0))) for j in range(3)] for i in range(3)]
    wt=[[[red(W0*(-s.Rational(28,27)*g[i]*g[j]*g[k]/I0**3+s.Rational(4,9)*(H[i][j]*g[k]+H[i][k]*g[j]+H[j][k]*g[i])/I0**2-J[i][j][k]/(3*I0))) for k in range(3)] for j in range(3)] for i in range(3)]
    q=s.Matrix([at[x] for x in vs]);directions=[]
    for f,norm in zip(fs[:2],[4,10]):
        d=s.Matrix([conj(red(s.diff(f,x).subs(at)))/norm for x in vs]);directions +=[d,s.I*d]
    directions +=[q,s.I*q]
    # q perturbation = D*delta/sqrt2; cancel sqrt2 factors in the Hessian formula.
    out=s.zeros(6)
    for a in range(6):
        da=directions[a]
        first=red(sum(sum(wh[k][i]*da[i] for i in range(3))*conj(w[k]) for k in range(3))-9*sum(w[i]*da[i] for i in range(3)))
        assert red(first+conj(first))==0
        for b in range(a+1):
            db=directions[b]
            mixed=red(sum(sum(wh[k][i]*da[i] for i in range(3))*conj(sum(wh[k][j]*db[j] for j in range(3))) for k in range(3)))
            holo=red(sum(sum(wt[k][i][j]*da[i]*db[j] for i in range(3) for j in range(3))*conj(w[k]) for k in range(3))-9*sum(wh[i][j]*da[i]*db[j] for i in range(3) for j in range(3)))
            value=red((mixed+conj(mixed)+holo+conj(holo))/2);out[a,b]=out[b,a]=value
    out=out.applyfunc(s.simplify)
    target=s.diag(*([s.Rational(324,49)]*4+[s.Rational(648,7),s.Rational(486,7)]))
    target[0,2]=target[2,0]=s.Rational(162,49);target[1,3]=target[3,1]=-s.Rational(162,49)
    assert out==target
    expected=[s.Rational(162,49)]*2+[s.Rational(486,49)]*2+[s.Rational(648,7),s.Rational(486,7)]
    assert sorted(out.eigenvals().items(),key=lambda x:x[0])==sorted([(s.Rational(162,49),2),(s.Rational(486,49),2),(s.Rational(648,7),1),(s.Rational(486,7),1)])
    return {'status':'PASS_EXACT_CONDITIONAL_CONDENSATE_SIX_DIRECTION_MINIMUM',
      'declared_superpotential':'W=Lambda^9*(-729*I18)^(-1/3) on the branch W(T)>0',
      'potential':'V_F=sum|dW/dq_i|^2; V_AMSB=m*(sum q_i dW/dq_i-3W)+h.c.=-18m ReW',
      'homogeneity':'W degree-6 in q; phase transforms W->exp(-6i theta)W',
      'radial_potential':str(V),'stationary_radius':'r^8=14*Lambda^9/(3m)',
      'radial_mass_squared':str(h),'exact_slice_Hessian_matrix_over_m_squared':[[str(x) for x in row] for row in out.tolist()],'exact_slice_Hessian_eigenvalues_over_m_squared':list(map(str,expected)),
      'phase_mass_squared_over_m_squared':'486/7','angular_mass_squared_ratio':3,'angular_covariant_mixing':'Real gradI6/gradI12 directions have block(162/49)*[[2,1],[1,2]]; imaginary directions have the opposite offdiagonal sign. Equal-weight eigenvectors are slice scalar modes, not CKM/PMNS mixing.',
      'vacuum_energy':'-(72/7)*m*Lambda^9/r^6, before an independent gravitational counterterm',
      'exact_proof':'Exact cyclotomic restriction with sqrt3/i tangent coefficients and holomorphic derivatives through third order; all six canonical real gradients vanish and the full six-direction Hessian has positive exact eigenvalues and nonzero angular mixing blocks.',
      'physical_boundary':['Conditional condensate denominator/normalization and canonical slice Kähler metric.','Six real Cartan directions only; remaining global-flavor moduli, full Kähler corrections and gauge dynamics are not certified.','AMSB mass m and strong scale Lambda are inputs. Validity of published dynamics concerns small AMSB relative to Lambda, not a proved nonsupersymmetric continuation.','This is not the previous vectorlike EFT and its loop obstruction is not retracted.','The regular T vacuum still cannot retain the Standard Model gauge algebra, and the nonzero vacuum energy still needs a gravitational treatment.'],
      'prior_owners':['analysis/w33_20261001_degree18_phase_completion.py','analysis/w33_20261001_chiral_decuplet_and_condensate_bridge.py'],
      'primary_source':'https://arxiv.org/html/2505.07931v1'}
if __name__=='__main__':
    out=payload();(ROOT/'data/w33_20261001_condensate_cartan_potential.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
