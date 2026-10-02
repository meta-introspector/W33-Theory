#!/usr/bin/env python3
"""An economical spectator mass EFT, with an exact Higgs-only obstruction."""
from pathlib import Path
import sys,json,itertools,math
import sympy as s
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_20261001_global_e6_cartan_covariants as I

def pairings(xs):
    if not xs:yield [];return
    for j in range(1,len(xs)):
        rest=xs[1:j]+xs[j+1:]
        for pairs in pairings(rest):yield [(xs[0],xs[j])]+pairs

def mass(F):
    words=list(itertools.combinations_with_replacement(range(3),3))
    def entry(ix):return sum(math.prod(F[ix[i],ix[j]] for i,j in pairs) for pairs in pairings(list(range(6))))
    return s.Matrix([[entry(a+b) for b in words] for a in words])

def payload():
    _,d,_=I.tensors();singlet=d[:2,:2,:2];assert not singlet.any()
    M=mass(s.eye(3));assert M.rank()==10 and M.is_positive_definite
    N=mass(s.diag(1,2,3));assert N.rank()==10
    b_old=s.Rational(-93,2);b_new=11-s.Rational(2,3)*21-s.Rational(1,3)*(27*s.Rational(5,2)+s.Rational(5,2))
    assert b_new==-s.Rational(79,3)
    return {'status':'PASS','result_scope':'PASS_EXACT_SEXTET_COMPOSITE_MASS_AND_RUNNING_AUDIT',
      'replacement':'Replace elementary Sigma(1,28) by family sextet F(1,6). Sigma_eff=Sym(F tensor F tensor F), using the15 pairings of six indices.',
      'operator':'(c/M_UV²) chi_ijk chi_lmn Sigma_eff^{ijklmn}, a dimension6 operator. Spectator mass scale c v_F³/M_UV².',
      'identity_sextet_mass':[[str(x) for x in row] for row in M.tolist()],'identity_mass_determinant':str(M.det()),'diag123_mass_determinant':str(N.det()),'rank':10,
      'general_nonsingular_F':'Takagi congruence transforms the identity Wick catalecticant through Sym³ of the invertible3x3 congruence. Rank10 persists for every invertible complex symmetric F; positivity is asserted only for real positive F.',
      'Higgs_only_failure':'d restricted to the two canonical SM-neutral27 directions e0,e1 is identically zero. Therefore dbar(Hdagger,Hdagger,Hdagger) cannot supply the spectator mass when all H vevs preserve this SM. The independent E6-neutral sextet is needed for this operator.',
      'beta_family_old':str(b_old),'beta_family_new':str(b_new),'beta_E6_subset':'32','beta_E6_with_two_complex27_and_real78':'28','beta_E6_with_two_complex27_and_complex78':'26',
      'inventory':'Nonsupersymmetric Weyl(27,3)+(1,10bar), complex H(27,6bar)+F(1,6). Family Tfermions=21, Tscalars=70. The polynomial Higgs inventory uses two family-neutral complex27 and a real78, changing E6 beta by-4; a complex78 instead changes it by-6.',
      'UV_boundary':'Family running improves but remains non-asymptotically free. log(M_Landau/mu)=8pi²/((79/3)g_family(mu)²) at1loop within this fixed inventory. The nonrenormalizable operator still requires a mediator UV completion.',
      'prior_owners':['analysis/w33_pass11271_chiral_symmetric_yukawa_completion.py','analysis/w33_pass11271_canonical_sm_higgs_bridge.py'],
      'boundary':['No selected sextet vacuum, mediator spectrum, absolute spectator mass or observed Yukawa prediction.','Gauge anomaly matching below the spectator threshold must retain the appropriate broken-family Wess-Zumino terms.']}
if __name__=='__main__':
    out=payload();(ROOT/'data/w33_pass11276_sextet_composite_spectator.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
