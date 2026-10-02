#!/usr/bin/env python3
"""A narrow integer AF window ties the factor rank to the composite-Higgs inventory."""
from pathlib import Path
import json
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
def payload():
 r=s.symbols('r',integer=True,positive=True)
 # Three complex(27,2r) plus real adjoint SO(2r).
 aux=s.expand(s.Rational(11,3)*(2*r-2)-27-s.Rational(1,6)*(2*r-2));assert aux==7*r-34
 old=28-6*r;new=34-6*r
 oldgood=[n for n in range(3,100) if old.subs(r,n)>0 and aux.subs(r,n)>0]
 newgood=[n for n in range(3,100) if new.subs(r,n)>0 and aux.subs(r,n)>0]
 assert oldgood==[] and newgood==[5]
 return {'status':'PASS','result_scope':'PASS_EXACT_INTEGER_AF_WINDOW_FOR_THE_FIXED_THREE_FRAME_ARCHITECTURE',
 'auxiliary_group':'SO(2r), three complex(27,2r) matrices and one real antisymmetric adjoint; E6 trace conventions T27=3, vector T=1, C2(SO2r)=2r-2.',
 'beta_polynomials':{'auxiliary_SO2r':'7r-34','E6_elementary_H':'28-6r','E6_composite_H':'34-6r'},
 'analytic_window':'Auxiliary AF requires integer r>=5. Elementary-H E6 AF requires r<=4, so no intersection. Composite-H E6 AF requires r<=5, giving exactly r=5. This simultaneously matches the previous cubic selector rank5 lower bound.',
 'unique_r':5,'at_r5':{'E6':4,'SO10':1},
 'boundaries':['Uniqueness only within the specified three-frame plus real-adjoint architecture and fixed baseline inventory.','Additional mediators, scalar fields or Weyl matter change the inequalities. Family SU3 is not made AF by this argument.','This selects an allowable factor rank, not physical gauge couplings, Higgs scales or the Standard Model itself.'],
 'prior_owners':['analysis/w33_pass11280_factor_higgs_mediator.py','analysis/w33_pass11285_semisimple_factor_higgs.py']}
if __name__=='__main__':
 out=payload();(ROOT/'data/w33_pass11291_rank_five_AF_window.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
