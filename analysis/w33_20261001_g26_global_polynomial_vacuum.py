#!/usr/bin/env python3
"""Exact global polynomial vacuum selector, using the prior G26 invariants.

Prior owners: Pass11256 invariant potentials, Pass11269 invariant dictionary,
20261001 balanced barrier T-magic stationary point and unitary kinetic form.
New: a reduced 72-point complete intersection selects the entire known orbit
as global minima without logarithmic singularities or fitted angular targets.
"""
from pathlib import Path
import sys,json,argparse
import sympy as s
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11255_11259_g26_common as G
OUT=ROOT/'data/w33_20261001_g26_global_polynomial_vacuum.json'

def payload():
    x,y,z=G.VARS;u6,u12,u18=G.invariants()
    A,B=s.symbols('A B')
    def cubes(f):
        out=0
        for ex,c in s.Poly(f,x,y,z).terms():
            assert all(n%3==0 for n in ex)
            out+=c*A**(ex[0]//3)*B**(ex[1]//3) # chart z=1
        return s.expand(out)
    p,q=cubes(u6),cubes(u12)
    gb=s.groebner([p,q],B,A)
    linear,elim=[g.as_expr() for g in gb.polys]
    assert s.degree(linear,B)==1 and s.degree(elim,A)==8
    assert s.diff(linear,B)==1539
    assert s.gcd(elim,s.diff(elim,A))==1 and elim.subs(A,0)==1
    # No chart-at-infinity or vanishing coordinate roots.
    assert s.gcd(p.subs(B,0),q.subs(B,0))==1
    inf6=s.Poly(u6.subs(z,0).subs(y,1),x)
    inf12=s.Poly(u12.subs(z,0).subs(y,1),x)
    assert s.gcd(inf6,inf12).as_expr()==1
    assert u6.subs({x:1,y:0,z:0})==1
    assert s.gcd(u6,u12)==1
    # All 72 intersections avoid the reflection divisor.
    u9=(A-B)*(A-1)*(B-1)
    assert list(s.groebner([p,q,u9],B,A))==[s.Integer(1)]
    # Exact verification at the T ray in Q(zeta9).
    r=s.symbols('r');mod=s.Poly(s.cyclotomic_poly(9,r),r)
    red=lambda f:s.rem(s.Poly(s.expand(f),r),mod).as_expr()
    at={x:1,y:r,z:r**8}
    assert red(u6.subs(at))==0 and red(u12.subs(at))==0
    assert red(u18.subs(at))==-27
    jac=s.Matrix([[s.diff(f,v) for v in (x,y,z)] for f in (u6,u12,u18)])
    dj=red(jac.det().subs(at));assert dj!=0
    angular_minor=red(s.det(jac[:2,:2]).subs(at));assert angular_minor!=0
    rows=jac[:2,:].applyfunc(lambda f:red(f.subs(at)))
    gram=(rows*rows.applyfunc(lambda f:red(f.subs(r,r**8))).T).applyfunc(red)
    assert gram==s.diag(3888,17714700)
    assert [gram[i,i]/3**(d-1) for i,d in enumerate([6,12])]==[16,100]
    return {
      'schema':'w33.20261001.g26-global-polynomial-vacuum.v1',
      'status':'PASS_EXACT_72_REDUCED_GLOBAL_POLYNOMIAL_VACUUM_RAYS',
      'potential':{
        'formula':'V(v)=lambda*(v^dag v-r0^2)^2 + alpha*|u6(v)|^2 + beta*|u12(v)|^2',
        'parameters':'lambda,alpha,beta,r0 > 0 (appropriate powers of a cutoff for dimensionful fields)',
        'degrees':[4,12,24],
        'global_minimum_value':0,
        'unit_T_angular_gradient_Gram':[[16,0],[0,100]],
        'real_Hessian_spectrum':'0, 8*lambda*r0^2, 32*alpha*r0^10 (twice), 200*beta*r0^22 (twice)',
        'Hessian_scope':'Euclidean real kinetic coordinates of the standard unitary representation; coefficients are not physical predictions',
        'minimizers':'72 projective rays at norm r0; each retains a common U(1) phase',
        'angular_coefficient_independence':True,
        'proof':'All summands are nonnegative and T has u6=u12=0. The reduced complete intersection consists of 72 rays. Radial square fixes the norm. The angular Hessian is 2 Re(J^dag diag(alpha,beta) J), positive on the projective tangent; no global numerical optimization is used.'
      },
      'exact_complete_intersection':{
        'chart':'z=1; A=x^3, B=y^3',
        'groebner_basis':list(map(str,[linear,elim])),
        'eliminant_factorization':str(s.factor(elim)),
        'eliminant_discriminant':str(s.discriminant(elim,A)),
        'cube_ratio_solutions':8,'cube_lifts_per_solution':9,
        'projective_solutions':72,'all_reduced':True,
        'no_coordinate_zero':True,'no_points_at_infinity':True,
        'no_u9_zero':True,'u18_at_T_raw':'-27',
        'full_invariant_jacobian_at_T':str(dj),
        'angular_jacobian_minor_at_T':str(angular_minor)
      },
      'orbit_identification':{
        'argument':'The prior standard G26 invariant ring is C[u6,u12,u18]. Every nonzero common zero has nonzero u18 and can be scaled to the T value. Fibers of a finite-group invariant quotient are single group orbits. Therefore these 72 rays are the projective G26/Clifford T orbit, independently of numerical orbit counting.',
        'prior_numeric_crosscheck':'data/w33_20261001_balanced_g26_tmagic_vacuum.json'
      },
      'physical_boundary':[
        'A constructed scalar effective potential, not an action uniquely derived from W33.',
        'Degree-12 and degree-24 interactions are nonrenormalizable in four dimensions; a cutoff or UV mechanism is required.',
        'Positive coefficients and radius are external inputs; angular vacuum location is independent of their values, masses are not.',
        'The U(1) phase is flat unless gauged or lifted; no CP/time orientation is selected.',
        'A classical scalar vacuum is not a prepared qutrit state or a fault-tolerant magic factory.',
        'The chosen zero of V is not a prediction for the cosmological constant; additive vacuum energy and quantum corrections are unconstrained.'
      ],
      'prior_owners':['analysis/w33_pass11256_g26_invariant_potential.py','analysis/w33_pass11269_g26_qutrit_dictionary.py','analysis/w33_20261001_balanced_g26_tmagic_vacuum.py'],
      'checks':{'coprime':True,'eight_simple_cube_solutions':True,'all_9_lifts_unramified':True,'T_zero_exact':True,'T_regular_exact':True,'global_sum_of_squares':True}
    }

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--check',action='store_true');a=ap.parse_args()
    p=payload();txt=json.dumps(p,indent=2,sort_keys=True)+'\n'
    if a.check:assert OUT.read_text()==txt,'certificate drift'
    else:OUT.write_text(txt)
    print(p['status'])
if __name__=='__main__':main()
