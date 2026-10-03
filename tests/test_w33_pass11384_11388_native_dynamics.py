"""Independent native invariant, hard-class, graph, bulk and frame controls."""
from pathlib import Path
import sys
import json
import numpy as np
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11384_11388_native_dynamics as m


def test_native_cp_minimum_is_not_the_old_gauge_real_T_orbit():
    d=m.native_cp_vacuum()
    assert d['positive_normal_modes']==84 and d['gauge_orbit_rank']==78
    assert min(d['normal_eigenvalues'])>.02
    x,y=s.symbols('x y',real=True)
    V=(x*x+y*y)**2-(x*x+y*y)+x*x+x/2
    p={x:-s.Rational(1,4),y:s.sqrt(7)/4}
    assert V.subs(p)==-s.Rational(5,16)
    assert all(s.simplify(s.diff(V,z).subs(p))==0 for z in [x,y])
    H=s.hessian(V,[x,y]).subs(p)
    assert H.det()>0 and H[0,0]>0
    assert s.im((-1+s.I*s.sqrt(7))/4)!=0
    # The remaining common-phase Z6 leaves I6 fixed; it cannot identify the
    # two nonreal invariant values as the old real-I18 T selector did.
    assert s.simplify(s.exp(6*s.I*s.pi/3))==1


def test_compact_inputs_and_analytic_native_gradient_are_self_contained():
    from scipy.linalg import expm
    import w33_20261001_degree18_phase_completion as D
    B,t=m.native_tensors();H=m.compact_generators(B)
    data=json.loads((ROOT/'data/w33_pass11384_11388_native_dynamics.json').read_text())['sections']['11384_native_cp']
    phi=np.array(data['full_field_real'])+1j*np.array(data['full_field_imag'])
    old=D.I.tensors;D.I.tensors=lambda:t
    try:
        a,b,ga,gb=D.I.evaluate(phi);c,gc=D.evaluate(phi)
        rotation=expm(.03j*H[5]+.04j*H[-2]);rotated=rotation@phi
        ar,br,_,_=D.I.evaluate(rotated);cr,_=D.evaluate(rotated)
        assert max(abs(np.array([ar,br,cr])-np.array([a,b,c])))<1e-8
        rng=np.random.default_rng(81);v=rng.normal(size=81)+1j*rng.normal(size=81);v/=np.linalg.norm(v)
        h=2e-6
        ap,bp,_,_=D.I.evaluate(phi+h*v);am,bm,_,_=D.I.evaluate(phi-h*v)
        assert abs((ap-am)/(2*h)-ga@v)<1e-6
        assert abs((bp-bm)/(2*h)-gb@v)<1e-5
    finally:D.I.tensors=old


def test_hard_matching_classes_do_not_conflate_algebra_and_commutant():
    d=m.hard_matching_classes()
    assert (d['unrestricted_symmetric_normal_dimension'],d['central_projector_preserving_dimension'],
            d['matrices_commuting_with_cut_algebra_symmetric_dimension'],d['invariant_under_full_orthogonal_commutant_dimension'])==(105,36,29,7)
    # In the doubled irrep, preserving the protecting multiplicity rotations
    # means A tensor I, whereas commuting with every cut means I tensor B.
    A=np.array([[1.,2.],[2.,3.]]); B=np.array([[2.,0.],[0.,-1.]])
    assert np.linalg.norm(np.kron(A,np.eye(2))@np.kron(np.eye(2),B)-np.kron(np.eye(2),B)@np.kron(A,np.eye(2)))==0
    assert d['central_breaking_residual']>.1


def test_hub_graph_is_a_tree_without_erasing_native_matter_cycles():
    d=m.hub_gravity_architecture()
    assert d['spin2_interactions']==d['spin2_fields']-1==80
    assert d['conditional_spin2_polarizations']==402
    assert sum(d['native_internal_matter_spectrum'].values())==80
    assert 'g0 only' in d['matter_action']
    assert 'tree' in d['boundary']


def test_bulk_relaxed_wall_has_an_exact_negative_direction():
    d=m.scalar_wall_relaxation()
    H=s.Matrix([[s.Rational(x) for x in row] for row in d['Hessian']])
    assert H.det()==-s.Rational(13312,3675)
    direction=s.Matrix([1,s.Rational(5,4)])
    assert (direction.T*H*direction)[0]==-s.Rational(100,147)
    assert min(d['bulk_relaxed_Hessian_eigenvalues'])<-.28
    # This is not the old one-junction frozen-bulk displacement: q, C and c
    # change to keep flux, pole regularity and both bulk Ricci equations.
    u,v=s.symbols('u v',positive=True)
    S,q,c,C,F=m.relaxed_wall_action(u,v)
    assert s.simplify(q/c*(1/u-1/v)-1)==0
    ff=s.lambdify((u,v),(q,c,F),'numpy')
    for L in [.999,1.001]:
        qv,cv,fw=ff(L,5*L/4)
        assert cv>0 and fw>0 and abs(qv/cv*(1/L-4/(5*L))-1)<1e-12


def test_radial_vector_positivity_does_not_imply_scalar_wall_stability():
    import w33_pass11379_11383_full_frontier as prior
    assert prior.wall_factorization()['status']=='PASS'
    assert prior.fixed_action_scale()['second_derivative_at_L1']=='4'
    d=m.scalar_wall_relaxation()
    assert s.Rational(d['negative_direction_curvature'])<0


def test_relaxed_action_matches_direct_four_dimensional_quadrature():
    from scipy.integrate import quad
    u,v=s.symbols('u v',positive=True)
    action,q,c,C,_=m.relaxed_wall_action(u,v)
    evaluate=s.lambdify((u,v),(action,q,c,C),'numpy')
    for pole,wall in [(1.,1.25),(.999,1.25),(1.001,1.251)]:
        expected,qv,cv,constant=evaluate(pole,wall)
        F=lambda b:constant/b-qv*qv/(2*b*b)-16/45
        Fp=lambda b:-constant/b**2+qv*qv/b**3
        Fpp=lambda b:2*constant/b**3-3*qv*qv/b**4
        R=lambda b:-Fpp(b)-4*Fp(b)/b-2*F(b)/b**2
        bulk=quad(lambda b:2*b*b/cv*(-R(b)/2+qv*qv/(2*b**4)+(16/45)/b**2),pole,wall)[0]
        area=np.sqrt(F(wall))*wall*wall/cv
        extrinsic=np.sqrt(F(wall))*(Fp(wall)/(2*F(wall))+2/wall)
        direct=bulk-2*area*extrinsic+(64/75)*area
        assert abs(direct-expected)<1e-12


def test_induced_gravity_requires_a_different_stationarity_condition():
    d=m.native_dimensional_transmutation()
    assert d['trace_L_squared']==1600
    z,B,A,xi=s.symbols('z B A xi',positive=True)
    V=z**4*(A+B*s.log(z*z))
    stationary_log=-s.Rational(1,2)-A/B
    assert s.simplify(s.diff(V,z).subs(s.log(z*z),stationary_log))==0
    VE=V/(xi*xi*z**4)
    assert s.simplify(s.diff(VE,z))==2*B/(xi*xi*z)
    assert s.diff(VE,z).is_positive


def test_packet_records_five_executed_directions_with_boundaries():
    d=json.loads((ROOT/'data/w33_pass11384_11388_native_dynamics.json').read_text())
    assert d['status']=='PASS' and len(d['sections'])==5
    assert all(v['status']=='PASS' and v['boundary'] for v in d['sections'].values())
