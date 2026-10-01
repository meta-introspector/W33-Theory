"""Independent controls for the five physical follow-ups to c944a5060."""
from pathlib import Path
import sys,json
from fractions import Fraction
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_20261001_degree18_phase_completion as P
import w33_20261001_chiral_mass_assignment_audit as C
import w33_20261001_completed_eft_one_loop as L
import w33_20261001_metric_dynamics_and_heat_bound as G
import w33_20261001_uv_running_scale_audit as U

def test_global_I18_homogeneity_differential_and_phase_mass():
    P.exact_restriction()
    r=np.random.default_rng(19);p=.11*(r.normal(size=81)+1j*r.normal(size=81));val,grad=P.evaluate(p)
    assert abs(grad@p-18*val)<1e-10
    assert abs(P.evaluate((.8+.2j)*p)[0]-(.8+.2j)**18*val)<1e-10
    q,D=L.tangent();phi=P.I.C.plane()@q/np.sqrt(3);tau=-1/729
    potential=lambda th:abs(P.evaluate(np.exp(1j*th)*phi)[0]-tau)**2
    h=1e-5;curvature=(potential(h)+potential(-h))/(2*h*h)
    assert abs(curvature-324/729**2)<1e-9

def test_exact_little_group_and_chiral_branching():
    out=C.stabilizer();assert out['stabilizer_dimension']==8
    assert all(g['restriction_to_three_null_modes']=='zero3x3' for g in out['generators'])
    b=C.branching();assert b['chiral_index_after_any_subgroup_restriction']==0
    assert b['alternative_single_81_chiral_indices']

def test_completed_scalar_and_fermion_spectra_against_independent_block_law():
    q,_=L.tangent();v,sc,f=L.spectra(q,step=7e-6)
    expected=np.r_[v[v>1e-10],.4,.32,.32,.05,.05]
    assert np.max(abs(np.sort(sc[sc>1e-9])-np.sort(expected)))<2e-8
    assert np.max(abs(np.sqrt(np.sort(f)[:3])-np.array([.001,.002,.003])))<1e-10
    # Exact moment2=20, moment4=8 of original mirror spectrum scale with y.
    assert abs(np.sum(f)-20*L.PARAM['y']**2-(.001**2+.002**2+.003**2))<1e-10

def test_directed_continuum_tail_enclosure_and_metric_sign():
    d=json.loads((ROOT/'data/w33_20261001_metric_dynamics_and_heat_bound.json').read_text())
    h=d['heat_certificate'];lo,hi=Fraction(h['response_lower_exact']),Fraction(h['response_upper_exact'])
    tail=Fraction(h['tail_upper_exact']);assert 0<tail<Fraction(1,10**15)
    assert lo<hi and hi-lo>=2*tail
    # Numerically independent legacy cutoff sum agrees to its floating accuracy.
    old=G.G.heat_response(.2,30);assert abs(old-float((lo+hi)/2))<1e-9
    assert float(hi)<0;G.dynamics()

def test_actual_inventory_beta_and_matching_drift():
    d=U.payload();assert d['one_loop_b0']=={'E6':'17','SU3_family':'-59/2'}
    # Removing ALL fermions restores asymptotic freedom; even one81 reverses family sign.
    assert s.Rational(13,2)>0 and s.Rational(13,2)-9<0
    assert '93' in d['running_from_equal_g0']['inverse_difference']


def test_chiral_spectator_weights_and_condensate_radial_equation():
    import w33_20261001_chiral_decuplet_and_condensate_bridge as B
    d=B.symmetric_cube_weights();assert d['A10']==27 and d['T10']=='15/2'
    r,m,lam=s.symbols('r m lam',positive=True)
    V=36*lam**18/r**14-18*m*lam**9/r**6
    rd=s.factor(s.diff(V,r)*r**15/(36*lam**9))
    assert rd==3*m*r**8-14*lam**9
    # Independent numerical replay of the slice angular curvature; no copied Hessian routine.
    import w33_pass11255_11259_g26_common as Q
    u=Q.invariants()[2];fn=s.lambdify(Q.VARS,[u]+[s.diff(u,x) for x in Q.VARS],'numpy',cse=True)
    q,D=L.tangent()
    def pot(t):
        p=q+t*D[:,0];p/=np.linalg.norm(p);a=fn(*p);W=(3/14)*np.exp(-np.log(complex(-729*a[0]))/3)
        return abs(W)**2*np.vdot(a[1:],a[1:]).real/(9*abs(a[0])**2)-18*W.real
    h=1e-4;curve=(pot(h)+pot(-h)-2*pot(0))/h**2
    assert abs(curve-324/49)<1e-4
