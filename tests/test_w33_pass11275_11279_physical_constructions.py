"""Second-realization, propagator and cochain regression controls."""
from pathlib import Path
import sys,json
import numpy as np
import sympy as s
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11275_polynomial_sm_higgs as H
import w33_pass11276_sextet_composite_spectator as S
import w33_pass11277_radial_momentum_pole as P
import w33_pass11278_symplectic_ring_geometry as G
import w33_pass11279_history_flux_obstruction as F

def cert(n,name):return json.loads((ROOT/'data'/f'w33_pass{n}_{name}.json').read_text())

def test_polynomial_gauge_invariance_and_exact_normal_rank():
    J,T=H.jacobians();assert H.rank_mod(J,103)==120 and H.rank_mod(T,103)==66
    assert np.array_equal(J@T,np.zeros((len(J),78),int))
    basis,v,A,d=H.setup();N=d[:,:,0].T@d[:,:,1]
    assert H.rank_mod(N,103)==5 and np.array_equal((A@A+A-6*np.eye(27))@N,np.zeros((27,27)))
    U=expm(.11j*basis[10]+.17j*basis[23])
    assert H.potential(U@v,U@A@U.conj().T)<1e-14
    delta=np.zeros_like(v);delta[0,0]=.001
    assert H.potential(v+delta,A)>1e-6
    assert all(np.max(abs(z))<1e-12 for z in H.lifted_constraints(U@v,U@A@U.conj().T,U@N@U.conj().T,U@A@N@U.conj().T))
    # Perturbing an independently added field violates its defining constraint.
    X=N.astype(complex);X[0,0]+=.01
    assert np.linalg.norm(H.lifted_constraints(v,A,X,A@N)[0])>.009

def test_sextet_nonsingular_congruence_and_singular_control():
    Q=s.Matrix([[1,1,0],[0,1,1],[1,0,2]]);F=Q*Q.T
    assert F.det()!=0 and S.mass(F).rank()==10 and S.mass(F).is_positive_definite
    assert S.mass(s.diag(1,1,0)).rank()<10
    x=cert(11276,'sextet_composite_spectator')
    assert s.Rational(x['beta_family_new'])==-s.Rational(79,3)
    assert x['beta_E6_with_two_complex27_and_real78']=='28'

def test_bubble_analytic_value_derivative_and_pole():
    # Direct integration vs closed form below threshold, independently derived.
    for z in (.01,.5,1,2,3):
        closed=2-2*np.sqrt(4/z-1)*np.arctan(1/np.sqrt(4/z-1))
        assert abs(P.bubble(z)-closed)<1e-12
    assert abs(P.bubble(1e-5)/1e-5-1/6)<2e-7
    c=cert(11277,'radial_momentum_pole');M=c['curvature_mass_squared'];a=c['bubble_prefactor'];p=c['one_loop_Dyson_root_squared']
    assert abs(p-M+a*P.bubble(p/M))<1e-13
    assert abs(c['one_loop_pole_squared_strict']-(M-a*P.bubble(1)))<1e-14
    # This belongs to the declared scalar truncation; the full EFT remains excluded.
    assert 'Ten other' in c['excluded'][0] and 0<c['residue_at_Dyson_root']<1

def test_ring_projective_normalization_and_neighbor_fibers():
    p=G.points(9);first=p[0]
    same=[v for v in p if np.all((v-first)%3==0)]
    assert len(same)==27
    assert sum(int(first@G.J@v)%9==0 for v in same)==9
    assert {tuple(2*v%9) for v in p}!=set(map(tuple,p)) # coordinates differ; rays agree.
    def normal(v):
        i=next(i for i,x in enumerate(v) if x%3);return tuple(v*pow(int(v[i]),-1,9)%9)
    assert {normal(2*v%9) for v in p}==set(map(tuple,p))

def test_history_boundary_rank_and_action_variation():
    _,_,_,D=F.complex_data();assert np.array_equal(D.T@D,4*np.eye(40,dtype=int))
    # All circle periods retain injective D3 tensor I; not an accident of3ticks.
    for ticks in (2,4):
        dc=np.zeros((ticks,ticks),int)
        for i in range(ticks):dc[i,i]=-1;dc[(i+1)%ticks,i]=1
        top=np.vstack([np.kron(D,np.eye(ticks,dtype=int)),-np.kron(np.eye(40,dtype=int),dc)])
        assert H.rank_mod(top,103)==40*ticks
        rng=np.random.default_rng(ticks);a=rng.integers(-3,4,size=top.shape[0]);lam=rng.integers(-3,4,size=top.shape[1])
        assert lam@(top.T@a)==(top@lam)@a
        assert np.any(top@np.ones(40*ticks,dtype=int))
