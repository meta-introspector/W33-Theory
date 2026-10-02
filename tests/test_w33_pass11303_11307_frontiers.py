"""Independent controls of the five model changes, including their open boundaries."""
from pathlib import Path
import sys,json
import numpy as np
import mpmath as mp
import sympy as s
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11303_scalar_source_and_uv_barrier as U
import w33_pass11304_singlet_phase_probes as F
import w33_pass11305_vector_matrix_and_total_jet_audit as V
import w33_pass11306_relational_graph_clock_constraints as G
import w33_pass11307_quantized_membrane_radial_mode as M

def cert(n):return json.loads(next((ROOT/'data').glob(f'w33_pass{n}_*.json')).read_text())

def test_single_scalar_can_match_two_sources_without_two_mediators():
    # The scalar propagator acts on a source, not a factorized pair of vertices.
    h,J1,J2,y,mu1,mu2,m=s.symbols('h J1 J2 y mu1 mu2 m',nonzero=True)
    source=mu1*J1+mu2*J2
    potential=m*m*h*h+2*h*source+source*source/(m*m)
    sol=s.solve(s.diff(potential,h),h)[0]
    assert s.expand(potential.subs(h,sol))==0
    assert s.expand(y*sol)==-y*mu1*J1/m**2-y*mu2*J2/m**2
    assert s.diag(mu1,mu2).det()!=0

def test_portal_independent_Riccati_certificate():
    x,y=s.symbols('x y');z=x+s.Rational(3,10)*y
    f=124*x*x+s.Rational(224,5)*x*y+s.Rational(159,50)*y*y-104*x+18
    g=24*x*y+s.Rational(127,5)*y*y-104*y+60
    low=s.Rational(16580,159)*z*z-104*z+36
    assert s.factor(f+s.Rational(3,10)*g-low)==(56*x-15*y)**2/159
    assert s.Rational(23,32)*s.Rational(1,10)+s.Rational(9,32)*s.Rational(73,90)==s.Rational(3,10)
    rng=np.random.default_rng(54)
    for _ in range(4):
        A=rng.normal(size=(54,54));A=(A+A.T)/2
        B=rng.normal(size=(54,7));C=rng.normal(size=(7,7));C=(C+C.T)/2
        H=np.block([[A,B],[B.T,C]])
        assert np.trace(H@H)>=np.trace(A@A)

def test_adjoint_encoding_in_a_second_gauge_realization():
    import w33_pass11285_semisimple_factor_higgs as H
    v,A,u,w,K,S=U.H.economical_reference();L=S@K
    X=np.zeros((10,10));X[1,8]=1;X[8,1]=-1;R=expm(.17*X)
    E=expm(.013j*H.F.H.setup()[0][0]);A=E@A@E.conj().T;u=E@u@R.T
    K=R@K@R.T;L=R@L@R.T;S2=-K@L
    assert max(np.max(abs(L+L.T)),np.max(abs(L@L+K@L+6*np.eye(10))),np.max(abs(A@u-1j*u@L)),abs(np.trace(K@L)))<1e-10
    assert np.max(abs(S2-S2.T))<1e-12 and abs(np.trace(S2))<1e-12
    assert np.max(abs(S2@S2+S2-6*np.eye(10)))<1e-10

def test_added_SO10_Weyls_support_a_bounded_AF_subsector_ray():
    d=cert(11303)['Yukawa_matter_escape'];rows=d['bounded_fixed_rays']
    r=next(x for x in rows if max(x['quartic_ratio_stability_eigenvalues'])<0)
    x,y,R=r['lambda1_over_g_squared'],r['lambda2_over_g_squared'],r['Yukawa_squared_over_g_squared']
    b=8/3;B=120-2*b-8*R
    assert abs((41/5)*R-27+b)<1e-12
    assert abs(124*x*x+44.8*x*y+3.18*y*y-B*x+18)<1e-10
    assert abs(24*x*y+25.4*y*y-B*y+60-8*R*R)<1e-10
    assert x+73/90*y>0 and x+y/10>0
    for t in (0.,10.,1e4):
        g2=.1/(1+b*.1*t/(8*np.pi**2))
        assert all(a>=0 for a in [g2,R*g2]) and x*g2+y*g2*73/90>0

def test_family_probe_mass_spectrum_is_covariant_and_phase_sensitive():
    u,d=F.fields(np.zeros(4));u2,d2=F.fields(np.array([.1,-.07,.13,.04]))
    assert max(np.max(abs(u@u.conj().T-u2@u2.conj().T)),np.max(abs(d@d.conj().T-d2@d2.conj().T)))<1e-12
    X=np.array([[0,.2+.1j,.17j],[-.2+.1j,0,.1], [.17j,-.1,0]])
    R=expm(X);assert np.max(abs(R.conj().T@R-np.eye(3)))<1e-12
    def spec(u,d):
        A=.2*(u+d).conj();mass=np.block([[A,3*np.eye(3)],[3*np.eye(3),np.zeros((3,3))]])
        return np.linalg.eigvalsh(mass.conj().T@mass)
    assert np.max(abs(spec(u,d)-spec(R@u@R.T,R@d@R.T)))<1e-10
    assert np.linalg.norm(spec(u,d)-spec(u2,d2))>1e-4
    h=1e-4;J=np.column_stack([(F.moments(np.eye(4)[i]*h)-F.moments(-np.eye(4)[i]*h))/(2*h) for i in range(4)])
    J/=np.linalg.norm(J,axis=1)[:,None];assert np.linalg.svd(J,compute_uv=False)[-1]>.01

def test_full_vector_cut_replays_radial_width_and_is_PSD():
    import w33_pass11290_gauged_family_decay as Old
    sm,O,light,cs,ev,kappa,d=V.vector_channels();t=.3
    VV=sum((V.vector_rho(c,t) for c in cs if c['kind']=='VV'),start=np.zeros((22,22)))
    expected=sum(Old.width(x,t)*np.sqrt(t) for x in ev)
    assert abs(VV[-1,-1]-expected)<1e-12
    full=sum((V.vector_rho(c,t) for c in cs),start=np.zeros((22,22)))
    assert min(np.linalg.eigvalsh(full))>-1e-11
    assert cert(11305)['Goldstone_VS_coupling_norm']>0
    assert sum(r['multiplicity'] for r in cert(11305)['rows'])==14

def test_total_scalar_IR_coefficient_is_not_a_Gram_kernel():
    H,C,M0,Y=V.light_data();w,O=np.linalg.eigh(H);P=O[:,w<1e-8]@O[:,w<1e-8].T
    CG=np.array([P@c@P for c in C]);IR=1e-8*np.einsum('aij,bji->ab',CG,CG)/(32*np.pi**2)
    d=cert(11305)['total_jet_audit'];assert d['IR_rank']==11
    assert np.max(abs(IR-np.array(d['Goldstone_IR_log_Hessian_coefficient'])))<1e-12
    assert np.linalg.norm(CG)>1 and all(r['jet_error']<1e-35 for r in d['regulated_counterfunctions'])

def test_relational_observable_survives_two_clock_orders():
    J,A,edges=G.native_generators();i,j=sorted(edges[0]);t=np.zeros(80);dt=.013;ds=.019
    _,B0=G.connection(t,[i,j],A,J);z=np.random.default_rng(6).normal(size=160)
    z1=expm(dt*B0[i])@z;t[i]=dt;_,B1=G.connection(t,[i,j],A,J);z1=expm(ds*B1[j])@z1;t[j]=ds
    U,_=G.connection(t,[i,j],A,J);assert np.max(abs((-J@U.T@J)@z1-z))<1e-11
    z2=expm(ds*A[j])@z;z2=expm(dt*A[i])@z2
    assert np.max(abs(z1-z2))<1e-11
    assert np.linalg.norm((expm(ds*A[j])@expm(dt*A[i])-expm(dt*A[i])@expm(ds*A[j]))@z)>1e-6
    assert cert(11306)['phase_space_count']['physical_real_phase_dimension']==160

def test_quantized_membrane_action_matches_flux_and_has_radial_mode():
    x,N,Qmem,Q,Qhat,mu,_=M.saddle();H=M.hessian(x,Qmem,Q,Qhat,mu)
    assert abs(mp.mpf('.4')*Qmem-2*mp.pi*N)<mp.mpf('1e-40')
    assert max(abs(v) for v in M.gradient(x,Qmem,Q,Qhat,mu))<mp.mpf('1e-35')
    curvature=H[0,0]-(H[0,1:4]*H[1:4,1:4]**-1*H[1:4,0])[0]
    # Independently re-solve the global conditions away from the stationary radius.
    def value(R):
        aux=mp.findroot(lambda L,G,f:tuple(M.gradient([R,L,G,f],Qmem,Q,Qhat,mu)[1:]),tuple(x[1:]),tol=mp.mpf('1e-40'))
        return M.action([R,*aux],Qmem,Q,Qhat,mu)
    h=mp.mpf('0.00001');finite=(value(x[0]+h)-2*value(x[0])+value(x[0]-h))/h**2
    assert curvature<0 and abs(finite-curvature)<mp.mpf('.0001')
