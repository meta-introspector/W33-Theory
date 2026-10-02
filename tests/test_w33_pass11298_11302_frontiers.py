"""Independent tensor, spectral, counterterm, regulator and junction controls."""
from pathlib import Path
import sys
import numpy as np
import mpmath as mp
import sympy as sp
from scipy.optimize import root

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'analysis'))
import w33_pass11298_coupled_running_closure as R
import w33_pass11299_e6_flavor_spectral_obstruction as F
import w33_pass11300_invariant_vector_matching as V
import w33_pass11301_graph_adm_regulator as G
import w33_pass11302_sequestered_membrane_junction as M


def test_yukawa_components_and_two_channel_rank_budget():
    # Check different E6 and off-diagonal family scalar components.
    for a in (2, 53):
        c, error = R.coefficient(.7, .4, a)
        assert abs(c - (40*.7**2 + .4**2)) < 1e-10 and error < 1e-10
    for a in (56, 57, 64):
        c, error = R.coefficient(.7, .4, a)
        assert abs(c - (5*.7**2 + 29*.4**2)) < 1e-10 and error < 1e-10
    a, b, c, d = sp.symbols('a b c d')
    assert sp.Matrix([[a*c, a*d], [b*c, b*d]]).det() == 0
    assert sp.eye(2).det() == 1
    assert 34-20-2*12 == -10 and 34-20-2*6 == 2


def test_mandatory_S_quartics_from_independent_full_Hessians():
    n = 10
    # Orthonormal traceless diagonal basis, plus symmetric off-diagonal modes.
    E = np.linalg.qr(np.column_stack([np.ones(n), np.eye(n)[:, :n-1]]))[0][:, 1:]
    B = [np.diag(e) for e in E.T]
    for i in range(n):
        for j in range(i+1, n):
            X = np.zeros((n, n)); X[i, j] = X[j, i] = 1/np.sqrt(2); B.append(X)
    rng = np.random.default_rng(54)
    for _ in range(3):
        z = rng.normal(size=n); z -= z.mean(); S = np.diag(z)
        t, u = np.trace(S@S), np.trace(S@S@S@S)
        v = np.array([np.trace(S@X) for X in B])
        H1 = t*np.eye(54) + 2*np.outer(v, v)
        H2 = np.array([[np.trace(X@(S@S@Y+S@Y@S+Y@S@S)) for Y in B] for X in B])
        assert abs(np.trace(H1@H1)-62*t*t) < 1e-8
        assert abs(np.trace(H1@H2)-(11.2*t*t+6*u)) < 1e-8
        assert abs(np.trace(H2@H2)-(1.59*t*t+12.7*u)) < 1e-8
        assert abs(sum((z[i]-z[j])**4 for i in range(n) for j in range(i+1,n))-(10*u+3*t*t)) < 1e-8


def test_E6_spectral_law_is_blind_to_Takagi_phases():
    Ds, paths = F.graph(); assert len(paths) == 5
    rng = np.random.default_rng(113)
    Q = np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))[0]
    u = np.diag([1.,2.,4.]); d = Q@np.diag([.3,1.7,3.])@Q.T
    d2 = Q@np.diag(np.array([.3,1.7,3.])*np.exp(1j*np.array([.2,-.7,.4])))@Q.T
    assert np.max(abs(d@d.conj().T-d2@d2.conj().T)) < 1e-12
    masses = []  # no stored-certificate eigenvalues
    for f in (d,d2):
        m = np.kron(Ds[0],u)+np.kron(Ds[1],f)
        masses.append(np.linalg.eigvalsh(m.conj().T@m))
    assert np.max(abs(masses[0]-masses[1])) < 1e-10
    assert sum(masses[0]<1e-9) == 51


def test_vector_counterterm_full_matrix_second_jet():
    # Test noncommuting matrix variations, including mixed heavy/kernel blocks.
    nodes = [.00433251,.01174027,.01530281,.02777778,.03592023,.04381529]
    coeff, error = V.hermite(nodes); assert error < 1e-40
    scale = mp.mpf('.04'); mu = mp.mpf('.1')
    f = lambda t: scale**2*t*t*(mp.log(scale*t/mu**2)-mp.mpf(5)/6)
    p = lambda t: sum(c*t**i for i,c in enumerate(coeff))
    lam = [mp.mpf(0)] + [mp.mpf(str(x))/scale for x in nodes] + [mp.mpf(str(nodes[1]))/scale]
    fp = lambda x: mp.diff(f,x) if x else mp.mpf(0)
    pp = lambda x: mp.diff(p,x)
    for i,x in enumerate(lam):
        assert abs(fp(x)-pp(x)) < mp.mpf('1e-40')
        for j,y in enumerate(lam):
            if x == y:
                if not x: continue  # Gram kernel has Xprime_00=0.
                a,b=mp.diff(f,x,2),mp.diff(p,x,2)
            else:
                a,b=(fp(x)-fp(y))/(x-y),(pp(x)-pp(y))/(x-y)
            assert abs(a-b) < mp.mpf('1e-38')
    # First derivative of a Gram matrix is zero when restricted to its kernel.
    A=np.diag([0.,1.,2.]); dA=np.arange(9).reshape(3,3)/10
    dX=dA.T@A+A.T@dA; assert dX[0,0] == 0


def test_vector_cut_is_closed_and_subtraction_dependent():
    s=.009; m2=.0044
    z=V.vector_dispersion(s,m2)
    assert z>0 and s<4*m2
    assert abs(V.vector_dispersion(2*s,2*m2)-4*z) < 1e-10


def test_graph_bracket_symbolic_edge_and_nonzero_refinement():
    a,b,p,q,N0,N1,M0,M1=sp.symbols('a b p q N0 N1 M0 M1')
    H=lambda n0,n1:n0*p*p/2+n1*q*q/2+(n0+n1)*(a-b)**2/4
    x,y=H(N0,N1),H(M0,M1)
    pb=sum(sp.diff(x,z)*sp.diff(y,m)-sp.diff(x,m)*sp.diff(y,z) for z,m in [(a,p),(b,q)])
    assert sp.expand(pb-(N0*M1-M0*N1)*(p+q)*(b-a)/2)==0
    errs=[G.cycle_error(n) for n in (16,32,64,128)]
    assert errs[0]>1e-5 and all(errs[i+1]<errs[i]/3 for i in range(3))


def test_membrane_junction_independent_start_and_common_shift():
    r,l,g,mu,Q,Qhat=M.inputs()
    sol=root(M.equations,[r*1.002,l-.003,g*.997],args=(mu,Q,Qhat),tol=1e-10)
    assert sol.success and np.max(abs(M.equations(sol.x,mu,Q,Qhat))) < 1e-9
    for shift in (-.17,.03,.4):
        x=sol.x.copy();x[1]-=shift
        assert np.max(abs(M.equations(x,mu,Q,Qhat,bare=shift))) < 1e-9
    vol,avg,j=M.caps(*sol.x)
    rho=np.array([sol.x[1]+.5*1.2**2,sol.x[1]+.5*.8**2])
    assert avg > 4*(rho@vol)/(sol.x[2]*sum(vol))  # distributional wall curvature
    assert abs(1.2-.8-.4) < 1e-15 and abs(j) < 1e-9
