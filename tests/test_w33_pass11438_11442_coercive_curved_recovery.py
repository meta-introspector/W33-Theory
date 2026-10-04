"""Independent scope, geometry, RG and quantum-channel controls."""
import sys
import json
import hashlib
from pathlib import Path
from itertools import combinations
import numpy as np
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11438_11442_coercive_curved_recovery as N


def packet():
    return json.loads(N.OUT.read_text())


def test_source_certificates():
    for name,digest in packet()['source_sha256'].items():
        canonical=json.dumps(json.loads((ROOT/name).read_text()),sort_keys=True,separators=(',',':')).encode()
        assert hashlib.sha256(canonical).hexdigest()==digest


def test_runaway_exact_homogeneity_and_second_scaling():
    p=N.prior()['native_alignment'];phi=N.M.P.read(p['phi']);psi=N.M.P.read(p['psi'])
    for t in [.3,.7,1.4]:
        assert abs(N.M.cross_pin(t*phi,psi)[3]*t**9-p['CP_flux'])<1e-10
        assert abs(N.M.cross_pin(phi,t*psi)[3]/t**9-p['CP_flux'])<1e-10
    rows=packet()['coercive_alignment']['scans']
    assert rows[-1]['old_total']<-100 and rows[-1]['native']>-.4


def test_bounded_cross_random_fields_and_origin():
    rng=np.random.default_rng(11438)
    for scales in [(0,0),(1e-8,1e4),(1,1),(100,100)]:
        fields=[s*(rng.normal(size=(27,3))+1j*rng.normal(size=(27,3))) for s in scales]
        C,q=N.bounded_cross(*fields)
        assert np.linalg.norm(C,2)<=1+1e-12 and abs(q)<=24
    assert packet()['coercive_alignment']['relative_orbit_dimension']==70


def test_bounded_map_covariance_and_CP():
    p=N.prior()['native_alignment'];phi=N.M.P.read(p['phi']);psi=N.M.P.read(p['psi'])
    h=np.array([[.1,.2j,.3],[-.2j,-.2,.1],[.3,.1,.1]])
    U=expm(1j*h);C,q=N.bounded_cross(phi,psi)
    transformed,qp=N.bounded_cross(phi@U.T,psi@U.T)
    assert np.linalg.norm(transformed-U@C@U.conj().T)<1e-13
    assert abs(q-qp)<1e-14 and abs(q+N.bounded_cross(phi.conj(),psi.conj())[1])<1e-14


def test_finite324_native_action_gradient_covariance_and_hessian():
    import w33_pass11438_finite_native_model as F
    p=packet()['coercive_alignment']['finite_field_witness'];x=np.array(p['coordinates'])
    H=np.array(p['full324_Hessian']);rng=np.random.default_rng(11443)
    with F.native_context():
        value,g=F.fun(x)
        assert abs(value-p['action_value'])<1e-11 and np.linalg.norm(g)<5e-6
        assert value < -5/8 and p['global_CP_flux_existence']
        assert np.linalg.norm(g-np.array(p['gradient']))<1e-11
        for _ in range(3):
            d=rng.normal(size=324);d/=np.linalg.norm(d);h=1e-6
            difference=(F.fun(x+h*d)[0]-F.fun(x-h*d)[0])/(2*h)
            assert abs(difference-g@d)<1e-7
        phi=x[:81]+1j*x[81:162];psi=x[162:243]+1j*x[243:]
        U=expm(.03j*np.einsum('a,aij->ij',rng.normal(size=86),F.H))
        assert abs(F.fun(F.pack(U@phi,U@psi))[0]-value)<1e-10
        assert abs(F.fun(F.pack(phi.conj(),psi.conj()))[0]-value)<1e-10
        # The full matrix is independently replayed along random directions.
        for _ in range(2):
            d=rng.normal(size=324);d/=np.linalg.norm(d);h=1e-5
            assert np.linalg.norm((F.fun(x+h*d)[1]-F.fun(x-h*d)[1])/(2*h)-H@d)<1e-4
    ev=np.array(p['normal_eigenvalues'])
    assert len(ev)==238 and sum(ev>1e-3)==234 and max(abs(ev[:4]))<1e-5
    assert p['gauge_orbit_rank']==86 and p['gauge_Hessian_residual']<1e-7
    assert p['CP_flux']>1e-4 and abs(p['I6_phi'][1])>.1


def test_flux_periodicity_all_charges_and_size_boundary():
    p=packet()['admissible_measure'];L=p['minimum_L_for_unit_flux'];U0,U1=N.uniform_flux_links(L)
    for q in p['one_family_charges']:
        a,b=U0**q,U1**q
        plaquette=a*np.roll(b,-1,axis=0)*np.roll(a.conj(),-1,axis=1)*b.conj()
        assert np.max(abs(plaquette-np.exp(2j*np.pi*q/L**2)))<1e-12
        assert np.max(abs(np.angle(plaquette)))<1/30
    assert 6*2*np.pi/(L-1)**2>=1/30
    assert p['odd_absolute_charge_multiplicities']=={'1':6,'3':2}


def test_weak_patch_admissibility_under_gauge_transform():
    p=N.prior()['local_chiral_measure'];A=1e-4*np.array(p['background']).reshape((3,)*4+(4,))
    alpha=np.random.default_rng(11439).normal(size=(3,)*4)*.01
    B=A.copy()
    for mu in range(4):B[...,mu]+=alpha-np.roll(alpha,-1,axis=mu)
    assert np.max(abs(N.plaquettes(A,3)-N.plaquettes(B,3)))<1e-16


def spherical_area(lengths,tri,k):
    # Independently use spherical cosine law, rather than inverse-Gram dihedrals.
    a,b,c=[np.sqrt(k)*lengths[tri[i],tri[j]] for i,j in [(1,2),(0,2),(0,1)]]
    angles=[np.arccos(np.clip((np.cos(x)-np.cos(y)*np.cos(z))/(np.sin(y)*np.sin(z)),-1,1))
            for x,y,z in [(a,b,c),(b,a,c),(c,a,b)]]
    return (sum(angles)-np.pi)/k


def test_all_radial_equations_with_independent_area_derivatives():
    p=packet()['curved_stationarity'];bd=np.array(p['boundary_lengths']);k=p['sectional_curvature']
    for scan in p['scans']:
        L,fine,deficits=N.spherical_refinement(bd,k,np.array(scan['weights']))
        co=N.M.refined_geometry(bd)[0]
        gradients=[];h=1e-5
        for i,s in enumerate(co):
            for v in s:
                plus=L.copy();minus=L.copy();plus[6+i,v]=plus[v,6+i]=L[6+i,v]+h
                minus[6+i,v]=minus[v,6+i]=L[6+i,v]-h
                g=sum(d*(spherical_area(plus,t,k)-spherical_area(minus,t,k))/(2*h)
                      for t,d in deficits.items() if 6+i in t and v in t)
                gradients.append(g)
        assert np.linalg.norm(gradients)<1e-7
        assert max(abs(d) for d in deficits.values())<1e-9


def test_spherical_subdivision_is_coordinate_independent():
    p=packet()['curved_stationarity'];bd=np.array(p['boundary_lengths']);k=p['sectional_curvature']
    rng=np.random.default_rng(11440);Q=np.linalg.qr(rng.normal(size=(5,5)))[0]
    for i,s in enumerate(N.M.refined_geometry(bd)[0]):
        G=np.cos(np.sqrt(k)*bd[np.ix_(s,s)]);X=np.linalg.cholesky(G);w=np.array(p['scans'][1]['weights'][i])
        c=w@X;c/=np.linalg.norm(c);Y=X@Q;d=w@Y;d/=np.linalg.norm(d)
        assert np.linalg.norm(X@c-Y@d)<1e-12
        assert np.linalg.norm(X@X.T-G)<1e-12


def test_renormalized_force_curvature_and_scale():
    p=packet()['renormalized_slice']
    assert abs(p['residual_force'])<.01 and abs(p['measured_radial_curvature']-100)<.05
    assert max(r['scale_spread'] for r in p['scans'])<1e-7
    assert max(abs(r['supertrace_polynomial_error']) for r in p['scans'])<1e-4
    # c0 is an input: changing it changes the vacuum energy without changing force.
    c=np.array(p['counterterm_coefficients']);phi=p['imposed_stationary_phi']
    assert abs(np.array([1,phi**2,phi**4])@np.array([3.,0,0])-3)<1e-12
    assert abs(c[0])>1


def test_curved_complete_hessian_and_center_null_space():
    p=packet()['curved_stationarity'];bd=np.array(p['boundary_lengths']);k=p['sectional_curvature']
    controls=p['full15_radial_Hessian_controls'];small=[max(abs(np.array(c['eigenvalues'])[:12])) for c in controls]
    assert 20<small[0]/small[1]<30
    H=np.array(controls[-1]['Hessian']);weights=np.ones((3,5))/5;h=1e-6;T=[]
    coarse=N.M.refined_geometry(bd)[0]
    def radial(W):
        L,_,_=N.spherical_refinement(bd,k,W)
        return np.concatenate([L[6+i,list(s)] for i,s in enumerate(coarse)])
    for star in range(3):
        for j in range(4):
            shift=np.zeros((3,5));shift[star,j]=h;shift[star,4]=-h
            T.append((radial(weights+shift)-radial(weights-shift))/(2*h))
    T=np.array(T).T
    assert np.linalg.matrix_rank(T,tol=1e-6)==12
    assert np.linalg.norm(H@T)/(np.linalg.norm(H)*np.linalg.norm(T))<1e-6
    assert min(np.linalg.eigvalsh(H)[-3:])>1000


def test_exact_erasure_KL_and_nonselective_density_recovery():
    V=N.steane_frame();rho=np.array([[.4,.2+.1j],[.2-.1j,.6]])
    for sites in [(0,),(4,),(0,6),(2,5)]:
        Ks=N.reset_kraus(sites);r=len(Ks);out=np.zeros((128,128),complex)
        for i,K in enumerate(Ks):
            for j,L in enumerate(Ks):
                assert np.linalg.norm(V.conj().T@K.conj().T@L@V-(np.eye(2)/r if i==j else np.zeros((2,2))))<1e-12
        noisy=sum(K@V@rho@V.conj().T@K.conj().T for K in Ks)
        Rs=[np.sqrt(r)*V@(K@V).conj().T for K in Ks]
        out=sum(R@noisy@R.conj().T for R in Rs)
        assert np.linalg.norm(out-V@rho@V.conj().T)<1e-12
    assert len(packet()['leakage_recovery']['exact_recovery_supports'])==28
    assert len(packet()['leakage_recovery']['uncorrectable_three_erasure_supports'])==7


def test_routes_and_swap_spectator_restoration():
    p=packet()['leakage_recovery'];A=N.M.P.load()['A']
    for path,cost in zip(p['native_graph_paths'],p['CNOT_cost_per_routed_syndrome_or_flag_pair']):
        assert all(A[u,v] for u,v in zip(path,path[1:])) and cost==6*(len(path)-2)+1
        # Exhaust all classical basis states on the actual path; SWAP routing
        # must restore every spectator and perform exactly the remote CNOT.
        n=len(path)
        for bits in range(1<<n):
            x=[(bits>>i)&1 for i in range(n)];initial=x.copy()
            for i in range(n-2):x[i],x[i+1]=x[i+1],x[i]
            x[-1]^=x[-2]
            for i in reversed(range(n-2)):x[i],x[i+1]=x[i+1],x[i]
            initial[-1]^=initial[0]
            assert x==initial


if __name__=='__main__':
    for name,fn in list(globals().items()):
        if name.startswith('test_'):
            fn();print(name,'PASS',flush=True)
