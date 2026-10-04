"""Independent controls on claimed scope, channels, geometry and running."""
import sys,json,hashlib
from pathlib import Path
from itertools import combinations,product
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11443_11447_preserved_gauge_lorentz_recovery as P


def packet():return json.loads(P.OUT.read_text())


def test_source_hashes():
    for name,h in packet()['source_sha256'].items():
        raw=json.dumps(json.loads((ROOT/name).read_text()),sort_keys=True,separators=(',',':')).encode()
        assert hashlib.sha256(raw).hexdigest()==h


def test_stationary_preserved_gauge_and_second_group_basis():
    import w33_pass11438_finite_native_model as F
    p=packet()['preserved_gauge'];x=np.array(p['coordinates']);a=x[:81]+1j*x[81:162];b=x[162:243]+1j*x[243:]
    K,B=P.gauge_fixed_space();rng=np.random.default_rng(11443);R=np.linalg.qr(rng.normal(size=(8,8)))[0];K2=np.einsum('ab,bij->aij',R,K)
    assert np.linalg.norm(K2@a)<1e-9 and np.linalg.norm(K2@b)<1e-9
    with F.native_context():
        v,g=F.fun(x)
    assert abs(v-p['value'])<1e-10 and np.linalg.norm(g)<5e-6
    assert v < -2.94 and p['common_gauge_rank']==78
    assert p['negative_normal_count']==0 # Resolution threshold1e-3; soft modes unresolved.
    assert len(p['normal_eigenvalues'])==246


def test_relaxation_invalidates_line_stiffness_inference():
    rows=packet()['relaxed_soft'];assert len(rows)==8
    assert all(abs(r['relaxed_over_h4'])<.005 and r['raw_over_h4']>5 for r in rows)


def test_measure_cocycle_is_insufficient_and_flux_is_nonzero():
    p=packet()['measure'];assert p['orbit_cocycle_error']<1e-11
    assert p['negative_control']['linear']==0 and p['negative_control']['cubic']!=0
    for r in p['charged_fluxes']:
        assert abs(abs(r['flux'])-r['charge'])<1e-10 and r['max_plaquette']<1/30


def test_lorentzian_elimination_and_off_shell_boundary_fields():
    X=np.array(packet()['lorentzian']['coordinates']);rng=np.random.default_rng(11445);f=rng.normal(size=5)
    K,v=P.scalar_matrix(X,0);weights=np.array([.1,.15,.2,.25,.3]);S,vol,pivot=P.scalar_refinement(X,weights,0)
    assert abs(f@(S-K)@f)<1e-10 and abs(vol-v)<1e-10
    assert sum(np.linalg.eigvalsh(K)<-1e-9)==1
    Sm,_,_=P.scalar_refinement(X,weights,.2);Km,_=P.scalar_matrix(X,.2)
    assert np.linalg.norm(Sm-Km)>1e-6


def test_running_normalization_and_pole():
    p=packet()['running'];assert p['heavy_triplet_species']==966 and p['total_triplet_species']==972
    assert p['high_energy_b0']==11-2*972/3==-637
    for r in p['examples']:
        assert abs(1/r['supplied_g']**2-637*r['log_Landau_ratio']/(8*np.pi**2))<1e-11


def pauli(xmask,zmask):
    j=np.arange(128);A=np.zeros((128,128),complex);A[j^xmask,j]=1-2*np.array([(int(k&zmask).bit_count()%2) for k in j]);return A


def test_actual_steane_erasure_circuit_branches():
    V=P.N.steane_frame();rng=np.random.default_rng(11447);a=rng.normal(size=2)+1j*rng.normal(size=2);a/=np.linalg.norm(a);psi=V@a
    # Each mask's complete Pauli erasure channel has orthogonal syndrome images.
    # The operational syndrome projector is E V Vdag E† and correction E†.
    for sites in combinations(range(7),2):
        branches=[]
        for labels in product(range(4),repeat=2):
            xm=sum((label&1)<<site for site,label in zip(sites,labels));zm=sum((label>>1)<<site for site,label in zip(sites,labels));E=pauli(xm,zm);branches.append(E)
        images=np.hstack([E@V for E in branches]);assert np.linalg.norm(images.conj().T@images-np.eye(32))<1e-11
        recovered=sum(E.conj().T@(E@V)@((E@V).conj().T@(E@psi))/16 for E in branches)
        assert np.linalg.norm(recovered-psi)<1e-11
    sites=(0,3);rho=np.outer(psi,psi.conj());reset=sum(K@rho@K.conj().T for K in P.N.reset_kraus(sites))
    twirled=sum(E@reset@E.conj().T for E in [pauli(sum((l&1)<<s for s,l in zip(sites,labels)),sum((l>>1)<<s for s,l in zip(sites,labels))) for labels in product(range(4),repeat=2)])/16
    erased=sum(E@rho@E.conj().T for E in [pauli(sum((l&1)<<s for s,l in zip(sites,labels)),sum((l>>1)<<s for s,l in zip(sites,labels))) for labels in product(range(4),repeat=2)])/16
    assert np.linalg.norm(twirled-erased)<1e-11


def test_clock_bounds_and_optimal_route():
    p=packet()['recovery'];clock=p['native_clock'];assert clock['joint_dark_dimension']==2 and clock['logical_clock_residual']<1e-12
    for row in clock['zero_phase_detection']:assert row['max_false_accept']<=row['envelope']+1e-14
    r=p['route'];assert r['baseline_routed_CNOTs']==r['lower_bound']==288 and r['conditional_routed_CNOTs']==372
    assert len(set(r['data']+[r['hub'],r['flag']]))==9
