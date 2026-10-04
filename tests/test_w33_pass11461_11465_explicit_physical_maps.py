import sys,json,hashlib
from pathlib import Path
from itertools import product
import numpy as np
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11461_11465_explicit_physical_maps as M


def packet():return json.loads(M.OUT.read_text())


def test_source_binding():
    for n,h in packet()['source_sha256'].items():
        assert hashlib.sha256(json.dumps(json.loads((ROOT/n).read_text()),sort_keys=True,separators=(',',':')).encode()).hexdigest()==h


def test_native_charge_algebra_and_anomalies():
    p=packet()['charges'];K=M.dec(p['colour']);L=M.dec(p['left']);R=M.dec(p['right']);W=M.dec(p['weak']);Y=M.dec(p['hypercharge']);CW=M.dec(p['weak_Casimir'])
    for A,B in [(K,L),(K,R),(L,R)]:
        assert max(np.linalg.norm(a@b-b@a) for a in A for b in B)<1e-9
    for A in [K,L,R,W]:
        flat=A.reshape(len(A),-1).T
        for a in A:
            for b in A:
                v=(1j*(a@b-b@a)).ravel();assert np.linalg.norm(v-flat@np.linalg.lstsq(flat,v,rcond=None)[0])<1e-9
    assert np.linalg.norm(CW-sum(a@a for a in W))<1e-10
    assert max(np.linalg.norm(Y@a-a@Y) for a in list(K)+list(W))<1e-9
    assert abs(np.trace(Y))<1e-10 and abs(np.trace(Y@Y@Y))<1e-10
    assert abs(np.trace(Y@sum(a@a for a in K)))<1e-10 and abs(np.trace(Y@CW))<1e-10
    assert abs(np.trace(Y+np.eye(27)))>20 # Deliberately bad charge shift.


def test_casimir_branching_survives_another_basis():
    p=packet()['charges'];Y=M.dec(p['hypercharge']);C=sum(a@a for a in M.dec(p['colour']));W=M.dec(p['weak_Casimir'])
    rng=np.random.default_rng(11461);U=np.linalg.qr(rng.normal(size=(27,27))+1j*rng.normal(size=(27,27)))[0]
    A=3*C+5*W+Y;assert np.max(abs(np.linalg.eigvalsh(A)-np.linalg.eigvalsh(U.conj().T@A@U)))<1e-10
    assert sum(r['states'] for r in p['branching'])==27
    assert {(tuple(r['key']),r['states']) for r in p['branching']}=={((1,1,1),6),((1,0,-2),3),((1,0,-4),3),((1,0,2),6),((0,1,3),2),((0,1,-3),4),((0,0,6),1),((0,0,0),2)}


def test_photon_obstruction_on_actual_candidate():
    p=packet()['charges'];x=np.array(M.N.previous()['preserved_gauge']['coordinates']);phi=(x[:81]+1j*x[81:162]).reshape(27,3);psi=(x[162:243]+1j*x[243:]).reshape(27,3);Q=M.dec(p['electric_charge'])
    measured=np.linalg.norm(Q@phi)**2+np.linalg.norm(Q@psi)**2
    assert measured>1e-8 and abs(measured-p['photon_mass_Gram'])<1e-10
    K=M.dec(p['colour']);assert max(np.linalg.norm(a@phi)+np.linalg.norm(a@psi) for a in K)<1e-9
    assert min(p['weak_Y_gauge_Gram_eigenvalues'])>1e-8


def test_flux_fourier_blocks_equal_full_operator():
    # Independent direct4D matrix at smallerL verifies exact block assembly.
    L=4;t2=.31;t3=.23;gamma,g5=M.N.P.N.M.P.gamma_matrices();u,v=M.N.P.N.uniform_flux_links(L)
    sites=list(product(range(L),range(L),range(2),range(2)));lookup={s:i for i,s in enumerate(sites)};n=len(sites);W=3*np.eye(4*n,dtype=complex)
    for mu in range(4):
        T=np.zeros((n,n),complex)
        for i,s in enumerate(sites):
            q=list(s);q[mu]=(q[mu]+1)%([L,L,2,2][mu]);link=[u[s[0],s[1]],v[s[0],s[1]],np.exp(1j*t2/2),np.exp(1j*t3/2)][mu]
            if mu==3 and s[3]==1:link*=-1
            T[i,lookup[tuple(q)]]=link
        W-=.5*(np.kron(T,np.eye(4)-gamma[mu])+np.kron(T.conj().T,np.eye(4)+gamma[mu]))
    H=np.kron(np.eye(n),g5)@W;w,E=np.linalg.eigh(H);P=E[:,w>0]@E[:,w>0].conj().T;assembled=np.zeros_like(P)
    for p2,p3 in product([0.,np.pi],[np.pi/2,3*np.pi/2]):
        F,g,index=M.flux_frame(t2,t3,p2,p3,L);embed=np.zeros((4*n,4*L*L),complex)
        for i,(x,y,z,t) in enumerate(sites):
            for a in range(4):embed[4*i+a,4*(x*L+y)+a]=np.exp(1j*(p2*z+p3*t))/2
        B=embed@F;assembled+=B@B.conj().T
    assert np.linalg.norm(assembled-P)<1e-10


def test_flux_chart_and_rectangle_have_same_orientation():
    p=packet()['flux'];assert abs(p['magnetic_flux']-1)<1e-10 and p['maximum_plaquette']<1/30
    for r in p['rows']:
        assert r['Wilson_gap']>.1 and r['minimum_overlap_singular_value']>.9
        assert abs(r['curvatures'][0]-r['curvatures'][1])<1e-3
        assert abs(r['rectangle_phase_over_area']-r['curvatures'][1])<.01
    assert 'not a proved spatially local' in p['scope']


def test_gravity_schlafli_under_independent_direction():
    p=packet()['gravity'];X=np.array(p['coordinates']);d=np.random.default_rng(11463).normal(size=X.shape);h=1e-5
    A,T,V=M.lorentz_terms(X);tp=M.lorentz_terms(X+h*d)[1];tm=M.lorentz_terms(X-h*d)[1]
    assert abs(A@(tp-tm)/(2*h))<1e-6
    assert np.linalg.norm(M.dec(p['combined_force'])-M.dec(p['gravity_force'])-np.array(p['matter_force']))<1e-10
    assert p['boost_error']<1e-6 and p['translation_force_error']<1e-5
    assert all(r['type']=='timelike' for r in p['causal_hinges'])


def test_retained_determinant_slice_is_scoped_finite_eft():
    p=packet()['joint'];assert np.linalg.norm(p['force'])<1e-3 and np.linalg.norm(p['half_step_force'])<1e-3
    assert min(np.linalg.eigvalsh(p['Hessian']))>50
    assert p['retained_heavy_species']==972 and p['high_energy_b0']==-637
    assert p['finite_cutoffs'][0]['inverse_bare_g2']>0 and p['finite_cutoffs'][-1]['inverse_bare_g2']<0
    assert 'not UV completion' in p['scope']


def test_reset_pulses_and_control_obstruction():
    p=packet()['reset'];target=M.N.reset_target();U=np.eye(144,dtype=complex);X=np.array([[0,1],[1,0]],complex)
    pulse=expm(1j*np.pi*np.eye(2)/2)@expm(-1j*np.pi*X/2)
    for a,b in p['transpositions']:U[[a,b]]=pulse@U[[a,b]]
    assert np.linalg.norm(U-np.eye(144)[:,target['permutation']])<1e-10
    P=np.kron(np.diag([1.,1.]+[0.]*10),np.eye(12));assert np.linalg.norm(P@U-U@P)>1
    assert p['pulse_count']==2*len(p['transpositions'])


def test_ideal_relay_reuse_does_not_inject_endpoint_faults():
    p=packet()['reset'];assert len(p['relay_reuse'])==105 and p['relay_residues']==72
    assert all(r['endpoint_difference']==0 for r in p['relay_reuse'])
    # The actual seven-gate unitary is CNOT endpoints tensor relay identity,
    # checked on all basis states, including both initial relay values.
    for bits in product([0,1],repeat=3):
        b=list(bits)
        for c,t in [(0,1),(1,0),(0,1),(1,2),(0,1),(1,0),(0,1)]:b[t]^=b[c]
        expected=list(bits);expected[2]^=expected[0];assert b==expected


def test_cubic_colour_orientation_is_central_and_chiral_charge_pattern():
    p=packet()['charges'];C=M.dec(p['colour_cubic_Casimir']);K=M.dec(p['colour'])
    assert max(np.linalg.norm(C@k-k@C) for k in K)<1e-10
    assert np.max(abs(np.linalg.eigvalsh(C)-np.repeat([-1.,0.,1.],9)))<1e-9
    for r in p['branching']:
        colour,weak,y=r['key'];expected=0 if colour==0 else (1 if y in [1,-2] else -1)
        assert r['colour_orientation']==expected


def test_photon_preserving_search_checks_full_transverse_hessian():
    p=packet()['photon_candidate'];assert p['gradient_norm']<5e-6 and p['photon_mass_Gram']<1e-14
    assert len(p['normal_eigenvalues'])==324-p['common_gauge_rank']
    assert p['lowest_direction_replay_error']<1e-3
    assert p['fixed_complex_dimension']==15
    # The candidate classification must agree with its independently replayed
    # negative direction if a resolved transverse negative mode is present.
    if p['negative_normal_count']:
        assert p['normal_eigenvalues'][0]<-1e-3
        assert p['lowest_direction_energy_scans'][-1]['symmetric_energy_change']<0


def test_vacuum_rank_sensitivity_is_not_an_exact_symmetry_claim():
    p=packet()['rank_sensitivity']
    assert len(p['rows'])==3
    for r in p['rows']:
        sv=np.array(r['singular_values']);assert len(sv)==86
        assert [int(sum(sv>row['threshold'])) for row in r['threshold_ranks']]==[row['rank'] for row in r['threshold_ranks']]
        assert r['threshold_ranks'][-1]['rank']==58
        assert r['threshold_ranks'][0]['rank']>r['threshold_ranks'][-1]['rank']
    assert 'No exact limiting stationary stabilizer' in p['scope']
