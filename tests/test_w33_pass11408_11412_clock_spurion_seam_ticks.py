"""Independent native-map, covariance, perturbation and counted-unitary replays."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
from w33_pass11408_11412_clock_spurion_seam_ticks import read
from w33_pass11403_11407_native_bimodule_collective_gates import gellmann,extended
from w33_pass11390_11394_context_matter_clock import load
def packet():return json.loads((ROOT/'data/w33_pass11408_11412_clock_spurion_seam_ticks.json').read_text())
def previous():return json.loads((ROOT/'data/w33_pass11403_11407_native_bimodule_collective_gates.json').read_text())

def test_source_manifest_and_global_alignment_bound():
    p=packet();old=previous();key=next(iter(p['source_sha256']))
    assert hashlib.sha256(json.dumps(old,sort_keys=True,separators=(',',':')).encode()).hexdigest()==p['source_sha256'][key]
    D=np.diag(p['competing_clocks']['supplied_orbit_radii']);rng=np.random.default_rng(801)
    for row in p['competing_clocks']['scans']:
        B=read(row['source']);U=read(row['orientation']);H=U@D@U.conj().T
        assert np.linalg.norm(H-read(row['selected_family_operator']))<1e-12
        best=-np.trace(H@B).real
        for _ in range(40):
            X=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3));V=np.linalg.qr(X)[0]
            assert -np.trace(V@D@V.conj().T@B).real>=best-1e-12

def test_alignment_normal_Hessian_and_CP_conjugation():
    p=packet()['competing_clocks'];d=np.array(p['supplied_orbit_radii']);D=np.diag(d)
    for row in p['scans'][::2]:
        U=read(row['orientation']);B=read(row['source']);curv=[]
        for i,j in [(0,1),(0,2),(1,2)]:
            for imaginary in [False,True]:
                T=np.zeros((3,3),complex);T[i,j]=1j if imaginary else 1;T[j,i]=T[i,j].conjugate()
                def V(t):
                    R=U@expm(1j*t*T);return -np.trace(R@D@R.conj().T@B).real
                h=.0001;curv.append((V(h)+V(-h)-2*V(0))/h**2)
        assert np.linalg.norm(np.array(curv)-row['orbit_Hessian_eigenvalues'])<1e-6
        J=np.imag(U[0,0]*U[1,1]*U[0,1].conjugate()*U[1,0].conjugate())
        C=U.conjugate();Jc=np.imag(C[0,0]*C[1,1]*C[0,1].conjugate()*C[1,0].conjugate())
        assert abs(J+Jc)<1e-15 and abs(J)>1e-15
    r=p['scans'][0];assert max(r['angles_degrees'])<.1

def test_competing_clock_cubic_J_and_small_mixing_ratio():
    p=packet()['competing_clocks'];d=np.array(p['supplied_orbit_radii']);F=read(previous()['colour_safe_bridge']['other_native_clock_mixing'])
    source=F@np.diag(d)@F.conj().T;off=[abs(source[i,j])for i,j in [(0,1),(0,2),(1,2)]]
    assert max(off)-min(off)<1e-12
    for row in p['scans']:
        b=np.array(row['source_eigenvalues']);gap=lambda v:(v[1]-v[0])*(v[2]-v[0])*(v[2]-v[1])
        assert abs(abs(row['J'])-abs(row['r']**3*p['near_identity_J_over_r3']*gap(d)/gap(b)))<1e-12
    tiny=p['scans'][0]['angles_degrees'];assert abs(tiny[2]/tiny[1]-(d[2]-d[1])/(d[2]-d[0]))<.001
    assert p['small_mixing_ratio_obstruction']['theta13_over_theta23_limit']>.999
    r,z=sp.symbols('r z');char=sp.sympify(p['exact_characteristic_polynomial']);disc=sp.sympify(p['positive_discriminant_decomposition'])
    assert sp.expand(sp.discriminant(char,z)-disc)==0
    # Verify the characteristic coefficients on the actual native clock source.
    for row in p['scans']:
        coeff=np.array([float(x)for x in sp.Poly(char.subs(r,row['r']),z).all_coeffs()])
        assert np.linalg.norm(np.poly(read(row['source']))-coeff)<1e-11

def test_full_quark_Yukawa_tensor_Ward_identity():
    # Build the entire18 x9 x2 coupling tensor instead of testing one field value.
    Y=np.kron(np.array([[1,2j,0],[3,4,5],[1j,0,6]],complex),np.eye(3))
    T=np.zeros((18,9,2),complex)
    for c in range(9):
        for w in range(2):T[2*c+w,:,w]=Y[:,c]
    pauli=[np.array([[0,1],[1,0]])/2,np.array([[0,-1j],[1j,0]])/2,np.diag([1,-1])/2]
    gauges=[(np.kron(np.eye(3),G),np.zeros((2,2)))for G in gellmann()]+[(np.zeros((9,9)),G)for G in pauli]
    for Gc,Gw in gauges:
        Rq=np.kron(Gc,np.eye(2))-np.kron(np.eye(9),Gw.T);Ru=-Gc.T
        variation=np.einsum('abs,ac->cbs',T,Rq)+np.einsum('abs,bc->acs',T,Ru)+np.einsum('abs,st->abt',T,Gw)
        assert np.linalg.norm(variation)<1e-12

def test_actual_neutrino_character_Majorana_and_charge_condition():
    p=packet()['yukawa_arrows'];old=previous();s=old['fermion_bimodule']['native_slot_frames'];L=read(s['lepton']);R=read(s['right_up'])
    operators=[R@L[:,i:i+1].conj().T for i in range(3)];_,gs,_=extended(load());M=np.array(p['H27_invariant_symmetric_Majorana_mask'])
    for g in gs:
        rho=np.array([[np.trace(A.conj().T@g@B@g.conj().T)for B in operators]for A in operators])
        assert np.linalg.norm(rho.T@M@rho-M)<1e-10
    assert np.linalg.matrix_rank(M)==2 and np.array_equal(M@M,np.diag([0,1,1]))
    assert p['conditional_hypercharges']=={'Q':'1/6','u_conjugate':'-2/3','d_conjugate':'1/3','L':'-1/2','e_conjugate':'1','nu_conjugate':'0'}
    # Gauge invariance fails when the old B-L ambiguity is retained.
    wrong_charge=.1;assert np.linalg.norm(np.exp(2j*wrong_charge)*M-M)>.1

def test_messenger_exact_native_resolvent_and_graph_paths():
    p=packet()['messenger_locality'];A=np.array(p['native_adjacency']);indices=p['endpoint_indices']
    # Exact symbolic inverse from the three native eigenvalues, not float fitting.
    t=sp.symbols('t');den=1+2*t-8*t*t;alpha=(1+2*t)/den;beta=t/den;gamma=4*t*t/(den*(1-12*t))
    R=alpha*sp.eye(40)+beta*sp.Matrix(A)+gamma*sp.ones(40)
    assert all(sp.cancel(x)==0 for x in (sp.eye(40)-t*sp.Matrix(A))*R-sp.eye(40))
    paths=np.eye(40,dtype=np.int64);first=[]
    for degree in range(5):
        first.append(paths[indices,0].tolist());paths=paths@A
    assert np.array_equal(np.array(first).T,p['path_count_coefficients'])
    assert [next(j for j,v in enumerate(row)if v)for row in np.array(first).T]==p['distances']
    assert p['minimum_link_spurion_degrees']==[0,2,4]

def test_messenger_spurion_covariance_and_shortcut_failure():
    p=packet()['messenger_locality'];A=np.array(p['native_adjacency']);pin=np.zeros((40,40));pin[0,0]=1
    eta=.003;K=np.eye(40)-eta*A;R=np.linalg.inv(K);rng=np.random.default_rng(802)
    G,H=[np.diag(np.exp(1j*rng.normal(size=40)))for _ in range(2)]
    left=np.linalg.inv(G@K@G.conj().T)@(G@pin@H.conj().T)@np.linalg.inv(H@K@H.conj().T)
    assert np.linalg.norm(left-G@R@pin@R@H.conj().T)<1e-12
    i=p['endpoint_indices'][2];closed=(R@pin@R)[i,i]
    shortcut=pin.copy();shortcut[i,i]=.001
    assert (R@shortcut@R)[i,i]>10000*closed

def test_common_pin_rankone_and_three_copy_Schur_matching():
    p=packet()['messenger_locality'];A=np.array(p['native_adjacency']);ids=p['endpoint_indices'];K=np.eye(40)-.01*A
    P=np.zeros((40,40));P[0,0]=1;R=np.linalg.inv(K);family=.1*np.outer(R[ids,0],R[0,ids])
    assert np.linalg.matrix_rank(family,tol=1e-12)==1 and np.linalg.norm(family-read(p['common_pin_three_family_mass']))<1e-12
    repaired=read(p['repaired_three_copy_family_mass']);assert np.linalg.matrix_rank(repaired)==3
    # Direct heavy equations with each family source, not an inverse formula.
    for j,i in enumerate(ids):
        heavy=np.block([[K,.1*P],[np.zeros((40,40)),K]]);source=np.eye(80)[:,40+i];wave=np.linalg.solve(heavy,source)
        assert abs(-wave[i]-repaired[j,j])<1e-12
    assert 'not constructed' in p['joint_alignment_boundary']

def test_seam_maps_in_both_native_and_second_FCC_bases():
    p=packet()['coherent_seams'];F=sp.Matrix(p['FCC_primitive']);R=sp.Matrix(p['rational_rotation']).applyfunc(sp.Rational)
    M=sp.Matrix(p['coincidence_left_periods']);N=sp.Matrix(p['coincidence_right_periods']);B=sp.Matrix(p['native_to_second_FCC_basis'])
    native=sp.Matrix(json.loads((ROOT/'data/w33_pass11389_parabolic_spatial_cover.json').read_text())['line']['fcc_roots']).applyfunc(sp.Rational)
    assert native*B==F and abs(B.det())==1 and F*M==R*F*N and abs(M.det())==abs(N.det())==5
    assert R.T*R==sp.eye(3)and R.det()==1
    # Count all mod5 residue classes. Exactly25 of125 survive, so index5 is exact.
    from itertools import product
    Q=F.inv()*R.inv()*F
    count=sum(all(v.q==1 for v in Q*sp.Matrix(t))for t in product(range(5),repeat=3));assert count==25

def test_nonlinear_bracket_defect_and_refinement():
    p=packet()['coherent_seams'];rows=p['nonlinear_transport_scans']
    assert all(abs(rows[i+1]['closure_defect'])<abs(rows[i]['closure_defect'])/12 for i in range(3))
    n=64;x=2*np.pi*np.arange(n)/n;f=np.exp(4j*x);xi=np.exp(2j*x);eta=np.exp(3j*x)
    derivative=lambda f:(np.roll(f,-1)-np.roll(f,1))/2
    lhs=xi*derivative(eta*derivative(f))-eta*derivative(xi*derivative(f))
    rhs=(xi*derivative(eta)-eta*derivative(xi))*derivative(f)
    assert abs(np.linalg.norm(lhs-rhs)/np.sqrt(n)-p['discrete_plane_wave_bracket_defect'])<1e-12

def test_counted_ticks_full160_Hadamard_direct_power():
    c=load();D=np.array(c['D'],float);Hp=D[:40].T@D[:40];Hl=D[40:].T@D[40:];old=previous();Q=read(old['local_gates']['logical_basis']);a,b=old['local_gates']['gate_pair']
    ideal=read(old['local_gates']['ideal_qubit_gate']);row=packet()['finite_ticks']['Hadamard_scans'][-1];dt=row['dt'];u=row['strength'];state=Q.copy()
    clock=expm(-.5j*dt*Hp)@expm(-1j*dt*Hl)@expm(-.5j*dt*Hp)
    for e,angle,N in zip([a,b,a,b,a],old['local_gates']['five_pulse_phase_angles'],row['ticks_per_pulse']):
        kick=np.ones(160,complex);kick[e]=np.exp(-.5j*dt*np.copysign(u,angle));S=kick[:,None]*clock*kick[None,:]
        # Integer binary power of actual native160 operator; no producer eigendecomposition.
        state=np.linalg.matrix_power(S,N)@state
    error=np.linalg.norm(state-Q@ideal,ord=2)
    assert abs(error-row['total_encoded_error'])<2e-8 and error<row['dressing_and_rounding_bound']

def test_counted_native_pi_over_8_gate():
    c=load();D=np.array(c['D'],float);Hp=D[:40].T@D[:40];Hl=D[40:].T@D[40:];old=previous();Q=read(old['local_gates']['logical_basis']);a=old['local_gates']['gate_pair'][0]
    row=packet()['finite_ticks']['pi_over_8_gate_scans'][-1];dt,u=row['dt'],row['strength']
    clock=expm(-.5j*dt*Hp)@expm(-1j*dt*Hl)@expm(-.5j*dt*Hp)
    kick=np.ones(160,complex);kick[a]=np.exp(-.5j*dt*u);S=kick[:,None]*clock*kick[None,:]
    state=np.linalg.matrix_power(S,row['ticks'])@Q;target=Q@np.diag([np.exp(-1j*np.pi/4),1]);error=np.linalg.norm(state-target,ord=2)
    assert abs(error-row['total_encoded_error'])<2e-8 and error<row['dressing_and_rounding_bound']

def test_entangler_joint_frame_and_exact_counted_power():
    p=packet()['finite_ticks'];c=load();D=np.array(c['D'],float);Hp=D[:40].T@D[:40];Hl=D[40:].T@D[40:];B=read(p['entangler_star_frame']);edge=np.eye(160)[:,previous()['local_gates']['gate_pair'][0]]
    assert np.linalg.norm(B.conj().T@B-np.eye(B.shape[1]))<1e-10
    Kp=B.conj().T@Hp@B;Kl=B.conj().T@Hl@B
    assert max(np.linalg.norm(H@B-B@K)for H,K in [(Hp,Kp),(Hl,Kl)])<1e-10
    assert np.linalg.norm(B@B.conj().T@edge-edge)<1e-10
    # Build zero-mode projector independently from native incidence SVD.
    _,s,V=np.linalg.svd(D,full_matrices=True);P=V[79:].conj().T@V[79:]
    v=B.conj().T@(P@edge/np.sqrt(P[0,0]));logical=np.kron(v,v);f=B.conj().T@edge;delta=np.kron(f,f)
    row=p['entangler_scans'][-1];dt,g=row['dt'],row['coupling'];clock=expm(-.5j*dt*Kp)@expm(-1j*dt*Kl)@expm(-.5j*dt*Kp)
    kick=np.eye(len(delta))+(np.exp(-.5j*dt*g)-1)*np.outer(delta,delta.conj());S=kick@np.kron(clock,clock)@kick
    state=np.linalg.matrix_power(S,row['ticks'])@logical;err=np.linalg.norm(state+logical)
    assert abs(err-row['total_encoded_error'])<2e-8 and err<row['dressing_and_rounding_bound']

def test_floquet_budget_scopes_and_bright_gap():
    p=packet()['finite_ticks'];old=previous();c=load();D=np.array(c['D'],float);Hp=D[:40].T@D[:40];Hl=D[40:].T@D[40:]
    for dt in [.04,.02,.01]:
        U=expm(-.5j*dt*Hp)@expm(-1j*dt*Hl)@expm(-.5j*dt*Hp);angles=np.angle(np.linalg.eigvals(U))
        assert sum(abs(angles)<1e-9)==81 and min(abs(a)for a in angles if abs(a)>1e-9)>.01*dt
    for row in p['Hadamard_scans']:
        assert sum(row['ticks_per_pulse'])==row['total_ticks'] and row['total_encoded_error']<.02
    assert all(row['conditional_concurrence']>.999 for row in p['entangler_scans'])
    assert 'no fault tolerance' in p['adversarial_noise_bound']
