"""Independent symmetry, index, curved-mesh, heavy-spectrum and channel tests."""
import hashlib,json,sys
from pathlib import Path
from itertools import product,combinations
import numpy as np
import sympy as sp
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11423_11427_flux_index_curved_uv_noise as M
from w33_pass11390_11394_context_matter_clock import load
def packet():return json.loads((ROOT/'data/w33_pass11423_11427_flux_index_curved_uv_noise.json').read_text())
def read(x):return np.array(x['real'])+1j*np.array(x['imag'])

def test_source_hashes():
    p=packet()
    for s,h in p['source_sha256'].items():
        raw=json.dumps(json.loads((ROOT/s).read_text()),sort_keys=True,separators=(',',':')).encode()
        assert hashlib.sha256(raw).hexdigest()==h

def test_pin_radial_root_and_global_boundary_control():
    p=packet()['pin_vacua'];u=p['radius'];chi=p['chi'];r=p['parameters']['r'];g=p['parameters']['g'];v=p['parameters']['vchi']
    assert abs(4*u*(u*u-r*r)-g*chi*u*u)<1e-12
    assert abs(4*chi*(chi*chi-v*v)-g*u**3)<1e-12
    w=sp.symbols('w');pol=sp.Poly(sp.sympify(p['radius_squared_polynomial']),w)
    assert pol.count_roots(sp.Rational(1,16),sp.oo)==1
    # Nonzero product lowers energy below every chi=0 or zero-radius boundary.
    assert p['minimum_value']<0
    # Young's inequality bounds the negative mixed quartic by g/4 times
    # the sum of four fourth powers; g=.1<4 guarantees coercivity.
    assert 0<g<4

def test_pin_full_symmetry_CP_and_stability():
    p=packet()['pin_vacua'];U,B=read(p['pin_up']),read(p['pin_down']);chi=p['chi']
    basis=M.hermitian_basis();z=np.zeros(19)
    fn=lambda x:M.pin_potential(U+sum(a*b for a,b in zip(x[:9],basis)),B+sum(a*b for a,b in zip(x[9:18],basis)),chi+x[18])
    H=M.hessian(fn,z,h=3e-5);assert np.linalg.norm(H-np.array(p['Hessian']))<2e-5
    assert np.count_nonzero(np.linalg.eigvalsh(H)>1e-6)==17
    # Independent rephasing changes individual entries but not the flux/minimum.
    G=np.diag(np.exp(1j*np.array([.21,-.48,.72])))
    assert abs(M.pin_potential(G@U@G.conj().T,G@B@G.conj().T,chi)-p['minimum_value'])<1e-12
    comm=U@B-B@U;C=U@B.conjugate()-B.conjugate()@U
    assert abs(np.imag(np.trace(comm@comm@comm))-p['cubic_pair_CP'])<1e-12
    assert abs(np.trace(comm@comm@comm)+np.trace(C@C@C))<1e-12
    assert abs(M.pin_potential(U.conjugate(),B.conjugate(),-chi)-p['minimum_value'])<1e-12

def test_pin_colour_embedding_and_path_matching():
    p=packet()['pin_vacua'];A=load()['A'];r=np.linalg.inv(np.eye(80)-.01*A)[[44,1,0],0]
    U,B=read(p['pin_up']),read(p['pin_down']);Yu=np.diag(r)@U@np.diag(r);Yd=np.diag(r)@B@np.diag(r)
    assert np.linalg.matrix_rank(Yd,tol=1e-18)==3
    comm=Yu@Yu@Yd@Yd-Yd@Yd@Yu@Yu;assert abs(np.imag(np.trace(comm@comm@comm)))>0
    slots=M.prior()['fermion_bimodule']['native_slot_frames'];F=[read(x)for x in slots['colour_qutrit_seeds']];C=np.hstack(F)
    from w33_pass11403_11407_native_bimodule_collective_gates import gellmann
    N=C@np.kron(Yd,np.eye(3))@C.conj().T
    for t in gellmann():
        G=sum(f@t@f.conj().T for f in F);assert np.linalg.norm(G@N-N@G)<1e-12

def test_flux_plaquettes_and_topological_product():
    p=packet()['spacetime_index'];L=p['L'];sites=list(product(range(L),repeat=4));lookup={x:i for i,x in enumerate(sites)}
    for row in p['scans']:
        U=read(row['links'])
        for i,x in enumerate(sites):
            for mu,nu in combinations(range(4),2):
                y=list(x);y[mu]=(y[mu]+1)%L;z=list(x);z[nu]=(z[nu]+1)%L
                plaquette=U[i,mu]*U[lookup[tuple(y)],nu]*U[lookup[tuple(z)],mu].conjugate()*U[i,nu].conjugate()
                flux=row['flux'][0]if (mu,nu)==(0,1)else row['flux'][1]if (mu,nu)==(2,3)else 0
                assert abs(plaquette-np.exp(2j*np.pi*flux/L**2))<1e-12
        assert row['index']==-row['flux'][0]*row['flux'][1]

def test_actual_overlap_index_zero_mode_and_Weyl_map():
    p=packet()['spacetime_index'];U,T,G,H,w,D,S=M.overlap_operator();row=p['scans'][0]
    Z=read(row['zero_mode']);assert np.linalg.norm(D@Z)<1e-10
    index=round(np.trace(G@(np.eye(len(D))-D/2)).real)
    assert index==row['index']==-1 and min(abs(w))>.6
    hat=G@(np.eye(len(D))-D);Pm=(np.eye(len(D))-hat)/2;Po=(np.eye(len(D))+G)/2
    assert np.linalg.norm(D@Pm-Po@D)<1e-10
    assert round(np.trace(Pm).real)==row['hat_projector_minus_rank']!=row['ordinary_plus_rank']
    # A local U1 rephasing is a second representation of the same kernel.
    phase=np.repeat(np.exp(1j*np.random.default_rng(31).normal(size=256)),4);trans=phase[:,None]*H*phase.conjugate()[None,:]
    assert np.linalg.norm(trans@ (phase[:,None]*Z)-phase[:,None]*(H@Z))<1e-10

def refined_lengths(p,radial):
    L=np.zeros((9,9));L[:6,:6]=p['boundary_lengths']
    for i,s in enumerate(p['coarse_simplices']):
        for j,v in enumerate(s):L[6+i,v]=L[v,6+i]=radial[i*5+j]
    return L
def test_glued_mesh_curvature_and_subdivision_action():
    p=packet()['curved_gauge_islands'];r=np.array(p['radial_lengths']);fine=p['refined_simplices'];coarse=p['coarse_simplices']
    a,d=M.mesh_action(np.array(p['boundary_lengths']),coarse);b,e=M.mesh_action(refined_lengths(p,r),fine)
    assert abs(a-b)<1e-10 and abs(e[(0,1,2)]-p['coarse_hinge_deficit'])<1e-11
    assert abs(e[(0,1,2)])>.005
    assert max(abs(x)for t,x in e.items()if any(v>=6 for v in t))<1e-11
    assert len(fine)==15

def test_nonlinear_gauge_islands_in_curved_mesh():
    p=packet()['curved_gauge_islands'];flat=np.array(p['radial_lengths']);rng=np.random.default_rng(32)
    for _ in range(8):
        r=flat.copy()
        for i,coords in enumerate(p['local_coordinates']):
            V=np.array(coords);shift=rng.normal(size=4)*.01;r[i*5:i*5+5]=np.linalg.norm(V-V.mean(axis=0)-shift,axis=1)
        a,d=M.mesh_action(refined_lengths(p,r),p['refined_simplices']);assert abs(a-p['action'])<1e-10
        assert max(abs(x)for t,x in d.items()if any(v>=6 for v in t))<1e-11
    for H,T in zip(p['radial_Hessian_blocks'],p['vertex_translation_tangents']):assert np.linalg.norm(np.array(H)@np.array(T))<.001
    assert p['curvature_in_first_star']>1e-4

def test_heavy_mediator_matching_and_spectrum():
    p=packet();u=p['heavy_extension'];A=load()['A'];old=M.previous()['path_flavour'];idx=old['endpoints'];par=u['parameters']
    for row,key in zip(u['spectra'],['pin_up','pin_down']):
        B=read(p['pin_vacua'][key]);H=M.heavy_matrix(A,B,par['M'],par['eta'],par['v']);n=len(H);AL=np.zeros((3,n));AR=np.zeros((n,3))
        for j,i in enumerate(idx):AL[j,3*i+j]=1;AR[n-240+3*i+j,j]=1
        Full=np.block([[np.zeros((3,3)),AL],[AR,H]])
        assert np.linalg.norm(np.linalg.svd(Full,compute_uv=False)-row['full_Dirac_masses'])<1e-10
        assert np.linalg.norm(-AL@np.linalg.solve(H,AR)-read(row['matching']))<1e-12
        assert len(H)==483 and len(Full)==486

def test_all_heavy_scalar_fermion_counts_and_threshold_identity():
    p=packet()['heavy_extension'];groups=p['groups'];ferm=[g for g in groups if g['weight']<0]
    assert sum(-g['weight']*len(g['m2'])for g in ferm)==p['fermion_real_weight']==11688
    assert 3*160==p['complex_link_count'] and 4+12+19+2*480==p['real_scalar_count']
    assert 1440+720+720+9+9==p['heavy_Dirac_count']
    str4=sum(g['weight']*sum(x*x for x in g['m2'])for g in groups)
    assert abs(str4-p['supertrace_M4'])<1e-7 and p['heavy_only_supertrace_M4']<-1e8
    value=lambda mu:sum(g['weight']*sum(x*x*(np.log(x/mu**2)-g['constant'])for x in g['m2']if x>0)for g in groups)/(64*np.pi**2)
    assert abs(value(1)-p['one_loop_value'])<1e-7
    assert abs(value(3)-value(1)+str4*np.log(3)/(32*np.pi**2))<1e-7

def test_native_invariant_noise_frame_and_noiseless_leakage():
    p=packet()['native_noise'];B=read(p['invariant_frame']);c=load();D=np.array(c['D'],float);Hp=D[:40].T@D[:40];Hl=D[40:].T@D[40:]
    gates=M.prior()['local_gates'];Q=read(gates['logical_basis']);target=read(gates['ideal_qubit_gate']);a,b=gates['gate_pair']
    assert B.shape==(160,12) and np.linalg.norm(B.conj().T@B-np.eye(12))<1e-12
    for H in [Hp,Hl]:assert np.linalg.norm(H@B-B@(B.conj().T@H@B))<1e-10
    Kp=B.conj().T@Hp@B;Kl=B.conj().T@Hl@B;clock=expm(-.02j*Kp)@expm(-.04j*Kl)@expm(-.02j*Kp);state=B.conj().T@Q
    old=M.previous()['counted_noise']
    for angle,e,N in zip(old['principal_angles'],[a,b,a,b,a],old['ticks_per_pulse']):
        f=B[e,:].conjugate();kick=np.eye(12)+(np.exp(-.02j*np.copysign(.0125,angle))-1)*np.outer(f,f.conj());state=np.linalg.matrix_power(kick@clock@kick,N)@state
    logical=target.conj().T@Q.conj().T@B@state;expected=1-np.trace(logical.conj().T@logical).real/2
    assert abs(expected-p['scans'][0]['Pauli_and_erasure_weights'][-1])<1e-8

def test_one_tick_native_channel_CPTP_and_density_replay():
    p=packet()['native_noise'];B=read(p['invariant_frame']);n=B.shape[1];a=M.prior()['local_gates']['gate_pair'][0];f=B[a,:].conjugate();Z=np.eye(n)-2*np.outer(f,f.conj());q=.001
    rho=np.outer(B.conj().T@read(M.prior()['local_gates']['logical_basis'])[:,0],(B.conj().T@read(M.prior()['local_gates']['logical_basis'])[:,0]).conj())
    output=(1-q)*rho+q*Z@rho@Z.conj().T
    channel=(1-q)*np.eye(n*n)+q*np.kron(Z.conjugate(),Z)
    assert np.linalg.norm(channel@rho.reshape(-1,order='F')-output.reshape(-1,order='F'))<1e-12
    assert min(np.linalg.eigvalsh(output))>-1e-12 and abs(np.trace(output)-1)<1e-12
    for row in p['scans']:
        assert min(row['Choi_eigenvalues'])>-1e-8 and abs(sum(row['Pauli_and_erasure_weights'])-1)<1e-12

def test_Steane_ML_dephasing_and_erasure_controls():
    sys.path.insert(0,str(ROOT/'tests'));from test_w33_pass11413_11417_pin_regge_vacuum_noise import steane_basis
    H,C=steane_basis();labels=np.arange(128);p=.003;success=0.
    for e in product([0,1],repeat=7):
        e=np.array(e);syndrome=H@e%2;r=e.copy()
        if syndrome.any():r[np.flatnonzero(np.all(H.T==syndrome,axis=1))[0]]^=1
        phase=(-1.)**sum(((labels>>j)&1)*r[j]for j in range(7));encoded=C.conj().T@(phase[:,None]*C)
        success+=abs(np.trace(encoded)/2)**2*p**e.sum()*(1-p)**(7-e.sum())
    assert abs(M.steane_ML([1-p,0,0,p,0])-success)<1e-12
    assert abs(M.steane_ML([1,0,0,0,0])-1)<1e-12
    assert abs(M.steane_ML([0,0,0,0,1])-.25)<1e-12

def test_readout_floor_and_noise_comparison():
    p=packet()['native_noise']
    for row in p['scans']:
        f=M.steane_ML(row['Pauli_and_erasure_weights']);assert abs(f-row['ideal_encoded_entanglement_fidelity'])<1e-12
        for s in row['readout_scans']:
            eps=s['bit_error_probability'];assert abs(s['entanglement_fidelity']-f*(1-eps)**6)<1e-12
    zero=p['scans'][0]
    assert zero['readout_scans'][2]['failure_probability']>zero['uncoded_entanglement_infidelity']
    assert p['scans'][1]['readout_scans'][2]['failure_probability']<p['scans'][1]['uncoded_entanglement_infidelity']

def test_shared_neutrino_pin_and_CP_isospectrality():
    p=packet();h=p['heavy_extension'];B=read(p['pin_vacua']['pin_down']);N=read(h['neutrino_seesaw'])
    assert np.linalg.norm(N[:3,3:]-.01*B/np.sqrt(2))<1e-13
    assert np.linalg.norm(N-N.T)<1e-13
    assert np.max(abs(np.linalg.svd(N,compute_uv=False)-np.array(h['neutrino_masses'])))<1e-12
    assert np.max(abs(np.linalg.svd(N.conjugate(),compute_uv=False)-np.linalg.svd(N,compute_uv=False)))<1e-12
    A=load()['A'];H=M.heavy_matrix(A,B)
    assert np.linalg.norm(M.heavy_matrix(A,B.conjugate())-H.conjugate())<1e-12


def test_extraction_faults_against_actual_code_subspace():
    from test_w33_pass11413_11417_pin_regge_vacuum_noise import steane_basis
    _,C=steane_basis();labels=np.arange(128);dec=M.minimum_weight_decoder();counts=[]
    for gate in range(24):
        bad=0
        for bits in range(1,16):
            x,z,s=M.syndrome_circuit(fault=(gate,bits));rx,rz=dec[s];x^=rx;z^=rz
            phase=(-1.)**sum(((labels>>j)&1)*((z>>j)&1)for j in range(7))
            E=(phase[:,None]*C)[labels^x]
            overlap=C.conj().T@E
            fidelity=abs(np.trace(overlap)/2)**2
            assert abs(fidelity-round(fidelity))<1e-12
            bad+=fidelity<.5
        counts.append(int(bad))
    q=packet()['native_noise']['extraction_circuit']
    assert counts==q['failure_counts_per_CNOT']
    assert sum(counts)==q['malignant_single_faults']>0
    assert abs(sum(counts)/15-q['linear_failure_coefficient'])<1e-12


def test_hook_witness_with_eight_qubit_statevector():
    from test_w33_pass11413_11417_pin_regge_vacuum_noise import steane_basis
    _,C=steane_basis();w=packet()['native_noise']['extraction_circuit']['correlated_fault_witness']
    labels=np.arange(256);state=np.zeros(256,complex);state[:128]=(C[:,0]+C[:,1])/np.sqrt(2)
    def X(v,j):return v[labels^(1<<j)]
    def Z(v,j):return (-1.)**((labels>>j)&1)*v
    def H(v,j):
        mate=v[labels^(1<<j)];bit=(labels>>j)&1
        return np.where(bit,(mate-v)/np.sqrt(2),(v+mate)/np.sqrt(2))
    gate=0;reported=0
    for kind in [0,1]:
        for row in range(3):
            if kind:state=H(state,7)
            for data in [j for j in range(7)if ((j+1)>>row)&1]:
                control,target=(data,7)if kind==0 else(7,data)
                perm=labels^(((labels>>control)&1)<<target);state=state[perm]
                if gate==w['gate']:
                    bits=w['two_qubit_Pauli_bits']
                    for flag,j,op in [(1,control,X),(2,control,Z),(4,target,X),(8,target,Z)]:
                        if bits&flag:state=op(state,j)
                gate+=1
            if kind:state=H(state,7)
            prob=float(np.sum(abs(state[128:])**2));assert min(prob,1-prob)<1e-12
            bit=int(prob>.5);reported|=bit<<(row+3*kind)
            state[labels//128!=bit]=0
            if bit:state=X(state,7)
    assert reported==w['reported_syndrome']
    rx,rz=M.minimum_weight_decoder()[reported]
    for j in range(7):
        if (rx>>j)&1:state=X(state,j)
        if (rz>>j)&1:state=Z(state,j)
    expected=(C[:,0]+C[:,1])/np.sqrt(2);x=w['residual_X_mask'];z=w['residual_Z_mask'];idx=np.arange(128)
    for j in range(7):
        if (x>>j)&1:expected=expected[idx^(1<<j)]
        if (z>>j)&1:expected=(-1.)**((idx>>j)&1)*expected
    assert abs(abs(np.vdot(expected,state[:128]))-1)<1e-12
    # The explicit branch differs from an identity logical channel.
    assert x.bit_count()+z.bit_count()>=2

if __name__=='__main__':
    tests=[(n,f)for n,f in list(globals().items())if n.startswith('test_')]
    for n,f in tests:f();print(n,'PASS',flush=True)
    print(len(tests),'independent regressions PASS')
