"""Independent realizations and physical-scope regressions for11428-11432."""
import sys,json,hashlib
from pathlib import Path
from itertools import product
import numpy as np
import sympy as sp
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11428_11432_native_alignment_measure_coarse_flags as M
P=M.P

def packet():return json.loads((ROOT/'data/w33_pass11428_11432_native_alignment_measure_coarse_flags.json').read_text())
def read(x):return np.array(x['real'])+1j*np.array(x['imag'])

def test_canonical_source_hashes():
    for name,want in packet()['source_sha256'].items():
        raw=json.dumps(json.loads((ROOT/name).read_text()),sort_keys=True,separators=(',',':')).encode();assert hashlib.sha256(raw).hexdigest()==want

def test_exact_native_Cartan_contraction_algebra():
    p=packet()['native_alignment'];E=np.array(p['integer_Cartan_embedding']);tc=np.array(p['integer_cubic_coefficients']);rng=np.random.default_rng(901)
    X=np.roll(np.eye(3),1,axis=0);Z=np.diag(np.exp(2j*np.pi*np.arange(3)/3));_,(_,d,_)=M.native_tensors()
    for _ in range(5):
        q=rng.normal(size=3)+1j*rng.normal(size=3);phi=np.einsum('aix,x->ai',E,q)/np.sqrt(3);G=phi.T@phi.conj();T=np.einsum('abc,ai,bj,ck->ijk',d,phi,phi,phi)
        assert np.linalg.norm(G-np.trace(G)/3*np.eye(3))<1e-12
        assert np.linalg.norm(T-np.einsum('ijkxyz,x,y,z->ijk',tc,q,q,q)/(3*np.sqrt(3)))<1e-11
        for F in [X,Z]:assert np.linalg.norm(np.einsum('ia,jb,kc,abc->ijk',F,F,F,T)-T)<1e-11
    x=sp.Matrix(X.astype(int));omega=(-1+sp.I*sp.sqrt(3))/2;z=sp.diag(1,omega,omega**2)
    constraints=sp.Matrix.vstack(sp.kronecker_product(sp.eye(3),x)-sp.kronecker_product(x.T,sp.eye(3)),sp.kronecker_product(sp.eye(3),z)-sp.kronecker_product(z.T,sp.eye(3)))
    assert constraints.rank()==8

def test_two_field_native_invariants_and_D_flatness():
    import w33_20261001_degree18_phase_completion as D
    p=packet()['native_alignment'];basis,tensors=M.native_tensors();H=M.compact_generators(basis);original=D.I.tensors;D.I.tensors=lambda:tensors
    try:
        for key in ['phi','psi']:
            field=read(p[key]).ravel();i6,i12,_,_=D.I.evaluate(field);i18,_=D.evaluate(field)
            a=(-1+1j*np.sqrt(7))/4
            assert max(abs(i6-a),abs(i12-a*a),abs(i18-a**3))<1e-8
            assert max(abs(np.einsum('i,aij,j->a',field.conj(),H,field)))<1e-11
    finally:D.I.tensors=original

def test_native_pin_covariance_CP_and_positive_rank():
    p=packet()['native_alignment'];phi,psi=read(p['phi']),read(p['psi']);C,U,B,q=M.cross_pin(phi,psi);rng=np.random.default_rng(902)
    h=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3));h=(h+h.conj().T)/2;h-=np.trace(h)/3*np.eye(3);F=expm(.3j*h)
    _,u,b,q2=M.cross_pin(phi@F.T,psi@F.T)
    assert max(np.linalg.norm(u-F@U@F.conj().T),np.linalg.norm(b-F@B@F.conj().T))<1e-12
    assert abs(q-q2)<1e-12 and q>.01
    assert min(np.linalg.eigvalsh(U))>=1-1e-10 and min(np.linalg.eigvalsh(B))>=1-1e-10
    qc=M.cross_pin(phi.conjugate(),psi.conjugate())[3];assert abs(q+qc)<1e-12
    assert abs(-p['chi']*q-(+p['chi'])*qc)<1e-12

def test_native_pin_CP_vs_spectral_Jarlskog():
    p=packet()['native_alignment'];U,B=read(p['pin_up']),read(p['pin_down']);a,A=np.linalg.eigh(U);b,V=np.linalg.eigh(B);mix=A.conj().T@V
    J=np.imag(mix[0,0]*mix[1,1]*mix[0,1].conjugate()*mix[1,0].conjugate());disc=lambda x:np.prod([x[j]-x[i]for i in range(3)for j in range(i+1,3)])
    assert abs(abs(6*J*disc(a)*disc(b))-abs(p['CP_flux']))<1e-12
    assert np.all(np.diff(p['flux_history'])>=-1e-12)

def test_native_alignment_actual_colour_embedding():
    p=packet()['native_alignment'];A=P.load()['A'];r=np.linalg.inv(np.eye(80)-.08*A)[P.previous()['path_flavour']['endpoints'],0];Y=np.diag(r)@read(p['pin_down'])@np.diag(r)
    frames=[P.read(f)for f in P.prior()['fermion_bimodule']['native_slot_frames']['colour_qutrit_seeds']];C=np.hstack(frames)
    from w33_pass11403_11407_native_bimodule_collective_gates import gellmann
    native=C@np.kron(Y,np.eye(3))@C.conj().T
    assert np.linalg.matrix_rank(native,tol=1e-15)==9
    for t in gellmann():
        gauge=sum(F@t@F.conj().T for F in frames);assert np.linalg.norm(native@gauge-gauge@native)<1e-12
    assert p['actual_colour_Ward_residual']<1e-11


def test_endpoint_spurion_common_family_realization():
    p=packet()['native_alignment'];U,B=read(p['pin_up']),read(p['pin_down']);A=P.load()['A'];D=np.diag(np.linalg.inv(np.eye(80)-.08*A)[P.previous()['path_flavour']['endpoints'],0])
    h=np.array([[.1,.2+.1j,.3],[.2-.1j,-.2,.1j],[.3,-.1j,.1]]);F=expm(.4j*h);dp=F@D@F.conj().T;u=F@U@F.conj().T;b=F@B@F.conj().T
    Yu=D@U@D;Yd=D@B@D;up=dp@u@dp;down=dp@b@dp
    assert np.linalg.norm(up-F@Yu@F.conj().T)<1e-12 and np.linalg.norm(down-F@Yd@F.conj().T)<1e-12
    comm=lambda a,b:a@a@b@b-b@b@a@a
    W=comm(Yu,Yd);Wp=comm(up,down)
    # This invariant is~1e-23: a rotated dense double matrix loses its
    # cancellations. Verify the same transformation at80-digit precision.
    with M.mp.workdps(80):
        matrix=lambda a:M.mp.matrix([[M.mp.mpc(str(v.real),str(v.imag))for v in row]for row in np.asarray(a,dtype=complex)])
        um,bm,dm,hm=map(matrix,[U,B,D,h]);fm=M.mp.expm(M.mp.mpf('.4')*1j*hm);dnew=fm*dm*fm.H
        ua=dm*um*dm;ba=dm*bm*dm;ub=dnew*(fm*um*fm.H)*dnew;bb=dnew*(fm*bm*fm.H)*dnew
        wa=ua*ua*ba*ba-ba*ba*ua*ua;wb=ub*ub*bb*bb-bb*bb*ub*ub
        trace=lambda a:sum(a[i,i]for i in range(3))
        assert abs(trace(wa*wa*wa)-trace(wb*wb*wb))<M.mp.mpf('1e-65')
        assert abs(trace(wa*wa*wa))>M.mp.mpf('1e-24')
    # A fixed endpoint source is additional data, not a free gauge transform.
    assert np.linalg.norm(D@u@D-F@Yu@F.conj().T)>1e-4


def test_charged_overlap_projector_derivative():
    p=packet()['local_chiral_measure'];b=np.array(p['background']);v=np.array(p['tangents']);eps=.005;h=1e-5
    proj,ds,_,gap,_,_=M.weak_overlap(1,eps,b,v);assert gap>.1
    plus=M.weak_overlap(1,1.,eps*b+h*v[0],[])[0];minus=M.weak_overlap(1,1.,eps*b-h*v[0],[])[0]
    assert np.linalg.norm((plus-minus)/(2*h)-ds[0])/(1+np.linalg.norm(ds[0]))<1e-7
    assert np.linalg.norm(proj@proj-proj)<1e-11
    assert np.linalg.norm(proj@ds[0]@proj)<1e-11

def test_Weyl_projector_gauge_covariance():
    p=packet()['local_chiral_measure'];b=np.array(p['background']);sites=list(product(range(3),repeat=4));lookup={x:i for i,x in enumerate(sites)};theta=np.random.default_rng(903).normal(size=81)*.2
    gauge=np.empty((81,4))
    for i,x in enumerate(sites):
        for mu in range(4):
            y=list(x);y[mu]=(y[mu]+1)%3;gauge[i,mu]=theta[i]-theta[lookup[tuple(y)]]
    a=M.weak_overlap(2,.005,b,[])[0];trans=M.weak_overlap(2,1.,.005*b+gauge,[])[0];G=np.diag(np.repeat(np.exp(2j*theta),4))
    assert np.linalg.norm(trans-G@a@G.conj().T)<1e-10


def test_anomaly_cancellation_and_finite_curvature_boundary():
    p=packet()['local_chiral_measure'];q=np.array(p['integer_hypercharges']);d=np.array(p['multiplicities']);assert np.dot(q,d)==np.dot(q**3,d)==0
    rows=p['scans']
    for r in rows:
        assert abs(np.dot(d,r['curvatures'])-r['multiplet_curvature'])<1e-12
        assert r['minimum_Wilson_gap']>.1
    # Linear response cancels, finite lattice bundle curvature need not.
    assert abs(rows[-1]['multiplet_curvature'])>1e-6
    ratios=[rows[i]['multiplet_curvature']/rows[i+1]['multiplet_curvature']for i in range(2)]
    assert all(abs(x-8)<.1 for x in ratios)

def test_local_frame_holonomy_and_curvature():
    p=packet()['local_chiral_measure'];b=np.array(p['background']);v=np.array(p['tangents'])
    cur=M.weak_overlap(1,1.,.005*b+.0005*(v[0]+v[1]),v)[2]
    assert abs(abs(p['patch_loop_phase']/p['loop_area'])-abs(cur))<.06*abs(cur)+1e-8
    assert p['local_frame_rank']==162

def test_Regge_volume_gradient_and_remaining_equations():
    p=packet()['coarse_constraints'];bd=np.array(p['boundary_lengths']);rad=np.array(p['stationary_radials']);act,g,_,_=M.regge_volume_gradient(bd,rad,p['Lambda']);h=1e-6
    for j in [0,3,7,12]:
        step=np.eye(15)[j]*h
        fd=(M.regge_volume_gradient(bd,rad+step,p['Lambda'])[0]-M.regge_volume_gradient(bd,rad-step,p['Lambda'])[0])/(2*h)
        assert abs(fd-g[j])<2e-5
    assert np.linalg.norm(g)>1e-7
    for r in p['scans']:assert np.linalg.norm(r['normal_gradient'])<1e-7
    assert p['scans'][-1]['local_curvature']>p['scans'][1]['local_curvature']>0

def test_stationary_coarse_Schur_response():
    p=packet()['coarse_constraints'];H=np.array(p['unreduced_Hessian']);schur=H[0,0]-H[0,1:]@np.linalg.solve(H[1:,1:],H[1:,0])
    assert abs(schur-p['boundary_Schur_Hessian'])<1e-10
    bd=np.array(p['boundary_lengths']);h=.001;vals=[]
    for x in [h,0,-h]:
        b=bd.copy();b[0,3]+=x;b[3,0]=b[0,3];vals.append(M.normal_stationary(b,p['Lambda'])[0])
    fd=(vals[0]+vals[2]-2*vals[1])/h**2;assert abs(fd-schur)<.05*(1+abs(fd))
    assert np.linalg.norm(H[1:,0])>0

def test_neutrino_relative_phase_lift_and_gauge_control():
    p=packet()['radiative_vacuum'];old=M.previous();B=read(old['pin_vacua']['pin_down']);S=P.read(P.previous()['majorana_vacuum']['Majorana_source']);a=np.array(p['phase_angles'])
    H=np.array(p['phase_Hessian']);G=np.array(p['phase_kinetic_metric']);from scipy.linalg import eigh
    assert min(eigh(H,G,eigvals_only=True))>0
    assert np.max(abs(eigh(H,G,eigvals_only=True)-p['canonical_phase_mass_squared']))<1e-15
    with M.mp.workdps(45):
        center=M.neutrino_phase_potential(B,S,a);h=.001
        for delta in [np.array([h,0]),np.array([0,h])]:assert M.neutrino_phase_potential(B,S,a+delta)>center
        assert abs(float(M.neutrino_phase_potential(B.conjugate(),S,-a)-center))<1e-25
    Q=np.diag(np.exp(1j*np.array([.31,-.23,0])));Y=.01*B/np.sqrt(2);N=np.block([[np.zeros((3,3)),Y],[Y.T,S]]);J=np.block([[Q,np.zeros((3,3))],[np.zeros((3,3)),Q.conjugate()]])
    Yp=Q@Y@Q.conj().T;Sp=Q.conjugate()@S@Q.conj().T;Np=np.block([[np.zeros((3,3)),Yp],[Yp.T,Sp]])
    assert np.linalg.norm(Np-J@N@J.T)<1e-12
    assert np.max(abs(np.linalg.svd(Np,compute_uv=False)-np.linalg.svd(N,compute_uv=False)))<1e-12

def test_link_node_gradient_rank_and_radial_force():
    p=packet()['radiative_vacuum'];D=np.array(P.load()['D'],float);D[40:]*=-1
    assert np.linalg.matrix_rank(D)==79 and 3*79==p['anchored_node_phase_count']
    assert 3*(160-79)==p['remaining_link_cycle_phase_count']
    assert max(p['node_phase_spectral_residuals'])<1e-11
    r=p['radial_link_threshold_scans'];fd=(-3*r[0]['total']+4*r[1]['total']-r[2]['total'])/.001
    assert abs(fd-p['one_sided_radial_force'])<1e-7 and fd<0
    for s in r:
        assert abs(s['link_scalar']-480*sum(x*x*(np.log(x)-1.5)for x in [2*.1*(3*s['phi']**2-.8**2),s['link_phase_mass_squared']]if x>0)/(64*np.pi**2))<1e-10

def class_representatives():
    reps={}
    for x,z in product(range(128),repeat=2):
        c=M.pauli_class(x,z);w=(x|z).bit_count()
        if c not in reps or w<(reps[c][0]|reps[c][1]).bit_count():reps[c]=(x,z)
    return reps

def test_all_flag_faults_against_actual_Steane_states():
    sys.path.insert(0,str(ROOT/'tests'));from test_w33_pass11413_11417_pin_regge_vacuum_noise import steane_basis
    H,C=steane_basis();decoder,single,cases=M.flag_decoder();reps=class_representatives();labels=np.arange(128)
    for e,h,exact in cases:
        residual=e^decoder[h];assert residual==0 if exact else residual in single
        x,z=reps[e];cx,cz=reps[decoder[h]];x^=cx;z^=cz
        # Recover using the actual code's parity-check columns, rather than
        # trusting the producer's syndrome/logical-label classification.
        sx=H@np.array([(x>>j)&1 for j in range(7)])%2;sz=H@np.array([(z>>j)&1 for j in range(7)])%2
        if sx.any():x^=1<<int(np.flatnonzero(np.all(H.T==sx,axis=1))[0])
        if sz.any():z^=1<<int(np.flatnonzero(np.all(H.T==sz,axis=1))[0])
        phase=(-1.)**sum(((labels>>j)&1)*((z>>j)&1)for j in range(7));state=(phase[:,None]*C)[labels^x];assert abs(abs(np.trace(C.conj().T@state))/2-1)<1e-12
    p=packet()['flagged_native_correction'];assert len(cases)==p['single_fault_cases']==1714
    assert len(decoder)==p['history_count']==476


def test_flag_hook_using_nine_qubit_statevector():
    from test_w33_pass11413_11417_pin_regge_vacuum_noise import steane_basis
    _,C=steane_basis();idx=np.arange(512);state=np.zeros(512,complex);state[:128]=(C[:,0]+np.exp(.37j)*C[:,1])/np.sqrt(2)
    def X(v,j):return v[idx^(1<<j)]
    def Z(v,j):return (-1.)**((idx>>j)&1)*v
    def H(v,j):
        mate=v[idx^(1<<j)];return np.where((idx>>j)&1,(mate-v)/np.sqrt(2),(v+mate)/np.sqrt(2))
    def measure(v,j):
        probability=float(np.sum(abs(v[(idx>>j)&1==1])**2));assert min(probability,1-probability)<1e-11
        bit=int(probability>.5);v[((idx>>j)&1)!=bit]=0
        return (X(v,j)if bit else v),bit
    gate=0;sy=[];flags=[]
    def round_checks(flagged):
        nonlocal state,gate
        out=0;fl=0
        for kind in [0,1]:
            for row in range(3):
                if kind:state=H(state,7)
                if not kind and flagged:state=H(state,8)
                dg=[(j,7)if kind==0 else(7,j)for j in range(7)if((j+1)>>row)&1];fg=(8,7)if kind==0 else(7,8);gates=[dg[0],fg,dg[1],dg[2],fg,dg[3]]if flagged else dg
                for c,t in gates:
                    state=state[idx^(((idx>>c)&1)<<t)]
                    if gate==2:state=Z(state,t)
                    gate+=1
                if kind:state=H(state,7)
                if not kind and flagged:state=H(state,8)
                state,b=measure(state,7);out|=b<<(row+3*kind)
                if flagged:state,b=measure(state,8);fl|=b<<(row+3*kind)
        return out,fl
    for _ in range(3):a,b=round_checks(True);sy.append(a);flags.append(b)
    extra=round_checks(False)[0]if any(flags)else -1
    x,z,h,_=M.flagged_circuit(faults={2:8});assert h==tuple(sy+flags+[extra]) and any(flags)
    expected=(C[:,0]+np.exp(.37j)*C[:,1])/np.sqrt(2);labels=np.arange(128)
    for j in range(7):
        if(x>>j)&1:expected=expected[labels^(1<<j)]
        if(z>>j)&1:expected=(-1.)**((labels>>j)&1)*expected
    assert abs(abs(np.vdot(expected,state[:128]))-1)<1e-11


def test_native_CNOT_carrier_and_after_gate_channel():
    p=packet()['flagged_native_correction']['native_CNOT'];U=read(p['native_CNOT_unitary']);F=read(p['logical_frame']);target=read(p['logical_target'])
    assert np.linalg.norm(U.conj().T@U-np.eye(144))<1e-7 and np.linalg.norm(F.conj().T@F-np.eye(4))<1e-10
    error=F.conj().T@U@F@target.conj().T;weights=[]
    X=np.array([[0,1],[1,0]]);Z=np.diag([1,-1])
    for bits in range(16):
        Pa=np.kron(np.linalg.matrix_power(X,bool(bits&1))@np.linalg.matrix_power(Z,bool(bits&2)),np.linalg.matrix_power(X,bool(bits&4))@np.linalg.matrix_power(Z,bool(bits&8)))
        bell=Pa.reshape(-1,order='F')/2;vec=error.reshape(-1,order='F')/2;weights.append(abs(np.vdot(bell,vec))**2)
    assert np.max(abs(np.array(weights)-p['scans'][0]['Pauli_weights']))<1e-12
    assert abs(1-sum(weights)-p['scans'][0]['flagged_leakage'])<1e-12
    assert p['CNOT_total_ticks']==2*p['H_ticks']+2*p['Z_ticks']+p['CZ_ticks']


def test_multifault_intervals_and_leakage_policy():
    p=packet()['flagged_native_correction'];assert p['single_fault_uncorrectable_cases']==0
    for s,g in zip(p['scans'],p['native_CNOT']['scans']):
        k,n=s['uncorrectable_residual_count'],s['samples'];assert k/n==s['uncorrectable_residual_estimate']
        lo,hi=s['Wilson95'];assert 0<=lo<=k/n<=hi<=1
        assert abs(s['worst_case_no_leakage_acceptance']-(1-g['flagged_leakage'])**132)<1e-12
        assert 'accepted' in p['leakage_policy'] and n==16384

if __name__=='__main__':
    tests=[(n,f)for n,f in list(globals().items())if n.startswith('test_')]
    for n,f in tests:f();print(n,'PASS',flush=True)
    print(len(tests),'independent regressions PASS')
