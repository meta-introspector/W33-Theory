"""Independent matching, geometry, action, spectrum and QEC checks."""
import hashlib,json,sys
from pathlib import Path
from itertools import product,combinations
import numpy as np
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
from w33_pass11390_11394_context_matter_clock import load
def read(x):return np.array(x['real'])+1j*np.array(x['imag'])
def packet():return json.loads((ROOT/'data/w33_pass11413_11417_pin_regge_vacuum_noise.json').read_text())

def test_source_hashes():
    p=packet()
    for s,h in p['source_sha256'].items():
        value=json.loads((ROOT/s).read_text());canonical=json.dumps(value,sort_keys=True,separators=(',',':')).encode()
        assert hashlib.sha256(canonical).hexdigest()==h
    assert all(p[k]['status']=='PASS'for k in ['path_flavour','majorana_vacuum','regge_patch','vacuum_inventory','counted_noise'])

def test_full_messenger_matching_and_rank():
    f=packet()['path_flavour'];A=load()['A'];n=len(A);P=np.zeros_like(A);P[0,0]=1
    row=f['scans'][1];K=np.eye(n)-row['eta']*A;B=read(f['pin_down']);k=np.kron(K,np.eye(3))
    M=np.block([[k,np.kron(P,B)],[np.zeros_like(k),k]])
    src=np.zeros((2*n*3,3));dst=np.zeros((3,2*n*3))
    for j,i in enumerate(f['endpoints']):src[n*3+i*3+j,j]=1;dst[j,i*3+j]=1
    Y=-dst@np.linalg.solve(M,src)
    assert np.linalg.norm(Y-read(row['Yd']))<1e-12 and np.linalg.matrix_rank(Y)==3
    # Replacing the family fibre by a scalar common pin reproduces old rank1.
    r=np.linalg.inv(K)[f['endpoints'],0];assert np.linalg.matrix_rank(np.outer(r,r))==1

def test_path_selection_coefficients_and_colour_safe_basis():
    f=packet()['path_flavour'];A=load()['A']
    for i,d in zip(f['endpoints'],f['distances']):
        for n in range(d):assert np.linalg.matrix_power(A,n)[i,0]==0
        assert np.linalg.matrix_power(A,d)[i,0]==1
    B=read(f['pin_down']);tail=np.linalg.det(B[1:,1:]);lead=[np.linalg.det(B)/tail,tail/B[2,2],B[2,2]]
    assert np.linalg.norm(np.real(lead)-f['asymptotic_mass_coefficients'])<1e-12
    row=f['scans'][-1];eta=row['eta']
    assert np.max(abs(np.array(row['down_masses'])/eta**np.array([6,4,0])/np.real(lead)-1))<.002
    # Re-express family fibres in a second orthonormal basis; eigenvalues agree.
    U=np.linalg.qr(np.random.default_rng(12).normal(size=(3,3))+1j*np.eye(3))[0]
    Y=read(row['Yd']);assert np.linalg.norm(np.linalg.eigvalsh(U@Y@U.conj().T)-row['down_masses'])<1e-12

def test_mixing_powers_and_CP_triangle_identity():
    f=packet()['path_flavour']
    for row in f['scans']:
        Y=read(row['Yd']);w,V=np.linalg.eigh(Y)
        J=np.imag(V[0,0]*V[1,1]*V[0,1].conjugate()*V[1,0].conjugate())
        VT=V.T;JT=np.imag(VT[0,0]*VT[1,1]*VT[0,1].conjugate()*VT[1,0].conjugate())
        assert abs(J-JT)<1e-15  # Parallel11371: flavour CP is transpose-even.
        tri=np.imag(Y[0,1]*Y[1,2]*Y[2,0]);den=np.prod([w[j]-w[i]for i,j in combinations(range(3),2)])
        assert abs(abs(J)-abs(tri/den))<1e-12 and abs(J)>0
        _,Vc=np.linalg.eigh(Y.conjugate());Jc=np.imag(Vc[0,0]*Vc[1,1]*Vc[0,1].conjugate()*Vc[1,0].conjugate())
        assert abs(J+Jc)<1e-12
    r=f['scans'][-1];assert max(abs(np.array(r['angles_over_eta_powers'])/f['asymptotic_angle_coefficients']-1))<.01

def test_exact_radial_recurrence_and_critical_hierarchy_loss():
    import sympy as sp
    t=sp.symbols('t');den=(1-16*t*t)*(1-6*t*t)
    r=[(1-18*t*t+36*t**4)/den,t*(1-15*t*t)/den,t*t*(1-12*t*t)/den,t**3/den,4*t**4/den]
    equations=[r[0]-4*t*r[1]-1,r[1]-t*(r[0]+3*r[2]),r[2]-t*(r[1]+3*r[3]),r[3]-t*(r[2]+3*r[4]),r[4]-4*t*r[3]]
    assert all(sp.factor(e)==0 for e in equations)
    assert all(sp.limit(x/r[0],t,sp.Rational(1,4),dir='-')==1 for x in r)
    f=packet()['path_flavour'];assert f['distance_shell_sizes']==[1,4,12,36,27]
    assert min(f['critical_scans'][-1]['endpoint_to_pin_ratios'])>.95
    assert min(f['critical_scans'][-2]['endpoint_to_pin_ratios'])<.1

def test_majorana_polynomial_Hessian_and_minimum():
    m=packet()['majorana_vacuum'];S0=read(m['Majorana_source']);E=[read(b)for b in m['symmetric_field_basis']];p=m['parameters']
    def V(x):
        delta=sum(a*b for a,b in zip(x,E));S=S0+delta;n=np.linalg.norm(delta)**2
        return (p['t']+p['kappa']*p['v']**2/2)*n+p['lambda_S']*np.linalg.norm(S@S.conj().T-S0@S0.conj().T)**2
    h=1e-5;A=np.eye(12)*h;H=np.zeros((12,12));z=np.zeros(12)
    for i in range(12):
        H[i,i]=(V(A[i])+V(-A[i])-2*V(z))/h**2
        for j in range(i):H[i,j]=H[j,i]=(V(A[i]+A[j])-V(A[i]-A[j])-V(-A[i]+A[j])+V(-A[i]-A[j]))/(4*h*h)
    assert np.linalg.norm(H-np.array(m['scalar_Hessian']))<1e-7 and min(np.linalg.eigvalsh(H))>0
    for x in np.random.default_rng(13).normal(size=(20,12)):assert V(x)>0

def test_full_seesaw_and_balanced_chirality():
    m=packet()['majorana_vacuum'];N=read(m['seesaw_matrix']);Y=read(m['neutrino_Dirac']);M=read(m['Majorana_source'])
    assert np.linalg.norm(N-N.T)<1e-12 and np.linalg.matrix_rank(N,tol=1e-18)==6
    expected=-np.linalg.det(Y/np.sqrt(2))**2
    assert abs(np.linalg.det(N)-expected)<3e-14*abs(expected)
    approx=-(Y/np.sqrt(2))@np.linalg.inv(M)@(Y.T/np.sqrt(2))
    assert np.max(abs(np.sort(np.linalg.svd(approx,compute_uv=False))/np.sort(m['neutrino_masses'])[:3]-1))<.001
    D=np.block([[np.zeros_like(N),N],[N.conj().T,np.zeros_like(N)]])
    assert np.linalg.norm(np.linalg.eigvalsh(D)+np.linalg.eigvalsh(D)[::-1])<1e-12
    assert m['balanced_grading_index']==0

def test_Regge_angles_against_cartesian_normals():
    from w33_pass11413_11417_pin_regge_vacuum_noise import simplex_geometry
    p=packet()['regge_patch'];V=np.vstack([p['interior_vertex'],p['boundary_vertices']]);L=np.linalg.norm(V[:,None]-V[None,:],axis=2)
    for omit in range(1,6):
        inds=[0]+[i for i in range(1,6)if i!=omit];coords=V[inds];normals=[]
        for a in range(5):
            face=np.delete(coords,a,axis=0);_,_,vh=np.linalg.svd(face[1:]-face[0]);normal=vh[-1]
            if np.dot(normal,coords[a]-face[0])>0:normal=-normal
            normals.append(normal)
        angles=simplex_geometry(inds,L)
        for a,b in combinations(range(5),2):
            tri=tuple(sorted(set(inds)-{inds[a],inds[b]}));independent=np.arccos(np.clip(-np.dot(normals[a],normals[b]),-1,1))
            assert abs(angles[tri]-independent)<1e-12

def test_Regge_nonlinear_flat_branch_and_Hessian_convergence():
    from w33_pass11413_11417_pin_regge_vacuum_noise import regge_action,hessian
    p=packet()['regge_patch'];V=np.array(p['boundary_vertices']);b=np.array(p['boundary_lengths']);c=np.array(p['interior_vertex']);r=np.array(p['radial_lengths'])
    for dc in np.random.default_rng(14).normal(size=(20,4))*.02:
        s,d=regge_action(np.linalg.norm(V-c-dc,axis=1),b)
        assert abs(s-p['flat_action'])<1e-10 and max(abs(d))<1e-12
    fn=lambda r:regge_action(r,b)[0];H=hessian(fn,r,h=3e-5);T=np.array(p['vertex_translation_tangents'])
    coarse=hessian(fn,r,h=8e-5)
    assert np.linalg.norm(H@T)<np.linalg.norm(coarse@T)
    assert np.linalg.norm(H@T)<.002 and np.linalg.norm(H-np.array(p['flat_Hessian']))<.002
    assert np.count_nonzero(abs(np.linalg.eigvalsh(H))>.02)==1
    assert p['curved_offshell_scans'][-1]['max_deficit']>1e-4

def test_complete_declared_vacuum_inventory_and_scale():
    full=packet();p=full['vacuum_inventory'];rows=p['inventory'];ferm=[r for r in rows if r['spin']in ['Dirac','Majorana']]
    lookup={r['name']:r for r in rows};v=full['majorana_vacuum']['parameters']['v'];f=full['path_flavour']['scans'][1]
    for name,Y in [('up',read(f['Yu'])),('down',read(f['Yd'])),('charged_lepton',.03*np.diag([.01,.1,1.]))]:
        expected=np.linalg.svd(v*Y/np.sqrt(2),compute_uv=False)**2
        assert np.linalg.norm(expected-[lookup[f'{name}{i}']['m2']for i in range(3)])<1e-12
    # Gauge masses from the doublet kinetic term, not a stored formula.
    g,gp=p['input_gauge_couplings'];h=np.array([0,v/np.sqrt(2)])
    T=[g*np.array([[0,1],[1,0]])/2,g*np.array([[0,-1j],[1j,0]])/2,g*np.diag([1,-1])/2,gp*np.eye(2)/2]
    gram=2*np.real(np.array([[np.vdot(a@h,b@h)for b in T]for a in T]))
    eigen=np.linalg.eigvalsh(gram)
    assert np.linalg.norm(eigen-[0,lookup['W_pair']['m2'],lookup['W_pair']['m2'],lookup['Z']['m2']])<1e-12
    assert sum(-r['degeneracy']for r in ferm)==96
    assert len([r for r in rows if r['spin']=='real_scalar'])==13
    assert sum(r['degeneracy']for r in rows if r['spin']=='vector')==9
    str4=sum(r['degeneracy']*r['m2']**2 for r in rows)
    value=lambda mu:sum(r['degeneracy']*r['m2']**2*(np.log(r['m2']/mu**2)-r['c'])for r in rows)/(64*np.pi**2)
    assert abs(str4-p['supertrace_M4'])<1e-10 and abs(value(1)-p['one_loop_value'])<1e-12
    assert abs((value(np.exp(.001))-value(np.exp(-.001)))/.002+str4/(32*np.pi**2))<1e-9
    assert abs(str4)>1

def test_principal_phase_native_counted_gate_replay():
    p=packet()['counted_noise'];old=json.loads((ROOT/'data/w33_pass11403_11407_native_bimodule_collective_gates.json').read_text())['local_gates']
    Q=read(old['logical_basis']);target=read(old['ideal_qubit_gate']);D=np.array(load()['D'],float)
    clock=expm(-.02j*D[:40].T@D[:40])@expm(-.04j*D[40:].T@D[40:])@expm(-.02j*D[:40].T@D[:40]);a,b=old['gate_pair']
    state=Q.copy()
    assert max(abs(np.exp(1j*np.array(p['principal_angles']))-np.exp(1j*np.array(p['original_angles']))))<1e-14
    for angle,e,N in zip(p['principal_angles'],[a,b,a,b,a],p['ticks_per_pulse']):
        kick=np.eye(160,dtype=complex);kick[e,e]=np.exp(-.5j*.04*np.copysign(.0125,angle));state=np.linalg.matrix_power(kick@clock@kick,N)@state
    err=np.linalg.norm(state-Q@target,ord=2)
    assert abs(err-p['coherent_encoded_error'])<1e-8 and err<p['coherent_error_bound']
    assert p['total_ticks']<p['previous_ticks']/4

def steane_basis():
    H=np.array([[((j>>i)&1)for j in range(1,8)]for i in range(3)])
    words=[tuple(np.array(x)@H%2)for x in product([0,1],repeat=3)];C=np.zeros((128,2),complex)
    for word in words:
        i=sum(a<<j for j,a in enumerate(word));C[i,0]=1/np.sqrt(8);C[i^127,1]=1/np.sqrt(8)
    return H,C

def test_Steane_quantum_single_error_Knill_Laflamme():
    _,C=steane_basis();states=[C];labels=np.arange(128)
    for j in range(7):
        z=1-2*((labels>>j)&1);states.extend([C[labels^(1<<j)],z[:,None]*C,1j*z[:,None]*C[labels^(1<<j)]])
    for i,X in enumerate(states):
        for j,Y in enumerate(states):
            E=X.conj().T@Y;assert np.linalg.norm(E-(np.eye(2)if i==j else np.zeros((2,2))))<1e-12

def test_Steane_exhaustive_dephasing_recovery():
    H,C=steane_basis();labels=np.arange(128);counts=np.zeros(8,int)
    for e in product([0,1],repeat=7):
        e=np.array(e);syndrome=H@e%2;flip=np.zeros(7,int)
        if syndrome.any():flip[np.flatnonzero(np.all(H.T==syndrome,axis=1))[0]]=1
        residual=e^flip;phase=(-1.)**sum(((labels>>j)&1)*residual[j]for j in range(7));logical=C.conj().T@(phase[:,None]*C)
        if abs(logical[0,0]-logical[1,1])>1:counts[e.sum()]+=1
    p=packet()['counted_noise'];assert counts.tolist()==p['Steane_failure_weight_counts']
    for row in p['noise_scans']:
        q=row['effective_Z_probability'];f=sum(n*q**w*(1-q)**(7-w)for w,n in enumerate(counts))
        assert abs(f-row['encoded_logical_Z_probability'])<1e-15 and f<q

if __name__=='__main__':
    tests=[(n,f)for n,f in list(globals().items())if n.startswith('test_')]
    for n,f in tests:f();print(n,'PASS',flush=True)
    print(len(tests),'independent regressions PASS')
