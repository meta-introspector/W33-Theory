"""Independent field equations, exchange symmetry, Ward and flagged-channel tests."""
import sys,json,hashlib
from pathlib import Path
import numpy as np
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11481_11485_fibers_metric_yukawa_noisy as M
cert=lambda:json.loads(M.OUT.read_text())
def test_source_binding_and_section_scope():
 c=cert()
 for p,h in c['source_sha256'].items():assert hashlib.sha256(json.dumps(M.read(p),sort_keys=True,separators=(',',':')).encode()).hexdigest()==h
 for key in ['fibers','metric','yukawa','ward','recovery']:assert c[key]['status']=='PASS' and len(c[key]['scope'])>150

def test_finite_fibers_on_actual_native_potential():
 import w33_pass11438_finite_native_model as F
 old=M.read(M.OLD)['vacuum'];p=old['polydisc'];A=M.dec(p['A_at_vacuum']);B=M.dec(p['B_at_vacuum']);E=M.dec(old['isotope']['native_tripotents']);H=M.dec(p['compressed_moment_generators'])
 import w33_pass11476_11480_cubic_geometry_correlated as P
 base=P.polydisc_value(A,B,H)
 with F.native_context():
  for row in cert()['fibers']['rows']:
   U=M.dec(row['unitary']);a=U@A;b=U@B
   assert np.linalg.norm(U.conj().T@U-np.eye(3))<1e-12 and abs(np.linalg.det(U)-1)<1e-12
   for C in [A,B]:assert np.linalg.norm(np.diag(U@C@C.conj().T@U.conj().T-C@C.conj().T))<1e-10
   x=F.pack((E@a).reshape(-1),(E@b).reshape(-1));v,g=F.fun(x)
   assert abs(v-base)<1e-9 and abs(np.linalg.norm(g)-row['native_gradient_norm'])<1e-10
 assert cert()['fibers']['physical_kernel_singular_values'][2]<1e-9
 assert max(z['native_gradient_norm'] for z in cert()['fibers']['rows'])>1e-9

def test_metric_hamilton_equations_by_independent_energy_derivatives():
 r=cert()['metric'];y=np.array(r['initial_canonical_state']);edges=np.array(r['edges']);h=np.array(r['harmonic_displacements']);energy,flow=M.metric_hamiltonian(y,edges,h);rng=np.random.default_rng(98)
 # Frobenius matrix pairing accounts for off-diagonal symmetric coordinates.
 for offset,partner,sign in [(0,9,-1),(9,0,1)]:
  d=rng.normal(size=(3,3));d=(d+d.T)/2;delta=np.zeros(len(y));delta[offset:offset+9]=d.ravel();step=1e-6
  fd=(M.metric_hamiltonian(y+step*delta,edges,h)[0]-M.metric_hamiltonian(y-step*delta,edges,h)[0])/(2*step)
  assert abs(fd-sign*np.dot(flow[partner:partner+9],d.ravel()))<1e-7
 for offset,partner,sign in [(18,98,-1),(98,18,1)]:
  d=rng.normal(size=80);delta=np.zeros(len(y));delta[offset:offset+80]=d;step=1e-7
  fd=(M.metric_hamiltonian(y+step*delta,edges,h)[0]-M.metric_hamiltonian(y-step*delta,edges,h)[0])/(2*step)
  assert abs(fd-sign*np.dot(flow[partner:partner+80],d))<1e-6
 for row in r['Hamiltonian_history']:
  assert abs(M.metric_hamiltonian(np.array(row['state']),edges,h)[0])<1e-7
 assert np.linalg.norm(np.array(r['Hamiltonian_history'][-1]['state'])[:9]-y[:9])>1e-4

def test_metric_covariance_on_an_independent_basis():
 r=cert()['metric'];G=np.array(r['metric']);dg=np.array(r['constrained_metric_velocity']);u=np.array(r['scalar']);du=np.array(r['scalar_velocity']);edges=np.array(r['edges']);h=np.array(r['harmonic_displacements']);J=np.array([[1,.2,.1],[0,1.3,0],[.1,0,.8]]);Ji=np.linalg.inv(J)
 f=M.metric_action(G,dg,.8,u,du,edges,h)[0]
 assert abs(f-M.metric_action(Ji.T@G@Ji,Ji.T@dg@Ji,.8,u,du,edges,h@J.T,cell=np.linalg.det(J))[0])<1e-10
 assert r['positive_kinetic_lapse_force']<0 and r['lapse_constraint_residual']<1e-10

def test_yukawa_exchange_and_native_mass_map():
 from w33_pass11384_11388_native_dynamics import native_tensors
 d=np.array(native_tensors()[1][1]);r=cert()['yukawa'];eps=np.array(r['family_epsilon']);H=M.dec(r['enlarged_Higgs']);Y=M.dec(r['allowed_mass_matrix'])
 assert np.array_equal(d,d.transpose(1,0,2)) and np.array_equal(eps,-eps.transpose(1,0,2))
 # Rebuild by explicit 3x3 blocks, without producer combined-index einsum.
 blocks=np.block([[sum(d[:,:,k]*H[k,i,j] for k in range(27)) for j in range(3)] for i in range(3)])
 order=np.arange(81).reshape(3,27).T.ravel();assert np.linalg.norm(blocks[np.ix_(order,order)]-Y)<1e-12
 assert np.linalg.norm(Y-Y.T)<1e-12
 assert np.linalg.norm(np.linalg.svd(Y,compute_uv=False)-r['allowed_singular_values'])<1e-10
 assert 'Pass11271 owns' in r['scope']

def test_actual_anomaly_derivative_and_restricted_flows():
 import w33_pass11428_11432_native_alignment_measure_coarse_flags as Q
 import w33_pass11476_11480_cubic_geometry_correlated as P
 r=cert()['ward'];A=np.array(r['background']);J=np.array(r['anomaly_Jacobian']);B=P.gauge_boundary(2);d=np.random.default_rng(72).normal(size=64);direct=np.zeros(16)
 for q,m in zip([1,-4,2,-3,6],[6,3,3,2,1]):
  _,ds,*_=Q.weak_overlap(q,1.,A.reshape(16,4),d.reshape(1,16,4),L=2)
  direct+=m*q*np.trace(ds[0].reshape(16,4,16,4),axis1=1,axis2=3).diagonal().real
 assert np.linalg.norm(direct-J@d)<1e-11
 for row in r['rows']:
  flow=np.array(row['flows']);assert abs(max(np.linalg.norm(B.T@flow+J,axis=0))-row['maximum_divergence_error'])<1e-14
 assert np.linalg.norm(J@B)<1e-10
 assert r['rows'][0]['maximum_divergence_error']>=r['rows'][-1]['maximum_divergence_error']

def test_noisy_readout_erasure_channels_and_MAP_objective():
 r=cert()['recovery'];priors=np.array(r['syndrome_priors']);bits=np.array([[(s>>b)&1 for b in range(6)] for s in range(64)]);dist=np.sum(bits[:,None,:]!=bits[None,:,:],axis=2)
 for row in r['rows']:
  error=row['readout_bit_error'];C=error**dist*(1-error)**(6-dist);guess=np.array(row['guess']);success=sum(C[t,guess[t]]*priors[guess[t]] for t in range(64));J=M.dec(row['flagged_Choi'])
  assert abs(success-row['success_probability'])<1e-10
  assert np.linalg.norm(np.einsum('iaja->ij',J.reshape(2,3,2,3))-np.eye(2)/2)<1e-10
  assert min(np.linalg.eigvalsh(J))>-1e-10
  if row['adaptive']:assert abs(success-np.max(C*priors[None,:],axis=1).sum())<1e-10
 for e in [.001,.01,.05]:
  rows=[z for z in r['rows'] if z['readout_bit_error']==e];assert rows[1]['success_probability']>=rows[0]['success_probability']-1e-12

def test_native_leakage_cell_and_echo():
 import w33_pass11471_11475_frames_currents_matching as Q
 n=Q.P.native();old=n.N.P.N.M.P;B=old.read(n.N.P.N.M.previous()['native_noise']['invariant_frame']);logical=B.conj().T@old.read(old.prior()['local_gates']['logical_basis']);r=cert()['recovery']['native_leakage'];mask=np.isin(np.arange(160),r['block']);h=B.conj().T@(mask[:,None]*B);U=expm(1j*r['phase']*h);P=logical@logical.conj().T
 assert np.linalg.norm(U-M.dec(r['cell_unitary']))<1e-12
 assert abs(np.linalg.norm((np.eye(12)-P)@U@logical)**2/2-r['mean_leakage'])<1e-12
 assert np.linalg.norm(U.conj().T@U@logical-logical)<1e-12

def test_exact_generic_rank_control():
 import sympy as sp
 r=cert()['fibers']['exact_generic_control'];C=sp.Matrix([[sp.sympify(x) for x in row] for row in r['constraint_Jacobian']]);N=sp.Matrix([[sp.sympify(x) for x in row] for row in r['kernel']])
 assert C.rank()==4 and N.rank()==4 and C*N==sp.zeros(4,4)
 assert M.exact_fiber_control()['exact_ranks']==[4,10,12]

def test_portal_loop_feedback_on_finite_fibers():
 c=cert();r=c['yukawa'];d=M.dec(r['native_cubic']);S=M.dec(r['spurion']);old=M.read(M.OLD)['vacuum'];E=M.dec(old['isotope']['native_tripotents']);A=M.dec(old['polydisc']['A_at_vacuum']);mu=10
 for row,loop in zip(c['fibers']['rows'],r['fiber_portal_rows']):
  field=E@M.dec(row['unitary'])@A;H=np.einsum('ck,kij->cij',field,S);Y=np.einsum('abc,cij->aibj',d,H).reshape(81,81);mass=np.linalg.svd(Y,compute_uv=False);mass=mass[mass>1e-12];value=-sum(x**4*(np.log(x*x/mu**2)-1.5) for x in mass)/(32*np.pi**2)
  assert abs(value-r['supplied_Weyl_CW_value']-loop['Weyl_loop_shift'])<1e-12

def test_soft_restoring_decomposition_against_full_native_gradient():
 import w33_pass11438_finite_native_model as F
 old=M.read(M.OLD)['vacuum'];p=old['polydisc'];A=M.dec(p['A_at_vacuum']);B=M.dec(p['B_at_vacuum']);E=M.dec(old['isotope']['native_tripotents']);r=cert()['fibers']['restoring_decomposition'];gram=np.array(r['Gram_term']);H=np.array(r['generator_path_Hessian']);w,V=np.linalg.eigh(gram)
 assert sum(w>1e-10)==4 and sum(w>1e-3)==2
 assert np.linalg.norm(H-gram-np.array(r['moment_curvature_term']))<1e-12
 def derivative(t,L):
  U=expm(1j*t*L);a=U@A;b=U@B;x=F.pack((E@a).reshape(-1),(E@b).reshape(-1));dx=F.pack((E@(1j*L@a)).reshape(-1),(E@(1j*L@b)).reshape(-1));return F.fun(x)[1]@dx
 with F.native_context():
  for index in [4,5,6,7]:
   c=V[:,index];L=np.einsum('a,aij->ij',c,M.su3());h=1e-4;fd=(derivative(h,L)-derivative(-h,L))/(2*h);expected=c@H@c
   assert abs(fd-expected)<(1e-9 if index<6 else 1e-7)
