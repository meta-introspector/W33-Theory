"""Independent exact algebra, locality/stability, loop and channel checks."""
import sys,json,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11471_11475_frames_currents_matching as M

def cert():return json.loads(M.OUT.read_text())

def test_canonical_source_binding():
 for path,h in cert()['source_sha256'].items():
  assert hashlib.sha256(json.dumps(json.loads((ROOT/path).read_text()),sort_keys=True,separators=(',',':')).encode()).hexdigest()==h

def test_integer_frame_constraints_and_closure():
 r=cert()['frame'];K=np.array(r['integer_generators'],dtype=np.int64);Z=np.array(r['integer_generator_coefficients']);B,t=M.__dict__['P'].native().N.P.N.M.native_tensors()
 assert np.array_equal(K,np.einsum('ab,bij->aij',Z,B))
 assert not np.any(K[:,:,r['frame_indices']])
 # Rank unchanged when all brackets are adjoined: independent closure replay.
 f=K.reshape(28,-1).T;Q=np.linalg.qr(f)[0]
 for a in K:
  comm=(a@K-K@a).reshape(28,-1).T
  assert np.linalg.norm(comm-Q@(Q.T@comm))<1e-9
 assert np.linalg.matrix_rank(f)==28

def test_radial_connection_integrability():
 r=cert()['current'];assert r['integrability_error']<2e-7
 assert r['admissible_plaquette_angle_bound']<1/30
 replay=M.current();assert abs(replay['curvature']-r['curvature'])<1e-12
 assert sum(q*m for q,m in zip(r['charges'],r['multiplicities']))==0
 assert sum(q**3*m for q,m in zip(r['charges'],r['multiplicities']))==0
 # Exact homotopy identity on an independent nonlinear curvature control.
 import sympy as s
 x,y,t=s.symbols('x y t');F=1+x+x*y+y*y
 I=s.integrate(t*F.subs({x:t*x,y:t*y}),(t,0,1))
 assert s.expand(s.diff(x*I,x)-s.diff(-y*I,y)-F)==0

def test_wave_spectrum_stability_and_energy():
 r=cert()['wave'];A=np.array(r['adjacency']);L=4*np.eye(80)-A;dt=r['dt'];eig=np.linalg.eigvalsh(L)
 assert np.max(eig)<4/dt**2
 q=np.array(r['initial_data']);prev=q.copy();en=[]
 for _ in range(r['steps']):
  nxt=2*q-prev-dt*dt*L@q;en.append(np.sum((nxt-q)**2)/(2*dt**2)+q@L@nxt/2);prev,q=q,nxt
 assert np.max(abs(q-np.array(r['final_data'])))<1e-12 and np.ptp(en)<1e-12
 # A strict spectral stability bound, with a deliberately unstable control.
 assert r['unstable_growth_factor']>1 and r['unstable_dt']>r['CFL_dt_limit']

def test_loop_matching_independent_series():
 import mpmath as mp
 r=cert()['matching'];mp.mp.dps=35;p=mp.mpf(str(r['phi']));mu=mp.mpf(r['renormalization_scale']);Nc=r['color_multiplicity'];n=r['multiplicities']
 f=lambda z:-Nc/(16*mp.pi**2)*sum(nn*m**4*(mp.log(m*m/(mu*mu))-mp.mpf('1.5')) for nn,m in zip(n,[10-mp.sqrt(6)*z,mp.mpf(10),10+mp.sqrt(6)*z]))
 for k,stored in enumerate(r['CW_value_and_first_four_derivatives']):assert abs(float(mp.diff(f,p,k))-stored)<1e-7
 for row in r['polarization']:
  Q=mp.mpf(str(row['Euclidean_momentum']));v=sum(nn*mp.quad(lambda x:x*(1-x)*mp.log(1+Q*Q*x*(1-x)/(m*m)),[0,1]) for nn,m in zip(n,r['masses']))/(4*mp.pi**2)
  assert abs(float(v)-row['subtracted_polarization'])<1e-12
 assert sum(n)==876

def test_decoder_corrects_every_single_pauli_without_twirl():
 V,R=M.recovery_rows();I=np.eye(2);X=np.array([[0,1],[1,0]]);Z=np.diag([1,-1]);Y=1j*X@Z
 for wire in range(7):
  for a in [X,Y,Z]:
   E=np.array([[1.]])
   for j in range(7):E=np.kron(E,a if j==wire else I)
   amplitudes=R@E@V
   # Every recovery amplitude proportional to identity implies identity channel.
   scalar=np.trace(amplitudes,axis1=1,axis2=2)/2
   assert np.linalg.norm(amplitudes-scalar[:,None,None]*I)<1e-12
   assert abs(np.sum(abs(scalar)**2)-1)<1e-12

def test_dissipative_logical_channels_and_relay_dependence():
 r=cert()['decoded'];rows=r['rows'];omega=np.array([1,0,0,1])/np.sqrt(2)
 for row in rows:
  C=M.P.dec(row['logical_choi']);assert min(np.linalg.eigvalsh(C))>-1e-12
  assert np.linalg.norm(np.einsum('iaja->ij',C.reshape(2,2,2,2))-np.eye(2)/2)<1e-12
  assert abs(1-(omega@C@omega).real-row['logical_entanglement_infidelity'])<1e-12
 for i in range(0,len(rows),2):assert np.linalg.norm(M.P.dec(rows[i]['logical_choi'])-M.P.dec(rows[i+1]['logical_choi']))>1e-5
 # Realistic route noise can exceed the useful range of the ideal decoder.
 assert rows[-1]['logical_entanglement_infidelity']>rows[-1]['physical_entanglement_infidelity']

def test_native_mask_controls_and_triality_dictionary():
 r=cert();f=r['frame'];K=np.array(f['integer_generators']);blocks=f['invariant_coordinate_blocks']
 assert sorted(map(len,blocks))==[1,1,1,8,8,8]
 for b in blocks:
  outside=sorted(set(range(27))-set(b));assert not np.any(K[:,b][:,:,outside])
 assert sorted(t['count'] for t in f['cubic_block_patterns'])==[1,4,4,4,32]
 c=r['decoded']['controls'];assert sorted(v for b in c['edge_phase_masks'] for v in b)==list(range(160))
 assert c['cell_invariance_error']<1e-10 and c['maximum_logical_leakage_commutator']>.5
 E=c['diagonal_energies'];g=[E[j]-E[i] for i,j in c['connected_transition_edges']];assert len(g)==len(set(g))
 seen={0}
 for _ in range(12):
  for i,j in c['connected_transition_edges']:
   if i in seen or j in seen:seen.update([i,j])
 assert len(seen)==12

 # Independently replay the physical160-mode masks, not just compressed flags.
 N=M.P.native().N.P.N.M;P=N.P;B=P.read(N.previous()['native_noise']['invariant_frame']);Q=B.conj().T@P.read(P.prior()['local_gates']['logical_basis'])
 logical=Q@Q.conj().T;mix=[]
 for b in c['edge_phase_masks']:
  mask=np.zeros(160);mask[b]=1;A=B.conj().T@(mask[:,None]*B)
  assert np.linalg.norm(mask[:,None]*B-B@A)<1e-10
  mix.append(np.linalg.norm(A@logical-logical@A))
 assert max(mix)>.5


def test_exact_native_to_clock_albert_intertwiner():
 r=cert()['frame'];b=r['albert_bridge'];K=np.array(r['integer_generators'],dtype=np.int64);T=M.P.dec(b['native_from_clock_basis']);Ti2=M.P.dec(b['twice_inverse']);Pi=np.array(b['twice_Jordan_product'],dtype=np.int64);R=M.P.dec(b['twice_conjugated_generators'])
 assert np.array_equal(Ti2@T,2*np.eye(27))
 assert np.array_equal(R,np.array([Ti2@a@T for a in K]))
 for C in R:
  for A in [C.real.astype(np.int64),C.imag.astype(np.int64)]:
   left=np.einsum('ka,ija->ijk',A,Pi)
   right=np.einsum('ai,ajk->ijk',A,Pi)+np.einsum('aj,iak->ijk',A,Pi)
   assert np.array_equal(left,right)
 assert not np.any(R[:,:,b['clock_idempotent_indices']])
 # Verify the stored product really is the prior clock object, rather than
 # a product invented to accommodate the new generators.
 import w33_pass10950_clock_albert_lorentz_spinor as C
 J=C.build_clock_albert();assert np.array_equal(np.array(J['prod'],float)*2,Pi)

 # Independent determinant reconstruction from trace powers of actual Jordan
 # left multiplication, with exact denominator clearing.
 B,tensors=M.P.native().N.P.N.M.native_tensors();signed=tensors[1]
 tr=np.einsum('ijj->i',Pi);assert not np.any(tr%18);tr=tr//18
 gram=np.einsum('k,ijk->ij',tr,Pi);cube=np.einsum('l,ial,jka->ijk',tr,Pi,Pi)
 N12=2*np.einsum('i,j,k->ijk',tr,tr,tr)-(np.einsum('i,jk->ijk',tr,gram)+np.einsum('j,ik->ijk',tr,gram)+np.einsum('k,ij->ijk',tr,gram))+cube
 native12=2*np.einsum('abc,ai,bj,ck->ijk',signed,T,T,T,optimize=True)
 assert np.array_equal(native12,N12) and b['cubic_coefficient_residual']==0
