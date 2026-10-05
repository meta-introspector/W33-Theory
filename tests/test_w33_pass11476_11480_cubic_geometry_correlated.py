"""Independent carrier identities, native derivatives and channel controls."""
import sys,json,hashlib
from pathlib import Path
from itertools import product
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11476_11480_cubic_geometry_correlated as M

def cert():return json.loads(M.OUT.read_text())

def test_source_bindings():
 for path,h in cert()['source_sha256'].items():
  assert hashlib.sha256(json.dumps(json.loads((ROOT/path).read_text()),sort_keys=True,separators=(',',':')).encode()).hexdigest()==h

def test_isotope_product_and_conjugation_from_independent_operator_formula():
 r=cert()['vacuum']['isotope'];c=M.dec(r['complete_tripotent']);p=np.array(M.read('data/w33_pass11471_11475_frames_currents_matching.json')['frame']['albert_bridge']['twice_Jordan_product'])/2
 L=lambda x:np.einsum('i,ijk->kj',x,p)
 Lbar=L(c.conj());pc=np.einsum('ai,ajk->ijk',Lbar,p)+np.einsum('aj,aik->ijk',Lbar,p)-np.einsum('ka,ija->ijk',Lbar,p)
 assert np.linalg.norm(pc-M.dec(r['isotope_product']))<1e-12
 cc=np.einsum('i,j,ijk->k',c,c,p);star=2*L(c)@L(c)-L(cc)
 assert np.linalg.norm(star-M.dec(r['isotope_conjugation']))<1e-12
 assert np.linalg.norm(star@star.conj()-np.eye(27))<1e-12
 assert np.linalg.norm(np.einsum('ka,ija->ijk',star,pc.conj())-np.einsum('ai,bj,abk->ijk',star,star,pc))<1e-12
 trace=np.einsum('ijj->i',pc)/9;metric=np.einsum('k,ajk,ai->ij',trace,pc,star)
 assert min(np.linalg.eigvalsh(metric))>.9
 E=M.dec(r['primitive_tripotents'])
 for i,j in product(range(3),repeat=2):
  z=np.einsum('a,b,abk->k',E[:,i],E[:,j],pc)
  assert np.linalg.norm(z-(E[:,i] if i==j else 0))<1e-12
 assert np.linalg.norm(np.einsum('i,ijk->kj',c,pc)-np.eye(27))<1e-12
 for A in M.dec(r['native_conjugated_generators']):
  assert np.linalg.norm(A@c)<1e-11
  assert np.linalg.norm(np.einsum('ka,ija->ijk',A,pc)-np.einsum('ai,ajk->ijk',A,pc)-np.einsum('aj,iak->ijk',A,pc))<1e-11
 # Negative control: the original binary product did not close this space.
 assert cert()['vacuum']['Jordan_closure_residual']>.5

def test_isotope_cubic_against_original_signed_native_triads():
 from w33_pass11384_11388_native_dynamics import native_tensors
 r=cert()['vacuum']['isotope'];p=M.dec(r['isotope_product']);t=np.einsum('ijj->i',p)/9
 b=M.read('data/w33_pass11471_11475_frames_currents_matching.json')['frame']['albert_bridge'];T=M.dec(b['native_from_clock_basis']);tri=native_tensors()[1][0]
 rng=np.random.default_rng(11480)
 for y in rng.normal(size=(8,27))+1j*rng.normal(size=(8,27)):
  native=T@y;N=sum(s*native[a]*native[b]*native[c] for a,b,c,s in tri)
  yy=np.einsum('i,j,ijk->k',y,y,p);yyy=np.einsum('i,j,ijk->k',y,yy,p)
  new=((t@y)**3-3*(t@y)*(t@yy)+2*t@yyy)/6
  assert abs(new-N)<1e-10

def test_mixed_sextic_derivatives_at_a_third_step():
 import w33_pass11438_finite_native_model as F
 old=M.read('data/w33_pass11466_11470_stationary_ports_channels.json');x=np.array(old['vacuum']['coordinates']);S=M.dec(old['soft']['soft_frame']).real;r=cert()['vacuum']
 p=x[:81]+1j*x[81:162];s=x[162:243]+1j*x[243:];dp=S[:81]+1j*S[81:162];ds=S[162:243]+1j*S[243:];expected=M.dec(r['analytic_Jacobian']);h=2.5e-5;actual=[]
 with F.native_context():
  for a in M.dec(r['mixed_sextic_probes']):
   actual.append([(F.D.I.evaluate(p+a*s+h*(dp[:,j]+a*ds[:,j]))[0]-F.D.I.evaluate(p+a*s-h*(dp[:,j]+a*ds[:,j]))[0])/(2*h) for j in range(4)])
 assert np.linalg.norm(np.array(actual)-expected)<1e-8
 assert np.linalg.norm(np.array(actual)-expected)<r['finite_difference_errors'][0]

def test_native_polydisc_reduction_on_independent_field_samples():
 import w33_pass11438_finite_native_model as F
 r=cert()['vacuum']['polydisc'];E=M.dec(r['native_embedding']);H=M.dec(r['compressed_moment_generators']);rng=np.random.default_rng(37)
 assert r['exact_canonical_constants']==[27,729,0]
 assert np.linalg.norm(E.conj().T@E-np.eye(9))<1e-12
 with F.native_context():
  for size in [.001,.007]:
   A=M.dec(r['A_at_vacuum'])+size*(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
   B=M.dec(r['B_at_vacuum'])+size*(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
   x=F.pack(E@A.reshape(-1),E@B.reshape(-1))
   assert abs(F.fun(x)[0]-M.polydisc_value(A,B,H))<1e-9
   a,b,*_=F.D.I.evaluate(E@A.reshape(-1));c,_=F.D.evaluate(E@A.reshape(-1));d=27*np.linalg.det(A)**2
   assert abs(a-d)<1e-10 and abs(b-d*d)<1e-10 and abs(c)<1e-10

def test_determinant_signature_and_cover_map():
 b,p,n=M.determinant_tensor();r=cert()['geometry'];ids=b['clock_idempotent_indices']
 # Independent central differences of cubic gradients recover the Hessian.
 grad=lambda y:3*np.einsum('ijk,j,k->i',n,y,y)
 for row in r['frame_controls']:
  y=np.zeros(27);y[ids]=row['frame_values'];h=1e-4
  H=np.column_stack([(grad(y+h*d)-grad(y-h*d))/(2*h) for d in np.eye(27)])
  assert np.linalg.norm(H-np.array(row['Hessian']))<1e-10
 assert r['frame_controls'][0]['inertia']==[1,26,0]
 assert r['frame_controls'][0]['log_Hessian_inertia']==[0,27]
 assert r['frame_controls'][2]['inertia']==[10,17,0]
 assert r['frame_controls'][3]['inertia'][2]==9
 E=np.array(r['selected_4D_embedding']);H=np.array(r['frame_controls'][0]['Hessian'])
 assert np.linalg.norm(E.T@H@E-np.diag([1,-1,-1,-1]))<1e-12
 original=M.read('data/w33_pass11389_parabolic_spatial_cover.json')['line']
 edges=np.array(r['cover_edges']);degree=np.bincount(edges.ravel());parameters=r['graph_parameters']
 assert len(np.unique(edges))==parameters['v']==80 and len(edges)==parameters['edges']==160 and np.all(degree==parameters['k']) and parameters['k']==4
 assert np.array_equal(r['integer_voltage'],np.array(original['integer_voltage'],int))
 for a,b in zip(r['acoustic_scans'][::2],r['acoustic_scans'][1::2]):assert abs(a['lowest_band']/b['lowest_band']-4)<1e-3

def test_native_loop_gradient_and_Hessian_by_high_precision_difference():
 import mpmath as mp
 mp.mp.dps=45;r=cert()['loops'];old=M.read('data/w33_pass11466_11470_stationary_ports_channels.json');x=np.array(old['vacuum']['coordinates']);rad=mp.mpf(str(r['radius']))
 def f(z):
  phi=mp.mpf('.81')*mp.sqrt(sum(a*a for a in z))/rad
  masses=[10-mp.sqrt(6)*phi,mp.mpf(10),10+mp.sqrt(6)*phi]
  return -3/(16*mp.pi**2)*sum(k*a**4*(mp.log(a*a/100)-mp.mpf('1.5')) for k,a in zip([268,340,268],masses))
 S=M.dec(old['soft']['soft_frame']).real;g=np.array(r['loop_gradient']);H=np.array(r['loop_Hessian']);xx=list(map(lambda a:mp.mpf(str(a)),x));h=mp.mpf('1e-5')
 for d in [x/np.linalg.norm(x),S[:,0],S[:,3]]:
  dd=list(map(lambda a:mp.mpf(str(a)),d));plus=[a+h*b for a,b in zip(xx,dd)];minus=[a-h*b for a,b in zip(xx,dd)]
  fd=float((f(plus)-f(minus))/(2*h));second=float((f(plus)-2*f(xx)+f(minus))/h**2)
  assert abs(fd-g@d)<2e-6 and abs(second-d@H@d)<2e-5
 assert r['raw_total_gradient_norm']>1000 and r['native_tree_gradient_norm']<1e-8

def test_Ward_and_mixed_curvature_from_actual_overlap_derivatives():
 import w33_pass11428_11432_native_alignment_measure_coarse_flags as Q
 r=cert()['ward'];B=np.array(r['gauge_boundary']);j=np.array(r['current']);a=np.array(r['anomaly_density']);dirs=np.array(r['directions']);A=np.array(r['background']).reshape(16,4)
 assert np.linalg.norm(B.T@j+a)<1e-11
 bound=max(M.plaquette_bound(A.reshape(-1)+step*d) for d in dirs for step in [0.,-1e-4,1e-4])
 assert bound==r['admissible_plaquette_bound'] and bound<1/30
 assert np.linalg.norm(B.T@np.array(r['Coulomb_horizontal_projector']))<1e-12
 actual=np.zeros((3,3));da=np.zeros(16)
 for q,m in zip([1,-4,2,-3,6],[6,3,3,2,1]):
  P,ds,_,_,_,_=Q.weak_overlap(q,1.,A,dirs.reshape(3,16,4),L=2)
  for i,k in product(range(3),repeat=2):actual[i,k]+=m*(1j*np.trace(P@(ds[i]@ds[k]-ds[k]@ds[i]))).real
  da+=m*q*np.trace(ds[0].reshape(16,4,16,4),axis1=1,axis2=3).diagonal().real
 assert np.linalg.norm(actual-np.array(r['actual_curvature']))<1e-12
 # Moment-map identity F(eta,B omega)=-d_eta a. omega delta atsite0.
 assert abs(actual[0,2]+da[0])<1e-12
 assert np.linalg.norm(np.array(r['finite_difference_curls'][-1])-actual)<1e-11
 assert np.linalg.norm(actual)>1e-10

def test_Pauli_transform_and_Choi_conventions_independently():
 rng=np.random.default_rng(3);A=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));A=A+A.conj().T
 c=M.pauli_coefficients(A);reconstructed=sum(c[i,j]*np.kron(M.PAULI[i],M.PAULI[j]) for i,j in product(range(4),repeat=2))
 assert np.linalg.norm(A-reconstructed)<1e-12
 # A genuinely nonunital two-qubit CPTP channel, evaluated directly by Kraus.
 p=.17;ks=[np.diag([1,np.sqrt(1-p)]),np.array([[0,np.sqrt(p)],[0,0]])];U=np.diag([1,1,1,-1]);K=[U@np.kron(a,b) for a,b in product(ks,repeat=2)]
 C=sum(np.outer(k.T.reshape(-1),k.T.reshape(-1).conj()) for k in K)/4;T=M.pair_transfer(C)
 basis=[np.kron(a,b) for a,b in product(M.PAULI,repeat=2)]
 direct=np.array([[sum(np.trace(a@k@b@k.conj().T).real for k in K)/4 for b in basis] for a in basis])
 assert np.linalg.norm(T-direct)<1e-12 and np.linalg.norm(M.transfer_choi(T)-C)<1e-12

def test_correlated_recovery_against_dense_independent_marginal_control():
 import w33_pass11471_11475_frames_currents_matching as Q
 r=cert()['recovery'];row=r['rows'][0];T=np.array(row['physical_transfer']);V,R=Q.recovery_rows();locals=[]
 for t in [T[np.ix_([0,4,8,12],[0,4,8,12])],T[:4,:4]]:
  superop=np.einsum('pq,pab,qji->abij',t,M.PAULI,M.PAULI)/2;logical=np.zeros((4,4))
  for b in range(4):
   state=Q.apply_product(superop,V@M.PAULI[b]@V.conj().T);out=np.einsum('sab,bc,sdc->ad',R,state,R.conj())
   for a in range(4):logical[a,b]=np.trace(M.PAULI[a]@out).real/2
  locals.append(logical)
 independent=np.kron(*locals);correlated=np.array(row['decoded_transfer'])
 assert abs(1-np.trace(independent)/16-row['independent_marginal_infidelity'])<1e-12
 assert abs(np.linalg.norm(independent-correlated)-row['independent_channel_difference'])<1e-12
 assert np.linalg.norm(independent-correlated)>.04
 for row in r['rows']:
  C=M.dec(row['decoded_Choi']);R=np.array(row['decoded_transfer']);omega=np.eye(4).reshape(-1)/2
  assert min(np.linalg.eigvalsh(C))>-1e-12 and np.linalg.norm(np.einsum('iaja->ij',C.reshape(4,4,4,4))-np.eye(4)/4)<1e-12
  assert abs(1-(omega@C@omega).real-row['decoded_entanglement_infidelity'])<1e-12
  assert np.linalg.norm(M.transfer_choi(R)-C)<1e-12
 assert r['shared_fresh_difference']>1e-3 and r['fresh_twirl_difference']>1e-3
