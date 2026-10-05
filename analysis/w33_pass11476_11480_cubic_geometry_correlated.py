"""Five constructed diagnostics following11471, with inherited results cited.

11471 owns the full signed-cubic carrier map;11389 owns the rank-three cover;
11466 owns the stationary point and relay channel;11475 owns the Steane decoder.
No regulator, renormalization convention or geometry choice predicts a TOE.
"""
import json, hashlib
from pathlib import Path
from itertools import product
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/w33_pass11476_11480_cubic_geometry_correlated.json'
read=lambda p:json.loads((ROOT/p).read_text())
dec=lambda x:np.array(x['real'])+1j*np.array(x['imag'])
enc=lambda x:dict(real=np.asarray(x).real.tolist(),imag=np.asarray(x).imag.tolist())
PAULI=np.array([np.eye(2),[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)

def pauli_coefficients(A):
 """All coefficients Tr(P A)/2^n, without allocating the full Pauli basis."""
 n=int(round(np.log2(len(A))));axes=sum(([i,n+i] for i in range(n)),[])
 tensor=A.reshape((2,)*(2*n)).transpose(axes).reshape((4,)*n)
 transform=PAULI.transpose(0,2,1).reshape(4,4)/2
 for _ in range(n):tensor=np.tensordot(tensor,transform.T,axes=(0,0))
 return tensor.real

def sparse_terms(A):
 c=pauli_coefficients(A);ix=np.array(np.where(abs(c)>1e-10)).T
 return ix,c[tuple(ix.T)]

def decoder_terms():
 import w33_pass11471_11475_frames_currents_matching as Q
 V,R=Q.recovery_rows()
 inputs=[sparse_terms(V@p@V.conj().T) for p in PAULI]
 outputs=[sparse_terms(np.einsum('sai,ab,sbj->ij',R.conj(),p,R)) for p in PAULI]
 return inputs,outputs

def pair_transfer(choi):
 C=choi.reshape(4,4,4,4)*4
 basis=np.array([np.kron(a,b) for a,b in product(PAULI,repeat=2)])
 # Normalized Choi convention: E(|i><j|)_{ab}=4 C_{ia,jb}.
 images=np.einsum('pij,iajb->pab',basis,C)
 return np.einsum('qab,pba->qp',basis,images).real/4

def contract_channel(T,inputs,outputs):
 """Exact finite Pauli expansion, retaining all inter-block channel entries.

 Decoder observable coefficients and encoder coefficients have different
 normalizations; Tr on two seven-qubit blocks contributes2^14/4=4096.
 """
 out=np.zeros((16,16))
 for a,b,c,d in product(range(4),repeat=4):
  ia,ca=outputs[a];ib,cb=outputs[b];ic,cc=inputs[c];id_,cd=inputs[d]
  # Batch output words; vectorize the two input stabilizer sums.
  input_pair=4*ic[:,None,:]+id_[None,:,:]
  val=0.
  for wa,xa in zip(ia,ca):
   opair=4*wa[None,:]+ib
   weights=np.ones((len(ib),len(ic),len(id_)))
   for j in range(7):weights*=T[opair[:,j,None,None],input_pair[None,:,:,j]]
   val+=xa*np.einsum('b,c,d,bcd->',cb,cc,cd,weights)
  out[4*a+b,4*c+d]=4096*val
 return out

def transfer_choi(T):
 basis=np.array([np.kron(a,b) for a,b in product(PAULI,repeat=2)])
 # E(P_j)=sum_i T_ij P_i; J=(1/16)sum_j P_j^T tensor E(P_j).
 return sum(T[i,j]*np.kron(basis[j].T,basis[i]) for i,j in product(range(16),repeat=2))/16

def correlated_recovery():
 import w33_pass11471_11475_frames_currents_matching as Q
 inputs,outputs=decoder_terms();U=np.zeros((4,4))
 for j in range(4):U[(j//2)*2+((j%2)^(j//2)),j]=1
 invert=np.kron(np.eye(4),U.T);rows=[];transfers={}
 for p,relay in product([.01,.03,.12],range(2)):
  C=Q.P.endpoint_channel(np.diag([1-relay,relay]),'damping',p)
  C=invert@C@invert.T;T=pair_transfer(C);R=contract_channel(T,inputs,outputs);J=transfer_choi(R)
  tp=np.linalg.norm(R[0]-np.eye(16)[0]);w=np.linalg.eigvalsh(J)
  # Deliberately wrong independent-marginal control, not an approximation label.
  T1=T[np.ix_([0,4,8,12],[0,4,8,12])];T2=T[:4,:4]
  independent=contract_channel(np.kron(T1,T2),inputs,outputs)
  assert tp<1e-10 and min(w)>-1e-10
  transfers[p,relay]=R
  rows.append(dict(p=p,relay=relay,physical_transfer=T.tolist(),decoded_transfer=R.tolist(),decoded_Choi=enc(J),bare_entanglement_infidelity=float(1-np.trace(T)/16),decoded_entanglement_infidelity=float(1-np.trace(R)/16),independent_marginal_infidelity=float(1-np.trace(independent)/16),independent_channel_difference=float(np.linalg.norm(R-independent)),TP_error=float(tp),minimum_Choi_eigenvalue=float(min(w))))
  print('correlated row',p,relay,rows[-1]['decoded_entanglement_infidelity'],flush=True)
 # Shared classical relay label across all seven positions versus fresh labels.
 p=.03;mix=(transfers[p,0]+transfers[p,1])/2
 C=sum(Q.P.endpoint_channel(np.diag([1-r,r]),'damping',p) for r in range(2))/2
 T=pair_transfer(invert@C@invert.T);fresh=contract_channel(T,inputs,outputs)
 twirl=contract_channel(np.diag(np.diag(T)),inputs,outputs)
 ideal=contract_channel(np.eye(16),inputs,outputs)
 assert np.linalg.norm(ideal-np.eye(16))<1e-10
 return dict(status='PASS',input_term_counts=[len(t[0]) for t in inputs],output_term_counts=[len(t[0]) for t in outputs],rows=rows,ideal_error=float(np.linalg.norm(ideal-np.eye(16))),shared_classical_relay_transfer=mix.tolist(),fresh_classical_relay_transfer=fresh.tolist(),shared_fresh_difference=float(np.linalg.norm(mix-fresh)),fresh_twirl_difference=float(np.linalg.norm(fresh-twirl)),scope='Full untwirled channel on two Steane blocks under seven independent routed endpoint channels, ideal encode/readout/recovery and ideal endpoint CNOT inversion. Endpoint correlations retained. Shared classical relay preparation additionally correlates positions through a stored latent bit. Not correlated quantum baths, noisy correction, physical leakage hardware or a fault-tolerance threshold.')

def determinant_tensor():
 bridge=read('data/w33_pass11471_11475_frames_currents_matching.json')['frame']['albert_bridge']
 p=np.array(bridge['twice_Jordan_product'],float)/2;t=np.einsum('ijj->i',p)/9
 g=np.einsum('k,ijk->ij',t,p);h=np.einsum('l,ial,jka->ijk',t,p,p)
 n=(np.einsum('i,j,k->ijk',t,t,t)-3*(np.einsum('i,jk->ijk',t,g)+np.einsum('j,ik->ijk',t,g)+np.einsum('k,ij->ijk',t,g))/3+2*h)/6
 return bridge,p,n

def vacuum_isotope(Z,K,bridge,p,n):
 """A numerical vacuum-adapted carrier; unitary isotope is classical theory."""
 t=np.einsum('ijj->i',p)/9;g=np.einsum('k,ijk->ij',t,p)
 def mul(a,b):return np.einsum('i,j,ijk->k',a,b,p)
 def triple(a,b,c):return mul(mul(a,b.conj()),c)+mul(mul(c,b.conj()),a)-mul(mul(a,c),b.conj())
 Q=np.linalg.qr(Z)[0];alltr=np.array([triple(a,b,c) for a,b,c in product(Z.T,repeat=3)]).T
 closure=np.linalg.norm(alltr-Q@(Q.conj().T@alltr))
 z=Z@np.array([1,2j,3+.7j]);D=np.column_stack([triple(z,z,a) for a in Z.T]);small=Z.conj().T@g@D
 ew,U=np.linalg.eigh((small+small.conj().T)/2);E=Z@U;c=E.sum(axis=1)
 norm=np.einsum('ijk,i,j,k',n,c,c,c);E[:,0]*=np.exp(-1j*np.angle(norm));c=E.sum(axis=1)
 basic=np.eye(27);pc=np.array([[triple(a,c,b) for b in basic] for a in basic])
 star=np.column_stack([triple(c,a,c) for a in basic])
 complete=np.linalg.norm(np.column_stack([triple(c,c,a) for a in basic])-basic)
 primitive=max(np.linalg.norm(triple(a,a,a)-a) for a in E.T)
 orthogonal=max(np.linalg.norm(triple(a,a,b)) for i,a in enumerate(E.T) for j,b in enumerate(E.T) if i!=j)
 Ti=dec(bridge['twice_inverse'])/2;T=dec(bridge['native_from_clock_basis']);R=np.array([Ti@(1j*a)@T for a in K])
 deriv=lambda A:np.einsum('ka,ija->ijk',A,pc)-np.einsum('ai,ajk->ijk',A,pc)-np.einsum('aj,iak->ijk',A,pc)
 deriv_error=max(np.linalg.norm(deriv(a)) for a in R)
 # The entire isotope determinant, not a basis-dimension or frame-only check.
 tc=np.einsum('ijj->i',pc)/9;gc=np.einsum('k,ijk->ij',tc,pc);hc=np.einsum('l,ial,jka->ijk',tc,pc,pc)
 nc=(np.einsum('i,j,k->ijk',tc,tc,tc)-(np.einsum('i,jk->ijk',tc,gc)+np.einsum('j,ik->ijk',tc,gc)+np.einsum('k,ij->ijk',tc,gc))+2*hc)/6
 norm_error=np.linalg.norm(nc-n);star_error=np.linalg.norm(star@star.conj()-basic)
 star_product_error=np.linalg.norm(np.einsum('ka,ija->ijk',star,pc.conj())-np.einsum('ai,bj,abk->ijk',star,star,pc))
 hermitian_trace=np.einsum('k,ajk,ai->ij',tc,pc,star);metric_error=np.linalg.norm(hermitian_trace-g)
 rng=np.random.default_rng(11476);jordan=[]
 for x in rng.normal(size=(5,27))+1j*rng.normal(size=(5,27)):
  xx=np.einsum('i,j,ijk->k',x,x,pc);Lx=np.einsum('i,ijk->kj',x,pc);Lxx=np.einsum('i,ijk->kj',xx,pc);jordan.append(float(np.linalg.norm(Lx@Lxx-Lxx@Lx)))
 assert closure<1e-12 and complete<1e-12 and primitive<1e-12 and orthogonal<1e-12 and deriv_error<1e-11 and norm_error<1e-12 and star_error<1e-12 and star_product_error<1e-12 and metric_error<1e-12 and max(jordan)<1e-10
 return dict(status='PASS',triple_closure_residual=float(closure),diagonalizer_spectrum=ew.tolist(),primitive_tripotents=enc(E),complete_tripotent=enc(c),native_tripotents=enc(T@E),isotope_product=enc(pc),isotope_conjugation=enc(star),native_conjugated_generators=enc(R),primitive_residual=float(primitive),orthogonality_residual=float(orthogonal),completeness_residual=float(complete),involution_residual=float(star_error),involution_product_residual=float(star_product_error),positive_trace_metric_residual=float(metric_error),trace_metric_eigenvalues=np.linalg.eigvalsh(hermitian_trace).tolist(),derivation_residual=float(deriv_error),full_cubic_tensor_residual=float(norm_error),Jordan_identity_controls=jordan,scope='Numerical actual-vacuum frame in the inherited Hermitian Jordan triple, then explicit unitary isotope x circle_c y={x,c,y}, star_c x={c,x,c}. All28 native numerical stabilizer generators derive this product and the entire cubic tensor is unchanged. Classical isotope theory retained; exact/interval vacuum certificate, full compact E6 conjugator and physical real-form/dynamics remain open.')

def polydisc_value(A,B,compressed):
 """36-real-field value of the supplied native potential on this frame."""
 def one(C):
  z=C.reshape(-1);mom=np.einsum('i,aij,j->a',z.conj(),compressed,z).real;u=27*np.linalg.det(C)**2
  return .5*mom@mom+abs(u)**4-abs(u)**2+u.real**2+u.real/2+abs(u)**6
 norm=np.linalg.norm(A)**2+np.linalg.norm(B)**2;cross=A.T@B.conj()/(norm/2+.01)
 U=2*np.eye(3)+(cross+cross.conj().T)/2;V=np.eye(3)+cross@cross.conj().T;W=U@V-V@U
 a=27*np.linalg.det(A)**2;chi=a.imag
 return float(one(A)+one(B)-1e4*chi/(1+chi*chi)*np.imag(np.trace(W@W@W))+1e-6*norm**2)

def polydisc_reduction(isotope,F,x):
 E=dec(isotope['native_tripotents']);embedding=np.kron(E,np.eye(3));H=np.array([embedding.conj().T@a@embedding for a in F.H])
 phi=(x[:81]+1j*x[81:162]).reshape(27,3);psi=(x[162:243]+1j*x[243:]);A=E.conj().T@phi;B=(embedding.conj().T@psi).reshape(3,3)
 rank=int(np.linalg.matrix_rank(np.c_[H.real.reshape(86,-1),H.imag.reshape(86,-1)],tol=1e-9))
 # SL3 family invariance plus homogeneity reduces a polynomial on invertible
 # matrices to its identity value; polynomial continuity covers detA=0.
 # Canonical frame arithmetic is exact at the integer point, not a random fit.
 from w33_pass11384_11388_native_dynamics import native_tensors
 tri=native_tensors()[1][0][0];canonical=np.zeros((27,3));canonical[tri[0],0]=1;canonical[tri[1],1]=1;canonical[tri[2],2]=tri[3]
 with F.native_context():
  u,v,*_=F.D.I.evaluate(canonical.reshape(-1));w,_=F.D.evaluate(canonical.reshape(-1))
 import sympy as sp
 with F.native_context():
  tensor=np.einsum('ijk,ia,jb,kc->abc',F.D.I.tensors()[1],canonical,canonical,canonical)
 keys,variables,S,T,_=F.D.cubic_circuit();sub=dict(zip(variables,[sp.Integer(int(tensor[k])) for k in keys]));rawT=T.subs(sub);rawS=S.subs(sub)
 exact18=(19683*rawT-sp.Rational(48384,5)*27**3+sp.Rational(13824,5)*27*729)/1492992
 assert u==27 and v==729 and exact18==0 and abs(w)<1e-10
 rng=np.random.default_rng(11476);rows=[]
 with F.native_context():
  for size in [0.,.003,.01]:
   C=A+size*(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)));D=B+size*(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
   full=F.pack((E@C).reshape(-1),(E@D).reshape(-1));value=F.fun(full)[0];small=polydisc_value(C,D,H)
   inv=F.D.I.evaluate((E@C).reshape(-1));i18=F.D.evaluate((E@C).reshape(-1))[0]
   errors=[abs(inv[0]-27*np.linalg.det(C)**2),abs(inv[1]-(27*np.linalg.det(C)**2)**2),abs(i18)]
   rows.append(dict(amplitude=size,A=enc(C),B=enc(D),native_value=value,reduced_value=small,error=abs(value-small),invariant_errors=list(map(float,errors))))
 assert np.linalg.norm(E.conj().T@E-np.eye(3))<1e-12 and max(row['error'] for row in rows)<1e-9
 return dict(status='PASS',native_embedding=enc(embedding),compressed_moment_generators=enc(H),compressed_real_span_rank=rank,A_at_vacuum=enc(A),B_at_vacuum=enc(B),canonical_frame=list(map(int,tri)),exact_canonical_constants=[27,729,0],exact_Aronhold_S=int(rawS),exact_Hessian_T=int(rawT),canonical_float_I18_residual=float(abs(w)),invariant_laws=['I6(E A)=27 det(A)^2','I12(E A)=729 det(A)^4','I18(E A)=0'],rows=rows,scope='Classical SL3 invariance/homogeneity proves these restrictions on the canonical integer frame from exact constants27,729,0; actual numerical tripotent-frame controls replay them. Supplied324-field potential reduces to two3x3 complex matrices and compressed moments on the fixed stratum. No full transverse/global stability, exact numerical-frame algebraicity or physical vacuum prediction.')

def vacuum_invariants():
 from scipy.linalg import null_space
 import w33_pass11438_finite_native_model as F
 src=read('data/w33_pass11466_11470_stationary_ports_channels.json');x=np.array(src['vacuum']['coordinates']);S=dec(src['soft']['soft_frame']).real
 bridge,p,n=determinant_tensor();Ti=dec(bridge['twice_inverse'])/2
 K=dec(src['vacuum']['generators']);W=null_space(K.reshape(-1,27),rcond=1e-9);Z=Ti@W;Q=np.linalg.qr(Z)[0]
 unit=np.zeros(27);unit[bridge['clock_idempotent_indices']]=1
 products=np.einsum('ia,jb,ijk->kab',Z,Z,p);res=products-np.einsum('ki,ij,jab->kab',Q,Q.conj().T,products)
 # Gauge-invariant mixed native sextics, with analytic gradients and two FD steps.
 phi=x[:81]+1j*x[81:162];psi=x[162:243]+1j*x[243:]
 dphi=S[:81]+1j*S[81:162];dpsi=S[162:243]+1j*S[243:]
 probes=[1,-1,1j,2,-2,2j];analytic=[];fds=[]
 with F.native_context():
  for a in probes:
   grad=F.D.I.evaluate(phi+a*psi)[2];analytic.append(grad@(dphi+a*dpsi))
  for step in [1e-4,5e-5]:
   row=[]
   for a in probes:
    row.append([(F.D.I.evaluate(phi+a*psi+step*(dphi[:,j]+a*dpsi[:,j]))[0]-F.D.I.evaluate(phi+a*psi-step*(dphi[:,j]+a*dpsi[:,j]))[0])/(2*step) for j in range(4)])
   fds.append(np.array(row))
 D=np.array(analytic);stack=np.r_[D.real,D.imag];sv=np.linalg.svd(stack,compute_uv=False)
 # A vanishing linear probe is inconclusive on the enhanced-symmetry stratum.
 nonlinear=[]
 with F.native_context():
  for a in probes:
   f0=F.D.I.evaluate(phi+a*psi)[0];rows=[]
   for amp in [.03,.06]:
    rows.append([complex(F.D.I.evaluate(phi+a*psi+amp*(dphi[:,j]+a*dpsi[:,j]))[0]-f0) for j in range(4)])
   nonlinear.append(enc(np.array(rows)))
 assert W.shape==(27,3) and np.linalg.norm(D-fds[-1])<2e-5
 isotope=vacuum_isotope(Z,K,bridge,p,n);polydisc=polydisc_reduction(isotope,F,x)
 return dict(status='PASS',isotope=isotope,polydisc=polydisc,native_fixed_basis=enc(W),clock_fixed_basis=enc(Z),clock_unit_distance=float(np.linalg.norm(unit-Q@(Q.conj().T@unit))),Jordan_closure_residual=float(np.linalg.norm(res)),clock_reality_fraction=float(np.linalg.norm(Z.imag)/np.linalg.norm(Z)),mixed_sextic_probes=enc(np.array(probes)),analytic_Jacobian=enc(D),finite_difference_Jacobians=[enc(d) for d in fds],finite_difference_errors=[float(np.linalg.norm(D-d)) for d in fds],singular_values=sv.tolist(),nonlinear_changes=nonlinear,scope='Actual numerical native fixed space transported by11471 exact cubic map. Ordinary unit membership/Jordan closure fail; the Hermitian triple instead supplies a verified numerical vacuum-adapted isotope. Analytic and independent finite-difference mixed-sextic probes distinguish numerical sensitivity from exact moduli; no exact stationary point, moduli theorem or mass prediction.')

def determinant_cover():
 bridge,p,n=determinant_tensor();ids=bridge['clock_idempotent_indices'];e=np.zeros(27);e[ids]=1
 rows=[]
 for abc in [[1,1,1],[1,2,3],[1,-1,1],[1,0,1]]:
  y=np.zeros(27);y[ids]=abc;N=np.einsum('ijk,i,j,k',n,y,y,y);grad=3*np.einsum('ijk,j,k->i',n,y,y);H=6*np.einsum('ijk,k->ij',n,y);ev=np.linalg.eigvalsh(H)
  logH=H/N-np.outer(grad,grad)/N**2 if N!=0 else np.zeros_like(H)
  rows.append(dict(frame_values=abc,determinant=float(N),Hessian=H.tolist(),inertia=[int(sum(ev>1e-9)),int(sum(ev < -1e-9)),int(sum(abs(ev)<=1e-9))],log_Hessian_inertia=[int(sum(np.linalg.eigvalsh(logH)>1e-9)),int(sum(np.linalg.eigvalsh(logH)<-1e-9))] if N!=0 else None))
 assert rows[0]['inertia']==[1,26,0] and rows[1]['inertia']==[1,26,0] and rows[0]['log_Hessian_inertia']==[0,27]
 H=np.array(rows[0]['Hessian']);w,V=np.linalg.eigh(H);time=V[:,-1]/np.sqrt(w[-1]);space=V[:,:3]/np.sqrt(-w[:3]);embedding=np.column_stack([time,space]);metric=embedding.T@H@embedding
 from fractions import Fraction
 cover=read('data/w33_pass11389_parabolic_spatial_cover.json')['line'];edges=np.array(cover['edges']);voltage=np.array(cover['integer_voltage'],int);harmonic=np.array([[float(Fraction(a)) for a in row] for row in cover['fcc_harmonic_displacements']])
 vertices=np.unique(edges);degree=np.bincount(edges.ravel());assert np.array_equal(vertices,np.arange(80)) and np.all(degree==4)
 graph_parameters=dict(v=len(vertices),k=int(degree[0]),edges=len(edges),meaning='Actual four-regular80-vertex Levi base of11389 voltage cover; not SRG parameters or a physical interpretation of code dimensions.')
 def bands(k):
  lap=np.diag(degree).astype(complex)
  for (a,b),v in zip(edges,voltage):lap[a,b]-=np.exp(1j*k@v);lap[b,a]-=np.exp(-1j*k@v)
  return np.linalg.eigvalsh(lap)
 scans=[]
 for direction in np.eye(3):
  for step in [.01,.005]:scans.append(dict(momentum=(step*direction).tolist(),lowest_band=float(bands(step*direction)[0])))
 metric_error=float(np.linalg.norm(metric-np.diag([1,-1,-1,-1])))
 assert metric_error<1e-12 and all(r['lowest_band']>0 for r in scans)
 return dict(status='PASS',graph_parameters=graph_parameters,frame_controls=rows,selected_4D_embedding=embedding.tolist(),selected_metric=metric.tolist(),metric_error=metric_error,cover_edges=edges.tolist(),integer_voltage=voltage.tolist(),harmonic_FCC_displacements=harmonic.tolist(),acoustic_scans=scans,scope='Transported native cubic gives a Lorentz-signature Hessian on positive diagonal Albert elements, whereas log determinant gives a negative-definite Hessian. A chosen normalized four-plane is paired with11389 actual three-period cover. Choice of Hessian, positive real form, time and spatial embedding are inputs; signature does not derive spacetime, Einstein equations or backreaction. Known Jordan/cone theory and11389 cover retain ownership.')

def native_loop_feedback():
 import w33_pass11438_finite_native_model as F
 import mpmath as mp
 src=read('data/w33_pass11466_11470_stationary_ports_channels.json');x=np.array(src['vacuum']['coordinates']);S=dec(src['soft']['soft_frame']).real;r=np.linalg.norm(x);p0=.81;mu=10.;counts=[268,340,268]
 def V(phi):
  masses=[10-mp.sqrt(6)*phi,mp.mpf(10),10+mp.sqrt(6)*phi]
  return -3/(16*mp.pi**2)*sum(c*m**4*(mp.log(m*m/mu**2)-mp.mpf('1.5')) for c,m in zip(counts,masses))
 with mp.workdps(40):v=[float(mp.diff(V,mp.mpf('.81'),j)) for j in range(5)]
 # Explicit gauge-invariant radius port: supplied normalization, not derived Yukawa map.
 dphi=p0*x/r**2;Hphi=p0*(np.eye(len(x))/r**2-np.outer(x,x)/r**4)
 loopgrad=v[1]*dphi;loopH=v[2]*np.outer(dphi,dphi)+v[1]*Hphi
 with F.native_context():value,g=F.fun(x)
 rows=[]
 for amplitude in [-.05,-.01,.01,.05]:
  phi=p0*(1+amplitude)
  with mp.workdps(40):raw=float(V(mp.mpf(str(phi)))-V(mp.mpf('.81')))
  sub=raw-v[1]*(phi-p0)-v[2]*(phi-p0)**2/2
  with F.native_context():tree=F.fun((1+amplitude)*x)[0]-value
  rows.append(dict(radial_amplitude=amplitude,phi=phi,tree_change=tree,loop_change=raw,tadpole_mass_subtracted_loop_change=sub))
 assert np.linalg.norm(loopgrad)>100*np.linalg.norm(g)
 return dict(status='PASS',radius=r,phi_at_vacuum=p0,phi_rule='phi(x)=.81*||x||/||x_vac||',native_tree_gradient_norm=float(np.linalg.norm(g)),raw_matched_loop_gradient_norm=float(np.linalg.norm(loopgrad)),raw_total_gradient_norm=float(np.linalg.norm(g+loopgrad)),loop_gradient=loopgrad.tolist(),loop_Hessian=loopH.tolist(),soft_loop_Hessian=(S.T@loopH@S).tolist(),CW_derivatives=v,radial_rows=rows,renormalization_control='Subtract V(phi0), Vprime(phi0)*(phi-phi0), Vsecond(phi0)*(phi-phi0)^2/2: native stationary point and Hessian are restored by construction, but higher loop terms remain.',scope='Actual876 omitted species coupled to the actual324-field native potential through an explicitly supplied gauge-invariant radius port. Large tadpole destroys the prior tree stationary point unless matching conditions are changed. Finite tadpole/mass subtraction is a distinct convention, not a prediction. No observed masses, cosmological constant or unique quantum vacuum.')

def gauge_boundary(L):
 sites=list(product(range(L),repeat=4));ix={s:i for i,s in enumerate(sites)};B=np.zeros((4*len(sites),len(sites)))
 for i,s in enumerate(sites):
  for mu in range(4):
   t=list(s);t[mu]=(t[mu]+1)%L;j=ix[tuple(t)];B[4*i+mu,i]=1;B[4*i+mu,j]-=1
 return B

def plaquette_bound(field,L=2):
 sites=list(product(range(L),repeat=4));ix={s:i for i,s in enumerate(sites)};A=np.asarray(field).reshape(len(sites),4);bound=0.
 for i,s in enumerate(sites):
  for mu,nu in product(range(4),repeat=2):
   if mu>=nu:continue
   sm=list(s);sn=list(s);sm[mu]=(sm[mu]+1)%L;sn[nu]=(sn[nu]+1)%L
   angle=A[i,mu]+A[ix[tuple(sm)],nu]-A[ix[tuple(sn)],mu]-A[i,nu];bound=max(bound,6*abs(angle))
 return float(bound)

def ward_connection():
 import w33_pass11428_11432_native_alignment_measure_coarse_flags as M
 from numpy.polynomial.legendre import leggauss
 L=2;B=gauge_boundary(L);inv=np.linalg.pinv(B.T@B);horizontal=np.eye(64)-B@inv@B.T
 rng=np.random.default_rng(11479);A=.0005*rng.normal(size=64);charges=[1,-4,2,-3,6];mult=[6,3,3,2,1]
 basis=np.eye(64).reshape(64,16,4);nodes,weights=leggauss(8);nodes=(nodes+1)/2;weights=weights/2
 def local_data(field,tangents):
  anomaly=np.zeros(16);curv=np.zeros((len(tangents),len(tangents)));gap=10.
  for q,m in zip(charges,mult):
   P,ds,_,g,_,_=M.weak_overlap(q,1.,field.reshape(16,4),np.array(tangents).reshape(-1,16,4),L=2);gap=min(gap,g)
   anomaly+=m*q*np.trace(P.reshape(16,4,16,4),axis1=1,axis2=3).diagonal().real
   for i,j in product(range(len(ds)),repeat=2):curv[i,j]+=m*(1j*np.trace(P@(ds[i]@ds[j]-ds[j]@ds[i]))).real
  return anomaly,curv,gap
 def connection(field):
  h=horizontal@field;radial=np.zeros(64);integrated_anomaly=np.zeros(16)
  for t,w in zip(nodes,weights):
   # Only radial versus each link derivative is needed, not65^2 curvatures.
   for q,m in zip(charges,mult):
    P,ds,_,_,_,_=M.weak_overlap(q,1.,(t*h).reshape(16,4),np.r_[h.reshape(1,16,4),basis],L=2)
    radial+=w*t*m*np.array([(1j*np.trace(P@(ds[0]@d-d@ds[0]))).real for d in ds[1:]])
    integrated_anomaly+=w*m*q*np.trace(ds[0].reshape(16,4,16,4),axis1=1,axis2=3).diagonal().real
  anomaly,_,gap=local_data(h,[])
  # a(0)=0; integrating its analytic derivative avoids subtracting order-one
  # diagonal traces to obtain order1e-12 density, especially in FD curls.
  return horizontal@radial-B@inv@integrated_anomaly,integrated_anomaly,gap,float(np.linalg.norm(anomaly-integrated_anomaly))
 j,a,gap,density_error=connection(A);ward=B.T@j+a
 eta=horizontal[:,0];zeta=horizontal[:,7];omega=np.zeros(16);omega[0]=1;vertical=B@omega
 dirs=[eta,zeta,vertical];_,F,_=local_data(A,dirs);curls=[]
 for step in [1e-4,5e-5]:
  curl=np.zeros((3,3));derivatives=np.column_stack([(connection(A+step*d)[0]-connection(A-step*d)[0])/(2*step) for d in dirs])
  for i,k in product(range(3),repeat=2):curl[i,k]=dirs[k]@derivatives[:,i]-dirs[i]@derivatives[:,k]
  curls.append(curl)
 gauge_error=np.linalg.norm(connection(A+.0003*vertical)[0]-j)
 bound=max(plaquette_bound(A+step*d) for d in dirs for step in [0.,-1e-4,1e-4]);assert bound<1/30
 kernels=[]
 for size in [2,3,4]:
  bb=gauge_boundary(size);ii=np.linalg.pinv(bb.T@bb);kernel=bb@ii;col=kernel[:,0];distance=[]
  sites=list(product(range(size),repeat=4))
  for d in range(2*size+1):
   ids=[4*i+mu for i,s in enumerate(sites) for mu in range(4) if sum(min(a,size-a) for a in s)==d]
   if ids:distance.append(dict(distance=d,maximum=float(np.max(abs(col[ids])))))
  kernels.append(dict(L=size,distance_profile=distance))
 assert np.linalg.norm(ward)<1e-11 and np.linalg.norm(curl-F)<1e-11 and gauge_error<1e-10 and density_error<1e-11
 return dict(status='PASS',L=L,background=A.tolist(),gauge_boundary=B.tolist(),Coulomb_horizontal_projector=horizontal.tolist(),current=j.tolist(),anomaly_density=a.tolist(),admissible_plaquette_bound=bound,direct_density_replay_error=density_error,Ward_residual=float(np.linalg.norm(ward)),directions=np.array(dirs).tolist(),actual_curvature=F.tolist(),finite_difference_curls=[c.tolist() for c in curls],finite_difference_steps=[1e-4,5e-5],curl_error=float(np.linalg.norm(curl-F)),gauge_invariance_error=float(gauge_error),minimum_Wilson_gap=gap,Coulomb_kernel_profiles=kernels,scope='Finite contractible trivial-flux Coulomb-slice connection has exact Ward divergence, gauge covariance and curvature integrability, with independent derivative checks. Hodge projection uses massless inverse graph Laplacian: spatial locality is not established and cannot be inferred from this finite construction. Luscher local reconstruction, flux-sector gluing, non-Abelian completion and continuum remain open; anomalous multiplets can also admit nonlocal patch constructions.')

def run():
 results={}
 for name,fn in [('vacuum',vacuum_invariants),('geometry',determinant_cover),('loops',native_loop_feedback),('ward',ward_connection),('recovery',correlated_recovery)]:
  results[name]=fn();print(name,'PASS',flush=True)
  Path('/tmp/w33_11476_partial.json').write_text(json.dumps(results))
 files=['data/w33_pass11471_11475_frames_currents_matching.json','data/w33_pass11466_11470_stationary_ports_channels.json','data/w33_pass11389_parabolic_spatial_cover.json']
 results.update(status='PASS',passes=list(range(11476,11481)),reservation='394364e1c',source_sha256={p:hashlib.sha256(json.dumps(read(p),sort_keys=True,separators=(',',':')).encode()).hexdigest() for p in files})
 OUT.write_text(json.dumps(results,indent=2)+'\n');return results
if __name__=='__main__':run()
