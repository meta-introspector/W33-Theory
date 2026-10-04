"""Five constructed follow-ups to11466; finite controls are not a TOE closure.

Prior owners:11384 native signed E6;10956 exact clock Albert Spin8;
11389 periodic spatial cover and4045 supplied wave action;
11428 overlap derivative;11466 mediator reduction and relay noise;
11443 erasure decoder. Imported methods are cited in the companion report.
"""
import json,hashlib
from pathlib import Path
from itertools import product
import numpy as np
from scipy.linalg import null_space
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/w33_pass11471_11475_frames_currents_matching.json'
import w33_pass11466_11470_stationary_ports_channels as P
enc=P.enc

def source():return json.loads(P.OUT.read_text())

def albert_bridge(K):
 import w33_pass10950_clock_albert_lorentz_spinor as C
 J=C.build_clock_albert();T=np.zeros((27,27),complex)
 # A searched signed witness, verified below rather than inferred from labels.
 flips={1,4,6,7,26}
 for j,(kind,a,b) in enumerate(J['labels']):
  if kind=='fixed':T[a,j]=1
  else:
   signs=[sg for triple,sg in J['dsign'].items() if a in triple and b in triple and any(g in triple for g in J['frame'])];assert len(signs)==1;sg=-signs[0]
   if kind=='pair_re':T[a,j]=1;T[b,j]=sg
   else:T[a,j]=1j;T[b,j]=-1j*sg
 for i in flips:T[i]*=-1
 Ti2=2*np.linalg.inv(T);Pi2=np.array(J['prod'],float)*2
 # All values are Gaussian integers after these denominator clearings.
 assert np.array_equal(Ti2.real,np.rint(Ti2.real)) and np.array_equal(Ti2.imag,np.rint(Ti2.imag))
 Pi2=np.rint(Pi2).astype(np.int64);assert np.array_equal(Pi2/2,np.array(J['prod'],float))
 assert np.array_equal(Ti2@T,2*np.eye(27))
 R2=np.array([Ti2@a@T for a in K])
 assert np.array_equal(R2.real,np.rint(R2.real)) and np.array_equal(R2.imag,np.rint(R2.imag))
 def deriv(A):return np.einsum('ka,ija->ijk',A,Pi2)-np.einsum('ai,ajk->ijk',A,Pi2)-np.einsum('aj,iak->ijk',A,Pi2)
 for A in R2:
  assert not np.any(deriv(A.real.astype(np.int64))) and not np.any(deriv(A.imag.astype(np.int64)))
 # Compare the entire signed native determinant, not just D4 dimensions.
 # For Jordan trace t, 12N=2t^3-3t*Tr(x^2)*2+Tr(x^3)*4.
 from w33_pass11384_11388_native_dynamics import native_tensors
 signed=native_tensors()[1][1]
 raw=np.einsum('ijj->i',Pi2);assert not np.any(raw%18);trace=raw//18
 gram2=np.einsum('k,ijk->ij',trace,Pi2)
 cube4=np.einsum('l,ial,jka->ijk',trace,Pi2,Pi2)
 assert np.array_equal(cube4,cube4.transpose(1,0,2)) and np.array_equal(cube4,cube4.transpose(2,1,0))
 norm12=2*np.einsum('i,j,k->ijk',trace,trace,trace)-(np.einsum('i,jk->ijk',trace,gram2)+np.einsum('j,ik->ijk',trace,gram2)+np.einsum('k,ij->ijk',trace,gram2))+cube4
 # signed is six times the native cubic coefficient tensor; T is Gaussian
 # integral and intermediates are far below the exact integer range of doubles.
 mapped=2*np.einsum('abc,ai,bj,ck->ijk',signed,T,T,T,optimize=True)
 assert not np.any(mapped.imag) and np.array_equal(mapped.real,norm12)
 compact=np.r_[K-K.transpose(0,2,1),1j*(K+K.transpose(0,2,1))];compactR2=np.array([Ti2@a@T for a in compact])
 assert not np.any(compactR2.imag) and np.linalg.matrix_rank(compactR2.real.reshape(56,-1))==28
 # Both signed coordinate frame and actual Albert idempotents are fixed.
 idempotents=[J['idx'][g] for g in J['G']]
 assert not np.any(R2[:,:,idempotents])
 return dict(status='PASS',clock_frame=list(J['frame']),clock_basis_labels=[list(x) for x in J['labels']],sign_flips=sorted(flips),native_from_clock_basis=enc(T),twice_inverse=enc(Ti2),twice_Jordan_product=Pi2.tolist(),twice_conjugated_generators=enc(R2),compact_real_dimension=28,clock_idempotent_indices=idempotents,cubic_equality='12*N_native(T*y)=12*N_Albert(y)',cubic_coefficient_residual=0,scope='Explicit Gaussian-rational carrier map from the canonical signed-native compact frame algebra to the actual clock Albert frame derivations. Denominator-cleared product derivation identities are checked with integer arithmetic; compact images are real and span28. The full native signed cubic also equals the clock Albert determinant coefficientwise after denominator clearing. Physical real-form selection and a conjugator from the numerical stationary vacuum remain open.')

def frame():
 import sympy as sp
 from w33_pass11384_11388_native_dynamics import native_tensors
 B,tensors=native_tensors();tri=tensors[0];indices=tri[0][:3]
 # A complex E6 stabilizer of three actual signed-cubic coordinate vectors.
 J=np.eye(27,dtype=int)[:,indices];C=np.array([b@J for b in B]).reshape(78,-1).T
 ns=sp.polys.matrices.DomainMatrix.from_Matrix(sp.Matrix(C.tolist())).nullspace().to_Matrix()
 Z=np.array(ns,int);K=np.einsum('ab,bij->aij',Z,B)
 assert len(K)==28 and not np.any(K@J)
 # The subspace is transpose-closed, hence defines a compact real form.
 flat=K.reshape(28,-1).T;Q,_=np.linalg.qr(flat);transpose_error=np.linalg.norm(K.transpose(0,2,1).reshape(28,-1).T-Q@(Q.T@K.transpose(0,2,1).reshape(28,-1).T))
 cas=sum(k.T@k for k in K);ev=np.linalg.eigvalsh(cas)
 assert transpose_error<1e-10 and sum(ev<1e-9)==3
 # Exact invariant support blocks and the native cubic triality contraction.
 adj=np.any(K!=0,axis=0);unseen=set(range(27));blocks=[]
 while unseen:
  seed=min(unseen);block=[seed];unseen.remove(seed)
  for i in block:
   for j in np.flatnonzero(adj[i]):
    if int(j) in unseen:unseen.remove(int(j));block.append(int(j))
  blocks.append(sorted(block))
 assert sorted(map(len,blocks))==[1,1,1,8,8,8]
 from collections import Counter
 labels={v:i for i,b in enumerate(blocks) for v in b}
 patterns=Counter(tuple(sorted(labels[x] for x in t[:3])) for t in tri)
 assert sorted(patterns.values())==[1,4,4,4,32]
 # Probe the current four soft modes with gauge-invariant observables.
 r=source();x=np.array(r['vacuum']['coordinates']);soft=P.dec(r['soft']['soft_frame']).real
 def obs(y):
  a=(y[:81]+1j*y[81:162]).reshape(27,3);b=(y[162:243]+1j*y[243:]).reshape(27,3)
  C=a.T@b.conj();aa=a.conj().T@a;bb=b.conj().T@b
  return np.r_[np.linalg.eigvalsh(aa),np.linalg.eigvalsh(bb),np.linalg.svd(C,compute_uv=False),np.trace(C).real,np.trace(C).imag]
 h=1e-5;D=np.column_stack([(obs(x+h*d)-obs(x-h*d))/(2*h) for d in soft.T]);s=np.linalg.svd(D,compute_uv=False)
 return dict(status='PASS',frame_indices=list(map(int,indices)),frame_cubic_sign=int(tri[0][3]),albert_bridge=albert_bridge(K),integer_generator_coefficients=Z.tolist(),integer_generators=K.tolist(),constraint_rank=int(np.linalg.matrix_rank(C.astype(float))),dimension=28,transpose_closure_error=float(transpose_error),fixed_dimension=3,invariant_coordinate_blocks=blocks,cubic_block_patterns=[dict(blocks=list(k),count=v) for k,v in sorted(patterns.items())],soft_observable_Jacobian=D.tolist(),soft_observable_singular_values=s.tolist(),soft_observable_rank=int(sum(s>1e-6)),scope='Exact native complex E6 kernel fixing a nonzero signed-cubic coordinate frame, with transpose closure checked numerically and three fixed vectors. Closure of the kernel follows from the parent Lie algebra and fixing equations; the separately stored Gaussian-rational map identifies this canonical compact frame algebra with clock Albert frame derivations; no numerical-vacuum conjugator is asserted. Soft observable Jacobian tests whether candidate modes change gauge invariants; it is not an exact moduli theorem.')

def current():
 import w33_pass11428_11432_native_alignment_measure_coarse_flags as M
 from numpy.polynomial.legendre import leggauss
 rng=np.random.default_rng(11472);L=2;background=.0005*rng.normal(size=(16,4));f=np.zeros((2,16,4));f[0,0,0]=1;f[1,1,1]=1
 charges=[1,-4,2,-3,6];mult=[6,3,3,2,1]
 def curvature(x,y):
  vals=[];gaps=[]
  for q in charges:
   _,_,c,g,_,_=M.weak_overlap(q,1.,background+x*f[0]+y*f[1],f,L=L);vals.append(c);gaps.append(g)
  return float(np.dot(mult,vals)),min(gaps)
 nodes,weights=leggauss(8);nodes=(nodes+1)/2;weights=weights/2
 def j(x,y):
  integral=sum(w*t*curvature(t*x,t*y)[0] for t,w in zip(nodes,weights));return np.array([-y,x])*integral
 x,y=.001,-.001;h=2e-5;c,g=curvature(x,y);J=j(x,y)
 curl=((j(x+h,y)[1]-j(x-h,y)[1])-(j(x,y+h)[0]-j(x,y-h)[0]))/(2*h)
 # The complete radial patch and derivative endpoints stay in the same
 # admissible trivial-flux sector: every plaquette angle is a linear form.
 sites=list(product(range(L),repeat=4));ix={s:i for i,s in enumerate(sites)};bound=0.
 for ax,ay in [(0.,0.),(x+h,y),(x-h,y),(x,y+h),(x,y-h)]:
  field=background+ax*f[0]+ay*f[1]
  for i,site in enumerate(sites):
   for mu in range(4):
    for nu in range(mu+1,4):
     sm=list(site);sn=list(site);sm[mu]=(sm[mu]+1)%L;sn[nu]=(sn[nu]+1)%L
     angle=field[i,mu]+field[ix[tuple(sm)],nu]-field[ix[tuple(sn)],mu]-field[i,nu]
     bound=max(bound,max(abs(q*angle) for q in charges))
 assert abs(curl-c)<2e-7 and g>.5 and bound<1/30
 return dict(status='PASS',L=L,charges=charges,multiplicities=mult,background=background.tolist(),tangents=f.tolist(),point=[x,y],current=J.tolist(),curvature=c,finite_difference_curl=float(curl),integrability_error=float(abs(curl-c)),minimum_Wilson_gap=g,admissible_plaquette_angle_bound=bound,quadrature_order=8,formula='j_x=-y integral_0^1 t F(tx,ty)dt; j_y=x integral_0^1 t F(tx,ty)dt',scope='Constructed radial-homotopy measure connection on a two-link contractible patch of the anomaly-free hypercharge multiplet. Its curl equals the actual overlap determinant curvature. This current is local in parameter space, not established spatially local or gauge covariant; neither anomalous nor anomaly-free bundle curvature alone supplies the Luscher reconstruction conditions. Nonzero flux sectors and non-Abelian completion remain open.')

def wave():
 M=P.native();A=np.array(M.N.P.N.M.P.load()['A'],int);L=4*np.eye(80)-A
 eig=np.linalg.eigvalsh(L);dt=.5;N=100;rng=np.random.default_rng(11473);q=rng.normal(size=80);q/=np.linalg.norm(q);prev=q.copy();history=[q.copy()];energies=[]
 for _ in range(N):
  nxt=2*q-prev-dt*dt*L@q
  energies.append(float(np.dot((nxt-q)/dt,(nxt-q)/dt)/2+q@L@nxt/2));prev,q=q,nxt;history.append(q.copy())
 # Exact discrete action stationarity for every interior site/time, not a
 # homogeneous mode ansatz. Speed of support is <=one Levi edge per tick.
 H=np.array(history);force=(H[2:]-2*H[1:-1]+H[:-2])/dt**2+H[1:-1]@L
 probe=np.zeros(80);probe[0]=1;last=probe.copy();supports=[]
 reach=np.eye(80,dtype=int)[0].astype(bool)
 for t in range(1,5):
  nxt=2*probe-last-dt*dt*L@probe;reach|=(A@reach.astype(int))>0
  assert not np.any(abs(nxt[~reach])>1e-12)
  supports.append(int(sum(abs(nxt)>1e-12)));last,probe=probe,nxt
 unstable_dt=.8;growth=np.roots([1,unstable_dt**2*8-2,1]);assert max(abs(growth))>1
 assert max(abs(force.ravel()))<1e-12 and np.ptp(energies)<1e-12
 return dict(status='PASS',adjacency=A.tolist(),dt=dt,steps=N,initial_data=H[0].tolist(),final_data=H[-1].tolist(),Laplacian_eigenvalues=eig.tolist(),dispersion='4 sin^2(omega dt/2)=dt^2 lambda',CFL_dt_limit=float(2/np.sqrt(8)),maximum_Euler_Lagrange_residual=float(max(abs(force.ravel()))),discrete_energy_drift=float(np.ptp(energies)),support_sizes=supports,unstable_dt=unstable_dt,unstable_growth_factor=float(max(abs(growth))),scope='Explicit W33 Levi vertices map to spatial field sites, Levi incidence edges to spatial action terms, and a discrete time product supplies inhomogeneous local hyperbolic dynamics. Kinetic sign, time coordinate and coefficient are supplied. Finite graph dispersion and a causal support cone are verified; no continuum Lorentz symmetry, emergent3-space, Regge geometry map or gravitational backreaction is derived. This is an independent diagnostic alternative to homogeneous frusta.')

def matching():
 import mpmath as mp
 from numpy.polynomial.legendre import leggauss
 phi=.81;mu=10.;Nc=3;mass=np.array([10-np.sqrt(6)*phi,10,10+np.sqrt(6)*phi]);n=np.array([268,340,268]);T=.5
 def V(z):
  m=[10-mp.sqrt(6)*z,mp.mpf(10),10+mp.sqrt(6)*z]
  return -Nc/(16*mp.pi**2)*sum(int(nn)*a**4*(mp.log(a*a/mu**2)-mp.mpf('1.5')) for nn,a in zip(n,m))
 with mp.workdps(40):
  p=mp.mpf(str(phi));derivatives=[float(mp.diff(V,p,k)) for k in range(5)];coeff=[float(mp.diff(V,mp.mpf(0),k)/mp.factorial(k)) for k in [0,2,4,6]]
  h=mp.mpf('1e-5');fd=float((V(p+h)-V(p-h))/(2*h))
 nodes,w=leggauss(32);x=(nodes+1)/2;w=w/2;rows=[]
 for Q in [.01,.1,1.]:
  exact=T/(2*np.pi**2)*sum(nn*np.dot(w,x*(1-x)*np.log1p(Q*Q*x*(1-x)/(m*m))) for nn,m in zip(n,mass))
  expansion=T/(2*np.pi**2)*(Q*Q/30*np.sum(n/mass**2)-Q**4/280*np.sum(n/mass**4))
  rows.append(dict(Euclidean_momentum=Q,subtracted_polarization=float(exact),matched_Q2_Q4=float(expansion),error=float(abs(exact-expansion))))
 assert abs(fd-derivatives[1])<1e-5 and rows[0]['error']<1e-12 and np.all(mass>0)
 return dict(status='PASS',discarded_Dirac_species=int(sum(n)),color_multiplicity=Nc,phi=phi,renormalization_scale=mu,masses=mass.tolist(),multiplicities=n.tolist(),CW_value_and_first_four_derivatives=derivatives,CW_even_Taylor_coefficients=coeff,CW_derivative_replay_error=abs(fd-derivatives[1]),Dynkin_index=T,polarization=rows,scope='One-loop MSbar Coleman-Weinberg potential and subtracted Euclidean gauge polarization of the876 discarded triplet Dirac species in the actual port reduction. Supplied dimensionless mass units and phi background; normalized vacuum-polarization convention stated in the report. Scalar counterterms and UV boundary conditions remain free. These matched operators are required to equate the96-state description to the original theory below heavy thresholds, not a new UV completion or mass/gravity prediction.')

def recovery_rows():
 H=np.array([[(j>>k)&1 for j in range(1,8)] for k in range(3)]);bits=np.array([[(j>>k)&1 for k in range(7)] for j in range(128)])
 words={sum(int(a)<<i for i,a in enumerate((np.array(v)@H)%2)) for v in product(range(2),repeat=3)}
 V=np.zeros((128,2));V[list(words),0]=1/np.sqrt(8);V[[w^127 for w in words],1]=1/np.sqrt(8);R=[]
 # Six-bit CSS syndrome; each nonzero three-bit value is a Hamming column.
 for sz,sx in product(range(8),repeat=2):
  z=np.zeros(7,int);x=z.copy()
  if sx:z[sx-1]=1
  if sz:x[sz-1]=1
  # row Vdag E^-1; product convention E=X^x Z^z, phases cancel in channel.
  perm=np.arange(128)^sum(int(a)<<i for i,a in enumerate(x));phase=(-1.)**(bits@z)
  R.append((V.T[:,perm]*phase[None,:]))
 R=np.array(R);assert np.linalg.norm(np.einsum('sab,sac->bc',R.conj(),R)-np.eye(128))<1e-12
 return V,R

def marginal_noise(relay,p):
 C=P.endpoint_channel(relay,'damping',p);U=np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
 S=np.zeros((2,2,2,2),complex)
 for i,j in product(range(2),repeat=2):
  rho=sum(2*C[(2*i+k)*4:(2*i+k+1)*4,(2*j+k)*4:(2*j+k+1)*4] for k in range(2));rho=U.T@rho@U
  S[:,:,i,j]=np.einsum('ikjk->ij',rho.reshape(2,2,2,2))
 return S

def apply_product(S,rho):
 t=rho.reshape([2]*14)
 for wire in range(7):
  out=np.tensordot(S,t,axes=([2,3],[wire,wire+7]));labels=[wire,wire+7]+[i for i in range(14) if i not in [wire,wire+7]];t=out.transpose(np.argsort(labels))
 return t.reshape(128,128)

def leakage_controls():
 # Native phase masks constant on support components preserve the isolated
 # cell exactly at the stored-frame tolerance, unlike generic edge phases.
 M=P.native();old=M.N.P.N.M.P;B=old.read(M.N.P.N.M.previous()['native_noise']['invariant_frame']);Q=B.conj().T@old.read(old.prior()['local_gates']['logical_basis']);D=np.array(old.load()['D'],float)
 C=B@B.conj().T;adj=abs(C)>1e-8;unseen=set(range(160));blocks=[]
 while unseen:
  seed=min(unseen);block=[seed];unseen.remove(seed)
  for i in block:
   for j in np.flatnonzero(adj[i]):
    if int(j) in unseen:unseen.remove(int(j));block.append(int(j))
  blocks.append(sorted(block))
 masks=np.zeros((len(blocks),160))
 for j,b in enumerate(blocks):masks[j,b]=1
 H=np.array([B.conj().T@(m[:,None]*B) for m in masks]);U=np.array([B[min(b)].conj()/np.linalg.norm(B[min(b)]) for b in blocks]).T
 assert U.shape==(12,12) and np.linalg.norm(U.conj().T@U-np.eye(12))<1e-10
 Hp=B.conj().T@(D[:40].T@D[:40])@B;Hl=B.conj().T@(D[40:].T@D[40:])@B
 A=U.conj().T@(Hp+.123*Hl)@U;edges=np.argwhere(np.triu(abs(A)>1e-9,1));seen={0}
 for _ in range(12):
  for i,j in edges:
   if int(i) in seen or int(j) in seen:seen.update([int(i),int(j)])
 assert len(seen)==12
 # Powers of two give pairwise distinct positive energy differences: each
 # is a binary string of ones times a unique power of two.
 energies=[2**i for i in range(12)];gaps=[energies[j]-energies[i] for i,j in edges]
 assert len(set(gaps))==len(gaps)
 cellerr=max(np.linalg.norm(m[:,None]*B-B@h) for m,h in zip(masks,H));logical=Q@Q.conj().T;mix=max(np.linalg.norm(h@logical-logical@h) for h in H)
 assert cellerr<1e-10 and mix>.5
 return dict(status='PASS',edge_phase_masks=[list(map(int,b)) for b in blocks],mask_sizes=list(map(len,blocks)),cell_invariance_error=float(cellerr),maximum_logical_leakage_commutator=float(mix),control_basis=enc(U),connected_transition_edges=edges.tolist(),drift_combination=[1.,.123],diagonal_energies=energies,distinct_addressed_gaps=len(gaps),scope='Explicit160-edge phase masks preserve the actual12-cell and mix its logical/leakage sectors. Combined with native point/line Hamiltonians, connected transitions and distinct diagonal gaps yield full one-cell u12 controllability by polynomial commutator filtering. Independent grouped phase addressing is a supplied new control resource; experimental precision, two-cell entangling implementation, environment reset and fault-tolerant synthesis remain open.')

def decoded():
 V,R=recovery_rows();rows=[]
 for p in [.01,.03,.12]:
  for label,relay in [('0',np.diag([1.,0.])),('1',np.diag([0.,1.]))]:
   S=marginal_noise(relay,p);choi=np.zeros((4,4),complex)
   for i,j in product(range(2),repeat=2):
    rho=apply_product(S,np.outer(V[:,i],V[:,j]));out=np.einsum('sab,bc,sdc->ad',R,rho,R.conj());choi[2*i:2*i+2,2*j:2*j+2]=out/2
   omega=np.array([1,0,0,1])/np.sqrt(2);cp=float(min(np.linalg.eigvalsh((choi+choi.conj().T)/2)));tp=float(np.linalg.norm(np.einsum('iaja->ij',choi.reshape(2,2,2,2))-np.eye(2)/2));physical=S.transpose(2,0,3,1).reshape(4,4)/2
   rows.append(dict(probability=p,relay=label,logical_choi=enc(choi),logical_entanglement_infidelity=float(1-(omega@choi@omega).real),physical_entanglement_infidelity=float(1-(omega@physical@omega).real),minimum_Choi_eigenvalue=cp,TP_error=tp))
 assert min(r['minimum_Choi_eigenvalue'] for r in rows)>-1e-12 and max(r['TP_error'] for r in rows)<1e-12
 return dict(status='PASS',controls=leakage_controls(),rows=rows,recovery_completeness_error=float(np.linalg.norm(np.einsum('sab,sac->bc',R.conj(),R)-np.eye(128))),scope='Full Steane encode/channel/six-syndrome recovery/decode under seven independent copies of an actual dissipative relay marginal. Endpoint CNOT is inverted before tracing a maximally mixed second endpoint to define noise; correlations between the two logical blocks are discarded by this declared marginal model. Every Kraus branch and multiple faults are included, without Pauli twirling. Readout and correction are ideal; no full two-block correlated gate recovery, leakage-control hardware realization or threshold is claimed.')

def run():
 r={}
 for name,fn in [('frame',frame),('current',current),('wave',wave),('matching',matching),('decoded',decoded)]:
  r[name]=fn();print(name,r[name]['status'],flush=True)
 paths=['data/w33_pass11384_native_inputs.json',str(P.OUT.relative_to(ROOT))]
 r.update(status='PASS',passes=list(range(11471,11476)),reservation='c630384c2',source_sha256={p:hashlib.sha256(json.dumps(json.loads((ROOT/p).read_text()),sort_keys=True,separators=(',',':')).encode()).hexdigest() for p in paths})
 OUT.write_text(json.dumps(r,indent=2)+'\n');return r
if __name__=='__main__':run()
