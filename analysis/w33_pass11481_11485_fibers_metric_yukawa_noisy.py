"""Constructed follow-ups, with classical and repository ownership retained."""
import json,hashlib
from pathlib import Path
from itertools import product
import numpy as np
from scipy.linalg import expm,null_space
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/w33_pass11481_11485_fibers_metric_yukawa_noisy.json'
read=lambda p:json.loads((ROOT/p).read_text())
enc=lambda a:dict(real=np.asarray(a).real.tolist(),imag=np.asarray(a).imag.tolist())
dec=lambda a:np.array(a['real'])+1j*np.array(a['imag'])
OLD='data/w33_pass11476_11480_cubic_geometry_correlated.json'
def su3():
 out=[np.diag([1,-1,0])/np.sqrt(2),np.diag([1,1,-2])/np.sqrt(6)]
 for i,j in [(0,1),(0,2),(1,2)]:
  a=np.zeros((3,3),complex);a[i,j]=a[j,i]=1/np.sqrt(2);out.append(a)
  a=np.zeros((3,3),complex);a[i,j]=1j/np.sqrt(2);a[j,i]=-1j/np.sqrt(2);out.append(a)
 return np.array(out)
def exact_fiber_control():
 import sympy as sp
 A=sp.Matrix([[1,sp.I,0],[0,2,1],[1,0,3]]);B=sp.Matrix([[2,1,sp.I],[sp.I,1,0],[1,2,2]])
 L=[sp.diag(1,-1,0),sp.diag(1,1,-2)]
 for i,j in [(0,1),(0,2),(1,2)]:
  x=sp.zeros(3);x[i,j]=x[j,i]=1;L.append(x);x=sp.zeros(3);x[i,j]=sp.I;x[j,i]=-sp.I;L.append(x)
 grams=[A*A.conjugate().T,B*B.conjugate().T]
 C=sp.Matrix.hstack(*[sp.Matrix([sp.re((sp.I*(l*g-g*l))[i,i]) for g in grams for i in range(2)]) for l in L])
 def pack(a,b):return sp.Matrix([sp.re(z) for z in a]+[sp.im(z) for z in a]+[sp.re(z) for z in b]+[sp.im(z) for z in b])
 D=sp.Matrix.hstack(*[pack(sp.I*l*A,sp.I*l*B) for l in L]);gauge=sp.Matrix.hstack(*([pack(sp.I*l*A,sp.I*l*B) for l in L[:2]]+[pack(sp.I*A*l,sp.I*B*l) for l in L]))
 N=sp.Matrix.hstack(*C.nullspace());ranks=[C.rank(),gauge.rank(),gauge.row_join(D*N).rank()]
 assert ranks==[4,10,12] and C*N==sp.zeros(4,4)
 return dict(A=[[str(x) for x in row] for row in A.tolist()],B=[[str(x) for x in row] for row in B.tolist()],constraint_Jacobian=[[str(x) for x in row] for row in C.tolist()],kernel=[[str(x) for x in row] for row in N.tolist()],exact_ranks=ranks,physical_fiber_dimension=2,scope='Exact Gaussian-integer control proves nonempty generic rank-four locus and two fiber tangents modulo the ten internal gauge directions; it is not the numerical stationary vacuum.')

def fiber_restoring(A,B,H):
 L=su3();ops=np.array([np.kron(l,np.eye(3)) for l in L]);gram=np.zeros((8,8));curvature=gram.copy();sensitive=[]
 for a in [A.ravel(),B.ravel()]:
  mu=np.einsum('i,aij,j->a',a.conj(),H,a).real;da=1j*(ops@a);dm=2*np.einsum('xi,aij,j->xa',da.conj(),H,a).real
  dd=np.array([[-(ops[i]@ops[j]+ops[j]@ops[i])@a/2 for j in range(8)] for i in range(8)])
  d2=2*(np.einsum('xi,aij,yj->xya',da.conj(),H,da)+np.einsum('i,aij,xyj->xya',a.conj(),H,dd)).real
  gram+=dm@dm.T;curvature+=np.einsum('a,xya->xy',mu,d2);sensitive.append(float(np.linalg.norm(mu[np.linalg.norm(dm,axis=0)>1e-9])))
 eigen=np.linalg.eigvalsh(gram)
 assert sum(eigen>1e-10)==4 and sum(eigen>1e-3)==2
 return dict(Gram_term=gram.tolist(),moment_curvature_term=curvature.tolist(),generator_path_Hessian=(gram+curvature).tolist(),Gram_eigenvalues=eigen.tolist(),path_eigenvalues=np.linalg.eigvalsh(gram+curvature).tolist(),sensitive_moment_norms=sensitive,moment_curvature_norm=float(np.linalg.norm(curvature)),scope='Analytic second derivative along common-left SU3 exponential paths. Four nonzero Gram restoring values: two order0.4, two order8e-9. Its four-dimensional kernel includes two left-diagonal gauge directions and two physical potential-value fibers. Moment-curvature remainder6.13e-10 and numerical stationarity prevent an exact physical-mass claim; this explains the additional weak pair without declaring four exact moduli.')

def fibers():
 import w33_pass11476_11480_cubic_geometry_correlated as M
 import w33_pass11438_finite_native_model as F
 from scipy.optimize import least_squares
 r=read(OLD)['vacuum'];p=r['polydisc'];A=dec(p['A_at_vacuum']);B=dec(p['B_at_vacuum']);E=dec(r['isotope']['native_tripotents']);H=dec(p['compressed_moment_generators']);S=dec(read('data/w33_pass11466_11470_stationary_ports_channels.json')['soft']['soft_frame']).real
 L=su3();grams=[A@A.conj().T,B@B.conj().T]
 C=np.array([np.concatenate([np.diag(1j*(l@g-g@l)).real[:2] for g in grams]) for l in L]).T
 U,s,Vh=np.linalg.svd(C,full_matrices=True);N=Vh[4:].T;normal=Vh[:4].T
 pack=lambda a,b:F.pack((E@a).reshape(-1),(E@b).reshape(-1))
 D=np.column_stack([pack(1j*l@A,1j*l@B) for l in L])
 gauge=np.column_stack([pack(1j*l@A,1j*l@B) for l in L[:2]]+[pack(1j*A@l,1j*B@l) for l in L])
 q=np.linalg.qr(gauge)[0];physical=D@N-q@(q.T@D@N);uu,ss,vv=np.linalg.svd(physical,full_matrices=False)
 free=N@vv[:2].T;rows=[];base=M.polydisc_value(A,B,H);x0=pack(A,B)
 def constraint(U):return np.concatenate([np.diag(U@g@U.conj().T-g).real[:2] for g in grams])
 with F.native_context():
  for axis,t in product(range(2),[-.003,-.001,.001,.003]):
   def unit(y):return expm(1j*np.einsum('a,aij->ij',t*free[:,axis]+normal@y,L))
   opt=least_squares(lambda y:constraint(unit(y)),np.zeros(4),xtol=1e-14,gtol=1e-14,ftol=1e-14,max_nfev=500)
   u=unit(opt.x);a=u@A;b=u@B;x=pack(a,b);v,g=F.fun(x);chord=x-x0;physical_chord=chord-q@(q.T@chord)
   rows.append(dict(axis=axis,t=t,correction=opt.x.tolist(),unitary=enc(u),constraint_residual=float(np.linalg.norm(constraint(u))),reduced_energy_shift=M.polydisc_value(a,b,H)-base,native_energy_shift=v-base,native_gradient_norm=float(np.linalg.norm(g)),nongauge_chord_norm=float(np.linalg.norm(physical_chord)),soft_tangent_remainder=float(np.linalg.norm((np.eye(324)-S@S.T)@(uu[:,axis])))))
 assert min(s)>1e-5 and sum(ss>1e-7)==2 and max(z['constraint_residual'] for z in rows)<1e-10 and max(abs(z['native_energy_shift']) for z in rows)<1e-8
 return dict(status='PASS',restoring_decomposition=fiber_restoring(A,B,H),exact_generic_control=exact_fiber_control(),constraint_Jacobian=C.tolist(),constraint_singular_values=s.tolist(),physical_kernel_singular_values=ss.tolist(),free_coefficients=free.tolist(),normal_coefficients=normal.tolist(),rows=rows,scope='Two generic physical common-left level-set directions are integrated numerically while preserving both left-Gram diagonals. The implicit-function theorem applies at any exact pair with this rank-four Jacobian; stored numerical vacuum has no interval stationarity certificate. These are flat potential-value fibers, not a proof that all fiber points are stationary or that all four Hessian soft modes are exact moduli. Simultaneous hollowisation is prior art of Damm/Fassbender.')

def metric_action(G,dG,N,u,du,edges,h,cell=1.,coupling=1.):
 vol=cell*np.sqrt(np.linalg.det(G));Q=np.linalg.solve(G,dG)
 scalar_kin=.5*vol*np.dot(du,du);metric_kin=coupling*vol*(np.trace(Q@Q)-np.trace(Q)**2)/8
 weights=vol/np.einsum('ei,ij,ej->e',h,G,h)
 potential=.5*np.sum(weights*(u[edges[:,1]]-u[edges[:,0]])**2)
 return float((scalar_kin+metric_kin)/N-N*potential),float(scalar_kin+metric_kin),float(potential)
def metric_hamiltonian(state,edges,h):
 G=state[:9].reshape(3,3);P=state[9:18].reshape(3,3);u=state[18:98];p=state[98:];vol=np.sqrt(np.linalg.det(G));inv=np.linalg.inv(G);GP=G@P;t=np.trace(GP)
 Tg=2/vol*(np.trace(GP@GP)-t*t/2);Ts=np.dot(p,p)/(2*vol);length=np.einsum('ei,ij,ej->e',h,G,h);w=vol/length;delta=u[edges[:,1]]-u[edges[:,0]];V=.5*np.dot(w,delta**2)
 Vg=.5*V*inv-.5*vol*np.einsum('e,ei,ej->ij',delta**2/length**2,h,h)
 dG=4/vol*G@P@G-2/vol*t*G;dP=.5*(Tg+Ts)*inv-4/vol*P@G@P+2/vol*t*P-Vg
 force=np.zeros(80);np.add.at(force,edges[:,0],w*delta);np.add.at(force,edges[:,1],-w*delta)
 return float(Tg+Ts+V),np.r_[dG.ravel(),dP.ravel(),p/vol,force]

def metric():
 from fractions import Fraction
 r=read(OLD)['geometry'];edges=np.array(r['cover_edges']);h=np.array(r['harmonic_FCC_displacements']);rng=np.random.default_rng(11482)
 G=np.array([[1.2,.03,.02],[.03,.9,-.01],[.02,-.01,1.1]]);dG=np.diag([.02,-.02,0.]);u=rng.normal(size=80)*.01;du=rng.normal(size=80)*.02
 _,K,V=metric_action(G,dG,1.,u,du,edges,h);N=np.sqrt(K/V);value=metric_action(G,dG,N,u,du,edges,h)[0]
 J=np.array([[1.1,.1,0],[0,.9,.08],[.02,0,1.2]]);Ji=np.linalg.inv(J);Gp=Ji.T@G@Ji;dGp=Ji.T@dG@Ji;hp=h@J.T;vp=np.linalg.det(J)
 coordinate=metric_action(Gp,dGp,N,u,du,edges,hp,cell=vp)[0];alpha=1.7
 time=alpha*metric_action(G,dG/alpha,N/alpha,u,du/alpha,edges,h)[0]
 step=1e-7;fd=(metric_action(G,dG,N+step,u,du,edges,h)[0]-metric_action(G,dG,N-step,u,du,edges,h)[0])/(2*step)
 lapse=-K/N**2-V # For L=K/N-NV, positive scalar K makes lapse equation impossible.
 # Correct finite constraint diagnosis: positive shear and scalar energy cannot
 # solve lapse without negative conformal kinetic energy (or extra terms).
 from scipy.optimize import brentq
 f=lambda a:sum(metric_action(G,dG+a*G,1.,u,du,edges,h)[1:])
 a=brentq(f,0,10);action,K2,V2=metric_action(G,dG+a*G,1.,u,du,edges,h)
 fd2=(metric_action(G,dG+a*G,1.+step,u,du,edges,h)[0]-metric_action(G,dG+a*G,1.-step,u,du,edges,h)[0])/(2*step)
 assert abs(coordinate-value)<1e-12 and abs(time-value)<1e-12 and abs(fd-lapse)<1e-8 and abs(K2+V2)<1e-10 and abs(fd2)<1e-8
 from scipy.integrate import solve_ivp
 velocity=dG+a*G;vol=np.sqrt(np.linalg.det(G));inv=np.linalg.inv(G);Q=inv@velocity;P=vol/4*(inv@velocity@inv-np.trace(Q)*inv)
 state=np.r_[G.ravel(),P.ravel(),u,vol*du];sol=solve_ivp(lambda t,y:metric_hamiltonian(y,edges,h)[1],[0,.002],state,t_eval=np.linspace(0,.002,5),rtol=1e-10,atol=1e-12)
 history=[dict(time=float(t),constraint=metric_hamiltonian(y,edges,h)[0],minimum_metric_eigenvalue=float(min(np.linalg.eigvalsh(y[:9].reshape(3,3)))),state=y.tolist()) for t,y in zip(sol.t,sol.y.T)]
 assert sol.success and max(abs(z['constraint']) for z in history)<1e-7 and min(z['minimum_metric_eigenvalue'] for z in history)>0
 return dict(status='PASS',Hamiltonian_history=history,initial_canonical_state=state.tolist(),graph_parameters=r['graph_parameters'],edges=edges.tolist(),harmonic_displacements=h.tolist(),metric=G.tolist(),metric_velocity=dG.tolist(),scalar=u.tolist(),scalar_velocity=du.tolist(),coordinate_transform=J.tolist(),coordinate_covariance_residual=abs(coordinate-value),time_reparametrization_residual=abs(time-value),positive_kinetic_lapse_force=lapse,positive_kinetic_FD_force=fd,conformal_velocity=a,constrained_metric_velocity=(dG+a*G).tolist(),constrained_kinetic=K2,potential=V2,lapse_constraint_residual=abs(K2+V2),lapse_FD_residual=abs(fd2),scope='Explicit supplied homogeneous dynamical SPD metric and lapse, with coupled Hamiltonian scalar/metric evolution on actual periodic-cover edges, with DeWitt kinetic term and scalar spring action. Affine spatial-coordinate and time-reparametrization covariance are verified. Lapse variation gives K+N^2 V=0: positive scalar/shear kinetics alone fail; a negative conformal direction satisfies this finite constraint. No spatial diffeomorphism invariance, intrinsic Einstein action, local gravitational constraints, physical G or continuum derivation.')

def yukawa(fiber_rows):
 from w33_pass11384_11388_native_dynamics import native_tensors
 _,tensors=native_tensors();tri=tensors[0];d=np.array(tensors[1],complex)
 # tensors[1] is the native cubic determinant coefficient tensor.
 eps=np.zeros((3,3,3),int)
 for p in [(0,1,2),(1,2,0),(2,0,1)]:eps[p]=1;eps[p[1],p[0],p[2]]=-1
 T=np.einsum('abc,ijk->aibjck',d,eps).reshape(81,81,81)
 antisym=np.linalg.norm(T+T.transpose(1,0,2));sym=np.linalg.norm((T+T.transpose(1,0,2))/2)
 old=read(OLD);E=dec(old['vacuum']['isotope']['native_tripotents']);A=dec(old['vacuum']['polydisc']['A_at_vacuum']);fields=E@A
 zero=np.einsum('ijk,k->ij',(T+T.transpose(1,0,2))/2,fields.reshape(-1))
 # A scalar H_c^{ij} in (27,bar6) supplies the missing symmetric family tensor.
 # Build from existing scalar via S^{ij}_k symmetric in family indices: an
 # explicit additional spurion, not a consequence of the original field.
 spurion=np.zeros((3,3,3),complex);spurion[0]=np.diag([1.,0,0]);spurion[1]=np.diag([0,1.,0]);spurion[2]=np.diag([0,0,1.])
 H=np.einsum('ck,kij->cij',fields,spurion);Y=np.einsum('abc,cij->aibj',d,H).reshape(81,81)
 eig=np.linalg.svd(Y,compute_uv=False);mu=10.;m=eig[eig>1e-12]
 # Supplied Weyl MSbar one-loop value (two spin degrees), all81 masses kept.
 cw=-np.sum(m**4*(np.log(m*m/mu**2)-1.5))/(32*np.pi**2)
 portal_rows=[]
 for row in fiber_rows:
  U=dec(row['unitary']);HH=np.einsum('ck,kij->cij',E@U@A,spurion);YY=np.einsum('abc,cij->aibj',d,HH).reshape(81,81);mm=np.linalg.svd(YY,compute_uv=False);positive=mm[mm>1e-12];value=-np.sum(positive**4*(np.log(positive**2/mu**2)-1.5))/(32*np.pi**2)
  portal_rows.append(dict(axis=row['axis'],t=row['t'],tree_shift=row['native_energy_shift'],Weyl_loop_shift=float(value-cw)))
 loop_derivatives=[]
 for axis,step in product(range(2),[.001,.003]):
  plus=next(z['Weyl_loop_shift'] for z in portal_rows if z['axis']==axis and z['t']==step);minus=next(z['Weyl_loop_shift'] for z in portal_rows if z['axis']==axis and z['t']==-step)
  loop_derivatives.append(dict(axis=axis,step=step,first=(plus-minus)/(2*step),second=(plus+minus)/step**2))
 assert antisym<1e-12 and sym<1e-12 and np.linalg.norm(Y-Y.T)<1e-12 and np.linalg.norm(zero)<1e-12
 return dict(status='PASS',fiber_loop_derivatives=loop_derivatives,fiber_portal_rows=portal_rows,native_cubic=enc(d),family_epsilon=eps.tolist(),exchange_symmetrization_norm=float(sym),forbidden_mass_matrix_norm=float(np.linalg.norm(zero)),spurion=enc(spurion),enlarged_Higgs=enc(H),allowed_mass_matrix=enc(Y),allowed_singular_values=eig.tolist(),rank=int(sum(eig>1e-9)),supplied_Weyl_CW_value=float(cw),scope='One identical left-Weyl species in (27,3) with a scalar (27,3): symmetric E6 cubic times antisymmetric family epsilon has zero symmetric fermion-bilinear projection, so this renormalizable self-Yukawa vanishes. Distinct Weyl species, larger E6 tensors or family bar6 Higgs evade this. Explicit (27,bar6) spurion construction produces a mass matrix; spurion, coupling and MSbar scale are supplied; the map from the original triplet to bar6 is not an intertwiner without transforming its additional spurion. No SM projection, observed masses or unique matching conditions; Pass11271 owns the identical-species zero and (27,bar6) repair; Pass11281 owns the sextet selection obstruction. This extends their map to the actual11466 numerical vacuum, not a new Yukawa selection theorem.')

def local_ward():
 import w33_pass11428_11432_native_alignment_measure_coarse_flags as Q
 import w33_pass11476_11480_cubic_geometry_correlated as M
 L=2;n=16;B=M.gauge_boundary(L);rng=np.random.default_rng(11484);A=.0005*rng.normal(size=64);directions=np.eye(64).reshape(64,16,4);J=np.zeros((16,64));a=np.zeros(16)
 for q,m in zip([1,-4,2,-3,6],[6,3,3,2,1]):
  P,ds,*_=Q.weak_overlap(q,1.,A.reshape(16,4),directions,L=2)
  a+=m*q*np.trace(P.reshape(16,4,16,4),axis1=1,axis2=3).diagonal().real
  J+=m*q*np.trace(np.asarray(ds).reshape(64,16,4,16,4),axis1=2,axis2=4).diagonal(axis1=1,axis2=2).real.T
 sites=np.array(list(product(range(2),repeat=4)));rows=[]
 for radius in range(5):
  flows=np.zeros((64,64));res=[];tails=[]
  for ell in range(64):
   root=ell//4;distance=np.sum(abs(sites-sites[root]),axis=1)
   keep=np.array([j for j in range(64) if max(distance[j//4],distance[np.flatnonzero(B[j]<0)[0]])<=radius])
   if len(keep):flows[keep,ell]=np.linalg.lstsq(B[keep].T,-J[:,ell],rcond=None)[0]
   res.append(np.linalg.norm(B.T@flows[:,ell]+J[:,ell]));tails.append(np.linalg.norm(J[distance>radius,ell]))
  rows.append(dict(radius=radius,maximum_divergence_error=float(max(res)),maximum_omitted_anomaly_derivative=float(max(tails)),flow_norm=float(np.linalg.norm(flows)),gauge_derivative_defect=float(np.linalg.norm(flows@B)),flows=flows.tolist()))
 assert np.linalg.norm(J.sum(axis=0))<1e-10 and np.linalg.norm(J@B)<1e-10 and rows[-1]['maximum_divergence_error']<1e-10
 return dict(status='PASS',L=L,background=A.tolist(),anomaly_density=a.tolist(),anomaly_Jacobian=J.tolist(),gauge_invariance_of_anomaly_derivative=float(np.linalg.norm(J@B)),rows=rows,scope='Explicit support-restricted flow solves the differentiated Ward equation increasingly over finite balls. Full-radius divergence is numerical roundoff; gauge derivative defect is separately retained. This is a construction attempt at the local anomaly-divergence primitive, not a local chiral measure: L2 cannot bound infinite-volume exponential tails, flows need not be gauge-invariant, and curvature/axial symmetry/flux sectors remain open. Luscher equation5.8 requires a local gauge-invariant anomaly primitive before current reconstruction; previous Coulomb connection does not supply it.')

def noisy_recovery():
 import w33_pass11471_11475_frames_currents_matching as Q
 import w33_pass11476_11480_cubic_geometry_correlated as M
 V,R=Q.recovery_rows();states=[];rows=[];PA=M.PAULI
 S=Q.marginal_noise(np.diag([1.,0.]),.01)
 images=[Q.apply_product(S,V@p@V.T) for p in PA]
 syndrome_images=np.array([[r@rho@r.conj().T for rho in images] for r in R])
 priors=np.trace(syndrome_images[:,0],axis1=1,axis2=2).real/2
 bits=np.array([[(s>>b)&1 for b in range(6)] for s in range(64)])
 for error in [0.,.001,.01,.05]:
  dist=np.sum(bits[:,None,:]!=bits[None,:,:],axis=2);confusion=error**dist*(1-error)**(6-dist)
  # A final ideal code-space verification accepts only correctly identified
  # true syndromes. Rejected branches become a flagged erasure, not renormalized.
  for adaptive in [False,True]:
   guess=np.argmax(confusion*priors[None,:],axis=1) if adaptive else np.arange(64)
   acceptance=np.array([sum(confusion[t,s] for t in range(64) if guess[t]==s) for s in range(64)])
   logical=np.einsum('s,spij->pij',acceptance,syndrome_images);out=np.zeros((4,3,3),complex);out[:,:2,:2]=logical
   for k in range(4):out[k,2,2]=2*(k==0)-np.trace(logical[k])
   J=sum(np.kron(PA[k].T,out[k]) for k in range(4))/4;ev=np.linalg.eigvalsh(J)
   success=float(np.trace(logical[0]).real/2);target=np.zeros(6);target[0]=target[4]=1/np.sqrt(2);fidelity=float(np.vdot(target,J@target).real)
   assert min(ev)>-1e-10 and abs(np.trace(J)-1)<1e-10
   rows.append(dict(readout_bit_error=error,adaptive=adaptive,guess=guess.tolist(),acceptance=acceptance.tolist(),success_probability=success,unconditional_entanglement_fidelity=fidelity,minimum_Choi_eigenvalue=float(min(ev)),flagged_Choi=enc(J)))
 # Native12-cell phase controls provide a separate explicit leakage witness.
 ctrl=Q.leakage_controls();native=Q.P.native();old=native.N.P.N.M.P
 B=old.read(native.N.P.N.M.previous()['native_noise']['invariant_frame']);logical=B.conj().T@old.read(old.prior()['local_gates']['logical_basis']);P=logical@logical.conj().T
 masks=ctrl['edge_phase_masks'];best=None
 for block in masks:
  h=B.conj().T@(np.isin(np.arange(160),block)[:,None]*B)
  u=expm(1j*.1*h);leak=(np.eye(12)-P)@u@logical;rate=np.linalg.norm(leak)**2/2
  if best is None or rate>best['mean_leakage']:best=dict(block=block,phase=.1,cell_unitary=enc(u),mean_leakage=float(rate),logical_subchannel=enc(logical.conj().T@u@logical),echo_residual=float(np.linalg.norm(u.conj().T@u@logical-logical)))
 assert best['mean_leakage']>1e-5
 return dict(status='PASS',syndrome_priors=priors.tolist(),rows=rows,native_leakage=best,scope='Untwirled native relay marginal on a full Steane block, explicit six-bit readout confusion, MAP syndrome inference and ideal final code verification. Failed branches are retained as erasures; MAP maximizes acceptance under maximally mixed logical prior, not necessarily unconditional quantum fidelity. Correlated two-block recovery remains owned by11480 and is not approximated here as independent marginals. Native grouped edge phase creates quantified logical leakage inside12-cell; noiseless inverse echo is a declared control. No noisy verification circuit, all-fault correlated recovery or hardware threshold.')

def run():
 results={}
 for name,fn in [('fibers',fibers),('metric',metric),('yukawa',lambda:yukawa(results['fibers']['rows'])),('ward',local_ward),('recovery',noisy_recovery)]:
  results[name]=fn();print(name,'PASS',flush=True);Path('/tmp/w33_11481_partial.json').write_text(json.dumps(results))
 sources=[OLD,'data/w33_pass11466_11470_stationary_ports_channels.json','data/w33_pass11471_11475_frames_currents_matching.json'];results.update(status='PASS',passes=list(range(11481,11486)),reservation='ef90a0dc7',source_sha256={p:hashlib.sha256(json.dumps(read(p),sort_keys=True,separators=(',',':')).encode()).hexdigest() for p in sources})
 OUT.write_text(json.dumps(results,indent=2)+'\n');return results
if __name__=='__main__':run()
