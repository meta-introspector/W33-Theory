"""Constructed five-frontier controls; imported physical methods retain ownership.

11461 owns native fields/charges,11428 owns quark network,11443 owns scalar
Lorentz map,10956 owns a different exact Albert Spin8. This packet does not
identify the two carriers, derive couplings or complete a chiral measure.
"""
import json,hashlib
from pathlib import Path
from itertools import product
import numpy as np
from scipy.linalg import block_diag,null_space
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/w33_pass11466_11470_stationary_ports_channels.json'
enc=lambda a:dict(real=np.asarray(a).real.tolist(),imag=np.asarray(a).imag.tolist())
dec=lambda a:np.array(a['real'])+1j*np.array(a['imag'])
def prior():return json.loads((ROOT/'data/w33_pass11461_11465_explicit_physical_maps.json').read_text())
def native():
 import w33_pass11461_11465_explicit_physical_maps as M
 return M

def vacuum():
 import w33_pass11438_finite_native_model as F
 old=np.array(prior()['photon_candidate']['coordinates'])
 def orbit(x):
  a=x[:81]+1j*x[81:162];b=x[162:243]+1j*x[243:]
  return np.concatenate([(1j*(F.H@a)).real,(1j*(F.H@a)).imag,(1j*(F.H@b)).real,(1j*(F.H@b)).imag],axis=1).T
 z=np.linalg.svd(orbit(old),full_matrices=True)[2][58:].T
 K=np.einsum('ab,aij->bij',z,F.H);w,E=np.linalg.eigh(sum(k@k for k in K));E=E[:,w<1e-9]
 assert E.shape==(81,9)
 R=np.block([[E.real,-E.imag],[E.imag,E.real]]);T=block_diag(R,R);y=T.T@old;rows=[]
 with F.native_context():
  for it in range(3):
   v,g=F.fun(T@y);rows.append(dict(value=v,gradient_norm=float(np.linalg.norm(g))))
   h=1e-5;H=np.column_stack([T.T@(F.fun(T@(y+h*d))[1]-F.fun(T@(y-h*d))[1])/(2*h) for d in np.eye(36)]);H=(H+H.T)/2
   ew,V=np.linalg.eigh(H);y-=V@np.divide(V.T@(T.T@g),ew,out=np.zeros(36),where=abs(ew)>1e-3)
  x=T@y;v,g=F.fun(x);O=orbit(x);sv=np.linalg.svd(O,compute_uv=False)
  h=2e-5;H=np.column_stack([(F.fun(x+h*d)[1]-F.fun(x-h*d)[1])/(2*h) for d in np.eye(324)]);H=(H+H.T)/2
  Z=np.linalg.svd(O,full_matrices=True)[0][:,58:];ew,V=np.linalg.eigh(Z.T@H@Z)
  d=Z@V[:,0];check=(F.fun(x+1e-5*d)[1]-F.fun(x-1e-5*d)[1])/2e-5
 K=K[:,::3,::3];G=np.einsum('aij,bji->ab',K,K).real;gw,gv=np.linalg.eigh(G);K=np.einsum('ab,aij->bij',gv/np.sqrt(gw),K)
 flat=K.reshape(28,-1).T;ads=[];err=[]
 for a in K:
  bracket=np.array([(1j*(a@b-b@a)).ravel() for b in K]).T;c=np.linalg.lstsq(flat,bracket,rcond=None)[0];err.append(np.linalg.norm(bracket-flat@c));ads.append(c.real)
 ads=np.array(ads);cartan=null_space(np.einsum('a,aij->ij',np.arange(1,29),ads),rcond=1e-9);A=np.einsum('ab,aij->bij',cartan,ads)
 rootew,U=np.linalg.eigh(1j*np.einsum('a,aij->ij',np.arange(1,len(cartan.T)+1),A));roots=np.array([[np.vdot(U[:,j],1j*t@U[:,j]).real for t in A] for j in range(28)])
 lengths=np.sum(roots**2,axis=1);nz=lengths>1e-9
 assert cartan.shape[1]==4 and sum(nz)==24 and np.ptp(lengths[nz])<1e-9 and max(err)<1e-10
 assert np.linalg.norm(g)<1e-6 and sum(sv>1e-10)==58 and sum(ew < -1e-3)==0
 return dict(status='PASS',coordinates=x.tolist(),generators=enc(K),fixed_frame=enc(E),value=v,gradient_norm=float(np.linalg.norm(g)),Newton_history=rows,orbit_singular_values=sv.tolist(),threshold_ranks=[int(sum(sv>t)) for t in [1e-8,1e-9,1e-10]],closure_error=float(max(err)),Cartan_dimension=4,roots=roots.tolist(),root_squared_lengths=lengths.tolist(),representation_Casimir=np.linalg.eigvalsh(sum(k@k for k in K)).tolist(),normal_eigenvalues=ew.tolist(),lowest_Hessian_replay_error=float(np.linalg.norm(check-H@d)),scope='Constructed numerical native D4 stratum: compact rank4,24 equal-length roots,28 generators,3 fixed coordinates in27 and9 in81. Gauge rank58 is stable across three cutoffs after projection. Lie closure is checked at roundoff, not symbolically or by intervals; stationary residual and finite-difference Hessian remain numerical. No exact stationary symmetry or relaxed quartic/global stability theorem.10956 owns a separate exact Albert Spin8; no carrier intertwiner is supplied.')


def soft_relaxation(vac):
 import w33_pass11438_finite_native_model as F
 from scipy.optimize import minimize
 x=np.array(vac['coordinates']);E=dec(vac['fixed_frame']);R=np.block([[E.real,-E.imag],[E.imag,E.real]]);T=block_diag(R,R)
 a=x[:81]+1j*x[81:162];b=x[162:243]+1j*x[243:]
 O=np.concatenate([(1j*(F.H@a)).real,(1j*(F.H@a)).imag,(1j*(F.H@b)).real,(1j*(F.H@b)).imag],axis=1).T
 U,sv,Vh=np.linalg.svd(T.T@O,full_matrices=True);rank=int(sum(sv>1e-9));Q=T@U[:,rank:]
 with F.native_context():
  value,g=F.fun(x);h=1e-5;H=np.column_stack([Q.T@(F.fun(x+h*d)[1]-F.fun(x-h*d)[1])/(2*h) for d in Q.T]);H=(H+H.T)/2
  w,V=np.linalg.eigh(H);soft=Q@V[:,:4];hard=Q@V[:,4:];D=hard/np.sqrt(w[4:]);rows=[]
  rng=np.random.default_rng(11466);directions=list(np.eye(4))+[c/np.linalg.norm(c) for c in rng.normal(size=(2,4))]
  for j,c in enumerate(directions):
   d=soft@c
   for t in [.06,.12]:
    for sign in [-1,1]:
     def obj(y):
      vv,gg=F.fun(x+sign*t*d+D@y);return vv-value,D.T@gg
     opt=minimize(obj,np.zeros(D.shape[1]),jac=True,method='L-BFGS-B',options=dict(maxiter=120,gtol=1e-9,ftol=1e-15,maxls=40))
     xx=x+sign*t*d+D@opt.x;vv,gg=F.fun(xx)
     rows.append(dict(direction=j,coefficients=c.tolist(),amplitude=t,sign=sign,energy_change=float(vv-value),hard_gradient_norm=float(np.linalg.norm(hard.T@gg)),soft_gradient_norm=float(np.linalg.norm(soft.T@gg)),iterations=int(opt.nit),coordinates=xx.tolist()))
 assert rank==10 and len(w)==26 and all(np.isfinite(r['energy_change']) for r in rows)
 return dict(status='PASS',reduced_gauge_rank=rank,reduced_normal_eigenvalues=w.tolist(),soft_frame=enc(soft),rows=rows,scope='Line-searched massive relaxation at fixed four-soft-mode amplitudes inside the native D4 fixed space, including mixed directions. An unprotected frozen-Hessian iteration can diverge and is not used as evidence. Finite energies and hard/soft residuals are reported, not replaced by solver success flags. No uniform quartic positivity or exact flat-moduli theorem.')

def flux_sectors():
 M=native();L=34;u,v=M.N.P.N.uniform_flux_links(L)
 plaq=lambda a,b:a*np.roll(b,-1,axis=0)*np.roll(a.conj(),-1,axis=1)*b.conj()
 charges=sorted(set(r['key'][2] for r in prior()['charges']['branching']));rows=[]
 for q in charges:
  p=plaq(u**q,v**q);rows.append(dict(charge=q,flux=float(np.angle(p).sum()/2/np.pi),maximum_plaquette=float(np.max(abs(1-p)))))
 assert all(r['maximum_plaquette']<1/30 and abs(r['flux']-r['charge'])<1e-9 for r in rows)
 path=[]
 for t in np.linspace(0,1,201):
  p=plaq(np.exp(1j*t*np.angle(u)),np.exp(1j*t*np.angle(v)))
  path.append(dict(t=float(t),flux=float(np.angle(p).sum()/2/np.pi),maximum_plaquette=float(np.max(abs(1-p)))))
 assert abs(path[0]['flux'])<1e-9 and abs(path[-1]['flux']-1)<1e-9 and max(r['maximum_plaquette'] for r in path)>1.9
 return dict(status='PASS',L=L,charges=rows,path=path,scope='All native integer6Y charges on admissible nonzero flux links. Known Luscher topology: disconnected admissible flux sectors have no continuous transitions within this configuration space. The explicit link interpolation crosses inadmissible fields. A local chiral current and independent sector phase conventions remain unbuilt; a Berry chart cannot supply them. This corrects the requested transition target, not prior certificates.')

def frustum(a,b,h,p,q,G=1.,coupling=1.):
 import mpmath as mp
 d=b-a;area=(a+b)*mp.sqrt(h*h-d*d/2)/2
 gravity=6*(a*a-b*b)*mp.asinh(d/mp.sqrt(4*h*h-d*d))+12*area*mp.asin(d*d/(4*h*h-d*d))
 matter=(a+b)**3*(q-p)**2/(16*h)
 return gravity/(8*mp.pi*G)+coupling*matter

def cosmology():
 import mpmath as mp
 from scipy.optimize import least_squares
 with mp.workdps(40):
  def action(y):
   s,h0,h1,p,q=y
   return frustum(mp.mpf(1),s,h0,mp.mpf(0),p)+frustum(s,mp.mpf('1.2'),h1,p,q)
  def equations(z):
   s,h1,p,q=z;y=list(map(lambda x:mp.mpf(float(x)),[s,.3,h1,p,q]))
   if min(float(y[1]**2-(y[0]-1)**2/2),float(y[2]**2-(mp.mpf('1.2')-y[0])**2/2))<=0:return np.ones(4)*100
   return np.array([float(mp.diff(lambda t:action(y[:j]+[t]+y[j+1:]),y[j])) for j in [0,1,2,3]])
  sol=least_squares(equations,[1.095,.2,.045,.09],bounds=([1.001,.15,-1,.001],[1.199,3,1,2]),xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=150)
  s,h1,p,q=sol.x;y=list(map(lambda x:mp.mpf(float(x)),[s,.3,h1,p,q]));force=equations(sol.x)
  # Boundary q is solved for compatibility, then fixed for all variation tests.
  check=[]
  for j in [0,1,2,3]:
   h=mp.mpf('1e-5');yp=y.copy();ym=y.copy();yp[j]+=h;ym[j]-=h;check.append(float((action(yp)-action(ym))/(2*h)))
  w0=(1+s)**3/(16*.3);w1=(s+1.2)**3/(16*h1);currents=[w0*p,w1*(q-p)]
  deficits=[4*np.arcsin((s-1)**2/(4*.3**2-(s-1)**2)),4*np.arcsin((1.2-s)**2/(4*h1**2-(1.2-s)**2))]
 assert max(abs(force))<1e-10 and max(abs(np.array(check)-force))<1e-7 and abs(currents[0]-currents[1])<1e-10 and min(deficits)>0
 return dict(status='PASS',interior_scale=s,heights=[.3,h1],scalar_values=[0.,p,q],boundary_scales=[1.,1.2],forces=force.tolist(),finite_difference_forces=check,scalar_current=currents,timelike_trapezoid_deficits=deficits,causal_margins=[.3**2-(s-1)**2/2,h1**2-(1.2-s)**2/2],G=1.,matter_coupling=1.,Lambda=0.,scope='Imported Jercher-Steinhaus regular Lorentzian two-frustum action, equations19/31/44. Both lapse equations, interior geometry equation and scalar equation solved simultaneously, with nonzero timelike hinge deficits. Boundary scalar separation is determined for compatibility and then fixed; first lapse .3 parametrizes the search but its equation is enforced. Homogeneous global variations, supplied G and coupling, Lambda0. No W33-to-frustum map, local inhomogeneous solution, observed constants or continuum derivation.')

def mediators():
 import sympy as sp
 M=native();A=np.array(M.N.P.N.M.P.load()['A'],int);I=np.eye(80,dtype=int)
 assert not np.any(A@(A@A-6*I)@(A@A-16*I))
 basis=[];exact=[];rows=[]
 for j in [44,1,0]:
  J=I[:,[0,j]];S=np.concatenate([np.linalg.matrix_power(A,k)@J for k in range(5)],axis=1);X=sp.Matrix(S.tolist());cols=X.columnspace();Q=sp.Matrix.hstack(*cols);T=(Q.T*Q).inv()*Q.T*sp.Matrix(A.tolist())*Q
  assert sp.Matrix(A.tolist())*Q==Q*T
  B=np.linalg.qr(np.array(Q,float))[0];basis.append(B);exact.append(dict(endpoint=j,integer_basis=np.array(Q,int).tolist(),rational_restriction=[[str(x) for x in row] for row in T.tolist()],characteristic_polynomial=str(T.charpoly().as_expr())))
  rows.append(dict(endpoint=j,dimension=Q.cols,invariance_error=float(np.linalg.norm(A@B-B@(B.T@A@B)))))
 assert [r['dimension'] for r in rows]==[8,8,5]
 V=np.zeros((240,21));start=0
 for family,B in enumerate(basis):
  V[family::3,start:start+B.shape[1]]=B;start+=B.shape[1]
 W=block_diag(np.eye(3),V,np.eye(3),V);Pin=M.N.P.N.M.P.read(M.N.P.N.M.previous()['pin_vacua']['pin_down']);scans=[]
 for phi in [.2,.81,1.5]:
  K=np.kron(10*np.eye(80)-phi*A,np.eye(3));left=np.zeros((240,3));left[:3]=np.eye(3)/np.sqrt(2);right=np.zeros((3,240),complex);right[:,:3]=Pin
  H=np.block([[K,left,np.zeros_like(K)],[np.zeros((3,240)),10*np.eye(3),right],[np.zeros_like(K),np.zeros((240,3)),K]])
  AL=np.zeros((3,483));AR=np.zeros((483,3))
  for j,i in enumerate([44,1,0]):AL[j,i*3+j]=1;AR[243+i*3+j,j]=1
  full=np.block([[np.zeros((3,3)),AL],[AR,H]]);small=W.conj().T@full@W
  assert small.shape==(48,48)
  invariant=float(np.linalg.norm(full@W-W@small));schur=-AL@np.linalg.solve(H,AR);schursmall=small[:3,:3]-small[:3,3:]@np.linalg.solve(small[3:,3:],small[3:,:3])
  dark=np.r_[np.repeat(10-phi*np.sqrt(6),134),np.repeat(10,170),np.repeat(10+phi*np.sqrt(6),134)]
  spectra=np.sort(np.r_[np.linalg.svd(small,compute_uv=False),dark]);fullspec=np.sort(np.linalg.svd(full,compute_uv=False))
  logfactor=float(np.linalg.slogdet(full)[1]-np.linalg.slogdet(small)[1]-np.log(dark).sum())
  scans.append(dict(phi=phi,invariance_error=invariant,Schur_error=float(np.linalg.norm(schur-schursmall)),spectrum_error=float(np.max(abs(spectra-fullspec))),logdet_factorization_error=logfactor,reduced_masses=np.linalg.svd(small,compute_uv=False).tolist()))
 assert max(r['spectrum_error'] for r in scans)<1e-10 and max(r['Schur_error'] for r in scans)<1e-10
 return dict(status='PASS',minimal_polynomial='x*(x^2-6)*(x^2-16)',exact_bases=exact,rows=rows,scans=scans,per_sector_sizes=[486,48],two_sector_species=[972,96],two_sector_b0=[-637.,-53.],discarded_species=876,discarded_det_both_sectors='10^340*(100-6*phi^2)^268',discarded_logdet_force_at_081=float(-3216*.81/(100-6*.81**2)),scope='Exact port-invariant reduction of the actual signed-pin486-state mass network. Integer Krylov bases certify dimensions8,8,5 and rational restrictions; SVD replay retains all discarded singular values. Minimal per-family all-port realization by controllability and symmetric observability; this is not a minimal physical UV theory. Original Grassmann determinant and876 gauge-charged species remain. A new theory retaining only96 fermions changes the determinant and has b0=-53, still not asymptotically free; no UV completion.')

def cnot(c,t):
 U=np.zeros((8,8))
 for j in range(8):
  bits=[(j>>k)&1 for k in [2,1,0]];bits[t]^=bits[c];i=4*bits[0]+2*bits[1]+bits[2];U[i,j]=1
 return U
OPS=[(0,1),(1,0),(0,1),(1,2),(0,1),(1,0),(0,1)]
I2=np.eye(2);X=np.array([[0,1],[1,0]]);Z=np.diag([1,-1]);Y=1j*X@Z;PAULI=[I2,X,Y,Z]
def embed_local(A,j):return np.kron(np.kron(A if j==0 else I2,A if j==1 else I2),A if j==2 else I2)
def endpoint_channel(relay,kind,p=.12):
 out=np.zeros((16,16),complex)
 for i,j in product(range(4),repeat=2):
  rho=np.zeros((8,8),complex)
  for r,s in product(range(2),repeat=2):rho[(i//2)*4+r*2+i%2,(j//2)*4+s*2+j%2]=relay[r,s]
  for c,t in OPS:
   U=cnot(c,t);rho=U@rho@U.T
   if kind=='pauli':
    noise=np.zeros_like(rho)
    for a,b in product(range(4),repeat=2):
     weight=1-p if a==b==0 else p/15;P=embed_local(PAULI[a],c)@embed_local(PAULI[b],t);noise+=weight*P@rho@P.conj().T
    rho=noise
   elif kind=='damping':
    for wire in [c,t]:
     K=[embed_local(np.diag([1,np.sqrt(1-p)]),wire),embed_local(np.array([[0,np.sqrt(p)],[0,0]]),wire)]
     rho=sum(k@rho@k.conj().T for k in K)
  reduced=np.einsum('arbcrd->abcd',rho.reshape(2,2,2,2,2,2)).reshape(4,4)
  out[i*4:(i+1)*4,j*4:(j+1)*4]=reduced/4
 return out

def relay_channels():
 # Each complete fault trajectory equals P_end tensor P_relay after the
 # ideal endpoint CNOT tensor I. Partial trace removes P_relay even on an
 # initially correlated state; only the endpoint marginal survives.
 zero=np.diag([1.,0.]);one=np.diag([0.,1.]);plus=np.ones((2,2))/2;rows=[]
 for kind in ['pauli','damping']:
  C0=endpoint_channel(zero,kind);C1=endpoint_channel(one,kind);Cp=endpoint_channel(plus,kind)
  tp=np.einsum('iaja->ij',C0.reshape(4,4,4,4));w=np.linalg.eigvalsh((C0+C0.conj().T)/2)
  rows.append(dict(noise=kind,probability=.12,choi_zero=enc(C0),choi_one=enc(C1),choi_plus=enc(Cp),zero_one_difference=float(np.linalg.norm(C0-C1)),zero_plus_difference=float(np.linalg.norm(C0-Cp)),minimum_Choi_eigenvalue=float(min(w)),TP_error=float(np.linalg.norm(tp-np.eye(4)/4))))
 assert rows[0]['zero_one_difference']<1e-12 and rows[0]['zero_plus_difference']<1e-12 and rows[1]['zero_one_difference']>1e-3
 assert max(r['TP_error'] for r in rows)<1e-12 and min(r['minimum_Choi_eigenvalue'] for r in rows)>-1e-12
 return dict(status='PASS',operations=[list(x) for x in OPS],rows=rows,scope='Full noisy endpoint Choi channels after every routed CNOT, allowing simultaneous faults. State-independent stochastic Pauli trajectory theorem holds for arbitrary classical correlations in their probabilities and initially correlated data/relay states: endpoint output depends only on data marginal. Amplitude damping provides a complete CPTP counterexample to relay-state independence. Physical leakage-control realization, full decoder and threshold remain open; no claim that all noise can be twirled exactly.')

def run():
 results={}
 for name,fn in [('vacuum',vacuum),('flux',flux_sectors),('gravity',cosmology),('mediators',mediators),('relay',relay_channels)]:
  results[name]=fn();print(name,results[name]['status'],flush=True)
 results['soft']=soft_relaxation(results['vacuum']);print('soft',results['soft']['status'],flush=True)
 files=['data/w33_pass11461_11465_explicit_physical_maps.json','data/w33_pass11428_11432_native_alignment_measure_coarse_flags.json']
 results.update(status='PASS',passes=list(range(11466,11471)),reservation='90957fb96',source_sha256={p:hashlib.sha256(json.dumps(json.loads((ROOT/p).read_text()),sort_keys=True,separators=(',',':')).encode()).hexdigest() for p in files})
 OUT.write_text(json.dumps(results,indent=2)+'\n');return results
if __name__=='__main__':run()
