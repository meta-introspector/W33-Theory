import json,sys,hashlib
from pathlib import Path
from itertools import product
import numpy as np
import sympy as sp
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11466_11470_stationary_ports_channels as P

def packet():return json.loads(P.OUT.read_text())

def test_source_binding():
 for name,digest in packet()['source_sha256'].items():
  raw=(ROOT/name).read_text()
  for text in [raw,raw.replace('\n','\r\n')]:assert hashlib.sha256(json.dumps(json.loads(text),sort_keys=True,separators=(',',':')).encode()).hexdigest()==digest

def test_constructed_D4_in_another_basis():
 v=packet()['vacuum'];K=P.dec(v['generators']);E=P.dec(v['fixed_frame']);rng=np.random.default_rng(11466)
 U=np.linalg.qr(rng.normal(size=(27,27))+1j*rng.normal(size=(27,27)))[0];L=np.array([U.conj().T@k@U for k in K]);flat=L.reshape(28,-1).T
 for a in L[::4]:
  for b in L[::3]:
   c=(1j*(a@b-b@a)).ravel();assert np.linalg.norm(c-flat@np.linalg.lstsq(flat,c,rcond=None)[0])<1e-10
 assert np.linalg.norm(E.conj().T@E-np.eye(9))<1e-10
 assert max(np.linalg.norm(np.kron(k,np.eye(3))@E) for k in K)<1e-10
 roots=np.array(v['roots']);roots=roots[np.sum(roots**2,axis=1)>1e-9];assert len(roots)==24
 G=roots@roots.T;assert np.max(np.min(abs(G[:,:,None]-np.array([-1/3,-1/6,0,1/6,1/3])),axis=2))<1e-10
 assert v['threshold_ranks']==[58,58,58]
 assert len(v['normal_eigenvalues'])==266 and sum(x>1e-3 for x in v['normal_eigenvalues'])==262
 assert 'not symbolically' in v['scope']

def test_stationarity_and_gradient_against_scalar_difference():
 import w33_pass11438_finite_native_model as F
 v=packet()['vacuum'];x=np.array(v['coordinates']);d=np.random.default_rng(11466).normal(size=324);d/=np.linalg.norm(d)
 with F.native_context():
  value,g=F.fun(x);h=1e-5;fd=(F.fun(x+h*d)[0]-F.fun(x-h*d)[0])/(2*h)
 assert np.linalg.norm(g)<1e-6 and abs(fd-g@d)<1e-6 and abs(value-v['value'])<1e-10


def test_all_charges_and_topological_negative_control():
 M=P.native();p=packet()['flux'];u,v=M.N.P.N.uniform_flux_links(p['L']);rng=np.random.default_rng(11467);g=np.exp(1j*rng.normal(size=u.shape))
 u2=g*u*np.roll(g.conj(),-1,axis=0);v2=g*v*np.roll(g.conj(),-1,axis=1)
 for row in p['charges']:
  q=row['charge'];a=u2**q;b=v2**q;pl=a*np.roll(b,-1,axis=0)*np.roll(a.conj(),-1,axis=1)*b.conj()
  assert abs(np.angle(pl).sum()/2/np.pi-q)<1e-9 and np.max(abs(1-pl))<1/30
 assert max(r['maximum_plaquette'] for r in p['path'])>1.9
 assert p['path'][0]['flux']==0 and abs(p['path'][-1]['flux']-1)<1e-9


def test_four_interior_equations_and_conserved_matter_current():
 import mpmath as mp
 p=packet()['gravity'];s=p['interior_scale'];h0,h1=p['heights'];_,f,q=p['scalar_values']
 # Independent Schlaefli height equation, area derivative from equation5.
 for a,b,h,x,y in [(1,s,h0,0,f),(s,1.2,h1,f,q)]:
  d=b-a;root=np.sqrt(h*h-d*d/2);area_h=h*(a+b)/(2*root)
  equation=area_h*np.arcsin(d*d/(4*h*h-d*d))-np.pi*(a+b)**3*(y-x)**2/(24*h*h)
  assert abs(equation)<1e-10
 assert abs(p['scalar_current'][0]-p['scalar_current'][1])<1e-10
 assert max(abs(x) for x in p['finite_difference_forces'])<1e-7 and min(p['causal_margins'])>0
 assert min(p['timelike_trapezoid_deficits'])>.01
 with mp.workdps(35):
  def a(t):return P.frustum(mp.mpf(1),t,mp.mpf(h0),mp.mpf(0),mp.mpf(f))+P.frustum(t,mp.mpf('1.2'),mp.mpf(h1),mp.mpf(f),mp.mpf(q))
  assert abs(float(mp.diff(a,mp.mpf(s))))<1e-10
  assert abs(float(mp.diff(a,mp.mpf(s)+mp.mpf('.001'))))>1e-5


def test_exact_Krylov_restrictions_and_unseen_port_parameter():
 p=packet()['mediators'];M=P.native();A=np.array(M.N.P.N.M.P.load()['A'],int);Am=sp.Matrix(A.tolist());rng=np.random.default_rng(11469)
 for row in p['exact_bases']:
  B=sp.Matrix(row['integer_basis']);T=sp.Matrix([[sp.Rational(x) for x in r] for r in row['rational_restriction']]);assert Am*B==B*T and B.rank()==B.cols
  Q=np.linalg.qr(np.array(B,float))[0];rot=np.linalg.qr(rng.normal(size=(Q.shape[1],Q.shape[1])))[0];Q=Q@rot
  J=np.eye(80)[:,[0,row['endpoint']]];small=Q.T@J;phi=.437
  assert np.linalg.norm(J.T@np.linalg.solve(10*np.eye(80)-phi*A,J)-small.T@np.linalg.solve(10*np.eye(Q.shape[1])-phi*Q.T@A@Q,small))<1e-10
 assert p['two_sector_species']==[972,96] and p['two_sector_b0']==[-637.,-53.]
 assert p['discarded_species']==876


def test_discarded_determinant_cannot_be_silently_removed():
 p=packet()['mediators'];h=1e-5;phi=.81
 logdark=lambda t:340*np.log(10)+268*np.log(100-6*t*t)
 derivative=(logdark(phi+h)-logdark(phi-h))/(2*h)
 assert abs(derivative-p['discarded_logdet_force_at_081'])<1e-7 and abs(derivative)>20
 assert max(abs(r['logdet_factorization_error']) for r in p['scans'])<1e-9
 assert max(r['spectrum_error'] for r in p['scans'])<1e-10


def test_complete_channels_and_nonpauli_counterexample():
 p=packet()['relay']
 for row in p['rows']:
  for key in ['choi_zero','choi_one','choi_plus']:
   C=P.dec(row[key]);assert np.linalg.norm(C-C.conj().T)<1e-10 and min(np.linalg.eigvalsh(C))>-1e-10
   assert np.linalg.norm(np.einsum('iaja->ij',C.reshape(4,4,4,4))-np.eye(4)/4)<1e-10
 assert p['rows'][0]['zero_one_difference']<1e-12 and p['rows'][0]['zero_plus_difference']<1e-12
 assert p['rows'][1]['zero_one_difference']>.001
 # Independent ideal route, then arbitrary initially entangled relay/data.
 U=np.eye(8)
 for c,t in P.OPS:U=P.cnot(c,t)@U
 assert np.linalg.norm(U-P.cnot(0,2))==0
 rng=np.random.default_rng(11470);psi=rng.normal(size=8)+1j*rng.normal(size=8);psi/=np.linalg.norm(psi);rho=np.outer(psi,psi.conj())
 partial=lambda r:np.einsum('arbcrd->abcd',r.reshape(2,2,2,2,2,2)).reshape(4,4)
 base=partial(U@rho@U.T)
 for a,b,c in product(range(4),repeat=3):
  E=np.kron(np.kron(P.PAULI[a],P.PAULI[b]),P.PAULI[c]);got=partial(E@U@rho@U.T@E.conj().T);D=np.kron(P.PAULI[a],P.PAULI[c]);assert np.linalg.norm(got-D@base@D.conj().T)<1e-12


def test_soft_relaxation_replays_actual_energy_and_soft_constraints():
 import w33_pass11438_finite_native_model as F
 p=packet();soft=p['soft'];B=P.dec(soft['soft_frame']);x=np.array(p['vacuum']['coordinates'])
 assert soft['reduced_gauge_rank']==10 and len(soft['reduced_normal_eigenvalues'])==26
 assert len(soft['rows'])==24 and np.linalg.norm(B.T@B-np.eye(4))<1e-10
 with F.native_context():
  base=F.fun(x)[0]
  for row in soft['rows']:
   xx=np.array(row['coordinates']);v,g=F.fun(xx);expected=row['sign']*row['amplitude']*np.array(row['coefficients'])
   assert np.linalg.norm(B.T@(xx-x)-expected)<1e-9
   assert abs(v-base-row['energy_change'])<1e-10
   assert abs(np.linalg.norm(B.T@g)-row['soft_gradient_norm'])<1e-9
 assert 'No uniform quartic positivity' in soft['scope']
