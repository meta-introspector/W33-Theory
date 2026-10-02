"""Full24-real-field invariant flavor bifurcation and spontaneous relative CP."""
from pathlib import Path
import sys,json
from scipy.linalg import expm, null_space
ROOT=Path(__file__).resolve().parents[1]
import numpy as np
from scipy.optimize import minimize
basis=[]
for i in range(3):
 for j in range(i,3):
  x=np.zeros((3,3));x[i,j]=x[j,i]=1 if i==j else 1/np.sqrt(2);basis.append(x)
basis=np.array(basis)
def fields(x):
 U=np.einsum('a,aij->ij',x[:6]+1j*x[6:12],basis);D=np.einsum('a,aij->ij',x[12:18]+1j*x[18:],basis);return U,D

def val(x,sig=.2):
 U,D=fields(x);A=U@U.conj().T;B=D@D.conj().T;q=np.trace(A).real;r=np.trace(B).real;Z=U.conj().T@D
 t=-q-1.05*r+.18*(q+r)**2+.07*(np.trace(A@A)+np.trace(B@B)).real+.4*np.trace(A@B).real+.08*abs(np.trace(Z))**2+sig*np.trace(Z@Z).real
 def ad(a):return np.column_stack([np.cross(a[1],a[2]),np.cross(a[2],a[0]),np.cross(a[0],a[1])])
 t+=(-.06*np.linalg.det(U)-.045*np.linalg.det(D)+.035*np.trace(ad(U)@D)+.027*np.trace(ad(D)@U)).real
 for k in [1.,2.]:
  M=np.block([[.2*(U+k*D).conj(),3*np.eye(3)],[3*np.eye(3),np.zeros((3,3))]]);w=np.linalg.eigvalsh(M.conj().T@M);t-=2*np.sum(w*w*(np.log(w/9)-1.5))/(64*np.pi**2)
 return float(t)

WITNESSES=[[0.3, 3.252908321592342, 0.0005089871508377615, -0.0006760246975690599, -0.0074492408983739536, -0.019012654909658457, -0.006221624514609614, 0.004690892594114998, -0.010189554783976667, 0.005222713680277627, -0.027282552255090124, -0.010842765763493934, -0.003346800403769369, 0.008425817893100495, -0.02084197903206308, 0.027678649392724693, 0.3050340090065354, 0.7785521403493184, 0.2547748752559109, -0.19209057543007824, 0.41726381452424194, -0.21386638480544143, 1.1171871269039486, 0.44399254326339965, 0.13705381087850224, -0.34503639378122], [0.37, 3.2517785893945, 0.09122764879932714, -0.018414976887751573, 0.135642512318932, 0.18647994444564897, 0.13388774627623426, -0.03466502622839455, 0.005090994783236085, -0.24168286748924897, 0.04622878825591571, -0.05757401333963392, 0.11137897142514419, -0.14295970350251602, -0.42306416679565745, -0.03894309770170429, -0.007047997036917913, -0.6312758923101476, 0.5682562063995358, 0.6814569659311915, 0.7899095693569516, -0.13423613394704598, 0.04294973032662418, -0.4321420876097745, 0.34045272903920976, 0.3590161079159839], [0.39, 3.24747765575671, -0.1727285341509043, -0.028913014048361493, -0.1170644454039104, -0.06887554556042362, -0.031469990975872435, -0.08696265555046348, 0.102989333183986, -0.2851609544226903, 0.37086806083383395, -0.1992309570069998, -0.1837200778235536, -0.09697274002772864, 0.5425979993211476, -0.05862779335692185, -0.23807122466232344, 0.496303173651864, 0.028722115473080536, 0.763219598917467, 0.5071752304960091, 0.45875187711754783, 0.4002391087196791, -0.5815396775921094, 0.34032612881286844, 0.0525451813846653]]

def pack(U,D):
 return np.r_[np.einsum('aij,ij->a',basis,U).real,np.einsum('aij,ij->a',basis,U).imag,np.einsum('aij,ij->a',basis,D).real,np.einsum('aij,ij->a',basis,D).imag]
def gradient(x,sig,h=2e-5):
 return np.array([(val(x+e,sig)-val(x-e,sig))/(2*h) for e in np.eye(24)*h])
def hessian(x,sig,h=4e-4):
 H=np.column_stack([(gradient(x+e,sig)-gradient(x-e,sig))/(2*h) for e in np.eye(24)*h]);return (H+H.T)/2

def generators():
 out=[]
 for i in range(3):
  for j in range(i+1,3):
   T=np.zeros((3,3),complex);T[i,j]=T[j,i]=1;out.append(T/2)
   T=np.zeros((3,3),complex);T[i,j]=-1j;T[j,i]=1j;out.append(T/2)
 out.extend([np.diag([1,-1,0])/2,np.diag([1,1,-2])/(2*np.sqrt(3))]);return out

def audit(sig,x):
 U,D=fields(x);A=U@U.conj().T;B=D@D.conj().T
 orbit=np.column_stack([pack(1j*(T@U+U@T.T),1j*(T@D+D@T.T)) for T in generators()]);rank=np.linalg.matrix_rank(orbit,tol=1e-5);normal=null_space(orbit.T,rcond=1e-5)
 H=hessian(x,sig);H2=hessian(x,sig,2e-4);ev=np.linalg.eigvalsh(normal.T@H@normal);err=np.max(abs(H-H2));assert min(ev)>1e-3 and np.linalg.norm(gradient(x,sig))<2e-6
 assert rank==(5 if sig<.36 else 7)
 R=expm(1j*(.21*generators()[0]-.13*generators()[4]+.19*generators()[-1]));rot=pack(R@U@R.T,R@D@R.T)
 inv_error=abs(val(rot,sig)-val(x,sig));cp=pack(U.conj(),D.conj());cp_error=abs(val(cp,sig)-val(x,sig));assert max(inv_error,cp_error)<1e-10
 return {'sigma':sig,'coordinates':x.tolist(),'potential':val(x,sig),'gradient_norm':float(np.linalg.norm(gradient(x,sig))),'Hessian':H.tolist(),'Hessian_step_control_max_difference':float(err),'orbit_rank':int(rank),'normal_Hessian_eigenvalues':ev.tolist(),'Goldstone_residual':float(np.linalg.norm(H@orbit)),'U_Gram_eigenvalues':np.linalg.eigvalsh(A).tolist(),'D_Gram_eigenvalues':np.linalg.eigvalsh(B).tolist(),'Im_trace_UdaggerD':float(np.trace(U.conj().T@D).imag),'Im_det_U':float(np.linalg.det(U).imag),'Im_det_D':float(np.linalg.det(D).imag),'commutator_norm':float(np.linalg.norm(A@B-B@A)),'Jarlskog_commutator_cubic':float(np.trace(np.linalg.matrix_power(A@B-B@A,3)).imag),'SU3_potential_covariance_error':inv_error,'CP_conjugate_potential_error':cp_error}

def reproduce_discovery(starts=12):
 """Re-run the original two independent seed11321 grids; no witness initialization."""
 rows=[]
 for grid in [[.2,.3,.35],[.37,.39,.399]]:
  rng=np.random.default_rng(11321)
  for sig in grid:
   best=None
   for _ in range(starts):
    z=minimize(lambda x:val(x,sig),rng.normal(size=24)*.3,method='BFGS',options={'gtol':1e-7,'maxiter':400})
    if best is None or z.fun<best.fun:best=z
   rows.append([sig,float(best.fun),*best.x.tolist()])
 return rows

def payload():
 rows=[audit(float(row[0]),np.array(row[2:])) for row in WITNESSES]
 assert abs(rows[0]['Im_trace_UdaggerD'])<1e-5 and all(abs(r['Im_trace_UdaggerD'])>.01 for r in rows[1:])
 def symfield(z):return pack((z[0]+1j*z[1])*np.eye(3),(z[2]+1j*z[3])*np.eye(3))
 sym=minimize(lambda z:val(symfield(z),.37),[-.0215,0,.932,0],method='BFGS',options={'gtol':1e-8});symx=symfield(sym.x);syev=np.linalg.eigvalsh(hessian(symx,.37));assert sum(syev<-1e-3)>=5
 return {'status':'PASS','result_scope':'PASS_FULL_FIELD_LOCAL_FLAVOR_BIFURCATION_AND_RELATIVE_CP_NOT_CKM','fields':'Two complex symmetric3x3 fields U,D transform by congruence R U R^T,R D R^T for familySU3. All24 real Frobenius coordinates fluctuate. All coefficients are real.','tree_potential':'-q-1.05r+.18(q+r)^2+.07(TrA²+TrB²)+.4TrAB+.08|Tr(UdagD)|²+sigma ReTr[(UdagD)²]+Re[-.06detU-.045detD+.035Tr(adjU D)+.027Tr(adjD U)], A=UUdag,B=DDdag,q=TrA,r=TrB.','probe_loops':'The actual11304 two vectorlike Weyl singlet probes, kappa1,2, mass3 and Yukawa.2; CW=-2 sum masses^4[ln(masses²/9)-3/2]/(64pi²). No chosen singular values, angles or CP roots.','tree_boundedness':'|Tr[(UdagD)²]|<=TrAB. For |sigma|<.4 the mixed quartics are nonnegative; .18(q+r)^2 remains strictly coercive. Cubics do not change tree boundedness. Fermion CW alone is not a global high-field completion.','discovery':'12 deterministic-seed11321 random24-field BFGS starts for each sampled sigma. These stored witnesses replay stationarity, finite-difference step controls, symmetry covariance and every normal Hessian mode; no global-minimum or certified interval claim.','vacua':rows,'symmetric_branch_control_sigma037':{'coordinates':symx.tolist(),'full_Hessian_eigenvalues':syev.tolist(),'negative_modes_below_minus001':int(sum(syev<-1e-3))},'interpretation':'The symmetric CP-conserving SU3/SO3 branch develops traceless instabilities. Sampled stable sigma.37,.39 branches split singular values2+1, break familySU3 to a one-dimensional stabilizer and have nonzero CP-odd invariants ImTr(UdagD). Their CP conjugates have equal energy. This is spontaneous relative CP with real input couplings.','scope':['The Gram matrices remain commuting with a degenerate doublet. The Jarlskog commutator invariant vanishes; no observed CKM mixing or full three-generation hierarchy is obtained.','SU3 is a family symmetry of this supplied EFT; if gauged, its Goldstones are eaten and gauge loops must be included.','Positive normal Hessians are numerical local witnesses, not full UV stability, measured masses or a cosmological-constant prediction.','The transition is bracketed by sampled points, not a determined critical coupling or exhaustive phase diagram.'],'prior_owners':['analysis/w33_pass11304_singlet_phase_probes.py','analysis/w33_pass11314_untargeted_phase_vacuum.py']}
if __name__=='__main__':
 if '--search' in sys.argv:
  print(json.dumps(reproduce_discovery(),indent=2));raise SystemExit
 d=payload();(ROOT/'data/w33_pass11321_full_flavor_bifurcation.json').write_text(json.dumps(d,indent=2)+'\n');print(d['status'],[(r['sigma'],r['orbit_rank'],r['Im_trace_UdaggerD']) for r in d['vacua']])
