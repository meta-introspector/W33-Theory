"""CP-even angular flavor operator: noncommuting full24-field local vacuum, not observed CKM."""
from pathlib import Path
import sys,json
import numpy as np
from scipy.linalg import null_space,expm
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11321_full_flavor_bifurcation as F
J_NORMALIZATION=143327232

def invariants(x):
 U,D=F.fields(x);A=U@U.conj().T;B=D@D.conj().T;R=np.trace(A+B).real;C=A@B-B@A;J=np.trace(C@C@C).imag
 return U,D,A,B,R,C,J

def val(x,eta):
 U,D,A,B,R,C,J=invariants(x)
 return F.val(x,.39)-eta*J*J/(.1+R)**12

def gradient(x,eta,h=1e-5):return np.array([(val(x+e,eta)-val(x-e,eta))/(2*h) for e in np.eye(24)*h])
def hessian(x,eta,h=2e-4):
 H=np.column_stack([(gradient(x+e,eta)-gradient(x-e,eta))/(2*h) for e in np.eye(24)*h]);return (H+H.T)/2

def audit(row):
 eta=float(row[0]);x=np.array(row[2:]);U,D,A,B,R,C,J=invariants(x);orbit=np.column_stack([F.pack(1j*(T@U+U@T.T),1j*(T@D+D@T.T)) for T in F.generators()]);rank=np.linalg.matrix_rank(orbit,tol=1e-5);N=null_space(orbit.T,rcond=1e-5);H=hessian(x,eta);H2=hessian(x,eta,1e-4);ev=np.linalg.eigvalsh(N.T@H@N);assert min(ev)>.01 and np.linalg.norm(gradient(x,eta))<1e-5
 a,VA=np.linalg.eigh(A);b,VB=np.linalg.eigh(B);V=VA.conj().T@VB;delta=lambda v:(v[1]-v[0])*(v[2]-v[0])*(v[2]-v[1]);jckm=np.imag(V[0,0]*V[1,1]*V[0,1].conj()*V[1,0].conj());assert abs(abs(J)-abs(6*delta(a)*delta(b)*jckm))<1e-10
 Rg=expm(1j*(.23*F.generators()[0]-.14*F.generators()[5]));xr=F.pack(Rg@U@Rg.T,Rg@D@Rg.T);cp=F.pack(U.conj(),D.conj());cov=abs(val(xr,eta)-val(x,eta));cperror=abs(val(cp,eta)-val(x,eta));assert max(cov,cperror)<1e-10
 return {'eta_coefficient':eta,'normalized_angular_strength':eta/J_NORMALIZATION,'coordinates':x.tolist(),'potential':val(x,eta),'gradient_norm':float(np.linalg.norm(gradient(x,eta))),'normal_Hessian_eigenvalues':ev.tolist(),'Hessian_step_control_error':float(np.max(abs(H-H2))),'family_orbit_rank':int(rank),'U_Gram_eigenvalues':a.tolist(),'D_Gram_eigenvalues':b.tolist(),'commutator_norm':float(np.linalg.norm(C)),'Im_trace_commutator_cubed':float(J),'mixing_moduli_squared':(abs(V)**2).tolist(),'CKM_type_Jarlskog':float(jckm),'covariance_error':cov,'CP_conjugate_energy_error':cperror,'CP_conjugate_J':float(invariants(cp)[-1])}

def reproduce_discovery():
 from scipy.optimize import minimize
 rng=np.random.default_rng(113260);rows=[]
 for eta in [1e9,1e10]:
  best=None
  for _ in range(12):
   z=minimize(lambda x:val(x,eta),rng.normal(size=24)*.3,method='BFGS',options={'gtol':1e-7,'maxiter':600})
   if best is None or z.fun<best.fun:best=z
  rows.append([eta,float(best.fun),*best.x.tolist()])
 return rows

WITNESSES=[[10000000.0, 3.247477655756729, -0.3214648539179008, -0.11339520679329826, 0.026485699960930663, -0.19688866211909922, -0.0743588296822027, 0.18600863927867267, 0.08021772346576536, 0.17012831855679272, 0.16245203215315768, -0.16612498461188593, -0.03757149897126624, 0.28881984654033593, -0.2615113828410017, 0.43081270890186724, -0.9497716715721071, -0.43538263889232043, -0.6502006812275948, 0.08147348368887536, 0.0027605152639468434, -0.3683560988433444, -0.4070394867142176, 0.42196931618700984, -0.11200417588002637, -0.04274628011225897], [1000000000.0, -1.3053317729148048, -0.13865376117162237, -0.14191415903744137, -0.09736409411792124, -0.44236122009864437, 0.14029033070507343, 0.12965419461720462, -0.41772531555845166, 0.978399103791248, 0.47936222665045414, -0.14192862020083158, -0.2548074724318767, -0.12651621765087134, 0.4097330041064832, 0.6118150823037374, 0.13408898757668525, -0.1596622117748167, 0.12357466031483812, -0.37380167810585047, -0.33443385361926165, -0.2083428442854047, 0.6026247641606349, 0.4762576108223942, 0.5024283371828637, -0.09798271230811645]]

def payload():
 rows=[audit(r) for r in WITNESSES];strong=rows[1];assert strong['family_orbit_rank']==8 and abs(strong['CKM_type_Jarlskog'])>.09
 return {'status':'PASS','result_scope':'PASS_ENGINEERED_CP_EVEN_FULL_FIELD_NONCOMMUTING_LOCAL_FLAVOR_VACUUM_NOT_OBSERVED_MIXING','base':'Actual11321 full24-field flavor potential and two singlet probe loops at sigma.39.','added_operator':'-eta [ImTr([UUdag,DDdag]³)]²/(.1+Tr(UUdag+DDdag))^12. CP-even, familySU3 invariant, smooth at the origin, with supplied cutoff. No observed mass, angle, phase or CP root is used. This explicitly rewards noncommuting chirality rather than deriving the operator from W33.','bounded_angular_operator':'For positive3x3 A,B, |ImTr[A,B]³|=6|Delta(A)Delta(B)Jmix|. At fixed traceq, |Delta(A)|<=q³/(6sqrt3), and |Jmix|<=1/(6sqrt3). Thus J²<=R^12/143327232, R=q+r. The added rational term is bounded below by-eta/143327232 and cannot spoil tree coercivity. Fermion CW still needs a high-field UV completion.','vacua':rows,'interpretation':'The eta1e9 local vacuum has all three distinct singular values in each field, nonzero Gram commutator, nonzero CKM-type invariant and positive16 normal modes after removing the8 family orbit directions. Its conjugate reverses CP invariants at equal energy. The small light eigenvalues are numerical local outputs of chosen couplings, not observed masses.','scope':['This engineered high-order/rational EFT operator proves a structural escape from the previous commuting2+1 branch, not a prediction or a renormalizable weak-coupling UV theory.','Mixing moduli are nearly1/3 and |Jmix| approximately.0962, far from observed quark mixing. Do not fit or label them physical CKM values.','The normalized angular strength is about6.98 for the useful witness; scalar radiative control of this added sector and origin of its coefficient remain open.','Local finite-difference Hessians and multistart searches are numerical witnesses, not global or interval-certified vacua.'],'prior_owners':['analysis/w33_pass11321_full_flavor_bifurcation.py','analysis/w33_pass11304_singlet_phase_probes.py'],'primary_sources':['https://doi.org/10.1103/PhysRevLett.55.1039']}
if __name__=='__main__':
 if '--search' in sys.argv:print(json.dumps(reproduce_discovery(),indent=2));raise SystemExit
 d=payload();(ROOT/'data/w33_pass11326_noncommuting_flavor_operator.json').write_text(json.dumps(d,indent=2)+'\n');print(d['status'],[(r['eta_coefficient'],r['CKM_type_Jarlskog']) for r in d['vacua']])
