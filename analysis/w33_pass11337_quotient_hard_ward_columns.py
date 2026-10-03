"""Actual quotient family vector fields fix hard zero-momentum Ward columns."""
from pathlib import Path
import sys,json
import numpy as np
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
from w33_pass11327_nonlinear_quotient_mass_maps import Chart
import w33_pass11322_hard_soft_invariant_matching as H

def orbit_coordinates(C,F,theta):
 q=expm(1j*theta*F)@C.p;orbit=np.einsum('aij,j->ia',C.H,C.p);U=np.linalg.svd(orbit,full_matrices=False)[0][:,:70]
 for _ in range(10):
  r=U.conj().T@(q-C.p)
  if max(abs(r))<2e-13:break
  action=np.einsum('aij,j->ia',C.H,q);t=-np.linalg.pinv(U.conj().T@action,rcond=1e-10)@r;q=expm(np.einsum('a,aij->ij',t,C.H))@q
 else:raise RuntimeError('holomorphic slice failed')
 z=C.N.conj().T@(q-C.p);assert max(abs(q-C.p-C.N@z))<1e-11
 return np.sqrt(2)*np.r_[z.real,z.imag]
def hard_value(maps):
 m=.01;sw=np.linalg.eigvalsh(m*m*maps['scalar']);fw=np.linalg.eigvalsh(m*m*maps['Weyl'].conj().T@maps['Weyl']);vw=np.linalg.eigvalsh(maps['vector'])
 def f(w,c):
  w=w[w>1e-6];return np.sum(w*w*(np.log(w/.01)-c))
 return float((f(sw,1.5)-2*f(fw,1.5)+3*f(vw,5/6))/(64*np.pi**2))
def payload():
 print('chart',flush=True);C=Chart();g=np.array(json.loads((ROOT/'data/w33_pass11322_hard_soft_invariant_matching.json').read_text())['actual_unregulated_first_jet']);F=H.Q.Q.L.generators()[78:];orbit=np.einsum('aij,j->ia',C.H,C.p);T=[];J=[]
 for f in F:
  v=1j*f@C.p;a=np.linalg.pinv(orbit,rcond=1e-10)@v;L=1j*f-np.einsum('a,aij->ij',a,C.H);t=C.N.conj().T@v;j=C.N.conj().T@L@C.N;T.append(np.sqrt(2)*np.r_[t.real,t.imag]);J.append(np.block([[j.real,-j.imag],[j.imag,j.real]]))
 T=np.array(T).T;J=np.array(J);assert np.linalg.matrix_rank(T,tol=1e-9)==8
 K=np.column_stack([-j.T@g for j in J]);ward1=float(max(abs(T.T@g)));compat=float(max(abs(T.T@K-K.T@T).flat));assert ward1<1e-10 and compat<1e-8
 A=K@np.linalg.pinv(T,rcond=1e-10);P=T@np.linalg.pinv(T,rcond=1e-10);W=A+A.T-P@A;W=(W+W.T)/2;res=float(max(abs(W@T-K).flat));assert res<1e-8
 f=F[0];eps=.001;xp=orbit_coordinates(C,f,eps);xm=orbit_coordinates(C,f,-eps);t=T[:,0];j=J[0];tangent_error=float(np.linalg.norm((xp-xm)/(2*eps)-t));acceleration_error=float(np.linalg.norm((xp+xm)/(eps*eps)-j@t));assert tangent_error<1e-6 and acceleration_error<1e-5
 print('Ward/orbit',ward1,compat,res,tangent_error,acceleration_error,flush=True);maps=C.mass_maps(np.zeros(22));v0=hard_value(maps);print('reference hard',v0,flush=True);vp=hard_value(C.mass_maps(xp));vm=hard_value(C.mass_maps(xm));err=max(abs(vp-v0),abs(vm-v0));assert err<1e-8
 return {'status':'PASS','scope':'Explicit family-vector-field and hard first-jet Ward determination of eight zero-momentum columns, with actual nonlinear orbit/mass-map control; not complete hard1PI normal block or total resummation.','family_fields':'Holomorphic gauge compensation L=iF-Ha, a=(Hp)^+ iFp. Atp, t=Ndag iFp, DT=Ndag L N. Convert to22 real coordinates. Finite complex-E6 slice retraction independently checks t and DT.t.','hard_gradient':g.tolist(),'Goldstone_tangents':T.tolist(),'hard_Ward_columns':K.tolist(),'symmetric_minimal_completion':W.tolist(),'Ward_equation':'Hhard T=-DT^T grad(Vhard). Therefore hard Goldstone mass columns do not generally vanish before tadpole matching. An invariant matching functional cancelling hard first/second jets cancels these columns consistently.','first_Ward_residual':ward1,'column_symmetry_compatibility':compat,'completion_column_residual':res,'finite_orbit_tangent_error':tangent_error,'finite_orbit_acceleration_error':acceleration_error,'hard_potential_on_family_orbit':[v0,vp,vm],'hard_family_orbit_error':err,'boundary':['The displayed minimal symmetric completion fixes Goldstone columns and sets the free normal-normal block by convention; it is not that block calculated from diagrams.','The hard nonlocal cuts are prior11287/11295; full analytic normal matching, heavy-vector/gaugino UV thresholds and higher-loop Goldstone subtractions remain uncomputed.','Numerical chart and derivatives retain their finite-difference accuracy. Globalfamily remains global; no gauged-family Goldstone count is substituted.'],'prior_owners':['analysis/w33_pass11327_nonlinear_quotient_mass_maps.py','analysis/w33_pass11322_hard_soft_invariant_matching.py','analysis/w33_pass11295_ward_matched_modulus_poles.py'],'primary_sources':['https://arxiv.org/abs/hep-ph/9603341','https://arxiv.org/abs/1609.06977']}
if __name__=='__main__':
 d=payload();(ROOT/'data/w33_pass11337_quotient_hard_ward_columns.json').write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['hard_family_orbit_error'],d['completion_column_residual'])
