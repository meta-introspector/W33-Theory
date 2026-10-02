"""Numerical local D-flat quotient chart, metric, and covariant scalar/Weyl mass maps."""
from pathlib import Path
import sys,json
import numpy as np
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11287_quotient_self_energy_matrix as Q

class Chart:
 def __init__(self,p=None,N=None):
  p0,N0,self.H,*_=Q.geometry();self.p=p0 if p is None else p;self.N=N0 if N is None else N;self.S=np.column_stack([self.N,1j*self.N])/np.sqrt(2)
 def retract(self,x):
  q=self.p+self.S@x
  for _ in range(10):
   act=np.einsum('aij,j->ai',self.H,q);mu=np.real(np.einsum('i,ai->a',q.conj(),act))
   if max(abs(mu))<2e-13:return q
   R=(act.conj()@act.T).real;t=-.5*np.linalg.pinv(R,rcond=1e-10)@mu;q=expm(np.einsum('a,aij->ij',t,self.H))@q
  raise RuntimeError('D-flat projection failed')
 def metric(self,x,h=2e-5):
  q=self.retract(x);J=np.column_stack([(self.retract(x+e)-self.retract(x-e))/(2*h) for e in np.eye(22)*h]);act=np.einsum('aij,j->ia',self.H,q);U=np.linalg.svd(act,full_matrices=False)[0][:,:70];Jh=J-U@(U.conj().T@J);G=2*np.real(Jh.conj().T@Jh)
  assert np.linalg.matrix_rank(act,tol=1e-7)==70 and np.linalg.eigvalsh(G)[0]>.9
  return G
 def potential(self,x):return Q.Q.potential(self.retract(x))
 def vector_gram(self,x,g=.5):
  hp=g*np.einsum('aij,j->ai',self.H,self.retract(x));return 2*np.real(hp.conj()@hp.T)
 def mass_maps(self,x,h=.0004,gh=.0002):
  E0=np.eye(22);G=self.metric(x);dG=np.array([(self.metric(x+gh*e)-self.metric(x-gh*e))/(2*gh) for e in E0]);gi=np.linalg.inv(G)
  Gamma=np.einsum('kl,ilj->kij',gi,dG)/2+np.einsum('kl,jli->kij',gi,dG)/2-np.einsum('kl,lij->kij',gi,dG)/2
  v0=self.potential(x);grad=np.array([(self.potential(x+h*e)-self.potential(x-h*e))/(2*h) for e in E0]);Hess=np.zeros((22,22))
  for i in range(22):
   ei=h*E0[i];Hess[i,i]=(self.potential(x+ei)+self.potential(x-ei)-2*v0)/(h*h)
   for j in range(i):
    ej=h*E0[j];Hess[i,j]=Hess[j,i]=(self.potential(x+ei+ej)+self.potential(x-ei-ej)-self.potential(x+ei-ej)-self.potential(x-ei+ej))/(4*h*h)
  cov=Hess-np.einsum('kij,k->ij',Gamma,grad);w,O=np.linalg.eigh(G);E=(O/np.sqrt(w))@O.T;Xs=E.T@cov@E
  # The quotient class has holomorphic coordinates z=(x_real+i*x_imag)/sqrt2.
  # Gauge-invariant W can be differentiated on the holomorphic affine representative.
  phi=self.p+self.S@x;W,ambient=Q.Q.superpotential(phi);w1=self.N.T@ambient
  W2=np.column_stack([(self.N.T@Q.Q.superpotential(phi+h*self.N[:,j])[1]-self.N.T@Q.Q.superpotential(phi-h*self.N[:,j])[1])/(2*h) for j in range(11)]);W2=(W2+W2.T)/2
  K=G[:11,:11]-1j*G[:11,11:];dK=dG[:,:11,:11]-1j*dG[:,:11,11:];dz=(dK[:11]-1j*dK[11:])/np.sqrt(2);cg=np.einsum('kl,ilj->kij',np.linalg.inv(K),dz);Mw=W2-np.einsum('kij,k->ij',cg,w1);Mw=(Mw+Mw.T)/2
  wk,Uk=np.linalg.eigh(K);Ek=(Uk/np.sqrt(wk))@Uk.conj().T;Mw=Ek.T@Mw@Ek
  return {'scalar':Xs,'Weyl':Mw,'metric':G,'Christoffel':Gamma,'gradient':grad,'vector':self.vector_gram(x),'connection_Hessian_correction_norm':float(np.linalg.norm(np.einsum('kij,k->ij',Gamma,grad))),'complex_metric_Hermiticity_error':float(np.max(abs(K-K.conj().T))),'Kahler_block_error':float(max(np.max(abs(G[:11,:11]-G[11:,11:])),np.max(abs(G[:11,11:]+G[:11,11:].T))))}

def payload():
 C=Chart();x0=np.zeros(22);maps=C.mass_maps(x0);expected=json.loads((ROOT/'data/w33_pass11287_quotient_self_energy_matrix.json').read_text());H=np.array(expected['Hessian']);M=np.array(expected['Weyl_mass_real'])+1j*np.array(expected['Weyl_mass_imag']);serror=float(np.max(abs(np.linalg.eigvalsh(maps['scalar'])-np.linalg.eigvalsh(H))));ferror=float(np.max(abs(np.linalg.svd(maps['Weyl'],compute_uv=False)-np.linalg.svd(M,compute_uv=False))));assert serror<.005 and ferror<.001
 rng=np.random.default_rng(11327);x=rng.normal(size=22)*.001;q=C.retract(x);act=np.einsum('aij,j->ai',C.H,q);mu=np.real(np.einsum('i,ai->a',q.conj(),act));G=C.metric(x);generic_maps=C.mass_maps(x)
 cubic=np.array(expected["cubic_tensor"]);Y=np.array(expected["Yukawa_real"])+1j*np.array(expected["Yukawa_imag"]);jet_rows=[]
 for scale,actual in [(1.,generic_maps),(.5,C.mass_maps(.5*x))]:
  dx=scale*x;serr=float(np.linalg.norm(actual["scalar"]-H-np.einsum("a,aij->ij",dx,cubic)));ferr=float(np.linalg.norm(actual["Weyl"]-M-np.einsum("a,aij->ij",dx,Y)));jet_rows.append({"scale":scale,"scalar_linear_jet_residual":serr,"Weyl_linear_jet_residual":ferr})
 assert jet_rows[1]["scalar_linear_jet_residual"]<.4*jet_rows[0]["scalar_linear_jet_residual"] and jet_rows[1]["Weyl_linear_jet_residual"]<.4*jet_rows[0]["Weyl_linear_jet_residual"]
 # Change the ambient family frame, not just a mass eigenbasis.
 F=Q.Q.L.generators()[78:];U=expm(1j*(.19*F[0]-.11*F[3]+.07*F[7]));D=Chart(U@C.p,U@C.N);qd=D.retract(x);cov_error=float(np.max(abs(qd-U@q)));metric_error=float(np.max(abs(D.metric(x)-G)));potential_error=abs(D.potential(x)-C.potential(x));vector_error=float(np.max(abs(np.linalg.eigvalsh(D.vector_gram(x))-np.linalg.eigvalsh(C.vector_gram(x)))));assert max(cov_error,metric_error,potential_error,vector_error)<1e-8
 return {'status':'PASS','result_scope':'PASS_NUMERICAL_NONLINEAR_LOCAL_QUOTIENT_METRIC_AND_COVARIANT_MASS_MAPS_NOT_FULL_RESUMMATION','chart':'phi0=p+N z; iterate exp(t_a H_a)phi0 to solve all actual E6 moment maps. This extends prior11287 directional retraction to all22 coordinates. Rank70 and quotient dimension11 are checked locally.','metric':'G_ab=2Re(J_a† P_horizontal J_b), J=dphi_Dflat/dx, P=I-Uorbit Uorbit†. Christoffel from metric derivatives; scalar map is G^-1/2 [d²V-Gamma.dV] G^-1/2.','Weyl_map':'Wcov_ij=W_ij-GammaK^k_ij W_k on holomorphic quotient coordinates; normalize with the inverse square root of the Hermitian Kahler metric. Vector map is the actual78-generator Gram on the D-flat representative.','reference_scalar_spectral_error':serror,'reference_Weyl_singular_value_error':ferror,'reference_metric':maps['metric'].tolist(),'reference_scalar_mass_matrix':maps['scalar'].tolist(),'reference_Weyl_real':maps['Weyl'].real.tolist(),'reference_Weyl_imag':maps['Weyl'].imag.tolist(),'reference_connection_norm':float(np.linalg.norm(maps['Christoffel'])),'prior_22_field_linear_jet_controls':jet_rows,'generic_coordinates':x.tolist(),'generic_moment_residual':float(max(abs(mu))),'generic_metric_eigenvalues':np.linalg.eigvalsh(G).tolist(),'generic_covariant_scalar_eigenvalues':np.linalg.eigvalsh(generic_maps['scalar']).tolist(),'generic_Weyl_singular_values':np.linalg.svd(generic_maps['Weyl'],compute_uv=False).tolist(),'generic_connection_norm':float(np.linalg.norm(generic_maps['Christoffel'])),'generic_connection_Hessian_correction_norm':generic_maps['connection_Hessian_correction_norm'],'generic_Kahler_block_error':generic_maps['Kahler_block_error'],'family_frame_covariance':{'Dflat_map_error':cov_error,'metric_error':metric_error,'potential_error':potential_error,'vector_spectrum_error':vector_error},'hard_matching_bridge':'The nonlinear mass-map functions named by11322 now have a concrete numerical local implementation. Covariant scalar and Weyl masses use the quotient connection rather than arbitrary frozen affine Hessians. The prior hard spectral counterfunction theorem can be applied on a gap-preserving chart with its numerical accuracy stated.','resummation_boundary':['This builds the previously missing local maps; it does not compute full hard1PI Ward kernels, higher-loop double-counting subtractions or total observed poles.','The prior11272 radial Ward-resummed result remains a separate established restriction. Its radial conclusion is not extended to all22 fields without the hard kernel.','Chart and derivatives are finite-difference numerical witnesses, not directed intervals, global atlas or complete E6 UV threshold matching.'],'prior_owners':['analysis/w33_pass11270_full_condensate_moduli.py','analysis/w33_pass11287_quotient_self_energy_matrix.py','analysis/w33_pass11322_hard_soft_invariant_matching.py','analysis/w33_pass11272_resummed_condensate_radial_control.py'],'primary_sources':['https://arxiv.org/abs/1609.06977']}
if __name__=='__main__':
 d=payload();(ROOT/'data/w33_pass11327_nonlinear_quotient_mass_maps.json').write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['reference_scalar_spectral_error'],d['reference_Weyl_singular_value_error'])
