"""Actual native pair-vielbein boosts on a collinear ADM slice; no generic ghost claim."""
from pathlib import Path
import sys,json
import numpy as np
from scipy.optimize import root
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
from w33_pass11289_cycle_gram_gluing_flux import cycle_basis
from w33_pass11323_native_cycle_lapse_hessian import solve_currents

def solve_boosts(B,p):
 # Fix the common Lorentz gauge eta[0]=0. Convex graph cosh gives a unique root.
 D=B[1:];f=lambda z:D@np.sinh(D.T@z)-p[1:];jac=lambda z:(D*np.cosh(D.T@z))@D.T
 z=root(f,np.zeros(79),jac=jac,tol=1e-11);eta=np.r_[0.,z.x];j=np.sinh(B.T@eta);assert max(abs(B@j-p))<1e-10
 return eta,j

def payload():
 import sympy as s
 ni,nj,si,sj,t=s.symbols('ni nj si sj t',real=True);ei=s.Matrix([[ni,0],[si,1]]);ej=s.Matrix([[s.cosh(t),s.sinh(t)],[s.sinh(t),s.cosh(t)]])*s.Matrix([[nj,0],[sj,1]])
 pair=s.simplify((ni*s.trace(ei.inv()*ej)+nj*s.trace(ej.inv()*ei))/2);target=(ni+nj)*s.cosh(t)+(sj-si)*s.sinh(t);assert s.simplify(pair-target)==0
 B,C=cycle_basis();rng=np.random.default_rng(11338);p=rng.normal(size=80)*.15;p-=p.mean();N=1+rng.uniform(-.2,.2,80);eta,j=solve_boosts(B,p);c=np.cosh(B.T@eta);L=abs(B).T@N
 A=(B*c)@B.T;shift=-np.linalg.pinv(A,rcond=1e-11)@(B@(L*j));lorentz=B@(L*j+c*(B.T@shift));assert max(abs(lorentz))<1e-10
 def energy(N):return float((abs(B).T@N)@c)
 grad=abs(B)@c;h=1e-4;direction=rng.normal(size=80);fd=(energy(N+h*direction)+energy(N-h*direction)-2*energy(N))/(h*h);assert abs(fd)<1e-4
 jm,Hm,res=solve_currents(B,C,N,p);rank=int(np.linalg.matrix_rank(Hm,tol=1e-10));assert rank==78
 return {'status':'PASS','scope':'Exact nonlinear lapse-linearity on one-spatial-direction collinear boost slice of an actual native pair-vielbein action, not generic4D ghost freedom.','pair_action':'Half of det(ei)Tr(ei^-1 ej)+det(ej)Tr(ej^-1 ei), with two identity spectator directions, gives (Ni+Nj)cosh(eta_j-eta_i)+(sj-si)sinh(eta_j-eta_i)+(Ni+Nj). The omitted last term is lapse-linear. Add -p.s.','shift_constraint':'B sinh(B.T eta)=p. This convex-gradient graph equation fixes79 relative boosts independently of lapses. No separate edge boosts are introduced.','Lorentz_constraint':'B[(|B|.T N)*sinh(B.T eta)+cosh(B.T eta)*(B.T s)]=0 then fixes79 relative shifts linearly in N.','reduced_action':'Sum_e (Ni+Nj)cosh((B.T eta)_e); exactly linear in every lapse. Thus the full80x80 reduced lapse Hessian vanishes on this slice despite81 graph cycles.','boosts':eta.tolist(),'currents':j.tolist(),'shift_constraint_residual':float(max(abs(B@j-p))),'Lorentz_constraint_residual':float(max(abs(lorentz))),'cycle_boost_residual':float(max(abs(C.T@np.arcsinh(j)))),'lapse_gradient':grad.tolist(),'finite_difference_lapse_curvature':fd,'metric_slice_comparison':{'lapse_rank':rank,'cycle_boost_residual':float(max(abs(C.T@np.arcsinh(jm))))},'boundary':['The independent transverse boosts, spatial rotations and all generic spatial coframes remain untested. A slice with lapse-linear reduced action does not prove all secondary/tertiary constraints.','This is not a retraction of11323: integrating independent pair metric square roots and imposing one Lorentz frame per vertex define different cycle theories. The boost-cycle compatibility is an explicit difference.','No supplied spacetime, Planck scale or cosmological constant is derived.'],'prior_owners':['analysis/w33_pass11323_native_cycle_lapse_hessian.py','analysis/w33_pass11328_determinant_multivielbein_scope.py'],'primary_sources':['https://arxiv.org/abs/1410.7774','https://arxiv.org/html/2510.03014v2']}
if __name__=='__main__':
 d=payload();(ROOT/'data/w33_pass11338_native_boost_constraint_slice.json').write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['cycle_boost_residual'],d['metric_slice_comparison'])
