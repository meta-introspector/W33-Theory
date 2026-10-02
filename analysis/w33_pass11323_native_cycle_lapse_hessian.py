"""Nonlinear ADM shift elimination on the actual80-site160-edge metric theory graph."""
from pathlib import Path
import sys,json
import numpy as np
from scipy.optimize import root
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
from w33_pass11289_cycle_gram_gluing_flux import cycle_basis

def solve_currents(B,C,N,p):
 U=abs(B);L=U.T@N;assert min(L)>0
 j0=B.T@np.linalg.pinv(B@B.T)@p
 if C.shape[1]==0:return j0,np.zeros((len(N),len(N))),0.
 fun=lambda c:C.T@(L*(j0+C@c)/np.sqrt(1+(j0+C@c)**2))
 jac=lambda c:C.T@np.diag(L/(1+(j0+C@c)**2)**1.5)@C
 q=root(fun,np.zeros(C.shape[1]),jac=jac,tol=1e-11);res=float(max(abs(fun(q.x))));assert res<1e-10
 j=j0+C@q.x;G=U@((j/np.sqrt(1+j*j))[:,None]*C);H=-G@np.linalg.solve(jac(q.x),G.T)
 return j,H,res

def flag_orbit():
 from w33_pass11278_symplectic_ring_geometry import points,J
 from w33_pass11279_history_flux_obstruction import complex_data
 P=points(3);_,_,lines,_=complex_data();lookup={tuple(v):i for i,v in enumerate(P)};llook={tuple(sorted(t)):i for i,t in enumerate(lines)}
 flags=[(p,i) for i,t in enumerate(lines) for p in t];flook={x:i for i,x in enumerate(flags)};perms=[]
 for v in P:
  T=(np.eye(4,dtype=int)+np.outer(v,v@J))%3;assert not np.any((T.T@J@T-J)%3)
  perm=[]
  for w in P:
   z=T@w%3;ix=next(i for i,x in enumerate(z) if x);z=z*pow(int(z[ix]),-1,3)%3;perm.append(lookup[tuple(z)])
  lp=[llook[tuple(sorted(perm[x] for x in t))] for t in lines];perms.append([flook[(perm[p],lp[i])] for p,i in flags])
 orbit={0};queue=[0]
 for x in queue:
  for perm in perms:
   y=perm[x]
   if y not in orbit:orbit.add(y);queue.append(y)
 assert len(orbit)==160
 return {'symplectic_transvections':40,'native_flag_orbit_size':160,'invariant_tree_obstruction':'These native symplectic automorphisms are edge-transitive. An invariant edge subset is empty or all160 edges; no79-edge spanning tree preserves this full subgroup.'}

def payload():
 B,C=cycle_basis();rng=np.random.default_rng(11323);p=rng.normal(size=80)*.15;p-=np.mean(p);N=1+rng.uniform(-.2,.2,80);j,H,res=solve_currents(B,C,N,p);ev=np.linalg.eigvalsh(H);rank=int(sum(abs(ev)>1e-10));assert rank==78 and max(ev)<1e-10
 adj=[[] for _ in range(80)]
 for e in range(160):
  i,k=np.where(B[:,e])[0];adj[i].append((k,e));adj[k].append((i,e))
 seen={0};queue=[0];tree=[]
 for i in queue:
  for k,e in adj[i]:
   if k not in seen:seen.add(k);queue.append(k);tree.append(e)
 Bt=B[:,tree];jt,Ht,_=solve_currents(Bt,np.zeros((79,0)),N,p);assert Ht.shape==(80,80) and len(tree)==79
 shift=(abs(B).T@N)*j/np.sqrt(1+j*j);assert max(abs(C.T@shift))<1e-10
 checker=np.r_[np.ones(40),-np.ones(40)];assert max(abs(H@N))<1e-10 and max(abs(H@checker))<1e-10
 return {'status':'PASS','result_scope':'PASS_NONLINEAR_NATIVE_METRIC_CYCLE_CONSTRAINT_FAILURE_WITNESS','metric_slice':'Each site has spatial metricI3, lapseN_i>0 and shift s_i in one physical spatial direction. The beta1 dRGT pair potential is w_e[sqrt((N_i+N_j)²-(s_i-s_j)²)+2N_i], principal real square-root branch. Add the ADM momentum source sum p_i s_i. Linear lapse terms are omitted from the Hessian only.','reduced_potential':'Let B be oriented incidence, C an81-column native integral cycle basis, j=j0+C c with B j=p. Ured(N)=min_c sum_e (N_i+N_j)*sqrt(1+j_e²). Stationarity C^T[(N_i+N_j)j_e/sqrt(1+j_e²)]=0 reconstructs an exact shift potential.','exact_Hessian_formula':'H_N=-G A^-1 G^T; G=|B| diag(j/sqrt(1+j²)) C, A=C^T diag((N_i+N_j)/(1+j²)^(3/2)) C. It is negative semidefinite and nonzero when G is nonzero.','lapse_Hessian_rank':rank,'lapse_Hessian_eigenvalues':ev.tolist(),'lapse_Hessian':H.tolist(),'lapses':N.tolist(),'momentum_source':p.tolist(),'edge_currents':j.tolist(),'shift_cycle_residual':float(max(abs(C.T@shift))),'stationary_residual':res,'equal_lapse_control':{'rank':int(np.linalg.matrix_rank(solve_currents(B,C,np.ones(80),p)[1],tol=1e-10))},'symmetry_obstruction':flag_orbit(),'tree_control':{'edges':79,'cycle_dimension':0,'reduced_lapse_Hessian_rank':0},'null_vectors':'N by degree-one homogeneity; point/line checkerboard because every edge joins the two parts.','scope':['This explicit nonlinear metric pair extension fails the lapse-linearity test after actual shift elimination on the native graph. It cannot inherit the linear11316 constraint count by assumption.','The numerical rank78 is for the stored generic source. A full Dirac count or every-source rank theorem is not claimed.','A chosen tree admits the usual ghost-free pair construction but changes the native graph and its symmetry. Multivielbein/higher interactions require a separate Lorentz and constraint audit.'],'prior_owners':['analysis/w33_pass11316_native_edge_fierz_pauli.py','analysis/w33_pass11289_cycle_gram_gluing_flux.py'],'primary_sources':['https://arxiv.org/abs/1410.7774','https://arxiv.org/abs/1109.3515']}
if __name__=='__main__':
 d=payload();(ROOT/'data/w33_pass11323_native_cycle_lapse_hessian.json').write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['lapse_Hessian_rank'])
