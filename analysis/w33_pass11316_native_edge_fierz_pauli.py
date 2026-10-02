"""Explicit linear spin-two field architecture on the native Levi graph."""
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
from w33_pass11289_cycle_gram_gluing_flux import cycle_basis
ETA=np.diag([-1.,1.,1.,1.])
B=[]
for i in range(4):
 for j in range(i,4):
  h=np.zeros((4,4));h[i,j]=h[j,i]=1 if i==j else 1/np.sqrt(2);B.append(h)
B=np.array(B)
def einstein(h,k):
 ku=ETA@k;k2=k@ku;tr=np.trace(ETA@h);hk=h@ku
 return (k2*h+np.outer(k,k)*tr-np.outer(k,hk)-np.outer(hk,k)-ETA*(k2*tr-ku@h@ku))/2

def operator(k,m2=0):
 return np.column_stack([np.einsum('aij,ij->a',B,einstein(h,k)+m2*(h-ETA*np.trace(ETA@h))/2) for h in B])
def gauge(k):
 return np.column_stack([np.einsum('aij,ij->a',B,np.outer(k,v)+np.outer(v,k)) for v in np.eye(4)])
def payload():
 inc,_=cycle_basis();L=inc@inc.T;w=np.linalg.eigvalsh(L.astype(float));w[abs(w)<1e-10]=0;unique=[]
 for x in [0,4-np.sqrt(6),4,4+np.sqrt(6),8]:unique.append({'mass_squared_over_m0_squared':float(x),'multiplicity':int(sum(abs(w-x)<1e-8))})
 assert [r['multiplicity'] for r in unique]==[1,24,30,24,1]
 k=np.array([1.,0,0,1]);E=operator(k);G=gauge(k);err=float(max(abs(E@G).flat));assert err<1e-12
 assert np.linalg.matrix_rank(E,tol=1e-10)==4 and np.linalg.matrix_rank(G,tol=1e-10)==4
 plus=np.diag([0.,1.,-1.,0.])/np.sqrt(2);cross=np.zeros((4,4));cross[1,2]=cross[2,1]=1/np.sqrt(2)
 pol=np.column_stack([np.einsum('aij,ij->a',B,h) for h in [plus,cross]])
 assert max(abs(E@pol).flat)<1e-12 and np.linalg.matrix_rank(np.column_stack([G,pol]),tol=1e-10)==6
 km=np.array([np.sqrt(2),0.,0.,1.]);M=operator(km,1.);assert np.linalg.matrix_rank(M,tol=1e-10)==5
 return {'status':'PASS','result_scope':'PASS_CONSTRUCTED_LINEAR_NATIVE_EDGE_SPIN_TWO_ARCHITECTURE','action':'Supply4D Minkowski base and80 symmetric tensor fields h_i. Sum standard positive-sign linearized Einstein actions. Add local edge potential -(m0²/4)sum_edges[(h_i-h_j)_munu(h_i-h_j)^munu-(trace(h_i-h_j))²]. Internal graph is the actual80-site160-edge Levi graph.','native_graph_spectrum':unique,'massless_equation_rank':4,'massless_gauge_rank':4,'massless_solution_dimension':6,'massless_physical_polarizations':2,'massive_equation_rank':5,'massive_physical_polarizations_per_graph_mode':5,'total_linear_polarizations':2+79*5,'prior_tensor_control':{'source':'manuscripts/parts/PART_CCLXIX_GRAVITON_EMERGENCE.md','symmetric_traceless_dimension':299,'natural_additive_tensor_Laplacian_eigenvalue':20,'natural_additive_zero_modes':0,'scope':'The old24-dimensional adjacency eigenspace has point-Laplacian eigenvalue10. Its symmetric traceless tensor square has eigenvalue20 under L tensorI+I tensorL. A different gauge operator needs an explicit construction; an internal tensor dimension does not establish Lorentz helicity.'},'Bianchi_gauge_error':err,'mass_gap_over_m0':float(np.sqrt(4-np.sqrt(6))),'explicit_massless_polarizations':[plus.tolist(),cross.tolist()],'boundary':['The4D Lorentzian base and tensor action are added inputs, not an emergence theorem from the finite graph.','Uniform graph mode carries two massless polarizations;79 massive modes carry five each. Below the declared graph gap only the uniform mode propagates.','Fierz-Pauli tuning removes the extra scalar at quadratic order; nonlinear ghost freedom, cosmological stability and interacting constraint closure are not established.','No Newton coefficient, light-cone refinement, observed graviton scale or full nonlinear Einstein limit is derived.'],'prior_owners':['analysis/w33_pass11306_relational_graph_clock_constraints.py','analysis/w33_pass11288_metric_constraint_audit.py','analysis/w33_pass11289_cycle_gram_gluing_flux.py'],'primary_sources':['https://arxiv.org/abs/1109.3515']}
if __name__=='__main__':
 d=payload();(ROOT/'data/w33_pass11316_native_edge_fierz_pauli.json').write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['native_graph_spectrum'])
