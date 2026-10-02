"""A constructed multi-clock first-class completion of the native quadratic graph."""
from pathlib import Path
import sys,json
import numpy as np
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
from w33_pass11289_cycle_gram_gluing_flux import cycle_basis

def native_generators():
    inc,_=cycle_basis();edges=[tuple(np.where(inc[:,j])[0]) for j in range(160)];n=80
    J=np.block([[np.zeros((n,n)),np.eye(n)],[-np.eye(n),np.zeros((n,n))]])
    K=np.zeros((n,2*n,2*n))
    for i in range(n):K[i,i,i]=1;K[i,n+i,n+i]=1
    for i,j in edges:
        v=np.zeros(2*n);v[i]=1;v[j]=-1
        K[i]+=.5*np.outer(v,v);K[j]+=.5*np.outer(v,v)
    return J,np.einsum('ij,ajk->aik',J,K),edges

def connection(t,indices,A,J):
    U=np.eye(160);out={}
    for i in range(80):
        if i in indices:out[i]=U@A[i]@(-J@U.T@J)
        if t[i]:U=U@expm(t[i]*A[i])
    return U,out

def payload():
    J,A,edges=native_generators();a,b=edges[0];other=next(i for i,j in edges if j==b and i!=a);ix=list(map(int,sorted([a,b,other])))
    t=np.zeros(80);t[ix]=[.017,.023,.031];U,B=connection(t,ix,A,J)
    sym=float(np.max(abs(U.T@J@U-J)));assert sym<1e-11
    h=1e-5;errors=[]
    for i in ix:
        for j in ix:
            if i>=j:continue
            ei=np.zeros(80);ei[i]=h;ej=np.zeros(80);ej[j]=h
            _,bp=connection(t+ei,ix,A,J);_,bm=connection(t-ei,ix,A,J);di=(bp[j]-bm[j])/(2*h)
            _,bp=connection(t+ej,ix,A,J);_,bm=connection(t-ej,ix,A,J);dj=(bp[i]-bm[i])/(2*h)
            err=float(np.max(abs(di-dj-(B[i]@B[j]-B[j]@B[i]))));assert err<1e-8;errors.append(err)
    i,j=min(a,b),max(a,b);u=.017;v=.031
    E=expm(u*A[i]);F=expm(v*A[j]);dressed=expm(v*(E@A[j]@(-J@E.T@J)))@E
    assert np.max(abs(dressed-E@F))<1e-11
    naive=float(np.linalg.norm(E@F-F@E));assert naive>1e-6
    # Exact long magic circuits can remain reversible: T has finite order9.
    z=np.exp(2j*np.pi/9);T=np.diag([1,z,z**-1]);period=float(np.max(abs(np.linalg.matrix_power(T,9)-np.eye(3))));assert period<1e-12
    return {'status':'PASS','result_scope':'PASS_NATIVE_GRAPH_RELATIONAL_CLOCK_FIRST_CLASS_COMPLETION',
      'native_graph':{'vertices':80,'edges':160},'clock_constraints':'Add80 canonical(t_i,P_i) pairs. U(t)=ordered product exp(t_i A_i), A_i=J K_i of the native quadratic site Hamiltonians. B_i=(partial_i U)U^-1. Define C_i=P_i+(1/2)z^T(-J B_i)z.',
      'exact_closure':'partial_i B_j-partial_j B_i-[B_i,B_j]=0 by right Maurer-Cartan. Hence{C_i,C_j}=0 exactly. The constraint Jacobian contains the80x80 identity in the P columns; rank80. BRST charge Omega=sum c_i C_i is nilpotent without higher ghost terms.',
      'relational_observable':'z0=U(t)^-1 z is invariant along every constraint flow and supplies80 physical canonical pairs. The dressed connection transports the old edge-current bracket into explicit clock dependence; it does not set the old bracket to zero.',
      'phase_space_count':{'total_real_phase_dimension':320,'first_class_constraints':80,'physical_real_phase_dimension':160},
      'finite_controls':{'vertices':ix,'flatness_errors':errors,'symplectic_error':sym,'undressed_path_difference':naive,'dressed_path_difference':float(np.max(abs(dressed-E@F)))},
      'parallel_intake_control':{'T_order9_identity_error':period,'boundary':'The sampled decay in11312 is probabilistic, not every sufficiently long word: all-Clifford-identity words with9m cubic gates equalI. Haar measure zero alone does not prove the random-walk limit or its geometric rate. The11309 affine-subspace claim must be read against its evolving certificate, not the stale hyperplane docstring.'},
      'boundaries':['This is a parametrized oscillator theory with an explicit first-class algebra, not dynamical spacetime or two graviton polarizations.','Quadratic site potential is used; the quartic11301 scalar interaction is not included in the matrix witness.','Clock ordering and normalizations are additional choices; no Lorentzian3D refinement or local Einstein action is derived.'],
      'prior_owners':['analysis/w33_pass11301_graph_adm_regulator.py','analysis/w33_pass11296_covariant_matter_constraints.py','analysis/PASS11312_DEPTH_LAW.md','analysis/PASS11309_ONE_GATE_LAW.md'],
      'primary_sources':['https://arxiv.org/abs/1011.2463']}
if __name__=='__main__':
    out=payload();(ROOT/'data/w33_pass11306_relational_graph_clock_constraints.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
