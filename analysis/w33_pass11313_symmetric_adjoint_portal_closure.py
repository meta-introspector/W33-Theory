"""Complete quartic beta operator for the real SO10 54+45 restriction."""
from pathlib import Path
import sys,json
import numpy as np
from scipy.optimize import root
from fractions import Fraction
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))

def bases(n=10):
    d=np.linalg.qr(np.column_stack([np.ones(n),np.eye(n)[:,:n-1]]))[0][:,1:]
    S=[np.diag(x) for x in d.T];K=[]
    for i in range(n):
        for j in range(i+1,n):
            a=np.zeros((n,n));a[i,j]=a[j,i]=1/np.sqrt(2);S.append(a)
            b=np.zeros((n,n));b[i,j]=1/np.sqrt(2);b[j,i]=-1/np.sqrt(2);K.append(b)
    return np.array(S),np.array(K)
BS,BK=bases();LABELS=['Sdouble','Ssingle','Kdouble','Ksingle','norm_portal','S2K2','SKSK']

def invariants(S,K):
    q=np.trace(S@S);r=-np.trace(K@K)
    return np.array([q*q/4,np.trace(S@S@S@S)/4,r*r/4,np.trace(K@K@K@K)/4,q*r/2,-np.trace(S@S@K@K)/2,-np.trace(S@K@S@K)/2])

def trpair(a,b):return np.einsum('aij,bji->ab',a,b)

def hessians(S,K):
    nS=len(BS);nK=len(BK);N=nS+nK;out=np.zeros((7,N,N));ss=np.s_[:nS,:nS];kk=np.s_[nS:,nS:]
    q=np.trace(S@S);r=-np.trace(K@K);a=2*np.einsum('ij,aji->a',S,BS);b=-2*np.einsum('ij,aji->a',K,BK)
    out[0][ss]=np.outer(a,a)/2+q*np.eye(nS);out[2][kk]=np.outer(b,b)/2+r*np.eye(nK)
    out[1][ss]=trpair(BS,S@S@BS+S@BS@S+BS@S@S)
    out[3][kk]=trpair(BK,K@K@BK+K@BK@K+BK@K@K)
    out[4][ss]=r*np.eye(nS);out[4][kk]=q*np.eye(nK);out[4,:nS,nS:]=np.outer(a,b)/2
    out[5][ss]=-trpair(BS,BS@K@K+K@K@BS)/2
    out[5][kk]=-trpair(BK,BK@S@S+S@S@BK)/2
    out[5,:nS,nS:]=-trpair(S@BS+BS@S,K@BK+BK@K)/2
    out[6][ss]=-trpair(BS,K@BS@K);out[6][kk]=-trpair(BK,S@BK@S)
    out[6,:nS,nS:]=-trpair(BS,BK@S@K+K@S@BK)
    for h in out:h[nS:,:nS]=h[:nS,nS:].T
    return out

def vector_gram(S,K):
    AS=BK@S-S@BK;AK=BK@K-K@BK
    return np.einsum('aij,bij->ab',AS,AS)+np.einsum('aij,bij->ab',AK,AK)

def sample(rng):
    S=np.einsum('a,aij->ij',rng.normal(size=54)/np.sqrt(54),BS)*rng.uniform(.4,1.5)
    K=np.einsum('a,aij->ij',rng.normal(size=45)/np.sqrt(45),BK)*rng.uniform(.4,1.5)
    return S,K

def coefficients():
    rng=np.random.default_rng(11313);X=[];ys=[];gs=[]
    for _ in range(48):
        S,K=sample(rng);H=hessians(S,K);X.append(invariants(S,K));ys.append(np.einsum('aij,bji->ab',H,H)/2);G=vector_gram(S,K);gs.append(1.5*np.trace(G@G))
    X=np.array(X);raw=np.linalg.lstsq(X,np.array(ys).reshape(48,49),rcond=None)[0].reshape(7,7,7);graw=np.linalg.lstsq(X,gs,rcond=None)[0]
    def rat(v):return float(Fraction(float(v)).limit_denominator(1800))
    C=np.vectorize(rat)(raw);D=np.vectorize(rat)(graw);error=0.
    for _ in range(12):
        S,K=sample(rng);H=hessians(S,K);I=invariants(S,K);G=vector_gram(S,K)
        error=max(error,float(np.max(abs(np.einsum('a,aij->ij',I,C)-np.einsum('aij,bji->ab',H,H)/2))),abs(I@D-1.5*np.trace(G@G)))
    assert error<1e-9
    return C,D,error

def exact_coefficients():
    import inspect
    import sympy as s
    SB=[];KB=[];norm=[]
    for j in range(1,10):
     SB.append(np.diag([1]*j+[-j]+[0]*(9-j)));norm.append(j*(j+1))
    for i in range(10):
     for j in range(i+1,10):
      a=np.zeros((10,10),int);a[i,j]=a[j,i]=1;SB.append(a);norm.append(2)
      b=np.zeros((10,10),int);b[i,j]=1;b[j,i]=-1;KB.append(b)
    SB=np.array(SB);KB=np.array(KB);weights=np.array([2520//x for x in norm]+[1260]*45,dtype=np.int64)
    text=inspect.getsource(hessians).replace('q*np.eye(nS)','q*np.diag(normS)').replace('r*np.eye(nK)','r*np.diag(normK)').replace('r*np.eye(nS)','r*np.diag(normS)').replace('q*np.eye(nK)','q*np.diag(normK)')
    ns={'np':np,'BS':SB,'BK':KB,'normS':np.array(norm),'normK':np.array([2]*45),'trpair':trpair};exec(text,ns)
    rng=np.random.default_rng(11313);xs=[];ts=[];ds=[]
    for b in range(8):
     S=rng.integers(-1,2,(10,10));S=np.triu(S)+np.triu(S,1).T;S[-1,-1]=-np.trace(S[:-1,:-1]);K=rng.integers(-1,2,(10,10));K=np.triu(K,1)-np.triu(K,1).T
     I4=4*invariants(S,K);assert np.max(abs(I4-np.rint(I4)))==0;xs.append(list(map(int,I4)))
     h2=2*ns['hessians'](S,K);assert np.max(abs(h2-np.rint(h2)))==0;h2=np.array(h2,np.int64)
     assert int(np.max(abs(h2)))**2*int(max(weights))**2*99**2<np.iinfo(np.int64).max
     T=np.einsum('aij,bij,i,j->ab',h2,h2,weights,weights,optimize=True)
     ts.append([s.Rational(int(x),8*2520**2) for x in T.flat])
     X=KB@S-S@KB;Y=KB@K-K@KB;V=np.einsum('aij,bij->ab',X,X)+np.einsum('aij,bij->ab',Y,Y);ds.append(s.Rational(3*int(np.sum(V*V)),8))
    X=s.Matrix(xs[:7])/4;det=X.det();raw=X.inv()*s.Matrix(ts[:7]);dg=X.inv()*s.Matrix(ds[:7]);assert s.Matrix([xs[7]])/4*raw==s.Matrix([ts[7]])
    assert (s.Matrix([xs[7]])/4*dg)[0]==ds[7]
    
    return np.array(raw,float).reshape(7,7,7),np.array(dg,float).ravel(), {"evaluation_matrix_determinant":str(det),"scalar_tensor":[[str(v) for v in row] for row in raw.tolist()],"gauge_source":[str(v) for v in dg],"proof":"Seven independent integer backgrounds determine all seven degree-four invariants; exact rational coordinate metrics, integer contractions and an eighth exact replay. Float projections on twelve further backgrounds are independent controls."}

def beta(v,C,D,R=365/123,b=8/3,active=2):
    waves=np.array([120-4*active*R]*2+[96]*2+[108-2*active*R]*3)-2*b
    z=np.einsum('aij,i,j->a',C,v,v)+D-waves*v;z[1]-=4*active*R*R
    return z

def dual_yukawa_beta(v,C,D,active=2,total=8):
    b=8-2*total/3;R=(27-b)/(9.7+active)
    wave=np.array([120-4*active*R]*2+[96-4*active*R]*2+[108-4*active*R]*3)-2*b
    z=np.einsum('aij,i,j->a',C,v,v)+D-wave*v
    z[[1,3,5,6]]-=active*R**2*np.array([4,4,8,4])
    return z

def yukawa_contraction_control():
    eps=np.array([[0.,1.],[-1.,0.]])
    def coefficient(ys,yk,ix):
        Y=np.array([np.kron(np.eye(2),s)*ys for s in BS]+[np.kron(eps,k)*yk for k in BK])
        a=Y[ix];square=np.einsum('aji,ajk->ik',Y,Y)
        val=(square.T@a+a@square)/2+2*np.einsum('aij,jk,akl->il',Y,a.T,Y)
        val+=np.einsum('a,aij->ij',np.einsum('ij,aji->a',a,Y),Y)
        return float(np.sum(val*a)/np.sum(a*a))
    rows=[coefficient(1,0,0),coefficient(0,1,54),coefficient(1,1,0),coefficient(1,1,54)]
    assert np.max(abs(np.array(rows)-[8.2,7.5,11.7,11.7]))<1e-10
    rng=np.random.default_rng(45);S,K=sample(rng);M=np.kron(np.eye(2),S)+np.kron(eps,K)
    exact=2*(np.trace(S@S@S@S)+np.trace(K@K@K@K)-4*np.trace(S@S@K@K)-2*np.trace(S@K@S@K))
    error=abs(np.trace(M@M@M@M)-exact);assert error<1e-10
    return {'contractions':rows,'box_trace_error':float(error),'equal_Yukawa_AF_ratio':'730/351','scope':'Allowed antisymmetric flavor tensor for the same two active vector Weyls; Hermitian99 realification, not an SU10 gauge group. Equal Yukawas preserved at one loop; no full portal AF ray asserted.'}

def exact_portal_barriers(C):
    import sympy as s
    def rational(x):return s.Rational(str(Fraction(float(x)).limit_denominator(1800)))
    rows=[]
    for name,index,r,w,c in [('S_only_Yukawa_K_barrier',2,s.Rational(1,5),s.Rational(272,3),s.Rational(102,5)),('equal_dual_Yukawa_S_barrier',0,s.Rational(3,10),s.Rational(34408,351),s.Rational(1052092,41067))]:
        A=s.Matrix([[rational(C[index,i,j])+r*rational(C[index+1,i,j]) for j in range(7)] for i in range(7)])
        ix=[i for i in range(7) if any(A[i,j]!=0 for j in range(7))];sub=A.extract(ix,ix)
        minors=[s.factor(sub[:j,:j].det()) for j in range(1,len(ix)+1)];assert all(x>0 for x in minors)
        Q=A.extract([index,index+1],[index,index+1]);v=s.Matrix([1,r]);k=s.factor(1/(v.T*Q.inv()*v)[0]);minimum=s.factor(c-w*w/(4*k));assert minimum>0
        rows.append({'case':name,'functional_weights':[str(1),str(r)],'active_quartic_indices':ix,'positive_matrix':[[str(x) for x in row] for row in sub.tolist()],'leading_principal_minors':[str(x) for x in minors],'Riccati_coefficient':str(k),'linear_coefficient':str(w),'constant':str(c),'global_lower_bound':str(minimum),'global_lower_bound_numeric':float(minimum),'identity':'dz/dtau >= k*z²-w*z+c; all three mixed-quartic contributions are a positive definite quadratic form. Negative discriminant gives no real fixed ray and finite-tau blowup if one-loop evolution persists.'})
    # Actual backgrounds realize the functional trace ratios, excluding a fictitious invariant evaluation.
    a=s.sqrt((5+s.sqrt(5))/20);b=s.sqrt((5-s.sqrt(5))/20);S=s.diag(a,-a,b,-b,*([0]*6))
    ka=s.sqrt((10+s.sqrt(10))/60);kb=s.sqrt((5-s.sqrt(10))/30);K=s.zeros(10)
    for j,z in enumerate([ka,ka,kb]):K[2*j,2*j+1]=z;K[2*j+1,2*j]=-z
    assert s.simplify(s.trace(S*S))==1 and s.simplify(s.trace(S**4))==s.Rational(3,10)
    assert s.simplify(-s.trace(K*K))==1 and s.simplify(s.trace(K**4))==s.Rational(1,5)
    return {'barriers':rows,'backgrounds':'Exact unit-normalized S and K backgrounds have fourth trace3/10 and1/5. Positivity is checked by exact rational Sylvester minors after rational beta projection; coefficients are independently replayed on12 backgrounds.','scope':'b=8/3 and specified eight-Weyl Yukawa texture. Other omitted scalar Hessian blocks add nonnegative background TrH² contributions. Additional Yukawas or changed charged inventory can alter the barrier; this is not a no-go for every SO10 matter model.'}

def bounded_control(v):
    a=min(v[0]+v[1]/10,v[0]+v[1]*73/90)/4
    b=min(v[2]+v[3]/10,v[2]+v[3]/2)/4
    c=v[4]/2+min(0,v[5]-abs(v[6]))/4
    return a,b,c,c+2*np.sqrt(max(0,a*b))

def fixed_rays(C,D):
    rng=np.random.default_rng(113130);sols=[]
    for _ in range(240):
        x=np.array([.32,-.12,.2,.1,0,0,0])+rng.normal(size=7)*.2
        q=root(lambda z:beta(z,C,D),x,tol=1e-10)
        if not q.success or max(abs(beta(q.x,C,D)))>1e-7 or any(np.linalg.norm(q.x-z)<1e-5 for z in sols):continue
        sols.append(q.x)
    rows=[]
    for v in sols:
        a,b,c,margin=bounded_control(v);h=1e-5;J=np.column_stack([(beta(v+h*np.eye(7)[i],C,D)-beta(v-h*np.eye(7)[i],C,D))/(2*h) for i in range(7)])
        rows.append({'ratios':v.tolist(),'sufficient_boundedness':[a,b,c,margin],'certified_bounded_by_trace_inequalities':bool(a>0 and b>0 and margin>0),'quartic_UV_eigenvalues_real':np.linalg.eigvals(J).real.tolist(),'residual':float(max(abs(beta(v,C,D))))})
    return rows

def payload():
    C,D,error=coefficients();Ce,De,exact=exact_coefficients();assert np.max(abs(C-Ce))<1e-12 and np.max(abs(D-De))<1e-12;C,D=Ce,De;np.savez('/tmp/w33_11313_portal_coefficients.npz',C=C,D=D)
    rows=fixed_rays(C,D)
    return {'status':'PASS','result_scope':'PASS_COMPLETE_SEVEN_QUARTIC_SO10_54_45_BETA_RESTRICTION',
      'fields':'Canonical real symmetric-tracelessS54 and antisymmetricK45. SO generators are unit trace-norm antisymmetric matrices.',
      'quartics':LABELS,'basis':'[q²/4,TrS4/4,r²/4,TrK4/4,qr/2,-TrS²K²/2,-TrSKSK/2], q=TrS²,r=-TrK².',
      'exact_projection_certificate':exact,'scalar_beta_tensor':C.tolist(),'gauge_quartic_source':D.tolist(),'independent_projection_error':error,
      'operator':'16pi² betaV = (1/2)Tr(H_V²)+(3/2)Tr(X_vector²)-3g² sum_fields C2 phi.dV + gammaS S.dV -nactive*y4 TrS4. C2S=10,C2K=8,gammaS=nactive*y². Ratio flow adds2b times each quartic.',
      'gauge_budget':'b=8/3 imports the full11303 field budget; this7-quartic restriction does not remove charged frames from that budget.',
      'dual_yukawa_control':yukawa_contraction_control(),'exact_portal_barriers':exact_portal_barriers(C),'fixed_ray_search':rows,'search_scope':'240 deterministic root starts. Discovered roots are numerical witnesses; absence of any other root is not proved.',
      'full_inventory_boundary':'Seven invariants close for the54+45 two-field quartic restriction, not for allE6/frame/sextet/mediator portals. Gauge sources explicitly generate these mixed portals at zero; setting all omitted frame portals to zero does not close the full inventory.',
      'prior_owners':['analysis/w33_pass11303_scalar_source_and_uv_barrier.py','analysis/w33_pass11298_coupled_running_closure.py'],'primary_sources':['https://arxiv.org/abs/hep-ph/0211440','https://arxiv.org/abs/1705.00751']}
if __name__=='__main__':
 d=payload();(ROOT/'data/w33_pass11313_symmetric_adjoint_portal_closure.json').write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(d['status'],d['gauge_quartic_source']);print([(r['certified_bounded_by_trace_inequalities'],r['ratios']) for r in d['fixed_ray_search']])
