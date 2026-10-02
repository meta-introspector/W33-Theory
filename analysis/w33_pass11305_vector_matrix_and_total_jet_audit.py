"""All light-field physical vector cuts and a separate regulated total-jet audit."""
from pathlib import Path
import sys,json
import numpy as np
import mpmath as mp
from scipy.integrate import quad
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11287_quotient_self_energy_matrix as Q
import w33_pass11295_ward_matched_modulus_poles as P

def light_data():
    d=json.loads((ROOT/'data/w33_pass11287_quotient_self_energy_matrix.json').read_text())
    return np.array(d['Hessian']),np.array(d['cubic_tensor']),np.array(d['Weyl_mass_real'])+1j*np.array(d['Weyl_mass_imag']),np.array(d['Yukawa_real'])+1j*np.array(d['Yukawa_imag'])

def vector_channels(g=.5):
    p,N,H,*_=Q.geometry();S=np.column_stack([N,1j*N])/np.sqrt(2)
    hp=g*np.einsum('aij,j->ai',H,p);hs=g*np.einsum('aij,jk->aik',H,S)
    X=2*np.real(hp.conj()@hp.T);ev,R=np.linalg.eigh(X);use=ev>1e-10;ev=ev[use];R=R[:,use]
    derivatives=2*np.real(np.einsum('aix,bi->xab',hs.conj(),hp)+np.einsum('ai,bix->xab',hp.conj(),hs))
    H0,C,M,Y=light_data();sm,fm,O,groups=Q.cut_groups(H0,C,M,Y)
    kappa=np.einsum('xi,xab,ac,bd->icd',O,derivatives,R,R)
    current=2*g*np.imag(np.einsum('ix,aij,jy->axy',S.conj(),H,S))
    current=np.einsum('ac,axy,xi,yj->cij',R,current,O,O)
    assert np.max(abs(current+current.transpose(0,2,1)))<1e-10
    vc={}
    def add(kind,a,b,outer):
        key=(kind,round(float(a),10),round(float(b),10))
        if key not in vc:vc[key]=dict(kind=kind,a=float(a),b=float(b),weight=np.zeros((22,22)))
        vc[key]['weight']+=outer
    for a in range(70):
        for b in range(a,70):
            v=kappa[:,a,b];add('VV',ev[a],ev[b],np.outer(v,v)/(32*np.pi)*(1 if a==b else 2))
        for j in range(22):
            v=current[a,:,j];add('VS',ev[a],sm[j],np.outer(v,v)/(16*np.pi*ev[a]))
    channels=[c for c in vc.values() if np.linalg.norm(c['weight'])>1e-20]
    return sm,O,groups,channels,ev,kappa,derivatives

def coefficients(c):
    a,b=c['a'],c['b']
    if c['kind']=='VV':return np.array([2+(a+b)**2/(4*a*b),-(a+b)/(2*a*b),1/(4*a*b)])
    return np.array([(a-b)**2,-2*(a+b),1.])

def vector_rho(c,t):
    a,b=c['a'],c['b'];th=(np.sqrt(a)+np.sqrt(b))**2
    if t<=th:return np.zeros((22,22))
    beta=np.sqrt(max(0,1-2*(a+b)/t+(a-b)**2/t**2));v=coefficients(c)
    return beta*(v[0]+v[1]*t+v[2]*t*t)*c['weight']

def vector_dispersion(c,s):
    a,b=c['a'],c['b'];th=(np.sqrt(a)+np.sqrt(b))**2
    assert 0<s<th  # All nonzero probe channels are below threshold here.
    v=coefficients(c)
    def f(t):
        beta=np.sqrt(max(0,1-2*(a+b)/t+(a-b)**2/t**2))
        return beta*(v[0]+v[1]*t+v[2]*t*t)/(t**3*(t-s))
    return s**3/np.pi*quad(f,th,np.inf,epsabs=1e-9,limit=200)[0]*c['weight']

def regulated_hermite(nodes,delta,constant):
    mp.mp.dps=90;nodes=[mp.mpf(str(x)) for x in nodes];d=mp.mpf(str(delta));mu2=mp.mpf('.01')
    f=lambda x:(x+d)**2*(mp.log((x+d)/mu2)-mp.mpf(str(constant)))
    n=2*len(nodes);A=mp.matrix(n);b=mp.matrix(n,1)
    for i,x in enumerate(nodes):
        for j in range(n):A[2*i,j]=x**j;A[2*i+1,j]=j*x**(j-1) if j else 0
        b[2*i]=mp.diff(f,x);b[2*i+1]=mp.diff(f,x,2)
    c=mp.lu_solve(A,b);error=float(max(abs(x) for x in A*c-b))
    assert error<1e-35
    return [str(mp.mpf(0))]+[str(c[i]/(i+1)) for i in range(n)],error

def total_jet_audit(ev,derivatives):
    H,C,M,Y=light_data();m=.01;w,O=np.linalg.eigh(H);w[abs(w)<1e-8]=0;sm=m*m*w
    P0=O[:,w<1e-8]@O[:,w<1e-8].T;CG=np.einsum('ij,ajk,kl->ail',P0,C,P0)
    IR=m**4*np.einsum('aij,bji->ab',CG,CG)/(32*np.pi**2)
    eigen=np.linalg.eigvalsh(IR);assert min(eigen)>-1e-14
    W=M.conj().T@M*m*m;fw,FU=np.linalg.eigh(W)
    _,_,Hg,*_=Q.geometry();p,N=Q.Q.reference();acts=.5*np.einsum('aij,j->ai',Hg,p)
    X=2*np.real(acts.conj()@acts.T);vw,VU=np.linalg.eigh(X);vw[abs(vw)<1e-10]=0
    delta=1e-8;fp=lambda x,c:2*(x+delta)*(np.log((x+delta)/.01)-c+.5)
    gs=np.einsum('aij,ik,jk,k->a',m*m*C,O,O,fp(sm,1.5))
    fd=np.array([m*m*(M.conj().T@y+y.conj().T@M) for y in Y])
    gf=np.einsum('aij,ik,jk,k->a',fd,FU.conj(),FU,fp(fw,1.5)).real
    gv=np.einsum('aij,ik,jk,k->a',derivatives,VU,VU,fp(vw,5/6))
    grad=(gs-2*gf+3*gv)/(64*np.pi**2)
    rows=[]
    for kind,masses,const in [('scalar',sm,1.5),('Weyl',fw,1.5),('vector',vw,5/6)]:
        nodes=sorted(set(round(float(x),10) for x in masses));coef,error=regulated_hermite(nodes,delta,const)
        rows.append({'kind':kind,'nodes':nodes,'polynomial':coef,'jet_error':error})
    return {'regulator':delta,'total_CW_tadpole':grad.tolist(),'tadpole_norm':float(np.linalg.norm(grad)),
      'Goldstone_IR_log_Hessian_coefficient':IR.tolist(),'IR_rank':int(sum(eigen>1e-12)),'IR_eigenvalues':eigen.tolist(),
      'regulated_counterfunctions':rows,
      'functional':'At fixeddelta, use -[Tr Ps(Xs)-2Tr Pf(M†M)+3Tr Pv(Xv)]/(64pi²). Each Pprime,Psecond matches the regulated CW function at all vacuum spectral nodes. This gives invariant spectral first/second-jet matching for the named covariant mass maps.',
      'IR_boundary':'The scalar Goldstone block has nonzero first mass variation, unlike the vector Gram kernel. Its logdelta Hessian coefficient has rank11. No finite polynomial gives an unregulated C² match at the tree vacuum: Ward-consistent Goldstone resummation and UV matching remain required.',
      'scope':'Actual total first-jet contractions plus formal spectral second-jet matching at a regulator. Scalar quartic mass derivatives and a nonlinear background implementation are not numerically built here. This is separate from the spacelike dispersion matching below, not a single finished UV renormalization scheme.'}

def payload():
    sm,O,light,vc,ev,kappa,derivatives=vector_channels();C,_,P0,_=P.ward_matching(light,sm)
    # VV threshold safety is global; VS matrices may have low thresholds only in low external blocks.
    rows=[];gold_error=0.
    for t in sorted(set(round(float(x),12) for x in sm if x>1e-10)):
        ix=np.where(abs(sm-t)<1e-10)[0];V=np.zeros((22,22));rho=np.zeros((22,22));active=0
        for c in vc:
            block=c['weight'][np.ix_(ix,ix)]
            if np.linalg.norm(block)<1e-18:continue
            th=(np.sqrt(c['a'])+np.sqrt(c['b']))**2
            assert t<th
            V+=vector_dispersion(c,t);rho+=vector_rho(c,t);active+=1
        S=P.sigma(light,t,C)+V;poles=t-np.linalg.eigvals(S[np.ix_(ix,ix)])
        assert max(poles.imag)<1e-10 and np.max(abs(rho))<1e-10
        rows.append({'tree_mass_squared':t,'multiplicity':len(ix),'pole_squared_real':poles.real.tolist(),'pole_squared_imag':poles.imag.tolist(),'vector_block_norm':float(np.linalg.norm(V[np.ix_(ix,ix)])),'nonzero_vector_channels':active})
    for c in vc:gold_error=max(gold_error,float(np.linalg.norm(c['weight']@P0)))
    vv_gold=max(float(np.linalg.norm(c['weight']@P0)) for c in vc if c['kind']=='VV')
    return {'status':'PASS','result_scope':'PASS_FULL_LIGHT_FIELD_VECTOR_CUT_MATRIX_AND_REGULATED_TOTAL_JET_AUDIT',
      'VV_and_VS_channel_groups':len(vc),'massive_vectors':len(ev),'Goldstone_VV_projection_error':vv_gold,'Goldstone_VS_coupling_norm':gold_error,'rows':rows,
      'vector_vertices':'Kappa_iab=dX_ab/dq_i from actual E6 vector Gram matrix; derivative Vh_ih_j vertex T_aij=2g Im(S_i†H_a S_j), projected onto70 massive vectors and22 light scalar eigenmodes.',
      'spectral_formulas':'VV: beta[2+(t-a-b)²/(4ab)]kappa_i kappa_j/(32pi), twice for unequal vector labels. VS: beta[(t-a-b)²-4ab]T_i T_j/(16pi a). Each channel is a positive outer product.',
      'matching':'Add three-times-zero-subtracted physicalVV/VS dispersions to prior11295 scalar/Weyl spacelike-matched bubbles. This computes conditional leading14 massive block poles; all nonzero channels are below their on-shell thresholds. VV and VS Goldstone couplings need not vanish: quotient Goldstone directions can include E6 compensators rotating the vector basis. The three-subtracted vector contribution hasSigma(0)=0 as a matching condition.',
      'total_jet_audit':total_jet_audit(ev,derivatives),
      'boundaries':['Earlier condensate EFT and projected physical vector cuts; not the enlarged11303 Higgs spectrum.','Gauginos/heavy chiral thresholds, contact/ghost analytic pieces, kinetic curvature corrections and complete UV matching remain open.','The regulator jet audit and conditional dispersion pole scheme are distinct calculations; neither predicts observed masses or the CC.'],
      'prior_owners':['analysis/w33_pass11287_quotient_self_energy_matrix.py','analysis/w33_pass11295_ward_matched_modulus_poles.py','analysis/w33_pass11300_invariant_vector_matching.py'],
      'primary_sources':['https://arxiv.org/abs/1609.06977','https://arxiv.org/abs/1910.02094','https://arxiv.org/abs/1808.07615']}
if __name__=='__main__':
    out=payload();(ROOT/'data/w33_pass11305_vector_matrix_and_total_jet_audit.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
