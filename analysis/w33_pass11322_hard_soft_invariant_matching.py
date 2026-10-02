"""Finite invariant hard-sector jets; retain the actual rank-eleven soft IR boundary."""
from pathlib import Path
import json,sys
import numpy as np
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
from w33_pass11305_vector_matrix_and_total_jet_audit import light_data
import w33_pass11287_quotient_self_energy_matrix as Q

def hermite(nodes,constant,soft_zero):
    mp.mp.dps=100;nodes=[mp.mpf(str(x)) for x in nodes]
    if soft_zero:nodes=[mp.mpf(0)]+nodes
    n=2*len(nodes);A=mp.matrix(n);b=mp.matrix(n,1)
    for i,x in enumerate(nodes):
        for j in range(n):A[2*i,j]=x**j;A[2*i+1,j]=j*x**(j-1) if j else 0
        b[2*i]=2*x*(mp.log(x/mp.mpf('.01'))-constant+mp.mpf('.5')) if x else 0
        b[2*i+1]=2*mp.log(x/mp.mpf('.01'))-2*constant+3 if x else 0
    c=mp.lu_solve(A,b);error=float(max(abs(v) for v in A*c-b));assert error<1e-50
    return [str(mp.mpf(0))]+[str(c[i]/(i+1)) for i in range(n)],error

def unique(x):return sorted(set(round(float(y),12) for y in x if y>1e-10))

def payload():
    H,C,M,Y=light_data();m=.01;sw,O=np.linalg.eigh(m*m*H);sw[abs(sw)<1e-10]=0
    W=m*m*M.conj().T@M;fw,FU=np.linalg.eigh(W)
    p,N,Hg,*_=Q.geometry();S=np.column_stack([N,1j*N])/np.sqrt(2);hp=.5*np.einsum('aij,j->ai',Hg,p);hs=.5*np.einsum('aij,jk->aik',Hg,S)
    X=2*np.real(hp.conj()@hp.T);vw,VU=np.linalg.eigh(X);vw[abs(vw)<1e-10]=0
    vd=2*np.real(np.einsum('aix,bi->xab',hs.conj(),hp)+np.einsum('ai,bix->xab',hp.conj(),hs))
    fd=m*m*np.array([M.conj().T@y+y.conj().T@M for y in Y])
    fp=lambda w,c:np.array([2*x*(np.log(x/.01)-c+.5) if x>1e-10 else 0 for x in w])
    gs=np.einsum('aij,ik,jk,k->a',m*m*C,O,O,fp(sw,1.5));gf=np.einsum('aij,ik,jk,k->a',fd,FU.conj(),FU,fp(fw,1.5)).real
    gv=np.einsum('aij,ik,jk,k->a',vd,VU,VU,fp(vw,5/6));grad=(gs-2*gf+3*gv)/(64*np.pi**2)
    P0=O[:,sw==0]@O[:,sw==0].T;CG=np.einsum('ij,ajk,kl->ail',P0,m*m*C,P0);IR=np.einsum('aij,bji->ab',CG,CG)/(32*np.pi**2);rank=int(sum(np.linalg.eigvalsh(IR)>1e-12));assert rank==11
    rows=[];jet_error=0.
    for kind,w,U,D,c,soft in [('scalar',sw,O,m*m*C,1.5,True),('Weyl',fw,FU,fd,1.5,False),('vector_hard',vw,VU,vd,5/6,True)]:
        mp.mp.dps=100;c=mp.mpf(5)/6 if kind=="vector_hard" else mp.mpf("1.5")
        nodes=unique(w);coef,err=hermite(nodes,c,soft)
        # Fréchet Hessian of Tr f(X) uses divided differences of fprime,
        # including hard/soft cross blocks, plus Tr fprime(X) Xsecond.
        # Jet agreement therefore cancels every unknown Xsecond contribution too.
        allnodes=([0.] if soft else [])+nodes
        def poly1(x):return sum(mp.mpf(coef[j])*j*mp.mpf(str(x))**(j-1) for j in range(1,len(coef)))
        def poly2(x):return sum(mp.mpf(coef[j])*j*(j-1)*mp.mpf(str(x))**(j-2) for j in range(2,len(coef)))
        f1=[mp.mpf(0) if x==0 else 2*mp.mpf(str(x))*(mp.log(mp.mpf(str(x))/mp.mpf('.01'))-mp.mpf(str(c))+mp.mpf('.5')) for x in allnodes]
        f2=[mp.mpf(0) if x==0 else 2*mp.log(mp.mpf(str(x))/mp.mpf('.01'))-2*mp.mpf(str(c))+3 for x in allnodes]
        for i,x in enumerate(allnodes):
            jet_error=max(jet_error,float(abs(poly1(x)-f1[i])),float(abs(poly2(x)-f2[i])))
            for j,y in enumerate(allnodes):
                if i!=j:jet_error=max(jet_error,float(abs((poly1(x)-poly1(y)-f1[i]+f1[j])/mp.mpf(str(x-y)))))
        rows.append({'kind':kind,'hard_nodes':nodes,'node_precision':'numerical mass representatives rounded to12 decimal places; interpolation precision applies at these nodes','soft_derivatives_set_to_zero':soft,'P_coefficients_ascending':coef,'Hermite_error':err,'degree':len(coef)-1})
    assert jet_error<1e-40
    return {'status':'PASS','result_scope':'PASS_REGULATOR_FREE_INVARIANT_HARD_TADPOLE_MATCHING_NOT_FULL_SOFT_RESUMMATION','mass_maps':'Xs=scalar covariant Hessian, W=M†M, Xv=actual E6 vector Gram. Their stored vacuum masses and first variations are used. Nonlinear scalar second variations are not numerically reconstructed.','functional':'Vct=-[Tr Ps(Xs)-2Tr Pf(W)+3Tr Pv(Xv)]/(64pi²). Ps and Pv match the positive gapped spectral branches and have Pprime(0)=Psecond(0)=0. Pf matches every positive Weyl node. Additive constants are unrestricted. This is an invariant higher-dimensional local EFT functional, not a renormalizable UV counterterm.','hard_soft_definition':'In a sufficiently small gap-preserving neighborhood, soft clusters are omitted from the hard potential via covariant spectral projectors. This defines a smooth hard function identically zero near the soft cluster; it is not the full unregulated CW function at zero. Vector soft contributions start at fourth field order and do not affect first/second field jets.','counterfunctions':rows,'maximum_first_second_divided_difference_error':jet_error,'actual_unregulated_first_jet':grad.tolist(),'first_jet_norm':float(np.linalg.norm(grad)),'soft_IR_log_coefficient_rank':rank,'soft_IR_log_eigenvalues':np.linalg.eigvalsh(IR).tolist(),'Ward_matching':'The hard CW first and second field jets cancel identically against Vct by the spectral Fréchet formula. Thus the supplied stationary tree vacuum has no residual one-loop hard tadpole in this chosen EFT matching prescription. This is input renormalization, not radiative vacuum selection.','remaining_obligations':['The rank11 scalar soft logarithmic Hessian survives; no regulator-free C² match of the full static CW potential is claimed.','A covariant nonlinear completion of the22-field scalar/background mass maps and quotient kinetic geometry is not built here. The invariant field-functional statement is conditional on that completion; the spectral counterfunctions, actual first-jet contraction and mass-jet identity are explicit.',
      'Full hard self-energy Ward resummation, two-loop double-counting control and total physical poles remain open.','UV completion, observed masses/mixing and absolute vacuum energy remain inputs or open problems.'],'prior_owners':['analysis/w33_pass11300_invariant_vector_matching.py','analysis/w33_pass11305_vector_matrix_and_total_jet_audit.py','analysis/w33_pass11315_regulator_free_goldstone_poles.py'],'primary_sources':['https://arxiv.org/abs/1406.2355','https://arxiv.org/abs/1609.06977']}
if __name__=='__main__':
    d=payload();(ROOT/'data/w33_pass11322_hard_soft_invariant_matching.json').write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['maximum_first_second_divided_difference_error'])
