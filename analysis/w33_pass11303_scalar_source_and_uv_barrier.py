"""One scalar supplies two channels; a portal-independent54-sector UV barrier."""
from pathlib import Path
import sys,json
import numpy as np
import sympy as s
from scipy.sparse import csr_matrix
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_20261001_global_e6_cartan_covariants as I
import w33_pass11293_yukawa_uv_budget as H

def scalar_yukawa(y):
    _,d,_=I.tensors();Ys=[]
    for a in range(27):
        for i in range(3):
            for j in range(i,3):
                T=np.zeros((3,3));T[i,j]=T[j,i]=1 if i==j else 1/np.sqrt(2)
                Y=csr_matrix(y*np.kron(d[:,:,a],T)/np.sqrt(2));Ys.extend([Y,1j*Y])
    W=sum((Y.conj().T@Y for Y in Ys),start=csr_matrix((81,81),dtype=complex))
    out=[]
    for a in (0,17,323):
        Y=Ys[a];B=(W.conj()@Y+Y@W)/2
        for Z in Ys:B+=2*Z@Y.conj().T@Z+Z*float(Z.conj().multiply(Y).sum().real)
        c=complex(Y.conj().multiply(B).sum()/Y.conj().multiply(Y).sum());err=float(np.sqrt(abs((B-c*Y).conj().multiply(B-c*Y).sum())))
        assert abs(c-25*y*y)<1e-10 and err<1e-10;out.append([a,c.real,err])
    return out

def uv_barrier():
    x,y=s.symbols('x y');r=s.Rational(3,10);z=x+r*y
    f1=124*x*x+s.Rational(224,5)*x*y+s.Rational(159,50)*y*y-104*x+18
    f2=24*x*y+s.Rational(127,5)*y*y-104*y+60
    A=s.Matrix([[124,26],[26,s.Rational(54,5)]]);L=s.Matrix([1,r]);a=1/(L.T*A.inv()*L)[0]
    remainder=s.factor(f1+r*f2-(a*z*z-104*z+36))
    assert A.det()>0 and a==s.Rational(16580,159)
    assert remainder==(56*x-15*y)**2/s.Integer(159)
    minimum=s.factor(36-104**2/(4*a));assert minimum==s.Rational(41736,4145)>0
    D=float(4*a*36-104**2);bound=2/np.sqrt(D)*(np.pi/2-np.arctan((2*float(a)*.5-104)/np.sqrt(D)))
    return {'weighted_polynomial':str(s.factor(f1+r*f2)),'positive_square_remainder':str(remainder),'Riccati_coefficient':str(a),'strict_minimum':str(minimum),'example_z0':.5,'finite_rescaled_time_blowup_bound':float(bound),
      'proof':'For unit tracelessS, r=TrS4/(TrS²)² has controls1/10 and73/90. Weight23/32 and9/32 givesr=3/10. Other scalar blocks/mixed blocks add nonnegative Tr(Hfull²)-Tr(HSS²) to each background beta. With no SO10-charged Weyls, S has no Yukawa and no negative fermion-box term. Thus dz/dtau >= (16580/159)z²-104z+36 >0, tau=integral gSO² dt/(16pi²). Since bSO=8 makes tau divergent in the UV, the Riccati comparison blows up at finite tau for every finite initialz.',
      'scope':'A one-loop obstruction to complete weak-coupling asymptotic freedom of the unchanged11293 real54 inventory, including arbitrary scalar portals. Not a no-go for added SO10 Yukawas, a different representation, asymptotic safety, gravity corrections or strong UV completion.'}

def adjoint_reencoding():
    v,A,U,V,K,S=H.economical_reference();L=S@K
    assert np.max(abs(L+L.T))<1e-12
    err=max(np.max(abs(L@L+K@L+6*np.eye(10))),np.max(abs(A@U-1j*U@L)),abs(np.trace(K@L)))
    assert err<1e-12 and np.max(abs(-K@L-S))<1e-12
    return {'map':'L=S K is a second real SO10 adjoint45; inverse S=-K L on K²=-I and[K,L]=0. Replace S equations by L²+K L+6I=0, Tr(KL)=0 and AU=iUL. All constraints remain quadratic.',
      'reference_error':float(err),'local_kernel':'The inverse reconstructs the old symmetric tracelessS locally. The prior111-dimensional orbit kernel is retained; replacing54 by45 leaves1356 real alignment fields and1245 positive normal directions. This is local elimination using11293, not a global vacuum classification.',
      'gauge_beta_with_one_scalar_mediator':{'E6':8,'SO10':'26/3'},'boundary':'An explicit alternative representation that escapes the54-specific theorem; full two-adjoint/portal UV stability is not established.'}

def fermion_escape():
    # Real SO10 vector Weyls admit the symmetric54 Yukawa; E6 costs zero.
    n=10;E=np.linalg.qr(np.column_stack([np.ones(n),np.eye(n)[:,:n-1]]))[0][:,1:]
    basis=[np.diag(v) for v in E.T]
    for i in range(n):
        for j in range(i+1,n):
            X=np.zeros((n,n));X[i,j]=X[j,i]=1/np.sqrt(2);basis.append(X)
    Ys=[np.kron(np.eye(2),X) for X in basis];W=sum(Y.T@Y for Y in Ys);Y=Ys[12]
    beta=(W@Y+Y@W)/2+sum(2*Z@Y.T@Z+Z*np.trace(Z.T@Y) for Z in Ys)
    assert np.max(abs(beta-s.Rational(41,5).__float__()*Y))<1e-12
    x,y=s.symbols('x y');screen=[];candidates=[]
    for total in range(1,9):
        for active in range(1,total+1):
            b=8-s.Rational(2,3)*total;R=(27-b)/(s.Rational(31,5)+active);B=120-2*b-4*active*R
            f=124*x*x+s.Rational(224,5)*x*y+s.Rational(159,50)*y*y-B*x+18
            g=24*x*y+s.Rational(127,5)*y*y-B*y+60-4*active*R*R
            rr=s.resultant(f,g,x);count=int(s.Poly(rr,y).count_roots(-s.oo,s.oo))
            screen.append({'total_vector_Weyls':total,'active_Yukawas':active,'real_resultant_roots':count})
            if total==8 and active==2:
                for yy in s.nroots(rr,n=35,maxsteps=200):
                    if abs(s.im(yy))>1e-20:continue
                    yf=float(s.re(yy));xf=float((B*yy-s.Rational(127,5)*yy*yy-60+4*active*R*R)/(24*yy))
                    jac=np.array([[248*xf+44.8*yf-float(B),44.8*xf+6.36*yf],[24*yf,24*xf+50.8*yf-float(B)]])
                    minimum=xf+(yf/10 if yf>=0 else yf*73/90)
                    assert minimum>0
                    candidates.append({'lambda1_over_g_squared':xf,'lambda2_over_g_squared':yf,'Yukawa_squared_over_g_squared':float(R),'quartic_boundedness_margin':minimum,'quartic_ratio_stability_eigenvalues':np.linalg.eigvals(jac).real.tolist()})
    assert all(r['real_resultant_roots']==0 for r in screen if r['total_vector_Weyls']<8)
    assert any(max(c['quartic_ratio_stability_eigenvalues'])<0 for c in candidates)
    return {'new_fields':'Eight E6/family-singlet Weyls in the realSO10 vector10; two have y chi^T S chi/2, six have zeroS Yukawa. RealSO representations are anomaly-free. Masses and protection of spectator zeros require additional model input.',
      'beta_gauge':'bSO10=8-2*8/3=8/3; E6 beta with oneH remains8.',
      'Yukawa_beta':'16pi² beta_y=y[(31/5+nactive)y²-27gSO²]. Direct54-tensor contraction for nactive2 gives41/5.',
      'quartic_changes':'Add4*nactive*y²*lambda1,2 from scalar anomalous dimension and -4*nactive*y4 to beta_lambda2. These remove the no-Yukawa hypothesis of the Riccati barrier.',
      'screen':screen,'bounded_fixed_rays':candidates,
      'scope':'A constructed one-loop completely-AF gauge/Yukawa/two-quartic ray for the isolatedS/fermion block; the lower-lambda1 ray has two UV-attractive quartic ratios and a tuned Yukawa critical ratio. Six uncoupled spectator Weyls slow the gauge running. The full frame/family/mediator portal system is NOT solved, and added charged matter must be made phenomenologically acceptable.'}

def payload():
    _,d,_=I.tensors();rng=np.random.default_rng(11303);v=rng.normal(size=(27,2));F=rng.normal(size=(2,3,3));F=(F+F.transpose(0,2,1))/2
    mu=np.array([[1.,.2],[-.1,2.]]);M=7.;Y=.4
    source=np.einsum('ab,ia,bjk->ijk',mu,v,F);h=-source/M**2
    got=Y*np.einsum('abc,cij->aibj',d,h).reshape(81,81)
    target=-Y/M**2*sum((mu[a,b]*np.kron(np.einsum('abc,c->ab',d,v[:,a]),F[b]) for a in range(2) for b in range(2)),start=np.zeros((81,81)))
    assert np.max(abs(got-target))<1e-12
    return {'status':'PASS','result_scope':'PASS_TWO_CHANNEL_SCALAR_SOURCE_AND_PORTAL_INDEPENDENT_54_UV_BARRIER',
      'single_scalar_action':'H(27,bar6): V=||M H+J/M||², J=sum_ab mu_ab v_a F_bdagger; Yukawa y d(psi,psi,H). Here mu has mass dimension1. Eliminating H=-J/M² gives C_ab=-y mu_ab/M², whose rank can be2. The source square includes required cross portals.',
      'channel_rank':int(np.linalg.matrix_rank(mu)),'signed_E6_matching_error':float(np.max(abs(got-target))),
      'one_loop_beta':'16pi² beta_y=y(25y²-52gE²-8gF²), from324 real scalar Yukawa components on81 Weyls. One scalar retains bE6=8,bSO10=8; family gauging remains non-AF.',
      'tensor_controls':scalar_yukawa(.7),'UV_barrier':uv_barrier(),'alternative_adjoint_encoding':adjoint_reencoding(),'Yukawa_matter_escape':fermion_escape(),
      'prior_owners':['analysis/w33_pass11293_yukawa_uv_budget.py','analysis/w33_pass11298_coupled_running_closure.py'],
      'primary_sources':['https://arxiv.org/abs/hep-ph/0211440','https://arxiv.org/abs/1705.00751']}
if __name__=='__main__':
    out=payload();(ROOT/'data/w33_pass11303_scalar_source_and_uv_barrier.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
