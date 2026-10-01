#!/usr/bin/env python3
"""Numerical full D-flat normalization, signed mass pencil, exact design sum rules.

Prior: Pass11260/61 Cartan, 11269 SIC/MUB dictionary; dated global selector.
Numerical field normalization is kept separate from exact design identities
and outward interval signs for the declared mirror Coleman-Weinberg function.
"""
from pathlib import Path
import sys,json,argparse,math
import numpy as np
import sympy as s
from scipy.linalg import qr,expm
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11261_g26_cartan_coordinate_map as P
import w33_pass11269_g26_qutrit_dictionary as G
from w33_20261001_isovolume_dirac_gravity_response import compare_certificate
OUT=ROOT/'data/w33_20261001_dflat_cartan_mass_loop.json'

def generators():
    B=np.load(ROOT/'artifacts/e6_27rep_basis_export/E6_basis_78.npy').astype(complex)
    hs=np.array([(b+b.conj().T)/2 for b in B]+[(b-b.conj().T)/(2j) for b in B])
    real=np.vstack([hs.real.reshape(156,-1).T,hs.imag.reshape(156,-1).T])
    Q,R,_=qr(real,mode='economic',pivoting=True)
    rank=sum(abs(np.diag(R))>1e-10);assert rank==78
    H=(Q[:729,:rank].T+1j*Q[729:,:rank].T).reshape(rank,27,27)
    bs=B.reshape(78,-1).T
    err=float(np.max(abs(H.reshape(rank,-1).T-bs@np.linalg.lstsq(bs,H.reshape(rank,-1).T,rcond=1e-12)[0])))
    assert err<1e-10
    h=[]
    for i in range(3):
        for j in range(i+1,3):
            a=np.zeros((3,3),complex);a[i,j]=a[j,i]=1/np.sqrt(2);h.append(a)
            a=np.zeros((3,3),complex);a[i,j]=-1j/np.sqrt(2);a[j,i]=1j/np.sqrt(2);h.append(a)
    h +=[np.diag([1,-1,0])/np.sqrt(2),np.diag([1,1,-2])/np.sqrt(6)]
    return H,np.array(h),err

def balance(K,H,h):
    L=np.eye(27,dtype=complex);C=np.eye(3,dtype=complex)
    V=np.array(K[2],complex).reshape(27,3);history=[]
    def moments(v):
        return np.r_[np.einsum('aij,ji->a',H,v@v.conj().T).real,
                     np.einsum('aij,ji->a',h,v.conj().T@v).real]
    for it in range(30):
        m=moments(V);norm=float(np.linalg.norm(V)**2);err=float(np.linalg.norm(m))
        history.append({'iteration':it,'norm_squared':norm,'moment_norm':err})
        if err<1e-9:break
        T=np.r_[np.einsum('aij,jk->aik',H,V),np.einsum('ij,ajk->aik',V,h)].reshape(86,-1)
        Hess=4*np.real(T.conj()@T.T)
        direction=-np.linalg.pinv(Hess,rcond=1e-10)@(2*m)
        AE=np.einsum('a,aij->ij',direction[:78],H)
        AS=np.einsum('a,aij->ij',direction[78:],h);eta=1.
        for _ in range(20):
            ell=expm(eta*AE);cc=expm(eta*AS);nv=ell@V@cc
            if np.linalg.norm(nv)**2<norm+1e-12 and np.linalg.norm(moments(nv))<err:break
            eta*=.5
        else:raise AssertionError('line search failed')
        V=nv;L=ell@L;C=C@cc
    assert np.linalg.norm(moments(V))<1e-9
    Qnew=np.column_stack([(L@np.array(k,complex).reshape(27,3)@C).ravel() for k in K])
    allQ=Qnew.reshape(27,3,3);cross=[]
    for gen in H:cross.append(np.max(abs(np.einsum('ika,ij,jkb->ab',allQ.conj(),gen,allQ))))
    for gen in h:cross.append(np.max(abs(np.einsum('ika,kl,ilb->ab',allQ.conj(),gen,allQ))))
    err=float(max(cross));assert err<1e-9
    gram=Qnew.conj().T@Qnew;c=float(gram[2,2].real)
    target=c*np.diag([2**(2/3),2**(4/3),1])
    gerr=float(np.max(abs(gram-target)));assert gerr<1e-9
    return L,C,Qnew,history,err,gerr

def exact_designs():
    W=s.Symbol('W');mod=s.Poly(W*W+W+1,W)
    red=lambda f:s.rem(s.Poly(s.expand(f),W),mod).as_expr()
    sic=[(1,-W**k,0) for k in range(3)]+[(1,0,-W**k) for k in range(3)]+[(0,1,-W**k) for k in range(3)]
    mub=[(1,0,0),(0,1,0),(0,0,1)]+[(1,W**a,W**b) for a in range(3) for b in range(3)]
    swap=s.zeros(9)
    for i in range(3):
        for j in range(3):swap[3*i+j,3*j+i]=1
    for rays,count in [(sic,3),(mub,4)]:
        a=s.zeros(3);b=s.zeros(9)
        for ray in rays:
            v=s.Matrix(ray);vc=v.applyfunc(lambda x:x.subs(W,W**2))
            P0=(v*vc.T/sum(int(x!=0) for x in ray)).applyfunc(red)
            a+=P0;b+=s.kronecker_product(P0,P0)
        assert a.applyfunc(red)==count*s.eye(3)
        assert b.applyfunc(red)==s.Rational(count,4)*(s.eye(9)+swap)
    return {'first_moments':'SIC=3I, MUB=4I',
      'second_moments':'SIC=(3/4)(I+Swap), MUB=I+Swap',
      'mirror_mass_squared':'MUB: |<b|v>|^2 twice each; SIC: (2/3)|<s|v>|^2 six times each, plus 3 zeros (unit radius)',
      'trace_mass_squared':20,'trace_mass_fourth':8,
      'scale_independence':'For this mass law, the log(mu) and subtraction terms in sum m^4[log(m^2/mu^2)-3/2] are angular constants at fixed radius; quadratic/quartic mass moments do not select an angle.'}

def loop_saddle():
    iv=mp.iv;iv.dps=60
    z=iv.mpc(iv.cos(2*iv.pi/9),iv.sin(2*iv.pi/9))
    w=iv.mpc(iv.cos(2*iv.pi/3),iv.sin(2*iv.pi/3))
    conj=lambda x:iv.mpc(x.real,-x.imag)
    dot=lambda a,b:sum(conj(x)*y for x,y in zip(a,b))
    v=[1/iv.sqrt(3),z/iv.sqrt(3),conj(z)/iv.sqrt(3)]
    e1=[1/iv.sqrt(2),-z/iv.sqrt(2),0]
    e2=[1/iv.sqrt(6),z/iv.sqrt(6),-2*conj(z)/iv.sqrt(6)]
    E=[e1,[1j*x for x in e1],e2,[1j*x for x in e2]]
    sic=[(1,-w**k,0) for k in range(3)]+[(1,0,-w**k) for k in range(3)]+[(0,1,-w**k) for k in range(3)]
    mub=[(1,0,0),(0,1,0),(0,0,1)]+[(1,w**a,w**b) for a in range(3) for b in range(3)]
    H=[[iv.mpf(0) for _ in range(4)] for _ in range(4)]
    for rays,k,weight in [(sic,iv.mpf(2)/3,6),(mub,iv.mpf(1),2)]:
        for ray in rays:
            den=iv.sqrt(sum(int(x!=0) for x in ray));n=[x/den for x in ray]
            amp=dot(n,v);p=(conj(amp)*amp).real
            A=[dot(n,e) for e in E];dp=[2*(conj(amp)*a).real for a in A]
            mass=k*p;fp=mass*(2*iv.log(mass)+1);fpp=2*iv.log(mass)+3
            for i in range(4):
                for j in range(4):
                    dd=2*(conj(A[i])*A[j]).real-2*p*int(i==j)
                    H[i][j]+=weight*(fpp*k*k*dp[i]*dp[j]+fp*k*dd)
    def bounds(x):
        return [math.nextafter(float(x.a),-math.inf),math.nextafter(float(x.b),math.inf)]
    hn=np.array([[sum(bounds(x))/2 for x in row] for row in H])
    eig,U=np.linalg.eigh(hn);witness=[]
    # Fixed rational witnesses avoid basis drift inside degenerate eigenspaces.
    for direction in ([197,0,0,-1000],[0,-197,1000,0]):
        d=np.array(direction,dtype=int)
        ray=sum(H[i][j]*int(d[i])*int(d[j])/1000000 for i in range(4) for j in range(4))
        b=bounds(ray);witness.append({'direction_integer_thousandths':d.tolist(),'Rayleigh_bounds':b})
    assert witness[0]['Rayleigh_bounds'][1]<0 and witness[1]['Rayleigh_bounds'][0]>0
    return {'function':'F=sum_m multiplicity*m^4*log(m^2) for the declared mirror mass law',
      'Hessian_eigenvalues_numerical':eig.tolist(),
      'outward_interval_precision_digits':60,'sign_witnesses':witness,
      'conclusion':'T is a saddle for either overall bosonic or fermionic sign; this simplest shared mass-loop spectrum does not select the polynomial T minimum',
      'scope':'Interval sign proof is for the declared mirror loop function. Its identification with the normalized 81-field mass pencil is checked numerically below, not symbolically certified.'}

def selfadjoint_dirac(A):
    zero=np.zeros_like(A)
    return np.block([[zero,A],[A.conj().T,zero]])

def payload():
    K,_,_=P.slice_matrices();H,h,closure=generators()
    L,C,Q,history,cross,gerr=balance(K,H,h)
    old=P._load(ROOT/'analysis/w33_pass11218_e8_cubic_cartan_kernel.py','mass_loop_old')
    records=old.load_parent().build_records()
    def op(v):
        A=np.zeros((81,81),complex)
        for u,x,o,c in records:A[o,x]+=c*v[u]
        return A
    d=np.zeros((27,27,27))
    for u,x,o,c in records:
        if u%3==0 and x%3==1 and o%3==2:d[u//3,x//3,o//3]=c
    preservation=float(np.max(abs(np.einsum('abc,ai,bj,ck->ijk',d,L,L,L,optimize=True)-d)))
    assert preservation<1e-9 and abs(np.linalg.det(C)-1)<1e-9
    scale=float(np.vdot(Q[:,2],Q[:,2]).real)
    a=2**(1/3);t=np.exp(1j*np.pi/9);w=t**6
    coord=np.array([[a,a*a,1],[t**8*a,t**8*w*a*a,t**8*w*w],[-t*a,-t*w*w*a*a,-t*w]])
    S=np.array(G.num(G.SIC));S/=np.linalg.norm(S,axis=1)[:,None]
    B=np.array(G.num(G.STAB));B/=np.linalg.norm(B,axis=1)[:,None]
    samples=[np.array(x,complex) for x in [[1,0,0],[0,1,0],[0,0,1],[1,2,3],[1+.5j,-.4j,2-.2j]]]
    serr=0.;cov=0.;g81=np.kron(L,C.T);inv=np.linalg.inv(g81)
    for q in samples:
        M=op(Q@q);raw=op(np.array(K,complex).T@q)
        cov=max(cov,float(np.max(abs(M-inv.T@raw@inv))))
        val=np.sort(np.linalg.svd(M,compute_uv=False)**2);v=coord@q
        pred=np.sort(np.r_[np.repeat(scale/3*abs(B.conj()@v)**2,2),np.repeat(2*scale/9*abs(S.conj()@v)**2,6),[0,0,0]])
        serr=max(serr,float(np.max(abs(val-pred))))
    assert max(serr,cov)<1e-8
    val=np.linalg.svd(op(Q[:,2]/np.sqrt(scale)),compute_uv=False);groups=[]
    for x in val:
        if x<1e-9:continue
        if not groups or abs(x-groups[-1]['singular_value'])>1e-8:groups.append({'singular_value':float(x),'multiplicity':1})
        else:groups[-1]['multiplicity']+=1
    # A gauge-invariant nonholomorphic dimension-five mass spurion.
    phi=Q[:,2]/np.sqrt(scale);base=op(phi);B5=np.outer(phi.conj(),phi.conj());eps=.001
    cross5=float(np.max(abs(base.conj().T@B5)))
    gram5=float(np.max(abs((base+eps*B5).conj().T@(base+eps*B5)
                            -base.conj().T@base-eps**2*np.outer(phi,phi.conj()))))
    oldvals=np.linalg.svd(base,compute_uv=False);newvals=np.linalg.svd(base+eps*B5,compute_uv=False)
    assert np.linalg.matrix_rank(base,tol=1e-9)==78
    assert np.linalg.matrix_rank(base+eps*B5,tol=1e-9)==79
    assert max(cross5,gram5)<1e-10
    heavy_error=float(np.max(abs(oldvals[:78]-newvals[:78])))
    assert heavy_error<1e-10 and abs(newvals[78]-eps)<1e-10
    DF=selfadjoint_dirac(base+eps*B5);DF2=DF@DF
    F2=float(np.trace(DF2).real);F4=float(np.trace(DF2@DF2).real)
    assert max(abs(F2-(40+2*eps**2)),abs(F4-(16+2*eps**4)))<1e-10
    return {'schema':'w33.20261001.dflat-cartan-mass-loop.v1',
      'status':'PASS_NUMERICAL_FULL_DFLAT_CARTAN_WITNESS_EXACT_DESIGN_MOMENTS_AND_INTERVAL_LOOP_SADDLE',
      'normalization':{'assumption':'canonical Tr(V^dag V) positive kinetic norm; compact E6 x SU3 moment maps',
        'method':'noncompact Newton descent of norm, Hessian 4 Re<Ta V,Tb V>; pseudoinverse and line search',
        'history':history,'all_Cartan_cross_moment_max':cross,'metric_match_error':gerr,
        'balanced_Gram':'c*diag(2^(2/3),2^(4/3),1)', 'c_numerical':scale,
        'unitary_G26_pullback_relation':'balanced ambient Gram = (c/3) times the prior exact standard-unitary pullback',
        'Hermitian_E6_span_closure_error':closure,'E6_cubic_preservation_error':preservation,
        'E6_transform_real':L.real.tolist(),'E6_transform_imag':L.imag.tolist(),
        'SL3_transform_real':C.real.tolist(),'SL3_transform_imag':C.imag.tolist(),
        'scope':'floating witness, not an algebraic/interval-certified exact gauge transformation or a full F-flat physical vacuum'},
      'mass_interface':{'named_map':'L_Y = y chi_i^alpha epsilon_alpha_beta M(Phi)_ij psi_j^beta + h.c., M from signed d_E6 tensor epsilon_SL3',
        'species':'two distinct chiral fermion multiplets with canonical kinetic terms; a complete anomaly-free matter sector is additional',
        'alternating_tensor_boundary':'M is skew. The associated single-commuting-field scalar cubic vanishes, so this bracket is not a scalar Hessian or a one-field superpotential.',
        'sample_complex_directions_checked':len(samples),'signed_operator_covariance_error':cov,'mirror_spectrum_error':serr,
        'unit_radius_T_singular_value_groups':groups,'kernel_dimension':3,
        'scale':'masses = |y| r times the dimensionless singular values; y and r are not predicted'},
      'exact_mirror_model':exact_designs(),'loop_test':loop_saddle(),
      'mass_to_gravity_Dirac':{
        'map':'D_F = [[0,A],[A^dag,0]], A=y M(Phi)+(c5/Lambda) B5',
        'finite_Hilbert_dimension':162,
        'unit_radius_control_moments':[F2,F4],
        'self_adjoint_residual':float(np.max(abs(DF-DF.conj().T))),
        'moments_for_declared_mirror_model':'F2=40*|y|^2*r^2 + 2*|c5|^2*r^4/Lambda^2; F4=16*|y|^4*r^4 + 2*|c5|^4*r^8/Lambda^4',
        'heat_factor':'Theta_F(t)=2 Tr exp(-t A^dag A), a positive finite input to the isovolume Dirac response',
        'boundary':'162 is the specified self-adjoint doubling of the 81x81 mass map, not the 40/480 finite fixtures or an observed particle count. The EFT cutoff and gravitational heat cutoff are not identified.'},
      'dimension_five_light_mass':{
        'operator':'(c5/Lambda) (Phi^dag chi)^alpha epsilon_alpha_beta (Phi^dag psi)^beta + h.c.',
        'matrix':'B5=conj(Phi) conj(Phi)^T',
        'exact_identity':'M Phi=0 and M^T=-M imply M^dag B5=B5^dag M=0 and B5^dag B5=||Phi||^2 Phi Phi^dag',
        'consequence':'all old nonzero singular masses unchanged; one exact new mass |c5| r^2/Lambda; generic rank 78 -> 79, two zeros remain',
        'gauge_scope':'compact E6 x SU3 invariance follows from the invariant inner products Phi^dag chi and Phi^dag psi',
        'EFT_scope':'nonholomorphic nonsupersymmetric dimension-five operator; no holomorphic superpotential or orbifold selection-rule permission is assumed',
        'hierarchy':'relative to a heavy |y| r sigma, the light/heavy ratio is (|c5|/|y|)*(r/Lambda)/sigma; coefficients and scales are inputs',
        'angular_loop_effect':'the additional light mass is radial and cannot repair the T angular loop saddle',
        'control_epsilon':eps,'rank_before_after':[78,79],
        'cross_term_error':cross5,'squared_mass_identity_error':gram5,
        'old_heavy_singular_value_error':heavy_error,'new_light_singular_value':float(newvals[78])},
      'boundary':['The normalized field and spectral matching are numerical witnesses.','The moment sum rules and interval saddle signs concern the declared mirror mass model.','No observed particle assignments, Yukawa hierarchy, mixing, Newton scale, cosmological constant or complete TOE is derived.'],
      'prior_owners':['analysis/w33_pass11260_exact_semisimple_cartan.py','analysis/w33_pass11261_g26_cartan_coordinate_map.py','analysis/w33_pass11269_g26_qutrit_dictionary.py','analysis/w33_pass11267_tm1_selector.py','analysis/w33_20261001_g26_global_polynomial_vacuum.py'],
      'external_sources':['https://arxiv.org/abs/hep-th/9506098','https://doi.org/10.1103/PhysRevD.7.1888'],
      'checks':{'numerical_Dflat_plane':True,'numerical_kinetic_compatibility':True,'numerical_mass_pencil_control':True,'exact_design_sum_rules':True,'interval_saddle_signs':True}}

def main():
    a=argparse.ArgumentParser();a.add_argument('--check',action='store_true');args=a.parse_args()
    p=payload();txt=json.dumps(p,indent=2,sort_keys=True)+'\n'
    if args.check:assert compare_certificate(p,json.loads(OUT.read_text()))
    else:OUT.write_text(txt)
    print(p['status'])
if __name__=='__main__':main()
