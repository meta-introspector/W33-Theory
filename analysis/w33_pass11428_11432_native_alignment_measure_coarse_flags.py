"""Five scoped continuations of11423-11427. Prior11326 owns engineered noncommuting flavour;11384 owns native CP. Generic mechanisms retain ownership.
Native two-field stiff-orbit alignment; local overlap measure patch; curved
normal-stationary coarse action; actual threshold forces; flagged native gates.
"""
import json,hashlib
from pathlib import Path
from itertools import product,combinations
import numpy as np
from scipy.linalg import expm
from scipy.optimize import root_scalar,minimize
import mpmath as mp
import w33_pass11423_11427_flux_index_curved_uv_noise as P
from w33_pass11384_11388_native_dynamics import native_tensors,compact_generators
from w33_pass11408_11412_clock_spurion_seam_ticks import unitary_calibration
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/w33_pass11428_11432_native_alignment_measure_coarse_flags.json'
def previous():return json.loads((ROOT/'data/w33_pass11423_11427_flux_index_curved_uv_noise.json').read_text())

def cross_pin(phi,psi):
    scale=np.trace(phi.T@phi.conj()).real/3
    C=phi.T@psi.conj()/scale
    U=2*np.eye(3)+(C+C.conj().T)/2;B=np.eye(3)+C@C.conj().T;W=U@B-B@U
    return C,U,B,float(np.imag(np.trace(W@W@W)))

def native_alignment():
    old=json.loads((ROOT/'data/w33_pass11384_11388_native_dynamics.json').read_text())['sections']['11384_native_cp']
    phi=(np.array(old['full_field_real'])+1j*np.array(old['full_field_imag'])).reshape(27,3)
    basis,(_,d,_)=native_tensors();H=compact_generators(basis);he=H[:78,::3,::3]
    import w33_20261001_degree18_phase_completion as D
    Eint=np.array(D.I.C.plane(),dtype=int).reshape(27,3,3)
    gram_coeff=np.einsum('aix,ajy->ijxy',Eint,Eint)
    assert all(np.array_equal(gram_coeff[:,:,x,y],np.eye(3,dtype=int)*gram_coeff[0,0,x,y])for x,y in product(range(3),repeat=2))
    tc=np.einsum('abc,aix,bjy,ckz->ijkxyz',d,Eint,Eint,Eint,optimize=True)
    assert np.array_equal(tc,np.roll(tc,(1,1,1),axis=(0,1,2)))
    assert all(not np.any(tc[i,j,k])for i,j,k in product(range(3),repeat=3)if(i+j+k)%3)
    T=np.einsum('abc,ai,bj,ck->ijk',d,phi,phi,phi)
    G=phi.T@phi.conj();K=np.einsum('ikl,jkl->ij',T,T.conj());L=np.einsum('ikl,km,jml->ij',T,K.conj(),T.conj())
    X=np.roll(np.eye(3),1,axis=0);Z=np.diag(np.exp(2j*np.pi*np.arange(3)/3))
    sym=max(np.linalg.norm(np.einsum('ia,jb,kc,abc->ijk',F,F,F,T)-T)for F in [X,Z])
    assert sym<1e-12
    scalar=[float(np.linalg.norm(A-np.trace(A)/3*np.eye(3)))for A in [G,K,L]];assert max(scalar)<1e-12
    rng=np.random.default_rng(11428);h=np.einsum('a,aij->ij',rng.normal(size=78),he);psi=expm(2j*h)@phi
    chi=np.sqrt(7)/4;scale=np.trace(G).real/3
    def derivative(psi):
        C,U,B,q=cross_pin(phi,psi);dc=-1j*np.einsum('ai,haj->hij',phi,np.einsum('hab,bj->haj',he,psi).conj())/scale
        du=(dc+dc.conj().transpose(0,2,1))/2;db=dc@C.conj().T+C@dc.conj().transpose(0,2,1);W=U@B-B@U
        dw=du@B-B@du+U@db-db@U
        return q,3*np.imag(np.einsum('ij,hji->h',W@W,dw))
    history=[]
    for it in range(180):
        q,grad=derivative(psi);history.append(q);norm=np.linalg.norm(grad)
        if norm<1e-9:break
        generator=np.einsum('a,aij->ij',grad/norm,he);step=.4;accepted=False
        for _ in range(18):
            trial=expm(1j*step*generator)@psi
            if cross_pin(phi,trial)[3]>q+1e-4*step*norm:psi=trial;accepted=True;break
            step*=.5
        if not accepted:break
    C,U,B,q=cross_pin(phi,psi);assert q>0 and min(np.linalg.eigvalsh(U))>=1-1e-10 and min(np.linalg.eigvalsh(B))>=1-1e-10
    moments=max(np.max(abs(np.einsum('i,aij,j->a',v.conj(),H,v)))for v in [phi.ravel(),psi.ravel()]);assert moments<1e-11
    V=expm(1j*np.array([[.1,.2+.1j,.3],[.2-.1j,-.2,.1j],[.3,-.1j,.1]]));E=expm(.2j*h)
    Cp,Up,Bp,qp=cross_pin(E@phi@V.T,E@psi@V.T)
    covariance=max(np.linalg.norm(Ap-V@A@V.conj().T)for Ap,A in [(Cp,C),(Up,U),(Bp,B)]);assert covariance<1e-11
    assert abs(q-cross_pin(phi.conjugate(),psi.conjugate())[3]*-1)<1e-10
    r=np.linalg.inv(np.eye(80)-.08*P.load()['A'])[P.previous()['path_flavour']['endpoints'],0]
    Yu=np.diag(r)@U@np.diag(r);Yd=np.diag(r)@B@np.diag(r);W=Yu@Yu@Yd@Yd-Yd@Yd@Yu@Yu;cp=float(np.imag(np.trace(W@W@W)))
    assert abs(cp)>0 and np.linalg.matrix_rank(Yd,tol=1e-16)==3
    from w33_pass11403_11407_native_bimodule_collective_gates import gellmann
    slots=P.prior()['fermion_bimodule']['native_slot_frames'];frames=[P.read(x)for x in slots['colour_qutrit_seeds']];carrier=np.hstack(frames);native=carrier@np.kron(Yd,np.eye(3))@carrier.conj().T
    colour_residual=max(np.linalg.norm(native@g-g@native)for g in [sum(F@t@F.conj().T for F in frames)for t in gellmann()]);assert colour_residual<1e-11
    return dict(status='PASS',single_field_scalar_residuals=scalar,cubic_clock_invariance=sym,integer_Cartan_embedding=Eint.tolist(),integer_cubic_coefficients=tc.tolist(),exact_Cartan_Gram_and_XZ_cubic_checks=True,
        single_field_boundary='The named algebra generated from the native cubic T and its invariant contractions has irreducible X/Z symmetry, hence its Hermitian family endomorphisms are scalar. This is not a classification of every native covariant.',
        phi=P.enc(phi),psi=P.enc(psi),cross_matrix=P.enc(C),pin_up=P.enc(U),pin_down=P.enc(B),CP_flux=q,chi=chi,
        action='V_align=-kappa chi Im tr[U,B]^3; C=Phi^T Psi*/g0, g0=tr(Phi^T Phi*)/3; U=2I+Re_H C, B=I+CCdagger',
        regime='Two copies of the11384 native vacuum orbit, with stiff radial/normal fields; kappa>0 supplied. Phi and Psi transform under the same compact E6 x SU3, not independent gauge groups.',
        orbit_steps=len(history),flux_history=history,D_flat_residual=float(moments),covariance_residual=float(covariance),pin_eigenvalues=[np.linalg.eigvalsh(A).tolist()for A in [U,B]],matched_CP_cubic=cp,actual_colour_Ward_residual=float(colour_residual),
        existence='Compact relative E6 orbit and a positive CP-flux witness imply a global alignment minimum has nonzero flux; no unique minimizer, full324-field Hessian or physical coefficients are asserted.',
        scope='A constructed two-field gauge-covariant pin map and CP-even stiff-orbit alignment action bypass the tested single-field isotropy. The added field, stiffness, kappa, affine singlet offsets and endpoint embedding are inputs; finite-stiffness vacuum, full flavour fit and UV completion remain open.')

def weak_overlap(q,epsilon,background,tangents,L=3):
    sites=list(product(range(L),repeat=4));ix={x:i for i,x in enumerate(sites)};n=len(sites);gamma,g5=P.gamma_matrices();G=np.kron(np.eye(n),g5)
    W=3*np.eye(4*n,dtype=complex);der=[np.zeros_like(W)for _ in tangents]
    for mu in range(4):
        T=np.zeros((n,n),complex);dTs=[T.copy()for _ in tangents]
        for i,x in enumerate(sites):
            y=list(x);y[mu]=(y[mu]+1)%L;j=ix[tuple(y)];value=np.exp(1j*q*epsilon*background[i,mu])*(-1 if mu==3 and x[mu]==L-1 else 1)
            T[i,j]=value
            for a,f in enumerate(tangents):dTs[a][i,j]=1j*q*f[i,mu]*value
        W-=.5*(np.kron(T,np.eye(4)-gamma[mu])+np.kron(T.conj().T,np.eye(4)+gamma[mu]))
        for a,dT in enumerate(dTs):der[a]-=.5*(np.kron(dT,np.eye(4)-gamma[mu])+np.kron(dT.conj().T,np.eye(4)+gamma[mu]))
    H=G@W;w,V=np.linalg.eigh(H);sign=np.sign(w);Proj=(np.eye(4*n)+(V*sign)@V.conj().T)/2
    F=np.zeros((4*n,4*n));mask=sign[:,None]!=sign[None,:];F[mask]=(sign[:,None]-sign[None,:])[mask]/(w[:,None]-w[None,:])[mask]
    ds=[V@(F*(V.conj().T@(G@a)@V))@V.conj().T/2 for a in der]
    curvature=float((1j*np.trace(Proj@(ds[0]@ds[1]-ds[1]@ds[0]))).real)if len(ds)>=2 else None
    return Proj,ds,curvature,float(min(abs(w))),V[:,sign>0],G

def local_chiral_measure():
    rng=np.random.default_rng(11429);background=rng.normal(size=(81,4));tangents=rng.normal(size=(2,81,4));charges=[1,-4,2,-3,6,0];multiplicities=[6,3,3,2,1,1]
    assert sum(q*m for q,m in zip(charges,multiplicities))==0 and sum(q**3*m for q,m in zip(charges,multiplicities))==0
    rows=[]
    for eps in [.01,.005,.0025]:
        vals=[];gaps=[]
        for q in charges:
            _,_,cur,gap,_,_=weak_overlap(q,eps,background,tangents);vals.append(cur);gaps.append(gap)
        rows.append(dict(epsilon=eps,curvatures=vals,multiplet_curvature=float(np.dot(multiplicities,vals)),minimum_Wilson_gap=min(gaps)))
    # A numerical local frame on a contractible patch: polar transport from
    # a common reference frame. This works for anomalous multiplets too and
    # must not be promoted to a local gauge-invariant global measure current.
    ref,_,_,_,E,G=weak_overlap(1,.005,background,tangents)
    phase=0.;hol=np.eye(E.shape[1],dtype=complex);frames=[]
    for x,y in [(0,0),(.001,0),(.001,.001),(0,.001),(0,0)]:
        Proj,_,_,_,_,_=weak_overlap(1,1.,.005*background+x*tangents[0]+y*tangents[1],[])
        projected=Proj@E;gram=projected.conj().T@projected;w,U=np.linalg.eigh(gram);assert min(w)>.9
        frames.append(projected@(U/np.sqrt(w))@U.conj().T)
    for A,B in zip(frames,frames[1:]):
        u,s,vh=np.linalg.svd(A.conj().T@B);hol=hol@u@vh
    phase=float(np.angle(np.linalg.det(hol)));assert max(np.linalg.norm(F.conj().T@F-np.eye(F.shape[1]))for F in frames)<1e-10
    return dict(status='PASS',L=3,spinor_dimension=324,integer_hypercharges=charges,multiplicities=multiplicities,
        linear_anomaly=0,cubic_anomaly=0,background=background.tolist(),tangents=tangents.tolist(),scans=rows,
        local_frame_rank=E.shape[1],patch_loop_phase=phase,loop_area=1e-6,frame_rule='F(A)=P(A)E [Edagger P(A)E]^-1/2',
        measure_curvature='i Tr P[dP1,dP2]; dP from divided differences of Wilson sign',
        scope='Actual Weyl-bundle curvature and a smooth local section for a supplied anomaly-free hypercharge multiplet. The neutral species has no gauge curvature. Anomaly cancellation alone does not make finite-background measure curvature vanish. A global local gauge-covariant measure current, non-Abelian integrability, topology/mirror selection and physical spacetime remain open.')

def refined_geometry(boundary,centers=None):
    coarse=[tuple(s)for s in previous()['curved_gauge_islands']['coarse_simplices']];fine=[];coords=[];radii=[];normals=[]
    for i,s in enumerate(coarse):
        X=P.local_simplex_coordinates(boundary,s);c=X.mean(axis=0)+(np.zeros(4)if centers is None else centers[i]);r=np.linalg.norm(X-c,axis=1);T=(c-X)/r[:,None];n=np.linalg.svd(T.T)[2][-1]
        if n[np.argmax(abs(n))]<0:n=-n
        coords.append(X);radii.append(r);normals.append(n);fine.extend(tuple([6+i]+[v for v in s if v!=omit])for omit in s)
    return coarse,fine,coords,np.concatenate(radii),normals

def regge_volume_gradient(boundary,rad,lam):
    coarse,fine,_,_,_=refined_geometry(boundary);L=np.zeros((9,9));L[:6,:6]=boundary
    for i,s in enumerate(coarse):
        for j,v in enumerate(s):L[6+i,v]=L[v,6+i]=rad[5*i+j]
    act,d=P.mesh_action(L,fine);volume=0.;vg=np.zeros(15);grad=np.zeros(15)
    for i,s in enumerate(coarse):
        for omit in s:
            vs=[v for v in s if v!=omit];r=L[6+i,vs];G=(r[:,None]**2+r[None,:]**2-L[np.ix_(vs,vs)]**2)/2;vol=np.sqrt(np.linalg.det(G))/24;volume+=vol;g=vol*r*np.linalg.solve(G,np.ones(4))
            for v,a in zip(vs,g):vg[5*i+s.index(v)]+=a
        for j,v in enumerate(s):
            for o in s:
                if o==v:continue
                tri=tuple(sorted([6+i,v,o]));ell=L[6+i,v];opp=L[v,o];adj=L[6+i,o];semi=(ell+opp+adj)/2;area=np.sqrt(semi*(semi-ell)*(semi-opp)*(semi-adj))
                grad[5*i+j]+=ell*(opp**2+adj**2-ell**2)/(8*area)*d[tri]
    return float(act-lam*volume),grad-lam*vg,d,float(volume)

def normal_stationary(boundary,lam,centers=None):
    _,_,coords,flat,normals=refined_geometry(boundary,centers);rad=flat.copy();ts=[]
    for i,n in enumerate(normals):
        def fun(t):
            r=rad.copy();r[5*i:5*i+5]=flat[5*i:5*i+5]+t*n
            return float(n@regge_volume_gradient(boundary,r,lam)[1][5*i:5*i+5])
        t=root_scalar(fun,x0=0.,x1=1e-5,xtol=1e-12).root;ts.append(t);rad[5*i:5*i+5]+=t*n
    act,g,d,vol=regge_volume_gradient(boundary,rad,lam)
    return act,np.array(ts),rad,g,d,vol

def coarse_constraints():
    bd=np.array(previous()['curved_gauge_islands']['boundary_lengths']);rows=[]
    for lam in [0.,1e-4,5e-4,1e-3]:
        act,t,r,g,d,vol=normal_stationary(bd,lam);_,_,coords,flat,normals=refined_geometry(bd)
        ng=np.array([n@g[5*i:5*i+5]for i,n in enumerate(normals)])
        assert np.linalg.norm(ng)<1e-7
        rows.append(dict(Lambda=lam,action=act,normal_shifts=t.tolist(),normal_gradient=ng.tolist(),full_radial_gradient_norm=float(np.linalg.norm(g)),local_curvature=max(abs(v)for tri,v in d.items()if any(x>=6 for x in tri))))
    lam=.001;act,t,r,g,d,vol=normal_stationary(bd,lam)
    def unreduced(z):
        b=bd.copy();b[0,3]+=z[0];b[3,0]=b[0,3];_,_,_,flat,ns=refined_geometry(b)
        rad=flat+np.concatenate([z[i+1]*n for i,n in enumerate(ns)])
        return regge_volume_gradient(b,rad,lam)[0]
    z=np.r_[0,t];H=P.hessian(unreduced,z,h=4e-5);effective=float(H[0,0]-H[0,1:]@np.linalg.solve(H[1:,1:],H[1:,0]))
    def blocked(x,centers=None):
        b=bd.copy();b[0,3]+=x;b[3,0]=b[0,3];return normal_stationary(b,lam,centers)[0]
    h=.001;direct=(blocked(h)+blocked(-h)-2*blocked(0))/h**2
    assert abs(effective-direct)<.05*(1+abs(direct))
    # These are physical boundary responses of a named normal-stationary
    # slice. The remaining twelve center equations are not discarded as gauge.
    shifts=[]
    for j in range(4):
        plus=np.zeros((3,4));plus[0,j]=.01;minus=-plus
        shifts.append((blocked(0,plus)+blocked(0,minus)-2*act)/.01**2)
    assert rows[-1]['local_curvature']>1e-4 and rows[-1]['full_radial_gradient_norm']>1e-7
    return dict(status='PASS',action='S_Regge(boundary+radials)-Lambda sum4-simplex volumes; supplied Lambda',boundary_lengths=bd.tolist(),scans=rows,
        blocked_variables='One boundary length and three stationary radial-normal coordinates; twelve centroid coordinates remain explicit, not assumed gauge.',
        Lambda=lam,stationary_radials=r.tolist(),normal_shifts=t.tolist(),unreduced_Hessian=H.tolist(),boundary_Schur_Hessian=effective,boundary_direct_Hessian=float(direct),first_star_center_second_differences=shifts,
        scope='A constructed normal-stationary coarse action and Schur response on the native-period curved mesh. Nonzero remaining radial equations expose the missing curved-star constraints. This is saddle elimination on a specified slice, not integration over every fine variable, a perfect action, a full stationary gravity solution or Lorentzian spacetime.')

def neutrino_phase_potential(B,S,angles,precision=45):
    with mp.workdps(precision):
        bm=mp.matrix([[mp.mpc(float(z.real),float(z.imag))for z in row]for row in B]);sm=mp.matrix([[mp.mpc(float(z.real),float(z.imag))for z in row]for row in S]);G=mp.diag([mp.exp(1j*mp.mpf(str(angles[0]))),mp.exp(1j*mp.mpf(str(angles[1]))),1]);Y=mp.mpf('.01')/mp.sqrt(2)*G*bm*G.H;N=mp.matrix(6)
        for i in range(3):
            for j in range(3):N[i,j+3]=Y[i,j];N[j+3,i]=Y[i,j];N[i+3,j+3]=sm[i,j]
        w=mp.eighe(N.H*N,eigvals_only=True)
        return -sum(x*x*(mp.log(x)-mp.mpf('1.5'))for x in w)/(32*mp.pi**2)

def full_quark_masses(A,B,phi,phasesQ=None,phasesR=None):
    M=10.;n=len(A);Gq=np.diag(np.exp(1j*(np.zeros(n)if phasesQ is None else phasesQ)));Gr=np.diag(np.exp(1j*(np.zeros(n)if phasesR is None else phasesR)))
    kq=np.kron(M*np.eye(n)-phi*Gq@A@Gq.conj().T,np.eye(3));kr=np.kron(M*np.eye(n)-phi*Gr@A@Gr.conj().T,np.eye(3));left=np.zeros((n*3,3));left[:3]=np.eye(3)/np.sqrt(2);right=np.zeros((3,n*3),complex);right[:,:3]=B
    H=np.block([[kq,left,np.zeros_like(kq)],[np.zeros((3,n*3)),M*np.eye(3),right],[np.zeros_like(kq),np.zeros((n*3,3)),kr]])
    AL=np.zeros((3,len(H)));AR=np.zeros((len(H),3));idx=P.previous()['path_flavour']['endpoints']
    for j,i in enumerate(idx):AL[j,i*3+j]=1;AR[len(H)-240+i*3+j,j]=1
    return np.linalg.svd(np.block([[np.zeros((3,3)),AL],[AR,H]]),compute_uv=False)

@mp.workdps(45)
def radiative_vacuum():
    old=previous();B=P.read(old['pin_vacua']['pin_down']);U=P.read(old['pin_vacua']['pin_up']);S=P.read(P.previous()['majorana_vacuum']['Majorana_source']);A=P.load()['A'];base=neutrino_phase_potential(B,S,[0,0]);unit=mp.mpf('1e-10')
    fun=lambda x:float((neutrino_phase_potential(B,S,x)-base)/unit)
    grid=[(fun([x,y]),[x,y])for x,y in product(np.arange(6)*np.pi/6,repeat=2)];start=min(grid,key=lambda z:z[0])[1]
    opt=minimize(fun,start,method='BFGS',options=dict(gtol=1e-6,maxiter=80));angles=np.mod(opt.x,np.pi);h=mp.mpf('.0002');center=neutrino_phase_potential(B,S,angles);HH=np.zeros((2,2))
    for i in range(2):
        xp=angles.copy();xm=angles.copy();xp[i]+=float(h);xm[i]-=float(h)
        HH[i,i]=float((neutrino_phase_potential(B,S,xp)+neutrino_phase_potential(B,S,xm)-2*center)/(h*h))
    HH[0,1]=HH[1,0]=float((neutrino_phase_potential(B,S,angles+float(h))-neutrino_phase_potential(B,S,angles+float(h)*np.array([1,-1]))-neutrino_phase_potential(B,S,angles+float(h)*np.array([-1,1]))+neutrino_phase_potential(B,S,angles-float(h)))/(4*h*h))
    G=np.diag(np.exp(1j*np.r_[angles,0]));selected=G@B@G.conj().T;phase_tangents=[1j*(E@selected-selected@E)for E in [np.diag([1,0,0]),np.diag([0,1,0])]];metric=np.array([[np.trace(x@y).real for y in phase_tangents]for x in phase_tangents])
    from scipy.linalg import eigh
    physical=eigh(HH,metric,eigvals_only=True);assert min(physical)>0 and float(center-base)<=1e-17
    cp=neutrino_phase_potential(B.conjugate(),S,-angles);assert abs(float(cp-center))<1e-25
    probe=np.diag(np.exp(1j*np.array([.31,-.23,0.])));Sp=probe.conjugate()@S@probe.conj().T
    covariant=neutrino_phase_potential(B,Sp,[.31,-.23]);covariant_residual=abs(float(covariant-base));assert covariant_residual<1e-14
    # Node-gradient link phases, anchored at the pin, leave all fermion
    # singular values invariant. This is a determinant statement, not
    # all-orders protection of anomalous global transformations.
    rng=np.random.default_rng(11431);phases=rng.normal(size=(3,80));phases-=phases[:,0,None];spectral=[]
    for pin,k in [(U,1),(B,2)]:
        m=full_quark_masses(A,pin,.8);shift=full_quark_masses(A,pin,.8,phases[0],phases[k]);spectral.append(float(np.max(abs(m-shift))))
    assert max(spectral)<1e-11
    def threshold(phi):
        masses=[full_quark_masses(A,pin,phi)for pin in [U,B]]
        ferm=-12*sum(np.sum(m**4*(np.log(m*m)-1.5))for m in masses)/(64*np.pi**2)
        rad=2*.1*(3*phi**2-.8**2);phase=2*.1*(phi**2-.8**2)
        scalar=480*sum(x*x*(np.log(x)-1.5)for x in [rad,phase]if x>0)/(64*np.pi**2)
        tree=480*.1*(phi*phi-.8**2)**2
        return dict(phi=phi,fermion=float(ferm),link_scalar=float(scalar),tree=tree,total=float(ferm+scalar+tree),link_phase_mass_squared=phase)
    rows=[threshold(v)for v in [.8,.8005,.801,.802,.81]];step=.0005;force=(-3*rows[0]['total']+4*rows[1]['total']-rows[2]['total'])/(2*step)
    assert force<0
    return dict(status='PASS',phase_angles=angles.tolist(),phase_potential_shift=float(center-base),phase_Hessian=HH.tolist(),phase_kinetic_metric=metric.tolist(),canonical_phase_mass_squared=physical.tolist(),selected_down_pin=P.enc(selected),CP_branch_phase_threshold_residual=abs(float(cp-center)),phase_precision=45,covariant_Majorana_rotation_residual=covariant_residual,
        neutrino_phase_mechanism='Fixed real Majorana source breaks the isolated pin rephasing symmetry in this fermion determinant; quark singular values remain rephasing-invariant. Simultaneously transforming S as G* S Gdagger restores isospectrality; no gauged family Goldstone is lifted.',
        anchored_node_phase_count=237,remaining_link_cycle_phase_count=243,node_phase_profiles=phases.tolist(),node_phase_spectral_residuals=spectral,
        radial_link_threshold_scans=rows,one_sided_radial_force=float(force),
        radial_slice='Common link modulus phi>=.8; full up/down486-state fermion spectra and all480 link radial+480 phase scalar eigenvalues, plus tree link potential. Other scalar/gauge thresholds are constant on this declared slice.',
        scope='A numerically locally stable neutrino-loop lift of the two isolated pin phase directions and explicit joint link-modulus threshold forces.237 anchored node-gradient phases remain isospectral in the stated one-loop determinant, not protected against nonperturbative anomalies. The strong radial force means the old tree background is not a radiatively stationary vacuum. Full995-field loop stability, RG-improved action, counterterm selection and cosmological constant remain open.')

def pauli_class(x,z):
    return P.hamming_syndrome(x)+8*P.hamming_syndrome(z)+64*(x.bit_count()%2)+128*(z.bit_count()%2)

def flagged_circuit(initial=(0,0),faults=None,prep=None,read=None):
    faults={}if faults is None else faults;prep={}if prep is None else prep;read={}if read is None else read
    x,z=initial;out=[];flags=[];gate=0;check=0
    def cnot(c,t):
        nonlocal x,z,gate
        if(x>>c)&1:x^=1<<t
        if(z>>t)&1:z^=1<<c
        b=faults.get(gate,0)
        if b&1:x^=1<<c
        if b&2:z^=1<<c
        if b&4:x^=1<<t
        if b&8:z^=1<<t
        gate+=1
    def prepare(kind,flagged):
        nonlocal x,z
        x&=127;z&=127
        for which in ([0,1]if flagged else[0]):
            if prep.get((check,which),False):
                anc=7+which
                if(kind==0)==(which==0):x^=1<<anc
                else:z^=1<<anc
    for rnd in range(3):
        sy=0;fl=0
        for kind in [0,1]:
            for row in range(3):
                prepare(kind,True);data=[j for j in range(7)if((j+1)>>row)&1];dg=[(j,7)if kind==0 else(7,j)for j in data];fg=(8,7)if kind==0 else(7,8)
                for c,t in [dg[0],fg,dg[1],dg[2],fg,dg[3]]:cnot(c,t)
                a=((x if kind==0 else z)>>7)&1;f=((z if kind==0 else x)>>8)&1
                a^=bool(read.get((check,0),False));f^=bool(read.get((check,1),False));sy|=a<<(row+3*kind);fl|=f<<(row+3*kind);check+=1
        out.append(sy);flags.append(fl)
    extra=-1
    # Under the one-fault contract, a flag already locates that fault.
    # The conditional follow-up has no additional fault in that census.
    # Multi-fault simulation nevertheless includes its actual noisy gates.
    if any(flags):
        extra=0
        for kind in [0,1]:
            for row in range(3):
                prepare(kind,False)
                for data in [j for j in range(7)if((j+1)>>row)&1]:cnot(*((data,7)if kind==0 else(7,data)))
                a=((x if kind==0 else z)>>7)&1;a^=bool(read.get((check,0),False));extra|=a<<(row+3*kind);check+=1
    return x&127,z&127,tuple(out+flags+[extra]),gate

def flag_decoder():
    single={0}|{pauli_class(a<<j,b<<j)for j in range(7)for a,b in [(1,0),(0,1),(1,1)]};allowed={};cases=[]
    def add(result,exact=False):
        x,z,history,_=result;e=pauli_class(x,z);choices={e}if exact else{e^v for v in single};allowed[history]=allowed.get(history,set(range(256)))&choices;cases.append((e,history,exact))
    add(flagged_circuit(),True)
    for j in range(7):
        for a,b in [(1,0),(0,1),(1,1)]:add(flagged_circuit((a<<j,b<<j)),True)
    for g in range(108):
        for b in range(1,16):add(flagged_circuit(faults={g:b}))
    for k in range(18):
        for a in [0,1]:add(flagged_circuit(prep={(k,a):True}));add(flagged_circuit(read={(k,a):True}))
    assert all(allowed.values());decoder={h:min(v)for h,v in allowed.items()}
    assert all((e==decoder[h]if exact else e^decoder[h]in single)for e,h,exact in cases)
    return decoder,single,cases

def native_CNOT():
    old=previous();B=P.read(old['native_noise']['invariant_frame']);Q=B.conj().T@P.read(P.prior()['local_gates']['logical_basis']);D=np.array(P.load()['D'],dtype=float);Hp=B.conj().T@(D[:40].T@D[:40])@B;Hl=B.conj().T@(D[40:].T@D[40:])@B;n=len(Hp)
    a,b=P.prior()['local_gates']['gate_pair'];axes=[B[e].conjugate()for e in [a,b]];proj=[np.outer(f,f.conj())for f in axes];clock=expm(-.02j*Hp)@expm(-.04j*Hl)@expm(-.02j*Hp);H=np.eye(n,dtype=complex)
    for angle,axis,N in zip(P.previous()['counted_noise']['principal_angles'],[0,1,0,1,0],P.previous()['counted_noise']['ticks_per_pulse']):
        kick=np.eye(n)+(np.exp(-.02j*np.copysign(.0125,angle))-1)*proj[axis];H=np.linalg.matrix_power(kick@clock@kick,N)@H
    kick=np.eye(n)+(np.exp(-.02j*.0125)-1)*proj[0];Sz=kick@clock@kick;ang,V,Nz,_,_=unitary_calibration(Sz,Q[:,0],np.pi);Z=(V*np.exp(-1j*Nz*ang))@V.conj().T
    f=np.kron(axes[0],axes[0]);delta=np.outer(f,f.conj());kick=np.eye(n*n)+(np.exp(-.02j*.0125)-1)*delta;Scz=kick@np.kron(clock,clock)@kick;logical=np.kron(Q[:,0],Q[:,0]);ang,V,Ncz,_,_=unitary_calibration(Scz,logical,np.pi);CZ=(V*np.exp(-1j*Ncz*ang))@V.conj().T
    actual=np.kron(np.eye(n),H)@np.kron(Z,Z)@CZ@np.kron(np.eye(n),H);frame=np.kron(Q,Q);target=np.eye(4);target[[2,3]]=target[[3,2]]
    encoded=frame.conj().T@actual@frame;phase=np.vdot(target,encoded);phase/=abs(phase);error=float(np.linalg.norm(actual@frame-phase*frame@target,ord=2));leak=float(1-np.linalg.norm(encoded,'fro')**2/4)
    assert error<.1 and np.linalg.norm(actual.conj().T@actual-np.eye(n*n))<1e-7
    # The extra stochastic channel is native-edge phase noise AFTER the
    # compiled gate; it is not called noise at every internal tick.
    zedge=[np.eye(n)-2*p for p in proj];generators=[np.kron(z,np.eye(n))for z in zedge]+[np.kron(np.eye(n),z)for z in zedge]
    paulis=[]
    for bits in range(16):
        def one(x,z):
            X=np.array([[0,1],[1,0]]);Z=np.diag([1,-1]);return np.linalg.matrix_power(X,x)@np.linalg.matrix_power(Z,z)
        paulis.append(np.kron(one(bool(bits&1),bool(bits&2)),one(bool(bits&4),bool(bits&8))))
    scans=[]
    for q in [0.,1e-5,1e-4]:
        weights=np.zeros(16)
        for mask in range(16):
            state=actual@frame
            for j,g in enumerate(generators):
                if(mask>>j)&1:state=g@state
            probability=q**mask.bit_count()*(1-q)**(4-mask.bit_count());K=frame.conj().T@state@target.conj().T
            weights+=probability*np.array([abs(np.trace(Pa.conj().T@K)/4)**2 for Pa in paulis])
        erasure=max(0.,1-weights.sum());assert min(weights)>=0 and weights.sum()<=1+1e-7
        scans.append(dict(native_boundary_phase_flip_probability=q,Pauli_weights=weights.tolist(),flagged_leakage=erasure,identity_weight=float(weights[0])))
    return dict(status='PASS',invariant_cell_dimension=n,two_cell_dimension=n*n,native_CNOT_unitary=P.enc(actual),logical_frame=P.enc(frame),logical_target=P.enc(target),H_ticks=old['native_noise']['ticks'],Z_ticks=Nz,CZ_ticks=Ncz,CNOT_total_ticks=2*old['native_noise']['ticks']+2*Nz+Ncz,phase_aligned_operator_error=error,average_coherent_leakage=leak,scans=scans,
        compiler='CNOT = - (I tensor H) (Z tensor Z) CZ00 (I tensor H), up to global phases; native Z and CZ00 calibrated at pi',
        scope='Exact finite-penalty native CNOT compilation with actual carrier leakage and a specified gate-boundary edge-noise channel. Pauli twirling and ideal leakage detection are supplied. During-CZ/H/Z tick noise, correlated drift and physical readout remain open.')

def wilson_interval(k,n):
    z=1.96;den=1+z*z/n;mid=(k/n+z*z/(2*n))/den;half=z*np.sqrt(k/n*(1-k/n)/n+z*z/(4*n*n))/den;return [float(max(0,mid-half)),float(min(1,mid+half))]

def flagged_native_correction():
    decoder,single,cases=flag_decoder();native=native_CNOT();rng=np.random.default_rng(11432);rows=[];base=P.minimum_weight_decoder()
    for gate in native['scans']:
        weights=np.array(gate['Pauli_weights']);conditional=weights/weights.sum();cdf=np.cumsum(conditional);cdf[-1]=1.;n=16384;fail=0;branch=0;eps=1e-5
        for _ in range(n):
            draws=np.searchsorted(cdf,rng.random(132));faults={i:int(v)for i,v in enumerate(draws)if v};prep={(i,a):True for i in range(24)for a in [0,1]if rng.random()<eps};read={(i,a):True for i in range(24)for a in [0,1]if rng.random()<eps}
            x,z,h,count=flagged_circuit(faults=faults,prep=prep,read=read);branch+=count>108
            if h in decoder:correction=decoder[h]
            else:
                # Explicit fallback for histories outside the certified
                # single-fault set; not an optimized multi-fault decoder.
                sy=sum((((h[0]>>j)&1)+((h[1]>>j)&1)+((h[2]>>j)&1)>=2)<<j for j in range(6));rx,rz=base[sy];correction=pauli_class(rx,rz)
            fail+=pauli_class(x,z)^correction not in single
        rows.append(dict(native_boundary_phase_flip_probability=gate['native_boundary_phase_flip_probability'],samples=n,uncorrectable_residual_count=fail,uncorrectable_residual_estimate=fail/n,Wilson95=wilson_interval(fail,n),conditional_branch_fraction=branch/n,
                         accepted_Pauli_error_probability=float(1-conditional[0]),worst_case_no_leakage_acceptance=(1-gate['flagged_leakage'])**132,readout_and_preparation_error=eps))
    return dict(status='PASS',native_CNOT=native,decoder={','.join(map(str,h)):c for h,c in decoder.items()},single_fault_cases=len(cases),history_count=len(decoder),baseline_CNOTs=108,conditional_CNOTs=132,data_modules=7,reused_ancillas=2,single_fault_uncorrectable_cases=0,
        repair='Three full flagged rounds, then one unflagged full syndrome round if any flag fired. Under the single-fault contract the conditional branch is clean; multi-fault simulation includes its noisy gates.',
        initial_error_condition='All21 single input Paulis with no circuit faults are restored modulo stabilizers.',
        circuit_fault_condition='Any tested single CNOT Pauli, wrong-basis ancilla preparation or measurement-bit fault with clean input leaves a residual equivalent to at most one data Pauli.',
        scans=rows,leakage_policy='Ideal detection and rejection of any compiled-gate leakage; postselected accepted-Pauli simulations are reported separately from the explicit acceptance lower bound.',
        scope='A certified one-fault flagged circuit in the stated stochastic Pauli/preparation/readout model, attached to an actual noisy native CNOT model. Multi-fault rates are Monte Carlo with intervals and a declared fallback decoder. Rejecting leakage is not deterministic leakage correction. No scalable fault-tolerance threshold or laboratory hardware is proved.')

def produce():
    out=dict(status='PASS',passes='11428-11432',reservation='8e8bfdcb3')
    for name,fn in [('native_alignment',native_alignment),('local_chiral_measure',local_chiral_measure),('coarse_constraints',coarse_constraints),('radiative_vacuum',radiative_vacuum),('flagged_native_correction',flagged_native_correction)]:
        out[name]=fn();print(name,'PASS',flush=True)
    sources=['data/w33_pass11423_11427_flux_index_curved_uv_noise.json','data/w33_pass11384_11388_native_dynamics.json','data/w33_pass11384_native_inputs.json','data/w33_pass11413_11417_pin_regge_vacuum_noise.json','data/w33_pass11403_11407_native_bimodule_collective_gates.json']
    out['source_sha256']={s:hashlib.sha256(json.dumps(json.loads((ROOT/s).read_text()),sort_keys=True,separators=(',',':')).encode()).hexdigest()for s in sources};out['scope']='Five explicit models and controls; no common physically selected action, full TOE, measured parameters or global fault-tolerance threshold.'
    OUT.write_text(json.dumps(out,indent=2)+'\n');return out
if __name__=='__main__':produce()
