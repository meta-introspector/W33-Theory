"""Five finite constructions continuing17b4bb65b, with physical assumptions explicit.

Native E6, flux links, scalar Schur map, Gaussian kernel and reset target retain
ownership in11443--11455. Standard trinification and Regge methods are prior art.
"""
import hashlib
import json
from itertools import combinations, product
from pathlib import Path
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import minimize
import w33_pass11448_11455_subgroup_toron_auxiliary_control as N
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/w33_pass11461_11465_explicit_physical_maps.json'
enc=lambda a: {'real':np.asarray(a).real.tolist(),'imag':np.asarray(a).imag.tolist()}
dec=lambda a: np.array(a['real'])+1j*np.array(a['imag'])


def centralizer_split():
    import w33_pass11438_finite_native_model as F
    K,_=N.P.gauge_fixed_space();K=K[:,::3,::3];E=F.H[:78,::3,::3]
    comm=np.concatenate([np.array([(u@v-v@u).ravel() for v in E]).T for u in K])
    z=null_space(np.r_[comm.real,comm.imag],rcond=1e-9);C=np.einsum('ab,aij->bij',z,E)
    assert len(C)==16
    w,V=np.linalg.eigh(np.einsum('aij,bji->ab',C,C).real);C=np.einsum('ab,aij->bij',V/np.sqrt(w),C)
    flat=C.reshape(16,-1).T
    ads=np.array([np.linalg.lstsq(flat,np.array([(1j*(a@b-b@a)).ravel() for b in C]).T,rcond=None)[0].real for a in C])
    I=np.eye(16);equ=np.concatenate([np.kron(I,a)-np.kron(a.T,I) for a in ads])
    zz=null_space(equ,rcond=1e-8);assert zz.shape[1]==2
    centroid=sum((j+1)*v.reshape(16,16,order='F') for j,v in enumerate(zz.T));centroid=(centroid+centroid.T)/2
    ew,ev=np.linalg.eigh(centroid);assert ew[8]-ew[7]>.1 and np.ptp(ew[:8])<1e-8 and np.ptp(ew[8:])<1e-8
    L=np.einsum('ab,aij->bij',ev[:,:8],C);R=np.einsum('ab,aij->bij',ev[:,8:],C)
    assert max(np.linalg.norm(a@b-b@a) for a in L for b in R)<1e-9
    return K,L,R,ew


def cartan(A):
    g=np.einsum('a,aij->ij',np.arange(1,9),A)
    comm=np.array([(g@a-a@g).ravel() for a in A]).T
    z=null_space(np.r_[comm.real,comm.imag],rcond=1e-8);assert z.shape[1]==2
    return np.einsum('ab,aij->bij',z,A)


def charge_map(K=None,L=None,R=None,ew=None):
    if K is None:K,L,R,ew=centralizer_split()
    CK=sum(a@a for a in K);CL=sum(a@a for a in L);CR=sum(a@a for a in R)
    # Pick the colour-charged sector on which the requested factor acts.
    def sector(C):
        w,V=np.linalg.eigh(CK+2*C);return V[:,(w>1e-8)&(np.real(np.diag(V.conj().T@CK@V))>1e-8)&(np.real(np.diag(V.conj().T@C@V))>1e-8)]
    def assign(A,C,values):
        T=cartan(A);B=sector(C);assert B.shape[1]==9
        w,V=np.linalg.eigh(B.conj().T@(T[0]+np.sqrt(2)*T[1])@B);Q=B@V
        weights=np.array([[np.vdot(Q[:,j],t@Q[:,j]).real for t in T] for j in [0,3,6]])
        coeff=np.linalg.lstsq(weights,np.array(values),rcond=None)[0];assert np.linalg.norm(weights@coeff-values)<1e-9
        return np.einsum('a,aij->ij',coeff,T),T
    YL,TL=assign(L,CL,[1/6,1/6,-1/3]);YR,TR=assign(R,CR,[-2/3,1/3,1/3]);Y=YL+YR
    comm=np.array([(YL@a-a@YL).ravel() for a in L]).T
    z=null_space(np.r_[comm.real,comm.imag],rcond=1e-8);assert z.shape[1]==4
    coeff=np.einsum('aij,ji->a',L,YL).real
    z=z@null_space((coeff@z)[None,:]);W=np.einsum('ab,aij->bij',z,L);assert len(W)==3
    CW=sum(a@a for a in W);W*=np.sqrt(.75/np.linalg.eigvalsh(CW)[-1]);CW=sum(a@a for a in W)
    # A native cubic colour Casimir fixes conjugate orientation, which the
    # quadratic Casimir alone cannot do. The chosen colour-left sector is3.
    Bcolour=sector(CL);Pcolour=Bcolour@Bcolour.conj().T;C3=np.zeros((27,27),complex)
    for a,b,c in product(K,repeat=3):
        d=float(np.trace(Pcolour@a@(b@c+c@b)/2).real);C3+=d*a@b@c
    C3=(C3+C3.conj().T)/2;C3/=np.trace(Pcolour@C3).real/9
    # Simultaneous colour/weak Casimir and hypercharge eigenspaces.
    w,Q=np.linalg.eigh(3*CK+5*CW+Y);rows=[]
    for j in range(27):
        v=Q[:,j];c=float(np.vdot(v,CK@v).real);a=float(np.vdot(v,CW@v).real);y=float(np.vdot(v,Y@v).real)
        key=(int(c>1e-8),int(a>1e-8),int(round(6*y)))
        hit=next((r for r in rows if tuple(r['key'])==key),None)
        if hit:hit['states']+=1
        else:rows.append(dict(key=list(key),states=1,hypercharge=y,colour_orientation=int(round(np.vdot(v,C3@v).real))))
    rows.sort(key=lambda r:tuple(r['key']))
    expected={(1,1,1):6,(1,0,-2):3,(1,0,-4):3,(1,0,2):6,(0,1,3):2,(0,1,-3):4,(0,0,6):1,(0,0,0):2}
    assert {tuple(r['key']):r['states'] for r in rows}==expected
    anomalies=dict(gravity=float(np.trace(Y).real),cubic=float(np.trace(Y@Y@Y).real),colour=float(np.trace(Y@CK).real),weak=float(np.trace(Y@CW).real))
    assert max(abs(v) for v in anomalies.values())<1e-9
    residual=max(np.linalg.norm(Y@a-a@Y) for a in list(K)+list(W))
    T3,_=assign(L,CL,[.5,-.5,0.]);Qem=Y+T3
    x=np.array(N.previous()['preserved_gauge']['coordinates']);phi=(x[:81]+1j*x[81:162]).reshape(27,3);psi=(x[162:243]+1j*x[243:]).reshape(27,3)
    acts=np.array([np.r_[(a@phi).ravel(),(a@psi).ravel()] for a in list(W)+[Y]])
    mass=np.einsum('ai,bi->ab',acts.conj(),acts).real;massw=np.linalg.eigvalsh(mass)
    photon=float(np.linalg.norm(Qem@phi)**2+np.linalg.norm(Qem@psi)**2)
    assert min(massw)>1e-8 and photon>1e-8
    return dict(status='PASS',colour=enc(K),left=enc(L),right=enc(R),weak=enc(W),hypercharge=enc(Y),weak_Casimir=enc(CW),colour_cubic_Casimir=enc(C3),centroid_eigenvalues=ew.tolist(),branching=rows,anomalies=anomalies,commutation_error=float(residual),electric_charge=enc(Qem),weak_Y_gauge_Gram_eigenvalues=massw.tolist(),photon_mass_Gram=photon,scope='Native27 matrix embedding of a chosen SU3c x SU2L x U1Y inside the preserved su3 and its two su3 centralizer ideals. Branching is one SM family, a vectorlike coloured pair, a vectorlike weak pair and two neutral singlets. The cubic colour Casimir distinguishes3 from conjugate3 relative to a supplied orientation. Which ideal is weak, the hypercharge weight assignment, chirality and removal of extra multiplets are supplied choices, not derived by the vacuum; no unique physical embedding claim. All four weak/hypercharge directions and the chosen electromagnetic generator are broken by the preserved candidate; it is not an SM vacuum with a massless photon.')


def little_algebra(x):
    import w33_pass11438_finite_native_model as F
    a=x[:81]+1j*x[81:162];b=x[162:243]+1j*x[243:]
    orbit=np.concatenate([(1j*(F.H@a)).real,(1j*(F.H@a)).imag,(1j*(F.H@b)).real,(1j*(F.H@b)).imag],axis=1).T
    rank=int(np.linalg.matrix_rank(orbit,tol=1e-8));z=np.linalg.svd(orbit,full_matrices=True)[2][rank:].T
    K=np.einsum('ab,aij->bij',z,F.H);flat=K.reshape(len(K),-1).T;ads=[];closure=[]
    for u in K:
        v=np.array([(1j*(u@v-v@u)).ravel() for v in K]).T;c=np.linalg.lstsq(flat,v,rcond=None)[0];ads.append(c);closure.append(np.linalg.norm(v-flat@c))
    C=sum(a.conj().T@a for a in ads);w=np.linalg.eigvalsh(C)
    return dict(dimension=len(K),center_dimension=int(sum(w<1e-9)),family_component_norm=float(np.linalg.norm(z[78:])),closure_error=float(max(closure)),adjoint_Casimir_eigenvalues=w.tolist(),generators=enc(K[:,::3,::3]),representation_Casimir_eigenvalues=np.linalg.eigvalsh(sum(k@k for k in K[:,::3,::3])).tolist(),classification='At the recorded1e-8 cutoff, a center-free compact15-dimensional algebra containing the existing su3 is su4: the only other center-free compact15-dimensional possibility is su2^5, which cannot contain su3.')


def rank_sensitivity(photon_x):
    import w33_pass11438_finite_native_model as F
    old=json.loads(N.P.N.OUT.read_text())['coercive_alignment']['finite_field_witness']
    rows=[]
    for name,x in [('photon',photon_x),('colour',np.array(N.previous()['preserved_gauge']['coordinates'])),('older',np.array(old['coordinates']))]:
        a=x[:81]+1j*x[81:162];b=x[162:243]+1j*x[243:]
        orbit=np.concatenate([(1j*(F.H@a)).real,(1j*(F.H@a)).imag,(1j*(F.H@b)).real,(1j*(F.H@b)).imag],axis=1).T
        sv=np.linalg.svd(orbit,compute_uv=False)
        rows.append(dict(candidate=name,singular_values=sv.tolist(),threshold_ranks=[dict(threshold=t,rank=int(sum(sv>t))) for t in [1e-8,1e-7,1e-6,1e-5]]))
    return dict(rows=rows,scope='Recorded ranks at1e-8 reproduce prior numerical values, but many additional orbit singular values are at solver-noise scale. All three candidates have rank58 at1e-6. No exact limiting stationary stabilizer is certified; possible enhanced symmetry requires a constructed high-precision or symbolic stationary witness. This qualifies physical interpretations without changing the prior certificate values.')


def photon_preserving_candidate(charges,x=None,H=None):
    import w33_pass11438_finite_native_model as F
    from scipy.linalg import block_diag
    K=dec(charges['colour']);Q=dec(charges['electric_charge'])
    w,B=np.linalg.eigh(sum(a@a for a in K)+Q@Q);B=np.kron(B[:,w<1e-9],np.eye(3))
    assert B.shape[1]==15
    E=np.block([[B.real,-B.imag],[B.imag,B.real]]);T=block_diag(E,E)
    with F.native_context():
        if x is None:
            old=np.array(N.previous()['preserved_gauge']['coordinates']);start=T.T@old
            def obj(y):
                v,g=F.fun(T@y);return v,T.T@g
            sol=minimize(obj,start,jac=True,method='L-BFGS-B',options=dict(maxiter=400,gtol=1e-8,ftol=1e-14,maxls=30,maxcor=30));x=T@sol.x
        x=np.asarray(x);v,g=F.fun(x);a=x[:81]+1j*x[81:162];b=x[162:243]+1j*x[243:]
        orbit=np.concatenate([(1j*(F.H@a)).real,(1j*(F.H@a)).imag,(1j*(F.H@b)).real,(1j*(F.H@b)).imag],axis=1).T
        rank=int(np.linalg.matrix_rank(orbit,tol=1e-8))
        if H is None:
            h=2e-5;H=np.column_stack([(F.fun(x+h*d)[1]-F.fun(x-h*d)[1])/(2*h) for d in np.eye(324)]);H=(H+H.T)/2
        Z=np.linalg.svd(orbit,full_matrices=True)[0][:,rank:];eig,V=np.linalg.eigh(Z.T@H@Z);d=Z@V[:,0]
        replay=(F.fun(x+1e-5*d)[1]-F.fun(x-1e-5*d)[1])/2e-5
        difference=float(np.linalg.norm(replay-H@d))
        scan=[dict(amplitude=t,symmetric_energy_change=(F.fun(x+t*d)[0]+F.fun(x-t*d)[0])/2-v) for t in [.001,.003,.01]]
    photon=float(np.linalg.norm(Q@a.reshape(27,3))**2+np.linalg.norm(Q@b.reshape(27,3))**2)
    assert np.linalg.norm(g)<5e-6 and photon<1e-14
    return dict(status='PASS',little_algebra=little_algebra(x),coordinates=x.tolist(),value=float(v),gradient_norm=float(np.linalg.norm(g)),fixed_complex_dimension=B.shape[1],common_gauge_rank=rank,unbroken_dimension=86-rank,photon_mass_Gram=photon,CP_flux=N.P.N.bounded_cross(a.reshape(27,3),b.reshape(27,3))[1],normal_eigenvalues=eig.tolist(),negative_normal_count=int(sum(eig < -1e-3)),lowest_direction_replay_error=difference,lowest_direction_energy_scans=scan,scope='Stationary colour/electric-charge-preserving search in the15-complex-dimensional fixed subspace per field of the supplied324-field coercive action. Full normal Hessian and a smaller-step lowest-mode replay determine whether this candidate is a transverse saddle. It is a searched candidate, not a new selected Higgs potential, observed spectrum or global minimum; prior11271 owns the separate canonical SM-preserving configuration.')


def flux_frame(t2,t3,p2,p3,L=16):
    """Exact Fourier block of the4D overlap operator on16x16x2x2."""
    gamma,g5=N.P.N.M.P.gamma_matrices();u,v=N.P.N.uniform_flux_links(L);n=L*L
    W=(3-np.cos(p2+t2/2)-np.cos(p3+t3/2))*np.eye(4*n,dtype=complex)
    W+=1j*np.kron(np.eye(n),np.sin(p2+t2/2)*gamma[2]+np.sin(p3+t3/2)*gamma[3])
    for mu,link in enumerate([u,v]):
        T=np.zeros((n,n),complex)
        for x,y in product(range(L),repeat=2):
            target=((x+1)%L,y) if mu==0 else (x,(y+1)%L)
            T[x*L+y,target[0]*L+target[1]]=link[x,y]
        W-=.5*(np.kron(T,np.eye(4)-gamma[mu])+np.kron(T.conj().T,np.eye(4)+gamma[mu]))
    H=np.kron(np.eye(n),g5)@W;w,E=np.linalg.eigh(H)
    return E[:,w>0],float(min(abs(w))),int(sum(w<0)-sum(w>0))//2


def chart_frame(E,reference):
    u,s,vh=np.linalg.svd(reference.conj().T@E,full_matrices=False)
    assert min(s)>.1
    return E@vh.conj().T@u.conj().T,float(min(s))


def flux_transport():
    center=np.array([.31,.23]);h=.01;rows=[];L=16
    for p2,p3 in product([0.,np.pi],[np.pi/2,3*np.pi/2]):
        cache={}
        def get(x):
            key=tuple(np.round(x,12))
            if key not in cache:cache[key]=flux_frame(*x,p2,p3,L)
            return cache[key]
        E,g,index=get(center);P=E@E.conj().T
        curves=[]
        for step in [h,h/2]:
            d=[]
            for j in range(2):
                axis=np.eye(2)[j]*step;A=get(center+axis)[0];B=get(center-axis)[0]
                d.append((A@A.conj().T-B@B.conj().T)/(2*step))
            curves.append(float((-1j*np.trace(P@(d[0]@d[1]-d[1]@d[0]))).real))
        # Berry holonomy of the same positively oriented rectangle.
        corners=[center,center+[h,0],center+[h,h],center+[0,h],center]
        hol=np.eye(E.shape[1],dtype=complex);minimum=1.
        for a,b in zip(corners,corners[1:]):
            u,s,vh=np.linalg.svd(get(a)[0].conj().T@get(b)[0],full_matrices=False);hol=hol@u@vh;minimum=min(minimum,min(s))
        # Explicit determinant-line connection on a local Stiefel chart.
        point=center+np.array([h,h]);F,_=chart_frame(get(point)[0],E);current=[]
        for j in range(2):
            d=np.eye(2)[j]*h/2;plus,_=chart_frame(get(point+d)[0],E);minus,_=chart_frame(get(point-d)[0],E)
            current.append(float((-1j*np.trace(F.conj().T@(plus-minus)/h)).real))
        phase=float(np.angle(np.linalg.det(hol)))
        rows.append(dict(spectator_momenta=[p2,p3],Wilson_gap=min(v[1] for v in cache.values()),overlap_index=index,curvatures=curves,rectangle_phase=phase,rectangle_phase_over_area=phase/h**2,chart_connection=current,minimum_overlap_singular_value=float(minimum),sample_count=len(cache)))
    u,v=N.P.N.uniform_flux_links(L);plaquette=u*np.roll(v,-1,axis=0)*np.roll(u.conj(),-1,axis=1)*v.conj()
    maximum=float(max(abs(1-plaquette).ravel()));assert maximum<1/30
    assert min(r['Wilson_gap'] for r in rows)>.1
    total=sum(r['curvatures'][-1] for r in rows)
    return dict(status='PASS',lattice_shape=[L,L,2,2],magnetic_flux=float(np.angle(plaquette).sum()/(2*np.pi)),maximum_plaquette=maximum,toron_center=center.tolist(),rows=rows,total_curvature=total,total_rectangle_phase=sum(r['rectangle_phase'] for r in rows),scope='Exact spectator Fourier decomposition of one charged4D Weyl projector in an admissible nonzero magnetic-flux sector. Implements a local chart connection and polar transport on a two-toron slice. The connection is a field-space chart object, not a proved spatially local anomaly-cancelling measure current; all-SM-charge admissibility and global transition cocycles remain open. One magnetic plane has zero4D index and does not give a double-flux index witness.')


def lorentz_terms(X,epsilon=1e-10):
    """Upper Wick branch; the stored fixture has ten timelike triangles."""
    z=-1+1j*epsilon;g=np.diag([z,1,1,1]);cov=np.linalg.inv(np.c_[np.ones(5),X])[1:,:].T
    norms=cov@np.linalg.inv(g)@cov.T;roots=np.sqrt(np.diag(norms)+0j)
    area=[];angles=[]
    for i,j in combinations(range(5),2):
        tri=[k for k in range(5) if k not in [i,j]];d=X[tri[1:]]-X[tri[0]]
        area.append(.5*np.sqrt(np.linalg.det(d@g@d.T)+0j))
        angles.append(np.arccos(-norms[i,j]/(roots[i]*roots[j])))
    volume=np.sqrt(z)*abs(np.linalg.det(X[1:]-X[0]))/24
    return np.array(area),np.array(angles),volume


def coupled_lorentz():
    X=np.array(N.previous()['lorentzian']['coordinates']);f=np.array([.3,-.2,.7,.1,-.5]);weights=np.array([.1,.15,.2,.25,.3]);Newton=1.;Lambda=.001
    eta=np.diag([-1.,1,1,1]);causal=[]
    for tri in combinations(range(5),3):
        d=X[list(tri)[1:]]-X[tri[0]];w=np.linalg.eigvalsh(d@eta@d.T)
        causal.append(dict(triangle=list(tri),Gram_eigenvalues=w.tolist(),type='timelike' if min(w)<0<max(w) else 'spacelike'))
    assert all(r['type']=='timelike' for r in causal)
    def grav(Y,eps=1e-10):
        A,T,V=lorentz_terms(Y,eps);return -1j*(A@(np.pi-T)/Newton-Lambda*V)
    def matter(Y):return float(f@N.P.scalar_refinement(Y,weights,.2)[0]@f/2)
    a,t,_=lorentz_terms(X);force=np.zeros_like(X,dtype=complex);matterforce=np.zeros_like(X);errors=[];schlafli=[];h=1e-5
    for i,j in product(range(5),range(4)):
        d=np.zeros_like(X);d[i,j]=1
        force[i,j]=(grav(X+h*d)-grav(X-h*d))/(2*h)
        small=(grav(X+h*d/2)-grav(X-h*d/2))/h;errors.append(abs(force[i,j]-small))
        matterforce[i,j]=(matter(X+h*d)-matter(X-h*d))/(2*h)
        tp=lorentz_terms(X+h*d)[1];tm=lorentz_terms(X-h*d)[1];schlafli.append(abs(a@(tp-tm)/(2*h)))
    boost=np.eye(4);boost[:2,:2]=[[np.cosh(.4),np.sinh(.4)],[np.sinh(.4),np.cosh(.4)]]
    boosterror=abs(grav(X@boost.T)-grav(X));translation=max(abs(np.sum(force,axis=0)))
    assert max(schlafli)<1e-6 and max(errors)<1e-5 and boosterror<1e-6 and translation<1e-5
    return dict(status='PASS',coordinates=X.tolist(),causal_hinges=causal,supplied_Newton=Newton,supplied_Lambda=Lambda,gravitational_action=[float(grav(X).real),float(grav(X).imag)],matter_action=matter(X),gravity_force=enc(force),matter_force=matterforce.tolist(),combined_force=enc(force+matterforce),finite_difference_error=float(max(errors)),Schlafli_error=float(max(schlafli)),boost_error=float(boosterror),translation_force_error=float(translation),Wick_regulator_difference=float(abs(grav(X,1e-8)-grav(X))),scope='Single supplied Lorentzian4-simplex with ten timelike hinges, upper Wick branch and boundary Regge rotation angles; coupled to the prior stationary-eliminated scalar. Named complex action and all20 coordinate forces are computed, including a Schlaefli and boost check. No interior curvature, general causal-sector gluing, stationary coupled geometry, derived Newton coupling or cosmological constant follows; spacelike-hinge boost branches remain open.')


def nonlinear_joint():
    M=N.P.N.M;old=M.previous();A=M.P.load()['A'];pins=[M.P.read(old['pin_vacua'][key]) for key in ['pin_up','pin_down']];S=np.diag([1.,2.,3.]);angles=np.array(json.loads(N.OUT.read_text())['joint_slice']['pin_angles'])
    ren=json.loads(N.P.N.OUT.read_text())['renormalized_slice'];coef=np.array(ren['counterterm_coefficients'],dtype=np.longdouble)
    cache={}
    def link(phi):
        key=round(float(phi),12)
        if key not in cache:
            masses=[np.asarray(M.full_quark_masses(A,B,key),dtype=np.longdouble) for B in pins]
            r=np.longdouble('.2')*(3*np.longdouble(key)**2-np.longdouble('.8')**2);z=np.longdouble('.2')*(np.longdouble(key)**2-np.longdouble('.8')**2)
            assert z>0
            V=-12*sum(np.sum(m**4*(np.log(m*m)-np.longdouble('1.5'))) for m in masses)
            V+=480*sum(x*x*(np.log(x)-np.longdouble('1.5')) for x in [r,z]);V/=64*np.longdouble(str(np.pi))**2
            V+=480*np.longdouble('.1')*(np.longdouble(key)**2-np.longdouble('.8')**2)**2
            V+=np.array([1,np.longdouble(key)**2,np.longdouble(key)**4])@coef
            cache[key]=V
        return cache[key]
    base=link(.81);nubase=M.neutrino_phase_potential(pins[1],S,angles,precision=40)
    def potential(x):
        phi,s=x
        return float(link(phi)-base)+100*(s*s-1)**2+float(M.neutrino_phase_potential(pins[1],s*S,angles,precision=40)-nubase)+.1*(phi-.81)*(s-1)
    h=1e-4
    def gradient(x):return np.array([(potential(x+h*d)-potential(x-h*d))/(2*h) for d in np.eye(2)])
    sol=minimize(potential,[.81,1.0016],jac=gradient,method='L-BFGS-B',bounds=[(.805,.82),(.97,1.03)],options={'gtol':2e-5,'ftol':1e-14,'maxiter':60,'maxls':30})
    x=sol.x;g=gradient(x);H=np.column_stack([(gradient(x+h*d)-gradient(x-h*d))/(2*h) for d in np.eye(2)]);H=(H+H.T)/2
    small=np.array([(potential(x+h*d/2)-potential(x-h*d/2))/h for d in np.eye(2)])
    assert np.linalg.norm(g)<1e-3 and min(np.linalg.eigvalsh(H))>50
    # Explicit Wilsonian cutoff choices in the elementary-mediator reading.
    masses=np.sort(np.concatenate([M.full_quark_masses(A,B,x[0]) for B in pins]));cutoffs=[]
    for cutoff in [20.,100.,1000.]:
        points=np.r_[5.,masses[(masses>5)&(masses<cutoff)],cutoff];inv=25.
        for lo,hi in zip(points[:-1],points[1:]):
            nf=int(sum(masses<=np.sqrt(lo*hi)));inv+=(11-2*nf/3)*np.log(hi/lo)/(8*np.pi**2)
        cutoffs.append(dict(cutoff=cutoff,inverse_bare_g2=float(inv),positive_bare_kinetic_term=bool(inv>0)))
    return dict(status='PASS',coordinates=x.tolist(),force=g.tolist(),half_step_force=small.tolist(),Hessian=H.tolist(),Hessian_eigenvalues=np.linalg.eigvalsh(H).tolist(),potential_difference=potential(x),counterterms=ren['counterterm_coefficients'],supplied_portal=.1,retained_heavy_species=len(masses),high_energy_b0=11-2*len(masses)/3,finite_cutoffs=cutoffs,link_spectral_evaluations=len(cache),scope='Nonlinear two-variable link/Majorana vacuum slice using all486 singular values per quark sector and the actual6-state neutrino determinant, with retained heavy terms. Counterterms, pin orientations, neutrino angles, tree coefficients and low-scale coupling g(5)=.2 are supplied. A positive finite-cutoff bare kinetic coefficient is an EFT consistency condition, not UV completion or a continuum limit. Joint orientation variations and the full field Hessian remain open.')


def transpositions(perm):
    visited=set();swaps=[]
    for i in range(len(perm)):
        if i in visited:continue
        cycle=[];j=i
        while j not in visited:
            visited.add(j);cycle.append(j);j=int(perm[j])
        # Applying these swaps left-to-right realizes columns e_i -> e_perm[i].
        swaps.extend((cycle[0],j) for j in cycle[1:])
    return swaps


def reset_synthesis():
    p=N.reset_target();n=12;perm=np.array(p['permutation']);U=np.eye(n*n,dtype=complex);swaps=transpositions(perm)
    # Two pulses per transposition: exp(-i*pi*Xab/2), followed by
    # exp(+i*pi*(Pa+Pb)/2). Their product is exactly an unsigned swap.
    for a,b in swaps:
        U[[a,b]]=U[[b,a]]
    target=np.eye(n*n)[:,perm];error=np.linalg.norm(U-target);assert error<1e-12
    logical=np.diag([1.,1.]+[0.]*10);P=np.kron(logical,np.eye(n));obstruction=np.linalg.norm(P@target-target@P)
    assert obstruction>1
    # Reuse two SWAP-routed syndrome interactions; one initial relay Pauli
    # survives the first interaction and may propagate to data/second syndrome.
    ops1=[(0,2),(2,0),(0,2),(2,3),(0,2),(2,0),(0,2)]
    ops2=[(1,2),(2,1),(1,2),(2,4),(1,2),(2,1),(1,2)]
    def prop(x,z,ops):
        for c,t in ops:x[t]^=x[c];z[c]^=z[t]
        return x,z
    rows=[]
    for location,(c,t) in enumerate(ops1):
        for pa,pb in product(range(4),repeat=2):
            if pa==pb==0:continue
            x=np.zeros(5,dtype=int);z=x.copy()
            x[c]=pa&1;z[c]=pa>>1;x[t]=pb&1;z[t]=pb>>1
            x,z=prop(x,z,ops1[location+1:]);first=[i for i in range(5) if x[i] or z[i]]
            xx,zz=prop(x.copy(),z.copy(),ops2);xr=x.copy();zr=z.copy();xr[2]=zr[2]=0;xr,zr=prop(xr,zr,ops2)
            mask=[0,1,3,4];difference=int(np.sum(xx[mask]!=xr[mask])+np.sum(zz[mask]!=zr[mask]))
            rows.append(dict(location=location,paulis=[pa,pb],first_support=first,reused_support=[i for i in range(5) if xx[i] or zz[i]],fresh_reset_support=[i for i in range(5) if xr[i] or zr[i]],endpoint_difference=difference))
    # Prove route factorization on all32 computational basis states.
    def route(bits,ops):
        bits=bits.copy()
        for c,t in ops:bits[t]^=bits[c]
        return bits
    for bits in product([0,1],repeat=5):
        expected=list(bits);expected[3]^=expected[0]
        assert route(list(bits),ops1)==expected
        expected=list(bits);expected[4]^=expected[1]
        assert route(list(bits),ops2)==expected
    assert all(r['endpoint_difference']==0 for r in rows)
    return dict(status='PASS',transpositions=[list(s) for s in swaps],pulse_count=2*len(swaps),compiled_unitary_error=float(error),logical_sector_commutator_norm=float(obstruction),relay_reuse=rows,relay_residues=sum(2 in r['first_support'] for r in rows),endpoint_changes_due_to_relay_reset=sum(r['endpoint_difference'] for r in rows),scope='Complete algebraic two-level pulse decomposition of the stored144-state reset target. Joint leakage-sector Xab and diagonal pulses are additional supplied controls; all controls preserving the system logical projector cannot synthesize this reset. Reinitializing the environment is an open-system resource. All105 first-route Pauli faults are propagated through a second ideal route. Each ideal route is exactly CNOT on endpoints tensor identity on relay, so perfect relay reset changes no endpoint Pauli. Persistent relay residues alone do not prove harmful ideal reuse or a required reset cost; correlated or state-dependent noisy gates still require a full audit. No noisy decoder or threshold is established.')



def run():
    r=dict(status='PASS',passes=list(range(11461,11466)),reservation='3e9cc0867')
    for key,fun in [('charges',charge_map),('flux',flux_transport),('gravity',coupled_lorentz),('joint',nonlinear_joint),('reset',reset_synthesis)]:
        r[key]=fun();print(key,r[key]['status'],flush=True)
    r['photon_candidate']=photon_preserving_candidate(r['charges']);r['rank_sensitivity']=rank_sensitivity(np.array(r['photon_candidate']['coordinates']));print('photon_candidate',r['photon_candidate']['status'],flush=True)
    names=['data/w33_pass11448_11455_subgroup_toron_auxiliary_control.json','data/w33_pass11443_11447_preserved_gauge_lorentz_recovery.json','data/w33_pass11438_11442_coercive_curved_recovery.json']
    r['source_sha256']={n:hashlib.sha256(json.dumps(json.loads((ROOT/n).read_text()),sort_keys=True,separators=(',',':')).encode()).hexdigest() for n in names}
    OUT.write_text(json.dumps(r,indent=2)+'\n');return r

if __name__=='__main__':
    result=run();print({k:v['status'] for k,v in result.items() if isinstance(v,dict) and 'status' in v})
