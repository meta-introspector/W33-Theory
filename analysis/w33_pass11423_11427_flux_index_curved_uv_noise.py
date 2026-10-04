#!/usr/bin/env python3
"""Five constructed follow-ups to11413-11417, with explicit input boundaries.
Prior4084 overlap,11361/11374 phase selectors,11409 Yukawa arrows,
11415 Regge,11417 counted gates, PART_CCCCIV Steane. Generic methods prior.
"""
import hashlib,json
from pathlib import Path
from itertools import combinations,product
import numpy as np
import sympy as sp
from scipy.linalg import expm,block_diag
from scipy.optimize import brentq
from w33_pass11390_11394_context_matter_clock import load
from w33_pass11408_11412_clock_spurion_seam_ticks import prior,read
from w33_pass11413_11417_pin_regge_vacuum_noise import enc,hessian,simplex_geometry
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/w33_pass11423_11427_flux_index_curved_uv_noise.json'
def previous():return json.loads((ROOT/'data/w33_pass11413_11417_pin_regge_vacuum_noise.json').read_text())
def hermitian_basis():
    B=[np.diag(np.eye(3)[i]).astype(complex)for i in range(3)]
    for i,j in combinations(range(3),2):
        E=np.zeros((3,3),complex);E[i,j]=E[j,i]=1/np.sqrt(2);B.append(E)
        E=np.zeros((3,3),complex);E[i,j]=1j/np.sqrt(2);E[j,i]=-1j/np.sqrt(2);B.append(E)
    return B
def pin_potential(U,B,chi,a=3.,r=.25,vchi=.5,g=.1):
    poly=U@U@U-6*U@U+11*U-6*np.eye(3)
    up=np.linalg.norm(U-np.diag(np.diag(U)))**2+np.linalg.norm(poly)**2+(np.trace(U).real-6)**2+(np.trace(U@U).real-14)**2
    t=B[0,1]*B[1,2]*B[2,0]
    down=np.linalg.norm(np.diag(B)-a)**2+sum((abs(B[i,j])**2-r*r)**2 for i,j in [(0,1),(1,2),(2,0)])+(chi*chi-vchi*vchi)**2-g*chi*t.imag
    return float(up+down)
def pin_vacua():
    a,r,v,g=3.,.25,.5,.1
    # At a nonzero stationary point all three radii exceed r and are equal.
    # Eliminate chi, then certify the unique positive-radius root by Sturm.
    w=sp.symbols('w');rr,vv,gg=map(sp.Rational,['1/4','1/2','1/10'])
    poly=(256-gg**4)*w**3-(768*rr**2+16*vv**2*gg**2)*w**2+(768*rr**4+16*vv**2*gg**2*rr**2)*w-256*rr**6
    count=int(sp.Poly(poly,w).count_roots(rr**2,sp.oo));assert count==1
    u=np.sqrt(brentq(lambda x:float(poly.subs(w,x)),r*r+1e-10,1.))
    chi=4*(u*u-r*r)/(g*u);assert chi>v
    X=np.roll(np.eye(3),1,axis=0);U=np.diag([1.,2.,3.]);z=u*np.exp(-1j*np.pi/6);B=a*np.eye(3)+z*X+z.conjugate()*X.T
    basis=hermitian_basis();zero=np.zeros(19)
    def fn(x):return pin_potential(U+sum(q*e for q,e in zip(x[:9],basis)),B+sum(q*e for q,e in zip(x[9:18],basis)),chi+x[18])
    H=hessian(fn,zero,h=2e-5);ev=np.linalg.eigvalsh(H)
    assert np.count_nonzero(abs(ev)<1e-6)==2 and min(ev)>-1e-6
    rng=np.random.default_rng(11423);sym=[];omega=np.exp(2j*np.pi/3);Z=np.diag([1,omega,omega**2])
    for _ in range(10):
        A=sum(x*e for x,e in zip(rng.normal(size=9),basis));D=sum(x*e for x,e in zip(rng.normal(size=9),basis));c=rng.normal()
        e0=pin_potential(A,D,c)
        for G in [X,Z]:sym.append(abs(pin_potential(G@A@G.conj().T,G@D@G.conj().T,c)-e0))
        sym.append(abs(pin_potential(A.conjugate(),D.conjugate(),-c)-e0))
    assert max(sym)<1e-9
    comm=U@B-B@U;cp=float(np.imag(np.trace(comm@comm@comm)))
    assert abs(cp)>1e-3 and np.linalg.eigvalsh(B)[0]>0
    # Full pair invariant excludes a generalized CP mapping of this vacuum.
    commc=U@B.conjugate()-B.conjugate()@U
    assert abs(np.imag(np.trace(commc@commc@commc))+cp)<1e-12
    return dict(status='PASS',parameters=dict(a=a,r=r,vchi=v,g=g),pin_up=enc(U),pin_down=enc(B),chi=chi,radius=u,
        potential='Vup=||offdiag U||²+||U³-6U²+11U-6I||²+(tr U-6)²+(tr U²-14)²; Vdown=||diag B-a||²+sum_cycle(|Bij|²-r²)²+(chi²-vchi²)²-g chi Im(B01B12B20)',
        CP='U,B -> conjugates; chi -> -chi',Hessian=H.tolist(),Hessian_eigenvalues=ev.tolist(),minimum_value=fn(zero),
        radius_squared_polynomial=str(poly),positive_radius_root_count=count,triangle_flux=float(np.imag(B[0,1]*B[1,2]*B[2,0])),cubic_pair_CP=cp,symmetry_residual=float(max(sym)),
        scope='A native-family clock anisotropy and real-coupling potential select noncommuting full-rank pin minima with spontaneous opposite CP flux, without prescribing matrix entries or a phase. Spectral/radial scales, clock axis and couplings are supplied; two diagonal-rephasing flat directions remain. Up potential is degree6, a cutoff scalar EFT rather than a fundamental renormalizable UV completion. Analytic link path protection remains conditional; loop stability and measured parameters are open.')

def gamma_matrices():
    sx=np.array([[0,1],[1,0]]);sy=np.array([[0,-1j],[1j,0]]);sz=np.diag([1,-1]);I=np.eye(2)
    gamma=[np.kron(sx,s)for s in [sx,sy,sz]]+[np.kron(sy,I)]
    return gamma,gamma[0]@gamma[1]@gamma[2]@gamma[3]
def torus_links(L,flux=(1,1)):
    sites=list(product(range(L),repeat=4));lookup={x:i for i,x in enumerate(sites)};U=np.ones((len(sites),4),complex)
    for i,x in enumerate(sites):
        for first,m in [(0,flux[0]),(2,flux[1])]:
            U[i,first]=np.exp(-2j*np.pi*m*x[first+1]/L**2)
            if x[first+1]==L-1:U[i,first+1]=np.exp(2j*np.pi*m*x[first]/L)
    shifts=[]
    for mu in range(4):
        T=np.zeros((len(sites),len(sites)),complex)
        for i,x in enumerate(sites):
            y=list(x);y[mu]=(y[mu]+1)%L;T[i,lookup[tuple(y)]]=U[i,mu]
        shifts.append(T)
    return U,shifts
def overlap_operator(L=4,flux=(1,1),mass=1.):
    gamma,g5=gamma_matrices();U,T=torus_links(L,flux);n=L**4;G=np.kron(np.eye(n),g5)
    W=(4-mass)*np.eye(4*n,dtype=complex)
    for mu in range(4):W-=(np.kron(T[mu],np.eye(4)-gamma[mu])+np.kron(T[mu].conj().T,np.eye(4)+gamma[mu]))/2
    Hw=G@W;assert np.linalg.norm(Hw-Hw.conj().T)<1e-12
    w,V=np.linalg.eigh(Hw);sign=(V*np.sign(w))@V.conj().T;D=np.eye(4*n)+G@sign
    return U,T,G,Hw,w,D,sign
def spacetime_index():
    rows=[]
    for flux in [(1,1),(-1,1),(0,0)]:
        U,T,G,Hw,w,D,S=overlap_operator(flux=flux);index=int(round(-np.sign(w).sum()/2))
        GW=np.linalg.norm(G@D+D@G-D@G@D);hat=G@(np.eye(len(D))-D)
        intertwine=np.linalg.norm(D@hat+G@D);assert GW<1e-10 and intertwine<1e-10
        if flux==(0,0):assert index==0
        else:assert abs(index)==1
        row=dict(flux=list(flux),index=index,Wilson_gap=float(min(abs(w))),GW_residual=float(GW),Weyl_intertwining_residual=float(intertwine),
            hat_projector_plus_rank=int(round((len(D)+np.trace(hat).real)/2)),hat_projector_minus_rank=int(round((len(D)-np.trace(hat).real)/2)),
            ordinary_plus_rank=len(D)//2,links=enc(U),Wilson_sign_counts=[int((w<0).sum()),int((w>0).sum())])
        if flux==(1,1):
            vals,B=np.linalg.eigh(D.conj().T@D);zero=B[:,vals<1e-10];chir=np.linalg.eigvalsh(zero.conj().T@G@zero)
            assert len(chir)==1 and abs(chir[0]-index)<1e-10
            row.update(zero_mode=enc(zero),zero_mode_chirality=chir.tolist(),zero_residual=float(np.linalg.norm(D@zero)))
        rows.append(row)
    assert rows[0]['index']==-rows[1]['index']
    # A named extension of the native rank9 colour/family frame.
    slots=prior()['fermion_bimodule']['native_slot_frames'];C=np.hstack([read(x)for x in slots['colour_qutrit_seeds']]);assert np.linalg.norm(C.conj().T@C-np.eye(9))<1e-12
    return dict(status='PASS',L=4,spacetime_spinor_dimension=1024,Wilson_mass=1.,scans=rows,native_family_colour_rank=9,tensor_index=9*rows[0]['index'],
        extension='H_spacetime tensor native rank9 colour-family carrier C; Weyl action maps hat-P-minus to ordinary P-plus via overlap D',
        scope='Standard4D Euclidean overlap on a supplied compact torus and U1 flux has an actual nonzero background index and rectangular lattice Weyl projectors; attach the actual native colour/family frame. This is an explicit spacetime extension, not emergent spacetime, Lorentzian particle content, SM mirror removal or an anomaly-safe chiral gauge measure. Flux, lattice, spin matrices and projection choice remain inputs; background Dirac index alone does not select the SM.')

def mesh_action(L,simplices):
    counts={};angles={}
    for s in simplices:
        for tet in combinations(s,4):key=tuple(sorted(tet));counts[key]=counts.get(key,0)+1
        for tri,theta in simplex_geometry(list(s),L).items():angles[tri]=angles.get(tri,0)+theta
    boundary=set()
    for tet,count in counts.items():
        if count==1:boundary.update(tuple(sorted(t))for t in combinations(tet,3))
    action=0.;deficits={}
    for tri,total in angles.items():
        sides=[L[i,j]for i,j in combinations(tri,2)];a,b,c=sides;s=(a+b+c)/2;area=np.sqrt(s*(s-a)*(s-b)*(s-c))
        d=(np.pi if tri in boundary else 2*np.pi)-total;action+=area*d
        if tri not in boundary:deficits[tri]=d
    return float(action),deficits
def local_simplex_coordinates(L,s):
    o=s[0];vs=s[1:];G=np.array([[(L[o,i]**2+L[o,j]**2-L[i,j]**2)/2 for j in vs]for i in vs])
    X=np.linalg.cholesky(G);return np.vstack([np.zeros(4),X])
def curved_gauge_islands():
    old=previous()['regge_patch'];period=np.array(old['spatial_supercell_periods']);affine=block_diag(period,2.)
    V=np.array([[1,0,0,0],[0,1,0,0],[-1,-1,0,0],[0,0,1,0],[0,0,0,1],[0,0,-1,-1]])@affine.T
    boundary=np.linalg.norm(V[:,None]-V[None,:],axis=2);coarse=[(0,1,2,3,4),(0,1,2,4,5),(0,1,2,5,3)]
    boundary[0,3]*=1.002;boundary[3,0]=boundary[0,3]
    s0,def0=mesh_action(boundary,coarse);assert abs(def0[(0,1,2)])>1e-5
    coords=[local_simplex_coordinates(boundary,list(s))for s in coarse]
    radii=[np.linalg.norm(x-x.mean(axis=0),axis=1)for x in coords];flat=np.concatenate(radii)
    fine=[]
    for i,s in enumerate(coarse):fine.extend(tuple([6+i]+[v for v in s if v!=omitted])for omitted in s)
    def lengths(rad):
        L=np.zeros((9,9));L[:6,:6]=boundary
        for i,s in enumerate(coarse):
            for j,v in enumerate(s):L[6+i,v]=L[v,6+i]=rad[5*i+j]
        return L
    base,df=mesh_action(lengths(flat),fine);assert abs(base-s0)<1e-10
    assert abs(df[(0,1,2)]-def0[(0,1,2)])<1e-11
    assert max(abs(d)for t,d in df.items()if any(v>=6 for v in t))<1e-11
    tangents=[];tests=[];rng=np.random.default_rng(11425)
    for i,x in enumerate(coords):
        center=x.mean(axis=0);T=(center-x)/radii[i][:,None];tangents.append(T)
        for _ in range(4):
            r=flat.copy();shift=rng.normal(size=4)*.005;r[i*5:i*5+5]=np.linalg.norm(x-center-shift,axis=1)
            a,d=mesh_action(lengths(r),fine);tests.append(abs(a-base));assert abs(a-base)<1e-10
    # Refinement action separable in each simplex's radial variables.
    Hess=[];res=[]
    for i in range(3):
        def fn(q):
            r=flat.copy();r[5*i:5*i+5]=q;return mesh_action(lengths(r),fine)[0]
        # Schlaefli identity: angle derivatives cancel from the action gradient.
        # Radial edges only change triangles containing their interior vertex.
        def gradient(q):
            r=flat.copy();r[5*i:5*i+5]=q;L=lengths(r);_,d=mesh_action(L,fine);out=np.zeros(5)
            for j,v in enumerate(coarse[i]):
                for other in coarse[i]:
                    if other==v:continue
                    tri=tuple(sorted([6+i,v,other]));ell=q[j];opp=L[v,other];adj=L[6+i,other]
                    semi=(ell+opp+adj)/2;area=np.sqrt(semi*(semi-ell)*(semi-opp)*(semi-adj))
                    out[j]+=ell*(opp**2+adj**2-ell**2)/(8*area)*d[tri]
            return out
        step=1e-6;E=np.eye(5)*step
        H=np.column_stack([(gradient(radii[i]+e)-gradient(radii[i]-e))/(2*step)for e in E])
        assert np.linalg.norm(H-H.T)<.001
        H=(H+H.T)/2;Hess.append(H);res.append(float(np.linalg.norm(H@tangents[i])))
    assert max(res)<.01
    r=flat.copy();_,_,vh=np.linalg.svd(tangents[0].T);r[:5]+=.005*vh[-1]
    _,curved=mesh_action(lengths(r),fine);assert max(abs(d)for t,d in curved.items()if 6 in t)>1e-4
    return dict(status='PASS',boundary_vertices=V.tolist(),boundary_lengths=boundary.tolist(),coarse_simplices=coarse,refined_simplices=fine,
        coarse_internal_hinge=[0,1,2],coarse_hinge_deficit=def0[(0,1,2)],action=base,radial_lengths=flat.tolist(),
        local_coordinates=[x.tolist()for x in coords],vertex_translation_tangents=[x.tolist()for x in tangents],radial_Hessian_blocks=[h.tolist()for h in Hess],flat_star_residuals=res,
        vertex_translation_directions=12,nonlinear_action_residuals=tests,curvature_in_first_star=float(max(abs(d)for t,d in curved.items()if 6 in t)),
        scope='Three coarse4-simplices share a genuinely curved internal triangle and are refined into15 simplices. Three newly inserted flat vertex stars retain twelve nonlinear translation directions despite global curvature. Curvature inside a moved vertex star is a different obstruction. This standard subdivision/gauge mechanism is a native-period witness, not gauge symmetry for generic curved stars, a perfect action, Lorentzian gravity or the full Einstein constraint algebra.')

def heavy_matrix(A,B,M=10.,eta=.08,v=1.):
    n=len(A);k=np.kron(M*(np.eye(n)-eta*A),np.eye(3));left=np.zeros((n*3,3));left[:3]=v*np.eye(3)/np.sqrt(2);right=np.zeros((3,n*3),complex);right[:,:3]=B
    return np.block([[k,left,np.zeros_like(k)],[np.zeros((3,n*3)),M*np.eye(3),right],[np.zeros_like(k),np.zeros((n*3,3)),k]])
def heavy_extension(c,p):
    old=previous();A=c['A'];idx=old['path_flavour']['endpoints'];M,eta,v=10.,.08,1.;spectra=[];checks=[]
    for name,key in [('up','pin_up'),('down','pin_down')]:
        B=read(p[key]);H=heavy_matrix(A,B,M,eta,v);n=len(H);AL=np.zeros((3,n));AR=np.zeros((n,3))
        for j,i in enumerate(idx):AL[j,i*3+j]=1;AR[n-240+i*3+j,j]=1
        Full=np.block([[np.zeros((3,3)),AL],[AR,H]])
        invR=np.linalg.solve(H,AR);Y=-AL@invR;R=np.linalg.inv(np.eye(80)-eta*A);r=R[idx,0]
        expected=-v/np.sqrt(2)/M**3*np.diag(r)@B@np.diag(r)
        assert np.linalg.norm(Y-expected)<1e-12
        masses=np.linalg.svd(Full,compute_uv=False);heavy=np.linalg.svd(H,compute_uv=False)
        spectra.append(dict(name=name,full_dimension=len(Full),heavy_dimension=n,full_Dirac_masses=masses.tolist(),heavy_Dirac_masses=heavy.tolist(),matching=enc(Y),full_light_masses=masses[-3:].tolist()))
        checks.append(float(np.linalg.norm(Y-expected)))
    # Link phases are physical flat tree directions here; node symmetries are
    # formal/global, not extra gauged anomaly-free U1s.
    link_count=3*160;phi=eta*M;lambda_link=.1;link_mass2=4*lambda_link*phi*phi
    scalar=np.r_[old['majorana_vacuum']['scalar_mass_squared'],old['majorana_vacuum']['Higgs_mass_squared'],p['Hessian_eigenvalues']]
    scalar=np.maximum(scalar,0);groups=[]
    for s in spectra:groups.append(dict(name=s['name']+'_full',weight=-12,m2=(np.array(s['full_Dirac_masses'])**2).tolist(),constant=1.5))
    charged=np.array([.0003,.003,.03])/np.sqrt(2)
    # Use the dynamically selected pin in the lepton EFT too, rather than
    # retaining the previous packet's explicit complex Yukawa source.
    Ynu=.01*read(p['pin_down']);S0=read(old['majorana_vacuum']['Majorana_source'])
    seesaw=np.block([[np.zeros((3,3)),Ynu/np.sqrt(2)],[Ynu.T/np.sqrt(2),S0]])
    nu=np.linalg.svd(seesaw,compute_uv=False)
    assert np.linalg.norm(heavy_matrix(A,read(p['pin_down']).conjugate()).conjugate()-heavy_matrix(A,read(p['pin_down'])))<1e-12
    groups.extend([dict(name='charged_lepton',weight=-4,m2=(charged**2).tolist(),constant=1.5),dict(name='Majorana6',weight=-2,m2=(nu**2).tolist(),constant=1.5),dict(name='Higgs_S_pin',weight=1,m2=scalar.tolist(),constant=1.5),dict(name='link_radials480',weight=link_count,m2=[link_mass2],constant=1.5),dict(name='W',weight=6,m2=[.65**2/4],constant=5/6),dict(name='Z',weight=3,m2=[(.65**2+.35**2)/4],constant=5/6)])
    # Numerical pin zero modes are thresholded using certified Hessian tolerance.
    groups[4]['m2']=[x if x>1e-6 else 0. for x in groups[4]['m2']]
    def value(mu):
        return sum(g['weight']*sum(x*x*(np.log(x/mu**2)-g['constant'])for x in g['m2']if x>0)for g in groups)/(64*np.pi**2)
    str4=sum(g['weight']*sum(x*x for x in g['m2'])for g in groups);val=value(1.)
    assert abs(value(2)-val+str4*np.log(2)/(32*np.pi**2))<1e-7
    heavy_str4=sum(-12*sum(x**4 for x in s['heavy_Dirac_masses'])for s in spectra)
    return dict(status='PASS',parameters=dict(M=M,eta=eta,v=v,lambda_link=lambda_link,link_vev=phi),spectra=spectra,
        matching_residuals=checks,heavy_representations='Q:80*3 family Dirac(3,2,1/6); U,D:80*3 Dirac(3,1,2/3 or-1/3); F_u,F_d:3 each matching right singlet reps. All vectorlike.',
        heavy_Dirac_count=2898,full_Dirac_count=2919,Majorana_count=6,fermion_real_weight=11688,
        real_scalar_count=995,complex_link_count=480,link_radial_mass_squared=link_mass2,tree_zero_scalar_count=485,
        heavy_bridge='Q -- y_H Higgs --> F -- neutral family pin B --> right network; all dimension4 fermion vertices',
        neutrino_EFT='Ynu=.01 B_down with real coefficient; dimension5 lepton-pinning vertex, cutoff EFT',neutrino_seesaw=enc(seesaw),neutrino_masses=nu.tolist(),
        CP_branch_threshold_identity='All link/Higgs/attachment coefficients and S0 are real; conjugating the selected pin conjugates the full fermion mass matrices, so their singular spectra and one-loop vacuum thresholds coincide on the two CP branches.',
        groups=groups,supertrace_M4=str4,one_loop_value=val,heavy_only_supertrace_M4=heavy_str4,
        scale_identity='dV1/dln mu=-STr M4/(32pi²) at fixed masses',
        scope='Explicit finite heavy extension of the pin-matching low-energy model: vectorlike quark messengers, pin mediators,480 complex neutral link fields,19 neutral pin coordinates and prior Higgs/Majorana fields. All declared quadratic eigenvalues and heavy thresholds are counted; node phases are not gauged. Up pin potential remains a cutoff EFT, scalar flat directions require loop analysis, and vacuum counterterm is free. This is not a fundamental UV completion or a protected cosmological constant.')

def invariant_frame(c,seeds):
    D=np.array(c['D'],float);Hp=D[:40].T@D[:40];Hl=D[40:].T@D[40:];basis=[]
    def add(v):
        v=v.astype(complex).copy()
        for _ in range(2):
            for b in basis:v-=b*np.vdot(b,v)
        if np.linalg.norm(v)>1e-10:basis.append(v/np.linalg.norm(v))
    for v in seeds:add(v)
    i=0
    while i<len(basis):
        for H in [Hp,Hl]:add(H@basis[i])
        i+=1
    B=np.column_stack(basis);assert max(np.linalg.norm(H@B-B@(B.conj().T@H@B))for H in [Hp,Hl])<1e-10
    return B,Hp,Hl
def steane_ML(weights):
    # Exact joint syndrome/logical-class probabilities for all Pauli patterns,
    # conditional on each known erasure support. Erased sites = uniform Pauli.
    i,x,z=np.indices((1,128,128));x=x.ravel();z=z.ravel();syx=np.zeros(len(x),int);syz=syx.copy();px=syx.copy();pz=syx.copy();pairs=[]
    for j in range(7):
        xb=(x>>j)&1;zb=(z>>j)&1;syx^=(j+1)*xb;syz^=(j+1)*zb;px^=xb;pz^=zb;pairs.append(xb+2*zb)
    labels=4*(syx+8*syz)+px+2*pz
    w=np.array([weights[0],weights[1],weights[3],weights[2]]) # I,X,Z,Y
    success=total=0.
    for mask in range(128):
        prob=np.ones(len(x))
        for j in range(7):prob*=weights[4]/4 if (mask>>j)&1 else w[pairs[j]]
        groups=np.bincount(labels,weights=prob,minlength=256).reshape(64,4)
        success+=groups.max(axis=1).sum();total+=groups.sum()
    assert abs(total-1)<1e-10
    return float(success)
def hamming_syndrome(mask):
    out=0
    for j in range(7):
        if (mask>>j)&1:out^=j+1
    return out
def minimum_weight_decoder():
    dec={0:(0,0)}
    for weight in [1,2]:
        for support in combinations(range(7),weight):
            for types in product([(1,0),(0,1),(1,1)],repeat=weight):
                x=sum(a<<j for j,(a,b)in zip(support,types));z=sum(b<<j for j,(a,b)in zip(support,types))
                dec.setdefault(hamming_syndrome(x)+8*hamming_syndrome(z),(x,z))
    assert len(dec)==64;return dec
def syndrome_circuit(initial=(0,0),fault=None):
    # Seven data modules, one reused bare ancilla. Z checks then X checks.
    x,z=initial;syndrome=0;gate=0
    for kind in [0,1]:
        for row in range(3):
            x&=127;z&=127
            for data in [j for j in range(7)if ((j+1)>>row)&1]:
                control,target=(data,7)if kind==0 else (7,data)
                if (x>>control)&1:x^=1<<target
                if (z>>target)&1:z^=1<<control
                if fault is not None and gate==fault[0]:
                    bits=fault[1]
                    if bits&1:x^=1<<control
                    if bits&2:z^=1<<control
                    if bits&4:x^=1<<target
                    if bits&8:z^=1<<target
                gate+=1
            bit=((x if kind==0 else z)>>7)&1;syndrome|=bit<<(row+3*kind)
    return x&127,z&127,syndrome
def extraction_fault_audit():
    dec=minimum_weight_decoder();rows=[];bad=0;witness=None
    for gate in range(24):
        fails=0
        for bits in range(1,16):
            x,z,s=syndrome_circuit(fault=(gate,bits));rx,rz=dec[s];x^=rx;z^=rz
            failed=bool(hamming_syndrome(x)or hamming_syndrome(z)or x.bit_count()%2 or z.bit_count()%2)
            fails+=failed
            if x.bit_count()+z.bit_count()>=2 and failed and witness is None:
                witness=dict(gate=gate,two_qubit_Pauli_bits=bits,residual_X_mask=x,residual_Z_mask=z,reported_syndrome=s)
        rows.append(fails);bad+=fails
    assert bad>0
    for j in range(7):
        for a,b in [(1,0),(0,1),(1,1)]:
            x,z,s=syndrome_circuit((a<<j,b<<j));assert (x,z)==dec[s]
    return dict(CNOT_count=24,checks=6,data_qubits=7,reused_bare_ancillas=1,malignant_single_faults=bad,
        tested_faults=360,failure_counts_per_CNOT=rows,linear_failure_coefficient=bad/15,correlated_fault_witness=witness,
        scope='Explicit ideal Clifford extraction corrects all initial single data Paulis; uniform nonidentity two-qubit Pauli faults after each CNOT yield this first-order failure coefficient. Bare ancilla hooks make it non-fault-tolerant. Native CNOT compilation, leakage-aware measurement, repeated checks and full multi-fault channel remain open.')
def native_noise(c):
    p=prior()['local_gates'];old=previous()['counted_noise'];Q=read(p['logical_basis']);a,b=p['gate_pair'];target=read(p['ideal_qubit_gate']);ea,eb=np.eye(160)[:,[a,b]].T
    B,Hp,Hl=invariant_frame(c,[Q[:,0],Q[:,1],ea,eb]);Kp=B.conj().T@Hp@B;Kl=B.conj().T@Hl@B;n=B.shape[1];q0=B.conj().T@Q
    f=[B.conj().T@e for e in [ea,eb]];Ps=[np.outer(e,e.conj())for e in f];Zs=[np.eye(n)-2*P for P in Ps]
    clock=expm(-.02j*Kp)@expm(-.04j*Kl)@expm(-.02j*Kp);I=np.eye(n*n);rows=[]
    paulis=[np.eye(2),np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1])]
    bell=np.array([1,0,0,1])/np.sqrt(2)
    initial=np.column_stack([(q0[:,i,None]@q0[:,j,None].conj().T).reshape(-1,order='F')for i,j in product(range(2),repeat=2)])
    for probability in [0.,1e-8,1e-7,1e-6]:
        dephase=I.copy()
        for Z in Zs:dephase=((1-probability)*I+probability*np.kron(Z.conjugate(),Z))@dephase
        state=initial.copy()
        for angle,axis,N in zip(old['principal_angles'],[0,1,0,1,0],old['ticks_per_pulse']):
            kick=np.eye(n)+(np.exp(-.02j*np.copysign(.0125,angle))-1)*Ps[axis];S=kick@clock@kick
            channel=dephase@np.kron(S.conjugate(),S);state=np.linalg.matrix_power(channel,N)@state
        C=np.zeros((4,4),complex);trace=[]
        for k,(i,j)in enumerate(product(range(2),repeat=2)):
            density=state[:,k].reshape(n,n,order='F');logical=target.conj().T@q0.conj().T@density@q0@target
            C[2*i:2*i+2,2*j:2*j+2]=logical/2
            if i==j:trace.append(float(np.trace(density).real))
        assert max(abs(np.array(trace)-1))<1e-7
        weights=[]
        for P in paulis:
            v=np.kron(np.eye(2),P)@bell;weights.append(float(np.vdot(v,C@v).real))
        erasure=1-sum(weights);assert min(weights)>-1e-9 and erasure>-1e-8
        weights=[max(x,0.)for x in weights]+[max(erasure,0.)];weights=np.array(weights)/sum(weights)
        fidelity=steane_ML(weights);readout=[]
        for eps in [0.,1e-6,1e-5,1e-4]:readout.append(dict(bit_error_probability=eps,entanglement_fidelity=float(fidelity*(1-eps)**6),failure_probability=float(1-fidelity*(1-eps)**6)))
        break_even=1-(weights[0]/fidelity)**(1/6)if fidelity>weights[0]else 0.
        rows.append(dict(edge_flip_probability_per_tick=probability,Pauli_and_erasure_weights=weights.tolist(),uncoded_entanglement_infidelity=1-weights[0],ideal_encoded_entanglement_fidelity=fidelity,readout_scans=readout,one_shot_readout_break_even=float(break_even),
            Choi_eigenvalues=np.linalg.eigvalsh((C+C.conj().T)/2).tolist(),trace_residual=float(max(abs(np.array(trace)-1)))))
    return dict(status='PASS',native_invariant_dimension=n,invariant_frame=enc(B),ticks=old['total_ticks'],native_noise='After each actual tick, independent probability q of phase flip I-2|e><e| on both addressed native edges',
        scans=rows,syndrome_bits=6,extraction_circuit=extraction_fault_audit(),readout_factor='One-shot independent syndrome-bit errors leave nonzero residual syndrome unless all6 bits correct: F_e=(1-eps)^6 F_e(ideal readout)',
        scope='Exact counted native Hadamard density channel includes coherent control error and native edge-phase noise/leakage. Logical Pauli twirling plus ideal leakage flagging defines a stated CPTP Pauli/erasure model; exact maximum-likelihood Steane recovery counts all patterns and erasure supports. Six noisy classical readout bits and a separate bare-ancilla Clifford circuit single-fault audit are explicit. Native CNOT compilation, noisy preparation/readout, repeated extraction and multi-fault correlations remain open. No fault-tolerance threshold or physical correction machine follows.')

def produce():
    c=load();p=pin_vacua();print('pin_vacua PASS',flush=True)
    s=spacetime_index();print('spacetime_index PASS',flush=True)
    g=curved_gauge_islands();print('curved_gauge_islands PASS',flush=True)
    u=heavy_extension(c,p);print('heavy_extension PASS',flush=True)
    n=native_noise(c);print('native_noise PASS',flush=True)
    sources=['data/w33_pass11413_11417_pin_regge_vacuum_noise.json','data/w33_pass11403_11407_native_bimodule_collective_gates.json']
    hashes={q:hashlib.sha256(json.dumps(json.loads((ROOT/q).read_text()),sort_keys=True,separators=(',',':')).encode()).hexdigest()for q in sources}
    out=dict(status='PASS',passes='11423-11427',reservation='69e0bfcd1',source_sha256=hashes,source_hash_convention='canonical sorted compact JSON',pin_vacua=p,spacetime_index=s,curved_gauge_islands=g,heavy_extension=u,native_noise=n,
        scope='Five explicit conditional extensions; mathematical indices and numerical models do not constitute a physical TOE.')
    OUT.write_text(json.dumps(out,indent=2)+'\n');return out
if __name__=='__main__':produce()
