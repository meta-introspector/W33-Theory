#!/usr/bin/env python3
"""Conditional five-front packet. Prior11408-11412,11271/11274,BT986,
Pass235 (FN ansatz), PART_CCCCIV (Steane). Also prior2D triangle-count
action:analysis/w33_discrete_einstein_hilbert.py, distinct from this4D patch.
No observed-physics prediction.
"""
import json, hashlib
from pathlib import Path
from itertools import combinations, product
import numpy as np
from scipy.linalg import expm, block_diag
from w33_pass11390_11394_context_matter_clock import load
from w33_pass11408_11412_clock_spurion_seam_ticks import (
    prior, read, mixing, unitary_calibration)
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/w33_pass11413_11417_pin_regge_vacuum_noise.json'
def enc(a):
    a=np.asarray(a);return dict(real=a.real.tolist(),imag=a.imag.tolist())

def path_flavour(c):
    A=c['A'];n=len(A);dist=np.full(n,-1,int);dist[0]=0;queue=[0]
    for i in queue:
        for j in np.flatnonzero(A[i]):
            if dist[j]<0:dist[j]=dist[i]+1;queue.append(j)
    idx=[int(np.flatnonzero(dist==d)[0])for d in [3,2,0]]
    paths=[int(np.linalg.matrix_power(A,d)[i,0])for i,d in zip(idx,[3,2,0])]
    assert paths==[1,1,1]
    assert np.bincount(dist).tolist()==[1,4,12,36,27]
    X=np.roll(np.eye(3),1,axis=0);z=.2+.15j
    Bu=np.diag([1.,2.,3.]);Bd=np.diag([2.,3.,4.])+z*X+z.conjugate()*X.T
    assert np.linalg.eigvalsh(Bd)[0]>0
    P=np.zeros((n,n));P[0,0]=1;rows=[];critical=[]
    # Exact radial resolvent from intersection array{4,3,3,3;1,1,1,4}.
    for eta in [.04,.15,.225,.2499]:
        den=(1-16*eta**2)*(1-6*eta**2)
        radial=np.array([1-18*eta**2+36*eta**4,eta*(1-15*eta**2),eta**2*(1-12*eta**2),eta**3,4*eta**4])/den
        R=np.linalg.inv(np.eye(n)-eta*A)
        assert np.linalg.norm(R[:,0]-radial[dist])<1e-9
        critical.append(dict(eta=eta,radial_resolvent=radial.tolist(),endpoint_to_pin_ratios=(radial[[3,2,0]]/radial[0]).tolist()))
    for eta in [.08,.04,.02,.01,.005]:
        K=np.eye(n)-eta*A;R=np.linalg.inv(K);r=R[idx,0];L=np.diag(r)
        Yu=L@Bu@L;Yd=L@Bd@L
        su,U=np.linalg.eigh(Yu);sd,V=np.linalg.eigh(Yd);m=mixing(U.conj().T@V)
        ang=np.radians(m['angles_degrees']);comm=Yu@Yu@Yd@Yd-Yd@Yd@Yu@Yu
        cp=float(np.imag(np.trace(comm@comm@comm)))
        # Full480-dimensional heavy neutral-family block, not scalar-pin rank1.
        if eta==.04:
            k=np.kron(K,np.eye(3));bridge=np.kron(P,Bd)
            M=np.block([[k,bridge],[np.zeros_like(k),k]])
            right=np.zeros((2*n*3,3))
            for j,i in enumerate(idx):right[n*3+i*3+j,j]=1
            solution=np.linalg.solve(M,right)
            matched=-np.array([solution[i*3+j,:]for j,i in enumerate(idx)])
            assert np.linalg.norm(matched-Yd)<1e-12
            rng=np.random.default_rng(11413);G=np.diag(np.exp(1j*rng.normal(size=n)))
            gauge=np.linalg.inv(G@K@G.conj().T)
            assert np.linalg.norm(gauge-G@R@G.conj().T)<1e-12
        assert abs(np.linalg.det(Yd)-np.linalg.det(Bd)*np.prod(r)**2)<1e-20
        rows.append(dict(eta=eta,resolvent_endpoints=r.tolist(),Yu=enc(Yu),Yd=enc(Yd),
            up_masses=su.tolist(),down_masses=sd.tolist(),cubic_CP=cp,
            angles_over_eta_powers=(ang/np.array([eta,eta**2,eta**3])).tolist(),
            J_over_eta6=m['J']/eta**6,**m))
    powers=np.log(np.array(rows[-2]['down_masses'])/rows[-1]['down_masses'])/np.log(2)
    assert np.max(abs(powers-[6,4,0]))<.02
    schur=Bd[:2,:2]-np.outer(Bd[:2,2],Bd[2,:2])/Bd[2,2]
    mass_leading=[np.linalg.det(Bd)/np.linalg.det(Bd[1:,1:]),schur[1,1],Bd[2,2]]
    angle_leading=[abs(schur[0,1])/schur[1,1].real,abs(Bd[1,2])/Bd[2,2].real,abs(Bd[0,2])/Bd[2,2].real]
    J_leading=np.imag(Bd[0,1]*Bd[1,2]*Bd[2,0])/(Bd[2,2].real**2*schur[1,1].real)
    assert max(abs(np.array(rows[-1]['angles_over_eta_powers'])/angle_leading-1))<.01
    assert abs(abs(rows[-1]['J_over_eta6']/J_leading)-1)<.01
    # Native colour commutation without changing family/colour ordering.
    slots=prior()['fermion_bimodule']['native_slot_frames'];F=[read(x)for x in slots['colour_qutrit_seeds']];C=np.hstack(F)
    from w33_pass11403_11407_native_bimodule_collective_gates import gellmann
    native=C@np.kron(Yd,np.eye(3))@C.conj().T
    ward=max(np.linalg.norm(native@G-G@native)for G in [sum(f@t@f.conj().T for f in F)for t in gellmann()])
    assert ward<1e-12
    return dict(status='PASS',endpoints=idx,distances=[3,2,0],shortest_path_counts=paths,
        pin_up=enc(Bu),pin_down=enc(Bd),heavy_family_dimension=480,
        distance_shell_sizes=[1,4,12,36,27],radial_resolvent='[1-18eta²+36eta4,eta(1-15eta²),eta²(1-12eta²),eta³,4eta4]/[(1-16eta²)(1-6eta²)]',critical_scans=critical,
        critical_limit='As eta approaches1/4 from below, every radial endpoint/pin ratio tends to1: the Perron pole erases path hierarchy. Small-eta powers cannot be extrapolated as a Cabibbo fit.',
        matching_identity='Y=D_R B D_R; M=[[K tensor I3,P_pin tensor B],[0,K tensor I3]], K=I-eta A_Levi',
        mass_orders=[6,4,0],mixing_orders=[1,2,3],J_order=6,measured_mass_orders=powers.tolist(),
        asymptotic_mass_coefficients=np.real(mass_leading).tolist(),asymptotic_angle_coefficients=angle_leading,asymptotic_J_magnitude_coefficient=float(abs(J_leading)),
        colour_Ward_residual=float(ward),scans=rows,
        scope='Full-rank family fibre at the pin preserves link power d_i+d_j while allowing noncommuting up/down pin matrices. Graph, endpoint and pin choices, link spurions, eta and Higgs insertions are supplied. Perturbative analytic spurion rule, not anomaly-certified protection, unique vacuum or measured CKM prediction. Prior FN mechanism and11410 rank-one obstruction are cited.')

def symmetric_basis():
    B=[]
    for i,j in combinations(range(3),2):
        E=np.zeros((3,3),complex);E[i,j]=E[j,i]=1/np.sqrt(2);B.extend([E,1j*E])
    for i in range(3):
        E=np.zeros((3,3),complex);E[i,i]=1;B.extend([E,1j*E])
    return B

def majorana_vacuum(f):
    S0=np.diag([1.,2.,3.]).astype(complex);t,lam,kappa,v,lh=.3,.2,.1,1.,.13
    B=symmetric_basis();J=[E@S0.conj().T+S0@E.conj().T for E in B]
    H=2*(t+kappa*v*v/2)*np.eye(12)+2*lam*np.array([[np.real(np.vdot(a,b))for b in J]for a in J])
    assert np.linalg.eigvalsh(H)[0]>0
    rng=np.random.default_rng(11414);res=[]
    for _ in range(20):
        x=rng.normal(size=12)*.1;S=S0+sum(q*e for q,e in zip(x,B));h=rng.normal(size=4)
        delta=S-S0;norm=float(np.linalg.norm(delta)**2);hh=np.dot(h,h)/2
        pot=lh*(hh-v*v/2)**2+t*norm+lam*np.linalg.norm(S@S.conj().T-S0@S0.conj().T)**2+kappa*hh*norm
        assert pot>=0;res.append(float(pot))
    chars=[(2,0),(2,2),(1,1)];channels={}
    for i in range(3):
        for j in range(i,3):channels[f'{i}{j}']=[(-chars[i][a]-chars[j][a])%3 for a in range(2)]
    # Declared lepton attachments sit at the pin, unlike the quark endpoints.
    # This keeps the full seesaw numerically resolvable; no quark-lepton fit.
    Y=.01*read(f['pin_down']);N=np.block([[np.zeros((3,3)),v*Y/np.sqrt(2)],[v*Y.T/np.sqrt(2),S0]])
    masses=np.linalg.svd(N,compute_uv=False);assert masses[-1]>0
    # Any square finite Dirac block with balanced chirality has index zero.
    D=np.block([[np.zeros_like(N),N],[N.conj().T,np.zeros_like(N)]])
    gamma=np.diag([1]*6+[-1]*6)
    assert np.linalg.norm(D@gamma+gamma@D)<1e-12
    return dict(status='PASS',Majorana_source=enc(S0),symmetric_field_basis=[enc(b)for b in B],
        parameters=dict(t=t,lambda_S=lam,kappa=kappa,v=v,lambda_H=lh),
        potential='lambda_H (h^dag h-v²/2)²+t||S-S0||²+lambda_S||SS^dag-S0S0^dag||²+kappa h^dag h||S-S0||²',
        scalar_Hessian=H.tolist(),scalar_mass_squared=np.linalg.eigvalsh(H).tolist(),Higgs_mass_squared=2*lh*v*v,
        Goldstone_count=3,required_character_channels=channels,neutrino_Dirac=enc(Y),neutrino_attachment='all three lepton fibres at the pin; Ynu=.01 B_down is supplied, not the quark distance texture',
        seesaw_matrix=enc(N),neutrino_masses=masses.tolist(),balanced_grading_index=0,
        scope='An explicit globally nonnegative polynomial action has S=S0 and h^dag h=v²/2 as its full-rank minimum;12 singlet directions strictly stable and3 electroweak Goldstones. S0 is an explicit H27-breaking source and coefficients are supplied. Bare neutrality assumes h_charge=3q. This constructs a conditional vacuum, not a source-free hierarchy or physical chirality selector; finite balanced grading remains index zero.')

def simplex_geometry(vertices, lengths):
    # Gram at first vertex; barycentric gradients are outward normal negatives.
    o=vertices[0];vs=vertices[1:]
    G=np.array([[(lengths[o,i]**2+lengths[o,j]**2-lengths[i,j]**2)/2 for j in vs]for i in vs])
    assert np.linalg.eigvalsh(G)[0]>0
    inv=np.linalg.inv(G);grad=np.vstack([-np.ones(4),np.eye(4)])
    normals=grad@inv@grad.T
    angles={}
    for a,b in combinations(range(5),2):
        hinge=tuple(sorted(set(vertices)-{vertices[a],vertices[b]}))
        cos=-normals[a,b]/np.sqrt(normals[a,a]*normals[b,b])
        angles[hinge]=np.arccos(np.clip(cos,-1,1))
    return angles

def regge_action(radial,boundary):
    ell=np.zeros((6,6));ell[1:,1:]=boundary;ell[0,1:]=ell[1:,0]=radial
    angle={}
    for omitted in range(1,6):
        vertices=[0]+[j for j in range(1,6)if j!=omitted]
        for t,x in simplex_geometry(vertices,ell).items():angle[t]=angle.get(t,0)+x
    action=0.;deficits=[]
    for tri,total in angle.items():
        a,b,c=[ell[i,j]for i,j in combinations(tri,2)];s=(a+b+c)/2
        area=np.sqrt(max(0,s*(s-a)*(s-b)*(s-c)))
        deficit=(2*np.pi if 0 in tri else np.pi)-total
        action+=area*deficit
        if 0 in tri:deficits.append(deficit)
    return float(action),np.array(deficits)

def hessian(fn,x,h=1e-4):
    n=len(x);H=np.zeros((n,n));f=fn(x);E=np.eye(n)*h
    for i in range(n):
        H[i,i]=(fn(x+E[i])+fn(x-E[i])-2*f)/h**2
        for j in range(i):H[i,j]=H[j,i]=(fn(x+E[i]+E[j])-fn(x+E[i]-E[j])-fn(x-E[i]+E[j])+fn(x-E[i]-E[j]))/(4*h*h)
    return H

def regge_patch():
    old=json.loads((ROOT/'data/w33_pass11408_11412_clock_spurion_seam_ticks.json').read_text())['coherent_seams']
    F=np.array(old['FCC_primitive'],float);M=np.array(old['coincidence_left_periods'],float);period=F@M
    # Five boundary vertices: spatial native supercell edges and supplied time.
    V=np.zeros((5,4));V[1:4,:3]=period.T;V[4,3]=2.
    boundary=np.linalg.norm(V[:,None,:]-V[None,:,:],axis=2);center=V.mean(axis=0)
    radial=np.linalg.norm(V-center,axis=1);T=(center-V)/radial[:,None]
    fn=lambda r:regge_action(r,boundary)[0];base,deficit=regge_action(radial,boundary)
    assert max(abs(deficit))<1e-12
    moved=[]
    for delta in [[.01,0,0,0],[0,-.01,.01,.01],[.02,.01,-.01,-.02]]:
        r=np.linalg.norm(V-center-np.array(delta),axis=1);s,d=regge_action(r,boundary)
        moved.append(dict(shift=delta,action_change=s-base,max_deficit=float(max(abs(d)))))
        assert abs(s-base)<1e-10 and max(abs(d))<1e-12
    H=hessian(fn,radial,h=5e-5);res=np.linalg.norm(H@T);w=np.linalg.eigvalsh(H)
    # Finite differences are numerical diagnostics, not exact rank proof.
    assert res<.02 and np.count_nonzero(abs(w)>.02)==1
    # Normal to the4-dimensional exact flat-embedding manifold.
    _,_,vh=np.linalg.svd(T.T);normal=vh[-1]
    curved=[]
    for eps in [.002,.01,.03]:
        r=radial+eps*normal;s,d=regge_action(r,boundary);HH=hessian(fn,r)
        curved.append(dict(epsilon=eps,max_deficit=float(max(abs(d))),action_change=s-base,
            Hessian_eigenvalues=np.linalg.eigvalsh(HH).tolist(),flat_tangent_Hessian_residual=float(np.linalg.norm(HH@T))))
    return dict(status='PASS',spatial_supercell_periods=period.tolist(),boundary_vertices=V.tolist(),
        boundary_lengths=boundary.tolist(),interior_vertex=center.tolist(),radial_lengths=radial.tolist(),
        action='sum_internal A(2pi-sum theta)+sum_boundary A(pi-sum theta); Euclidean length Regge, G and Lambda supplied/omitted',
        flat_action=base,flat_Hessian=H.tolist(),flat_Hessian_eigenvalues=w.tolist(),
        vertex_translation_tangents=T.tolist(),finite_difference_step=5e-5,flat_tangent_residual=float(res),
        nonlinear_flat_moves=moved,curved_offshell_scans=curved,
        scope='Actual4D1-to5 simplex patch uses common native FCC supercell periods and a supplied time direction. Four finite vertex motions preserve the flat action; numerical Hessian has4 near-zero modes and1 curvature mode. Generic off-shell curvature lifts flat tangents. This is standard Regge symmetry on a new native-period witness, not Lorentzian dynamics, a closed nonlinear constraint algebra or ghost freedom.')

def vacuum_inventory(f,m):
    v=m['parameters']['v'];g,gp=.65,.35;mu=1.
    row=f['scans'][1];Yu=read(row['Yu']);Yd=read(row['Yd']);Ye=.03*np.diag([.01,.1,1.])
    groups=[]
    for name,Y,mult in [('up',Yu,3),('down',Yd,3),('charged_lepton',Ye,1)]:
        masses=np.linalg.svd(v*Y/np.sqrt(2),compute_uv=False)
        groups.extend(dict(name=f'{name}{i}',spin='Dirac',degeneracy=-4*mult,m2=float(x*x),c=1.5)for i,x in enumerate(masses))
    groups.extend(dict(name=f'neutrino{i}',spin='Majorana',degeneracy=-2,m2=x*x,c=1.5)for i,x in enumerate(m['neutrino_masses']))
    groups.extend(dict(name=f'S{i}',spin='real_scalar',degeneracy=1,m2=x,c=1.5)for i,x in enumerate(m['scalar_mass_squared']))
    groups.append(dict(name='Higgs',spin='real_scalar',degeneracy=1,m2=m['Higgs_mass_squared'],c=1.5))
    groups.extend([dict(name='W_pair',spin='vector',degeneracy=6,m2=g*g*v*v/4,c=5/6),dict(name='Z',spin='vector',degeneracy=3,m2=(g*g+gp*gp)*v*v/4,c=5/6)])
    for r in groups:r['CW']=r['degeneracy']*r['m2']**2*(np.log(r['m2']/mu**2)-r['c'])/(64*np.pi**2)
    str4=sum(r['degeneracy']*r['m2']**2 for r in groups);cw=sum(r['CW']for r in groups)
    assert sum(-r['degeneracy']for r in groups if r['spin']in ['Dirac','Majorana'])==96
    # Explicit RG identity, holding masses fixed: dV1/dln(mu)=-STr M4/(32pi²).
    cw2=sum(r['degeneracy']*r['m2']**2*(np.log(r['m2']/4)-r['c'])/(64*np.pi**2)for r in groups)
    assert abs(cw2-cw+str4*np.log(2)/(32*np.pi**2))<1e-12
    return dict(status='PASS',renormalization='Landau gauge MSbar, mu=1, fixed background and masses',
        input_gauge_couplings=[g,gp],inventory=groups,one_loop_value=cw,supertrace_M4=str4,
        fixed_mass_log_scale_derivative=-str4/(32*np.pi**2),real_scalar_count=16,physical_scalar_count=13,
        Weyl_count=48,gauge_generator_count=12,massive_vector_count=3,massless_vector_count=9,
        zero_mass_landau_terms='3 Goldstones, massless gluons/photon and Landau ghosts contribute zero to V1 at this stationary background',
        scope='Complete gauge/scalar/fermion determinant for the explicitly declared48-Weyl low-energy model and12-real-singlet/Higgs potential. Heavy messenger thresholds are integrated out and not supplied; their dynamical link scalars and UV spectrum are not counted. Nonzero supertrace and a freely renormalized vacuum counterterm leave the cosmological constant open. No quantum-stable vacuum or gravitational cancellation follows.')

def steane_dephasing(p):
    H=np.array([[((j>>i)&1)for j in range(1,8)]for i in range(3)])
    # Kernel H, syndrome minimum-weight recovery. Logical Z iff residual odd.
    counts=np.zeros(8,int)
    for pattern in product([0,1],repeat=7):
        e=np.array(pattern);s=H@e%2;r=e.copy()
        if s.any():r[np.flatnonzero(np.all(H.T==s,axis=1))[0]]^=1
        assert not (H@r%2).any()
        if r.sum()%2:counts[e.sum()]+=1
    fail=float(sum(n*p**w*(1-p)**(7-w)for w,n in enumerate(counts)))
    assert counts[0]==counts[1]==0 and counts[2]==21
    return fail,counts.tolist()

def counted_noise(c):
    p=prior();old=json.loads((ROOT/'data/w33_pass11408_11412_clock_spurion_seam_ticks.json').read_text())['finite_ticks']
    D=np.array(c['D'],float);Hp=D[:40].T@D[:40];Hl=D[40:].T@D[40:]
    Q=read(p['local_gates']['logical_basis']);a,b=p['local_gates']['gate_pair'];target=read(p['local_gates']['ideal_qubit_gate'])
    phases=np.array(p['local_gates']['five_pulse_phase_angles']);short=(phases+np.pi)%(2*np.pi)-np.pi
    dt,u=.04,.0125;clock=expm(-.5j*dt*Hp)@expm(-1j*dt*Hl)@expm(-.5j*dt*Hp)
    P=json.loads((ROOT/'data/w33_pass11395_11402_joint_atlas_spin_portal.json').read_text())['axis_bridge']['cycle_projector_160_numerator'];P=np.array(P,dtype=float)/160
    state=Q.copy();counts=[];bound=0
    for angle,e in zip(short,[a,b,a,b,a]):
        kick=np.ones(160,complex);kick[e]=np.exp(-.5j*dt*np.copysign(u,angle));S=kick[:,None]*clock*kick[None,:]
        axis=P[:,e]/np.sqrt(P[e,e]);theta,V,N,B,_=unitary_calibration(S,axis,angle)
        state=V@(np.exp(-1j*N*theta)[:,None]*(V.conj().T@state));counts.append(N);bound+=B
    err=float(np.linalg.norm(state-Q@target,ord=2));assert err<bound+1e-9
    long=next(r for r in old['Hadamard_scans']if r['dt']==dt and r['strength']==u)['total_ticks']
    N=sum(counts);assert N<long
    noise=[]
    for q in [1e-9,1e-8,1e-7,1e-6]:
        prob=-np.expm1(N*np.log1p(-2*q))/2;failure,coeffs=steane_dephasing(prob)
        noise.append(dict(probability_per_idle_tick=q,ticks=N,effective_Z_probability=float(prob),
            encoded_logical_Z_probability=failure,improves=bool(failure<prob)))
    return dict(status='PASS',original_angles=phases.tolist(),principal_angles=short.tolist(),
        ticks_per_pulse=counts,total_ticks=N,previous_ticks=long,overhead_ratio=long/N,
        coherent_encoded_error=err,coherent_error_bound=bound,Steane_failure_weight_counts=coeffs,noise_scans=noise,
        scope='Exact2pi projector-phase periodicity shortens the prior native5-pulse Hadamard, replayed on160 edges. QEC is an explicit independent post-gate/idle logical-Z stochastic model with perfect Steane syndrome/recovery; these Z channels commute, so tick aggregation is exact for idles. It is not native during-gate edge noise, leakage correction, a physical syndrome circuit or a fault-tolerance threshold. Prior PART_CCCCIV owns the Steane lift.')

def produce():
    c=load();f=path_flavour(c);print('path_flavour PASS',flush=True)
    m=majorana_vacuum(f);print('majorana_vacuum PASS',flush=True)
    r=regge_patch();print('regge_patch PASS',flush=True)
    v=vacuum_inventory(f,m);print('vacuum_inventory PASS',flush=True)
    q=counted_noise(c);print('counted_noise PASS',flush=True)
    sources=['data/w33_pass11403_11407_native_bimodule_collective_gates.json','data/w33_pass11408_11412_clock_spurion_seam_ticks.json','data/w33_pass11389_parabolic_spatial_cover.json','data/w33_pass11395_11402_joint_atlas_spin_portal.json']
    hashes={s:hashlib.sha256(json.dumps(json.loads((ROOT/s).read_text()),sort_keys=True,separators=(',',':')).encode()).hexdigest()for s in sources}
    out=dict(status='PASS',passes='11413-11417',reservation='77dc3d5fb',source_sha256=hashes,source_hash_convention='canonical sorted compact JSON',path_flavour=f,majorana_vacuum=m,regge_patch=r,vacuum_inventory=v,counted_noise=q,
        scope='Five conditional constructed witnesses; no complete TOE or measured-parameter derivation.')
    OUT.write_text(json.dumps(out,indent=2)+'\n');return out
if __name__=='__main__':produce()
