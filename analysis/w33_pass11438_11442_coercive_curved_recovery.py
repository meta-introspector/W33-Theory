"""Five explicit continuations, with changed actions and imported theorems named.

Prior owners:11384 native CP;11428 stiff alignment;11429 local measure;
11430 restricted gravity;11431 thresholds;11423 erasure ML;11432 flags.
"""
import json
import hashlib
from pathlib import Path
from itertools import combinations, product
import numpy as np
from scipy.linalg import block_diag
import w33_pass11428_11432_native_alignment_measure_coarse_flags as M

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/w33_pass11438_11442_coercive_curved_recovery.json'


def prior():
    return json.loads((ROOT / 'data/w33_pass11428_11432_native_alignment_measure_coarse_flags.json').read_text())


def bounded_cross(phi, psi, sigma=0.1):
    g = (np.linalg.norm(phi)**2 + np.linalg.norm(psi)**2)/2 + sigma**2
    C = phi.T @ psi.conj() / g
    U = 2*np.eye(3)+(C+C.conj().T)/2
    B = np.eye(3)+C@C.conj().T
    W = U@B-B@U
    return C, float(np.imag(np.trace(W@W@W)))


def finite_native_vacuum(candidate=None, hessian=None):
    import w33_pass11438_finite_native_model as F
    from scipy.optimize import minimize
    p=prior()['native_alignment'];phi=M.P.read(p['phi']).ravel();psi=M.P.read(p['psi']).ravel()
    with F.native_context():
        if candidate is None:
            sol=minimize(F.fun,F.pack(phi,psi),jac=True,method='L-BFGS-B',
                         options=dict(maxiter=400,gtol=1e-8,ftol=1e-14,maxls=30,maxcor=30))
            candidate=sol.x
        x=np.asarray(candidate);value,gradient=F.fun(x)
        assert np.linalg.norm(gradient)<5e-6
        assert value < -5/8 # Below every tilt-zero state: global minimizers have chi*q nonzero.
        if hessian is None:
            step=2e-5;eye=np.eye(324)
            hessian=np.column_stack([(F.fun(x+step*e)[1]-F.fun(x-step*e)[1])/(2*step) for e in eye])
        H=(hessian+hessian.T)/2
        phi=x[:81]+1j*x[81:162];psi=x[162:243]+1j*x[243:]
        ap=1j*(F.H@phi);aq=1j*(F.H@psi)
        orbit=np.concatenate([ap.real,ap.imag,aq.real,aq.imag],axis=1).T
        rank=int(np.linalg.matrix_rank(orbit,tol=1e-8));assert rank==86
        Q=np.linalg.svd(orbit,full_matrices=True)[0][:,rank:]
        eig,vec=np.linalg.eigh(Q.T@H@Q)
        # Independent smaller-step directional Hessian controls on the lowest modes.
        errors=[]
        for j in range(6):
            v=Q@vec[:,j];h=1e-5
            replay=(F.fun(x+h*v)[1]-F.fun(x-h*v)[1])/(2*h)
            errors.append(float(np.linalg.norm(replay-H@v)))
        ward=float(np.linalg.norm(H@orbit)/(np.linalg.norm(H)*np.linalg.norm(orbit)))
        a=F.native(phi)[2];b=F.native(psi)[2];C,q=bounded_cross(phi.reshape(27,3),psi.reshape(27,3))
        soft=[]
        for j in range(4):
            v=Q@vec[:,j];rows=[]
            for amplitude in [.01,.03,.1]:
                increase=(F.fun(x+amplitude*v)[0]+F.fun(x-amplitude*v)[0]-2*value)/2
                rows.append(dict(amplitude=amplitude,symmetric_energy_increase=float(increase)))
            soft.append(dict(mode=j,unrelaxed_line_scans=rows))
        return dict(status='PASS',coordinates=x.tolist(),full324_Hessian=H.tolist(),gradient=gradient.tolist(),
                    action_value=value,gradient_norm=float(np.linalg.norm(gradient)),gauge_orbit_rank=rank,
                    normal_eigenvalues=eig.tolist(),lowest_normal_direction_errors=errors,gauge_Hessian_residual=ward,
                    supplied_coefficients=dict(kappa=F.kappa,rho=F.rho,sigma=F.sigma),
                    I6_phi=[float(a.real),float(a.imag)],I6_psi=[float(b.real),float(b.imag)],CP_flux=q,
                    global_CP_flux_existence='Each native potential is at least -5/16; positive rho adds a nonnegative term. Every chi*q=0 state has V>=-5/8. This trial energy is lower, so every finite global minimizer has chi*q nonzero. CP is the named conjugation combined with compact gauge transformations, not a classification of other generalizedCP symmetries.',
                    positive_normal_mode_count=int(sum(eig>1e-3)),unresolved_soft_normal_mode_count=4,soft_direction_probes=soft,
                    scope='Numerical finite324-field stationary candidate for the separately supplied coercive action.234 normal modes are positive; four near-zero modes have unrelaxed positive approximately quartic line responses. These do not prove coupled quartic stability after massive-field relaxation. No isolated or global minimum identification, interval certificate, physical coefficient selection or observed flavor fit follows.')


def coercive_alignment():
    p = prior()['native_alignment']
    phi, psi = M.P.read(p['phi']), M.P.read(p['psi'])
    chi, q = p['chi'], p['CP_flux']
    a = (-1+1j*np.sqrt(7))/4
    rows = []
    for t in [1., .5, .25, .125, .0625]:
        i6 = t**6*a
        native = abs(i6)**4-abs(i6)**2+i6.real**2+i6.real/2-5/16
        qt = M.cross_pin(t*phi, psi)[3]
        assert abs(qt*t**9-q) < 1e-11
        _, qb = bounded_cross(t*phi, psi)
        rows.append(dict(t=t, native=native, old_alignment=-chi*q/t**3,
                         old_total=native-chi*q/t**3,
                         bounded_alignment=-i6.imag/(1+i6.imag**2)*qb))
    # Native scalar potentials are >=-5/16 per field by11384 squares.
    # ||C||<=sqrt(n1*n2)/((n1+n2)/2+sigma²)<=1.
    # ||[U,B]||<=2, |Im tr[U,B]^3|<=24 and |chi/(1+chi²)|<=1/2.
    C, qb = bounded_cross(phi, psi)
    assert np.linalg.norm(C, 2) <= 1 and qb > 0
    H = M.compact_generators(M.native_tensors()[0])
    acts = np.concatenate([np.einsum('aij,j->ai', H, phi.ravel()),
                           np.einsum('aij,j->ai', H, psi.ravel())], axis=1)
    common = np.concatenate([-acts.imag, acts.real], axis=1).T
    rank = int(np.linalg.matrix_rank(common, tol=1e-9))
    assert rank == 86
    return dict(status='PASS', scans=rows, common_gauge_orbit_rank=rank,
                separate_orbit_dimension=156, relative_orbit_dimension=156-rank,
                runaway='Phi=t Phi0, Psi=Psi0: ImI6=t^6 chi0, q=t^-9 q0, alignment=-kappa chi0 q0/t^3; native potentials stay finite.',
                finite_stiffness_old_action='Unbounded below for every finite stiffness and kappa>0. A global vacuum/Hessian cannot be supplied for this action.',
                replacement='g=(||Phi||F²+||Psi||F²)/2+sigma²; C=Phi^T Psi*/g; U=2I+Re_H C; B=I+CCdagger; V=sum Vnative-rk chi/(1+chi²) q +rho (||Phi||F²+||Psi||F²)², rk=kappa>0,rho>0,sigma>0 supplied.',
                lower_bound='V>=-5/8-12*kappa+rho (||Phi||F²+||Psi||F²)²',
                bounded_seed_flux=qb, normalized_cross=M.P.enc(C),
                result='The replacement is continuous, gauge invariant, CP even and coercive, hence has a finite global minimizer. A numerical finite-field stationary witness and complete Hessian are separately supplied below; global minimizing status is not established.70 relative orbit directions prevent treating the old84-normal single-field audit as a two-field stability result.')


def plaquettes(A, L):
    A = np.asarray(A).reshape((L,)*4+(4,))
    return np.array([A[...,mu]+np.roll(A[...,nu],-1,axis=mu)
                     -np.roll(A[...,mu],-1,axis=nu)-A[...,nu]
                     for mu,nu in combinations(range(4),2)])


def uniform_flux_links(L, m=1):
    x,y = np.indices((L,L)); b=2*np.pi*m/L**2
    U0=np.exp(-1j*b*y)
    U1=np.where(y==L-1,np.exp(1j*b*L*x),1+0j)
    return U0,U1


def admissible_measure():
    p=prior()['local_chiral_measure']; charges=p['integer_hypercharges']; mult=p['multiplicities']
    bound=1/30; epsilon=1e-4; bg=np.array(p['background']); F=epsilon*plaquettes(bg,3)
    maximum=max(abs(q)*np.max(abs(F)) for q in charges)
    assert maximum<bound
    odd={str(e):sum(m for q,m in zip(charges,mult) if abs(q)==e)
         for e in sorted(set(abs(q) for q in charges)) if e%2}
    assert all(m%2==0 for m in odd.values())
    minL=int(np.floor(np.sqrt(2*np.pi*max(abs(q) for q in charges)/bound)))+1
    U0,U1=uniform_flux_links(minL)
    pl=U0*np.roll(U1,-1,axis=0)*np.roll(U0.conj(),-1,axis=1)*U1.conj()
    flux=float(np.angle(pl).sum()/(2*np.pi))
    err=float(np.max(abs(pl-np.exp(2j*np.pi/minL**2))))
    assert minL==34 and abs(flux-1)<1e-12 and err<1e-12
    gaps=[]
    for q in charges:
        proj,_,_,gap,_,_=M.weak_overlap(q,epsilon,bg,[])
        gaps.append(dict(charge=q,gap=gap,rank=int(round(np.trace(proj).real))))
    return dict(status='PASS', one_family_charges=charges,multiplicities=mult,
                linear_anomaly=sum(q*m for q,m in zip(charges,mult)),cubic_anomaly=sum(q**3*m for q,m in zip(charges,mult)),
                odd_absolute_charge_multiplicities=odd, base_epsilon=epsilon,
                maximum_charged_plaquette_phase=float(maximum), locality_sufficient_bound=bound,
                charged_Wilson_gaps=gaps, minimum_L_for_unit_flux=minL,
                unit_flux_plaquette_residual=err, unit_flux=flux,
                unit_flux_maximum_charge_phase=6*2*np.pi/minL**2,
                three_site_topology='Under the common charged admissibility bound, L=3 cannot support nonzero magnetic flux: L²/(30*6)<2pi. The previous tiny-lattice patch cannot test sector gluing.',
                imported_theorem='Luscher hep-lat/9811032: cubic anomaly cancellation plus even multiplicity at each odd absolute charge gives the all-sector Abelian construction on sufficiently large admissible lattices. Bound |q|epsilon<1/30 and finite-size/localization hypotheses are separate. This is prior literature, not a new construction.',
                scope='Charge and admissibility hypotheses and a nontrivial periodic flux sector are explicit. L34 is only the flux/admissibility size bound, not a proof of the theorem finite-size/localization requirement. No implemented global measure current, non-Abelian measure or mirror-free native continuum theory is claimed.')


def spherical_dihedrals(lengths, vertices, kappa):
    G=np.cos(np.sqrt(kappa)*lengths[np.ix_(vertices,vertices)])
    assert np.linalg.eigvalsh(G)[0]>0
    inv=np.linalg.inv(G); angles={}
    for i,j in combinations(range(len(vertices)),2):
        hinge=tuple(sorted(v for k,v in enumerate(vertices) if k not in [i,j]))
        angles[hinge]=float(np.arccos(np.clip(-inv[i,j]/np.sqrt(inv[i,i]*inv[j,j]),-1,1)))
    return angles


def spherical_refinement(boundary, kappa, weights):
    coarse=M.refined_geometry(boundary)[0]; lengths=np.zeros((9,9));lengths[:6,:6]=boundary;fine=[]
    for i,s in enumerate(coarse):
        G=np.cos(np.sqrt(kappa)*boundary[np.ix_(s,s)])
        X=np.linalg.cholesky(G); c=weights[i]@X; c/=np.linalg.norm(c)
        r=np.arccos(np.clip(X@c,-1,1))/np.sqrt(kappa)
        lengths[6+i,list(s)]=r;lengths[list(s),6+i]=r
        fine.extend(tuple([6+i]+[v for v in s if v!=omit]) for omit in s)
    angles={}
    for s in fine:
        for t,a in spherical_dihedrals(lengths,list(s),kappa).items():
            angles[t]=angles.get(t,0)+a
    deficits={t:2*np.pi-a for t,a in angles.items() if any(v>=6 for v in t)}
    return lengths,fine,deficits


def curved_radial_hessian(lengths, fine, hinges, edges, kappa, h):
    def measures(L):
        angles={}
        for s in fine:
            for t,a in spherical_dihedrals(L,list(s),kappa).items():
                angles[t]=angles.get(t,0)+a
        areas=[]
        for t in hinges:
            G=np.cos(np.sqrt(kappa)*L[np.ix_(t,t)]); inv=np.linalg.inv(G)
            aa=sum(np.arccos(np.clip(-inv[i,j]/np.sqrt(inv[i,i]*inv[j,j]),-1,1)) for i,j in combinations(range(3),2))
            areas.append((aa-np.pi)/kappa)
        return np.array(areas),np.array([2*np.pi-angles[t] for t in hinges])
    JA=[];JD=[]
    for e in edges:
        plus=lengths.copy();minus=lengths.copy()
        plus[e]=plus[e[::-1]]=lengths[e]+h;minus[e]=minus[e[::-1]]=lengths[e]-h
        ap,dp=measures(plus);am,dm=measures(minus)
        JA.append((ap-am)/(2*h));JD.append((dp-dm)/(2*h))
    H=np.array(JA)@np.array(JD).T
    return (H+H.T)/2


def curved_stationarity():
    bd=np.array(M.previous()['curved_gauge_islands']['boundary_lengths']);k=.001
    weights=[np.ones((3,5))/5, np.array([[.1,.15,.2,.25,.3],[.3,.2,.1,.15,.25],[.18,.22,.25,.2,.15]])]
    scans=[]
    for W in weights:
        L,fine,d=spherical_refinement(bd,k,W)
        maximum=max(abs(x) for x in d.values());assert maximum<1e-9
        # The exact Schlaefli cancellation leaves dS=sum_h delta_h dA_h.
        # Boundary-only hinges have no radial derivatives; all remaining delta vanish.
        scans.append(dict(weights=W.tolist(),radials=[L[6+i,list(s)].tolist() for i,s in enumerate(M.refined_geometry(bd)[0])],
                          interior_hinge_count=len(d),maximum_interior_deficit=maximum,full_radial_stationarity=True))
    _,_,_,g,_,_=M.normal_stationary(bd,.001)
    L,fine,d=spherical_refinement(bd,k,weights[0])
    edges=[(6+i,v) for i,s in enumerate(M.refined_geometry(bd)[0]) for v in s]
    hessians=[]
    for h in [5e-5,1e-5]:
        H=curved_radial_hessian(L,fine,sorted(d),edges,k,h);ev=np.linalg.eigvalsh(H)
        assert min(ev[-3:])>1000 and max(abs(ev[:12]))/min(ev[-3:])<2e-6
        hessians.append(dict(step=h,Hessian=H.tolist(),eigenvalues=ev.tolist()))
    return dict(status='PASS',boundary_lengths=bd.tolist(),sectional_curvature=k,
                action='S_k=sum_internal A_k delta_k+sum_boundary A_k psi_k+3*k sum V4_k; geodesic curved simplices and boundary term.',
                cosmological_convention='Lambda=3*k in the Regge volume coefficient convention of11430; not a measured cosmological constant.',
                scans=scans,full15_radial_Hessian_controls=hessians,flat_discretization_full_radial_gradient_norm=float(np.linalg.norm(g)),
                theorem='For each geodesic1->5 subdivision, internal hinges tile the same constant-curvature simplex and have zero deficit. Sum A dtheta=3*k dV cancels the volume derivative, so all15 radial equations vanish. Each star has four independent center-displacement null directions;12 total. Action is independent of center placement on this subdivision family.',
                map='Gij=cos(sqrt(k)*ell_ij), X=Cholesky(G), c=normalize(sum wi Xi), ri=acos(Xi dot c)/sqrt(k). Positive barycentric weights specify the interior point.',
                prior='Bahr-Dittrich0907.4323 curved Regge/Schlaefli construction; not invented here.',
                scope='Exact local stationary elimination by replacing flat building blocks on the actual supplied boundary. This closes the radial/centroid equations of this subdivision action, not those of the old flat-simplex action, nor a globally perfect4D gravity action, Lorentzian dynamics or physicalCC selection. Volumes are named geometrically; no numerical4-volume integration or global boundary vacuum is asserted.')


def slice_threshold(phi, mu=1.):
    p=M.previous(); A=M.P.load()['A']; pins=[M.P.read(p['pin_vacua'][key]) for key in ['pin_up','pin_down']]
    masses=[M.full_quark_masses(A,B,phi) for B in pins]
    r=2*.1*(3*phi**2-.8**2);z=2*.1*(phi**2-.8**2)
    assert z>0
    st=-12*sum(np.sum(m**4) for m in masses)+480*(r*r+z*z)
    V=-12*sum(np.sum(m**4*(np.log(m*m/mu**2)-1.5)) for m in masses)
    V+=480*sum(x*x*(np.log(x/mu**2)-1.5) for x in [r,z])
    return float(V/(64*np.pi**2)),float(st),float(480*.1*(phi*phi-.8**2)**2)


def renormalized_slice():
    star=.81;h=.0002;K=100.; target_energy=0.
    rows={x:slice_threshold(x) for x in [star-2*h,star-h,star,star+h,star+2*h,.805,.82]}
    value=lambda x:rows[x][0]+rows[x][2]
    v0=value(star);d1=(value(star+h)-value(star-h))/(2*h)
    d2=(value(star+h)+value(star-h)-2*v0)/h**2
    # Three supplied renormalization conditions: energy, force, curvature.
    # ct=c0+c2 phi²+c4 phi⁴, even local radial slice polynomial.
    c4=(K-d2+d1/star)/(8*star**2);c2=-d1/(2*star)-2*c4*star**2
    c0=target_energy-v0-c2*star**2-c4*star**4
    coef=np.array([c0,c2,c4]);fitphi=np.array([.805,star,.82]);fit=np.array([[1,x*x,x**4] for x in fitphi])
    beta=np.linalg.solve(fit,np.array([rows[x][1] for x in fitphi])/(32*np.pi**2))
    replay=[]
    for x in [star-h,star,star+h]:
        base=rows[x];vec=np.array([1,x*x,x**4]);vals=[]
        for mu in [.5,1.,2.]:
            v,st,tree=slice_threshold(x,mu)
            vals.append(v+tree+vec@(coef+beta*np.log(mu)))
        replay.append(dict(phi=x,renormalized_potential_by_scale=vals,scale_spread=max(vals)-min(vals),supertrace_polynomial_error=float(vec@beta*32*np.pi**2-base[1])))
    assert max(r['scale_spread'] for r in replay)<1e-7
    f=(replay[2]['renormalized_potential_by_scale'][1]-replay[0]['renormalized_potential_by_scale'][1])/(2*h)
    curvature=(replay[2]['renormalized_potential_by_scale'][1]+replay[0]['renormalized_potential_by_scale'][1]-2*replay[1]['renormalized_potential_by_scale'][1])/h**2
    assert abs(f)<.01 and abs(curvature-K)<.05
    return dict(status='PASS',imposed_stationary_phi=star,imposed_radial_curvature=K,imposed_vacuum_energy=target_energy,
                counterterm_coefficients=coef.tolist(),counterterm_beta_coefficients=beta.tolist(),scans=replay,
                residual_force=float(f),measured_radial_curvature=float(curvature),
                action='Vtree(phi)+V1_full_quark_and_link(phi;mu)+c0(mu)+c2(mu)phi²+c4(mu)phi⁴',
                running='dc/dln(mu)=coefficients of STrM4(phi)/(32pi²); frozen spectral masses, one-loop scale cancellation on the declared slice.',
                scope='A self-consistent stationary stable radial slice with explicit one-loop counterterm running, constructed by supplied renormalization conditions. This does not derive phi, energy, curvature, fundamental coupling beta functions, joint pin/Majorana equations, full995-field Hessian or small cosmological constant. The scalar slice has phi>.8 to avoid tachyonic logarithms.')


def steane_frame():
    checks=[sum(1<<j for j in range(7) if ((j+1)>>row)&1) for row in range(3)]
    span={0}
    for c in checks:span|={x^c for x in list(span)}
    V=np.zeros((128,2),complex)
    for x in span:V[x,0]=1/np.sqrt(8);V[x^127,1]=1/np.sqrt(8)
    return V


def reset_kraus(sites):
    # Known erased qubits are replaced by |0>; outcomes are not postselected.
    mask=sum(1<<j for j in sites);out=[]
    for bits in product(range(2),repeat=len(sites)):
        chosen=sum(b<<j for b,j in zip(bits,sites));K=np.zeros((128,128),complex)
        for x in range(128):
            if x&mask==chosen:K[x&~mask,x]=1
        out.append(K)
    return out


def erasure_recovery(sites):
    V=steane_frame();Ks=reset_kraus(sites);r=len(Ks)
    # KL matrix is I/r for up to two known erasures in the distance-three code.
    W=np.hstack([np.sqrt(r)*K@V for K in Ks]);assert np.linalg.norm(W.conj().T@W-np.eye(2*r))<1e-12
    Rs=[np.sqrt(r)*V@(K@V).conj().T for K in Ks]
    complement=np.eye(128)-W@W.conj().T
    # Complete the recovery to a trace-preserving map on the entire128-space.
    ev,Q=np.linalg.eigh(complement)
    # Complement Kraus operators are |0L><v| for each orthonormal
    # complement vector. Sum their adjoint products without materializing128 maps.
    complement_frame=Q[:,ev>.5]
    tp=sum(R.conj().T@R for R in Rs)+complement_frame@complement_frame.conj().T
    assert np.linalg.norm(tp-np.eye(128))<1e-11
    operators=[V.conj().T@R@K@V for R in Rs for K in Ks]
    assert max(np.linalg.norm(complement_frame.conj().T@K@V) for K in Ks)<1e-12
    entanglement_fidelity=sum(abs(np.trace(O))**2 for O in operators)/4
    assert abs(entanglement_fidelity-1)<1e-12
    return float(entanglement_fidelity),float(np.linalg.norm(tp-np.eye(128)))


def leakage_recovery():
    rows=[]
    for size in [1,2]:
        for sites in combinations(range(7),size):
            fidelity,tp=erasure_recovery(sites)
            rows.append(dict(erased_sites=list(sites),entanglement_fidelity=fidelity,TP_residual=tp))
    V=steane_frame();noncorrectable=[]
    # Three erased sites include a minimum logical Pauli support; record the obstruction.
    for sites in combinations(range(7),3):
        mask=sum(1<<j for j in sites);Z=np.diag([(-1)**((x&mask).bit_count()) for x in range(128)])
        A=V.conj().T@Z@V
        if np.linalg.norm(A-np.trace(A)/2*np.eye(2))>.1:noncorrectable.append(list(sites))
    assert len(noncorrectable)==7
    # Native adjacency supplies routes without inventing a degree-eight hub.
    A=M.P.load()['A'];paths=[]
    from collections import deque
    def path(s,t):
        queue=deque([s]);prev={s:None}
        while queue:
            v=queue.popleft()
            if v==t:break
            for w in np.flatnonzero(A[v]):
                w=int(w)
                if w not in prev:prev[w]=v;queue.append(w)
        assert t in prev;out=[t]
        while out[-1]!=s:out.append(prev[out[-1]])
        return out[::-1]
    for data in range(7):paths.append(path(7,data))
    paths.append(path(7,8))
    costs=[6*(len(p)-2)+1 for p in paths] # swap forward, CNOT, swap back; each swap=3 CNOT.
    native_ticks=prior()['flagged_native_correction']['native_CNOT']['CNOT_total_ticks']
    return dict(status='PASS',exact_recovery_supports=rows,uncorrectable_three_erasure_supports=noncorrectable,
                recovery='Detect location, reset erased physical qubits to|0>, apply complete CPTP recovery Rs. Every Kraus branch is retained; no block rejection or postselection.',
                known_single_or_double_erasure_identity=True, native_graph_paths=paths,
                CNOT_cost_per_routed_syndrome_or_flag_pair=costs,native_ticks_per_CNOT=native_ticks,
                routing='An architecture of80 independent12-dimensional cell copies indexed by the native A80 graph; graph-neighbor joint gates are supplied couplings. SWAP moves control to the neighbor of target, CNOT, reverse SWAPs restore spectators. Serial scheduling avoids collisions; each SWAP uses three bidirectional CNOTs.',
                scope='Explicit deterministic recovery map and graph routing recipe, not a native implementation of the128-dimensional recovery channel or a realization of80 independent cells within one cover. Ideal location detection, reset, recovery synthesis, neighbor couplings and ancilla readout/preparation are additional primitives. Existing flags certify one Pauli fault, not arbitrary leakage propagation or routed multi-fault tolerance. No threshold is claimed.')


def produce(candidate=None, hessian=None):
    out=dict(status='PASS',passes='11438-11442',reservation='afa5402d3')
    for name,fn in [('coercive_alignment',coercive_alignment),('admissible_measure',admissible_measure),
                    ('curved_stationarity',curved_stationarity),('renormalized_slice',renormalized_slice),
                    ('leakage_recovery',leakage_recovery)]:
        out[name]=fn()
        if name=='coercive_alignment':out[name]['finite_field_witness']=finite_native_vacuum(candidate,hessian)
        print(name,'PASS',flush=True)
    sources=['data/w33_pass11428_11432_native_alignment_measure_coarse_flags.json',
             'data/w33_pass11423_11427_flux_index_curved_uv_noise.json']
    out['source_sha256']={s:hashlib.sha256(json.dumps(json.loads((ROOT/s).read_text()),sort_keys=True,separators=(',',':')).encode()).hexdigest() for s in sources}
    OUT.write_text(json.dumps(out,indent=2)+'\n')
    return out


if __name__=='__main__':
    produce()
