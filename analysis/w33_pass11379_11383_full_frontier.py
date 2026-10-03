"""Full ray normal stability, native cuts and rotations, and wound-wall proofs.

Exact statements and floating-point native tensor witnesses are labelled separately.
The actions and couplings remain supplied; no observed-physics fit is claimed.
"""
from pathlib import Path
import json
import sys
import numpy as np
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'analysis'))


def ray_data():
    z = s.symbols('z:18', real=True)
    X = s.Matrix(3, 3, lambda i,j: z[2*(3*i+j)] + s.I*z[2*(3*i+j)+1])
    G = X.conjugate().T*X
    B = G[0,1]*G[1,2]*G[2,0]
    f = s.Matrix([G[i,i]-1 for i in range(3)] +
                 [s.expand(G[i,j]*G[j,i])-t for i,j,t in
                  [(0,1,s.Rational(1,2)),(0,2,s.Rational(1,3)),(1,2,s.Rational(1,3))]] +
                 [s.expand_complex(B).as_real_imag()[0]-s.Rational(1,6)])
    x = s.Matrix.hstack(s.Matrix([1,0,0]),s.Matrix([1,1,0])/s.sqrt(2),
                       s.Matrix([1,s.I,1])/s.sqrt(3))
    subs = {z[2*k]:s.re(x[k]) for k in range(9)}
    subs.update({z[2*k+1]:s.im(x[k]) for k in range(9)})
    J = s.simplify(f.jacobian(z).subs(subs))
    return x, f, J


def full_ray_selector():
    x, f, J = ray_data()
    gram = s.simplify(J*J.T)
    minors = [s.factor(gram[:k,:k].det()) for k in range(1,8)]
    assert all(v > 0 for v in minors)
    def realvec(a):
        return s.Matrix([v for w in list(a) for v in (s.re(w),s.im(w))])
    generators = [s.I*s.eye(3)[:,i]*s.eye(3)[i,:] for i in range(3)]
    for i in range(3):
        for j in range(i+1,3):
            a = s.zeros(3); a[i,j]=1; a[j,i]=-1; generators.append(a)
            a = s.zeros(3); a[i,j]=s.I; a[j,i]=s.I; generators.append(a)
    tangents = [realvec(a*x) for a in generators]
    tangents += [realvec(x*s.diag(*[s.I if k==j else 0 for k in range(3)])) for j in range(3)]
    T = s.Matrix.hstack(*tangents)
    assert T.rank()==11 and J*T==s.zeros(7,12)
    G=s.simplify(x.conjugate().T*x); B=s.simplify(G[0,1]*G[1,2]*G[2,0])
    assert B==(1+s.I)/6 and G.det()==s.Rational(1,6)
    A=x*s.diag(1,2,4)*x.conjugate().T
    C=x*s.diag(2,3,5)*x.conjugate().T
    assert s.simplify(s.trace((A*C-C*A)**3)/s.I)==1
    return dict(status='PASS', arithmetic='exact', real_fields=18,
                constraint_rank=7, normal_Hessian_rank=7, symmetry_kernel_dimension=11,
                symmetry='U(3) times three ray U(1)s, divided by their common central U(1)',
                potential='sum(norm_i^2-1)^2 + (abs(g12)^2-1/2)^2 + (abs(g13)^2-1/3)^2 + (abs(g23)^2-1/3)^2 + (Re(g12*g23*g31)-1/6)^2',
                Jacobian=[[str(v) for v in row] for row in J.tolist()],
                Jacobian_Gram_leading_minors=list(map(str,minors)),
                zero_orbits=2, Bargmann_branches=['(1+i)/6','(1-i)/6'],
                CP_cubic_branches=[1,-1],
                proof='H=2 J^T J at a zero. Positive Gram minors give seven positive normals. All eleven kernel directions are symmetry tangents. Nonzero overlaps allow two phases to be removed; the loop has fixed modulus and real part, leaving precisely its two conjugate signs. Positive definite Gram fixes the frame up to U(3).',
                boundary='Coefficients, overlap targets and Yukawa weights are supplied. The degree-12 potential is not derived from native E6 dynamics or a UV completion.',
                prior_owners=['w33_pass11374_11378_five_frontier.py','w33_pass11326_noncommuting_flavor_operator.py'])


def orthobasis(matrices, tol=1e-9):
    basis=[]
    for a in matrices:
        v=a.copy()
        for _ in range(2):
            for b in basis: v-=np.vdot(b,v).real*b
        norm=np.linalg.norm(v)
        if norm>tol: basis.append(v/norm)
    return basis


def native_channel_span():
    import w33_pass11287_quotient_self_energy_matrix as q
    d=json.loads((ROOT/'data/w33_pass11287_quotient_self_energy_matrix.json').read_text())
    H=np.array(d['Hessian']); C=np.array(d['cubic_tensor'])
    M=np.array(d['Weyl_mass_real'])+1j*np.array(d['Weyl_mass_imag'])
    Y=np.array(d['Yukawa_real'])+1j*np.array(d['Yukawa_imag'])
    _,_,O,groups=q.cut_groups(H,C,M,Y)
    mats=[g[k][8:,8:] for g in groups for k in ['R','I'] if np.linalg.norm(g[k][8:,8:])>1e-13]
    normed=[a/np.linalg.norm(a) for a in mats]
    sv=np.linalg.svd(np.array([a.ravel() for a in normed]),compute_uv=False)
    basis=orthobasis(normed)
    assert len(basis)==7 and sv[6]>0.1 and sv[7]<1e-10
    algebra=orthobasis([np.eye(14)]+basis)
    for _ in range(5):
        old=len(algebra)
        algebra=orthobasis(algebra+[a@b for a in algebra for b in algebra])
        if len(algebra)==old: break
    assert len(algebra)==8
    comm=np.array([np.concatenate([(a@b-b@a).ravel() for b in basis]) for a in algebra]).T
    _,cs,vh=np.linalg.svd(comm,full_matrices=False)
    center_coeff=vh[cs<1e-9]
    assert len(center_coeff)==5
    centers=[sum(c*a for c,a in zip(v,algebra)) for v in center_coeff]
    rng=np.random.default_rng(11380)
    central=sum(c*a for c,a in zip(rng.normal(size=5),centers))
    w,U=np.linalg.eigh((central+central.T)/2)
    clusters=[]
    for i,v in enumerate(w):
        if not clusters or abs(v-w[clusters[-1][0]])>1e-7: clusters.append([i])
        else: clusters[-1].append(i)
    assert sorted(map(len,clusters))==[1,1,2,4,6]
    projectors=[U[:,ix]@U[:,ix].T for ix in clusters]
    residual=max(np.linalg.norm(p@a-a@p) for p in projectors for a in basis)
    assert residual<1e-9
    block_algebra_dimensions=[len(orthobasis([p@a@p for a in algebra])) for p in projectors]
    assert sorted(block_algebra_dimensions)==[1,1,1,1,4]
    # The noncommuting M2(R) block acts with multiplicity two on FOUR normal
    # coordinates. Construct its actual intertwiner, not just dimensions.
    active=block_algebra_dimensions.index(4)
    assert len(clusters[active])==4
    W=U[:,clusters[active]]
    A=sum(t*a for t,a in zip(rng.normal(size=7),basis))
    ew,EV=np.linalg.eigh(W.T@A@W)
    assert abs(ew[0]-ew[1])<1e-9 and abs(ew[2]-ew[3])<1e-9 and ew[2]-ew[1]>1e-3
    left=W@EV[:,:2]; right=W@EV[:,2:]
    chosen=max(basis,key=lambda a:np.linalg.norm(left.T@a@right))
    lu,lv,lr=np.linalg.svd(left.T@chosen@right)
    assert abs(lv[0]-lv[1])<1e-9
    intertwiner=np.column_stack([left@lu,right@lr.T])
    compressed=[]; compression_error=0.
    for a in basis:
        full=intertwiner.T@a@intertwiner
        little=np.array([[np.trace(full[2*i:2*i+2,2*j:2*j+2])/2 for j in range(2)] for i in range(2)])
        compressed.append(little.tolist())
        compression_error=max(compression_error,float(np.linalg.norm(full-np.kron(little,np.eye(2)))))
    assert compression_error<1e-9
    tree_normal=O[:,8:].T@H@O[:,8:]
    identity_error=np.linalg.norm(np.eye(14)-sum(np.vdot(a,np.eye(14)).real*a for a in basis))
    tree_error=np.linalg.norm(tree_normal-sum(np.vdot(a,tree_normal).real*a for a in basis))
    assert identity_error<1e-9 and tree_error<1e-9
    # Rotate all 22 external scalar coordinates, including their internal legs.
    R=np.linalg.qr(rng.normal(size=(22,22)))[0]
    Hr=R.T@H@R
    Cr=np.einsum('ia,jb,kc,ijk->abc',R,R,R,C,optimize=True)
    Yr=np.einsum('ia,ijk->ajk',R,Y)
    _,_,Or,gr=q.cut_groups(Hr,Cr,M,Yr)
    nr=Or[:,8:]; n=O[:,8:]; change=n.T@R@nr
    assert np.linalg.norm(change.T@change-np.eye(14))<1e-10
    mr=[change@g[k][8:,8:]@change.T for g in gr for k in ['R','I'] if np.linalg.norm(g[k][8:,8:])>1e-13]
    br=orthobasis([a/np.linalg.norm(a) for a in mr])
    err=max(np.linalg.norm(a-sum(np.vdot(b,a).real*b for b in basis)) for a in br)
    assert len(br)==7 and err<1e-9
    spectra={str(t):np.linalg.eigvalsh(sum(q.rho(g,t) for g in groups)[8:,8:]).tolist() for t in [.002,.02,.2]}
    assert all(min(v)>1e-9 for v in spectra.values())
    # Concrete absorptive-truncation inversion; no unknown real hard matching
    # is silently set equal to a measured propagator.
    absorptive=sum(q.rho(g,.02) for g in groups)[8:,8:]
    denominator=(.02+.003j)*np.eye(14)-.0001*tree_normal+1j*absorptive
    little_full=intertwiner.T@denominator@intertwiner
    little=np.array([[np.trace(little_full[2*i:2*i+2,2*j:2*j+2])/2 for j in range(2)] for i in range(2)])
    inv_small=intertwiner@np.kron(np.linalg.inv(little),np.eye(2))@intertwiner.T
    for i,p in enumerate(projectors):
        if i!=active: inv_small+=p/(np.trace(p@denominator)/np.trace(p))
    inverse_error=np.linalg.norm(inv_small@denominator-np.eye(14))
    assert inverse_error<1e-9
    return dict(status='PASS', arithmetic='floating-point native certificate, not an exact representation theorem',
                mass_pair_groups=len(groups), nonzero_group_matrices=len(mats), symmetric_span_dimension=7,
                generated_algebra_dimension=8, center_dimension=5, central_block_dimensions=list(map(len,clusters)),
                central_block_algebra_dimensions=block_algebra_dimensions,
                noncommuting_block_normal_dimension=4, noncommuting_block_irrep_dimension=2,
                noncommuting_block_multiplicity=2,
                noncommuting_block_intertwiner=intertwiner.tolist(),
                compressed_two_by_two_channel_matrices=compressed,
                compression_residual=compression_error,
                tree_Hessian_span_residual=float(tree_error), identity_span_residual=float(identity_error),
                absorptive_truncation_inverse_residual=float(inverse_error),
                inverse_control='D=(0.02+0.003i)I - 0.0001 H_normal + i rho(0.02), inverse by four scalar reciprocals and one 2x2 inverse tensor I2; unknown hard real matching is not included',
                normal_basis_in_original_coordinates=O[:,8:].tolist(),
                orthonormal_channel_basis=[a.tolist() for a in basis],
                central_projectors=[p.tolist() for p in projectors],
                normalized_span_singular_values=sv.tolist(), basis_rotation_residual=float(err),
                central_projector_commutator_residual=float(residual), absorptive_spectra=spectra,
                boundary='Grouped equal-mass cuts are invariant under internal degenerate-basis rotations. This compresses the native nonanalytic cut sector; unrestricted local hard counterterms still have 105 entries. No full renormalized pole spectrum follows.',
                prior_owners=['w33_pass11287_quotient_self_energy_matrix.py','w33_pass11337_quotient_hard_ward_columns.py'])


def modular_rank(a, p=1000003):
    a=np.array(a,dtype=np.int64)%p; row=0; pivots=[]
    for col in range(a.shape[1]):
        nz=np.flatnonzero(a[row:,col])
        if not len(nz): continue
        k=row+int(nz[0]); a[[row,k]]=a[[k,row]]
        a[row]=(a[row]*pow(int(a[row,col]),-1,p))%p
        for k in range(row+1,len(a)):
            if a[k,col]: a[k]=(a[k]-a[k,col]*a[row])%p
        pivots.append(col); row+=1
        if row==len(a): break
    return row,pivots


def native_rotation_data():
    from w33_pass11289_cycle_gram_gluing_flux import cycle_basis
    inc,cyc=cycle_basis(); edges=[tuple(np.flatnonzero(inc[:,k])) for k in range(160)]
    adj={}
    for k in np.flatnonzero(cyc[:,0]):
        i,j=edges[k]; adj.setdefault(i,[]).append(j); adj.setdefault(j,[]).append(i)
    walk=[min(adj)]
    while len(walk)<len(adj): walk.append(next(j for j in adj[walk[-1]] if j not in walk))
    anis=[(1,0),(1,1),(0,1),(-1,1),(-1,0),(-1,-1),(0,-1),(1,-1)]
    frames=[s.eye(2) for _ in range(80)]
    for k,(x,y) in zip(walk,anis): frames[k]=s.Matrix([[1+s.Rational(x,10),s.Rational(y,10)],[s.Rational(y,10),1-s.Rational(x,10)]])
    rot=s.Matrix([[0,-1],[1,0]])
    aa=[]; bb=[]
    for i,j in edges:
        a=frames[j]*frames[i].adjugate()
        aa.append(s.trace(a)); bb.append(s.trace(rot*a))
    return inc, edges, walk, frames, aa, bb


def native_rotational_lapse():
    from scipy.optimize import root
    inc,edges,walk,frames,aa,bb=native_rotation_data()
    B=s.Matrix(inc); D=B[1:,:]; U=B.applyfunc(abs)
    assert B*s.Matrix(bb)==s.zeros(80,1) and all(a>0 for a in aa)
    G=U*s.diag(*bb)*D.T/2
    assert G.rank()==6
    # Exact invertibility of the shift/boost Jacobian with vertex 0 fixed.
    S=np.zeros((160,160),dtype=np.int64)
    for i,j in edges:
        ei=10*frames[i]; ej=10*frames[j]
        # Coordinate shifts enter the time column as e_i*s_i. The common
        # edge weight is det(e_i)*e_j + det(e_j)*e_i, not a frame-coordinate
        # weight missing its final e_i factor. Scale 2000 clears 1/2 and
        # the cubic denominators of the coframes.
        di=ei.det()*ej+ej.det()*ei
        dj=di
        for a in range(2):
            for b in range(2):
                S[2*i+a,2*i+b]+=int(di[b,a]); S[2*i+a,2*j+b]-=int(di[b,a])
                S[2*j+a,2*i+b]-=int(dj[b,a]); S[2*j+a,2*j+b]+=int(dj[b,a])
    r,piv=modular_rank(S[2:,2:]); assert r==158
    dn=np.array(D,float); un=np.array(U,float); a=np.array(aa,float); b=np.array(bb,float)
    K=-(dn*a)@dn.T; gn=np.array(G,float)
    reduced=-gn@np.linalg.solve(K,gn.T)
    eig=np.linalg.eigvalsh(reduced)
    assert np.count_nonzero(eig>1e-10)==6 and min(eig)>-1e-12
    direction=np.zeros(80); direction[walk[1]]=1
    def energy(N):
        weights=un.T@N
        def station(theta):
            delta=dn.T@theta
            return dn@(weights*(-a*np.sin(delta)+b*np.cos(delta))/2)
        sol=root(station,np.zeros(79),jac=lambda th:(dn*(weights*(-a*np.cos(dn.T@th)-b*np.sin(dn.T@th))/2))@dn.T,tol=1e-11)
        assert np.linalg.norm(station(sol.x))<1e-10
        delta=dn.T@sol.x
        return float(weights@(a*np.cos(delta)+b*np.sin(delta))/2)
    N=np.ones(80); base=energy(N)
    predicted=float(direction@reduced@direction)
    fd=[(energy(N+h*direction)+energy(N-h*direction)-2*base)/h**2 for h in [.01,.005]]
    assert max(abs(v-predicted) for v in fd)<1e-7
    return dict(status='PASS', arithmetic='exact stationarity and rank, modular full-rank auxiliary Jacobian, numerical independent curvature control',
                vertices=80,edges=160, cycle_vertices=list(map(int,walk)),
                coframes=[[[str(v) for v in row] for row in e.tolist()] for e in frames],
                edge_a=list(map(str,aa)),edge_b=list(map(str,bb)),
                pair_action='(Ni+Nj)/2 * [a_ij cos(theta_j-theta_i)+b_ij sin(theta_j-theta_i)+det(ei)+det(ej)]',
                stationary_cycle_force='B b = 0 exactly', mixed_lapse_rotation_rank=6,
                reduced_lapse_Hessian_rank=6, reduced_lapse_nonzero_eigenvalues=eig[-6:].tolist(),
                shift_boost_Jacobian_integer_scale=2000, shift_boost_Jacobian_mod_prime=1000003,
                shift_boost_Jacobian_rank=r, shift_boost_Jacobian_pivots=piv,
                finite_difference_curvatures=fd, predicted_curvature=predicted,
                proof='K=-D diag(a) D^T is negative definite on the fixed common-rotation gauge. Eliminating rotations gives H_N=-G K^-1 G^T, positive semidefinite and of exact rank rank(G)=6. The shift/boost Jacobian is a graph Laplacian with positive definite weights (det(e_i)e_j+det(e_j)e_i)/2, hence invertible after common modes are fixed; the modular determinant replays this independently.',
                boundary='This disproves general lapse linearity for the actual native pair action beyond the collinear boost slice. It is not a complete Hamiltonian Dirac count or a count of six propagating ghosts.',
                prior_owners=['w33_pass11338_native_boost_constraint_slice.py','w33_pass11376_native_frame_flatness in w33_pass11374_11378_five_frontier.py'],
                literature='https://arxiv.org/abs/1505.01450')


def wall_factorization():
    b,q=s.symbols('b q',positive=True)
    bw=s.Rational(5,4); h=4*(q*q-q)/5; c=q/5
    rho=q*(4-3*q)/10; C=rho/3+q*q/2+h
    F=s.factor(C/b-rho*b*b/3-q*q/(2*b*b)-h)
    polynomial=s.factor(30*b*b*F/(q*(b-1)))
    assert s.factor(polynomial.subs(q,1))==15-b-b*b-b**3
    assert s.simplify(polynomial.subs(q,s.Rational(4,3))-(20-8*b))==0
    v0=-F/(q*c*b*b); chi0=1/b-1/bw
    ode=s.factor(-s.diff(b**4*s.diff(v0,b),b)+2*(h*b*b-q*q)*v0/F)
    assert ode==0 and s.simplify(s.diff(chi0,b)-q*c*v0/F)==0
    assert v0.subs(b,1)==0 and s.factor(s.diff(v0,b).subs(b,bw))==0
    assert chi0.subs(b,bw)==0
    # Local ground-state identity, independent of trial-function choice.
    v=s.Function('v')(b)
    lhs=b**4*s.diff(v,b)**2+2*(h*b*b-q*q)/F*v*v
    rhs=b**4*v0*v0*s.diff(v/v0,b)**2+s.diff(b**4*s.diff(v0,b)/v0*v*v,b)
    assert s.factor(lhs-rhs)==0
    da,V,dl=s.symbols('da V dl', real=True)
    assert s.expand(((da-dl)-(V-dl))**2-(da-V)**2)==0
    return dict(status='PASS', arithmetic='exact symbolic identities for 1<q<=4/3 (positive winding coefficient)',
                cap_F=str(F), ground_state_v=str(s.factor(v0)), ground_state_chi=str(chi0),
                cap_positivity='F=q(b-1)H/(30b^2); H is affine in q, and its endpoints 15-b-b^2-b^3 and 20-8b are both positive for 1<=b<=5/4',
                full_radial_form='integral [F chi_prime^2 + c^2 b^4 v_prime^2/2 + 2qc chi v_prime + h b^2 c^2 v^2/F] db',
                factorization='integral F(chi_prime-qc*v/F)^2 + c^2/2 integral b^4*v0^2*((v/v0)_prime)^2; boundary terms vanish for the stated regular parity domains',
                endpoint_conditions='pole v=0 with v=O(b-1); odd: chi(bw)=0,v_prime(bw)=0; even: chi_prime(bw)=0,v(bw)=0',
                odd_kernel_per_torus_axis=1, even_Wilson_kernel_per_torus_axis=1,
                radial_vector_kernel_total=4,
                axion_gauge='y_new=y+lambda, a_new=a-lambda, V_new=V-dlambda: da-V is invariant',
                physical_constant_axion_moduli=0,
                correction='The two constant axion shifts of Pass11377 are exact configuration symmetries but are torus-coordinate gauge transformations in the dynamical winding geometry with no fixed boundary frame.',
                boundary='Full radial Maxwell/metric vector block, not all torus/azimuthal harmonics, scalar breathing, or wall bending. Nonlinear lifting of its four vector zeros remains open. Fixed background axion positivity stays valid.',
                prior_owners=['w33_pass11339_coupled_maxwell_metric_modes.py','w33_pass11346_winding_wall_mixed_vector_modes.py','w33_pass11374_11378_five_frontier.py'])


def fixed_action_scale():
    b,L=s.symbols('b L',positive=True)
    q=s.Rational(4,3); h=s.Rational(16,45); c=s.Rational(4,15); bw=s.Rational(5,4)
    F=-8*(b-1)*(2*b-5)/(45*b*b); fw=F.subs(b,bw)
    T=4*s.sqrt(fw)/bw; area=s.sqrt(fw)*bw*bw/c
    Kwall=s.sqrt(fw)*(s.diff(F,b).subs(b,bw)/(2*fw)+2/bw)
    volume=2*b*b/c  # two caps, with angular factor (2pi)^3 divided out
    EH=s.integrate(-h/b**2*volume,(b,1,bw))
    axion=s.integrate(h/b**2*volume,(b,1,bw))
    maxwell=s.integrate(q*q/(2*b**4)*volume,(b,1,bw))
    ghy=-2*area*Kwall; wall=T*area
    action=s.factor(maxwell+(EH+axion+ghy)*L**2+wall*L**3)
    assert (EH+axion,maxwell,ghy,wall)==(0,s.Rational(4,3),-2,s.Rational(4,3))
    assert s.diff(action,L).subs(L,1)==0 and s.diff(action,L,2).subs(L,1)==4
    return dict(status='PASS', arithmetic='exact', scale_path='g(L)=L^2 g0 with coordinate flux and winding, K,h,T,rho=0 held fixed',
                convention='K=M_Pl^2 is the Einstein-Hilbert coefficient; repository kappa^2 uses this convention',
                Euclidean_action='-K/2 integral R - K sum_caps integral Kextr + 1/4 integral F^2 + h/2 sum_axes integral (dtheta)^2 + rho integral 1 + T integral_wall 1',
                angular_factor='(2*pi)^3', coefficients_without_angular_factor=dict(EH=str(EH),axion=str(axion),Maxwell=str(maxwell),GHY=str(ghy),wall=str(wall)),
                action_without_angular_factor=str(action), stationary_L='1', second_derivative_at_L1='4',
                physical_tension='T_phys=(64/75)*K_phys/L', physical_ratio='K_phys^2*R_phys/T_phys^2=5/8',
                boundary='Positive global uniform-scale curvature with fixed supplied couplings is not positivity of nonuniform conformal or wall modes. K_phys and T_phys remain inputs; their ratio can fix the classical scale but W33 does not yet predict them.',
                prior_owners=['w33_pass11341_winding_lifted_wall.py','w33_pass11374_11378_five_frontier.py'])


def payload():
    parts={}
    for name,f in [('11379_full_ray_selector',full_ray_selector),('11380_native_channel_span',native_channel_span),
                   ('11381_native_rotational_lapse',native_rotational_lapse),('11382_wall_factorization',wall_factorization),
                   ('11383_fixed_action_scale',fixed_action_scale)]:
        parts[name]=f(); print(name,parts[name]['status'],flush=True)
    return dict(status='PASS',sections=parts)


if __name__=='__main__':
    out=payload()
    (ROOT/'data/w33_pass11379_11383_full_frontier.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
