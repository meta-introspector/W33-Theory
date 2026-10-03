"""Five native-dynamics investigations; supplied couplings remain explicit.

No observed mass/CP/CC parameter is used. Exact identities are separated from
numerical native-field stability and imported tree-multigravity results.
"""
from pathlib import Path
import sys
import json
import itertools
import numpy as np
import sympy as s
from scipy.optimize import root

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analysis'))


def native_tensors():
    data=json.loads((ROOT/'data/w33_pass11384_native_inputs.json').read_text())
    B=np.zeros((78,27,27),int); d=np.zeros((27,27,27),int); eps=np.zeros((3,3,3),int)
    tri=data['signed_triads']
    for a,rows in enumerate(data['E6_sparse_basis']):
        for i,j,v in rows:B[a,i,j]=v
    for i,j,k,v in tri:
        for ix in itertools.permutations((i,j,k)):d[ix]=v
    for p in itertools.permutations(range(3)):
        eps[p]=(-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
    for b in B:
        for bb in [b,b.T]:
            assert not np.any(np.einsum('ai,ajk->ijk',bb,d)+np.einsum('aj,iak->ijk',bb,d)+np.einsum('ak,ija->ijk',bb,d))
    return B,(tri,d,eps)


def compact_generators(B):
    unique={}
    for b in B:
        for kind,a in [(0,b+b.T),(1,b-b.T)]:
            nz=np.flatnonzero(a)
            if not len(nz):continue
            a=a//int(np.gcd.reduce(abs(a.ravel()[nz])))
            if a.ravel()[nz[0]]<0:a=-a
            unique[kind,tuple(a.ravel())]=a*(1j if kind else 1)
    assert len(unique)==78
    he=np.array(list(unique.values()),complex)
    hf=[]
    for i in range(3):
        for j in range(i+1,3):
            a=np.zeros((3,3),complex);a[i,j]=a[j,i]=1;hf.append(a)
            a=np.zeros((3,3),complex);a[i,j]=-1j;a[j,i]=1j;hf.append(a)
    hf.extend([np.diag([1,-1,0]),np.diag([1,1,-2])])
    def orth(a):
        a=np.array(a);G=np.einsum('aij,bji->ab',a,a).real;w,U=np.linalg.eigh(G)
        assert min(w)>0
        return np.einsum('ab,bij->aij',(U/np.sqrt(w)).T,a)
    return np.array([np.kron(h,np.eye(3)) for h in orth(he)]+[np.kron(np.eye(27),h.T) for h in orth(hf)])


def native_cp_vacuum():
    import w33_pass11255_11259_g26_common as G
    import w33_20261001_degree18_phase_completion as D
    B,tensors=native_tensors(); H=compact_generators(B)
    E=D.I.C.plane()/np.sqrt(3)
    # The prior plane is D-flat for every complex q, not only the old T ray.
    assert max(np.max(abs(E.conj().T@h@E)) for h in H)<1e-13
    fs=G.invariants();calc=s.lambdify(G.VARS,fs,'numpy');jc=s.lambdify(G.VARS,s.Matrix(fs).jacobian(G.VARS),'numpy')
    a=(-1+1j*np.sqrt(7))/4;target=np.array([a,a*a,a**3])
    unpack=lambda x:x[:3]+1j*x[3:]
    def fun(x):
        z=np.array(calc(*unpack(x)))-target
        return np.r_[z.real,z.imag]
    def jac(x):
        J=np.array(jc(*unpack(x)));return np.block([[J.real,-J.imag],[J.imag,J.real]])
    rng=np.random.default_rng(11384);sol=None
    for _ in range(40):
        p=root(fun,rng.normal(size=6)*.6,jac=jac,tol=1e-11)
        if np.linalg.norm(fun(p.x))<1e-9:sol=p;break
    assert sol is not None
    # Newton refinement eliminates discovery tolerance from the witness.
    x=sol.x.copy()
    for _ in range(3):x-=np.linalg.solve(jac(x),fun(x))
    phi=E@unpack(x)
    original=D.I.tensors
    D.I.tensors=lambda:tensors
    try:
        i6,i12,g6,g12=D.I.evaluate(phi);i18,g18=D.evaluate(phi)
    finally:D.I.tensors=original
    assert max(abs(np.array([i6,i12,i18])-target))<1e-9
    acts=np.einsum('aij,j->ai',H,phi)
    mu=np.real(acts.conj()@phi);assert max(abs(mu))<1e-12
    Jmu=np.sqrt(2)*np.concatenate([acts.real,acts.imag],axis=1)
    def realJ(g):return np.array([np.r_[g.real,-g.imag],np.r_[g.imag,g.real]])/np.sqrt(2)
    J6=realJ(g6);J12=realJ(g12-2*i6*g6);J18=realJ(g18-3*i6*i6*g6)
    HF=8*np.outer([i6.real,i6.imag],[i6.real,i6.imag])+2*np.diag([1,0])
    HV=Jmu.T@Jmu+J6.T@HF@J6+2*J12.T@J12+2*J18.T@J18
    ev=np.linalg.eigvalsh(HV)
    orbit=np.sqrt(2)*np.concatenate([-acts.imag,acts.real],axis=1).T
    assert np.linalg.matrix_rank(orbit,tol=1e-9)==78
    assert sum(ev>1e-8)==84 and max(abs(ev[:78]))<1e-8
    ward=np.linalg.norm(HV@orbit)/(1+np.linalg.norm(HV)*np.linalg.norm(orbit))
    assert ward<1e-12
    # Global invariant minimum is derived by completion of squares in I6.
    xx,yy=s.symbols('x y',real=True)
    landau=(xx*xx+yy*yy)**2-(xx*xx+yy*yy)+xx*xx+xx/2
    squares=(xx*xx+yy*yy-s.Rational(1,2))**2+(xx+s.Rational(1,4))**2-s.Rational(5,16)
    assert s.expand(landau-squares)==0
    return dict(status='PASS',arithmetic='exact invariant minimization; numerical full162-field native witness',
                potential='1/2 sum_a mu_a^2 + |I6|^4-|I6|^2+(Re I6)^2+Re I6/2 + |I12-I6^2|^2 + |I18-I6^3|^2',
                invariant_minimum='-5/16',I6_branches=['(-1+i*sqrt(7))/4','(-1-i*sqrt(7))/4'],
                invariant_relations=['I12=I6^2','I18=I6^3'],
                CP_order='Im I6=+/-sqrt(7)/4, invariant under native compact E6 x SU3 and the remaining common-phase Z6',
                cartan_coordinates=x.tolist(),full_field_real=phi.real.tolist(),full_field_imag=phi.imag.tolist(),
                invariant_residual=float(max(abs(np.array([i6,i12,i18])-target))),moment_residual=float(max(abs(mu))),
                gauge_orbit_rank=78,positive_normal_modes=84,normal_eigenvalues=ev[78:].tolist(),
                gauge_Hessian_residual=float(ward),tensor_generator_checks=156,
                boundary='A native invariant CP-breaking completion with supplied real EFT couplings, not a derivation of those coefficients, observed Yukawas or a UV-complete model. Full-field Hessian is numerical, not interval-certified. Invariant CP branches are gauge-distinct; no claim excludes every additional generalized-CP symmetry outside the named group.',
                prior_owners=['w33_20261001_global_e6_cartan_covariants.py','w33_20261001_degree18_phase_completion.py','w33_pass11379_11383_full_frontier.py'])


def hard_matching_classes():
    import w33_pass11379_11383_full_frontier as P
    d=P.native_channel_span();dims=d['central_block_dimensions']
    basis=list(map(np.array,d['orthonormal_channel_basis']));Ps=list(map(np.array,d['central_projectors']))
    W=np.array(d['noncommuting_block_intertwiner'])
    assert sum(n*(n+1)//2 for n in dims)==36
    # Exact model of the double-commutant architecture:
    # scalar blocks R6,R2,R1,R1 and M2 acting on R2 tensor R2.
    assert 21+3+1+1+3==29 and 36+4+1+1+4==46
    rng=np.random.default_rng(11385)
    Q=np.eye(14)
    for p,n in zip(Ps,dims):
        if n==4:continue
        w,U=np.linalg.eigh(p);V=U[:,w>.5]
        R=np.linalg.qr(rng.normal(size=(n,n)))[0]
        Q+=V@(R-np.eye(n))@V.T
    R=np.linalg.qr(rng.normal(size=(2,2)))[0]
    Q+=W@(np.kron(np.eye(2),R)-np.eye(4))@W.T
    assert np.linalg.norm(Q.T@Q-np.eye(14))<1e-10
    assert max(np.linalg.norm(Q@a-a@Q) for a in basis)<1e-9
    # A normal counterterm can satisfy stationary Goldstone Ward columns yet
    # fail the central projectors; Ward completion alone does not protect cuts.
    u=np.linalg.eigh(Ps[0])[1][:,-1];v=np.linalg.eigh(Ps[1])[1][:,-1]
    counter=np.outer(u+v,u+v)
    obstruction=max(np.linalg.norm(p@counter-counter@p) for p in Ps)
    assert obstruction>.1
    return dict(status='PASS',arithmetic='exact dimensions of named real block models, numerical native embedding controls',
                unrestricted_symmetric_normal_dimension=105,central_projector_preserving_dimension=36,
                matrices_commuting_with_cut_algebra_symmetric_dimension=29,full_real_commutant_dimension=46,
                invariant_under_full_orthogonal_commutant_dimension=7,
                protecting_group='O(6) x O(2) x O(1) x O(1) x O(2); last factor acts on multiplicity of M2',
                group_control_commutator=float(max(np.linalg.norm(Q@a-a@Q) for a in basis)),
                Ward_compatible_normal_counterterm=counter.tolist(),central_breaking_residual=float(obstruction),
                boundary='This identifies sufficient extra UV symmetry and distinct weaker matching classes. The protecting group is a symmetry of the native cut matrices; its realization by the full interacting UV theory is not proved. Zero stationary Goldstone columns alone leave the105 normal entries free.',
                prior_owners=['w33_pass11337_quotient_hard_ward_columns.py','w33_pass11379_11383_full_frontier.py'])


def hub_gravity_architecture():
    from w33_pass11289_cycle_gram_gluing_flux import cycle_basis
    B,_=cycle_basis();L=B@B.T
    star=np.zeros((81,81),int)
    for j in range(1,81):
        star[0,0]+=1;star[j,j]+=1;star[0,j]=star[j,0]=-1
    assert np.linalg.matrix_rank(star)==80
    ev=np.linalg.eigvalsh(star)
    assert np.max(abs(ev-np.array([0]+[1]*79+[81])))<1e-10
    # A native permutation extends by fixing the hub. Exact invariance needs
    # only equal spoke coefficients; the native matter edges remain untouched.
    perm=np.r_[0,np.arange(80,0,-1)]
    assert np.array_equal(star[np.ix_(perm,perm)],star)
    vals=np.linalg.eigvalsh(L)
    expected=np.sort([0,8]+[4-np.sqrt(6)]*24+[4+np.sqrt(6)]*24+[4]*30)
    assert np.max(abs(vals-expected))<1e-10
    return dict(status='PASS',arithmetic='exact tree construction and matrix identities; constraint result imported from primary multimetric literature',
                spin2_fields=81,spin2_interactions=80,spin2_graph='one invariant hub g0 with80 equal Hassan-Rosen bimetric spokes; no cyclic spin2 interaction',
                action='sum_I K_I sqrt(-gI) R(gI)/2 - 2m^4 sum_i sqrt(-g0) sum_n beta_n e_n(sqrt(g0^-1 gi)) + S_matter(g0,chi)',
                matter_action='-1/2 integral sqrt(-g0)[sum_i (dchi_i)^2 + g^2 sigma^2 chi^T L_native chi]; all matter sees g0 only',
                star_species_Laplacian_spectrum={'0':1,'1':79,'81':1},
                native_internal_matter_spectrum={'0':1,'8':1,'4-sqrt(6)':24,'4+sqrt(6)':24,'4':30},
                conditional_spin2_polarizations=2+5*80,
                scope='Regular ghost-free tree-bimetric branch with individual EH terms, admissible square-root branch and hub-only minimal matter coupling. Select bimetric parameters with positive Fierz-Pauli coefficient for the displayed linear mass graph.',
                boundary='This is a changed architecture: native cycles are internal matter interactions, not80 cyclically coupled spacetime metrics. The hub, continuum, EH scales and couplings are supplied. The402 polarization count uses the published tree-multimetric theorem, not a new independent Dirac computation or a physical spectrum prediction.',
                prior_owners=['w33_pass11381_native_rotational_lapse in w33_pass11379_11383_full_frontier.py','w33_pass11296_covariant_matter_constraints.py','w33_pass11328_determinant_multivielbein_scope.py'],
                primary_source='https://arxiv.org/html/2604.07625v1')


def relaxed_wall_action(u,v):
    h=s.Rational(16,45);T=s.Rational(64,75)
    q=2*u*u*(1-u/v)+s.sqrt(4*u**4*(1-u/v)**2+2*h*u*u)
    c=q*q/(4*u**3)-h/(2*u);C=q*q/(2*u)+h*u
    F=C/v-q*q/(2*v*v)-h;Fp=-C/v**2+q*q/v**3
    action=q-(v*v*Fp+4*v*F)/c+T*s.sqrt(F)*v*v/c
    return action,q,c,C,F


def scalar_wall_relaxation():
    u,v=s.symbols('u v',positive=True)
    action,q,c,C,F=relaxed_wall_action(u,v);sub={u:1,v:s.Rational(5,4)}
    assert s.simplify(q/c*(1/u-1/v)-1)==0
    assert s.simplify(action.subs(sub))==s.Rational(2,3)
    grad=s.Matrix([s.simplify(s.diff(action,x).subs(sub)) for x in [u,v]])
    assert grad==s.zeros(2,1)
    H=s.Matrix(2,2,lambda i,j:s.simplify(s.diff(action,[u,v][i],[u,v][j]).subs(sub)))
    assert H==s.Matrix([[s.Rational(400,49),-s.Rational(1536,245)],[-s.Rational(1536,245),s.Rational(16064,3675)]])
    assert H.det()==-s.Rational(13312,3675)
    direction=s.Matrix([1,s.Rational(5,4)]);curvature=s.factor((direction.T*H*direction)[0]);assert curvature==-s.Rational(100,147)
    f=s.lambdify((u,v),action,'numpy');h=.0005
    fd=(f(1+h,1.25*(1+h))+f(1-h,1.25*(1-h))-2*f(1,1.25))/h**2
    assert abs(fd-float(curvature))<1e-5
    # Every nearby cap satisfies bulk field equations and pole smoothness.
    b=s.symbols('b',positive=True);FF=C/b-q*q/(2*b*b)-s.Rational(16,45)
    assert s.simplify(FF.subs(b,u))==0
    assert s.simplify(s.diff(FF,b).subs(b,u)-2*c)==0
    assert s.simplify(-s.diff(FF,b,2)/2-s.diff(FF,b)/b-q*q/(2*b**4))==0
    return dict(status='PASS',arithmetic='exact fixed-flux bulk-relaxed two-cap second variation and numerical difference control',
                fixed_inputs='K=1,h=16/45,T=64/75,rho=0, torus periods2pi and magnetic flux q/c*(1/u-1/v)=1',
                collective_coordinates='u=regular pole torus radius, v=wall torus radius; 0<u<v',
                q=str(q),c=str(c),C=str(C),wall_F=str(F),reduced_action_without_angular_factor=str(action),
                Hessian=[[str(x) for x in row] for row in H.tolist()],Hessian_determinant=str(H.det()),
                negative_test_direction=['1','5/4'],negative_direction_curvature=str(curvature),finite_difference_curvature=float(fd),
                bulk_relaxed_Hessian_eigenvalues=np.linalg.eigvalsh(np.array(H,float)).tolist(),
                scalar_reduction='S/(2pi)^3 = -2K b_pole^2 + 2 integral [-K(2bb_prime r_prime+r b_prime^2)/n + b^2 A_prime^2/(2nr) + nr(rho b^2+h)] dxi + T r_wall b_wall^2, including reduced GHY and smooth-pole term',
                conclusion='The Euclidean wall saddle is not a local minimum on this real fixed-flux regular-cap family, despite positive radial-vector and uniform-scale blocks.',
                boundary='This is a genuine bulk-relaxed Euclidean negative second variation, not a frozen-bulk junction artifact. It does not by itself fix the full Morse index, conformal integration contour, fluctuation determinant, nucleation rate or Lorentzian instability; nonconstant wall bending and other harmonics remain open.',
                prior_owners=['w33_pass11324_quantized_warped_history_wall.py','w33_pass11341_winding_lifted_wall.py','w33_pass11379_11383_full_frontier.py'])


def native_dimensional_transmutation():
    from w33_pass11289_cycle_gram_gluing_flux import cycle_basis
    B,_=cycle_basis();L=B@B.T
    assert int(np.trace(L))==320 and int(np.trace(L@L))==1600
    sigma,mu,g,A,xi=s.symbols('sigma mu g A xi',positive=True)
    beta=25*g**4/s.pi**2
    V=sigma**4*(A+beta*s.log(sigma*sigma/(mu*mu)))
    logstar=-s.Rational(1,2)-A/beta
    stationary=s.factor(s.diff(V,sigma).subs(s.log(sigma*sigma/(mu*mu)),logstar))
    assert stationary==0
    energy=s.factor(V.subs(s.log(sigma*sigma/(mu*mu)),logstar))
    curvature=s.factor(s.diff(V,sigma,2).subs(s.log(sigma*sigma/(mu*mu)),logstar))
    assert energy==-beta*sigma**4/2 and curvature==8*beta*sigma**2
    VE=V/(xi*xi*sigma**4)
    slope=s.factor(s.diff(VE,sigma));assert slope==2*beta/(xi*xi*sigma)
    return dict(status='PASS',arithmetic='exact native spectral traces and one-loop symbolic identities',
                model='80 real native-graph scalars chi with M_chi^2=g^2 sigma^2 L_native and a canonical scalar sigma; flat-space one-loop Landau/MS potential, ignoring further interactions and gravity loops',
                native_vertices=80,native_degree=4,trace_L=320,trace_L_squared=1600,
                Coleman_Weinberg_log_coefficient='B=25*g^4/pi^2',
                scale='v=mu*exp(-1/4-A/(2B)); A is a renormalized quartic boundary condition, not a W33 prediction',
                flat_space_radial_mass_squared='8Bv^2',flat_space_vacuum_energy='-Bv^4/2',
                mass_energy_relation='m_sigma^2=-16 V(v)/v^2',
                induced_gravity='If the only Einstein coefficient is F(sigma)=xi sigma^2 with constant xi>0, VE=V/F^2 and dVE/dsigma=2B/(xi^2 sigma)>0.',
                induced_gravity_conclusion='The flat-space Coleman-Weinberg minimum is not a constant-field stationary gravitational vacuum in this one-loop constant-xi induced-gravity truncation. Stable zero vacuum energy and absolute Newton scale are not obtained.',
                boundary='Running nonminimal coupling, curvature counterterms, additional fields, two-loop logarithms or a supplied bare EH coefficient can change this truncation. The80-scalar loop coefficient is an application of standard Coleman-Weinberg theory, not a new general no-go theorem.',
                prior_owners=['w33_pass11316_native_edge_fierz_pauli.py','w33_pass11379_11383_full_frontier.py'],
                primary_source='https://doi.org/10.1103/PhysRevD.7.1888')


def payload():
    result={}
    for name,fn in [('11384_native_cp',native_cp_vacuum),('11385_hard_classes',hard_matching_classes),
                    ('11386_hub_gravity',hub_gravity_architecture),('11387_relaxed_scalar_wall',scalar_wall_relaxation),
                    ('11388_native_transmutation',native_dimensional_transmutation)]:
        result[name]=fn();print(name,result[name]['status'],flush=True)
    return dict(status='PASS',sections=result)


if __name__=='__main__':
    out=payload()
    (ROOT/'data/w33_pass11384_11388_native_dynamics.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
