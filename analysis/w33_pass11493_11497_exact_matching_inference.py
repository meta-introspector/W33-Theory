"""Exact conditional identities and native numerical controls; no TOE claim.

Prior owners: Pass11481-11485, Pass11271/11276/11281/11293, BT870,
Pass11480 and the classical Matrix-Tree and Schur-complement identities.
"""
import hashlib
import json
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.linalg import expm
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
SOURCE = 'data/w33_pass11481_11485_fibers_metric_yukawa_noisy.json'
OLD = 'data/w33_pass11476_11480_cubic_geometry_correlated.json'
OUT = ROOT / 'data/w33_pass11493_11497_exact_matching_inference.json'


def read(name):
    return json.loads((ROOT / name).read_text())


def dec(z):
    return np.array(z['real']) + 1j * np.array(z['imag'])


def enc(z):
    z = np.asarray(z)
    return {'real': z.real.tolist(), 'imag': z.imag.tolist()}


def stationary():
    """Hollow both actual Grams, then relax only hard fixed-stratum modes."""
    import w33_pass11481_11485_fibers_metric_yukawa_noisy as P
    import w33_pass11438_finite_native_model as F
    r = read(OLD)['vacuum']
    A, B = (dec(r['polydisc'][k]) for k in ['A_at_vacuum', 'B_at_vacuum'])
    E = dec(r['isotope']['native_tripotents'])
    grams = [a @ a.conj().T for a in [A, B]]
    _, V = np.linalg.eigh(grams[0])
    fourier = np.exp(2j * np.pi * np.outer(np.arange(3), np.arange(3)) / 3) / np.sqrt(3)
    initial = fourier @ V.conj().T
    initial *= np.exp(-1j * np.angle(np.linalg.det(initial)) / 3)
    L = P.su3()
    target = np.concatenate([np.repeat(np.trace(g).real / 3, 2) for g in grams])
    def unit(z):
        return expm(1j * np.einsum('a,aij->ij', z, L)) @ initial
    def residual(z):
        U = unit(z)
        return np.concatenate([np.diag(U @ g @ U.conj().T).real[:2] for g in grams]) - target
    opt = least_squares(residual, np.zeros(8), gtol=1e-14, ftol=1e-14,
                        xtol=1e-14, max_nfev=500)
    U = unit(opt.x)
    a, b = U @ A, U @ B
    # Isometric real embedding of the two3x3 matrices in the324-real carrier.
    T = np.kron(E, np.eye(3))
    R = np.block([[T.real, -T.imag], [T.imag, T.real]])
    embedding = np.zeros((324, 36)); embedding[:162, :18] = R; embedding[162:, 18:] = R
    z = np.r_[a.ravel().real, a.ravel().imag, b.ravel().real, b.ravel().imag]
    rows = []
    with F.native_context():
        def evaluate(z):
            value, gradient = F.fun(embedding @ z)
            return value, embedding.T @ gradient, float(np.linalg.norm(gradient))
        value0, gradient0 = F.fun(F.pack((E @ A).ravel(), (E @ B).ravel()))
        for tick in range(4):
            value, g, full = evaluate(z)
            step = 2e-5
            H = np.column_stack([(evaluate(z + step * v)[1] - evaluate(z - step * v)[1]) / (2 * step)
                                 for v in np.eye(36)])
            H = (H + H.T) / 2
            eigen, frame = np.linalg.eigh(H)
            hard = eigen > 1e-3
            delta = -frame[:, hard] @ ((frame[:, hard].T @ g) / eigen[hard])
            rows.append(dict(iteration=tick, value=value, full_gradient_norm=full,
                             restricted_gradient_norm=float(np.linalg.norm(g)), hard_rank=int(sum(hard)),
                             hard_gradient_norm=float(np.linalg.norm(frame[:, hard].T @ g)),
                             soft_gradient_norm=float(np.linalg.norm(frame[:, ~hard].T @ g)),
                             Newton_step_norm=float(np.linalg.norm(delta)), eigenvalues=eigen.tolist()))
            # Fixed hard-mode projection, explicit energy-decreasing line search.
            if tick < 3:
                candidates = [(evaluate(z + scale * delta)[0], scale) for scale in [1., .5, .25, .125, 0.]]
                _, scale = min(candidates)
                z += scale * delta
    # Exact elementary completed-square control, not an exact vacuum.
    x, y, a0, b0, c0, j = sp.symbols('x y a b c j', real=True)
    Vtoy = a0*x*x/2 + b0*x*y + c0*y*y/2 + j*y
    yeq = -(b0*x+j)/c0
    reduced = sp.factor(Vtoy.subs(y, yeq))
    assert sp.simplify(reduced - (a0*x*x/2 - (b0*x+j)**2/(2*c0))) == 0
    assert np.linalg.norm(residual(opt.x)) < 1e-9
    return dict(status='PASS', hollowisation_residual=float(np.linalg.norm(residual(opt.x))),
                hollowising_unitary=enc(U), initial_value=float(value0),
                initial_gradient_norm=float(np.linalg.norm(gradient0)), history=rows,
                final_reduced_coordinates=z.tolist(), exact_heavy_elimination=str(reduced),
                exact_shape_Hessian='A-B C^{-1} B^T, when C is invertible on the specified heavy slice',
                scope='Numerical simultaneous hollowisation and hard-mode relaxation of the supplied fixed stratum. Exact Schur-complement identity is conditional, not an exact or transverse/global stationary-vacuum certificate.')


def metric_determinant():
    """Actual cover, exact cofactor and scale identity, numerical shape force."""
    r = read(OLD)['geometry']; edges = np.array(r['cover_edges'])
    h = np.array(r['harmonic_FCC_displacements']); n = 80
    incidence = np.zeros((160, n), int)
    incidence[np.arange(160), edges[:, 0]] = 1
    incidence[np.arange(160), edges[:, 1]] = -1
    # At G=I, each supplied displacement is an integer multiple of1/80.
    hi = np.rint(80*h).astype(int)
    assert np.max(abs(h-hi/80)) < 1e-12
    assert np.all(abs(hi)==abs(hi[:,0,None])) and np.all(hi[:,0]!=0)
    lengths = [int(v @ v) for v in hi]
    Lap = sp.zeros(n)
    for (u, v), length in zip(edges, lengths):
        w = sp.Rational(6400, length)
        for i, j, sign in [(u,u,1),(v,v,1),(u,v,-1),(v,u,-1)]:
            Lap[int(i),int(j)] += sign*w
    # DomainMatrix exact rational determinant avoids generic expression elimination.
    cofactor = Lap[:-1,:-1].det(method='domain-ge')
    assert cofactor > 0
    small = np.asarray(Lap[:-1,:-1], float)
    sign, logdet = np.linalg.slogdet(small)
    assert sign == 1 and abs(logdet-float(sp.log(cofactor))) < 1e-9
    inv = np.linalg.inv(small)
    resistance = np.einsum('ei,ij,ej->e', incidence[:,:-1], inv, incidence[:,:-1])
    w = 6400/np.array(lengths, float)
    probabilities = w*resistance
    # W=.5 log cofactor; d logw=.5 Tr(dG)-h dG h/(h h).
    force = .5*np.einsum('e,eij->ij', probabilities,
                        .5*np.eye(3)[None,:,:] - np.einsum('ei,ej->eij',hi,hi)/np.array(lengths)[:,None,None])
    def action(G):
        weights = np.sqrt(np.linalg.det(G))/np.einsum('ei,ij,ej->e',h,G,h)
        C = incidence[:,:-1].T @ (weights[:,None]*incidence[:,:-1])
        return np.linalg.slogdet(C)[1]/2
    D = np.array([[.3,.07,-.02],[.07,-.1,.04],[-.02,.04,.2]])
    fd = (action(np.eye(3)+1e-5*D)-action(np.eye(3)-1e-5*D))/2e-5
    assert abs(fd-np.sum(force*D)) < 1e-7
    assert abs(probabilities.sum()-79) < 1e-9 and abs(np.trace(force)-79/4) < 1e-9
    shear=np.diag([np.exp(.3),np.exp(-.3),1.])
    shear_shift=action(shear)-action(np.eye(3))
    predicted=-79/2*np.log((np.exp(.3)+np.exp(-.3)+1)/3)
    assert abs(shear_shift-predicted)<1e-10
    return dict(status='PASS',vertices=n,edges=edges.tolist(),displacements_times80=hi.tolist(),
                squared_integer_lengths=lengths,exact_reduced_Laplacian_determinant=str(cofactor),
                cofactor_logdet=float(logdet),edge_tree_probabilities=probabilities.tolist(),
                metric_force=force.tolist(),direction=D.tolist(),directional_derivative=float(fd),
                exact_dilation_law='W(lambda G)-W(G)=(79/4) log(lambda)',
                exact_diagonal_metric_law='W(diag(g1,g2,g3))-W(I)=(79/2) log(3 sqrt(g1 g2 g3)/(g1+g2+g3))',
                exact_unimodular_shear_Hessian='d²W(exp(t diag(a,b,c)))/dt² at0 = -(79/6)(a²+b²+c²), a+b+c=0',
                shear_control_shift=float(shear_shift),shear_control_parameter=.3,
                shear_scope='Exact for the actual equal-magnitude-component displacement realization and supplied grounded scalar measure. On positive diagonal detG=1, AM-GM makes I a strict maximum of W, and W runs to-minus-infinity under unbounded anisotropy. Different regulator/measure, added gravity or fermion terms can change the action. Affine coordinate change must transform both metric and displacements.',
                exact_tree_probability_sum=79,scope='Classical Gaussian integral and weighted Matrix-Tree on the actual supplied80-vertex cover. Grounded scalar measure supplies W=.5 log cofactor. It produces a metric force but nonzero uniform-scale derivative, no scale stationary point without other terms, and no intrinsic curvature or Einstein action. Grounding, regulator and measure are inputs. BT870 owns the spanning-tree gravity route; this is a different actual-cover operator.')


def portal():
    import w33_pass11384_11388_native_dynamics as N
    _, tensors = N.native_tensors(); d = np.array(tensors[1], complex)
    assert np.max(abs(d.imag)) == 0 and np.max(abs(d.real-np.rint(d.real))) == 0
    E = dec(read(OLD)['vacuum']['isotope']['native_tripotents'])
    fields = [E @ dec(read(OLD)['vacuum']['polydisc'][key]) for key in ['A_at_vacuum','B_at_vacuum']]
    J = [np.einsum('abc,ai,bj->cij',d.conj(),p.conj(),p.conj()) for p in fields]
    Y = [np.einsum('abc,cij->aibj',d,j).reshape(81,81) for j in J]
    grams = [y@y.conj().T for y in Y]; comm = grams[0]@grams[1]-grams[1]@grams[0]
    # Completion of the square for positive M², no imposed linear spurion.
    hr, hi, jr, ji, kr, ki, m = sp.symbols('hr hi jr ji kr ki M2', real=True)
    hc, jc, kc = hr+sp.I*hi, jr+sp.I*ji, kr+sp.I*ki
    v = m*(hr**2+hi**2)-2*sp.re(sp.conjugate(hc)*kc*jc)
    square = m*sp.expand_complex(abs(hc-kc*jc/m)**2) - (kr**2+ki**2)*(jr**2+ji**2)/m
    assert sp.simplify(v-square) == 0
    # Prior11276 owns cubic SM-neutral restriction; stronger required pair restriction is checked here.
    exact = sp.MutableDenseNDimArray(np.rint(d.real).astype(int).tolist())
    flat=sp.Matrix(np.rint(d.real).astype(int).reshape(27,729).tolist())
    contraction=flat*flat.T
    contraction_constant=contraction[0,0]
    assert contraction==contraction_constant*sp.eye(27)
    sm_source_entries = [exact[a,b,c] for a,b in product(range(2),repeat=2) for c in range(27)]
    sm_nonzero = [(a,b,c,str(exact[a,b,c])) for a,b in product(range(2),repeat=2) for c in range(27) if exact[a,b,c] != 0]
    # Do not infer pair contraction zero from the prior cubic restriction alone.
    U=expm(1j*np.diag([.3,-.1,-.2]))
    p=fields[0]; pp=p@U.T
    covariance=np.linalg.norm(np.einsum('abc,ai,bj->cij',d.conj(),pp.conj(),pp.conj())-
                              np.einsum('ik,ckl,jl->cij',U.conj(),J[0],U.conj()))
    assert covariance<1e-10 and all(np.linalg.norm(y-y.T)<1e-10 for y in Y)
    return dict(status='PASS',source_matrices=[enc(j) for j in J],induced_Yukawas=[enc(y) for y in Y],
                singular_values=[np.linalg.svd(y,compute_uv=False).tolist() for y in Y],
                family_covariance_residual=float(covariance),
                Gram_commutator_norm=float(np.linalg.norm(comm)),cubic_CP_invariant=float(np.trace(comm@comm@comm).imag),
                exact_SM_pair_source_nonzero_entries=sm_nonzero,
                SM_pair_source_zero=not any(sm_source_entries),
                exact_matching='H=kappa J/M²; DeltaV=-|kappa|² ||J||²/M²; Y_eff=y kappa d J/M²',
                operator_dimension=5,exact_cubic_contraction_constant=str(contraction_constant),
                sufficient_boundedness='lambda_Phi >= '+str(contraction_constant)+' |kappa|²/M² ensures V=lambda_Phi||Phi||^4+M²||H||²-2Re(kappa Hdagger J) >=0; sufficient, not necessarily sharp',
                scope='Declared dynamical H(27,bar6) with quadratic composite J=conj(d) conj(Phi) conj(Phi). Exact algebraic heavy-field elimination for M²>0, and leading derivative EFT matching; kinetic inverse expands as (M²-D²)^{-1}. Extra field, coupling and scale supplied. Canonical SM pair source explicitly checked, not inferred from prior three-singlet zero. No SM projection, selected vacuum, observed mixing or physical CP prediction. Prior11271 owns sextet repair;11276/11281 own related composite obstructions.')


def gaussian_CP_control():
    """Exact Gaussian-integer quadratic-source witness, not physical CKM."""
    d=dec(read(SOURCE)['yukawa']['native_cubic'])
    assert np.max(abs(d.imag))==0 and np.max(abs(d.real-np.rint(d.real)))==0
    d=np.rint(d.real).astype(np.int64).astype(object)
    tri=read(OLD)['vacuum']['polydisc']['canonical_frame']
    E=np.zeros((27,3),dtype=object)
    for i,k in enumerate(tri[:3]): E[k,i]=tri[3] if i==2 else 1
    A=np.array([[1,1j,0],[0,2,1],[1,0,3]],complex)
    B=np.array([[2,1,1j],[1j,1,0],[1,2,2]],complex)
    def multiply(a,b):
        ar,ai=a; br,bi=b
        return ar@br-ai@bi, ar@bi+ai@br
    grams=[]
    for field in [A,B]:
        fr=E@field.real.astype(np.int64).astype(object)
        fi=E@field.imag.astype(np.int64).astype(object)
        jr=np.zeros((27,3,3),dtype=object);ji=jr.copy()
        for a,b,c in np.argwhere(np.asarray(d,dtype=int)!=0):
            coefficient=d[a,b,c]
            jr[c]+=coefficient*(np.outer(fr[a],fr[b])-np.outer(fi[a],fi[b]))
            ji[c]-=coefficient*(np.outer(fr[a],fi[b])+np.outer(fi[a],fr[b]))
        yr=np.zeros((81,81),dtype=object);yi=yr.copy()
        for i,j in product(range(3),repeat=2):
            yr[i::3,j::3]=sum(d[:,:,k]*jr[k,i,j] for k in range(27))
            yi[i::3,j::3]=sum(d[:,:,k]*ji[k,i,j] for k in range(27))
        grams.append(multiply((yr,yi),(yr.T,-yi.T)))
    p=multiply(grams[0],grams[1]);q=multiply(grams[1],grams[0]);C=(p[0]-q[0],p[1]-q[1])
    cube=multiply(multiply(C,C),C);real=np.trace(cube[0]);imag=np.trace(cube[1])
    assert isinstance(imag,int) and real==0 and imag!=0
    return dict(status='PASS',A=enc(A),B=enc(B),canonical_frame=tri,
                exact_trace_commutator_cube_real=str(real),exact_trace_commutator_cube_imag=str(imag),
                algebra='All tensor contractions and matrix products after input split use arbitrary-precision Python Gaussian-integer pairs. ImTr[Yu Yu†,Yd Yd†]^3 is basis-invariant and flips under conjugation.',
                scope='Exact81-state quadratic-source witness proves this supplied two-source architecture permits a nonzero CP-like cubic invariant. It is not the canonical SM-neutral plane, an up/down SM decomposition, physical CKM CP, selected vacuum or observed phase. Prior flavor audits own the Jarlskog diagnostic.')


def SM_mixed_portal():
    """Nonzero SM-neutral source v F; explicit added dynamical sextet."""
    d=dec(read(SOURCE)['yukawa']['native_cubic'])
    assert np.max(abs(d.imag))==0 and np.max(abs(d.real-np.rint(d.real)))==0
    D=sp.Matrix(np.rint(d[:,:,0].real).astype(int).tolist())
    # e0 is in the prior11271 canonical SM-neutral plane, reused rather than selected anew.
    F1=sp.diag(1,2,3)
    F2=sp.Matrix([[1,1,sp.I],[1,2,1],[sp.I,1,3]])
    g1=F1*F1.conjugate().T;g2=(F2*F2.conjugate().T).applyfunc(sp.expand)
    C=(g1*g2-g2*g1).applyfunc(sp.expand)
    cp=sp.expand(sp.trace(C**3));e6_factor=sp.trace((D*D.T)**3)
    ranks=[D.rank()*f.rank() for f in [F1,F2]]
    assert ranks==[30,30] and cp==-21600*sp.I and e6_factor==10
    return dict(status='PASS',SM_neutral_vector_index=0,
                family_sextets=[[[str(x) for x in row] for row in f.tolist()] for f in [F1,F2]],
                exact_E6_mass_rank=D.rank(),exact_full_mass_ranks=ranks,
                exact_family_cubic_invariant=str(cp),exact_full_cubic_invariant=str(e6_factor*cp),
                matching='Declare family-neutral v(27,1), F(1,bar6), H(27,bar6), J_Cij=v_C Fij. H*=kappa v F/M²; DeltaV=-|kappa|²||v||²||F||²/M². Scalar Hdagger v F is dimension3 and renormalizable with dimension1 kappa.',
                sufficient_boundedness='lambda_cross >= |kappa|²/M² for a positive lambda_cross ||v||²||F||² term, plus nonnegative stabilizing self quartics.',
                exact_mass_factorization='Y=d(v) tensor F; rankY=rank d(v)*rankF; Tr[H1,H2]^3=Tr[(d(v)d(v)dagger)^3]*Tr[F1F1dagger,F2F2dagger]^3 for common v.',
                scope='Explicit dynamical mixed-source matching remains nonzero at the prior canonical SM-neutral v=e0, so the quadratic-triplet source obstruction is bypassed by adding an independent E6-neutral family sextet. Two independently matched H/F copies give the CP-like capacity control. Prior11293 already owns this scalar mediator matching and its positive-square completion;11271 owns canonical rank30 exotic mass interface;11281 owns mixed v/H sextet relation and family selection obstruction. This reuses that matched interface and adds exact mass/CP controls, not a new mediator or physical quark masses/CKM. Higgs/sextet selection, scales, extra inventory, anomaly matching and loop matching remain open.')


def source_norm_fibers():
    """Exact row-Gram identity plus numerical replay on actual frame/fibers."""
    # On the exact canonical diagonal frame, d(Ep,Eq)=remaining Er for p!=q.
    raw=read(SOURCE);d=dec(raw['yukawa']['native_cubic'])
    tri=read(OLD)['vacuum']['polydisc']['canonical_frame']
    E0=np.zeros((27,3),int)
    for i,k in enumerate(tri[:3]):E0[k,i]=tri[3] if i==2 else 1
    for p,q in product(range(3),repeat=2):
        image=np.einsum('abc,a,b->c',d,E0[:,p],E0[:,q])
        expected=np.zeros(27) if p==q else E0[:,3-p-q]
        assert np.array_equal(image,expected)
    # Exact arbitrary-row Gaussian symbolic identity; all18 real coordinates.
    xr=sp.symbols('x0:9',real=True);xi=sp.symbols('y0:9',real=True)
    A=sp.Matrix(3,3,[a+sp.I*b for a,b in zip(xr,xi)]);G=A*A.conjugate().T
    norm=0
    for p,q in [(0,1),(0,2),(1,2)]:
        for i,j in product(range(3),repeat=2):
            z=A[p,i]*A[q,j]+A[q,i]*A[p,j]
            norm+=z*sp.conjugate(z)
    formula=sp.trace(G)**2+sp.trace(G*G)-2*sum(G[i,i]**2 for i in range(3))
    assert sp.expand(norm-formula)==0
    old=read(OLD)['vacuum'];E=dec(old['isotope']['native_tripotents']);a=dec(old['polydisc']['A_at_vacuum'])
    def value(a):
        phi=E@a;j=np.einsum('abc,ai,bj->cij',d.conj(),phi.conj(),phi.conj())
        g=a@a.conj().T
        rhs=np.trace(g).real**2+np.trace(g@g).real-2*np.sum(np.diag(g).real**2)
        return float(np.linalg.norm(j)**2),float(rhs)
    base,rhs=value(a);controls=[]
    for row in raw['fibers']['rows']:
        direct,prediction=value(dec(row['unitary'])@a)
        controls.append(dict(axis=row['axis'],t=row['t'],source_norm_shift=direct-base,identity_error=direct-prediction))
    assert abs(base-rhs)<1e-10 and max(abs(z['source_norm_shift']) for z in controls)<1e-10
    return dict(status='PASS',symbolic_coordinates=18,exact_identity='||J||²=(TrG)²+Tr(G²)-2 sum_i Gii², G=A A†',
                native_source_norm_squared=base,native_formula_residual=base-rhs,controls=controls,
                exact_matched_tree_consequence='-c||J||²=constant+2c sum_i Gii² along common-left SU3, c=|kappa|²/M²>0. Its minimum is uniform Gram diagonal; it is constant on all diagonal-preserving fibers.',
                scope='Exact arbitrary-complex3x3 identity on the canonical signed integer frame; actual numerical tripotent frame and all eight prior fibers replay it. The quadratic-triplet portal supplies a hollowisation preference but does not lift the prior potential-value fibers at tree level. No exact native stationary vacuum, frame conjugator, physical masses or selected loop coefficients follows. Prior11476 owns the fixed stratum;11481 owns the fibers; standard cubic-adjoint/Gram identities retain their classical ownership.')


def ward_cut():
    r=read(SOURCE)['ward']; J=np.array(r['anomaly_Jacobian']); sites=np.array(list(product(range(2),repeat=4)))
    rows=[]
    for radius in range(5):
        bounds=[]
        for ell in range(64):
            inside=np.sum(abs(sites-sites[ell//4]),axis=1)<=radius
            # Incidence of the induced ball has image all zero-sum vectors on the connected ball.
            square=np.dot(J[~inside,ell],J[~inside,ell])+J[inside,ell].sum()**2/sum(inside)
            bounds.append(float(np.sqrt(square)))
        actual=r['rows'][radius]['maximum_divergence_error']
        assert abs(max(bounds)-actual)<1e-13
        rows.append(dict(radius=radius,exact_projection_formula_numeric=max(bounds),stored_least_squares_residual=actual))
    charges=[1,-4,2,-3,6]; multiplicities=[6,3,3,2,1]
    moments={str(k):sum(m*q**k for q,m in zip(charges,multiplicities)) for k in [1,3,5,7]}
    assert moments['1']==moments['3']==0
    return dict(status='PASS',rows=rows,charge_moments=moments,
                exact_minimum_squared_defect='||J_out||²+(sum J_inside)²/|inside| for each connected induced ball',
                exact_proof='Incidence image is the zero-sum subspace on the connected ball, and zero outside. Orthogonal projection subtracts the mean inside and leaves every outside entry.',
                scope='Exact graph-support obstruction explains numerical Ward defects, including unavoidable outside tails. Charge anomaly moments cancel exactly but higher moments need not. This does not contradict infinite-volume local cohomology: exponential locality allows tails. No local gauge-invariant nonlinear anomaly primitive or measure has been built. Prior11484 owns the support experiment; Luscher owns anomaly-free local reconstruction.')


def joint_inference(error=Fraction(1,20), physical=Fraction(1,100), verification=Fraction(1,1000)):
    """Exact correlated two-block control: at most one shared Pauli fault.

    This is a supplied CPTP Pauli ensemble, not the native untwirled relay.
    Both blocks share the same X/Y/Z error and location; each block corrects it.
    Independent bit readout errors, noisy binary final verification, full flags.
    """
    e,p,v=map(Fraction,[error,physical,verification])
    assert 0<e<1 and 0<p<1 and 0<=v<=1
    syndromes=[0]+list(range(1,8))+[s<<3 for s in range(1,8)]+[s|(s<<3) for s in range(1,8)]
    # Integer common denominator: exact Bayesian comparisons and probabilities.
    prior=[1-p]+[p/21]*21
    prior_den=sp.ilcm(*[a.denominator for a in prior]); prior_num=[int(a*prior_den) for a in prior]
    en,ed=e.numerator,e.denominator
    confusion=[[en**((r^s).bit_count())*(ed-en)**(6-(r^s).bit_count()) for s in syndromes] for r in range(64)]
    single=[max(range(22),key=lambda s:prior_num[s]*confusion[r][s]) for r in range(64)]
    joint_numerator=separate_numerator=raw_numerator=0
    counts=[0]*22
    for a,b in product(range(64),repeat=2):
        scores=[prior_num[s]*confusion[a][s]*confusion[b][s] for s in range(22)]
        guess=max(range(22),key=lambda s:scores[s]); counts[guess]+=1
        joint_numerator+=scores[guess]
        if single[a]==single[b]:separate_numerator+=scores[single[a]]
        if a==b and a in syndromes:raw_numerator+=scores[syndromes.index(a)]
    den=int(prior_den)*ed**12
    joint=Fraction(joint_numerator,den);separate=Fraction(separate_numerator,den);raw=Fraction(raw_numerator,den)
    assert joint>=separate>=raw
    # False-positive test branches occupy an explicit orthogonal leakage sector.
    # Outputs: logical4, accepted-leakage flag1, rejected-erasure flag1.
    good=joint*(1-v);bad=(1-joint)*v;erasure=1-good-bad
    assert good+bad+erasure==1
    # A concrete likelihood-ratio explanation for memory-enabled inference.
    prior_ratio=21*(1-p)/p
    witness=3 # X error at the column with Hamming weight2.
    single_ratio=prior_ratio*(e/(1-e))**2
    joint_ratio=prior_ratio*(e/(1-e))**4
    witness_scores=[prior_num[s]*confusion[witness][s]**2 for s in range(22)]
    witness_joint=syndromes[max(range(22),key=lambda s:witness_scores[s])]
    return dict(status='PASS',readout_error=str(e),physical_shared_fault_probability=str(p),verification_flip=str(v),
                common_syndromes=syndromes,prior_fractions=[str(a) for a in prior],joint_guess_counts=counts,
                exact_joint_success=str(joint),exact_separate_success=str(separate),exact_raw_success=str(raw),
                joint_success=float(joint),separate_success=float(separate),raw_success=float(raw),
                logical_weight=str(good),accepted_leakage_weight=str(bad),rejected_erasure_weight=str(erasure),
                two_report_witness=dict(reports=[witness,witness],separate_guess=syndromes[single[witness]],
                                        joint_guess=witness_joint,identity_vs_fault_single_odds=str(single_ratio),
                                        identity_vs_fault_joint_odds=str(joint_ratio)),
                exact_output_channel='rho -> good*rho direct_sum bad*Tr(rho)|L><L| direct_sum erasure*Tr(rho)|E><E|',
                scope='Exact rational two-Steane-block shared one-Pauli-fault CPTP control with12 noisy readout bits and a noisy final syndrome-mismatch test. Independent readout draws but correlated physical faults retained; MAP success optimized for supplied prior. False accepts retained as leakage flags, no postselection. Flag output coarse-grains erroneous-syndrome sectors by tracing out encoded logical information; it is a supplied CPTP postprocessor, not an isometry into a one-dimensional flag. Not native untwirled11480 channel, arbitrary multiple faults, leakage hardware or threshold.')


def run():
    result={}
    for name, fn in [('stationary',stationary),('metric',metric_determinant),('portal',portal),('ward',ward_cut),('joint',joint_inference),('exact_CP_control',gaussian_CP_control),('SM_mixed_portal',SM_mixed_portal),('source_norm_fibers',source_norm_fibers)]:
        result[name]=fn(); print(name,'PASS',flush=True)
        Path('/tmp/w33_11493_partial.json').write_text(json.dumps(result))
    result.update(status='PASS',passes=list(range(11493,11498)),reservation='2cda9139b',
                  source_sha256={p:hashlib.sha256(json.dumps(read(p),sort_keys=True,separators=(',',':')).encode()).hexdigest() for p in [SOURCE,OLD,'data/w33_pass11293_yukawa_uv_budget.json']})
    OUT.write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    run()
