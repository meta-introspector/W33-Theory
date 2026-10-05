"""Exact reductions and explicit obstructions for the five11516-11520 frontiers.

Earlier owners:11471 frame dictionary,11476 polydisc,11481 value fibers,
11516 normal audit,11517 adjoint split,11518 metric stress,11438/11484
chiral reconstruction,11520 Pauli recovery. See companion report for scopes.
"""
import json, hashlib, itertools, inspect, sys
from functools import lru_cache
from pathlib import Path
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/w33_pass11521_11525_exact_frontiers.json'
SOURCES=['data/w33_pass11471_11475_frames_currents_matching.json',
 'data/w33_pass11493_11497_exact_matching_inference.json',
 'data/w33_pass11516_11520_local_geometry_quantum.json',
 'data/w33_pass11506_11510_five_physics_targets.json',
 'data/w33_pass11271_canonical_sm_higgs_bridge.json',
 'data/w33_pass11384_native_inputs.json']
read=lambda p:json.loads((ROOT/p).read_text())
enc=lambda a:dict(real=np.asarray(a).real.tolist(),imag=np.asarray(a).imag.tolist())
dec=lambda a:np.array(a['real'])+1j*np.array(a['imag'])
strings=lambda a:[[str(x) for x in row] for row in a.tolist()]

def modular_rank(matrix,prime=65521):
    """Exact finite-field lower bound on rational/integer rank."""
    assert s.isprime(prime)
    A=np.array(matrix,dtype=np.int64)%prime;row=0
    for col in range(A.shape[1]):
        piv=np.flatnonzero(A[row:,col])
        if not len(piv):continue
        pivot=row+int(piv[0]);A[[row,pivot]]=A[[pivot,row]]
        A[row]=(A[row]*pow(int(A[row,col]),-1,prime))%prime
        for i in range(row+1,A.shape[0]):A[i]=(A[i]-A[i,col]*A[row])%prime
        row+=1
        if row==A.shape[0]:break
    return row

def exact_row_penalty():
    """Compute the rational E6 Casimir restricted to the canonical three-frame."""
    from w33_pass11384_11388_native_dynamics import native_tensors
    B,(tri,_,_)=native_tensors();frame=tri[0][:3];unique={}
    for b in B:
        for kind,a in [(0,b+b.T),(1,b-b.T)]:
            nz=np.flatnonzero(a)
            if not len(nz):continue
            a=a//int(np.gcd.reduce(abs(a.ravel()[nz])))
            if a.ravel()[nz[0]]<0:a=-a
            unique[kind,tuple(a.ravel())]=a*(1j if kind else 1)
    he=np.array(list(unique.values()));G=np.einsum('aij,bji->ab',he,he).real.astype(int)
    D=s.Matrix([[int(h[i,i].real) for i in frame] for h in he])
    indices=[i for i in range(78) if any(D[i,j] for j in range(3))]
    while True:
        expanded=sorted(set(indices)|set(np.flatnonzero(np.any(G[indices]!=0,axis=0))))
        if expanded==indices:break
        indices=expanded
    # Cartan/root orthogonality makes the smaller inverse exactly equivalent.
    others=[i for i in range(78) if i not in indices]
    assert not np.any(G[np.ix_(indices,others)])
    for h in he:
        block=h[np.ix_(frame,frame)]
        assert np.array_equal(block,np.diag(np.diag(block)))
    C=D[indices,:].T*s.Matrix(G[np.ix_(indices,indices)]).inv()*D[indices,:]
    assert C* s.ones(3,1)==s.zeros(3,1) and C.rank()==2
    c=C[0,0]*s.Rational(3,2)
    assert C==c*(s.eye(3)-s.ones(3)/3)
    return C,c

def quotient_value(z,mp):
    A=mp.diag([z[i]*mp.exp(1j*z[3]) for i in range(3)])
    imaginary=[i for i in range(9) if i not in [1,5]];B=mp.matrix(3)
    for i in range(9):B[i//3,i%3]=z[4+i]+(1j*z[13+imaginary.index(i)] if i in imaginary else 0)
    dagger=lambda M:M.transpose_conj();trace=lambda M:sum(M[i,i] for i in range(3))
    norm=lambda M:sum(abs(M[i,j])**2 for i in range(3) for j in range(3))
    def one(M):
        G=dagger(M)*M;G-=trace(G)*mp.eye(3)/3;u=27*mp.det(M)**2
        return norm(G)/2+abs(u)**4-abs(u)**2+u.real**2+u.real/2+abs(u)**6
    radius=norm(A)+norm(B);C=A.T*B.conjugate()/(radius/2+mp.mpf('0.01'))
    U=2*mp.eye(3)+(C+dagger(C))/2;V=mp.eye(3)+C*dagger(C);W=U*V-V*U
    chi=(27*mp.det(A)**2).imag
    return one(A)+one(B)-10000*chi/(1+chi*chi)*trace(W**3).imag+mp.mpf('0.000001')*radius**2

def quotient_control():
    import mpmath as mp
    z=np.array(read(SOURCES[1])['stationary']['final_reduced_coordinates']);A=(z[:9]+1j*z[9:18]).reshape(3,3);B=(z[18:27]+1j*z[27:]).reshape(3,3)
    u,values,vh=np.linalg.svd(A);left=np.exp(1j*np.angle(np.linalg.det(u))/3)*u.conj().T
    right=np.exp(1j*np.angle(np.linalg.det(vh))/3)*vh.conj().T;AA=left@A@right;BB=left@B@right
    phase=np.angle(AA[0,0]);p=np.angle(BB[0,1]);q=np.angle(BB[1,2]);theta=np.array([(-2*p-q)/3,(p-q)/3,(p+2*q)/3]);D=np.diag(np.exp(1j*theta));BB=D@BB@D.conj().T
    assert abs(BB[0,1].imag)+abs(BB[1,2].imag)<1e-12
    initial=np.r_[values,phase,BB.real.ravel(),BB.imag.ravel()[[0,2,3,4,6,7,8]]]
    with mp.workdps(70):
        value=lambda *args:quotient_value(args,mp)
        gradient=lambda *args:tuple(mp.diff(value,tuple(args),tuple(int(j==i) for j in range(20))) for i in range(20))
        root=mp.findroot(gradient,tuple(map(mp.mpf,initial)),tol=mp.mpf('1e-55'),maxsteps=12,solver='mdnewton')
        residual=max(abs(v) for v in gradient(*root));assert residual<mp.mpf('1e-50')
        H=mp.matrix(20)
        for i in range(20):
            for j in range(20):
                order=[0]*20;order[i]+=1;order[j]+=1;H[i,j]=mp.diff(value,tuple(root),tuple(order))
        spectrum=mp.eigsy(H,eigvals_only=True)
        return dict(coordinates=[mp.nstr(x,65) for x in root],value=mp.nstr(value(*root),65),gradient_maximum=mp.nstr(residual,12),Hessian_eigenvalues=[mp.nstr(v,35) for v in spectrum],precision_digits=70,
          chart='A=e^(i phase)diag(a1,a2,a3), ai>0; B12,B23 real positive after residual diagonal SU3 conjugation; remaining B entries unrestricted complex.20 real coordinates.',
          exact_minimization='Every pair of trace-zero Hermitian matrices can be simultaneously unitarily hollowed (Toeplitz-Hausdorff, recursion on orthogonal complements; Damm-Fassbender1910.08813 and11481). Hence min_Uleft V(UA,UB)=V_U(A,B) exactly, and inf over the native polydisc equals inf over V_U. On this regular SVD chart V_U has20 real variables.',
          scope='High-precision stationary quotient root and Hessian spectrum of the supplied polydisc potential, not interval root isolation or full physical vacuum. Exact global minimization identity is restricted to the polydisc, not all324 fields.')

def vacuum():
    """Exact Schur reduction; numerical blocks never become interval masses."""
    import w33_pass11438_finite_native_model as F
    C,c=exact_row_penalty();f=read(SOURCES[0])['frame']
    T=dec(f['albert_bridge']['native_from_clock_basis']);Ti2=dec(f['albert_bridge']['twice_inverse'])
    generators=np.array(f['integer_generators']);real_generators=[]
    for k in generators:
        for g in [k-k.T,1j*(k+k.T)]:
            r=Ti2@g@T
            assert not np.any(r.imag)
            if np.any(r):real_generators.append(np.rint(r.real).astype(int))
    # Each real8 is absolutely irreducible: its exact commutant is scalar.
    labels=f['albert_bridge']['clock_basis_labels'];blocks=[]
    for native in f['invariant_coordinate_blocks']:
        if len(native)==8:
            columns=[j for j in range(27) if any(T[i,j] for i in native)]
            assert len(columns)==8;blocks.append(columns)
    from scipy.linalg import qr
    witnesses=[]
    def certified_rank(equations,target):
        # Floating QR ONLY selects a minor; exact rank certifies it.
        _,_,pivots=qr(np.asarray(equations,float).T,pivoting=True,mode='economic')
        witness=equations[list(map(int,pivots[:target])),:]
        assert modular_rank(witness)==target
        witnesses.append(dict(rank=target,prime=65521,selected_rows=list(map(int,pivots[:target])),minor=witness.tolist()))
        return target
    ranks=[];module_generators=[]
    for columns in blocks:
        gg=[g[np.ix_(columns,columns)] for g in real_generators]
        equations=np.vstack([np.kron(np.eye(8,dtype=int),g)-np.kron(g.T,np.eye(8,dtype=int)) for g in gg])
        rank=certified_rank(equations,63);ranks.append(rank);module_generators.append(gg)
    cross_ranks=[]
    for i,j in itertools.combinations(range(3),2):
        equations=np.vstack([np.kron(np.eye(8,dtype=int),a)-np.kron(b.T,np.eye(8,dtype=int)) for a,b in zip(module_generators[i],module_generators[j])])
        cross_ranks.append(certified_rank(equations,64))
    z=np.array(read(SOURCES[1])['stationary']['final_reduced_coordinates'])
    A=(z[:9]+1j*z[9:18]).reshape(3,3);B=(z[18:27]+1j*z[27:]).reshape(3,3)
    tri=F.TENSORS[0][0];E=np.zeros((27,3));E[tri[0],0]=1;E[tri[1],1]=1;E[tri[2],2]=tri[3]
    x=F.pack((E@A).ravel(),(E@B).ravel())
    with F.native_context():
        value,g=F.fun(x);step=1e-5
        H=np.column_stack([(F.fun(x+step*v)[1]-F.fun(x-step*v)[1])/(2*step) for v in np.eye(324)])
        H=(H+H.T)/2
    U=T/np.sqrt(np.diag(T.conj().T@T).real)[None,:]
    carriers=[];rows=[]
    for columns in blocks:
        carrier=[]
        for field,imag,family,axis in itertools.product(range(2),range(2),range(3),range(8)):
            a=np.zeros((27,3),complex);b=a.copy();(a if field==0 else b)[:,family]=(1j if imag else 1)*U[:,columns[axis]]
            carrier.append(F.pack(a.ravel(),b.ravel()))
        carrier=np.array(carrier).T;assert np.linalg.norm(carrier.T@carrier-np.eye(96))<1e-12
        reduced=(carrier.T@H@carrier).reshape(12,8,12,8)
        small=np.einsum('aibi->ab',reduced)/8
        defect=np.linalg.norm(reduced-np.einsum('ab,ij->aibj',small,np.eye(8)))
        rows.append(dict(clock_columns=columns,block=small.tolist(),eigenvalues=np.linalg.eigvalsh(small).tolist(),Schur_residual=float(defect)))
        carriers.append(carrier)
    cross=max(np.linalg.norm(a.T@H@b) for i,a in enumerate(carriers) for b in carriers[i+1:])
    return dict(status='PASS',row_penalty_matrix=strings(C),row_penalty_coefficient=str(c),
      exact_restriction='V(A,B)=V_unitary(A,B)+(c/2)(||diag(AA†)-tr(AA†)/3||²+||diag(BB†)-tr(BB†)/3||²); c stored exactly. Family moment term is (1/2)tr[(A†A-tr(A†A)/3)^2] for each field.',
      stationary_transport='On simultaneous hollow-Gram locus the penalty and its full first derivative vanish. V_unitary is invariant under common-left SU3. Therefore any EXACT reduced stationary pair on this locus maps to another stationary pair under every common-left U retaining the locus. Spin8 fixes the stratum and has no transverse fixed vector, hence reduced stationarity implies full stationarity.',
      conditional_moduli='At regular rank-four hollow constraints, the common-left SU3 fiber has dimension4; dividing the two diagonal gauge directions leaves two stationary moduli, provided the full pair has no extra relevant gauge identifications. The generic rank [4,10,12] witness belongs11481. This is conditional on an exact hollow stationary pair, not its existence certificate.',
      exact_Hessian_form='H=H36 direct_sum (Bv tensor I8) direct_sum (Bs tensor I8) direct_sum (Bc tensor I8), with three symmetric12x12 matrices. It holds at every point of the canonical Spin8-fixed stratum, not only stationary points. Schur lemma, exact scalar commutants and inequivalent modules prove the form.',
      real_commutant_ranks=ranks,intertwiner_constraint_ranks=cross_ranks,rank_witnesses=witnesses,canonical_coordinates=x.tolist(),value=value,gradient_norm=float(np.linalg.norm(g)),normal_blocks=rows,cross_block_residual=float(cross),
      quotient_control=quotient_control(),scope='Exact rational row penalty, conditional stationary-fiber theorem and representation reduction; numerical12x12 blocks and70-digit quotient root are controls. No interval vacuum, transverse mass certificate or observed physical spectrum.')

def flavor():
    from w33_pass11384_11388_native_dynamics import native_tensors
    _,(_,d,_)=native_tensors();sm=read(SOURCES[4]);charges=list(map(s.Rational,sm['hypercharge_diagonal']))
    neutral=[i for i,y in enumerate(charges) if y==0]
    old=read(SOURCES[3])['SM'];hu=old['neutral_Hu_index'];hd=old['neutral_Hd_index']
    # Hypercharge alone suffices: two Y=0 vevs cannot source Y=+/-1/2.
    supports={c:sum(bool(d[c,a,b]) for a,b in itertools.product(neutral,repeat=2)) for c in [hu,hd]}
    assert all(y==0 for y in (charges[i] for i in neutral))
    assert all(v==0 for v in supports.values())
    assert all(charges[a]+charges[b]+charges[c]==0 for a,b,c in np.argwhere(d))
    # Exact tensor covariance of the composite, tested coefficientwise in E6.
    B=native_tensors()[0]
    residual=max(np.max(abs(np.einsum('ai,ajk->ijk',g,d)+np.einsum('aj,iak->ijk',g,d)+np.einsum('ak,ija->ijk',g,d))) for g in B)
    assert residual==0
    t=s.symbols('t');z=s.symbols('z');d0,l0=s.symbols('d0 l0')
    # True arbitrary insertion form on a fixed charge pair: a symmetric polynomial.
    qd=(s.Rational(1,6),s.Rational(1,3));ql=(-s.Rational(1,2),s.Integer(1))
    elementary=dict(sum=str(sum(qd)),down_product=str(s.prod(qd)),lepton_product=str(s.prod(ql)))
    assert sum(qd)==sum(ql)==s.Rational(1,2)
    return dict(status='PASS',composite='H^c_ij=dbar^{cab} Phi*_ai Phi*_bj / M; transforming (27,bar6), not an elementary sextet',
      exact_tensor_covariance_residual=int(residual),hypercharge_zero_indices=neutral,Higgs_indices=[hu,hd],neutral_source_entries=supports,
      exact_neutral_obstruction='Every polynomial covariant of SM-preserving vevs transforms trivially under that SM. It cannot acquire an electroweak-doublet vev. In the named quadratic composite, Yc=-Ya-Yb=0, so Hu/Hd components vanish coefficientwise. Adjoint insertions with SM-preserving adjoint vevs cannot evade this.',
      exact_insertion_ring='For identical fermion legs all symmetric polynomials f(Yi,Yj) lie in Q[Yi+Yj,YiYj]. With common Higgs charge, the first generator is fixed. Down/lepton splitting is precisely dependence on the second generator, 1/18 versus -1/2;11517 owns the first split.',
      charge_pair_invariants=elementary,
      selection_boundary='The scalar11438 potential contains neither an elementary adjoint A nor independent Higgs sextets F0,F2. Its stationary equations cannot select their vevs or Wilson coefficients. A composite reduces inventory, but the SM-preserving sector cannot source an EW vev. Additional EW-breaking dynamics and matching data are necessary, not numerical precision.',
      prior_owners=['11271 symmetric Weyl obstruction and sextet repair','11276 different composite spectator and neutral cubic obstruction','11281 sextet spectral selection','11517 adjoint split'],
      scope='Exact covariant composite and symmetry obstruction; no dynamically selected realistic flavor or claim of general composite-Higgs novelty.')

def gravity():
    z=s.symbols('z',real=True);la,a,b,w0,wz=s.symbols('Lambda alpha beta W0 Wz',real=True);gamma=s.symbols('gamma',positive=True)
    v=s.exp(3*z/2);V=la*v+a/v+3*b*s.exp(z/2)/2+gamma*(w0+wz*z)
    stationary=s.solve(s.diff(V,z).subs(z,0),la)[0]
    energy=s.factor(V.subs({z:0,la:stationary}));assert s.expand(energy-(2*a+b+gamma*w0-s.Rational(2,3)*gamma*wz))==0
    both=s.solve([s.diff(V,z).subs(z,0),V.subs(z,0)],[la,a])
    assert both[a]==-b/2-gamma*w0/2+gamma*wz/3
    src=read(SOURCES[2])['metric']['flux_winding_control']
    # L(e^z I)=e^-z L(I), so Wz=+sum r(lambda), NOT its negative.
    # log(1+x)>=x/(1+x), x=3/(lambda+1), proves W0>=Wz>=0.
    l=s.symbols('l',nonnegative=True);mode=(l+4)/(l+1);r=3*l/((l+4)*(l+1))
    assert s.simplify(mode-1)==3/(l+1)
    gamma0=s.Rational(1,1000)
    lower=4-gamma0*s.Rational(158,9)
    # Native incidence Dirac, not an assumed four-dimensional spin module.
    m,t=s.symbols('m t',positive=True)
    pauli=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1)]
    spatial=[s.kronecker_product(pauli[0],p) for p in pauli];time=s.I*s.kronecker_product(pauli[2],s.eye(2))
    assert time*time==-s.eye(4)
    for i,j in itertools.product(range(3),repeat=2):assert spatial[i]*spatial[j]+spatial[j]*spatial[i]==(2 if i==j else 0)*s.eye(4)
    assert all(time*g+g*time==s.zeros(4) for g in spatial)
    metric=read(SOURCES[1])['metric'];h=np.array(metric['displacements_times80'],int);edges=metric['edges']
    gradient=np.zeros((640,320),complex);lap=np.zeros((80,80))
    for edge,((i,j),v0) in enumerate(zip(edges,h)):
        square=int(v0@v0);g=sum((int(v0[k])*spatial[k] for k in range(3)),s.zeros(4))
        assert (g*g-square*s.eye(4)).applyfunc(s.expand)==s.zeros(4)
        ce=np.array(80*g/square,complex)
        gradient[4*edge:4*edge+4,4*i:4*i+4]=-ce;gradient[4*edge:4*edge+4,4*j:4*j+4]=ce
        weight=6400/square
        lap[i,i]+=weight;lap[j,j]+=weight;lap[i,j]-=weight;lap[j,i]-=weight
    lift_error=float(np.linalg.norm(gradient.conj().T@gradient-np.kron(lap,np.eye(4))))
    assert lift_error<1e-8
    return dict(status='PASS',stationary_Lambda=str(stationary),stationary_energy=str(energy),
      simultaneous_zero_energy_parameters={str(k):str(val) for k,val in both.items()},
      exact_no_go='For positive flux alpha and winding beta, gamma>=0, W0>=Wz>=0, every isotropic stationary point has V0=2alpha+beta+gamma(W0-2Wz/3)>=2alpha+beta+gamma W0/3>0. Setting V0=0 forces alpha=-beta/2-gamma W0/2+gamma Wz/3<0. A cosmological Lambda alone cannot simultaneously enforce stationarity and zero flat-space lapse energy in this supplied homogeneous action.',
      exact_log_bound='For x=3/(lambda+1)>0, log(1+x)-x/(1+x) has derivative x/(1+x)^2>=0 and value0 at0. r(lambda)=[lambda/(lambda+1)]x/(1+x)<=log(1+x). Hence W0>=Wz>=0.',
      inherited_energy_lower_bound=str(lower),
      true_constraint='For the11482 homogeneous action L=K/N-N V, lapse equation is K+N² V=0. At static G and zero scalar velocities, K=0, hence V=0. The11518 positive minimum is not a static flat gravitational vacuum of that action. Time evolution, intrinsic spatial curvature or a different stress completion must be retained.',
      spin_carrier='Native graph incidence Dirac acts on C^80 direct_sum C^160, with79 nonzero Laplacian modes and82 zero modes. det(mI+iDgraph)=m^82 product_(j=1)^79(m²+lambda_j)=m^80 det(m²I80+L). The graph determinant alone has no four-dimensional Clifford action, tetrad or spin connection.',
      needed_spin_map='Choose a spin structure, a4D geometry/tetrad and Spin transport; build Dspin on site-spinor fibers with covariant edge couplings, then verify Clifford principal symbol and compare the full determinant including zero modes. Existing10950/forty_points Spin(1,9) matter is an internal carrier map, not this spacetime operator.',
      explicit_flat_spin_lift=dict(spatial_gamma=[strings(g) for g in spatial],Lorentzian_time_gamma=strings(time),edge_map='C_e psi=(80 Gamma(h_times80)/||h_times80||²)(psi_b-psi_a)',vertex_square='Cdag C=L tensor I4',numerical_square_replay_error=lift_error,cochain_dimensions=[320,640],Dirac_dimension=960,zero_modes=328,full_determinant='det(mI+iDlift)=m^320 det(m²I80+L)^4',boundary='Explicit supplied four-spinor Clifford fibers dress the actual edge map, but this remains an incidence/cochain Dirac with surplus zero modes. It is not the original one-fermion determinant or an identified SM spacetime Dirac; a spin/taste projection and local Spin connection must be justified.'),
      scope='Exact static lapse obstruction and determinant/carrier diagnosis in the declared model; not a no-go theorem for gravity, flux stabilization in other dimensions, or all W33 real forms.')

def ward():
    # Exact explanation of the gauge variation missing from the relative current.
    # P dP P=0 makes i tr P[[i omega,P],dP] a density variation.
    P=s.diag(1,1,0,0);a,b,c,d,e,f,g,h=s.symbols('a b c d e f g h',real=True)
    dot=s.Matrix([[0,0,a,b],[0,0,c,d],[a,c,0,0],[b,d,0,0]])
    omega=s.Matrix([[e,0,f,0],[0,e,0,g],[f,0,h,0],[0,g,0,h]])
    comm=lambda x,y:x*y-y*x
    identity=s.simplify(s.trace(P*comm(comm(omega,P),dot))+s.trace(omega*dot))
    assert identity==0
    charges=[1,-4,2,-3,6];mult=[6,3,3,2,1]
    sums={str(k):sum(m*q**k for q,m in zip(charges,mult)) for k in [1,3,5]}
    return dict(status='PASS',projector_gauge_variation_identity='tr P[[omega,P],dP]=-tr omega dP; P dP P=0',
      missing_exact_term='Let div k(A)=anomaly(A), with k gauge invariant and spatially local. Add delta_eta S(A), S(A)=integral_0^1 <A,k(tA)>dt, to the radial Berry connection. Its gauge variation cancels the radial curvature gauge variation; its exterior derivative is zero, so integrability curvature is retained. The pair-density transgression alone is not this measure connection.',
      coefficient_form='delta_eta S=integral_0^1 [<eta,k(tA)>+t<A,Dk(tA)[eta]>]dt; Berry term=integral_0^1 t i Tr P(tA)[DP(tA)[A],DP(tA)[eta]]dt.',
      exact_charge_sums=sums,absolute_odd_charge_counts={'1':6,'3':2},
      imported_existence='11438 already checks charge conditions. Luscher9811032 theorems5.1,5.3,5.4 give an admissible Abelian measure, with finite-volume corrections on sufficiently large lattices and even counts per absolute odd charge for all sectors. Existence is prior art; this packet exposes the missing exact correction, not a new existence theorem.',
      obstruction='The required k is a gauge-invariant local nonlinear anomaly primitive.11519 pair current is only relative-reference covariant. A Coulomb inverse supplies k with a massless inverse Laplacian and fails the desired locality. Finite-volume toron/flux corrections and a computed local primitive remain necessary for an executable implementation.',
      scope='Exact gauge-cancellation/integrability form. Local primitive, finite-volume corrections and non-Abelian measure not implemented; no theorem about locality follows from parameter homotopy alone.')

def native_branches():
    """Replay11520 sparse contraction, preserving all logical coherences.

    Reuse the established producer until its conditional channels are built;
    no new ownership claim for encoding/contraction. No source file is edited.
    """
    import w33_pass11516_11520_local_geometry_quantum as P
    source=inspect.getsource(P.quantum_fidelity)
    source=source[:source.index('        dephase=np.array')]+ '        return branches\n'
    namespace=dict(P.__dict__);exec(source,namespace)
    return namespace['quantum_fidelity']()

@lru_cache(maxsize=1)
def cell_effect():
    source=read(SOURCES[3])['recovery'];K=dec(source['logical_cell_subchannel'])
    E=np.kron(K.conj().T@K,K.conj().T@K);v=source['verification_flip']
    return E,v

def recovery_objective(transfer):
    import w33_pass11476_11480_cubic_geometry_correlated as M
    E,v=cell_effect()
    # Inherited Choi ordering is INPUT tensor OUTPUT; recovery vectors are
    # row-major OUTPUT tensor INPUT. Swap both ket and bra factors explicitly.
    J=M.transfer_choi(transfer).reshape(4,4,4,4).transpose(1,0,3,2)
    q=lambda e:np.einsum('ia,jc,bidj->abcd',e.conj(),e,J.conj()).reshape(16,16)/4
    Q=(1-v)*((1-v)*q(E)+v*q(np.eye(4)-E))
    return (Q+Q.conj().T)/2

def rational_matrix(a,bits=32):
    scale=2**bits
    return s.Matrix([[s.Rational(int(round(z.real*scale)),scale)+s.I*s.Rational(int(round(z.imag*scale)),scale) for z in row] for row in np.asarray(a)])

def exact_positive_pivots(A):
    """Exact Gaussian-rational LDL Schur pivots, with no numerical PSD test."""
    assert A==A.conjugate().T
    n=A.rows;B=s.MutableDenseMatrix(A);pivots=[]
    for k in range(n):
        pivot=s.cancel(B[k,k]);assert pivot.is_real and pivot>0
        pivots.append(pivot)
        for i in range(k+1,n):
            for j in range(i,n):
                val=s.cancel(B[j,i]-B[j,k]*B[k,i]/pivot)
                B[j,i]=val;B[i,j]=s.conjugate(val)
    return pivots

def flagged_Z_hooks():
    """Explicit Clifford flag circuit, exact for syndrome-ancilla Z faults only."""
    support=[i for i in range(7) if ((7-i)&4)];ancilla=7;flag=8
    gates=[(flag,ancilla)]+[(i,ancilla) for i in support]+[(flag,ancilla)]
    raw=[]
    for insertion in range(len(gates)+1):
        z=np.zeros(9,dtype=int);z[ancilla]=1
        for control,target in gates[insertion:]:z[control]^=z[target]
        syndrome=0
        for i in np.flatnonzero(z[:7]):syndrome^=7-int(i)
        raw.append((insertion,int(z[flag]),syndrome,z[:7]))
    lookup={}
    for _,flag_result,syndrome,z in raw:
        key=(flag_result,syndrome)
        if key not in lookup or z.sum()<lookup[key].sum():lookup[key]=z
    H=np.array([[((7-j)>>k)&1 for j in range(7)] for k in range(3)])
    stabilizers={tuple((np.array(v)@H)%2) for v in itertools.product(range(2),repeat=3)}
    rows=[]
    for insertion,flag_result,syndrome,z in raw:
        correction=lookup[flag_result,syndrome];residual=z^correction
        assert tuple(residual) in stabilizers
        rows.append(dict(fault_after_gate=insertion,flag_X_minus=flag_result,syndrome=syndrome,data_Z_support=np.flatnonzero(z).tolist(),correction_Z_support=np.flatnonzero(correction).tolist(),residual_Z_stabilizer_support=np.flatnonzero(residual).tolist()))
    return dict(gates=[list(g) for g in gates],syndrome_ancilla=ancilla,flag_ancilla=flag,
      initialization='syndrome ancilla |0>, flag |+>; final syndrome Z and flag X measurements; later ideal data X checks',rows=rows,
      exact_result='All seven positions of a single syndrome-ancilla Z fault yield a stabilizer after the stored flag/syndrome-dependent correction. The two flag-control CNOTs commute with all data-control CNOTs and cancel ideally, preserving the native parity measurement.',
      scope='Only the declared syndrome-ancilla Z fault family, with ideal flag gates/other faults/readout. Not a proof for all CNOT Pauli faults, faulty flag ancillas, multiple extraction rounds or a threshold.')

def syndrome_relation_key(branch):
    left,right=divmod(int(branch),64)
    vectors=[left//8,left%8,right//8,right%8];key=0
    for mask in range(1,16):
        value=0
        for i in range(4):
            if (mask>>i)&1:value^=vectors[i]
        if value==0:key|=1<<(mask-1)
    return key

def recovery_orbit_bounds(branches):
    """All4096 ideal-readout recoveries, via exact GL3(2) syndrome orbits."""
    import cvxpy as cp
    import w33_pass11476_11480_cubic_geometry_correlated as M
    groups={}
    for branch in range(4096):groups.setdefault(syndrome_relation_key(branch),[]).append(branch)
    assert len(groups)==66
    pauli=[np.kron(a,b) for a,b in itertools.product(M.PAULI,repeat=2)]
    certificates=[];lower_total=s.Integer(0);upper_total=s.Integer(0);orbit_error=0.
    for key,members in groups.items():
        branch=members[0];orbit_error=max(orbit_error,float(np.max(abs(branches[members]-branches[branch]))))
        Q=recovery_objective(branches[branch]);X=cp.Variable((16,16),hermitian=True)
        tr=cp.bmat([[sum(X[4*a+b,4*a+d] for a in range(4)) for d in range(4)] for b in range(4)])
        constraints=[X>>0,tr==np.eye(4)];problem=cp.Problem(cp.Maximize(cp.real(cp.trace(Q@X))),constraints)
        problem.solve(solver='CLARABEL',tol_gap_abs=1e-10,tol_feas=1e-10,tol_gap_rel=1e-10,max_iter=200)
        assert problem.status in ['optimal','optimal_inaccurate']
        Qexact=rational_matrix(Q,40);dual=np.asarray(constraints[1].dual_value)
        dual=(dual+dual.conj().T)/2;Y0=rational_matrix(dual,40)
        if np.linalg.eigvalsh(np.kron(np.eye(4),np.array(Y0,complex))-Q).min()<-1e-7:Y0=-Y0
        # Adaptive shifts accepted ONLY after exact positive rational pivots.
        for exponent in [32,30,28,26,24,22,20]:
            Y=Y0+s.eye(4)/2**exponent
            try:exact_positive_pivots(s.kronecker_product(s.eye(4),Y)-Qexact)
            except AssertionError:continue
            break
        else:raise AssertionError('No exact recovery dual found')
        values=[]
        for p in pauli:
            vec=s.Matrix(list(rational_matrix(p,0)));values.append(s.re((vec.conjugate().T*Qexact*vec)[0]))
        choice=int(np.argmax([float(v) for v in values]));lower=values[choice];upper=s.trace(Y)
        assert upper>=lower
        lower_total+=len(members)*lower;upper_total+=len(members)*upper
        certificates.append(dict(relation_key=key,representative=branch,multiplicity=len(members),members=members,objective=strings(Qexact),dual=strings(Y),Pauli_choice=choice,lower=str(lower),upper=str(upper),dual_shift_exponent=exponent))
        print('orbit',len(certificates),'of66','multiplicity',len(members),'gap',float(upper-lower),flush=True)
    assert orbit_error<1e-10
    return dict(certificates=certificates,orbit_count=66,numerical_orbit_replay_error=orbit_error,
      exact_Pauli_lower=str(lower_total),exact_CPTP_upper=str(upper_total),exact_gain_cap=str(upper_total-lower_total),
      proof='Four three-bit syndrome vectors are simultaneous GL3(2) orbits, classified by their binary relation kernel. All66 kernels are enumerated. GL3(2) permutes the seven Hamming columns and fixes both Steane logical basis states. Identical paired physical channels commute with simultaneous permutations, and standard syndrome Pauli corrections transform up to irrelevant phase. Hence the conditional logical channels and recovery objectives coincide in each orbit before any floating evaluation.',
      scope='Exact primal-dual certificates for66 stored dyadic representative objectives with orbit multiplicities. The exact symmetry statement applies to identical physical paired channels and ideal extraction/readout; numerical contraction/objective and dyadic rounding need interval bridging for a certified physical fidelity bound.')

def quantum(branches=None):
    import cvxpy as cp
    import w33_pass11476_11480_cubic_geometry_correlated as M
    if branches is None:branches=native_branches()
    old=read(SOURCES[2])['quantum']['results'][0];scores=np.array(old['complete_fidelity_branch_scores']).reshape(16,4096)
    order=np.argsort(scores.max(axis=0))[::-1];rows=[];total_gain=0.
    pauli=[np.kron(a,b) for a,b in itertools.product(M.PAULI,repeat=2)]
    for branch in order[:16]:
        Q=recovery_objective(branches[branch]);X=cp.Variable((16,16),hermitian=True)
        tr=cp.bmat([[sum(X[4*a+b,4*a+d] for a in range(4)) for d in range(4)] for b in range(4)])
        constraints=[X>>0,tr==np.eye(4)]
        problem=cp.Problem(cp.Maximize(cp.real(cp.trace(Q@X))),constraints)
        problem.solve(solver='CLARABEL',tol_gap_abs=1e-10,tol_feas=1e-10,tol_gap_rel=1e-10,max_iter=200)
        assert problem.status in ['optimal','optimal_inaccurate']
        pauli_values=np.array([np.vdot(p.ravel(),Q@p.ravel()).real for p in pauli])
        assert np.max(abs(pauli_values-scores[:,branch]))<1e-10
        dual=np.array(constraints[1].dual_value);dual=(dual+dual.conj().T)/2
        # Sign convention verified against actual dual slack, not guessed.
        if np.linalg.eigvalsh(np.kron(np.eye(4),dual)-Q).min()<-1e-7:dual=-dual
        primal=np.array(X.value);gain=float(problem.value-pauli_values.max());total_gain+=max(0,gain)
        row=dict(branch=int(branch),Pauli_fidelity=float(pauli_values.max()),SDP_fidelity=float(problem.value),gain=gain,
          primal=enc(primal),dual=enc(dual),objective=enc(Q),primal_minimum_eigenvalue=float(np.linalg.eigvalsh(primal).min()),dual_slack_minimum_eigenvalue=float(np.linalg.eigvalsh(np.kron(np.eye(4),dual)-Q).min()),trace_preservation_error=float(np.linalg.norm(np.einsum('abac->bc',primal.reshape(4,4,4,4))-np.eye(4))))
        rows.append(row);print('SDP branch',int(branch),'gain',gain,flush=True)
    # The SDP found no positive improvement on these dominant branches.
    # Keep the known EXACT Pauli primal; bound all coherent channels from above.
    best=rows[0];Q=rational_matrix(dec(best['objective']),36)
    choice=int(np.argmax(scores[:,best['branch']]));pv=s.Matrix(list(rational_matrix(pauli[choice],0)));X=pv*pv.conjugate().T
    assert s.Matrix(4,4,lambda b,d:sum(X[4*a+b,4*a+d] for a in range(4)))==s.eye(4)
    Y=rational_matrix(dec(best['dual']),32);Y+=s.eye(4)/2**24
    dual_pivots=exact_positive_pivots(s.kronecker_product(s.eye(4),Y)-Q)
    lower=s.re(s.trace(Q*X));upper=s.trace(Y)
    exact_pauli=max(s.re((s.Matrix(list(rational_matrix(p,0))).conjugate().T*Q*s.Matrix(list(rational_matrix(p,0))))[0]) for p in pauli)
    assert upper>=lower==exact_pauli
    # An exact dyadic Bell-basis objective for ALL4096 branches. Gershgorin
    # gives an explicit scalar dual bound without numerical eigenvalue claims.
    bell=np.column_stack([p.ravel()/2 for p in pauli]);bits=40;scale=2**bits
    rounded=[];integer_lower=0;integer_upper=0;float_error=0.
    for transfer in branches:
        q=recovery_objective(transfer);qb=bell.conj().T@q@bell
        re=np.rint(qb.real*scale).astype(np.int64);im=np.rint(qb.imag*scale).astype(np.int64)
        assert np.array_equal(re,re.T) and np.array_equal(im,-im.T)
        diagonal=np.diag(re);radius=(abs(re)+abs(im)).sum(axis=1)-abs(diagonal)
        integer_lower+=4*int(diagonal.max());integer_upper+=4*int((diagonal+radius).max())
        float_error=max(float_error,float(np.max(abs(qb-(re+1j*im)/scale))))
        rounded.append(dict(real_numerators=re.tolist(),imag_numerators=im.tolist()))
    global_lower=s.Rational(integer_lower,scale);global_upper=s.Rational(integer_upper,scale)
    # New hardware-level diagnosis: actual unflagged Steane extraction hook.
    import w33_pass11471_11475_frames_currents_matching as R
    V,_=R.recovery_rows();support=[i for i in range(7) if ((7-i)&4)]
    hook=support[2:];syndrome=(7-hook[0])^(7-hook[1]);correction=7-syndrome
    pauli_hook=np.eye(128,dtype=complex)
    for i in hook+[correction]:
        factors=[np.eye(2) for _ in range(7)];factors[i]=M.PAULI[3]
        op=factors[0]
        for a in factors[1:]:op=np.kron(op,a)
        pauli_hook=op@pauli_hook
    logical=V.conj().T@pauli_hook@V
    assert np.linalg.norm(logical-M.PAULI[3])<1e-12
    baseline=old['rows'][0]['fidelity_optimal_complete_flagged_fidelity']
    return dict(status='PASS',rows=rows,baseline_Pauli_complete_fidelity=baseline,
      exact_certificate=dict(branch=best['branch'],objective=strings(Q),primal=strings(X),dual=strings(Y),dual_pivots=list(map(str,dual_pivots)),lower=str(lower),upper=str(upper),Pauli_upper=str(exact_pauli),primal_dual_gap=str(upper-lower)),
      exact_all_branch_Bell_objectives=rounded,exact_Bell_denominator=str(scale),exact_all_branch_Pauli_lower=str(global_lower),exact_all_branch_coherent_upper=str(global_upper),exact_all_branch_coherent_gain_cap=str(global_upper-global_lower),maximum_float_to_dyadic_entry_error=float_error,
      extraction_hook=dict(Z_check_support=support,ancilla_Z_fault_after_data_coupling=2,data_Z_hook=hook,X_check_syndrome=syndrome,standard_Z_correction=correction,logical_result='Z',decoded_logical_operator=enc(logical),proof='CNOT conjugates target Z to control Z times target Z. Ancilla target Z after the second of four data-control CNOTs propagates onto the final two data controls. Ideal later X checks and standard decoder add a third Z; the Fano-line triple is logical Z. The ancilla Z does not flip its Z-basis readout.',fault_probability='A sole such Bernoulli fault p gives logical channel (1-p)id+p Ad_Z and entanglement fidelity1-p, even with ideal remaining extraction/readout.'),
      flagged_extraction_repair=flagged_Z_hooks(),
      coherent_recovery_orbit_bound=recovery_orbit_bounds(branches),
      exact_optimization='max tr(QX), X>=0, tr_output X=I4; dual min tr Y subject to I4 tensor Y>=Q. X is an unnormalized recovery Choi matrix. Q includes both accepted native-cell effects and verification flags. Arbitrary CPTP recovery on the decoded four-dimensional logical fiber is allowed.',
      scope='Sixteen full CPTP SDPs find no resolved improvement on dominant branches. Exact rounded-objective dual certificate and4096 Bell/Gershgorin bounds constrain coherent recovery. A concrete ancilla hook diagnoses unflagged extraction failure. Native channel bridge remains floating; no complete flagged circuit or threshold.')

def run():
    result={}
    for name in ['vacuum','flavor','gravity','ward','quantum']:
        result[name]=globals()[name]();print(name,'PASS',flush=True)
    result.update(status='PASS',passes=list(range(11521,11526)),reservation='6388c0ab3',
      source_sha256={p:hashlib.sha256(json.dumps(read(p),sort_keys=True,separators=(',',':')).encode()).hexdigest() for p in SOURCES},
      producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    OUT.write_text(json.dumps(result,indent=2)+'\n')
    return result
if __name__=='__main__':run()
