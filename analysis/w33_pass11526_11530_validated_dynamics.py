"""Five physical-frontier investigations, with executable interval witnesses.

Owners:11521 quotient potential;11432 flagged protocol;11482 DeWitt history;
11438 native action;11484/11524 local-measure boundary. See companion report.
"""
import hashlib, itertools, json, sys
from pathlib import Path
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/w33_pass11526_11530_validated_dynamics.json'
read=lambda name:json.loads((ROOT/name).read_text())

class Jet:
    """Second real-coordinate jet with rigorous complex Arb ball coefficients."""
    n=20
    def __init__(self,v,g=None,h=None):
        from flint import acb
        self.v=acb(v)
        self.g=np.full(self.n,acb(0),object) if g is None else g
        self.h=np.full((self.n,self.n),acb(0),object) if h is None else h
    @classmethod
    def variable(cls,v,i):
        out=cls(v);out.g[i]=out.v*0+1;return out
    @staticmethod
    def coerce(other):return other if isinstance(other,Jet) else Jet(other)
    def __add__(self,other):
        b=self.coerce(other);return Jet(self.v+b.v,self.g+b.g,self.h+b.h)
    __radd__=__add__
    def __neg__(self):return Jet(-self.v,-self.g,-self.h)
    def __sub__(self,other):return self+-self.coerce(other)
    def __rsub__(self,other):return self.coerce(other)+-self
    def __mul__(self,other):
        b=self.coerce(other)
        return Jet(self.v*b.v,self.g*b.v+b.g*self.v,self.h*b.v+b.h*self.v+np.outer(self.g,b.g)+np.outer(b.g,self.g))
    __rmul__=__mul__
    def unary(self,v,first,second):
        return Jet(v,first*self.g,first*self.h+second*np.outer(self.g,self.g))
    def reciprocal(self):return self.unary(1/self.v,-1/self.v**2,2/self.v**3)
    def __truediv__(self,other):return self*self.coerce(other).reciprocal()
    def __rtruediv__(self,other):return self.coerce(other)*self.reciprocal()
    def __pow__(self,n):
        if n<0:return self.reciprocal()**(-n)
        if n==0:return Jet(1)
        if n==1:return self
        a=self**(n//2);return a*a*(self if n%2 else 1)
    def conj(self):return Jet(self.v.conjugate(),np.array([x.conjugate() for x in self.g],object),np.array([[x.conjugate() for x in row] for row in self.h],object))
    def real(self):return Jet(self.v.real,np.array([x.real for x in self.g],object),np.array([[x.real for x in row] for row in self.h],object))
    def imag(self):return Jet(self.v.imag,np.array([x.imag for x in self.g],object),np.array([[x.imag for x in row] for row in self.h],object))
    def exp(self):
        a=self.v.exp();return self.unary(a,a,a)

def jet_potential(coords):
    """Same declared potential as11521, no floating constants or differencing."""
    z=[Jet.variable(v,i) for i,v in enumerate(coords)];zero=Jet(0)
    A=[[z[i]*(1j*z[3]).exp() if i==j else zero for j in range(3)] for i in range(3)]
    imaginary=[i for i in range(9) if i not in (1,5)]
    B=[[z[4+3*i+j]+(1j*z[13+imaginary.index(3*i+j)] if 3*i+j in imaginary else 0) for j in range(3)] for i in range(3)]
    transpose=lambda X:list(map(list,zip(*X)))
    conj=lambda X:[[v.conj() for v in row] for row in X]
    dagger=lambda X:transpose(conj(X))
    mul=lambda X,Y:[[sum((X[i][k]*Y[k][j] for k in range(3)),zero) for j in range(3)] for i in range(3)]
    trace=lambda X:sum((X[i][i] for i in range(3)),zero)
    norm=lambda X:sum((v*v.conj() for row in X for v in row),zero).real()
    det=lambda X:sum(((-1 if sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))%2 else 1)*X[0][p[0]]*X[1][p[1]]*X[2][p[2]] for p in itertools.permutations(range(3))),zero)
    def one(X):
        G=mul(dagger(X),X);t=trace(G)/3;G=[[G[i][j]-(t if i==j else 0) for j in range(3)] for i in range(3)]
        u=27*det(X)**2;r=(u*u.conj()).real()
        return norm(G)/2+r**2-r+u.real()**2+u.real()/2+r**3
    radius=norm(A)+norm(B);den=radius/2+Jet(1)/100
    C=[[v/den for v in row] for row in mul(transpose(A),conj(B))];Cd=dagger(C)
    U=[[2*(i==j)+(C[i][j]+Cd[i][j])/2 for j in range(3)] for i in range(3)]
    CC=mul(C,Cd);V=[[int(i==j)+CC[i][j] for j in range(3)] for i in range(3)]
    UV,VU=mul(U,V),mul(V,U);W=[[UV[i][j]-VU[i][j] for j in range(3)] for i in range(3)]
    chi=(27*det(A)**2).imag();q=trace(mul(mul(W,W),W)).imag()
    return (one(A)+one(B)-10000*chi/(1+chi**2)*q+radius**2/1000000).real(),chi,q

def bound_string(x):return x.str(35)

def validated_vacuum():
    from flint import arb,ctx
    ctx.prec=256
    center=read('data/w33_pass11521_11525_exact_frontiers.json')['vacuum']['quotient_control']['coordinates']
    # Exact decimal rationals are enclosed by Arb; every operation rounds outwards.
    radius='1e-25';c=[arb(v) for v in center];box=[arb(v,radius) for v in center]
    f0,_,_=jet_potential(c);f,chi,q=jet_potential(box)
    H0=np.array([[float(x.real.mid()) for x in row] for row in f0.h])
    pre=np.linalg.inv(H0);R=[[arb(format(x,'.18g')) for x in row] for row in pre]
    H=[[x.real for x in row] for row in f.h];g=[x.real for x in f0.g]
    defect=[[int(i==j)-sum((R[i][k]*H[k][j] for k in range(20)),arb(0)) for j in range(20)] for i in range(20)]
    errors=[abs(sum((R[i][k]*g[k] for k in range(20)),arb(0)))+arb(radius)*sum((abs(x) for x in defect[i]),arb(0)) for i in range(20)]
    contraction=max(sum((abs(x) for x in row),arb(0)) for row in defect)
    assert contraction<arb('1e-8') and all(e<arb(radius)/100 for e in errors)
    # Interval LDL is a rigorous positivity proof throughout the enclosing box.
    L=[[arb(int(i==j)) for j in range(20)] for i in range(20)];pivots=[]
    for i in range(20):
        p=H[i][i]-sum((L[i][k]**2*pivots[k] for k in range(i)),arb(0));assert p>0;pivots.append(p)
        for j in range(i+1,20):L[j][i]=(H[j][i]-sum((L[j][k]*L[i][k]*pivots[k] for k in range(i)),arb(0)))/p
    assert all(box[i]>0 for i in (0,1,2,5,9)) and not chi.v.real.contains(0) and not q.v.real.contains(0)
    return dict(status='PASS',precision_bits=256,center=center,radius=radius,preconditioner=[[format(x,'.18g') for x in row] for row in pre],
        Krawczyk_contraction=bound_string(contraction),inclusion_errors=[bound_string(x) for x in errors],interval_LDL_pivots=[bound_string(x) for x in pivots],energy=bound_string(f.v.real),CP_chi=bound_string(chi.v.real),CP_commutator=bound_string(q.v.real),
        theorem='Strict Krawczyk inclusion and contraction prove a unique exact20-coordinate stationary root in this box. Interval LDL proves positive quotient Hessian throughout it. Nonzero CP-odd invariants exclude its CP-conserving orbit. Simultaneous hollowisation then gives an exact full324-field stationary representative by11521 fixed-stratum transport.',
        boundary='This proves existence and quotient stability for the declared potential, not full transverse stability, a global minimum, observed masses or a selected physical Standard Model vacuum.')

def electroweak_action():
    """An explicit positive portal; coefficients/inventory remain supplied."""
    chi,t=s.symbols('chi t',real=True);k=s.symbols('k',positive=True)
    # p=Hu^T epsilon Hd, S=|Hu|²+|Hd|²; C=|Hu†Hd|².
    # |p|²=(S²-D²)/4-C fixes all minima of this sum of squares.
    target=k*(1+s.I*chi);norm=s.sqrt(1+chi**2)
    source='k(1+i chi)';total=2*k*norm
    # The Gram-rank obstruction is in the actual signed E6 tensor, not a toy.
    from w33_pass11384_11388_native_dynamics import native_tensors
    d=native_tensors()[1][1]
    sm=read('data/w33_pass11506_11510_five_physics_targets.json')['SM']
    hu,hd=sm['neutral_Hu_index'],sm['neutral_Hd_index'];n=0
    assert d[hu,n,hd]==d[hd,n,hu]==1
    vs=s.Matrix(s.symbols('s0:3'));vh=s.Matrix(s.symbols('h0:3'))
    same=vs*vh.T+vh*vs.T;assert s.expand(same.det())==0
    # Two existing scalar copies need a genuinely mixed covariant.
    a=s.Matrix([1,0,0]);b=s.Matrix([0,1,0]);c=s.Matrix([0,0,1]);e=s.Matrix([2,1,3])
    mixed=(a*e.T+e*a.T+c*b.T+b*c.T)/2
    second=(a*s.Matrix([0,1,2*s.I*chi]).T+s.Matrix([0,1,2*s.I*chi])*a.T+s.Matrix([1,2,0])*b.T+b*s.Matrix([1,2,0]).T)/2
    assert mixed.det()==s.Rational(1,4) and mixed.rank()==3 and s.factor(second.det())==2*chi**2
    up=mixed*mixed.conjugate().T;down=second*second.conjugate().T;comm=up*down-down*up
    cp=s.factor(s.im(s.trace(comm**3)));assert cp==-s.Rational(69,2)*chi*(chi**2-19)
    pu,pd=up.charpoly(),down.charpoly();du=s.discriminant(pu.as_expr(),pu.gen);dd=s.factor(s.discriminant(pd.as_expr(),pd.gen))
    # A CP-even portal locks p to a CP-odd native order parameter. It does not
    # shift the native vacuum after minimizing the added EW fields.
    xx=s.symbols('x0:8',real=True);tt,la,ka,et,nu=s.symbols('t lambda kappa eta nu',positive=True)
    U=s.Matrix([xx[0]+s.I*xx[4],xx[1]+s.I*xx[5]]);D=s.Matrix([xx[2]+s.I*xx[6],xx[3]+s.I*xx[7]])
    uu=(U.conjugate().T*U)[0];dd=(D.conjugate().T*D)[0];p=U[0]*D[1]-U[1]*D[0];overlap=(U.conjugate().T*D)[0]
    potential=s.expand(la*(uu+dd-2*tt)**2+ka*(uu-dd)**2+et*overlap*s.conjugate(overlap)+nu*(p-tt)*s.conjugate(p-tt))
    hessian=s.hessian(potential,xx).subs(dict(zip(xx,[0,s.sqrt(tt),-s.sqrt(tt),0,0,0,0,0])))
    eigen=hessian.eigenvals();assert eigen[0]==3 and sum(eigen.values())==8
    return dict(status='PASS',potential='lambda(S-2k sqrt(1+chi²))²+kappa D²+eta |Hu†Hd|²+nu |Hu^T epsilon Hd-k(1+i chi)|²; all five coefficients positive',
      exact_minimum=dict(S=str(total),D='0',charged_overlap='0',p=str(target),value='0'),
      proof='All summands are nonnegative. Choose Hu=(0,h), Hd=(-k(1+i chi)/h,0), |h|²=k sqrt(1+chi²). This saturates every summand. CP sends chi->-chi and both doublets to their conjugates, so the joint action has real CP-even coefficients.',
      exact_real_coordinate_Hessian_spectrum={str(e):int(m) for e,m in eigen.items()},Hessian_parameter='t=k sqrt(1+chi²)>0; a unitary rephasing of Hd maps the complex target to t without changing the Hessian eigenvalues',
      backreaction='The minimized EW contribution is identically zero for every native chi. Thus this explicitly engineered positive portal preserves the exact native stationary solution; freezing chi in an uncompensated bilinear would instead shift it.',
      single_copy_rank_bound='On the restricted color/charge-preserving S,Hu,Hd ansatz, a self-composite has Hhu=s*hd*^T+hd*s*^T and Hhd=s*hu*^T+hu*s*^T, hence rank<=2. The actual signed tensor coefficient is1. A second SM singlet at native1 cannot help Hhd in this ansatz.',
      mixed_two_copy_covariant='H^c_ij=(dbar^{cab}/2M)(Phi*_ai Psi*_bj+Phi*_aj Psi*_bi). Two scalar copies already occur in the native324 model; this is distinct from the11522 self-composite.',
      mixed_full_rank_witness=[[str(v) for v in row] for row in mixed.tolist()],mixed_determinant=str(mixed.det()),mixed_second_witness=[[str(v) for v in row] for row in second.tolist()],second_determinant=str(2*chi**2),Gram_discriminants=[str(du),str(dd)],CP_commutator_cubic=str(cp),
      source_witness_scope='The integer family vectors and chi-dependent EW source are supplied ansatz data. They demonstrate two nondegenerate full-rank family matrices and CP transfer for chi!=0,chi²!=19, not dynamically selected CKM parameters. In particular the11526 native CP phase fails the SM symmetry screen; it cannot be directly identified with this hypothetical SM-compatible CP source.',
      healthy_auxiliary_completion='Define the unscaled dimension-two bilinear Jmix=M Hmix. Then M²||X-Jmix/M||², X in(27,bar6), selects X=Jmix/M. Expanding gives a positive mass, cubic X-Phi-Psi coupling and quartic counterterm; integrating X yields the allowed dimension-five identical-Weyl Yukawa. This adds a heavy sextet and its matching coefficients.',
      boundary='Exact EW minimum and rank escape in a supplied effective action/ansatz, not a full E6+EW gauge-compatible selected vacuum, measured electroweak scale, complete family-vector dynamics or observed masses/mixing. General 2HDM minimization and completing squares are established methods;11326 already owns an engineered noncommuting family vacuum.')

def symmetry_phase_audit():
    """Test physical gauge symmetry before interpreting any stationary root."""
    from w33_pass11384_11388_native_dynamics import native_tensors
    raw,(_,d,_)=native_tensors();unique={}
    for b in raw:
        for kind,a in [(0,b+b.T),(1,b-b.T)]:
            nz=np.flatnonzero(a)
            if not len(nz):continue
            a=a//int(np.gcd.reduce(abs(a.ravel()[nz])))
            if a.ravel()[nz[0]]<0:a=-a
            unique[kind,tuple(a.ravel())]=a*(1j if kind else 1)
    he=np.array(list(unique.values()));assert len(he)==78
    G=np.einsum('aij,bji->ab',he,he).real.astype(int)
    cartan=[i for i,a in enumerate(he) if np.array_equal(a,np.diag(np.diag(a)))];assert len(cartan)==6
    weights=s.Matrix([[int(he[j,i,i].real) for j in cartan] for i in range(27)])
    metric=s.Matrix(G[np.ix_(cartan,cartan)]).inv();weight_gram=weights*metric*weights.T
    assert set(weight_gram.diagonal())=={s.Rational(2,9)}
    # Two singlets of the canonical SM embedding cannot support these cubic
    # contractions: even d(P,P,-) is exactly zero on this subspace.
    assert not np.any(d[:,[0,1]][:,:,[0,1]])
    frame=[0,17,26];outside=[i for i in range(27) if i not in frame]
    blocks=he[:,outside,:][:,:,frame]
    R=s.Matrix(np.concatenate((blocks.real.reshape(78,-1),blocks.imag.reshape(78,-1)),axis=1).astype(int).T)
    kernel=R.nullspace();assert len(kernel)==30
    inside=s.Matrix([[int(a[i,i].real) for a in he] for i in frame])
    insidei=s.Matrix([[int(a[i,i].imag) for a in he] for i in frame])
    assert all(sum(a[i,i] for i in frame)==0 for a in he)
    # Cartan restrictions of Hermitian generators are real; root skew terms
    # have no diagonal. The frame normalizer acts by a two-dimensional torus.
    assert insidei==s.zeros(3,78)
    image=inside*s.Matrix.hstack(*kernel);assert image.rank()==2
    root=read('data/w33_pass11521_11525_exact_frontiers.json')['vacuum']['quotient_control']['coordinates']
    assert float(root[0])>float(root[1])>float(root[2])>0 and float(root[5])>0 and float(root[9])>0
    return dict(status='PASS',native_weights=[[str(x) for x in row] for row in weights.tolist()],Cartan_metric_inverse=[[str(x) for x in row] for row in metric.tolist()],weight_norm_squared='2/9',family_maximum='2/3',joint_maximum='8/9',
      frame_outside_constraints=[[int(x) for x in row] for row in R.tolist()],frame_normalizer_kernel=[[str(x) for x in v] for v in kernel],frame_normalizer_dimension=30,frame_torus_image_dimension=2,
      no_SM_singlet_stationary_point='The canonical SM singlets are native0,1 (owner11271). d(P,P,-)=0 there, so all native I6,I12,I18 and their first derivatives vanish. The CP coupling also vanishes with its first derivative. V restricted to both singlet-supported fields is one-half squared moments plus rho(||Phi||²+||Psi||²)², rho>0. Euler gives x.gradV=4V>0 for every nonzero pair. Thus the unchanged11438 action has no nonzero stationary point preserving this fixed SM embedding.',
      CP_root_stabilizer='Since Phi has full frame rank, any infinitesimal combined E6 x SU3 stabilizer must preserve that frame. Its compact frame normalizer is Spin8 plus a two-dimensional diagonal torus. For any common-left hollowising U, anti-Hermiticity of the compensating family generator forces Udag D U to commute with AA†. The three distinct certified singular values make it diagonal. Its commutator with B then vanishes only when all three diagonal entries agree, because B01 and B12 are nonzero. Trace zero forces D=0. Thus the combined stabilizer Lie algebra is Spin8.',
      Spin8_SM_obstruction='An injective su3 in so8 acts on real8 either as realification(3) plus1 plus1 or as the real adjoint8. The orthogonal centralizer is u1+so2 in the first case and zero in the second. Neither contains commuting su2. Therefore Spin8 cannot contain su3_c+su2_L, regardless of the SM embedding. The certified CP root is an exact stationary phase of the supplied model, not an SM-preserving vacuum; adding further vevs without changing this pair can only shrink its stabilizer.',
      alternate_coherent_selector='For one(27,3) field Phi let R=||Phi||² and mu be the same native E6/family moments. The exact bound ||mu||²<=8R²/9 follows by rotating the E6 moment into the Cartan and using the weight polytope with all27 weights of squared norm2/9, plus Tr(rho_F²)-R²/3<=2R²/3. Equality holds on compact coherent highest-weight products. Vnew=(R-1/2)²+kappa(8R²/9-||mu||²)+M²||Psi||², kappa,M²>0, has global zero-energy minima at R=1/2 on that coherent orbit, Psi=0. e0 tensor e1 is a representative; its E6 Spin10 stabilizer contains the standard SM chain. The potential is a supplied changed action, not a consequence of11438.',
      correct_next_phase='11271 already constructs e0,e1 plus a hypercharge adjoint with exactly the SM stabilizer. The coherent selector dynamically supplies a compatible first stage but does not select that full two-vector/adjoint orbit, measured couplings or CKM. CP/flavor dynamics must be rebuilt in an SM-compatible phase rather than transplanted from the Spin8 root.',
      boundary='Exact action-level and representation obstructions plus a globally minimized alternate supplied potential. Coherent-state moment maximization is classical. No change/retraction of11526 mathematical stationary existence, and no claim that Vnew is UV-selected, supersymmetric D-flat or a complete Standard Model.')

def coherent_parent_normal_certificate():
    """Exact full normal spectrum of the separately declared compatible phase."""
    from w33_pass11384_11388_native_dynamics import native_tensors
    raw=native_tensors()[0];unique={}
    for b in raw:
        for kind,a in [(0,b+b.T),(1,b-b.T)]:
            nz=np.flatnonzero(a)
            if not len(nz):continue
            a=a//int(np.gcd.reduce(abs(a.ravel()[nz])))
            if a.ravel()[nz[0]]<0:a=-a
            unique[kind,tuple(a.ravel())]=a*(1j if kind else 1)
    he=np.array(list(unique.values()));GE=s.Matrix(np.einsum('aij,bji->ab',he,he).real.astype(int));IE=GE.inv()
    family=[]
    for i,j in itertools.combinations(range(3),2):
        a=np.zeros((3,3),complex);a[i,j]=a[j,i]=1;family.append(a)
        a=np.zeros((3,3),complex);a[i,j]=-1j;a[j,i]=1j;family.append(a)
    family.extend([np.diag([1,-1,0]),np.diag([1,1,-2])]);family=np.array(family)
    GF=s.Matrix(np.einsum('aij,bji->ab',family,family).real.astype(int));IF=GF.inv();C=s.diag(IE,IF)
    cols=[np.kron(a[:,0],np.array([1,0,0])) for a in he]+[np.kron(np.array([1]+[0]*26),a[:,0]) for a in family]
    columns=np.array(cols).T;V=s.Matrix(np.concatenate([columns.real,columns.imag],axis=0).astype(int));purity_gradient=V*C*V.T
    # The moment at e0 is Cartan. Build the entire real N=sum mu_a M_a.
    diagE=s.Matrix([[int(a[i,i].real) for a in he] for i in range(27)]);diagF=s.Matrix([[int(a[i,i].real) for a in family] for i in range(3)])
    pairing=s.kronecker_product(diagE*IE*diagE.T,s.ones(3))+s.kronecker_product(s.ones(27),diagF*IF*diagF.T)
    native=s.diag(*list(pairing[:,0]),*list(pairing[:,0]));E=s.zeros(162);E[0,0]=1
    H=4*E+s.Rational(8,9)*(4*E+2*s.eye(162))-4*purity_gradient-2*native
    assert H==s.diag(*H.diagonal())
    from collections import Counter
    spectrum=Counter(H.diagonal());assert spectrum=={s.Integer(0):37,s.Rational(2,3):20,s.Rational(7,3):64,s.Rational(8,3):40,s.Integer(4):1}
    # All zero directions are actual combined compact gauge-orbit directions.
    tangents=s.Matrix(np.concatenate([-columns.imag,columns.real],axis=0).astype(int));assert tangents.rank()==37 and H*tangents==s.zeros(162,86)
    spectrum[s.Integer(2)]+=162
    return dict(status='PASS',changed_action='Vnew=(||Phi||²-1/2)²+(8||Phi||^4/9-||mu(Phi)||²)+||Psi||²',
      exact_stationary_representative='Phi=e0_native tensor e0_family /sqrt2, Psi=0',
      exact_full324_real_Hessian_spectrum={str(v):int(m) for v,m in spectrum.items()},Phi_Hessian_diagonal=list(map(str,H.diagonal())),
      purity_gradient_matrix=[[str(v) for v in row] for row in purity_gradient.tolist()],moment_action_diagonal=list(map(str,native.diagonal())),gauge_tangents=[[int(v) for v in row] for row in tangents.tolist()],gauge_orbit_rank=37,
      exact_stability='All287 normal real-coordinate Hessian eigenvalues are positive. The37 zeros are exactly the compact E6 x SU3 family orbit. Combined with the weight-polytope inequality, this is a global-minimum and full normal/Morse-Bott certificate for the changed supplied action, unlike the transverse-open CP root of11438.',
      boundary='This selected coherent phase preserves a Spin10-compatible E6 subgroup and a family SU2; it is not the final SM vacuum and has no CP-breaking invariant. The additive constant makes its potential minimum zero by construction, not a cosmological-constant solution. Supplied couplings, kinetic normalization and scale remain inputs. Additional SM/family breaking and CP dynamics require a separate compatible action.')

def coherent_two_stage_normal_certificate():
    """A second supplied action, with its own full normal certificate."""
    from w33_pass11384_11388_native_dynamics import native_tensors
    raw=native_tensors()[0];unique={}
    for b in raw:
     for kind,a in [(0,b+b.T),(1,b-b.T)]:
      nz=np.flatnonzero(a)
      if not len(nz):continue
      a=a//int(np.gcd.reduce(abs(a.ravel()[nz])))
      if a.ravel()[nz[0]]<0:a=-a
      unique[kind,tuple(a.ravel())]=a*(1j if kind else 1)
    he=np.array(list(unique.values()));IE=s.Matrix(np.einsum('aij,bji->ab',he,he).real.astype(int)).inv()
    family=[]
    for i,j in itertools.combinations(range(3),2):
     a=np.zeros((3,3),complex);a[i,j]=a[j,i]=1;family.append(a)
     a=np.zeros((3,3),complex);a[i,j]=-1j;a[j,i]=1j;family.append(a)
    family.extend([np.diag([1,-1,0]),np.diag([1,1,-2])]);family=np.array(family)
    IF=s.Matrix(np.einsum('aij,bji->ab',family,family).real.astype(int)).inv();C=s.diag(IE,IF)
    DE=s.Matrix([[int(a[i,i].real) for a in he] for i in range(27)]);DF=s.Matrix([[int(a[i,i].real) for a in family] for i in range(3)])
    pair=s.kronecker_product(DE*IE*DE.T,s.ones(3))+s.kronecker_product(s.ones(27),DF*IF*DF.T)
    Hs=[];tangents=[]
    for point in [0,1]:
     cols=np.array([np.kron(a[:,point],np.array([1,0,0])) for a in he]+[np.kron(np.eye(27)[point],a[:,0]) for a in family]).T
     V=s.Matrix(np.concatenate([cols.real,cols.imag],axis=0).astype(int));grad=V*C*V.T
     E=s.zeros(162);E[3*point,3*point]=1
     H=4*E+s.Rational(8,9)*(4*E+2*s.eye(162))-4*grad-2*s.diag(*list(pair[:,3*point]),*list(pair[:,3*point]))
     assert H==s.diag(*H.diagonal());Hs.append(H)
     tangents.append(s.Matrix(np.concatenate([-cols.imag,cols.real],axis=0).astype(int)))
    # DF_Phi for F=(KPhi-RPhi/18)Psi at Phi=e0/sqrt2,Psi=e1/sqrt2.
    inputE=np.array([np.kron(a[:,0],np.array([1,0,0])) for a in he]).T
    V=s.Matrix(np.concatenate([inputE.real,inputE.imag],axis=0).astype(int))
    outputE=np.array([np.kron(a[:,1],np.array([1,0,0])) for a in he]).T
    W=s.Matrix(outputE.real.astype(int))+s.I*s.Matrix(outputE.imag.astype(int))
    JP=W*IE*V.T;JP[3,0]-=s.Rational(1,18)
    M=DE*IE*DE[0,:].T
    JQ=s.Rational(1,2)*s.kronecker_product(s.diag(*list(M-s.ones(27,1)/18)),s.eye(3))
    JQ=JQ.row_join(s.I*JQ);J=JP.row_join(JQ)
    Jreal=J.applyfunc(s.re);Jimag=J.applyfunc(s.im);Jall=Jreal.col_join(Jimag)
    H=s.diag(*Hs)+2*Jall.T*Jall
    # Positive family density alignment Rphi Rpsi-tr(rhoFphi rhoFpsi).
    # Relative-family part of the quadratic form is (1/2) sum_f |deltaPhi_0f-deltaPsi_1f|².
    # All other off-family native components also have diagonal curvature1.
    A=s.zeros(4,324)
    for k,(f,imag) in enumerate(itertools.product([1,2],[False,True])):
     A[k,f+81*imag]=1;A[k,162+3+f+81*imag]=-1
    H+=s.diag(*([int(i%3!=0) for i in range(81)]*4))
    for row in range(4):
     ids=[j for j in range(324) if A[row,j]]
     H[ids[0],ids[1]]-=1;H[ids[1],ids[0]]-=1
    
    # Exact connected blocks provide positivity and complete eigenvalues.
    adj=[set() for _ in range(324)]
    for i in range(324):
     for j in range(i):
      if H[i,j]:adj[i].add(j);adj[j].add(i)
    seen=set();sizes=[];spectrum={};certblocks=[]
    for i in range(324):
     if i in seen:continue
     todo=[i];block=[];seen.add(i)
     while todo:
      a=todo.pop();block.append(a)
      for b in adj[a]-seen:seen.add(b);todo.append(b)
     sizes.append(len(block));small=H.extract(block,block);eig=small.eigenvals()
     certblocks.append(dict(indices=block,matrix=[[str(x) for x in row] for row in small.tolist()]))
     for e,m in eig.items():spectrum[e]=spectrum.get(e,0)+m
    T=s.Matrix.vstack(*tangents)
    assert T.rank()==58 and spectrum.get(0)==58 and H*T==s.zeros(324,86)
    assert sum(spectrum.values())==324 and all(e>=0 for e in spectrum)
    return dict(status='PASS',action='V2=(Rphi-1/2)^2+Dphi+(Rpsi-1/2)^2+Dpsi+||(Kphi-Rphi I/18)Psi||^2+Rphi Rpsi-tr(rhoFphi rhoFpsi); D=8R^2/9-||mu||^2',
     moment_map='Kphi=sum_ab (GE^-1)_ab tr(Phi^dag Ta Phi) Tb, E6 generators only',
     reference='Phi=e0_native tensor e0_family/sqrt2; Psi=e1_native tensor e0_family/sqrt2',
     full324_spectrum={str(e):m for e,m in spectrum.items()},single_field_Hessian_diagonals=[list(map(str,h.diagonal())) for h in Hs],selector_Jacobian=[[str(x) for x in row] for row in Jall.tolist()],gauge_tangents=[[int(x) for x in row] for row in T.tolist()],Hessian_blocks=certblocks,
     family_alignment_Hessian='Diagonal1 on all real/imaginary components with family index1 or2 in both fields; offdiagonal-1 between Phi(native0,f) and Psi(native1,f) in the matching real/imaginary components. Every other entry0.',
     proof='The coherent moment defects are globally nonnegative. PSD family densities obey tr(rhoPhi rhoPsi)<=tr(rhoPhi)tr(rhoPsi); equality aligns their rank-one family factors. At the displayed coherent pair the moment selector vanishes because native1 has E6 moment weight1/18. Thus V2>=0 and this pair is an exact global minimum. Rational connected Hessian blocks have58 zeros, all exactly the combined compact gauge orbit, and266 positive eigenvalues.',
     boundary='A separately supplied degree-six EFT action, not the first Vnew or unchanged11438. Its e0/e1 E6 stabilizer is24-dimensional and preserves the canonical SM gauge algebra (owner11271); additional SU5-to-SM breaking, physical CP/flavor/coefficients and UV completion remain open. Minimum0 is engineered, not a cosmological-constant derivation. Global orbit uniqueness is not claimed.')

def path_frames():
    """Exact first-jet directions on the actual harmonic periodic cover."""
    src=read('data/w33_pass11493_11497_exact_matching_inference.json')['metric'];edges=src['edges'];h=np.array(src['displacements_times80'],int)
    adj=[[] for _ in range(80)]
    for e,((a,b),v) in enumerate(zip(edges,h)):
        adj[a].append((b,v,e,1));adj[b].append((a,-v,e,-1))
    records=[];rank_histories=[]
    for start in range(80):
        walks=[(start,np.zeros(3,int),[],[start])];acc=[];ranks=[]
        for radius in range(1,4):
            nxt=[]
            for end,v,path,vertices in walks:
                for j,w,e,sgn in adj[end]:
                    if len(vertices)>1 and j==vertices[-2]:continue
                    nxt.append((j,v+w,path+[(e,sgn)],vertices+[j]))
            walks=nxt;acc+=walks;vectors=[v.tolist() for _,v,_,_ in acc if np.any(v)]
            ranks.append(s.Matrix(vectors).rank())
        rank_histories.append(ranks)
        G=sum((s.Matrix(v)*s.Matrix(v).T for _,v,_,_ in acc),s.zeros(3))
        assert G.det()>0
        inverse=G.inv();replay=sum((inverse*s.Matrix(v)*s.Matrix(v).T for _,v,_,_ in acc),s.zeros(3));assert replay==s.eye(3)
        records.append(dict(site=start,Gram=[[str(v) for v in row] for row in G.tolist()],determinant=str(G.det()),inverse=[[str(v) for v in row] for row in inverse.tolist()],walks=[dict(end=int(j),displacement=v.tolist(),edges=path) for j,v,path,_ in acc]))
    from collections import Counter
    census=[dict(Counter(r[k] for r in rank_histories)) for k in range(3)]
    return dict(status='PASS',rank_census_by_radius=census,site_rank_histories=rank_histories,sites=records,
      exact_map='For each site i, G_i=sum_p d_p d_p^T over nonbacktracking lifted walks of length<=3. nabla_i psi=sum_p (G_i^-1 d_p)(U_p psi_endpoint-psi_i). It reproduces every affine lifted first jet exactly at all80 sites.',
      local_spin_covariance='U_p is the ordered product of endpoint links and transforms S_i U_p S_endpoint^-1. Thus each reconstructed derivative transforms S_i; gamma_i transforms by conjugation. This is radius-three spatial locality with explicit transport, not a chosen globally trivial spin frame.',
      boundary='52 one-hop sites have rank1 in this supplied embedding. Three-hop first-jet completeness repairs that kinematic obstruction, but does not by itself build a Hermitian chiral spacetime Dirac, spin structure, Einstein constraints or continuum limit. Under another invertible spatial basis G_i transforms by congruence and all ranks/reconstruction identities remain unchanged.')

def coupled_spin_history():
    """Named site-spin operator and constrained homogeneous mean-field history."""
    from scipy.linalg import expm
    from scipy.integrate import solve_ivp
    src=read('data/w33_pass11493_11497_exact_matching_inference.json')['metric'];edges=src['edges'];h=np.array(src['displacements_times80'],float)
    # Alpha_i=PauliX tensor Pauli_i, beta=PauliZ tensor I.
    paulis=[np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]],complex),np.diag([1.,-1.])]
    gamma=[np.kron(paulis[0],v) for v in paulis];beta=np.kron(paulis[2],np.eye(2));B=np.kron(np.eye(80),beta)
    parent=list(range(80))
    def root(i):
        while parent[i]!=i:i=parent[i]
        return i
    chord=None
    for e,(a,b) in enumerate(edges):
        aa,bb=root(a),root(b)
        if aa==bb:chord=e;break
        parent[aa]=bb
    assert chord is not None
    link_generator=np.kron(np.eye(2),paulis[2])/2
    links=[expm(.2j*link_generator) if e==chord else np.eye(4) for e in range(len(edges))]
    lap=np.zeros((320,320),complex);kin=np.zeros((320,320),complex)
    for e,((a,b),v) in enumerate(zip(edges,h)):
        weight=6400/(v@v);g=sum(v[j]*gamma[j] for j in range(3));hop=-1j*80*g/(v@v)
        hop=-1j*40*(g@links[e]+links[e]@g)/(v@v)
        kin[4*a:4*a+4,4*b:4*b+4]+=hop;kin[4*b:4*b+4,4*a:4*a+4]+=hop.conj().T
        lap[4*a:4*a+4,4*a:4*a+4]+=weight*np.eye(4);lap[4*b:4*b+4,4*b:4*b+4]+=weight*np.eye(4)
        lap[4*a:4*a+4,4*b:4*b+4]-=weight*links[e];lap[4*b:4*b+4,4*a:4*a+4]-=weight*links[e].conj().T
    # This supplied Wilson term changes the inventory to320 site-spin states.
    K=kin+.01*B@lap;H=K+.1*B;assert np.linalg.norm(H-H.conj().T)<1e-10
    rng=np.random.default_rng(11528)
    rotations=[expm(1j*sum(float(v[j])*np.kron(np.eye(2),paulis[j])/2 for j in range(3))) for v in rng.normal(size=(80,3))]
    S=np.zeros((320,320),complex)
    for i,R in enumerate(rotations):S[4*i:4*i+4,4*i:4*i+4]=R
    # Rebuild transformed couplings rather than merely declaring H'=S H Sdag.
    Hg=np.zeros_like(H)
    for i in range(80):Hg[4*i:4*i+4,4*i:4*i+4]=H[4*i:4*i+4,4*i:4*i+4]
    for e,((a,b),v) in enumerate(zip(edges,h)):
        g=sum(v[j]*gamma[j] for j in range(3));ga=rotations[a]@g@rotations[a].conj().T;gb=rotations[b]@g@rotations[b].conj().T
        U=rotations[a]@links[e]@rotations[b].conj().T
        hop=-1j*40*(ga@U+U@gb)/(v@v)-.01*6400/(v@v)*beta@U
        Hg[4*a:4*a+4,4*b:4*b+4]+=hop;Hg[4*b:4*b+4,4*a:4*a+4]+=hop.conj().T
    cov=float(np.linalg.norm(Hg-S@H@S.conj().T));assert cov<1e-9
    # The inherited bosonic determinant uses the original scalar Laplacian,
    # distinct from the newly supplied connection Laplacian on spin fibers.
    scalar_lap=np.zeros((80,80))
    for (a,b),v in zip(edges,h):
        w=6400/(v@v);scalar_lap[a,a]+=w;scalar_lap[b,b]+=w;scalar_lap[a,b]-=w;scalar_lap[b,a]-=w
    eig=np.linalg.eigvalsh(scalar_lap);eig=np.maximum(eig,0);alpha=beta0=1.;gamma0=.001
    W=lambda a:np.sum(np.log((eig/a**2+4)/(eig/a**2+1)))
    Wp=lambda a:np.sum(6*eig/a**3/((eig/a**2+4)*(eig/a**2+1)))
    cosm=alpha-beta0/2-gamma0*Wp(1)/3
    V=lambda a:cosm*a**3+alpha/a**3+1.5*beta0*a+gamma0*W(a)
    Vp=lambda a:3*cosm*a**2-3*alpha/a**4+1.5*beta0+gamma0*Wp(a)
    psi=rng.normal(size=320)+1j*rng.normal(size=320)
    # Choose a positive beta-sector packet, so the declared expanding branch
    # has positive initial total potential; an arbitrary Dirac packet need not.
    psi.reshape(80,4)[:,2:]=0;psi/=np.linalg.norm(psi)
    Hspin=lambda a:K/a+.1*B
    energy=lambda a,p,psi:-p*p/(12*a)+V(a)+np.vdot(psi,Hspin(a)@psi).real
    initial_energy=V(1)+np.vdot(psi,Hspin(1)@psi).real;assert initial_energy>0
    p=-np.sqrt(12*initial_energy);start=np.r_[1.,p,psi.real,psi.imag]
    def rhs(time,y):
        a,p=y[:2];psi=y[2:322]+1j*y[322:];dpsi=-1j*(Hspin(a)@psi)
        force=-p*p/(12*a*a)-Vp(a)+np.vdot(psi,K@psi).real/(a*a)
        return np.r_[-p/(6*a),force,dpsi.real,dpsi.imag]
    sol=solve_ivp(rhs,[0,.0001],start,t_eval=np.linspace(0,.0001,5),rtol=1e-11,atol=1e-12);assert sol.success
    rows=[]
    for time,y in zip(sol.t,sol.y.T):
        a,p=y[:2];pp=y[2:322]+1j*y[322:]
        rows.append(dict(time=float(time),scale=float(a),scale_momentum=float(p),spinor_norm=float(np.vdot(pp,pp).real),Hamiltonian_constraint=float(energy(a,p,pp))))
    assert max(abs(r['Hamiltonian_constraint']) for r in rows)<1e-7 and max(abs(r['spinor_norm']-1) for r in rows)<1e-8
    return dict(status='PASS',site_spin_dimension=320,local_Spin_covariance_residual=cov,curved_chord=chord,fundamental_cycle_holonomy_trace=float(np.trace(links[chord]).real),initial_spinor=dict(real=psi.real.tolist(),imag=psi.imag.tolist()),initial_scale_momentum=float(p),history=rows,
      action='L=-3 a adot²/N-N V(a)+i psi† psidot-N psi† Hspin(a) psi; V contains the actual positive flux/winding/determinant inventory, Hspin(a)=K/a+m beta. This is the isotropic DeWitt reduction, not positive scalar kinetic substituted for gravity.',
      constraint='H=-p_a²/(12a)+V(a)+psi†Hspin(a)psi=0; adot=-p_a/(6a). The positive static energy is balanced by conformal momentum, giving a moving solution rather than a static-flat solution.',
      named_operator='K_ab=-i40[Gamma_a(h)U_ab+U_ab Gamma_b(h)]/||h||², K_ba=K_ab†, plus supplied Wilson r beta L_U; m=.1,r=.01. The first non-tree edge has Spin3 link exp(i .2 Sigma3/2), so a fundamental-cycle holonomy is nontrivial. One-hop principal incompleteness remains diagnosed by the exact path-frame audit.',
      boundary='Explicit numerical constrained homogeneous mean-field evolution with independent curved Spin3 links. The connection is supplied, not derived as a torsion-free Levi-Civita connection. The commuting complex spinor is a mean-field test, not quantized fermion statistics; coefficients are supplied. No local Einstein constraint algebra, derived spin structure, chiral SM projection or continuum Lorentz limit.11482 owns the general DeWitt dynamical mechanism.')

def nonlinear_anomaly_primitive():
    """Actual nonlinear gauge-invariant source-patch current, locality qualified."""
    import w33_pass11428_11432_native_alignment_measure_coarse_flags as Q
    import w33_pass11476_11480_cubic_geometry_correlated as M
    L=3;n=L**4;sites=list(itertools.product(range(L),repeat=4));ix={x:i for i,x in enumerate(sites)}
    B=M.gauge_boundary(L);A=np.zeros((n,4));A[0]=[.023,-.019,.017,.011];A[1,0]=-.009;A[4,2]=.007
    charges=[1,-4,2,-3,6];multiplicity=[6,3,3,2,1]
    def anomaly(A):
        out=np.zeros(n);gap=10.
        for q,m in zip(charges,multiplicity):
            P,_,_,g,_,_=Q.weak_overlap(q,1.,A,[],L=L);gap=min(gap,g)
            out+=m*q*np.trace(P.reshape(n,4,n,4),axis1=1,axis2=3).diagonal().real
        return out,gap
    # Rooted shortest coordinate paths produce a finite incidence right inverse.
    root=0;flows=np.zeros((n*4,n))
    for j,x in enumerate(sites):
        here=list(x)
        for mu in range(4):
            while here[mu]:
                before=ix[tuple(here)];step=-1 if here[mu]<=L/2 else 1;here[mu]=(here[mu]+step)%L;after=ix[tuple(here)]
                edge=next(e for e in range(n*4) if B[e,before]*B[e,after]<0)
                # The edge sign supplies +1 at the local departure vertex.
                flows[edge,j]=B[edge,before]
    assert np.array_equal(B.T@flows,np.eye(n)-np.eye(n)[:,[root]]@np.ones((1,n)))
    assert np.array_equal(np.abs(flows).sum(axis=0),[sum(min(v,L-v) for v in x) for x in sites])
    density,gap=anomaly(A);current=flows@density
    half,gaphalf=anomaly(A/2);nonlinear=np.linalg.norm(half-density/2)
    assert max(abs(density))>1e-9 and nonlinear>1e-9
    omega=np.linspace(-.003,.004,n);Ag=A.ravel()+B@omega;densityg,gapg=anomaly(Ag.reshape(n,4));currentg=flows@densityg
    ward=float(np.linalg.norm(B.T@current-density));cov=float(np.linalg.norm(current-currentg))
    assert abs(density.sum())<1e-10 and ward<1e-10 and cov<1e-10
    # Exact geometry, second construction: shortest-path routing gives a
    # source-centered exponential tail from a gapped finite-rank perturbation.
    return dict(status='PASS',L=L,background=A.tolist(),gauge_parameter=omega.tolist(),density=density.tolist(),current=current.tolist(),incidence=B.tolist(),routing=flows.tolist(),Ward_residual=ward,gauge_invariance_error=cov,minimum_Wilson_gap=min(gap,gapg,gaphalf),density_peak=float(max(abs(density))),nonlinear_half_amplitude_defect=float(nonlinear),
      nonlinear_formula='k(A)=T q(A), where T routes each density q(x) along a shortest path from the fixed source root. div(Tq)=q when sum q=0. q(A) is the actual charge-weighted nonlinear overlap density, not its linearization. No gauge-field reference is transformed or retained.',
      infinite_volume_source_patch_lemma='For a perturbation supported in a fixed bounded patch S with a uniformly gapped interpolation, resolvent locality gives |q(x)|<=C exp(-c dist(x,S)). Rooted shortest-path routing then bounds |k(e)| by a polynomial distance prefactor times the same exponential. Charge cancellation and unchanged index supply sum q=0. This proves a localized-source current tail, not the general functional locality condition for dense backgrounds.',
      retained_obstruction='T is geometry-global. The functional derivative T Dq for unrestricted dense backgrounds need not decay near its varied link; root choice is external. Gauge invariance and nonlinear Ward identity now hold, but spatial-current decay for bounded defects must not be promoted to a local chiral measure. A two-variable local Poincare construction and finite-volume toron/flux corrections remain required.',
      earlier_owners=['11509 nonlinear gauge-invariant finite Green primitive (w33_pass11506_11510_five_physics_targets.py)','11484 differentiated support flows','11479 Coulomb patch measure','11524 standard measure one-form'])

def coherent_flag_ports():
    """Lift the owned11432 fault census to a non-twirled port theorem."""
    import w33_pass11428_11432_native_alignment_measure_coarse_flags as Q
    decoder,single,cases=Q.flag_decoder();assert len(cases)==1714
    # A final ideal recovery is a reference decoding map, not another noisy
    # fault-tolerant hardware gadget. Enumerate the correctable residual classes.
    final={v&63:v for v in single};assert len(final)==22
    rows=[]
    for e,history,exact in cases:
        residual=e^decoder[history];assert residual in single
        assert residual^final[residual&63]==0
        rows.append(dict(error_class=int(e),history=list(history),correction=int(decoder[history]),residual_class=int(residual),terminal_correction=int(final[residual&63])))
    # All arbitrary two-qubit port Kraus operators lie in the 16-Pauli span.
    X=s.Matrix([[0,1],[1,0]]);Z=s.diag(1,-1)
    paulis=[s.kronecker_product(X**a*Z**b,X**c*Z**d) for a,b,c,d in itertools.product(range(2),repeat=4)]
    gram=s.Matrix([[s.trace(a.conjugate().T*b) for b in paulis] for a in paulis]);assert gram==4*s.eye(16)
    # Norm-bound the entire multi-fault remainder of a coherent port circuit.
    # At most132 gate slots; inactive conditional slots are controlled identities.
    envelopes=[]
    for delta in [s.Rational(1,10**5),s.Rational(1,10**4),s.Rational(1,10**3)]:
        bad=(1+delta)**132-1-132*delta;bound=bad**2
        envelopes.append(dict(delta=str(delta),bad_amplitude=str(bad),infidelity_bound=str(bound),float_bound=float(bound)))
    return dict(status='PASS',owner='11432: actual three flagged rounds and conditional clean followup',fault_cases=rows,case_count=len(rows),Pauli_operator_basis_rank=16,maximum_CNOT_slots=132,coherent_fault_envelopes=envelopes,
      instrument_span_theorem='After the fixed11432 history decoder and reference ideal terminal recovery, every Pauli basis fault at a given port induces a scalar logical identity in each fully recorded measurement branch. Linearity therefore corrects every CPTP channel at that single two-qubit port, including coherent superpositions and amplitude damping, without Pauli twirling. Ancilla preparation/readout basis faults have the same linear-span extension. Clean encoded input and ideal other operations are required.',
      coherent_many_port_bound='For unitary errors with ||F_j-I||<=delta at at most132 active CNOT ports, expand the measurement-dilated adaptive circuit into fault subsets. The zero/single-fault terms have scalar identity logical action after reference recovery. The remaining isometry error has norm<=b=(1+delta)^132-1-132delta. Projection orthogonal to the encoded Bell state annihilates all correctable terms; hence 1-F_e<=min(1,b²). This is O(delta^4) in infidelity, not O(delta^4) in diamond distance.',
      boundary='Logical-port fault theorem and rigorous coherent norm envelope for a declared ideal circuit. Actual native compiled gates have systematic errors/leakage at every port, so their measured residual norms and accepted carrier conditions must be matched before applying this bound. No arbitrary native-leakage correction, noisy terminal decoder, multi-fault threshold or new general QEC theorem is claimed.')

def produce():
    data={'status':'PASS','passes':'11526-11530','reservation':'057834b58'}
    for name,fn in [('vacuum',validated_vacuum),('symmetry_phase',symmetry_phase_audit),('coherent_parent',coherent_parent_normal_certificate),('coherent_two_stage',coherent_two_stage_normal_certificate),('electroweak',electroweak_action),('gravity_frames',path_frames),('spin_history',coupled_spin_history),('anomaly',nonlinear_anomaly_primitive),('quantum',coherent_flag_ports)]:
        data[name]=fn();print(name,'PASS',flush=True)
    sources=['data/w33_pass11521_11525_exact_frontiers.json','data/w33_pass11493_11497_exact_matching_inference.json','data/w33_pass11506_11510_five_physics_targets.json','data/w33_pass11428_11432_native_alignment_measure_coarse_flags.json']
    data['source_sha256']={p:hashlib.sha256(json.dumps(read(p),sort_keys=True,separators=(',',':')).encode()).hexdigest() for p in sources}
    data['producer_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    OUT.write_text(json.dumps(data,indent=2)+'\n');return data

if __name__=='__main__':produce()
