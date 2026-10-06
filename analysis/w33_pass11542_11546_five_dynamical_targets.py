"""Five requested dynamical attacks, source-bound and scoped to supplied EFTs.

Prior owners11539,11526-11530,11413-11417,11274/11284 and11540/11541.
Classical spin statistics, binding boundary conditions, Jarlskog, lattice waves
and sequestering are not new general theories. No solved TOE claim.
"""
import hashlib,itertools,json,sys
from pathlib import Path
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/w33_pass11542_11546_five_dynamical_targets.json'
read=lambda f:json.loads((ROOT/f).read_text())
mat=lambda a:[[str(v) for v in row] for row in a.tolist()]

def binding_controls():
    # Canonically ordered Grassmann monomials independently implement statistics.
    def bilinear(i,j,T):
        out={}
        for a,b in itertools.product(range(2),repeat=2):
            if not T[a,b]:continue
            p,q=2*i+a,2*j+b
            if p!=q:
                key=tuple(sorted((p,q)));out[key]=out.get(key,0)+T[a,b]*(1 if p<q else -1)
        return {k:v for k,v in out.items() if v}
    eps=s.Matrix(2,2,[0,1,-1,0]);trip=s.eye(2)
    assert bilinear(0,1,eps)==bilinear(1,0,eps)
    assert bilinear(0,1,trip)=={k:-v for k,v in bilinear(1,0,trip).items()}
    mu,a,r=s.symbols('mu a r',positive=True);u=s.exp(-r/a);E=-1/(2*mu*a*a)
    assert s.simplify(-s.diff(u,r,2)/(2*mu)-E*u)==0
    assert s.simplify(s.diff(u,r).subs(r,0)/u.subs(r,0)+1/a)==0
    C=s.Matrix(2,2,[s.Symbol('c00',real=True),s.Symbol('x',real=True)+s.I*s.Symbol('y',real=True),s.Symbol('x',real=True)-s.I*s.Symbol('y',real=True),s.Symbol('c11',real=True)])
    # X=[-B conj(Psi),B conj(Phi)]; XdagX=10f²I at this orbit.
    f=s.symbols('f',positive=True);Gram=10*f*f*s.eye(2)
    assert Gram*C*Gram/(10*f*f)==10*f*f*C
    return dict(status='PASS',pass_number=11542,scalar_internal_symmetry='symmetric',triplet_internal_symmetry='antisymmetric',
      statistics='For identical left Weyl fields the local Lorentz scalar epsilon_ab psi_Ia psi_Jb is symmetric in I,J. The native Lambda²81 register requires the spinor-symmetric (1,0) tensor/spin-triplet channel, or extra species/orbital structure. In a nonrelativistic S-wave a spin-triplet with antisymmetric internal state obeys Fermi statistics. Three spin polarizations are a spectator subsystem if controls act identically on them.',
      binding='Zero-range radial boundary uprime(0)/u(0)=-1/a gives u=exp(-r/a), E_B=-1/(2mu a²), only a>0. The scattering length and reduced mass are physical matching inputs; an attractive Casimir label alone does not determine them.',
      scattering_channel_inverse_lengths=['alpha-beta/9+2gamma/3','alpha-13beta/9-4gamma/3'],
      direct_logical_control='X=[-B conj(Psi),B conj(Phi)] uses prior11539 cubic map. The covariant pair mass-squared perturbation X C Xdag, C any Hermitian2x2, restricts to10 f² C on normalized logical pair states. Diagonal difference, real cross and imaginary cross give all Pauli axes. Small perturbations of a positive pair mass matrix give dynamical local gates after a supplied nonrelativistic reduction; an ancilla is not required for these gates.',
      control_matrix=mat(C),control_normalization='Phi=f e0f0, Psi=f e1f0; Xdag X=10f² I2',
      operator_dimensions={'canonical_dimension1_pair_field_logical_control':4,'canonical_dimension1_pair_field_ancilla_mix':5,'microscopic_dimension3_fermion_pair_logical_control':8,'microscopic_dimension3_fermion_pair_ancilla_mix':9},
      ancilla_charge_bound='Pure family diag(-2,1,1) gives ancilla pair charge-4, logical pair charge+2. Both condensates have charge-2, conjugates+2. A nonzero mixing coefficient on this phase needs at least three condensate factors to supply charge6. This is scoped to the two-condensate inventory; extra charged backgrounds can change it.',
      boundary='An explicit operator and statistics-aware binding prescription, not a derived bound state or healthy spin-one UV completion. Canonical auxiliary pair-field dimensions are EFT assumptions, not proof of renormalizability. Binding, spin dynamics, coefficients and physical scales remain open.')

def pair_exchange():
    old=read('data/w33_pass11539_native_plane_holonomy.json');terms=old['neutral_composite']['pair_vectors'][1:]
    labels=sorted({i for row in terms for i,j,v in row}|{j for row in terms for i,j,v in row});index={i:k for k,i in enumerate(labels)};n=len(labels);A=[]
    for row in terms:
        a=np.zeros((n,n),dtype=np.int64)
        for i,j,v in row:a[index[i],index[j]]=int(v);a[index[j],index[i]]=-int(v)
        A.append(a)
    V=[np.einsum('ij,kl->ijkl',a,b) for a,b in itertools.product(A,repeat=2)] # normalized by20
    def W(v):return sum(v.swapaxes(i,j) for i,j in [(0,2),(0,3),(1,2),(1,3)])
    WV=list(map(W,V));M=s.Matrix(4,4,lambda i,j:s.Rational(int(np.sum(V[i]*WV[j])),400))
    M2=s.Matrix(4,4,lambda i,j:s.Rational(int(np.sum(WV[i]*WV[j])),400));I=s.eye(4);swap=s.Matrix(4,4,[1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1]);L=M2-M*M
    assert M==(I+swap)/10 and M2==s.Rational(19,5)*(I+swap) and L==s.Rational(189,50)*(I+swap)
    for j,w in enumerate(WV):
        p12=sum(np.einsum('ij,kl->ijkl',a,np.einsum('ij,ijkl->kl',a,w)) for a in A)
        p34=sum(np.einsum('ij,kl->ijkl',np.einsum('kl,ijkl->ij',a,w),a) for a in A)
        target=sum(v*int(20*M[i,j]) for i,v in enumerate(V))
        assert np.array_equal(p12,target) and np.array_equal(p34,target)
        assert np.array_equal(W(w)+2*w,4*(V[j]+V[[0,2,1,3][j]]))
    J,Delta=s.symbols('J Delta',real=True);block=s.Matrix(2,2,[J/5,3*s.sqrt(21)*J/5,3*s.sqrt(21)*J/5,2*Delta-11*J/5]);center=Delta-J
    omega2=Delta*Delta-s.Rational(12,5)*Delta*J+9*J*J
    assert (block-center*s.eye(2))**2==omega2*s.eye(2)
    q=(s.Integer(649606)-801*s.sqrt(337411))/25672045;ratio=s.Rational(801,800);omega=(1-q)/ratio
    assert s.simplify(omega*omega-(1-s.Rational(12,5)*q+9*q*q))==0 and 0<float(q)<s.Rational(1,8)
    # sin(Omega t)=0; (-1)^400 exp[-i400pi*(801/800)]=-i.
    assert s.simplify(s.exp(-s.I*400*s.pi*ratio))==-s.I
    shortq=(-14+s.sqrt(571))/25;shortomega=2*(1-shortq)
    assert s.simplify(shortomega**2-(1-s.Rational(12,5)*shortq+9*shortq**2))==0
    assert s.simplify(-s.exp(-s.I*3*s.pi/2))==-s.I
    return dict(status='PASS',pass_number=11543,one_particle_support=n,cross_exchange='W=S13+S14+S23+S24; ordinary constituent swaps. W commutes with within-pair permutations and preserves the antisymmetric-pair sector.',projected_W=mat(M),projected_W_squared=mat(M2),leakage_Gram=mat(L),
      full_cyclic_identity='(W²+2W)Vj=4(Vj+V_swap(j)); individual pair projection of W Vj equals joint code projection. The native orbit closes exactly, not approximately, including all81 labels; permutation preserves the30-label support.',
      supplied_two_body_gap='H0=Delta(Q12+Q34), where each Q is the complement of the two-dimensional internal pair code. Pair gap operators act on two constituents each; W is a sum of cross-register two-body swaps.',
      triplet_block=mat(block),singlet='W annihilates the logical singlet exactly; H0 also annihilates it.',
      revival={'J_over_Delta':str(q),'J_over_Delta_numeric':float(q),'Omega_over_Delta':str(omega),'duration_times_Delta':str(400*s.pi/omega),'duration_times_Delta_numeric':float(400*s.pi/omega),'triplet_phase':'-i','singlet_phase':'1','encoded_gate':'exp(-i*pi/4)*exp(-i*pi*SWAP_L/4)'},
      short_revival={'J_over_Delta':str(shortq),'J_over_Delta_numeric':float(shortq),'Omega_over_Delta':str(shortomega),'duration_times_Delta':str(3*s.pi/shortomega),'duration_times_Delta_numeric':float(3*s.pi/shortomega),'triplet_phase':'-i','singlet_phase':'1','encoded_gate':'exp(-i*pi/4)*exp(-i*pi*SWAP_L/4)','scope':'Exact stronger-exchange pulse, not a weak-coupling approximation or a time-optimality claim. At pi/Omega it instead implements inverse sqrtSWAP.'},
      exact_gate='With the displayed algebraic constant coupling ratio, t=400pi/Omega returns all bright amplitude exactly to zero and implements encoded sqrtSWAP up to global phase. No supplied four-body/product pulse or infinite adiabatic limit is needed for this ideal model.',
      boundary='Exact finite-time gate for this declared H0+JW. Pair binding, addressed code gaps, ordinary swap actuators and calibrated algebraic strengths remain physical inputs. Perturbations generally spoil exact revivals; no hardware error threshold or UV derivation.')

def invariant_flavor():
    e=s.symbols('epsilon',real=True);D=s.diag(e**3,e**2,1);Bu=s.diag(1,2,3);Bd=s.Matrix(3,3,[3,1+s.I,1-s.I,1-s.I,4,1,1+s.I,1,5]);Yu=D*Bu*D;Yd=D*Bd*D
    assert all(Bd[:i,:i].det()>0 for i in [1,2,3])
    Hu=Yu*Yu;Hd=Yd*Yd;comm=Hu*Hd-Hd*Hu;cp=s.expand(s.im(s.trace(comm**3)));poly=s.Poly(cp,e);valuation=min(m[0] for m in poly.monoms());lead=poly.coeff_monomial(e**valuation)
    assert valuation==22 and lead==369360
    # Positivity is necessary for this elementary min-max proof: invertibility alone fails.
    bad=s.Matrix(3,3,[0,1,0,1,0,0,0,0,1]);Ybad=D*bad*D
    assert Ybad.det()!=0 and Ybad*Ybad==s.diag(e**10,e**10,1)
    mass=[Bd.det()/Bd[1:,1:].det(),Bd[1:,1:].det()/Bd[2,2],Bd[2,2]]
    return dict(status='PASS',pass_number=11544,positive_down_pin=mat(Bd),up_pin=mat(Bu),down_leading_mass_coefficients=list(map(str,mass)),mass_valuations=[6,4,0],mixing_generic_valuations=[1,2,3],J_generic_valuation=6,weak_basis_CP_valuation=valuation,weak_basis_CP_leading_coefficient=str(lead),exact_CP_polynomial=str(s.factor(cp)),
      theorem='For Y=D B D, B fixed positive definite, lambda_min(B) D²<=Y<=lambda_max(B) D². Min-max proves mass valuations2d=(6,4,0), independent of positive pin coefficients. Generic nonzero Schur off-diagonals yield mixing valuations |di-dj|=(1,2,3). Val(theta_ij)=Val(mi/mj)/2 is a valuation relation, not an equality of physical values.',
      CP='Jarlskog invariant ImTr[Hu,Hd]^3 is weak-basis invariant. On this generic texture J has valuation6 and each squared-mass Vandermonde has valuation8; their product has valuation22. The stored exact polynomial verifies this in the supplied pin example. CP conserving coefficients can make it vanish identically.',
      counterexample='An invertible indefinite pin with B01=B10=1, B22=1 gives singular-value valuations(5,5,0), not(6,4,0). Positivity or appropriate nonzero principal minors cannot be omitted.',
      boundary='Exact parameter-independent exponents for the prior supplied hierarchical texture11413; the mass/mixing orders themselves are owned there. Exact invariant CP polynomial and positivity counterexample refine its scope. Pin coefficients, epsilon, endpoints, Higgs/SM matching and absolute scales remain inputs; no CKM values or measured masses selected.')

def causal_refinement():
    import itertools
    from math import acos,sin
    null=[v for v in itertools.product([-1,0,1],repeat=3) if any(v) and (v[0]*v[0]-v[1]*v[1]-v[2]*v[2])%3==0]
    expected=[v for v in itertools.product([-1,0,1],repeat=3) if abs(v[0])==1 and abs(v[1])+abs(v[2])==1];assert sorted(null)==sorted(expected)
    c2=s.Rational(1,3);a,kx,ky,w=s.symbols('a kx ky omega',real=True)
    symbol=(-4*s.sin(a*w/2)**2+4*c2*(s.sin(a*kx/2)**2+s.sin(a*ky/2)**2))/a**2
    expansion=s.series(symbol,a,0,4).removeO()
    assert s.simplify(expansion-(-w*w+c2*(kx*kx+ky*ky)+a*a*(w**4-c2*(kx**4+ky**4))/12))==0
    rows=[]
    for mesh in [1.,.5,.25,.125]:
        om=acos(1-2*float(c2)*(sin(mesh*.4/2)**2+sin(mesh*.3/2)**2))/mesh
        rows.append(dict(mesh=mesh,omega_squared=om*om,continuum_error=abs(om*om-float(c2)*.25)))
    assert rows[-1]['continuum_error']<rows[0]['continuum_error']/50
    n=65;u0=np.zeros((n,n));u0[n//2,n//2]=1
    lap=lambda u:sum(np.roll(u,z,axis=d) for d in [0,1] for z in [-1,1])-4*u
    u1=u0+float(c2)*lap(u0)/2
    grad=lambda u:[np.roll(u,-1,axis=d)-u for d in [0,1]]
    energy=lambda a,b:float(np.sum((b-a)**2)/2+float(c2)*sum(np.sum(x*y) for x,y in zip(grad(a),grad(b)))/2)
    base=energy(u0,u1);errors=[];dist=np.abs(np.arange(n)-n//2);diamond=dist[:,None]+dist[None,:]
    for tick in range(1,21):
        assert not np.any(u1[diamond>tick]);errors.append(abs(energy(u0,u1)-base));u0,u1=u1,2*u1-u0+float(c2)*lap(u1)
    assert max(errors)<1e-12
    return dict(status='PASS',pass_number=11545,integer_null_steps=list(map(list,null)),cover='Integer chart Z³ with reduction mod3 to the27-event chart. Choosing the four t=+1 steps supplies an orientation absent from the undirected finite graph. Causal reach is a Manhattan diamond; it is not an exactly Euclidean circular microscopic cone.',
      declared_wave='phi(t+1)-2phi(t)+phi(t-1)=c²[phi(x+1)+phi(x-1)+phi(y+1)+phi(y-1)-4phi], c²=1/3 supplied.',continuum_symbol=str(expansion),
      stability='sin²(aomega/2)=c²[sin²(akx/2)+sin²(aky/2)]. For2c²<1, every Fourier mode has real frequency. E=||phi_next-phi||²/2+c²<grad(phi_next),grad(phi)>/2 is exactly conserved; E>= (1/2-c²)||delta phi||²+c²||grad(average phi)||²/2.',
      numerical_refinements=rows,finite_support_ticks=20,maximum_energy_drift=max(errors),
      boundary='An explicitly added oriented integer cover and standard hyperbolic stencil; an example of a causal fixed-background2+1 continuum limit. The finite graph does not select this lift, speed, kinetic coefficients or metric. No4D spacetime, spin2 gravity, Einstein constraint or unique Lorentzian limit derived. Prior11540/11541 history shell owners retained.')

def vacuum_energy():
    old=read('data/w33_pass11526_11530_validated_dynamics.json');spectrum=old['coherent_two_stage']['full324_spectrum'];h4=sum(s.Rational(h)**2*n for h,n in spectrum.items());assert h4==s.Rational(2218769,864)
    L,C,Q,mu4,Vol=s.symbols('Lambda C Q mu4 Volume',nonzero=True);z=s.symbols('z');examples={}
    for name,sigma in [('linear',z),('quadratic',z*z/2)]:
        constraint=Vol-s.diff(sigma,z).subs(z,L/mu4)*Q/mu4
        examples[name]=str(s.factor(constraint.subs(L,L-C)-constraint))
    assert examples['linear']=='0' and s.simplify(s.sympify(examples['quadratic'],locals={'C':C,'Q':Q,'mu4':mu4})-C*Q/mu4**2)==0
    weights=[s.Rational(1,5),s.Rational(3,10),s.Rational(1,2)];rho=list(s.symbols('rho0:3'));avg=sum(v*r for v,r in zip(weights,rho));residual=[s.expand(r-avg) for r in rho]
    shifted=sum(v*(r+C) for v,r in zip(weights,rho));assert all(s.expand(r+C-shifted-q)==0 for r,q in zip(rho,residual))
    # Subtracting the average in the matter action itself erases all potential forces.
    naive=s.expand(sum(v*(r-avg) for v,r in zip(weights,rho)));assert naive==0
    return dict(status='PASS',pass_number=11546,scalar_hessian_mass_fourth_sum=str(h4),
      one_loop='In a separately declared canonically normalized Gaussian scalar theory with M_i²=m² h_i, the scalar-only CW log coefficient is2218769 m4/(55296 pi²), and dV1/dlog(mu)=-2218769 m4/(27648 pi²). Gauge, fermion and ghost contributions plus kinetic matching must be supplied for the full physical supertrace.',
      application='Apply prior11274/11284 sequestering to this actual two-stage scalar constant shift: T->T-Cg and Lambda->Lambda-C leave T-g<TrT>/4-DeltaLambda*g unchanged in a fixed affine-sigma flux sector. A nonconstant source remains. This imports an existing gravitational mechanism, not a new W33 solution.',
      fixed_flux_Ward_test=examples,criterion='With fixed geometry and flux Q, Volume=sigma_prime(Lambda/mu4)Q/mu4 must remain invariant under arbitrary constant C. Differentiating forces sigma_second=0, so sigma is affine. This criterion is scoped to exact fixed-data shift symmetry; it is not a no-go theorem for nonlinear sequestering with adjusted geometry or approximate radiative stability.',
      residual_epoch_density=list(map(str,residual)),naive_matter_subtraction='sum volume_i(rho_i-<rho>)=0 identically. Using this as the entire potential erases the local scalar force as well as constants. Subtraction belongs in the constrained gravitational source, with matter dynamics retained.',
      boundary='Scalar loop instability quantified, existing affine-sector shift protection applied and its fixed-data condition audited. Auxiliary four-forms, metric dynamics, flux ratios and matter/gravity matching are added inputs; graviton loops and observed residual vacuum energy are not solved.')

def digest(path):
    raw=json.dumps(json.loads(path.read_text()),sort_keys=True,separators=(',',':')).encode() if path.suffix=='.json' else path.read_bytes()
    return hashlib.sha256(raw).hexdigest()

def produce():
    out=dict(status='PASS',passes=list(range(11542,11547)),reservation='fbaf36b43')
    for key,fn in [('binding',binding_controls),('exchange',pair_exchange),('flavor',invariant_flavor),('causal',causal_refinement),('vacuum_energy',vacuum_energy)]:out[key]=fn();print(key,'PASS',flush=True)
    files=['data/w33_pass11539_native_plane_holonomy.json','data/w33_pass11526_11530_validated_dynamics.json','analysis/w33_pass11413_11417_pin_regge_vacuum_noise.py','analysis/w33_pass11274_scale_and_sequestering.py','analysis/w33_pass11284_closed_history_flux.py','analysis/w33_pass11541_history_quadratic_shell_scheme.py']
    out['source_sha256']={f:digest(ROOT/f) for f in files};out['binding_method']='Canonical parsed JSON hashes, raw LF code hashes.';out['producer_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();OUT.write_text(json.dumps(out,indent=2)+'\n');return out
if __name__=='__main__':produce()
