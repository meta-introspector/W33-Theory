"""Native exceptional two-plane: information quotient, control and gauge audit.

Classical holonomic computation, coherent orbits and Casimir calculus are prior
art. This application binds the actual native27 to11526-11530's e0/e1 phase.
"""
import hashlib,itertools,json,sys
from functools import lru_cache
from pathlib import Path
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/w33_pass11539_native_plane_holonomy.json'
sys.path.insert(0,str(ROOT/'analysis'))

@lru_cache(None)
def native():
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
    mats=[s.SparseMatrix(a.real.astype(int).tolist())+s.I*s.SparseMatrix(a.imag.astype(int).tolist()) for a in he]
    G=s.Matrix(np.einsum('aij,bji->ab',he,he).real.astype(int));IE=G.inv()
    return he,mats,IE,d

mat=lambda A:[[str(x) for x in row] for row in A.tolist()]
complex_mat=lambda A:dict(real=np.asarray(A).real.tolist(),imag=np.asarray(A).imag.tolist())

def geometry():
    he,mats,IE,d=native();W=s.Matrix([[int(a[i,i].real) for a in he] for i in range(27)]);pair=W*IE*W.T
    K=s.diag(*list(pair[:,0]+pair[:,1]));I=s.eye(27);P=36*(K-I/9)*(K+I/18)*(K+2*I/9)
    assert P==s.diag(1,1,*([0]*25)) and not np.any(d[:,[0,1]][:,:,[0,1]])
    cols=np.array([np.r_[1j*a[:,0],1j*a[:,1]] for a in he]).T
    tangent=s.Matrix(np.r_[cols.real,cols.imag].astype(int));ranks=[];omegas=[]
    for a,b in [(1,1),(1,2)]:
        omega=s.Matrix((2*np.imag(a*he[:,:,0].conj()@he[:,:,0].T+b*he[:,:,1].conj()@he[:,:,1].T)).astype(int))
        diagonal=a*pair[:,0]+b*pair[:,1]
        central=[j for j,g in enumerate(he) if all(diagonal[k]==diagonal[l] for k,l in np.argwhere(g!=0))]
        assert omega.rank()==78-len(central)
        ranks.append(dict(radii_squared=[a,b],moment_centralizer_dimension=len(central),Berry_rank=omega.rank(),ordered_orbit_dimension=tangent.rank(),null_on_ordered_orbit=tangent.rank()-omega.rank()))
        omegas.append(omega)
    assert [x['Berry_rank'] for x in ranks]==[50,52] and tangent.rank()==54
    mixing=[j for j,g in enumerate(he) if g[0,1]!=0];assert len(mixing)==2
    return dict(status='PASS',native_joint_moment_diagonal=list(map(str,K.diagonal())),native_band_projector=mat(P),band_Hamiltonian='H0=5I/18-K, K=K(e0)+K(e1)',band_spectrum={'0':2,'1/6':10,'1/3':10,'1/2':5},gap='1/6',normalization='Control vectors e0,e1 have norm1. The certified scalar vacuum uses e0/sqrt2,e1/sqrt2, so its joint moment and corresponding gap are half these normalized values. No physical energy scale is selected.',projector_polynomial='P=36(K-I/9)(K+I/18)(K+2I/9)',rank_comparison=ranks,equal_Berry_matrix=mat(omegas[0]),E6_frame_tangent=mat(tangent),mixing_generator_indices=mixing,
      amplitude_tradeoff='On the two native basis states the energy difference for K=r0 K(e0)+r1 K(e1) is (r0-r1)/6. Unequal positive radii restore two Berry directions while splitting the logical band; no physical mass ratio is thereby selected.',
      scope='The ordered E6 orbit has dimension54; the equal-weight moment/Slater ray orbit has dimension50. Four missing frame directions are a quotient/fiber issue, not new scalar masses. In the gauged vacuum model all scalar gauge-orbit directions are redundant already. The first-order Berry model here is a separately declared control/quantization model, not the relativistic scalar kinetic action.')

def curvature():
    he,mats,IE,d=native();Q=s.diag(0,0,*([1]*25));pairs=[];vectors=[];sn,cs=s.symbols('sn cs',real=True)
    for a,b in itertools.combinations(range(78),2):
        F=(mats[a]*Q*mats[b]-mats[b]*Q*mats[a])[:2,:2];H=-s.I*F
        vector=[H[0,0],H[1,1],s.re(H[0,1]),s.im(H[0,1])]
        if s.Matrix(vectors+[vector]).rank()>len(vectors):
            assert mats[a][:2,:2]==mats[b][:2,:2]==s.zeros(2)
            assert mats[a]**3==mats[a] and mats[b]**3==mats[b]
            E=s.eye(27)+s.I*sn*mats[b]+(cs-1)*mats[b]**2
            rotated=(E.conjugate().T*mats[a]*E)[:2,:2].applyfunc(s.expand)
            pairs.append(dict(generators=[a,b],curvature=mat(F),Hermitian_coordinates=list(map(str,vector)),rectangle_rotated_generator=mat(rotated)))
            vectors.append(vector)
        if len(vectors)==4:break
    assert s.Matrix(vectors).rank()==4
    return dict(status='PASS',witnesses=pairs,curvature_rank=4,curvature_formula='F_ab=P(Ta Q Tb-Tb Q Ta)P',rectangle='U(s,t)=exp(i s Ta)exp(i t Tb), path (0,0)->(eps,0)->(eps,eps)->(0,eps)->(0,0). For these normal pairs A_s(t)=i P exp(-itTb)Ta exp(itTb)P and A_t=0; geometric gate=exp(+i eps P exp(-i eps Tb)Ta exp(i eps Tb)P).',
      theorem='Four curvature values span u(2). Ambrose-Singer therefore gives full connected U(2) holonomy of this rank-two bundle on the declared E6 control orbit. This is the standard holonomic-control theorem instantiated in the native representation, not a new general theorem.',
      boundary='Controllable anchored relative plane motion and an adiabatic gapped single-carrier band are premises. Pure simultaneous E6 gauge changes do not implement gates. Radiative protection, physical control actuators and errors are unproved.')

def exchange_operator():
    he,mats,IE,d=native();entries={}
    for a in range(78):
        for b in range(78):
            if IE[a,b]:
                for (i,j),x in mats[a].todok().items():
                    for (k,l),y in mats[b].todok().items():
                        ij=(27*i+k,27*j+l);entries[ij]=entries.get(ij,0)+IE[a,b]*x*y
    C=s.SparseMatrix(729,729,{ij:s.expand(v) for ij,v in entries.items() if v!=0});I=s.SparseMatrix(s.eye(729))
    swap=s.SparseMatrix(729,729,{(27*i+j,27*j+i):1 for i in range(27) for j in range(27)})
    assert s.Rational(9,2)*C*C+s.Rational(11,2)*C-s.Rational(4,9)*I==swap
    assert (C+13*I/9)*(C+I/9)*(C-2*I/9)==s.zeros(729)
    assert s.trace(C)==0 and s.trace(C*C)==78
    D=s.SparseMatrix(729,27,{(27*i+j,k):int(d[i,j,k]) for i,j,k in np.argwhere(d!=0)})
    assert D.conjugate().T*D==10*s.eye(27)
    assert 18*C==I+3*swap-3*D*D.conjugate().T
    X=s.Matrix(2,2,[0,1,1,0]);Y=s.Matrix(2,2,[0,-s.I,s.I,0]);Z=s.diag(1,-1);i2=s.eye(2)
    S=s.Matrix(4,4,[1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1]);generators=[s.I*s.kronecker_product(a,b) for a,b in [(X,i2),(Y,i2),(Z,i2),(i2,X),(i2,Y),(i2,Z)]]+[s.I*(S-s.eye(4)/2)]
    vector=lambda A:list(A.applyfunc(s.re))+list(A.applyfunc(s.im))
    basis=[];rows=[];steps=[]
    for A in generators:
        if s.Matrix(rows+[vector(A)]).rank()>len(rows):basis.append(A);rows.append(vector(A));steps.append('supplied local or exchange')
    while len(basis)<15:
        found=False
        for a,b in itertools.combinations(range(len(basis)),2):
            A=basis[a]*basis[b]-basis[b]*basis[a]
            if s.Matrix(rows+[vector(A)]).rank()>len(rows):
                basis.append(A);rows.append(vector(A));steps.append([a,b]);found=True;break
        assert found
    return dict(status='PASS',Casimir_entries=[[i,j,str(v)] for (i,j),v in sorted(C.todok().items())],Casimir_spectrum={'-13/9':27,'-1/9':351,'2/9':351},swap_polynomial='SWAP27=(9C²+11C)/2-4I/9, C=sum(G^-1)_ab Ta tensor Tb',encoded_swap=mat(S),su4_basis=list(map(mat,basis)),su4_generation_steps=steps,control_Lie_rank=15,cubic_fierz='D_(27i+j,k)=d_ijk; Ddag D=10I27; 18C=I+3SWAP27-3D Ddag; SWAP27=6C+D Ddag-I/3',
      exact_entangler='exp(-i pi SWAP27/4) preserves the aligned P tensor P band exactly and sends |01> to (|01>-i|10>)/sqrt2. Local U(2) holonomy plus tunable exchange has full su(4) control algebra.',
      gapped_pair='Hpair=Delta(Q tensor I+I tensor Q)+J SWAP27. For aligned registers Delta>2|J| gives a gap>=Delta-2|J| to the complement. All local loops return to the same base plane before exchange.',
      boundary='A declared exchange Hamiltonian built from native E6 invariant operators, not a derived interaction strength or hardware pulse. Two-register occupation, alignment, anchored local controls and accessible exchange are additional premises.')

def adiabatic_control(curv):
    from scipy.linalg import expm
    from scipy.integrate import solve_ivp
    he,*_=native();a,b=curv['witnesses'][1]['generators'];A=he[a];B=he[b];i=np.eye(27);Q=np.diag([0,0]+[1]*25);eps=.6
    rot=lambda T,x:i+1j*np.sin(x)*T+(np.cos(x)-1)*(T@T)
    K=(rot(B,eps).conj().T@A@rot(B,eps))[:2,:2];target=expm(1j*eps*K);rows=[]
    for total in [64.,256.,512.]:
        psi=i[:,:2].astype(complex);duration=total/4
        for leg in range(4):
            def rhs(time,y):
                u=time/duration
                st=[(eps*u,0),(eps,eps*u),(eps*(1-u),eps),(0,eps*(1-u))][leg]
                U=rot(A,st[0])@rot(B,st[1]);H=U@Q@U.conj().T
                return (-1j*H@y.reshape(27,2)).ravel()
            sol=solve_ivp(rhs,(0,duration),psi.ravel(),method='DOP853',rtol=1e-10,atol=1e-12);assert sol.success
            psi=sol.y[:,-1].reshape(27,2)
        rows.append(dict(duration=total,leakage_operator_norm_squared=float(np.linalg.norm(psi[2:,:],2)**2),geometric_operator_error=float(np.linalg.norm(psi[:2,:]-target,2)),norm_error=float(np.linalg.norm(psi.conj().T@psi-np.eye(2))),encoded_action=complex_mat(psi[:2,:])))
    assert rows[-1]['geometric_operator_error']<rows[0]['geometric_operator_error'] and max(r['norm_error'] for r in rows)<1e-7
    return dict(status='PASS',epsilon=eps,geometric_gate=complex_mat(target),history=rows,scope='Finite-time floating Schrödinger controls, with supplied Delta=1. Exact curvature/holonomy is geometric; finite-time gates have leakage and coherent error. These are convergence witnesses, not certified error bounds or hardware thresholds.')


@lru_cache(None)
def full_stabilizer():
    _,mats,_,_=native();family=[]
    for i,j in itertools.combinations(range(3),2):
        a=s.zeros(3);a[i,j]=a[j,i]=1;family.append(a)
        a=s.zeros(3);a[i,j]=-s.I;a[j,i]=s.I;family.append(a)
    family.extend([s.diag(1,-1,0),s.diag(1,1,-2)])
    full=[s.kronecker_product(a,s.eye(3)) for a in mats]+[s.kronecker_product(s.eye(27),a) for a in family]
    cart=[j for j,a in enumerate(full) if a==s.diag(*a.diagonal())]
    w=s.Matrix([[full[j][i,i] for j in cart] for i in range(81)])
    cs=s.Matrix.vstack(w[0,:],w[3,:]).nullspace();hc=s.Matrix.hstack(*cs);weights=w*hc
    roots=[j for j,a in enumerate(full) if j not in cart and not any(a[:,0]) and not any(a[:,3])]
    stabilizer=[full[j] for j in roots]+[sum((v[k]*full[j] for k,j in enumerate(cart)),s.zeros(81)) for v in cs]
    assert len(stabilizer)==28 and len(roots)==22 and len(cs)==6
    return full,stabilizer,weights,roots,cart,cs

def wedge_action(a,pairs):
    entries={}
    for col,(i,k) in enumerate(pairs):
        for r in range(81):
            if r!=k and a[r,i]:
                key=(tuple(sorted((r,k))),col);entries[key]=entries.get(key,0)+a[r,i]*(1 if r<k else -1)
            if r!=i and a[r,k]:
                key=(tuple(sorted((i,r))),col);entries[key]=entries.get(key,0)+a[r,k]*(1 if i<r else -1)
    return {k:s.expand(v) for k,v in entries.items() if v!=0}

def neutral_composite():
    full,stab,weights,roots,cart,cs=full_stabilizer()
    pairs=list(itertools.combinations(range(81),2))
    neutral=[p for p in pairs if all(weights[p[0],k]+weights[p[1],k]==0 for k in range(6))]
    blocks=[]
    for j in roots:
        entries=wedge_action(full[j],neutral);rows=sorted({p for p,c in entries})
        blocks.append(s.Matrix([[entries.get((p,c),0) for c in range(len(neutral))] for p in rows]))
    constraint=s.Matrix.vstack(*blocks);basis=constraint.nullspace()
    assert len(neutral)==21 and len(basis)==3
    singles=[i for i in range(81) if all(weights[i,k]==0 for k in range(6))]
    one=s.Matrix.vstack(*[a[:,singles] for a in stab]).nullspace();assert len(one)==2
    vectors=[];norms=[]
    for v in basis:
        sparse={neutral[i]:v[i] for i in range(len(neutral)) if v[i]};norm=(v.conjugate().T*v)[0];norms.append(int(norm))
        for a in stab:
            action=wedge_action(a,list(sparse));out={}
            for (pair,j),value in action.items():out[pair]=out.get(pair,0)+value*list(sparse.values())[j]
            assert all(s.expand(x)==0 for x in out.values())
        vectors.append([[i,j,str(c)] for (i,j),c in sparse.items()])
    assert norms==[1,10,10]
    # Independent dark-bundle curvature in the three singlet coordinates.
    P=s.diag(1,1,0);Q=s.eye(3)-P;normal=[]
    for i in range(2):
        a=s.zeros(3);a[i,2]=a[2,i]=1;normal.append(a)
        a=s.zeros(3);a[i,2]=-s.I;a[2,i]=s.I;normal.append(a)
    values=[]
    for a,b in itertools.combinations(normal,2):
        H=-s.I*(a*Q*b-b*Q*a)[:2,:2]
        values.append([H[0,0],H[1,1],s.re(H[0,1]),s.im(H[0,1])])
    assert s.Matrix(values).rank()==4
    return dict(status='PASS',carrier='Lambda²(27 tensor 3), wedge basis |i> wedge |j>, i=3*native+family',dimension=3240,stabilizer_dimension=28,root_generator_indices=roots,cartan_generator_indices=cart,cartan_kernel_basis=list(map(mat,cs)),zero_weight_pair_count=21,root_constraint_rank=constraint.rank(),full_neutral_one_particle_dimension=2,full_neutral_pair_dimension=3,pair_vectors=vectors,pair_norm_squared=norms,neutral_dark_curvature_coordinates=mat(s.Matrix(values)),neutral_dark_curvature_rank=4,
      exact_map='V has the three pair_vectors as columns divided by sqrt(pair_norm_squared), so Vdag V=I3. H(b)=Delta[I-VVdag+V b bdag Vdag], bdag b=1. Its spectrum is 0 twice and Delta 3238 times; Pdark=V(I3-b bdag)Vdag. It commutes with all28 unbroken generators. CP² bright-direction controls have full U2 dark holonomy.',
      covariant_extension='On the connected-stabilizer vacuum orbit cover G/H0, define H_field(g vacuum,b)=rho_wedge(g)H_ref(b)rho_wedge(g)dag. This is independent of g->gh because H_ref commutes with the connected stabilizer H0. Descent through any disconnected stabilizer component requires a further check. This explicitly specified orbit-effective operator is not a derived polynomial UV interaction.',
      entangler='For two pair registers in slots (1,2),(3,4), the whole-register SWAP is S13(81) S24(81), restricted to antisymmetric pairs. exp(-i pi S13 S24/4) is encoded sqrtSWAP. It is a supplied four-body/product Hamiltonian, not the sequence of two swaps. SWAP81=SWAP27 tensor SWAP3 after regrouping indices; SWAP3=CF+I/3 for trace-dual SU3 Casimir.',
      obstruction='The full unbroken gauge singlet space in native81 has dimension2. A rank-two band wholly inside it is constant. Its two-particle exterior sector has exactly three singlets and therefore supplies an internal neutral ancilla without charged excursions.',
      boundary='Pair binding, spatial/spin statistics, two-particle occupation, independently tunable relative bright controls, the four-body interaction and physical energy/error scales are supplied premises, not TOE predictions. No universal physical processor, mass spectrum, gravity or dynamical spacetime is claimed.')


def pair_resources(data):
    from collections import defaultdict
    he,_,_,d=native();pairs=list(itertools.combinations(range(81),2));pi={v:i for i,v in enumerate(pairs)}
    def epsilon(f,g,h):
        if len({f,g,h})<3:return 0
        return -1 if sum(a>b for a,b in itertools.combinations([f,g,h],2))%2 else 1
    entries={}
    for i,j in pairs:
        ni,f=divmod(i,3);nj,g=divmod(j,3)
        for h in range(3):
            e=epsilon(f,g,h)
            if e:
                for m in np.flatnonzero(d[ni,nj]):entries[pi[i,j],3*int(m)+h]=int(d[ni,nj,m])*e
    B=s.SparseMatrix(3240,81,entries);assert B.conjugate().T*B==10*s.eye(81)
    for a in he:
        assert not np.any(np.einsum('ai,ajk->ijk',a,d)+np.einsum('aj,iak->ijk',a,d)+np.einsum('ak,ija->ijk',a,d))
    vectors=[{(i,j):s.sympify(v) for i,j,v in terms} for terms in data['neutral_composite']['pair_vectors']]
    assert {(pairs[r]):v for (r,c),v in B.todok().items() if c==0}==vectors[2]
    assert {(pairs[r]):v for (r,c),v in B.todok().items() if c==3}=={k:-v for k,v in vectors[1].items()}
    C=defaultdict(list)
    for i,j,v in data['exchange']['Casimir_entries']:C[j].append((i,s.sympify(v)))
    eigs=[]
    for kind,expected in [('E6',[s.Rational(-1,9),s.Rational(-13,9),s.Rational(-13,9)]),('family',[s.Rational(2,3),s.Rational(-4,3),s.Rational(-4,3)])]:
        for v,eig in zip(vectors,expected):
            out=defaultdict(lambda:s.S.Zero)
            for (i,j),z in v.items():
                k,f=divmod(i,3);l,g=divmod(j,3)
                terms=[(3*a+f,3*b+g,z*c) for row,c in C[27*k+l] for a,b in [divmod(row,27)]] if kind=='E6' else [(3*k+g,3*l+f,z),(i,j,-z/3)]
                for a,b,c in terms:
                    if a!=b:out[tuple(sorted((a,b)))]+=c*(1 if a<b else -1)
            assert {k:s.expand(v) for k,v in out.items() if v!=0}=={k:eig*z for k,z in v.items()}
        eigs.append(list(map(str,expected)))
    canonical=json.loads((ROOT/'extracted_v13/W33-Theory-master/artifacts/canonical_su3_gauge_and_cubic.json').read_text())
    old={tuple(sorted(row['triple'])):int(row['sign']) for row in canonical['solution']['d_triples']}
    current={tuple(map(int,ix)):int(d[tuple(ix)]) for ix in np.argwhere(d!=0) if ix[0]<ix[1]<ix[2]}
    assert current==old
    return dict(status='PASS',prior_bracket_owner='analysis/w33_e6_cubic_hybrid81_transport.py; this B is the adjoint of the existing graded E8 matter bracket in the identical canonical cubic gauge.',canonical_cubic_equal=True,cubic_pair_map='B_(if wedge jg,mh)=d_ijm epsilon_fgh. Bdag B=10I81; B intertwines conjugate81 into Lambda²81.',map_entries=[[r,c,str(v)] for (r,c),v in sorted(B.todok().items())],native_pair_Casimir=eigs[0],family_pair_Casimir=eigs[1],
      covariant_frame='For unit normalized condensates u,v on the certified orbit: V=[u wedge v, -B conjugate(v)/sqrt10, B conjugate(u)/sqrt10]. At e0f0,e1f0 it is exactly the three normalized pair vectors. Each column transforms covariantly under the full gauge group, so H_ref(b) extended by this V is well defined even through disconnected stabilizer components that fix both labelled condensates.',
      sector_split='The ancilla is in (Lambda²27,Sym²3), dimension351*6; the two logical channels are in (conjugate27,conjugate3), dimension27*3. General decomposition Lambda²(27 tensor3)=(351 tensor6)+(conjugate27 tensor conjugate3)+(351prime tensor conjugate3) is multiplicity free. Full-gauge invariant operators are scalar on each irrep and cannot mix ancilla with logical channels. Condensate-referenced controls are essential.',
      boundary='The explicit native cubic map removes the orbit-section ambiguity, but neither arbitrary bright-direction actuator amplitudes nor binding nor physical energy scales follow from it. Normalized condensates introduce fixed-orbit normalization; the model is an effective interaction prescription, not a proved microscopic Lagrangian.')

def produce():
    data={'status':'PASS','pass':11539,'reservation':'b6a492645'}
    for name,fn in [('geometry',geometry),('curvature',curvature),('exchange',exchange_operator)]:data[name]=fn();print(name,'PASS',flush=True)
    data['neutral_composite']=neutral_composite();print('neutral_composite PASS',flush=True)
    data['pair_resources']=pair_resources(data);print('pair_resources PASS',flush=True)
    data['neutral_composite']['covariant_extension']=data['pair_resources']['covariant_frame']
    data['adiabatic']=adiabatic_control(data['curvature']);print('adiabatic PASS',flush=True)
    data['occupation_gauge_audit']=dict(status='PASS',occupation='One excitation in a two-dimensional band is a qubit. Filling both levels gives wedge²(C²), dimension1, and U acts by det(U); SU2 data disappear. Empty and doubly occupied sectors are not logical qubits.',gauge='D=partial-iA. Under a pure change U, A=-i dotU Udag and P=U P0 Udag. Then DP=0 and Udag D U=0: apparent frame Berry gates cancel against the transformed gauge connection. Nonzero control holonomy requires changing the condensate plane relative to a fixed connection/reference, not moving every object by one gauge transformation.',scope='No identification of bundle SU2 with weak isospin or gravity. No physical universal quantum processor follows until carrier occupation, relative controls, exchange and error budgets are instantiated.')
    sources=['analysis/w33_pass11384_11388_native_dynamics.py','data/w33_pass11384_native_inputs.json','analysis/w33_pass11526_11530_validated_dynamics.py','analysis/w33_e6_cubic_hybrid81_transport.py','extracted_v13/W33-Theory-master/artifacts/canonical_su3_gauge_and_cubic.json']
    def bound_bytes(path):
        return json.dumps(json.loads(path.read_text()),sort_keys=True,separators=(',',':')).encode() if path.suffix=='.json' else path.read_bytes()
    data['source_binding_method']='JSON inputs: canonical parsed JSON, sorted keys and compact separators; code: raw LF source bytes.'
    data['source_sha256']={p:hashlib.sha256(bound_bytes(ROOT/p)).hexdigest() for p in sources}
    data['producer_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    OUT.write_text(json.dumps(data,indent=2)+'\n');return data
if __name__=='__main__':produce()
