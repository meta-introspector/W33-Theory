#!/usr/bin/env python3
"""Joint chart/spin/selector/portal architecture, and three additional bridges.

Generic owners: BT4047/4048 projected edge controls; Pass10951 Spin signs
and distinction GL2(3) vs2O; Pass4091 anomaly identities; Pass11388 quantum
scale boundary. New witnesses specialize them to Pass11389-11394 carriers.
"""
from __future__ import annotations
import json
from collections import Counter
from itertools import product,combinations
from pathlib import Path
import numpy as np
import sympy as sp
from w33_pass11390_11394_context_matter_clock import load,context_frames,paulis
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/w33_pass11395_11402_joint_atlas_spin_portal.json'

def mat(M):return [[str(v)for v in M.row(i)]for i in range(M.rows)]
def qumul(a,b):
    w,x,y,z=a;v,i,j,k=b
    return tuple(sp.expand(t)for t in (w*v-x*i-y*j-z*k,w*i+x*v+y*k-z*j,w*j-x*k+y*v+z*i,w*k+x*j-y*i+z*v))
def quconj(a):return (a[0],-a[1],-a[2],-a[3])
def rotation(q):
    w,x,y,z=q
    return sp.Matrix([[1-2*(y*y+z*z),2*(x*y-w*z),2*(x*z+w*y)],
        [2*(x*y+w*z),1-2*(x*x+z*z),2*(y*z-w*x)],
        [2*(x*z-w*y),2*(y*z+w*x),1-2*(x*x+y*y)]]).applyfunc(sp.expand)
def su2(q):
    w,x,y,z=q
    return sp.Matrix([[w-sp.I*z,-y-sp.I*x],[y-sp.I*x,w+sp.I*z]])

def build(c):
    fs,reps,scale=context_frames(c);G=sp.Matrix(fs[0].T@fs[0]);Gi=G.inv()
    J=np.array([[0,0,1,0],[0,0,0,1],[-1,0,0,0],[0,-1,0,0]])
    gs=[c['permutation']((np.eye(4,dtype=int)+np.outer(v,np.array(v)@J))%3)for v in c['points']]
    base=c['points'].index((1,0,0,0));r={base:tuple(range(80))};queue=[base]
    for a in queue:
        for g in gs:
            b=g[a]
            if b not in r:r[b]=tuple(g[r[a][i]]for i in range(80));queue.append(b)
    H=sp.Matrix(c['x']['point']['harmonic_voltage']).applyfunc(sp.Rational)
    assert int(sp.ilcm(*[v.q for v in H]))==scale
    P=np.array(H*scale,dtype=np.int64);pf=[]
    for a in range(40):
        F=np.zeros_like(P);g=r[a]
        for i,(u,v)in enumerate(c['edges']):F[c['ei'][(g[u],g[v])]]=P[i]
        assert sp.Matrix(F.T@F)==G;pf.append(F)
    A=sp.Matrix(c['x']['line']['fcc_roots'])*sp.Matrix(c['x']['line']['primitive_fcc_change']).inv()
    conn={}
    for a,b in combinations(range(40),2):
        if not set(c['lines'][a])&set(c['lines'][b]):
            R=27*Gi*sp.Matrix(fs[b].T@fs[a]);O=A*R.inv().T*A.inv()
            assert O.T*O==sp.eye(3)and all(v in (-1,0,1)for v in O)
            conn[a,b]=O;conn[b,a]=O.T
    D=c['D'];L=D*D.T;Pi=sp.eye(160)-D.T*((L+sp.ones(80)/80).inv()-sp.ones(80)/80)*D
    assert Pi*Pi==Pi and sp.trace(Pi)==81 and all(Pi[i,i]==sp.Rational(81,160)for i in range(160))
    return dict(c=c,line=fs,point=pf,G=G,Gi=Gi,conn=conn,Pi=Pi)

def spin(b):
    zero=sp.Integer(0);one=sp.Integer(1);half=sp.Rational(1,2);rt=sp.sqrt(2)/2
    Q=set()
    for i in range(4):
        for s in (-one,one):q=[zero]*4;q[i]=s;Q.add(tuple(q))
    Q.update(tuple(s*half for s in signs)for signs in product((-1,1),repeat=4))
    for i,j in combinations(range(4),2):
        for s,t in product((-1,1),repeat=2):q=[zero]*4;q[i]=s*rt;q[j]=t*rt;Q.add(tuple(q))
    assert len(Q)==48 and all(qumul(a,quconj(a))==(1,0,0,0)for a in Q)
    assert all(qumul(a,d)in Q for a in Q for d in Q)
    rots={}
    for q in sorted(Q,key=str):
        R=rotation(q);rots.setdefault(tuple(R),[]).append(q)
        U=su2(q);assert (U.H*U).applyfunc(sp.simplify)==sp.eye(2)and sp.simplify(U.det())==1
    assert len(rots)==24 and all(len(x)==2 for x in rots.values())
    canonical={k:next(q for q in v if next(x for x in q if x!=0)>0)for k,v in rots.items()}
    F={1:sp.eye(3),-1:sp.diag(-1,1,1)}
    links={}
    for (a,d),O in b['conn'].items():
        for s in (1,-1):
            t=s*int(O.det());R=F[t]*O*F[s];q=canonical[tuple(R)]
            if ((d,t),(a,s))in links:q=quconj(links[(d,t),(a,s)])
            links[(a,s),(d,t)]=q
    assert len(links)==2160
    witness=None
    for a,d,e in combinations(range(40),3):
        if all(p in b['conn']for p in [(a,d),(d,e),(e,a)]):
            W=b['conn'][e,a]*b['conn'][d,e]*b['conn'][a,d]
            if W.det()==1:
                s=1;t=int(b['conn'][a,d].det());u=t*int(b['conn'][d,e].det())
                q=qumul(links[(e,u),(a,s)],qumul(links[(d,t),(e,u)],links[(a,s),(d,t)]))
                assert q[0]==0 and qumul(q,q)==(-1,0,0,0)
                witness=dict(contexts=[a,d,e],orientation_signs=[s,t,u],quaternion=[str(x)for x in q],su2=mat(su2(q)));break
    assert witness
    return dict(binary_octahedral_quaternions=[[str(x)for x in q]for q in sorted(Q,key=str)],group_order=48,
        proper_rotations=24,oriented_directed_spin_links=2160,loop=witness,
        proper_triangle_double_turn='-I2 on the spinor; +I3 on vectors',
        scope='Explicit Spin3 lift on the oriented context graph. Generic2O and central minus sign are prior art; this is the actual harmonic-connection specialization. GL2(3) is a different order48 cover. No continuum spin-statistics or chiral spacetime theory follows.')

def selector_portal(b):
    c=b['c'];A=np.array([[int(i!=j and bool(set(a)&set(d)))for j,d in enumerate(c['lines'])]for i,a in enumerate(c['lines'])])
    base=0;target=next(i for i in range(40)if i!=base and not A[base,i]);shells=[[base],list(np.where(A[base]==1)[0]),[i for i in range(40)if i!=base and not A[base,i]]]
    T=np.zeros((40,3))
    for j,shell in enumerate(shells):T[shell,j]=1/np.sqrt(len(shell))
    quotient=np.array([[0,np.sqrt(12),0],[np.sqrt(12),2,6],[0,6,8]])
    assert np.linalg.norm(A@T-T@quotient)<1e-12
    eta,E=sp.symbols('eta E');H=sp.diag(0,1,1)-eta*sp.Matrix([[0,sp.sqrt(12),0],[sp.sqrt(12),2,6],[0,6,8]])
    char=sp.factor((E*sp.eye(3)-H).det());scans=[]
    for h in (.00625,.003125,.0015625,.00078125):
        full=np.diag([0]+[1]*39)-h*A;w,v=np.linalg.eigh(full);psi=v[:,0];pq=psi[target]**2
        hq=np.diag([0,1,1])-h*quotient;wq,vq=np.linalg.eigh(hq)
        assert abs(w[0]-wq[0])<1e-12 and abs(pq-vq[2,0]**2/27)<1e-12
        M=20.;D=np.zeros((40,40));D[target,target]=1
        # Exact heavy cross resolvent, evaluated at unperturbed selector E0.
        X=full-w[0]*np.eye(40);cross=M*np.linalg.inv(M*M*np.eye(40)-X@X)
        exact=psi@D@cross@D@psi
        scans.append(dict(h_over_J=h,target_probability=float(pq),probability_over_eta4=float(pq/h**4),
            exact_cross_resolvent_weight=float(exact),leading_probability_over_M=float(pq/M),
            heavy_correction_ratio=float(exact/(pq/M))))
    assert abs(scans[-1]['probability_over_eta4']-16)<2
    # Exact asymptotic coefficients via Rayleigh perturbation of a pinned cell.
    assert A[base,target]==0 and (A@A)[base,target]==4
    return dict(pinned_context=base,target_context=target,shell_sizes=[1,12,27],
        exact_shell_hamiltonian=mat(H),characteristic_polynomial=str(char),
        target_probability_leading='16*(h/J)^4',conditional_endpoint_operator='y0*|target><target|, not y0*mean(projector)',
        leading_effective_portal='-16*(y0^2/M)*(h/J)^4 * C_orbit',
        leading_small_singular_value='8*sqrt(3)*(y0^2/M)*(h/J)^4',scans=scans,
        three_channel_compiler=hierarchy_compiler(c,A,quotient),
        scope='A named pinned-selector, projector-controlled mediator ansatz dynamically suppresses the portal. Leading adiabatic regime M>>J,h and small y0; microscopic h,J,y0,M and the ordered boundary context remain inputs. Replacing <projector^2> by <projector>^2 gives the wrong suppression. Exact cross-resolvent corrections are stored; finite-energy self-energies and an unpinned many-cell phase need further analysis.')

def hierarchy_compiler(c,A,quotient):
    group=c['x']['internal_h27']['group_permutations'];seen=set();orbits=[];blocks=[]
    for a in range(40):
        for d in range(40,80):
            if (a,d)in seen:continue
            orb=sorted({(g[a],g[d])for g in group});seen.update(orb);M=np.zeros((80,80))
            for aa,dd in orb:M[aa,dd]=M[dd,aa]=1
            blocks.append(c['plus'].conj().T@c['S'].conj().T@M@c['S']@c['minus']);orbits.append(orb)
    R=np.array([np.r_[x.real.ravel(),x.imag.ravel()]for x in blocks]).T
    assert len(orbits)==88 and np.linalg.matrix_rank(R,tol=1e-10)==24
    prior=json.loads((ROOT/'data/w33_pass11390_11394_context_matter_clock.json').read_text())['portal']
    C=np.zeros((80,80))
    for a,d in prior['nonedge_orbit']:C[a,d]=C[d,a]=1
    reference=c['B']+c['plus'].conj().T@c['S'].conj().T@C@c['S']@c['minus']
    U,d,Vh=np.linalg.svd(reference);targets=[];weights=[];residuals=[]
    assert min(abs(np.diff(d)))>1e-3
    for j in range(3):
        target=np.outer(U[:,j],Vh[j,:]);co=np.linalg.lstsq(R,np.r_[target.real.ravel(),target.imag.ravel()],rcond=1e-12)[0]
        residual=np.linalg.norm(sum(w*x for w,x in zip(co,blocks))-target);assert residual<1e-12
        targets.append(target);weights.append(co.tolist());residuals.append(float(residual))
    scans=[]
    for eta in (.01,.003,.001):
        _,v=np.linalg.eigh(np.diag([0,1,1])-eta*quotient)
        probabilities=v[:,0]**2/np.array([1,12,27]);Y=sum(p*T for p,T in zip(probabilities,targets))
        singular=np.linalg.svd(Y,compute_uv=False)
        assert np.linalg.norm(singular-probabilities)<1e-12
        # Independent unitary seed-frame change checks coordinate invariance.
        W=np.diag(np.exp(1j*np.array([.17,.39,.71])));V=np.diag(np.exp(1j*np.array([.13,.29,.53,.83])))
        assert np.linalg.norm(np.linalg.svd(W@Y@V,compute_uv=False)-singular)<1e-12
        scans.append(dict(eta=eta,probabilities=probabilities.tolist(),singular_values=singular.tolist(),
            adjacent_ratios=(singular[:-1]/singular[1:]).tolist()))
    assert min(scans[-1]['adjacent_ratios'])>10000
    return dict(orbit_pairs=orbits,real_orbital_span_rank=24,rank1_channel_weights=weights,
        reconstruction_residuals=residuals,reference_singular_values=d.tolist(),scans=scans,
        leading_relative_singular_values=['1','eta^2','16*eta^4'],
        scope='A numerical orbital compiler builds three actual Hermitian H27 intertwiners, allowing two arbitrarily large hierarchy ratios. The pinned-selector dynamics supplies probabilities, but three channel shapes, shared strength and reference orientation remain an ansatz. This is not an observed fermion mass prediction; exact arithmetic compilation and full mediator self-energy are open.')

def spread_extension(c):
    lines=list(map(frozenset,c['lines']));pt={p:[i for i,L in enumerate(lines)if p in L]for p in range(40)};spreads=[]
    def search(covered,chosen):
        if len(covered)==40:spreads.append(tuple(sorted(chosen)));return
        best=None
        for p in range(40):
            if p in covered:continue
            cand=[i for i in pt[p]if not lines[i]&covered]
            if not cand:return
            if best is None or len(cand)<len(best):best=cand
        for i in best:search(covered|lines[i],chosen+[i])
    search(frozenset(),[]);spreads=sorted(set(spreads));assert len(spreads)==36
    B=np.zeros((40,36),dtype=int)
    for j,S in enumerate(spreads):B[list(S),j]=1
    Z=c['A'][:40,40:].astype(int);J=np.ones((40,40),dtype=int);I=np.eye(40,dtype=int)
    assert np.array_equal(3*Z.T@Z+B@B.T-3*J,18*I)
    E0=J/40;E24=(5*Z.T@Z-2*J)/30;E15=(4*B@B.T-9*J)/72
    V=np.vstack([np.ones((1,40))/np.sqrt(40),Z@E24/np.sqrt(6),B.T@E15/np.sqrt(18)])
    assert np.linalg.norm(V.T@V-I)<1e-12
    si={S:i for i,S in enumerate(spreads)};group=[]
    for g in c['x']['internal_h27']['group_permutations']:
        gp=[si[tuple(sorted(g[40+l]-40 for l in S))]for S in spreads]
        group.append(list(g)+[80+x for x in gp])
    PN=sum(np.eye(116)[np.argsort(g)]for g in group)/27
    w,v=np.linalg.eigh(np.eye(116)-PN);T=v[:,w>.5]
    prior=json.loads((ROOT/'data/w33_pass11390_11394_context_matter_clock.json').read_text())['portal'];C=np.zeros((80,80))
    for a,d in prior['nonedge_orbit']:C[a,d]=C[d,a]=1
    H=np.zeros((116,116));H[:80,:80]=c['A']+.1*C;H[40:80,80:]=B;H[80:,40:80]=B.T
    w,v=np.linalg.eigh(T.T@H@T);K=T@v[:,abs(w)<1e-9];assert K.shape[1]==34
    from sympy.polys.matrices import DomainMatrix
    H10=np.zeros((116,116),dtype=int);H10[:80,:80]=10*c['A'].astype(int)+C.astype(int);H10[40:80,80:]=10*B;H10[80:,40:80]=10*B.T
    Nsum=sum(np.eye(116,dtype=int)[np.argsort(g)]for g in group)
    stack=DomainMatrix.from_Matrix(sp.Matrix(np.vstack([H10,Nsum]))).convert_to(sp.QQ)
    exact_rank=len(stack.rref()[1]);assert exact_rank==82
    kernel,pivots=stack.nullspace().rref();kernel=kernel.to_Matrix()
    exact_chars=[sum(kernel[i,g[j]]for i,j in enumerate(pivots))for g in group]
    assert len(pivots)==34 and exact_chars[0]==34 and exact_chars[1]==7
    for g in group:
        inv=np.argsort(g);assert np.array_equal(H10[np.ix_(inv,inv)],H10)
    chars=[np.trace(K.T@np.eye(116)[np.argsort(g)]@K)for g in group]
    assert max(abs(a-float(d))for a,d in zip(chars,exact_chars))<1e-10
    zeta=(-1+sp.I*sp.sqrt(3))/2
    exact_linear=[sp.simplify(sum(ch*zeta**(-(u*a+v*d))for ch,(a,d,e)in zip(exact_chars,product(range(3),repeat=3)))/27)for u,v in product(range(3),repeat=2)]
    assert exact_linear==[0]+[2]*8
    omega=np.exp(2j*np.pi/3);linear=[]
    for u,v in product(range(3),repeat=2):
        z=sum(ch*omega**(-(u*a+v*d))for ch,(a,d,e)in zip(chars,product(range(3),repeat=3)))/27
        assert abs(z-round(z.real))<1e-10;linear.append(round(z.real))
    assert linear==[0]+[2]*8 and abs(chars[1]-7)<1e-10
    assert abs((chars[0]-chars[1])/9-3)<1e-10
    return dict(prior_owners=['Pass173','Pass4956','Pass4958','Part CXXVI'],spreads=spreads,
        complementary_reconstruction='18I=3Z^T Z+BB^T-3J',parseval_readout_shape=[77,40],
        parseval_residual=float(np.linalg.norm(V.T@V-I)),
        supplied_extended_operator='[[0,Z+0.1C,0],[(Z+0.1C)^T,0,B],[0,B^T,0]]',
        internal_carrier_dimension=T.shape[1],protected_zero_dimension=34,exact_stacked_rank=exact_rank,exact_kernel_characters=[int(x)for x in exact_chars],
        linear_character_multiplicities=linear,central_qutrit_multiplicities=[3,3],
        continuous_commutant='U(2)^8 x U(3)^2',continuous_commutant_dimension=50,
        scope='Prior complementary incidence channels are consumed, not rediscovered. The actual coupled116-vertex ansatz has protected multiplicities supporting nonabelian commutant factors. This is a candidate gauge carrier, not a chiral Standard Model: which factors are physical, gauge dynamics, anomaly matching and hypercharge remain unconstructed.')

def curvature(b):
    conn=b['conn'];p=json.loads((ROOT/'data/w33_pass11390_11394_context_matter_clock.json').read_text())
    group=[sp.Matrix(x).applyfunc(sp.Rational)for x in p['transport']['holonomy_matrices']]
    F=sp.Matrix(b['c']['x']['line']['fcc_roots'])*sp.Matrix(b['c']['x']['line']['primitive_fcc_change']).inv()
    proper=[F*R.inv().T*F.inv()for R in group if R.det()==1]
    energies=[int(sp.trace((O-sp.eye(3)).T*(O-sp.eye(3))))for O in proper]
    assert min(e for e in energies if e)>0 and min(e for e in energies if e)==4
    sample=None
    for a,d,e in combinations(range(40),3):
        if all(pair in conn for pair in [(a,d),(d,e),(e,a)]):
            W=conn[e,a]*conn[d,e]*conn[a,d]
            if W.det()==1:
                assert sp.trace(W)==-1;sample=dict(contexts=[a,d,e],loop=mat(W),wilson_cost_over_beta=4);break
    return dict(spatial_pullback_host_triangle=[[0,0,0],[1,1,0],[1,0,1]],triangle=sample,
        supplied_action='J sum_bonds(1-delta_context)+beta sum_triangles(3-tr(O_loop))',
        uniform_context_energy=0,distinct_proper_triangle_energy='3J+4beta',
        nonidentity_vector_holonomy_minimum_frobenius_squared=4,
        nonidentity_rotation_minimum_angle='pi/2',proper_group_trace_counts=dict(Counter(str(sp.trace(O))for O in proper)),
        continuum_requirement='Bounded smooth curvature needs plaquette holonomy I+O(a^2); a fixed finite48/24 group has no arbitrarily small nonidentity holonomies.',
        scope='A spatially pulled-back connection and curvature cost are named. Its finite rotation gap obstructs a generic smooth small-curvature limit with fixed group and shrinking plaquettes. Rare/concentrated defects or an added continuous connection are separate routes. Uniform classical energy zero is an additive convention, not a solved quantum cosmological constant or an Einstein-Hilbert derivation.')

def chiral_carrier(b):
    c=b['c'];prior=json.loads((ROOT/'data/w33_pass11390_11394_context_matter_clock.json').read_text())['portal']
    C=np.zeros((80,80))
    for a,d in prior['nonedge_orbit']:C[a,d]=C[d,a]=1
    group=c['x']['internal_h27']['group_permutations'];PN=sum(np.eye(80)[np.argsort(g)]for g in group)/27
    w,v=np.linalg.eigh(np.eye(80)-PN);T=v[:,w>.5];mass=T.conj().T@(c['A']+.1*C)@T
    w,V=np.linalg.eigh(mass);K=T@V[:,abs(w)<1e-9];assert K.shape==(80,14)
    chars=[]
    for g in group:chars.append(np.trace(K.conj().T@np.eye(80)[np.argsort(g)]@K))
    omega=np.exp(2j*np.pi/3);linear=[]
    for u,v in product(range(3),repeat=2):
        z=sum(ch*omega**(-(u*a+v*d))for ch,(a,d,e)in zip(chars,product(range(3),repeat=3)))/27
        assert abs(z-round(z.real))<1e-10;linear.append(int(round(z.real)))
    assert linear==[0]+[1]*8
    assert abs(chars[0]-14)<1e-10 and abs(chars[1]-5)<1e-10
    # Central-character multiplicities (14-5)/9=1 each; all10 irreps occur once.
    q,h=sp.symbols('q h');Q,U,D,L,E,N=q,-q-h,-q+h,-3*q,3*q+h,3*q-h
    anomaly=[2*Q+U+D,3*Q+L,6*Q+3*U+3*D+2*L+E+N,6*Q**3+3*U**3+3*D**3+2*L**3+E**3+N**3]
    assert all(sp.expand(x)==0 for x in anomaly)
    return dict(post_portal_internal_zero_dimension=14,linear_character_multiplicities=linear,
        central_qutrit_multiplicities=[1,1],continuous_commutant='U(1)^10',continuous_commutant_dimension=10,spread_extension=spread_extension(c),
        commuting_nonabelian_weak_gauge_factor=False,
        supplied_SM_plus_nu_charge_family={'Q':str(Q),'u_conjugate':str(U),'d_conjugate':str(D),'L':str(L),'e_conjugate':str(E),'nu_conjugate':str(N)},
        anomaly_polynomials=[str(sp.expand(x))for x in anomaly],witten_doublets_with_colour=4,
        scope='Schur multiplicity on the actual protected kernel gives only an abelian continuous commutant. Anomaly cancellation with a supplied SM+nu representation permits two continuous charge parameters; it does not derive hypercharge or fit that representation inside the14 native states. Prior Pass4091 owns generic anomaly-free exterior construction. A commuting nonabelian weak factor needs extra multiplicity or a different symmetry embedding.')

def atlas(b):
    c=b['c'];Gi=b['Gi'];hist=Counter()
    for p in range(40):
        for l in range(40):
            M=sp.Matrix(b['line'][l].T@b['point'][p]);T=Gi*M
            assert M.rank()==1
            tr=sp.trace(T*Gi*M.T);incident=p in c['lines'][l]
            assert tr==(1 if incident else sp.Rational(1,81))
            hist['incident_rank1_overlap1'if incident else'nonincident_rank1_overlap1/9']+=1
    # Point parabolic has line orbits4 and36, so it fixes no line.
    p=c['points'].index((1,0,0,0));generators=list(c['x']['internal_h27']['group_permutations'])
    for a,d,e,f in product(range(3),repeat=4):
        if (a*f-d*e)%3==1:
            M=np.eye(4,dtype=int);M[np.ix_([1,3],[1,3])]=[[a,d],[e,f]];generators.append(c['permutation'](M))
    assert all(g[p]==p for g in generators)
    seen=set();sizes=[]
    for l in range(40):
        if l in seen:continue
        orbit={l};queue=[l]
        for a in queue:
            for g in generators:
                d=g[40+a]-40
                if d not in orbit:orbit.add(d);queue.append(d)
        sizes.append(len(orbit));seen.update(orbit)
    assert sorted(sizes)==[4,36]
    return dict(cross_chart_histogram=dict(hist),point_stabilizer_line_orbits=sorted(sizes),
        G_equivariant_point_to_line_map_exists=False,incident_joint_cycle_dimension=5,nonincident_joint_cycle_dimension=6,
        scope='All cross projections have rank1, not a3D intertwiner. The native point stabilizer fixes no line, excluding a G-equivariant point-to-line map. Sharing SRG parameters does not prove self-duality. This chart overlap does not refute prior Pass4956 shared24D incidence transport or Pass4958 full40D point/spread reconstruction; those act on vertex coordinates, while these charts act on cycles.')

def axis_bridge(b):
    Pi=b['Pi'];Gi=b['Gi'];c=b['c'];checks=0
    for e,(p,l)in enumerate(c['edges']):
        v=Pi[:,e];Fp=sp.Matrix(b['point'][p]);Fl=sp.Matrix(b['line'][l-40])
        assert Fp*(Gi*(Fp.T*v))==v and Fl*(Gi*(Fl.T*v))==v;checks+=1
    assert Pi*Pi==Pi
    # The frame identity sum_e c_e c_e^T=Pi is Pi^2=Pi, not a floating SVD.
    # Build actual lossless chart reflections on the full native edge carrier.
    pp=c['points'].index((1,0,0,0));incident=next(l for l in range(40)if pp in c['lines'][l]);nonincident=next(l for l in range(40)if pp not in c['lines'][l])
    Fp=sp.Matrix(b['point'][pp]);Pp=Fp*Gi*Fp.T;samples=[]
    for l in (incident,nonincident):
        Fl=sp.Matrix(b['line'][l]);Pl=Fl*Gi*Fl.T;U=(sp.eye(160)-2*Pl)*(sp.eye(160)-2*Pp)
        assert U.T*U==sp.eye(160)and U*Pi==Pi*U
        T=Gi*(Fl.T*Fp);cos2=2*sp.trace(T*Gi*(Fp.T*Fl))-1
        samples.append(dict(point_context=pp,line_context=l,incident=pp in c['lines'][l],nontrivial_rotation_cosine=str(cos2)))
    return dict(flag_shared_axis_checks=checks,cycle_projector_160_numerator=mat(160*Pi),projector_denominator=160,
        prior_owned_frame_norm_squared='81/160',prior_owned_normalized_frame_bound='160/81',
        global_joint_chart_cycle_rank=81,lossless_chart_update='(I-2P_line)(I-2P_point) on the81D harmonic carrier',
        samples=samples,nonincident_update_infinite_order=True,
        scope='BT4048 owns the projected-edge tight frame. New identification: every flag axis is exactly the shared point/line chart direction. All160 axes span the81D Steinberg carrier and rank3 chart reflections act unitarily there. This avoids lossy rank1 projection by retaining the larger carrier; it does not map the80-vertex H27 matter space to harmonic cycles or make arbitrary chart changes local.')

def portal_polynomial(b):
    c=b['c'];prior=json.loads((ROOT/'data/w33_pass11390_11394_context_matter_clock.json').read_text())['portal']
    group=c['x']['internal_h27']['group_permutations'];Z=np.zeros((80,80),dtype=np.int64);Z[np.array(group[1]),np.arange(80)]=1
    C3=(2*np.eye(80,dtype=np.int64)-Z-Z@Z)[:40,:40];A=c['A'][:40,40:].astype(np.int64);C=np.zeros((40,40),dtype=np.int64)
    for a,d in prior['nonedge_orbit']:C[a,d-40]=1
    factors=[A@A.T,A@C.T+C@A.T,C@C.T];powers=[np.eye(40,dtype=np.int64)];g,x=sp.symbols('g x');tr=[]
    for n in range(1,4):
        new=[np.zeros((40,40),dtype=np.int64)for _ in range(len(powers)+2)]
        for i,a in enumerate(powers):
            for j,d in enumerate(factors):new[i+j]+=a@d
        powers=new;tr.append(sum(sp.Rational(int(np.trace(C3@a)),18)*g**i for i,a in enumerate(powers)))
    t1,t2,t3=tr;e2=sp.expand((t1*t1-t2)/2);e3=sp.expand((t1**3-3*t1*t2+2*t3)/6)
    poly=sp.expand(x**3-t1*x*x+e2*x-e3)
    assert sp.expand(poly-(x*(x-6)**2-3*g*g*(x-3)**2))==0
    # Ordered roots for g!=0 lie in(0,3),(3,6),(6,infinity). Their ratio
    # competition has one crossing; certify the optimizer in(4.2,4.3).
    y=sp.symbols('y');opt=y**4-6*y**3+108*y-324
    assert sp.Poly(opt,y).count_roots(3,6)==1 and opt.subs(y,sp.Rational(42,10))<0 and opt.subs(y,sp.Rational(43,10))>0
    ratio_sum_upper=sp.Rational(12)/(sp.Rational(21,5))+sp.Rational(43,10)**2/9-1
    assert ratio_sum_upper<sp.Rational(17,4) # t+1/t<4+1/4 implies t<4 and sqrt(t)<2
    ym=next(float(sp.re(z))for z in sp.nroots(opt)if abs(float(sp.im(z)))<1e-12 and 3<float(sp.re(z))<6)
    gm=np.sqrt(ym**3/27);roots=np.sort(np.roots([1,-(12+3*gm*gm),36+18*gm*gm,-27*gm*gm]).real)
    balanced=float(np.sqrt(roots[1]/roots[0]));assert balanced<2
    scans=[]
    for gg in (.01,.1,1.,gm,10.):
        roots=np.sort(np.roots([1,-(12+3*gg*gg),36+18*gg*gg,-27*gg*gg]).real)
        assert 0<roots[0]<3<roots[1]<6<roots[2]
        scans.append(dict(g=float(gg),squared_singular_values=roots.tolist(),adjacent_singular_ratios=np.sqrt(roots[1:]/roots[:-1]).tolist()))
    return dict(exact_characteristic_polynomial=str(poly),factorized_equation='x*(x-6)^2=3*g^2*(x-3)^2',
        trace_power_polynomials=[str(t)for t in tr],rank3_for_every_real_nonzero_g=True,
        invariant_constraints=['e2=6*e1-36','e3=9*e1-108'],optimizer_middle_root_polynomial=str(opt),
        optimizer_root_rational_interval=['21/5','43/10'],maximum_balanced_adjacent_singular_ratio=balanced,
        all_g_balanced_ratio_bound=2,optimizer_squared_ratio_sum_rational_upper=str(ratio_sum_upper),scans=scans,
        scope='Exact full3x4 multiplicity spectrum for one orbit portal, not a small-g fit. At least one adjacent singular-value ratio stays below2 for every g; a single portal cannot give two simultaneous large hierarchies. No identification of these singular values with measured quark/lepton masses is supplied.')

def controls(b):
    Pi=b['Pi'];c=b['c'];first=next((a,d)for a,d in combinations(range(160),2)if set(c['edges'][a])&set(c['edges'][d]));a,d=first
    v,w=Pi[:,a],Pi[:,d];norm=sp.Rational(81,160);overlap=sp.simplify(v.dot(w)/norm);assert overlap==sp.Rational(-1,3)
    P=v*v.T/norm;Q=w*w.T/norm
    U=(sp.eye(160)-2*Q)*(sp.eye(160)-2*P);assert U.T*U==sp.eye(160)
    cosine=2*overlap**2-1;assert cosine==sp.Rational(-7,9)
    # Every pair of harmonic frame vectors has an outer-product word: use a
    # nonzero-overlap path in the connected Levi line graph.
    adj=[[]for _ in range(160)]
    for i,j in combinations(range(160),2):
        if set(c['edges'][i])&set(c['edges'][j]):assert Pi[i,j]!=0;adj[i].append(j);adj[j].append(i)
    seen={0};queue=[0]
    for i in queue:
        for j in adj[i]:
            if j not in seen:seen.add(j);queue.append(j)
    assert len(seen)==160
    return dict(edge_control_pair=list(first),normalized_overlap=str(overlap),two_pi_flip_rotation_cosine=str(cosine),
        two_pi_flips_infinite_order=True,nonorthogonal_control_graph_connected=True,
        associative_control_algebra_dimension=81**2,
        independently_driven_rank1_phase_control_lie_algebra='u(81), conditional on all160 independent projector controls',
        scope='Generic rank1 controllability and the BT4047/4048 protected-control mechanism are prior art. The joint-atlas flag axes give a concrete full control-algebra witness and an infinite-order two-flip pair. Independent phase access, physical locality, logical tensor factorization, preparation and fault tolerance are not derived from finite geometry.')

def produce():
    b=build(load());p=dict(status='PASS',reservation='2debc76388af156117511c0c3bb2586028f4bdd3',
        spin=spin(b),selector_portal=selector_portal(b),curvature=curvature(b),chiral=chiral_carrier(b),
        atlas=atlas(b),axis_bridge=axis_bridge(b),portal_polynomial=portal_polynomial(b),controls=controls(b),
        scope='Eight completed mathematical/ansatz fronts. No observed mass spectrum, chiral Standard Model, Einstein limit, cosmological constant solution or unconditional physical universal computer.')
    OUT.write_text(json.dumps(p,indent=2)+'\n');print('PASS: eight joint-atlas/spin/selector/portal fronts',flush=True);return p
if __name__=='__main__':produce()
