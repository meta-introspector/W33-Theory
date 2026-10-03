#!/usr/bin/env python3
"""Five explicit follow-ups to Pass11389, with theorem/ansatz boundaries.

Prior owners: BT861 Steinberg carrier; Pass11389 cover and H27 block;
Pass4084 overlap fermions; Pass3325/3348 Szegedy walks. Generic lattice
fermions, Potts ordering and quantum walks are not claimed as new.
"""
from __future__ import annotations
import json
from collections import Counter
from itertools import combinations, product
from pathlib import Path
import numpy as np
import sympy as sp
from w33_pass11389_parabolic_spatial_cover import geometry, bloch

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/w33_pass11390_11394_context_matter_clock.json'

def load():
    x=json.loads((ROOT/'data/w33_pass11389_parabolic_spatial_cover.json').read_text())
    points,lines,edges,ei,D,permutation=geometry()
    H=sp.Matrix(x['line']['harmonic_voltage']).applyfunc(sp.Rational)
    S=np.array(x['internal_h27']['qutrit_seed_real'])+1j*np.array(x['internal_h27']['qutrit_seed_imag'])
    gamma=np.diag([1]*40+[-1]*40)
    w,v=np.linalg.eigh(S.conj().T@gamma@S)
    plus,minus=v[:,w>.5],v[:,w<-.5]
    A=4*np.eye(80)-np.array(D*D.T,dtype=float)
    B=plus.conj().T@S.conj().T@A@S@minus
    return dict(x=x,points=points,lines=lines,edges=edges,ei=ei,D=D,
                permutation=permutation,H=H,S=S,plus=plus,minus=minus,A=A,B=B)

def context_frames(c):
    p,lines,perm=c['points'],c['lines'],c['permutation']
    J=np.array([[0,0,1,0],[0,0,0,1],[-1,0,0,0],[0,-1,0,0]])
    gens=[perm((np.eye(4,dtype=int)+np.outer(v,np.array(v)@J))%3)for v in p]
    base=next(i for i,ll in enumerate(lines)if all(p[a][2:]==(0,0)for a in ll))
    reps={base:tuple(range(80))};queue=[base]
    for a in queue:
        for g in gens:
            b=g[40+a]-40
            if b not in reps:
                reps[b]=tuple(g[reps[a][i]]for i in range(80));queue.append(b)
    assert len(reps)==40
    scale=int(sp.ilcm(*[v.q for v in c['H']]))
    H=np.array(c['H']*scale,dtype=np.int64)
    frames=[]
    for a in range(40):
        g=reps[a];F=np.zeros_like(H)
        for i,(v,w)in enumerate(c['edges']):F[c['ei'][(g[v],g[w])]]=H[i]
        assert np.array_equal(np.array(c['D'],dtype=int)@F,np.zeros((80,3),dtype=int))
        frames.append(F)
    return frames,reps,scale

def selector(c):
    lines=c['lines'];adj=np.array([[int(i!=j and bool(set(a)&set(b)))for j,b in enumerate(lines)]for i,a in enumerate(lines)])
    assert np.all(adj.sum(axis=1)==12)
    assert np.array_equal(adj@adj,8*np.eye(40,dtype=int)-2*adj+4*np.ones((40,40),dtype=int))
    ev=np.linalg.eigvalsh(adj)
    assert np.max(abs(ev-np.array([-4]*15+[2]*24+[12])))<1e-12
    # H_selector=-h A, h>0, has one uniform ground state, not 40 selected vacua.
    assert np.linalg.norm(adj@np.ones(40)-12*np.ones(40))==0
    # Supplied ferromagnet: diagonal J(1-delta_ab) on each spatial bond.
    energies=[int(i!=j)for i,j in product(range(40),repeat=2)]
    assert Counter(energies)=={0:40,1:1560}
    # A one-dimensional finite open wall witness: distinct fixed endpoints.
    cost=np.full(40,999,dtype=int);cost[0]=0
    for _ in range(5):cost=np.min(cost[:,None]+(1-np.eye(40,dtype=int)),axis=0)
    assert cost[1]==1 and cost[0]==0
    return dict(native_context_adjacency=adj.tolist(),spectrum={'12':1,'2':24,'-4':15},
        one_cell_hamiltonian='-h*A_context, h>0',one_cell_ground='unique uniform superposition',
        one_cell_gap='10*h',ordered_hamiltonian='-h sum_x A_context(x)+J sum_<xy> (1-delta(context_x,context_y))',
        classical_connected_host_vacua=40,wall_energy_per_mismatched_bond='J',
        finite_five_bond_wall_minimum_units_J=1,
        unordered_wall_types={'intersecting':240,'disjoint':540},
        scope='J and h are supplied. Classical h=0 ordering and a finite wall are proved; no finite-h phase diagram, spontaneous selection in a finite symmetric system, or gravitational defect dynamics is inferred.')

def transport(c):
    frames,reps,scale=context_frames(c);G=sp.Matrix(frames[0].T@frames[0]);Gi=G.inv()
    transitions={};ranks=Counter();nonintegral=0;checks=0
    for i,j in combinations(range(40),2):
        M=sp.Matrix(frames[j].T@frames[i]);T=Gi*M
        incident=bool(set(c['lines'][i])&set(c['lines'][j]))
        if incident:
            assert M.rank()==1 and sp.trace(T*Gi*M.T)==sp.Rational(1,9)
            ranks['intersecting_rank1']+=1
        else:
            R=27*T
            assert R.T*G*R==G and abs(R.det())==1
            transitions[i,j]=R;transitions[j,i]=R.inv()
            nonintegral+=int(any(v.q!=1 for v in R))
            ranks['disjoint_rank3']+=1
        checks+=1
    hist=Counter();example=None
    for a,b,d in combinations(range(40),3):
        if all(pair in transitions for pair in [(a,b),(b,d),(d,a)]):
            W=transitions[d,a]*transitions[b,d]*transitions[a,b]
            assert W*W==sp.eye(3)
            key=f'trace={sp.trace(W)},det={W.det()}'
            hist[key]+=1
            if example is None:example=dict(contexts=[a,b,d],matrix=[[str(v)for v in W.row(i)]for i in range(3)])
    assert hist=={'trace=-1,det=1':2160,'trace=1,det=-1':1080}
    # Section-induced permutation transition has identically trivial holonomy.
    def compose(a,b):return tuple(a[b[i]]for i in range(80))
    def inverse(a):return tuple(np.argsort(a))
    a,b,d=example['contexts']
    ab=compose(reps[b],inverse(reps[a]));bd=compose(reps[d],inverse(reps[b]));da=compose(reps[a],inverse(reps[d]))
    assert compose(da,compose(bd,ab))==tuple(range(80))
    # Unlike intersecting pairs, disjoint polar maps preserve the full integer
    # cohomology lattice. Deck vectors transform contragrediently by R^(-T).
    assert nonintegral==0
    trees={0:sp.eye(3)};queue=[0]
    for a in queue:
        for b in range(40):
            if (a,b)in transitions and b not in trees:
                trees[b]=transitions[a,b]*trees[a];queue.append(b)
    loops={tuple(trees[b].inv()*R*trees[a]):trees[b].inv()*R*trees[a]
           for (a,b),R in transitions.items()}
    group={tuple(sp.eye(3)):sp.eye(3)};queue=[sp.eye(3)]
    for a in queue:
        for g in loops.values():
            b=a*g
            if tuple(b)not in group:group[tuple(b)]=b;queue.append(b)
    assert len(group)==48
    F=sp.Matrix(c['x']['line']['fcc_roots']).applyfunc(sp.Rational)
    U=sp.Matrix(c['x']['line']['primitive_fcc_change']).applyfunc(sp.Rational)
    cartesian=F*U.inv()
    for R in group.values():
        O=cartesian*R.inv().T*cartesian.inv()
        assert all(v in (-1,0,1)for v in O) and O.T*O==sp.eye(3)
    # Orientation double cover: determinant signs are not a vertex gauge,
    # since 1080 triangles have negative determinant. Its 80 states connect.
    seen={(0,1)};queue=[(0,1)]
    for a,sign in queue:
        for b in range(40):
            if (a,b)in transitions:
                target=(b,sign*int(transitions[a,b].det()))
                if target not in seen:seen.add(target);queue.append(target)
    assert len(seen)==80
    # Actual PSp action on the family of selected harmonic spaces. Its oriented
    # lift is transitive even though no single rank3 kernel is G-invariant.
    J=np.array([[0,0,1,0],[0,0,0,1],[-1,0,0,0],[0,-1,0,0]])
    lifted=[]
    for v in c['points']:
        g=c['permutation']((np.eye(4,dtype=int)+np.outer(v,np.array(v)@J))%3)
        lift=[];maps={}
        for a in range(40):
            b=g[40+a]-40;F=np.zeros_like(frames[a])
            for i,(u,w)in enumerate(c['edges']):F[c['ei'][(g[u],g[w])]]=frames[a][i]
            R=Gi*sp.Matrix(frames[b].T@F)
            assert all(z.q==1 for z in R) and R.T*G*R==G
            assert np.array_equal(frames[b]@np.array(R,dtype=np.int64),F)
            maps[a]=R;flip=int(R.det()==-1);lift.extend([2*b+flip,2*b+(1^flip)])
        a,b=0,13;ga,gb=g[40+a]-40,g[40+b]-40
        assert maps[b]*transitions[a,b]==transitions[ga,gb]*maps[a]
        lifted.append(lift)
    orbit={0};queue=[0]
    for a in queue:
        for g in lifted:
            if g[a]not in orbit:orbit.add(g[a]);queue.append(g[a])
    assert len(orbit)==80
    return dict(all_pairs_checked=checks,overlap_ranks=dict(ranks),
        orthonormal_overlap_singular_values={'intersecting':['1/3','0','0'],'disjoint':['1/27']*3},
        exact_disjoint_transport='R_ji=27*(H_j^T H_j)^(-1)*H_j^T H_i',
        nonintegral_primitive_period_transports=nonintegral,
        disjoint_triangle_holonomy_histogram=dict(hist),exact_all_triangle_holonomy_squares_identity=True,
        example=example,section_permutation_holonomy='identity (telescoping)',
        holonomy_group_order=48,holonomy_group='full signed-permutation FCC point group O_h',
        holonomy_matrices=[[[str(v)for v in R.row(i)]for i in range(3)]for R in group.values()],
        integral_deck_transition='R^(-T)',orientation_double_cover_connected_states=80,
        orientation_preserving_holonomy_order=24,
        native_transvection_orientation_permutations=lifted,
        exact_native_harmonic_action_maps_checked=1600,
        native_PSp_oriented_context_orbit=80,
        scope='Disjoint-pair maps give an actual integral lattice bundle and finite O_h connection with nontrivial holonomy. They do not glue all intersecting charts, define a smooth spacetime connection or establish Einstein gravity. The orientation obstruction is removed by a connected double cover; a spin lift and physical spin dynamics remain unbuilt.')

def edge_orbit_blocks(c):
    gs=c['x']['internal_h27']['group_permutations'];seen=set();blocks=[];orbits=[]
    for i,(a,b)in enumerate(c['edges']):
        if i in seen:continue
        orb=sorted({c['ei'][(g[a],g[b])]for g in gs});seen.update(orb)
        M=np.zeros((80,80))
        for z in orb:
            aa,bb=c['edges'][z];M[aa,bb]=M[bb,aa]=1
        blocks.append(c['plus'].conj().T@c['S'].conj().T@M@c['S']@c['minus']);orbits.append(orb)
    return blocks,orbits

def portal(c):
    B=c['B'];U,s,Vh=np.linalg.svd(B,full_matrices=True);left=U[:,-1];right=Vh[2:].conj().T
    blocks,orbits=edge_orbit_blocks(c)
    residual=max(float(np.linalg.norm(left.conj()@v))for v in blocks)
    assert residual<1e-12
    # Stronger common-left-kernel witness: exact projector onto intersection of
    # kernels will be independently checked using integer group averaging in tests.
    samples=[]
    for seed in range(5):
        weights=np.random.default_rng(seed).normal(size=16)
        s2=np.linalg.svd(sum(w*v for w,v in zip(weights,blocks)),compute_uv=False)
        assert s2[-1]<1e-12;samples.append(s2.tolist())
    gs=c['x']['internal_h27']['group_permutations'];a,b=1,47
    Z=np.zeros((80,80),dtype=np.int64);Z[np.array(gs[1]),np.arange(80)]=1
    central3=2*np.eye(80,dtype=np.int64)-Z-Z@Z
    incidence=c['A'][:40,40:].astype(np.int64)
    # These are 18 times the exact common-left and native-right null projectors
    # across BOTH conjugate central qutrit sectors (ranks6 and12).
    P18=(6*np.eye(40,dtype=np.int64)-incidence@incidence.T)@central3[:40,:40]
    Q18=(6*np.eye(40,dtype=np.int64)-incidence.T@incidence)@central3[40:,40:]
    assert np.array_equal(P18@P18,18*P18) and np.trace(P18)==108
    assert np.array_equal(Q18@Q18,18*Q18) and np.trace(Q18)==216
    for orb in orbits:
        edgeblock=np.zeros((40,40),dtype=np.int64)
        for z in orb:
            aa,bb=c['edges'][z];edgeblock[aa,bb-40]=1
        assert not np.any(P18@edgeblock)
    assert (a,b)not in c['ei']
    orbit=sorted({(g[a],g[b])for g in gs});assert len(orbit)==27
    # Locate an actual three-edge Levi path to the nonincident endpoint.
    adjacency=[[]for _ in range(80)]
    for aa,bb in c['edges']:adjacency[aa].append(bb);adjacency[bb].append(aa)
    path=next([a,u,v,b]for u in adjacency[a]for v in adjacency[u]if b in adjacency[v])
    paths=sorted({tuple(g[z]for z in path)for g in gs});assert len(paths)==27
    M=np.zeros((80,80))
    for aa,bb in orbit:M[aa,bb]=M[bb,aa]=1
    C=c['plus'].conj().T@c['S'].conj().T@M@c['S']@c['minus']
    coefficient=float(np.linalg.norm(left.conj()@C@right))
    assert abs(coefficient-np.sqrt(3)/2)<1e-12
    portal_integer=M[:40,40:].astype(np.int64)
    assert np.array_equal(P18@portal_integer@Q18@portal_integer.T@P18,243*P18)
    # Two auxiliary bipartite modes per orbit path. Strong middle coupling M0,
    # weak endpoint y. Exact Schur complement at energy zero: -y^2/M0 * portal.
    y=.1;M0=np.sqrt(6);heavy=np.zeros((54,54));W=np.zeros((80,54))
    for j,(aa,u,v,bb)in enumerate(paths):
        heavy[2*j,2*j+1]=heavy[2*j+1,2*j]=M0
        W[aa,2*j]=y;W[bb,2*j+1]=y
    eff=-W@np.linalg.inv(heavy)@W.T
    assert np.linalg.norm(eff+y*y/M0*M)<1e-12
    assert np.max(abs(np.diag(eff)))<1e-12
    # Actual low-energy gaps; coefficient refers to small coupling, not exact
    # finite-coupling singular value. No arbitrary rank-one SVD perturbation used.
    scans=[]
    for y in (.2,.1,.05,.025):
        g=-y*y/M0;ss=np.linalg.svd(B+g*C,compute_uv=False)
        scans.append(dict(endpoint_coupling=y,heavy_coupling=M0,portal_coefficient=g,
                          effective_small_singular_value=float(ss[-1]),ratio_to_y_squared_over_M=float(ss[-1]/abs(g))))
    assert abs(scans[-1]['ratio_to_y_squared_over_M']-coefficient)<1e-6
    return dict(native_edge_orbit_count=16,common_left_kernel_numerical_residual=residual,
        exact_common_left_projector_numerator=P18.tolist(),exact_native_right_projector_numerator=Q18.tolist(),
        projector_denominator=18,common_left_projector_rank_both_conjugate_sectors=6,
        native_right_projector_rank_both_conjugate_sectors=12,
        exact_coefficient_identity='P_left * C * P_right * C^T * P_left = (3/4) P_left',
        arbitrary_native_edge_weight_rank_upper_bound=2,random_weight_singular_values=samples,
        nonedge_seed=[a,b],nonedge_orbit=orbit,virtual_path_seed=path,virtual_paths=paths,
        auxiliary_vertex_count=54,auxiliary_heavy_matrix='direct sum of [[0,M],[M,0]]',
        energy_zero_schur_portal='-y_left*y_right/M * orbit_incidence',
        leading_small_singular_value='sqrt(3)/2 * abs(y_left*y_right/M)',
        scans=scans,exact_schur_residual=float(np.linalg.norm(eff+y*y/M0*M)),
        scope='Native H27-equivariant incidence weights cannot lift this common kernel. An added 54-mode, three-edge-local mediator extension can lift it quadratically. The coefficient is computed, but y and M are inputs; the zero-energy Schur complement is an effective operator, not the exact finite-energy spectrum or a predicted observed mass.')

def paulis():
    return [np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]],complex),np.diag([1,-1]).astype(complex)]

def fermions(c):
    sig=paulis();alpha=[np.kron(sig[0],s)for s in sig];beta=np.kron(sig[2],np.eye(2))
    F=np.array(c['x']['line']['fcc_roots'],float);reciprocal=np.linalg.inv(F).T
    corners=[]
    for bits in product((0,1),repeat=3):
        charge=int(round(np.linalg.det(reciprocal)*np.prod([(-1)**b for b in bits])/abs(np.linalg.det(reciprocal))))
        corners.append(dict(bits=bits,weyl_charge=charge,wilson_mass_r1=2*sum(bits)))
    assert sum(x['weyl_charge']for x in corners)==0
    assert sum(x['wilson_mass_r1']==0 for x in corners)==1
    # CAR are defined on a Fock space; a graph eigenvector alone is not a fermion.
    Z=sig[2];lower=np.array([[0,1],[0,0]],complex);ann=[]
    for j in range(4):
        z=np.array([[1]],complex)
        for k in range(4):z=np.kron(z,Z if k<j else lower if k==j else np.eye(2))
        ann.append(z)
    residual=max(np.linalg.norm(a@b.conj().T+b.conj().T@a-(np.eye(16)if i==j else 0))for i,a in enumerate(ann)for j,b in enumerate(ann))
    residual=max(residual,max(np.linalg.norm(a@b+b@a)for a in ann for b in ann))
    assert residual==0
    h=np.array([.13,-.07,.11]);d=reciprocal@np.sin(h);r=sum(1-np.cos(h));D=sum(v*a for v,a in zip(d,alpha))+r*beta
    assert np.linalg.norm(D@D-(d@d+r*r)*np.eye(4))<1e-12
    # Tensor with the ACTUAL native seven-dimensional multiplicity operator,
    # not an unrelated list of proposed particle masses.
    a7=c['S'].conj().T@c['A']@c['S']
    internal_mass=a7+r*np.eye(7)
    H28=np.kron(sum(v*a for v,a in zip(d,alpha)),np.eye(7))+np.kron(beta,internal_mass)
    square_residual=float(np.linalg.norm(H28@H28-
        (np.eye(28)*(d@d)+np.kron(np.eye(4),internal_mass@internal_mass))))
    assert square_residual<1e-12 and np.linalg.norm(H28-H28.conj().T)<1e-12
    return dict(fcc_primitive_matrix=F.tolist(),naive_weyl_symbol='sigma dot (F^(-T) sin(theta))',
        corner_audit=corners,total_weyl_charge=0,
        wilson_dirac_symbol='alpha dot (F^(-T)sin(theta))+beta*r*sum_j(1-cos(theta_j))',
        wilson_light_external_dirac_species_per_massless_internal_carrier=1,
        spinor_dimension=4,qutrit_fiber_dimension=3,multiplicity_dimension=7,
        native_spinor_operator_dimension_per_central_sector=84,
        native_spinor_operator='I3 tensor [alpha dot d tensor I7 + beta tensor (S^dagger A_native S + Wilson_mass I7)]',
        native_spinor_square_residual=square_residual,
        native_zero_internal_qutrit_copies=3,after_generic_rank_lift_zero_internal_qutrit_copies=1,
        fock_modes_checked=4,fock_dimension=16,car_residual=float(residual),
        vectorlike_u1_cubic_anomaly_per_qutrit='3*(q^3-q^3)=0',vectorlike_u1_gravitational_anomaly_per_qutrit='3*(q-q)=0',
        scope='CAR quantization, a supplied spinor factor and Wilson hopping construct vectorlike lattice fermions. H27 dimension3 is not a derived SU(3) gauge charge; the graph grading index is not Lorentz chirality. No Standard Model charges or chiral anomaly-free spectrum is assigned. Generic Wilson/overlap precedent is Pass4084.')

def clock_matrices(c,k,voltage=None):
    edges=c['edges'];voltage=np.asarray(voltage if voltage is not None else c['x']['line']['integer_voltage'],float)
    # Arc index 2e is tail a and 2e+1 tail b; reversal carries opposite voltage.
    d=np.zeros((80,320),complex);S=np.zeros((320,320),complex)
    for i,((a,b),v)in enumerate(zip(edges,voltage)):
        d[a,2*i]=d[b,2*i+1]=.5
        S[2*i,2*i+1]=np.exp(1j*np.dot(v,k));S[2*i+1,2*i]=np.exp(-1j*np.dot(v,k))
    C=2*d.conj().T@d-np.eye(320);U=S@C
    return U,S,d

def clock(c):
    samples=[]
    for k in ([0,0,0],[.11,-.23,.31]):
        U,S,d=clock_matrices(c,k)
        L=bloch(c['edges'],c['x']['line']['integer_voltage'],k)
        P=np.eye(80)-L/4
        assert np.linalg.norm(S@S-np.eye(320))<1e-12
        assert np.linalg.norm(U.conj().T@U-np.eye(320))<1e-12
        assert np.linalg.norm(d@S@d.conj().T-P)<1e-12
        w,v=np.linalg.eigh(P);res=0.
        for j in range(80):
            if abs(w[j])>1-1e-8:continue
            z=np.exp(1j*np.arccos(w[j]));aa=d.conj().T@v[:,j];bb=S@aa
            psi=aa-z*bb
            res=max(res,float(np.linalg.norm(U@psi-z*psi)))
        assert res<1e-12
        samples.append(dict(momentum=k,unitarity_residual=float(np.linalg.norm(U.conj().T@U-np.eye(320))),spectral_mapping_residual=res))
    E=np.array(sp.Matrix(c['x']['line']['fcc_harmonic_displacements']).applyfunc(sp.Rational),float)
    maxstep=sp.simplify(max(sum(sp.Rational(z)**2 for z in row)for row in c['x']['line']['fcc_harmonic_displacements']))
    # Real-space locality: each coin remains at tail; shift crosses one edge.
    U,S,d=clock_matrices(c,[0,0,0]);adj=np.zeros((80,80),int)
    for a,b in c['edges']:adj[a,b]=adj[b,a]=1
    tails=[a for a,b in c['edges']for a in (a,b)]
    distance=np.full((80,80),999,int);distance[adj==1]=1;np.fill_diagonal(distance,0)
    for v in range(80):distance=np.minimum(distance,distance[:,v,None]+distance[None,v,:])
    power=np.eye(320)
    for ticks in range(1,5):
        power=U@power
        rows,cols=np.where(abs(power)>1e-12)
        assert all(distance[tails[a],tails[b]]<=ticks for a,b in zip(rows,cols))
    return dict(arc_dimension=320,coin='2*d^dagger*d-I, d(vertex, outgoing arc)=1/2',
        update='U(k)=S_voltage(k)*coin',spectral_mapping='cos(omega)=1-lambda_Levi(k)/4',
        acoustic_omega_squared_coefficient='27/6400',maximum_edge_displacement_squared=str(maxstep),
        strict_support_cone='after n ticks: at most n native edges, Euclidean displacement <= n*max_edge_length',
        finite_support_ticks_checked=4,bipartite_quasienergy_partner='Gamma-uniform state at omega=pi',
        samples=samples,scope='A local unitary discrete update replaces the supplied continuous wave evolution. The tick duration in seconds is still an input. The native internal flat bands remain flat; universal gates and the added FCC fiber hopping require additional controls. Szegedy/Grover spectral mapping is classical prior art, not a new theorem.')

def produce():
    c=load()
    packet=dict(status='PASS',reservation='86ebd9e685b5c9a2d1efa27b45bc80d48bee1d41',
        prior_owners=['Pass11389: native spatial cover and qutrit multiplicity block','BT861: irreducible Steinberg H1','Pass4084: overlap lattice fermions','Pass3325/3348: Szegedy walk compilation'],
        selector=selector(c),fermions=fermions(c),portal=portal(c),transport=transport(c),clock=clock(c),
        physical_boundary='Five mathematical constructions/witnesses. No complete TOE, observed mass spectrum, chiral Standard Model, dynamical Einstein gravity or derived dimensional constants.')
    OUT.write_text(json.dumps(packet,indent=2)+'\n');print(json.dumps({k:packet[k]for k in ['status','physical_boundary']}));return packet

if __name__=='__main__':produce()
