"""Constructed maps for five physics targets; supplied model, not a TOE.

Prior owners: 11271/11293 (SM/Higgs matching), 11476/11481/11493
(native frame, cover, Ward and decoding), and classical exact root isolation,
Schur complements, Gaussian determinants and graph Green functions.
"""
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/w33_pass11506_11510_five_physics_targets.json'
SOURCES = ['data/w33_pass11384_native_inputs.json',
           'data/w33_pass11271_canonical_sm_higgs_bridge.json',
           'data/w33_pass11476_11480_cubic_geometry_correlated.json',
           'data/w33_pass11481_11485_fibers_metric_yukawa_noisy.json',
           'data/w33_pass11493_11497_exact_matching_inference.json']


def read(p):
    return json.loads((ROOT/p).read_text())


def enc(a):
    a = np.asarray(a)
    return dict(real=a.real.tolist(), imag=a.imag.tolist())


def dec(a):
    return np.array(a['real'])+1j*np.array(a['imag'])


def tensors():
    src = read(SOURCES[0]); d = np.zeros((27,27,27), int)
    B = np.zeros((78,27,27), int)
    for i,j,k,v in src['signed_triads']:
        for ix in itertools.permutations((i,j,k)): d[ix] = v
    for a, rows in enumerate(src['E6_sparse_basis']):
        for i,j,v in rows: B[a,i,j] = v
    return B,d,src['signed_triads']


def balanced_vacuum():
    """Exact first jets on all81 complex coordinates and rational root isolation."""
    import w33_20261001_degree18_phase_completion as D
    import w33_pass11438_finite_native_model as F
    B,d,tri = tensors(); perms = list(itertools.permutations(range(3)))
    sign = lambda p: (-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
    E = sp.zeros(27,3)
    for a,k in enumerate(tri[0][:3]): E[k,a] = tri[0][3] if a==2 else 1
    P = E*sp.diag(1,1,sp.I)
    X = [sp.zeros(3) for _ in range(27)]
    Fc = np.empty((3,3,3), object); Fc.fill(sp.S.Zero)
    for i,j,k,v in tri:
        for a,b,c in itertools.permutations((i,j,k)):
            for u,w in itertools.product(range(3),repeat=2): X[a][u,w] += v*P[b,u]*P[c,w]
            for u,w,z in itertools.product(range(3),repeat=3): Fc[u,w,z] += v*P[a,u]*P[b,w]*P[c,z]
    def adj(b,c):
        ans = sp.zeros(3)
        for p,q in itertools.product(perms,repeat=2):
            i,k,m=p; j,l,n=q
            ans[i,j] += sign(p)*sign(q)*b[k,l]*c[m,n]
        return ans
    gx = [sp.zeros(3) for _ in range(27)]; raw6=0
    for i,j,k,v in tri:
        for a,b,c in [(i,j,k),(j,i,k),(k,i,j)]: gx[a] += 6*v*adj(X[b],X[c])
        raw6 += 6*v*sum(X[i][a,b]*adj(X[j],X[k])[a,b] for a,b in itertools.product(range(3),repeat=2))
    g6=sp.zeros(27,3)
    for i,j,k,v in tri:
        for a,b,c in itertools.permutations((i,j,k)):
            for u,w in itertools.product(range(3),repeat=2): g6[b,u] += sp.Rational(9,2)*v*gx[a][u,w]*P[c,w]
    a6=sp.Rational(9,4)*raw6
    with F.native_context():
        idx,sg=D.I.aronhold_terms()
        keys,z,S,T,_=D.cubic_circuit()
    gf=np.empty((3,3,3),object);gf.fill(sp.S.Zero);raw12=0
    for row,v in zip(idx,sg):
        vals=[Fc[tuple(k)] for k in row];raw12 += int(v)*sp.prod(vals)
        for j in range(4): gf[tuple(row[j])] += int(v)*sp.prod(vals[:j]+vals[j+1:])
    def pull(gf):
        ans=sp.zeros(27,3)
        for i,j,k,v in tri:
            for aa,bb,cc in itertools.permutations((i,j,k)):
                for u,w,z in itertools.product(range(3),repeat=3): ans[aa,u]+=3*v*gf[u,w,z]*P[bb,w]*P[cc,z]
        return ans
    a12=(152*a6*a6-3645*raw12)/32
    g12=(304*a6*g6-3645*pull(gf))/32
    subs=dict(zip(z,[Fc[k] for k in keys]))
    raw18=T.subs(subs);gt=np.empty((3,3,3),object);gt.fill(sp.S.Zero)
    for k,var in zip(keys,z):
        ps=set(itertools.permutations(k));val=sp.diff(T,var).subs(subs)/len(ps)
        for pp in ps:gt[pp]=val
    a18=(19683*raw18-sp.Rational(48384,5)*a6**3+sp.Rational(13824,5)*a6*a12)/1492992
    g18=(19683*pull(gt)-sp.Rational(145152,5)*a6*a6*g6+sp.Rational(13824,5)*(a12*g6+a6*g12))/1492992
    assert a6==-27 and a12==729 and a18==0
    assert g6==-54*P.conjugate() and g12==2*a6*g6 and g18==sp.zeros(27,3)
    # Moment zero without the arbitrary numerical orthonormalization of generators.
    assert all((P.conjugate().T*sp.Matrix(b)*P).trace()==0 for b in B)
    assert P.conjugate().T*P==sp.eye(3)
    r=sp.symbols('r');rho=sp.Rational(1,10**6)
    poly=sp.Poly(8*27**4*r**10-27*r+12*27**6*r**16+24*rho,r)
    intervals=sp.polys.polytools.intervals(poly,eps=sp.Rational(1,10**30))
    positive=[(ab,m) for ab,m in intervals if ab[0]>0]
    assert len(positive)==2 and all(m==1 for _,m in positive)
    lo,hi=positive[-1][0];assert poly.count_roots(lo,hi)==1
    derivative=poly.diff();assert derivative.count_roots(lo,hi)==0 and derivative.eval(lo)>0
    c=float(sp.sqrt((lo+hi)/2));pf=np.asarray(P,complex)*c;x=F.pack(pf.ravel(),pf.ravel())
    with F.native_context(): value,g=F.fun(x)
    assert np.linalg.norm(g)<1e-8
    prior=read(SOURCES[-1])['stationary']['history'][-1]['value']
    return dict(status='PASS',exact_template=[[str(v) for v in row] for row in P.tolist()],
                exact_first_jets=dict(I6=str(a6),I12=str(a12),I18=str(a18),
                    dI6='-54 conjugate(P)',dI12='2 I6 dI6',dI18='0 on all81 coordinates'),
                radial_polynomial=str(poly.as_expr()),positive_root_intervals=[list(map(str,ab)) for ab,_ in positive],
                selected_radius_squared_interval=[str(lo),str(hi)],selected_root_count=1,
                selected_radial_derivative_positive=True,radius=c,native_coordinates=x.tolist(),
                native_value=value,numerical_full_gradient_norm=float(np.linalg.norm(g)),
                prior_CP_candidate_value=prior,energy_gap_above_CP_candidate=value-prior,
                exact_stationarity_proof='Both fields cP have zero moments; cross C is scalar so commutator and its first derivative vanish. Homogeneity gives the three full holomorphic jets. The full324-gradient is c^3 times the stored radial polynomial times the real template, in both fields. Rational root isolation certifies a full stationary point and positive radial curvature.',
                scope='Alternate balanced stationary point of the SAME supplied native model; exact first jets and isolating interval, not a certificate of the lower-energy CP-breaking candidate, its transverse stability or a global vacuum. Gauge zeros are handled by exact moment identities; full gauge-fixed CP slice interval certification remains open.')


def sm_blocks():
    B,d,_=tensors(); sm=read(SOURCES[1]);y=list(map(sp.Rational,sm['hypercharge_diagonal']))
    heavy=sm['exotic_mass_support'];light=[i for i in range(27) if i not in heavy]
    # Identify the actual SU2 Cartan using its action on the light lepton doublet.
    k=next(k for k in sm['unbroken_root_basis_indices'] if np.any(np.diag(B[k]@B[k].T-B[k].T@B[k])[[15,16]]))
    t3=list(map(lambda z:sp.Rational(int(z),2),np.diag(B[k]@B[k].T-B[k].T@B[k])))
    hu=next(i for i in heavy if y[i]==sp.Rational(1,2) and y[i]+t3[i]==0)
    hd=next(i for i in heavy if y[i]==-sp.Rational(1,2) and y[i]+t3[i]==0)
    subsets=dict(up_Q=[i for i in light if y[i]+t3[i]==sp.Rational(2,3)],
                 down_Q=[i for i in light if y[i]==sp.Rational(1,6) and y[i]+t3[i]==-sp.Rational(1,3)],
                 Uc=[i for i in light if y[i]==-sp.Rational(2,3)],
                 Dc=[i for i in light if y[i]==sp.Rational(1,3)],
                 charged_L=[i for i in light if y[i]==-sp.Rational(1,2) and y[i]+t3[i]==-1],
                 Ec=[i for i in light if y[i]==1])
    blocks={name:sp.Matrix(d[:,:,h]).extract(subsets[a],subsets[b]) for name,a,b,h in
            [('up','up_Q','Uc',hu),('down','down_Q','Dc',hd),('lepton','charged_L','Ec',hd)]}
    assert blocks['up'].shape==blocks['down'].shape==(3,3)
    assert all(m*m.T==sp.eye(m.rows) for m in blocks.values())
    Fu=sp.diag(1,2,3);Fd=sp.Matrix([[1,1,sp.I],[1,2,1],[sp.I,1,3]])
    full=sp.kronecker_product(sp.Matrix(d[:,:,0]),sp.eye(3))
    perturb=sp.kronecker_product(sp.Matrix(d[:,:,hu]),Fu)+sp.kronecker_product(sp.Matrix(d[:,:,hd]),Fd)
    hh=[3*i+j for i in heavy for j in range(3)];ll=[3*i+j for i in light for j in range(3)]
    H0=full.extract(hh,hh);Ehh=perturb.extract(hh,hh);Elh=perturb.extract(ll,hh);Ell=perturb.extract(ll,ll)
    correction=-(Elh*H0.inv()*Elh.T)
    family_comm=Fu*Fu.conjugate().T*Fd*Fd.conjugate().T-Fd*Fd.conjugate().T*Fu*Fu.conjugate().T
    family_CP=sp.expand(sp.trace(family_comm**3))
    assert family_CP==-21600*sp.I
    assert all(y[i]+t3[i]+y[j]+t3[j]==0 for i,j in np.argwhere(np.array(d[:,:,hu]+d[:,:,hd])!=0))
    return dict(status='PASS',SU2_root=k,T3=list(map(str,t3)),neutral_Hu_index=hu,neutral_Hd_index=hd,
                heavy_E6_indices=heavy,light_E6_indices=light,SM_indices=subsets,
                exact_color_blocks={k:[[str(x) for x in row] for row in m.tolist()] for k,m in blocks.items()},
                Fu=[[str(x) for x in row] for row in Fu.tolist()],Fd=[[str(x) for x in row] for row in Fd.tolist()],
                heavy_rank=H0.rank(),light_first_order_rank=Ell.rank(),
                heavy_EW_block_nonzero=not Ehh.is_zero_matrix,
                light_second_order_rank=correction.rank(),
                exact_family_quark_CP_cubic=str(family_CP),exact_color_repeated_quark_CP_cubic=str(3*family_CP),
                second_order_nonzero=[[i,j,str(correction[i,j])] for i,j in itertools.product(range(51),repeat=2) if correction[i,j]],
                exact_matching='M_light=t E_ll - t² E_lh H0^{-1} E_hl + O(t³); exact Schur map t E_ll-t² E_lh(H0+t E_hh)^{-1}E_hl whenever the heavy block is invertible.',
                down_lepton_relation='Same family Fd multiplies signed orthogonal color down block and signed scalar charged-lepton block; their family singular values agree at leading order.',
                scope='Actual canonical SM sectors and neutral electroweak doublet components instantiated in signed E6 cubic, with supplied independently matched bar6 family tensors. Heavy Schur matching is explicit. Classical E6 branching and Schur formula are prior art;11271 owns SM/exotic interface and11293 mixed-source mediation. Family matrices, EW scales, Higgs selection and loop thresholds are supplied/open; no measured masses or realistic down/lepton fit is claimed.')


def metric_stabilization():
    """Massive vertex bosons plus one graph-Dirac Grassmann field, no grounding."""
    from scipy.optimize import minimize
    old=read(SOURCES[2])['geometry'];h=np.array(old['harmonic_FCC_displacements']);edges=np.array(old['cover_edges'])
    B=np.zeros((160,80));B[np.arange(160),edges[:,0]]=1;B[np.arange(160),edges[:,1]]=-1
    S=[np.diag([1,-1,0]),np.diag([1,1,-2])]
    for i,j in [(0,1),(0,2),(1,2)]:
        s=np.zeros((3,3));s[i,j]=s[j,i]=1;S.append(s)
    S=np.array(S)
    def operator(z):
        G=expm(np.einsum('a,aij->ij',z,S));w=1/np.einsum('ei,ij,ej->e',h,G,h)
        return G,B.T@(w[:,None]*B),w
    def action(z):
        _,L,_=operator(z)
        return np.linalg.slogdet(L+4*np.eye(80))[1]-np.linalg.slogdet(L+np.eye(80))[1]
    base=action(np.zeros(5));G,L,w=operator(np.zeros(5));eig=np.linalg.eigvalsh(L);eig[abs(eig)<1e-9]=0
    inv4=np.linalg.inv(L+4*np.eye(80));inv1=np.linalg.inv(L+np.eye(80));D=inv4-inv1
    length=np.einsum('ei,ei->e',h,h);r=np.einsum('ei,aij,ej->ae',h,S,h)/length
    Li=np.array([B.T@((-w*t)[:,None]*B) for t in r]);grad=np.einsum('ij,aji->a',D,Li)
    Hess=np.zeros((5,5))
    for i,j in itertools.product(range(5),repeat=2):
        u=np.einsum('ei,ij,ej->e',h,(S[i]@S[j]+S[j]@S[i])/2,h)/length
        Lij=B.T@((w*(2*r[i]*r[j]-u))[:,None]*B)
        Hess[i,j]=np.trace(D@Lij-inv4@Li[i]@inv4@Li[j]+inv1@Li[i]@inv1@Li[j])
    fd=np.array([(action(np.eye(5)[i]*1e-3)+action(-np.eye(5)[i]*1e-3)-2*base)/1e-6 for i in range(5)])
    assert max(abs(fd-np.diag(Hess)))<1e-5 and min(np.linalg.eigvalsh(Hess))>0 and np.linalg.norm(grad)<1e-8
    # Graph Dirac has79 singular-value pairs and82 zero modes.
    sqrtB=np.sqrt(w)[:,None]*B;Dirac=np.block([[np.zeros((80,80)),sqrtB.T],[sqrtB,np.zeros((160,160))]])
    dirac_logdet=np.linalg.slogdet(np.eye(240)+1j*Dirac)[1]
    assert abs(dirac_logdet-np.linalg.slogdet(np.eye(80)+L)[1])<1e-9
    certified=certify_metric_shape(S)
    slope=sum(v/(4+v)-v/(1+v) for v in eig)
    # Exact monotonicity proof is spectral and works for every positive eigenvalue.
    scans=[]
    for t in [.1,.5,1.,2.,4.]:
        z=np.array([t,0,0,0,0.]);f=3/(np.exp(t)+np.exp(-t)+1)
        prediction=sum(np.log((4+f*v)/(1+f*v))-np.log((4+v)/(1+v)) for v in eig)
        delta=action(z)-base;assert delta>0 and abs(delta-prediction)<1e-9
        scans.append(dict(t=t,shift=delta,spectral_prediction=float(prediction)))
    return dict(status='PASS',bosons=2,complex_graph_Dirac_fields=1,boson_mass_squared=4,fermion_mass=1,
                graph_Dirac_dimension=240,Dirac_zero_modes=82,ungrounded_Laplacian=L.tolist(),
                shape_basis=S.tolist(),shape_gradient=grad.tolist(),analytic_shape_Hessian=Hess.tolist(),
                shape_Hessian_eigenvalues=np.linalg.eigvalsh(Hess).tolist(),FD_diagonal=fd.tolist(),
                metric_reference=G.tolist(),Dirac_logdet_error=float(abs(dirac_logdet-np.linalg.slogdet(np.eye(80)+L)[1])),
                rational_shape_certificate=certified,
                diagonal_scans=scans,scale_spectral_slope=float(slope),
                exact_diagonal_proof='For detG=1 diagonal, L(G)=f L(I), f=3/TrG in(0,1]. Each nonzero lambda contributes log((4+f lambda)/(1+f lambda)), whose derivative is -3 lambda/((4+f lambda)(1+f lambda))<0. AM-GM therefore proves the unique diagonal shape minimum at I, and a finite positive barrier as f->0.',
                exact_Dirac_determinant='D=[[0,B^T sqrt(w)],[sqrt(w) B,0]]; det(m I+iD)=m^(160-80) det(m² I+L). Its log phase is zero for m>0.',
                scope='Supplied two-real-boson and one complex graph-Dirac Gaussian inventory on the ACTUAL80-vertex/160-edge graph with prior11494 rational displacements. Massive ungrounded measure differs from11494 grounded massless scalar. Exact diagonal global stabilization and full five-shape strict local minimum via graph symmetries plus rational inverse-residual/Sylvester certificate. This graph exterior-form fermion is not an identified physical SM/spin Dirac operator; masses and multiplicities are supplied, radial scale force remains nonzero, and intrinsic Einstein dynamics/CC are unbuilt.')


def certify_metric_shape(S):
    """Rational enclosure of Hessian; numerical inverse proposal never proves it."""
    import mpmath as mp
    import networkx as nx
    src=read(SOURCES[-1])['metric'];hi=np.array(src['displacements_times80'],dtype=int);edges=src['edges']
    den=3**9;scale=2**180;n=80
    B=np.zeros((160,n),dtype=object)
    for e,(a,b) in enumerate(edges):B[e,a]=1;B[e,b]=-1
    S=np.array(S,dtype=int);w=[sp.Rational(6400,int(v@v)) for v in hi]
    r=[[sp.Rational(int(v@s@v),int(v@v)) for v in hi] for s in S]
    def integer_laplacian(weights):
        ans=np.zeros((n,n),dtype=object)
        for (a,b),v in zip(edges,weights):
            vv=v*den;assert vv.q==1;vv=int(vv)
            ans[a,a]+=vv;ans[b,b]+=vv;ans[a,b]-=vv;ans[b,a]-=vv
        return ans
    N=integer_laplacian(w);Ni=[integer_laplacian([-x*y for x,y in zip(w,ri)]) for ri in r];Nij={}
    for i,j in itertools.product(range(5),repeat=2):
        vals=[]
        for e,v in enumerate(hi):
            u=sp.Rational(int(v@(S[i]@S[j]+S[j]@S[i])@v),2*int(v@v))
            vals.append(w[e]*(2*r[i][e]*r[j][e]-u))
        Nij[i,j]=integer_laplacian(vals)
    inverses=[];errors=[]
    for m2 in [4,1]:
        A=N+m2*den*np.eye(n,dtype=object)
        with mp.workdps(75):
            proposal=mp.inverse(mp.matrix([[mp.mpf(int(v))/den for v in row] for row in A]))
            Q=np.array([[int(mp.nint(v*scale)) for v in row] for row in proposal.tolist()],dtype=object)
        residual=den*scale*np.eye(n,dtype=object)-A@Q
        err=sp.Rational(n*max(abs(int(v)) for v in residual.ravel()),m2*den*scale)
        assert err<sp.Rational(1,10**40)
        inverses.append(Q);errors.append(err)
    Qb,Qf=inverses;eb,ef=errors
    Pb=[Qb@a for a in Ni];Pf=[Qf@a for a in Ni]
    trace_product=lambda a,b:sum(int(a[i,j])*int(b[j,i]) for i,j in itertools.product(range(n),repeat=2))
    norm=lambda a:sp.Rational(max(sum(abs(int(x)) for x in row) for row in a),den)
    H=sp.zeros(5);bound=sp.S.Zero
    for i,j in itertools.product(range(5),repeat=2):
        H[i,j]=sp.Rational(trace_product(Qb-Qf,Nij[i,j]),scale*den)+sp.Rational(trace_product(Pf[i],Pf[j])-trace_product(Pb[i],Pb[j]),(scale*den)**2)
        error=n*((eb+ef)*norm(Nij[i,j])+(eb*(sp.Rational(1,2)+eb)+ef*(2+ef))*norm(Ni[i])*norm(Ni[j]))
        bound=max(bound,error)
    lower=H-5*bound*sp.eye(5)
    minors=[lower[:i,:i].det(method='domain-ge') for i in range(1,6)]
    assert all(x>0 for x in minors)
    # Exact stationarity: two even sign flips negate every off-diagonal
    # coordinate derivative; diagonal traceless derivatives vanish edgewise.
    permutations=[];flips=[np.diag([1,-1,-1]),np.diag([-1,1,-1])]
    key=lambda v:tuple((v if v[0]>0 else -v).tolist())
    for flip in flips:
        g=nx.Graph();target=nx.Graph()
        for (a,b),v in zip(edges,hi):g.add_edge(a,b,color=key(v));target.add_edge(a,b,color=key(flip@v))
        matcher=nx.algorithms.isomorphism.GraphMatcher(g,target,edge_match=lambda a,b:a['color']==b['color'])
        assert matcher.is_isomorphic();perm=[matcher.mapping[i] for i in range(n)]
        assert sorted(perm)==list(range(n))
        assert all(g[a][b]['color']==target[perm[a]][perm[b]]['color'] for a,b in edges)
        permutations.append(perm)
    assert all(x==0 for row in r[:2] for x in row)
    return dict(denominator=den,inverse_quantization_scale=str(scale),
                inverse_numerators=[[[str(v) for v in row] for row in Q.tolist()] for Q in inverses],
                inverse_error_operator_norm_bounds=list(map(str,errors)),
                exact_approximate_Hessian=[[str(x) for x in row] for row in H.tolist()],
                Hessian_entry_error_bound=str(bound),positive_lower_matrix_leading_minors=list(map(str,minors)),
                sign_flip_vertex_permutations=permutations,sign_flips=[v.tolist() for v in flips],
                proof='For A=L+m²I>=m²I, ||A^-1-Q||2<=80 max|I-AQ|/m² exactly. Trace/submultiplicative bounds enclose each Hessian entry within epsilon. Htilde-5epsilon I has positive exact leading minors; the true Hessian is SPD. Inverse proposals may be approximate; all residuals, bounds and final positivity are rational integer checks. Exact colored-edge automorphisms force zero off-diagonal force, and equal-component displacements force zero diagonal traceless force.')


def exponential_ward():
    """Actual nonlinear weak-overlap anomaly and an explicit finite Green primitive."""
    import w33_pass11428_11432_native_alignment_measure_coarse_flags as Q
    import w33_pass11476_11480_cubic_geometry_correlated as M
    L=2;B=M.gauge_boundary(L);n=16;lap=B.T@B;mean=np.ones((n,n))/n
    # Integer spectrum0,4,8,12,16 gives an exact inverse on mean-zero vectors.
    inverse=sp.Matrix(np.rint(lap).astype(int).tolist())
    z=sp.symbols('lambda')
    assert sp.expand(inverse.charpoly(z).as_expr()-z*(z-4)**4*(z-8)**6*(z-12)**4*(z-16))==0
    projector=sp.ones(n)/n
    inv=(inverse+projector).inv()-projector
    assert inverse*inv==sp.eye(n)-projector
    K=-sp.Matrix(B.astype(int))*inv;Kf=np.array(K,float)
    src=read(SOURCES[3])['ward'];A=20*np.array(src['background']);J=np.array(src['anomaly_Jacobian'])
    rng=np.random.default_rng(11509);gauge=.001*rng.normal(size=n);dA=B@gauge
    def density(a):
        ans=np.zeros(n);gaps=[]
        for q,m in zip([1,-4,2,-3,6],[6,3,3,2,1]):
            P,_,_,gap,_,_=Q.weak_overlap(q,1.,a.reshape(n,4),[],L=L)
            ans+=m*q*np.trace(P.reshape(n,4,n,4),axis1=1,axis2=3).diagonal().real;gaps.append(gap)
        return ans,gaps
    a,gaps=density(A);ag,_=density(A+dA);flow=Kf@a
    # Free L2 antiperiodic Wilson gap is exactly1. Each directional
    # perturbation is block diagonal in (1 +/- gamma_mu)/2, hence bounded
    # by max_x |exp(i q A_mu)-1| <= |q| max_x |A_mu|.
    perturbation_bound=6*sum(sp.Rational(str(max(abs(A.reshape(n,4)[:,mu]))))+sp.Rational(1,10**14) for mu in range(4))
    assert perturbation_bound<1
    assert abs(sum(a))<1e-10 and np.linalg.norm(B.T@flow+a)<1e-10
    assert np.linalg.norm(ag-a)<1e-10 and np.linalg.norm(Kf@J@B)<1e-10
    P=np.eye(n)-lap/16;rho=.75;tail_rows=[]
    partial=np.zeros(64);term=a.copy()
    for tick in range(41):
        if tick in [0,1,2,4,8,16,32,40]:
            defect=np.linalg.norm(B.T@partial+a);bound=rho**tick*np.linalg.norm(a)
            assert defect<=bound+1e-11
            tail_rows.append(dict(terms=tick,Ward_defect=float(defect),proved_upper_bound=float(bound)))
        partial-=B@term/16;term=P@term
    torus_rates=[dict(L=ell,gap=4*np.sin(np.pi/ell)**2,rho=1-np.sin(np.pi/ell)**2/4) for ell in [2,3,4,8,16,32]]
    return dict(status='PASS',L=L,exact_Green_kernel=[[str(x) for x in row] for row in K.tolist()],
                background=A.tolist(),gauge_shift=dA.tolist(),anomaly_density=a.tolist(),primitive_flow=flow.tolist(),
                minimum_Wilson_gap=min(gaps),homotopy_gap_lower_bound=str(1-perturbation_bound),
                exact_Laplacian_spectrum={'0':1,'4':4,'8':6,'12':4,'16':1},
                anomaly_norm=float(np.linalg.norm(a)),Ward_residual=float(np.linalg.norm(B.T@flow+a)),
                gauge_invariance_residual=float(np.linalg.norm(ag-a)),prior_background_linearized_gauge_defect=float(np.linalg.norm(Kf@J@B)),
                local_series=tail_rows,box_gap_scaling=torus_rates,
                exact_primitive='k(A)=-B (B^T B)^+ a(A), sum a=0; equivalently -B/16 sum_{j>=0}(I-L/16)^j a(A).',
                exact_gauge_proof='Gauge transformation of Wilson sign projector is diagonal unitary conjugation; its onsite trace and a(A) are invariant. K is field-independent, so k(A) is gauge invariant. B^T K=-(I-mean).',
                locality_proof='Each j-th propagation term has graph radius j+1 relative to anomaly density, with Ward tail bounded by rho^N ||a||. Uniform exponential locality follows on a bounded-degree graph family only if its nonzero gap is uniformly positive and the anomaly functional is uniformly local.',
                scope='Nonlinear finite weak-field zero-index-sector primitive, with exact rational L2 kernel and verified gauge invariance. This is stronger than support-restricted differentiated flows, but NOT an infinite-volume local cohomological current: torus rho approaches1 as L grows, so this construction loses a uniform exponential bound. Flux sectors, integrability/curvature and complete chiral measure remain open. Luscher locality theorem and classical Green/Neumann identities are prior art.')


def syndrome_probabilities(T,words,signs):
    """Untwirled paired stabilizer expectation, not independent marginals."""
    q=np.empty((64,64));pair=4*words[:,None,:]+words[None,:,:]
    for a,b in itertools.product(range(64),repeat=2):
        weights=signs[:,None]*signs[None,:]
        for j in range(7):weights=weights*T[4*words[a,j]+words[b,j],pair[:,:,j]]
        q[a,b]=weights.sum()
    return q


def native_noisy_recovery():
    import w33_pass11471_11475_frames_currents_matching as Q
    import w33_pass11476_11480_cubic_geometry_correlated as M
    V,R=Q.recovery_rows();words,coeff=M.sparse_terms(V@V.T);assert len(words)==64
    signs=np.rint(64*coeff).astype(int);assert np.all(abs(signs)==1)
    # Tensor Pauli words follow reshape order: first wire is the highest bit.
    bits=np.array([[(j>>(6-k))&1 for k in range(7)] for j in range(128)])
    characters=np.empty((64,64),int)
    for j,word in enumerate(words):
        x=((word==1)|(word==2)).astype(int);z=((word==2)|(word==3)).astype(int)
        perm=np.arange(128)^sum(int(a)<<(6-i) for i,a in enumerate(x));phase=(1j)**sum(word==2)*(-1.)**(bits@z)
        characters[:,j]=np.rint(np.einsum('sai,sai,i->s',R[:,:,perm],R.conj(),phase).real/2).astype(int)
    assert np.array_equal(characters[0],signs) and np.array_equal(characters@characters.T,64*np.eye(64))
    row=next(z for z in read(SOURCES[2])['recovery']['rows'] if z['p']==.03 and z['relay']==0)
    T=np.array(row['physical_transfer']);q=syndrome_probabilities(T,words,signs)
    prior=characters@q@characters.T/4096
    assert min(prior.ravel())>-1e-10 and abs(prior.sum()-1)<1e-10
    prior=np.maximum(prior,0);prior/=prior.sum()
    qt=syndrome_probabilities(np.diag(np.diag(T)),words,signs);twirled=characters@qt@characters.T/4096
    bits6=np.array([[(j>>k)&1 for k in range(6)] for j in range(64)]);dist=np.sum(bits6[:,None]!=bits6[None,:],axis=2)
    rows=[]
    for error in [0.,.001,.01,.05]:
        confusion=error**dist*(1-error)**(6-dist)
        ga=np.argmax(confusion*prior.sum(axis=1)[None,:],axis=1)
        gb=np.argmax(confusion*prior.sum(axis=0)[None,:],axis=1)
        acceptance=np.zeros((64,64));joint=separate=0.
        for a,b in itertools.product(range(64),repeat=2):
            likelihood=confusion[a,:,None]*confusion[b,None,:];scores=prior*likelihood
            s,t=np.unravel_index(np.argmax(scores),scores.shape)
            acceptance[s,t]+=likelihood[s,t];joint+=scores[s,t];separate+=scores[ga[a],gb[b]]
        assert joint>=separate-1e-12 and np.max(acceptance)<=1+1e-12
        assert -1e-12<=joint<=1+1e-12 and -1e-12<=separate<=1+1e-12
        joint=float(np.clip(joint,0,1));separate=float(np.clip(separate,0,1))
        rows.append(dict(bit_error=error,joint_correct_syndrome_probability=float(joint),separate_correct_syndrome_probability=float(separate),acceptance=acceptance.tolist()))
    leak=read(SOURCES[3])['recovery']['native_leakage'];U=dec(leak['cell_unitary']);K=dec(leak['logical_subchannel'])
    effects=np.linalg.eigvalsh(np.kron(K.conj().T@K,K.conj().T@K));lmin=1-max(effects);lmax=1-min(effects)
    assert lmin>=-1e-10 and lmax>1e-4 and np.linalg.norm(U.conj().T@U-np.eye(12))<1e-10
    v=.001
    for row in rows:
        correct=row['joint_correct_syndrome_probability'];wrong=1-correct
        row.update(first_wrong_syndrome_false_accept=wrong*v,
                   correct_branch_weight=correct*(1-v),
                   native_cell_leakage_false_accept_bounds=[correct*(1-v)*v*max(0,lmin),correct*(1-v)*v*lmax])
    return dict(status='PASS',physical_noise=.03,relay=0,physical_transfer=T.tolist(),stabilizer_words=words.tolist(),
                syndrome_characters=characters.tolist(),syndrome_prior=prior.tolist(),rows=rows,
                untwirled_vs_twirled_prior_difference=float(np.linalg.norm(prior-twirled)),
                correlated_vs_marginal_prior_difference=float(np.linalg.norm(prior-np.outer(prior.sum(axis=1),prior.sum(axis=0)))),
                cell_unitary=enc(U),logical_cell_subchannel=enc(K),two_cell_leakage_probability_bounds=[float(lmin),float(lmax)],verification_flip=v,
                named_instrument='Encode V tensor V; apply seven copies of native two-qubit CPTP relay channel; measure actual stabilizer projectors R_s^dagger R_s tensor R_t^dagger R_t; draw12 noisy bits; joint MAP correction; projectively verify code membership, then flip the classical result with probability v. Keep false accepted wrong-syndrome density in its physical16380-dimensional orthogonal subspace. Correct syndrome branches decode using R_s tensor R_t, embed in actual two12-cells, apply U tensor U, verify actual logical-cell projector with another classical flip v. Retain its140-dimensional accepted leakage subspace and all rejected branches as a trace-preserving erasure flag.',
                exact_CP_TP_proof='Encoding is isometric, physical channel CPTP, stabilizer projectors resolve identity, corrections unitary, readout/flip matrices stochastic, native U unitary, final projective sectors orthogonal and complete. Their composition and flagged direct sum is CPTP without twirling or postselection.',
                scope='Native untwirled correlated syndrome prior and noisy joint inference computed, with a named complete verification/leakage instrument and state-independent actual-cell leakage bounds. MAP optimizes correct-syndrome probability for maximally mixed logical input, not quantum fidelity. Correct syndrome can retain uncorrectable logical errors. Syndrome extraction/correction/verification hardware faults beyond supplied readout flips remain unmodeled; no fault-tolerance threshold or physical gate implementation is claimed.')


def run():
    result={}
    for name,fn in [('vacuum',balanced_vacuum),('SM',sm_blocks),('metric',metric_stabilization),('ward',exponential_ward),('recovery',native_noisy_recovery)]:
        result[name]=fn();print(name,'PASS',flush=True)
        Path('/tmp/w33_11506_partial.json').write_text(json.dumps(result))
    result.update(status='PASS',passes=list(range(11506,11511)),reservation='09cc50896',
        source_sha256={p:hashlib.sha256(json.dumps(read(p),sort_keys=True,separators=(',',':')).encode()).hexdigest() for p in SOURCES})
    OUT.write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':run()
