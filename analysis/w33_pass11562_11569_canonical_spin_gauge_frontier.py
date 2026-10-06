#!/usr/bin/env python3
"""Passes 11562-11569: canonical spin, connection, cohomology, exceptional stabilizers,
axis/Wilson Dirac, heat trace, and winding.

All claims are finite/exact or controlled lattice-continuum calculations.  The script
writes one frozen JSON certificate per pass.
"""
from __future__ import annotations
import itertools, json, math, sys
from collections import Counter
from fractions import Fraction as Fr
from math import lcm
from pathlib import Path
import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_pass10950_clock_albert_lorentz_spinor as P50
import w33_pass11557_pin_equivariant_event_dirac as P57

def dump(passno,name,obj):
    obj=dict(obj)
    obj["pass"]=passno
    obj["schema"]=f"w33.pass{passno}.{name}.v1"
    (ROOT/"data"/f"PART_W33_PASS{passno}_{name.upper()}.json").write_text(
        json.dumps(obj,indent=2,sort_keys=True)+"\n",encoding="utf-8")

def parse_matrix(rows):
    return sp.Matrix([[sp.sympify(x) for x in row] for row in rows])

def rank_mod_p(A,p=3):
    A=np.array(A,dtype=np.int64)%p
    r=0
    for c in range(A.shape[1]):
        q=next((i for i in range(r,A.shape[0]) if A[i,c]%p),None)
        if q is None: continue
        if q!=r: A[[r,q]]=A[[q,r]]
        inv=pow(int(A[r,c]),-1,p)
        A[r]=(A[r]*inv)%p
        nz=np.flatnonzero(A[:,c])
        for i in nz:
            i=int(i)
            if i!=r:
                A[i]=(A[i]-A[i,c]*A[r])%p
        r+=1
        if r==A.shape[0]: break
    return r

def perm_matrix(p):
    n=len(p); M=np.zeros((n,n),dtype=np.int64)
    for i,j in enumerate(p): M[j,i]=1
    return M

def compose(p,q):
    return tuple(p[q[i]] for i in range(len(p)))

def signed_WD3():
    out=[]
    for pm in itertools.permutations(range(3)):
        for sg in itertools.product((1,2),repeat=3):
            prod=1
            for s in sg: prod*=1 if s==1 else -1
            if prod!=1: continue
            M=np.zeros((3,3),dtype=np.int64)
            for i,j in enumerate(pm): M[i,j]=sg[i]
            out.append(M%3)
    uniq={tuple(M.ravel()):M for M in out}
    return list(uniq.values())

def build_cartan_matrix(E):
    pairs=[(0,1),(0,2),(1,2)]
    unk=[(i,a,b) for i in range(3) for a,b in pairs]
    rows=[]
    for a in range(3):
        for i,j in pairs:
            row=[]
            for ii,aa,bb in unk:
                v=0
                if ii==i:
                    if a==aa:v+=E[bb,j]
                    if a==bb:v-=E[aa,j]
                if ii==j:
                    if a==aa:v-=E[bb,i]
                    if a==bb:v+=E[aa,i]
                row.append(v)
            rows.append(row)
    return sp.Matrix(rows),unk

def axis_operators(N=3):
    V=list(itertools.product(range(N),repeat=3)); idx={v:i for i,v in enumerate(V)}
    T=[]
    for a in range(3):
        M=np.zeros((N**3,N**3),dtype=np.int64)
        for i,x in enumerate(V):
            y=list(x); y[a]=(y[a]+1)%N
            M[i,idx[tuple(y)]]=1
        T.append(M)
    Ti=[np.linalg.matrix_power(M,N-1).astype(np.int64) for M in T]
    delta=[T[i]-Ti[i] for i in range(3)]
    return V,T,Ti,delta

def arrow_matrix_N3():
    V=list(itertools.product(range(3),repeat=3)); idx={v:i for i,v in enumerate(V)}
    D=np.zeros((27,27),dtype=np.int64)
    for i,x in enumerate(V):
        for j,y in enumerate(V):
            if i==j:continue
            d=tuple((y[k]-x[k])%3 for k in range(3))
            if all(d):
                eps=[1 if z==1 else -1 for z in d]
                D[i,j]=eps[0]*eps[1]*eps[2]
    return D

def heat_rows():
    rows=[]
    for t in [0.5,1.0]:
        target=4*(sum(math.exp(-t*p*p) for p in range(-60,61)))**3
        for N in [9,27,81]:
            a=2*math.pi/N
            k=2*math.pi*np.arange(N)/N
            s=np.sin(k); c=np.cos(k)
            q2=(s[:,None,None]**2+s[None,:,None]**2+s[None,None,:]**2)/(a*a)
            W=2*((1-c)[:,None,None]+(1-c)[None,:,None]+(1-c)[None,None,:])
            reg=q2+(W/a)**2
            hu=4*np.exp(-t*q2).sum()
            hr=4*np.exp(-t*reg).sum()
            rows.append(dict(t=t,N=N,target=target,
                             unregulated_ratio=float(hu/target),
                             wilson_ratio=float(hr/target)))
    return rows

def winding_numeric(m,r=1.0,N=41):
    k=(np.arange(N)+0.5)*2*np.pi/N
    x,y,z=np.meshgrid(k,k,k,indexing="ij")
    sx,sy,sz=np.sin(x),np.sin(y),np.sin(z)
    cx,cy,cz=np.cos(x),np.cos(y),np.cos(z)
    M=m+2*r*((1-cx)+(1-cy)+(1-cz))
    h=np.stack([sx,sy,sz,M],axis=-1)
    hx=np.stack([cx,np.zeros_like(x),np.zeros_like(x),2*r*sx],axis=-1)
    hy=np.stack([np.zeros_like(x),cy,np.zeros_like(x),2*r*sy],axis=-1)
    hz=np.stack([np.zeros_like(x),np.zeros_like(x),cz,2*r*sz],axis=-1)
    A=np.stack([h,hx,hy,hz],axis=-1)
    det=np.linalg.det(A)
    norm2=sx*sx+sy*sy+sz*sz+M*M
    dens=det/(norm2*norm2)
    return float(dens.mean()*(2*np.pi)**3/(2*np.pi**2))

def main():
    # ------------------------------------------------------------------ 11562
    gcert=json.loads((ROOT/"data/w33_pass10961_albert_clifford9_gammas.json").read_text())
    gamma=[parse_matrix(x) for x in gcert["gamma9"]]
    I16=sp.eye(16)
    for i in range(9):
        for j in range(9):
            assert gamma[i]*gamma[j]+gamma[j]*gamma[i]==(2 if i==j else 0)*I16
    sel=(0,1,2)
    products=[]
    for mask in range(8):
        M=I16
        for i in range(3):
            if (mask>>i)&1:M=M*gamma[i]
        products.append(M)
    algdim=sp.Matrix.hstack(*[M.reshape(256,1) for M in products]).rank()
    assert algdim==8
    biv=[]; pairs=[]
    for i,j in itertools.combinations(range(9),2):
        pairs.append((i,j));biv.append(gamma[i]*gamma[j])
    assert sp.Matrix.hstack(*[B.reshape(256,1) for B in biv]).rank()==36
    spin3=[gamma[0]*gamma[1],gamma[0]*gamma[2],gamma[1]*gamma[2]]
    blocks=[]
    for S in spin3:
        blocks.append(sp.Matrix.hstack(*[(B*S-S*B).reshape(256,1) for B in biv]))
    eq=sp.Matrix.vstack(*blocks)
    centdim=len(eq.nullspace())
    assert centdim==15
    plane_stab_dim=sum((i in sel)==(j in sel) for i,j in pairs)
    assert plane_stab_dim==18
    p50=json.loads((ROOT/"data/w33_pass10950_clock_albert_lorentz_spinor.json").read_text())
    p56=json.loads((ROOT/"data/w33_pass10956_albert_spin8_halfspin_normalizer.json").read_text())
    assert p50["lorentz"]["Der_c_dimension"]==36
    assert p56["spin8_stabilizer"]["dimension"]==28
    out62=dict(
      status="PASS_CL3_EMBEDDING_CLASSIFIED_BUT_NOT_CANONIC",
      clifford=dict(selected_gamma_indices=list(sel),generated_Cl3_algebra_dimension=algdim,
                    spin9_bivector_dimension=36),
      orbit=dict(coordinate_three_planes=math.comb(9,3),
                 spin9_plane_stabilizer_dimension=plane_stab_dim,
                 spin9_plane_orbit_dimension=36-plane_stab_dim,
                 centralizer_of_event_spin3_dimension=centdim,
                 centralizer_identification="so(6) ~= su(4)"),
      no_go="The committed Cl(9)/Spin(9) data alone do not select a unique event three-plane: Spin(9) acts transitively on oriented orthonormal 3-planes, with stabilizer Spin(3)xSpin(6).",
      frame_note="Pass10956 supplies a Spin(8) selector after a spatial axis is chosen, but that still leaves continuous plane choices; no existing Hesse/clock certificate is used here to assert a unique Cl3.",
      theorem="The Pass11558 choice of three Albert gamma directions is a valid Cl3 embedding, but it is not canonical under the parent Spin(9) symmetry. Its plane stabilizer has dimension 18 and its internal Spin(3) centralizer has dimension 15.",
      boundary="Selection theorem/no-go only; no fermion-family or spacetime-axis interpretation is derived.")
    dump(11562,"canonical_cl3_inside_albert_cl9",out62)

    # ------------------------------------------------------------------ 11563
    e=sp.symbols("e0:9")
    E=sp.Matrix(3,3,e)
    A,unk=build_cartan_matrix(E)
    detA=A.det(method="domain-ge")
    target=-2*E.det()**3
    assert sp.Poly(sp.expand(detA-target),*e).is_zero
    E0=sp.Matrix([[1,1,0],[0,1,1],[1,0,1]])
    A0,_=build_cartan_matrix(E0)
    cvec=sp.Matrix([2,-1,3,4,1,-2,-3,5,2])
    omega=A0.LUsolve(-cvec)
    assert A0*omega+cvec==sp.zeros(9,1)
    out63=dict(
      status="PASS_UNIQUE_DISCRETE_TORSION_FREE_METRIC_CONNECTION",
      local_system=dict(unknowns=9,equations=9,
                        unknown_layout="omega_i^{ab}, i=1..3, ab in {12,13,23}",
                        equation="C^a_ij + omega_i^a_b E^b_j - omega_j^a_b E^b_i = 0",
                        metric_compatibility="omega_i^{ab}=-omega_i^{ba}"),
      exact_determinant="det A(E) = -2 (det E)^3",
      uniqueness="For every nondegenerate 3x3 frame E, the local torsion-free metric-compatible spin connection is unique.",
      exact_sample=dict(frame=[[1,1,0],[0,1,1],[1,0,1]],det_frame=int(E0.det()),
                        solution=[str(sp.simplify(x)) for x in omega]),
      curvature_next="Once omega links are fixed kinematically, plaquette holonomy supplies a discrete curvature; no curvature action is selected by this pass.",
      theorem="The frame ambiguity left by Pass11559 is kinematically removable: the discrete Cartan equations have determinant -2(det E)^3, hence a unique local metric-compatible torsion-free connection for every invertible frame.",
      boundary="Kinematic Levi-Civita analogue only. No Einstein-Hilbert dynamics, Lorentzian signature, Newton constant, or frame equation of motion is derived.")
    dump(11563,"discrete_spin_connection",out63)

    # ------------------------------------------------------------------ 11564
    W=signed_WD3(); assert len(W)==24
    inv_eq=np.vstack([(M.T-np.eye(3,dtype=np.int64))%3 for M in W])
    h1T_fixed=3-rank_mod_p(inv_eq,3)
    assert h1T_fixed==0
    S4=list(itertools.permutations(range(4)))
    pidx={p:i for i,p in enumerate(S4)}
    # permutation-module 1-cocycles
    eqs=[]
    for g in S4:
        Pg=perm_matrix(g)%3
        for h in S4:
            gh=compose(g,h)
            for a in range(4):
                row=np.zeros(24*4,dtype=np.int64)
                row[pidx[gh]*4+a]+=1
                row[pidx[g]*4+a]-=1
                for b in range(4):
                    row[pidx[h]*4+b]-=Pg[a,b]
                eqs.append(row%3)
    zdim=24*4-rank_mod_p(np.array(eqs),3)
    cob=[]
    for g in S4:
        cob.append((perm_matrix(g)-np.eye(4,dtype=np.int64))%3)
    Bmap=np.vstack(cob)
    bdim=rank_mod_p(Bmap,3)
    fixedM=4-rank_mod_p(Bmap,3)
    assert (zdim,bdim,fixedM)==(3,3,1)
    # trivial-module H1(S4,F3)
    eqt=[]
    for g in S4:
        for h in S4:
            row=np.zeros(24,dtype=np.int64)
            row[pidx[compose(g,h)]]+=1;row[pidx[g]]-=1;row[pidx[h]]-=1
            eqt.append(row%3)
    h1W=24-rank_mod_p(np.array(eqt),3)
    assert h1W==0
    out64=dict(
      status="PASS_LHS_NO_GROUP_H1_ARROW_AND_RESIDUAL_FIXED_LINE",
      five_term=dict(H1_W_trivial_F3_dimension=h1W,
                     H1_T_F3_W_invariant_dimension=h1T_fixed,
                     conclusion_H1_affine_PSp_trivial_coefficients_dimension=0),
      direction_module=dict(module="M=F3^4 on the four positive temporal directions",
                            fixed_dimension=fixedM,
                            cocycle_dimension=zdim,coboundary_dimension=bdim,
                            H1_W_M_dimension=zdim-bdim,
                            fixed_vector=[1,1,1,1]),
      interpretation="The arrow is not an ordinary F3-valued group 1-cocycle/character of 3^3:W(D3). It lives as the unique H^0(W(D3),M) line of the translation-coinvariant direction module from Pass11560.",
      theorem="Low-degree LHS/inflation-restriction data force H^1(3^3:W(D3),F3)=0, while the four-direction coefficient module has one invariant line and no nontrivial H^1. The mod-3 arrow is therefore a residual coefficient-module invariant, not a hidden affine-group character.",
      boundary="This is low-degree group/coefficient-module cohomology; it is not a thermodynamic arrow, anomaly, or continuum Chern class.")
    dump(11564,"lhs_arrow_origin",out64)

    # ------------------------------------------------------------------ shared Albert trace computation for 11565
    J=P50.build_clock_albert()
    n,prodF=J["n"],J["prod"]
    den=lcm(*[x.denominator for x in prodF.flat])
    Pi=np.array([[[int(x*den) for x in prodF[i,j]] for j in range(n)] for i in range(n)],dtype=np.int64)
    Ls=np.stack([Pi[i].T for i in range(n)])
    tr=np.array([np.trace(Ls[i]) for i in range(n)],dtype=np.int64)
    comm_trace_max=0
    for i,j in itertools.combinations(range(n),2):
        C=Ls[i]@Ls[j]-Ls[j]@Ls[i]
        comm_trace_max=max(comm_trace_max,int(np.max(np.abs(tr@C))))
    assert comm_trace_max==0
    tl=P50.nullspace([tr.tolist()]); assert len(tl)==26
    tv=[]
    for v in tl:
        q=lcm(*[x.denominator for x in v])
        tv.append(np.array([int(x*q) for x in v],dtype=np.int64))
    def Lv(v): return np.tensordot(np.asarray(v,dtype=np.int64),Ls,axes=1)
    image_rows=[tr]+[tr@Lv(v) for v in tv]
    rank_trace_action=P50.rank([list(map(int,r)) for r in image_rows])
    assert rank_trace_action==27
    parent=json.loads((ROOT/"data/w33_pass10950_clock_albert_lorentz_spinor.json").read_text())
    assert parent["structure"]["Der_dimension"]==52 and parent["structure"]["str0_dimension"]==78
    out65=dict(
      status="PASS_HESSE_TRACE_CUBE_STABILIZER_IS_F4",
      executable_albert=dict(derivation_dimension=52,reduced_structure_dimension=78,
                             trace_covector_nonzero_entries=int(np.count_nonzero(tr)),
                             all_derivation_commutators_annihilate_trace=True,
                             traceless_multiplication_dimension=26,
                             rank_of_trace_plus_traceless_trace_images=rank_trace_action),
      stabilizer=dict(
        trace_line_in_E6_dimension=52,
        identification="Der(J3(O)) = compact F4 inside the executable E6(-26)",
        hesse_two_plane_stabilizer_dimension=52,
        reason="E6 preserves the Albert norm N. Infinitesimal preservation of span{N,T^3} forces preservation of the trace line T because N is not divisible by T^2; the executable trace-line stabilizer contains no nonzero traceless multiplication."),
      witness_trace_zero_nonzero_norm="diag(1,1,-2) has T=0 and N=-2, excluding an N component in 3 T^2 deltaT = aN+bT^3",
      theorem="The full continuous Lie-algebra stabilizer inside the committed E6(-26) of the Hesse polynomial-law two-plane <N,T^3> is exactly the 52-dimensional Albert derivation algebra F4.",
      boundary="This identifies the exceptional stabilizer; it does not derive a vacuum that selects T^3 or the Standard Model subgroup of F4.")
    dump(11565,"e6_hesse_stabilizer_f4",out65)

    # ------------------------------------------------------------------ 11566
    V,T,Ti,delta=axis_operators(3)
    Daxis=delta[0]@delta[1]@delta[2]
    Dold=arrow_matrix_N3()
    assert np.array_equal(Daxis,Dold)
    A1=sum(T[i]+Ti[i] for i in range(3))
    W1=6*np.eye(27,dtype=np.int64)-A1
    assert np.array_equal(-sum(d@d for d in delta),W1)
    gam,grade=P57.small_clifford()
    gn=[np.array(g.evalf(40),dtype=complex) for g in gam]
    g0=np.array(grade.evalf(40),dtype=complex)
    Qaxis=1j*sum(np.kron(gn[i],delta[i]) for i in range(3))
    assert np.linalg.norm(Qaxis-Qaxis.conj().T)<1e-10
    assert np.linalg.norm(Qaxis@Qaxis-np.kron(np.eye(4),W1))<1e-10
    corners=[]
    for bits in itertools.product((0,1),repeat=3):
        npi=sum(bits); corners.append(dict(bits=list(bits),chirality=(-1)**npi,W1=4*npi))
    assert Counter(c["chirality"] for c in corners)==Counter({1:4,-1:4})
    out66=dict(
      status="PASS_AXIS_DIRAC_AND_HAMMING_WILSON_REGULATOR",
      finite_N3=dict(axis_differences="delta_i=T_i-T_i^{-1} on the Hamming weight-1 orbital",
                     exact_arrow_factorization="D = delta_1 delta_2 delta_3",
                     axis_dirac="Q_axis = i sum Gamma_i tensor delta_i",
                     exact_square="Q_axis^2 = I4 tensor (6I-A1)",
                     Wilson_operator="W1=6I-A1 = -sum delta_i^2"),
      continuum=dict(
        axis_symbol="Q_axis/(2a): -sum_i Gamma_i sin(k_i)/a",
        unregulated_zero_set="8 Brillouin-zone corners {0,pi}^3",
        corner_data=corners,
        wilson_dirac="H_W = Q_axis/(2a) + Gamma_0 r W1/a",
        exact_square="H_W^2 = Q_axis^2/(4a^2) + r^2 W1^2/a^2",
        wilson_zero_set="origin only for r != 0 and zero bare mass"),
      correction_to_null_edge_refinement=dict(
        pass11559_local_principal_symbol_remains_valid=True,
        naive_body_diagonal_symbol="B_i ~ 8 i sin(k_i) prod_{j!=i} cos(k_j)",
        global_pathology="Besides the 8 corners, the naive body-diagonal lift has 12 nodal lines where at least two cosines vanish.",
        resolution="Refine the exact cubic arrow through its weight-1 factorization delta_1 delta_2 delta_3; use the weight-1 Hamming orbital for the first-order Dirac operator."),
      albert_weld="Replace the 4x4 gamma carrier by the Pass11558 Albert Peirce-16 Cl4 submodule; the finite operator is then four copies of the same axis Dirac module before further couplings.",
      theorem="The rank-5 Hamming fission contains its own canonical first-order refinement and Wilson regulator: D factors through the three weight-1 axis differences, Q_axis^2 is exactly the weight-1 Hamming Laplacian, and adding that same Laplacian as a graded Wilson term removes every Brillouin-zone zero except the origin.",
      boundary="A regulator/continuum operator theorem. The Wilson coefficient, bare mass/Hesse coupling, physical lattice spacing, and Lorentzian continuation remain inputs.")
    dump(11566,"axis_dirac_hamming_wilson",out66)

    # ------------------------------------------------------------------ 11567
    hrows=heat_rows()
    row81_t1=next(r for r in hrows if r["N"]==81 and r["t"]==1.0)
    assert abs(row81_t1["unregulated_ratio"]-8)<0.05
    assert abs(row81_t1["wilson_ratio"]-1)<0.03
    out67=dict(
      status="PASS_WILSON_SPECTRAL_ACTION_SINGLE_SPECIES_LIMIT",
      normalization="a=2pi/N on a periodic 2pi three-torus; Q_axis/(2a) has continuum symbol gamma.p",
      continuum_target="4 [sum_{n in Z} exp(-t n^2)]^3",
      rows=hrows,
      asymptotic_reading="The unregulated axis Dirac tends to eight continuum species; the Hamming-Wilson operator tends to one.",
      N81_t1=dict(unregulated_ratio=row81_t1["unregulated_ratio"],wilson_ratio=row81_t1["wilson_ratio"]),
      theorem="The heat trace distinguishes the lattice species count: after the canonical Hamming-Wilson term, the fixed-dimensional 3^r tower converges numerically to the one-Dirac-spinor heat trace rather than the eight-species unregulated limit.",
      boundary="Free flat-torus spectral-action test only. Gauge curvature, dynamical frame integration, renormalized couplings, and gravitational heat-kernel coefficients are not derived.")
    dump(11567,"wilson_spectral_action",out67)

    # ------------------------------------------------------------------ 11568
    samples=[]
    for m in [2,-2,-6,-10,-14]:
        exact=0
        for npi in range(4):
            w=math.comb(3,npi)*((-1)**npi)
            mass=m+4*npi
            exact += w*(1 if mass>0 else -1)
        exact=-exact/2
        num=winding_numeric(m,N=41)
        assert abs(num-exact)<2e-6
        samples.append(dict(m=m,corner_formula=int(exact),numerical_winding=num))
    phases=[
      dict(interval="m>0",winding=0),
      dict(interval="-4<m<0",winding=1),
      dict(interval="-8<m<-4",winding=-2),
      dict(interval="-12<m<-8",winding=1),
      dict(interval="m<-12",winding=0),
    ]
    out68=dict(
      status="PASS_INTEGER_WILSON_DIRAC_WINDING_PHASES",
      map="h(k)=(sin k1,sin k2,sin k3,m+2 sum_i(1-cos ki)); normalized h defines T3 -> S3 away from m=0,-4,-8,-12",
      corner_formula="nu = -1/2 sum_{n=0}^3 C(3,n)(-1)^n sign(m+4n)",
      phases=phases,
      direct_integral_samples=samples,
      relation_to_old_4plus4="The eight unregulated corners carry Jacobian signs (-1)^n, four positive and four negative; the Wilson mass converts that cancellation data into quantized winding phases.",
      theorem="The Hamming-Wilson Dirac family has a genuine integer degree. For unit Wilson coefficient its phases are 0,+1,-2,+1,0 across the five mass intervals.",
      boundary="Lattice topological invariant of the regulated free symbol. It is not yet an observed chiral index, anomaly coefficient, or Standard Model family count.")
    dump(11568,"wilson_dirac_winding",out68)

    # ------------------------------------------------------------------ 11569
    sm_dim=8+3+1
    assert sm_dim==12
    out69=dict(
      status="PASS_CURRENT_SELECTORS_STOP_AT_SPIN6_NOT_SM",
      exact_chain=[
        dict(selector="Hesse trace-cube line inside E6",algebra="f4",dimension=52,source="Pass11565"),
        dict(selector="primitive Albert idempotent / Peirce 16",algebra="spin9",dimension=36,source="Pass10950"),
        dict(selector="event Cl3 three-plane preserved setwise",algebra="spin3 + spin6",dimension=18,source="Pass11562"),
        dict(selector="commute with event spin3",algebra="spin6 ~= su4",dimension=15,source="Pass11562"),
      ],
      standard_model_target=dict(algebra="su3 + su2 + u1",dimension=sm_dim,
        repo_prior_art="Exploratory/documentary files quote the Todorov-Dubois-Violette intersection S(U2 x U3)=Spin9 intersect (SU3 x SU3)/Z3 inside F4."),
      no_go="The new Hesse+Peirce+event-Cl3 selectors by themselves leave a 15-dimensional Spin(6) centralizer, not the 12-dimensional Standard Model gauge algebra. Therefore the SM intersection cannot be claimed from the present selector set.",
      missing_selector="An independently constructed SU(3)xSU(3) / lepton-quark splitting (or equivalent exact Albert datum) must be welded to the executable F4 matrices before the Todorov intersection becomes a repository-derived gauge theorem.",
      outside_box_candidate="Spin6 ~= SU4 is precisely the Pati-Salam color-plus-lepton algebra; an additional exact complex-line/lepton-quark selector would naturally test SU4 -> U3, but weak SU2 still requires an independent internal selector rather than reusing spacetime Spin3.",
      theorem="The exceptional selector chain is now explicit and falsifiable: E6 -> F4 -> Spin9 -> Spin3 x Spin6, with internal centralizer Spin6. It does not yet equal the Standard Model algebra.",
      boundary="Gauge-subgroup selection result/no-go. No identification of event Spin3 with electroweak SU2 is made.")
    dump(11569,"exceptional_gauge_subgroup_search",out69)

    print(json.dumps({
      "status":"PASS_11562_11569_FRONTIER",
      "11562_centralizer":centdim,
      "11563_det":"-2(detE)^3",
      "11564_affine_H1":0,
      "11565_stabilizer":"F4 dim52",
      "11566_arrow_axis_factor":True,
      "11567_N81_t1":row81_t1,
      "11568_phases":[p["winding"] for p in phases],
      "11569_residual":"Spin6 dim15"
    },indent=2))

if __name__=="__main__":
    main()
