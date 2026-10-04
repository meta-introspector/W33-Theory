#!/usr/bin/env python3
"""Five follow-ups to11403-11407. Conditional actions, not observed physics.
Prior: Pass236/11326/11347 (CP), BT1041/1042,11271,11301 (derivations),
11389 (FCC),11403-11407 (actual frames and star primitives).
"""
from pathlib import Path
import hashlib,json
from itertools import combinations,product
import numpy as np
import sympy as sp
from scipy.linalg import expm,schur,block_diag
from sympy.matrices.normalforms import hermite_normal_form
from w33_pass11390_11394_context_matter_clock import load
from w33_pass11403_11407_native_bimodule_collective_gates import cm,gellmann
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/w33_pass11408_11412_clock_spurion_seam_ticks.json'
def read(d):return np.array(d['real'])+1j*np.array(d['imag'])
def prior():return json.loads((ROOT/'data/w33_pass11403_11407_native_bimodule_collective_gates.json').read_text())
def differences(w):return np.prod([w[j]-w[i]for i,j in combinations(range(3),2)])
def mixing(U):
    s13=abs(U[0,2])**2;s12=abs(U[0,1])**2/(1-s13);s23=abs(U[1,2])**2/(1-s13)
    J=np.imag(U[0,0]*U[1,1]*U[0,1].conjugate()*U[1,0].conjugate())
    den=np.sqrt(s12*(1-s12)*s23*(1-s23)*s13)*(1-s13)
    return dict(angles_degrees=(np.arcsin(np.sqrt([s12,s23,s13]))*180/np.pi).tolist(),J=float(J),sin_delta=float(J/den)if den>1e-16 else None)

def competing_clocks(p):
    F=read(p['colour_safe_bridge']['other_native_clock_mixing'])
    # Fixed distinct orbit radii: an orientation mechanism, not a mass fit.
    d=np.array([.01,.1,1.]);D=np.diag(d);other=F@D@F.conj().T;rows=[]
    JF=float(np.imag(F[0,0]*F[1,1]*F[0,1].conjugate()*F[1,0].conjugate()))
    r0,z=sp.symbols('r z');exactd=[sp.Rational(1,100),sp.Rational(1,10),sp.Integer(1)]
    a=sum(exactd);b0=sum(exactd[i]*exactd[j]for i,j in combinations(range(3),2));c0=sp.prod(exactd)
    characteristic=z**3-a*(1+r0)*z**2+(b0*(1+r0**2)+2*a*a*r0/3)*z-c0*(1+r0**3)-a*b0*(r0+r0**2)/3
    discriminant=sp.factor(sp.discriminant(characteristic,z))
    positive=sp.Rational(531441,250000000000)*(3025*(r0**6+1)+151959*r0**2*(r0-1)**2+208662*r0**3)
    assert sp.expand(discriminant-positive)==0
    Fc=[read(x)for x in p['fermion_bimodule']['native_slot_frames']['colour_qutrit_seeds']]
    C=np.hstack(Fc);Gs=[sum(A@T@A.conj().T for A in Fc)for T in gellmann()]
    for r in (.0001,.001,.01,.1,1.,10.):
        B=D+r*other;b,U=np.linalg.eigh(B);Hd=U@D@U.conj().T
        # von Neumann trace inequality certifies this global orbit minimum.
        curv=[2*(d[j]-d[i])*(b[j]-b[i])for i,j in combinations(range(3),2)]
        assert min(curv)>0 and np.linalg.norm(B@Hd-Hd@B)<1e-12
        m=mixing(U);triangle=np.imag(B[0,1]*B[1,2]*B[2,0])
        assert abs(abs(m['J'])-abs(triangle/differences(b)))<1e-12
        assert abs(abs(m['J'])-abs(r**3*JF*differences(d)/differences(b)))<1e-12
        Hu=D@D;Hds=U@(D@D)@U.conj().T;comm=Hu@Hds-Hds@Hu
        cp=np.imag(np.trace(comm@comm@comm));expected=6*m['J']*differences(d*d)**2
        assert abs(abs(cp)-abs(expected))<1e-12
        native=C@np.kron(Hd,np.eye(3))@C.conj().T
        colour_res=max(np.linalg.norm(native@G-G@native)for G in Gs);assert colour_res<1e-10
        rows.append(dict(r=r,source=cm(B),orientation=cm(U),selected_family_operator=cm(Hd),source_eigenvalues=b.tolist(),orbit_Hessian_eigenvalues=np.repeat(curv,2).tolist(),energy=float(-np.dot(d,b)),triangle_flux=float(triangle),cubic_commutator=float(cp),colour_commutator_residual=float(colour_res),**m))
    # A one-source Hermitian perturbation would SPOIL small masses if used as
    # the Yukawa itself. Restricting to an orbit keeps them fixed by assumption.
    h=p['colour_safe_bridge']['scans'][-1]['family_singular_values'];Bh=np.diag(h)+.001*F@np.diag(h)@F.conj().T
    hd=np.sort(h);near_identity_ratio=(hd[2]-hd[1])/(hd[2]-hd[0])
    return dict(status='PASS',action='V(U)=-Tr[U D U^dag (D+r F D F^dag)], U on U(3)/U(1)^3',supplied_orbit_radii=d.tolist(),scans=rows,
        native_hierarchy_additive_control=dict(original=h,additive_source_eigenvalues=np.linalg.eigvalsh(Bh).tolist()),
        exact_J_magnitude='|J(r)|=r^3 |J(F)| product_(i<j)(dj-di)/product_(i<j)(bj-bi)',
        near_identity_J_over_r3=JF,
        small_mixing_ratio_obstruction=dict(theta13_over_theta23_limit=float(near_identity_ratio),reason='A Hermitian function of the conjugate3-clock is circulant: all three offdiagonal magnitudes are equal. First-order angles are r|source_ij|/(dj-di); highly hierarchical spectra therefore force theta13 approximately theta23.',scope='small-r nondegenerate perturbation theory; not a no-go for additional clocks, spurions or strong mixing'),
        global_minimum='sorted eigenframes of B=D+r F D F^dag pair with sorted D; normal curvatures2(dj-di)(bj-bi), twice each',
        exact_characteristic_polynomial=str(characteristic),exact_discriminant=str(discriminant),positive_discriminant_decomposition=str(positive),
        all_positive_r_orientation_theorem='For supplied D=diag(1/100,1/10,1), the exact source discriminant is strictly positive for every r>=0; the orbit minimum is unique modulo diagonal phases and has six positive normal curvatures.',
        CP_scope='The chosen oriented complex clock sources carry explicit CP breaking; conjugating all sources reverses J. No spontaneous or observed CP prediction.',
        scope='Native competing clock sources choose a stable intermediate family orientation. r, fixed orbit spectra and pin/orientation are supplied; no measured CKM fit or dynamical generation of mass radii. Additive use of the source lifts light eigenvalues, so it is not interchangeable with an orbit alignment action.')

def yukawa_arrows(p):
    slots=p['fermion_bimodule']['native_slot_frames'];Fc=[read(x)for x in slots['colour_qutrit_seeds']];C=np.hstack(Fc)
    W,L,Ru,Rd=[read(slots[n])for n in ['weak','lepton','right_up','right_down']]
    col=[sum(A@T@A.conj().T for A in Fc)for T in gellmann()]
    weak=[W@T@W.conj().T/2 for T in [np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1])]]
    q,h=sp.symbols('q h');charge={k:sp.sympify(v)for k,v in p['fermion_bimodule']['charge_family'].items()}
    # Nonzero bare singlet Majorana coupling has charge2(3q-h).
    assert sp.solve(2*charge['nu_conjugate'],h)==[3*q]
    sm={k:str(v.subs(h,3*q).subs(q,sp.Rational(1,6)))for k,v in charge.items()}
    chars=[(2,0),(2,2),(1,1)]
    mask=np.array([[all((a+b)%3==0 for a,b in zip(x,y))for y in chars]for x in chars],int)
    assert np.linalg.matrix_rank(mask)==2
    rng=np.random.default_rng(11409);rand=lambda shape:rng.normal(size=shape)+1j*rng.normal(size=shape)
    Yf=rand((3,3));Yq=C@np.kron(Yf,np.eye(3))@C.conj().T;Yl=L@rand((3,3))@L.conj().T
    Q=C@rand((9,2))@W.conj().T;Ll=L@rand((3,2))@W.conj().T
    U=Ru@rand((1,9))@C.conj().T;D=Rd@rand((1,9))@C.conj().T
    E=Rd@rand((1,3))@L.conj().T;N=Ru@rand((1,3))@L.conj().T
    hv=rand((2,1));eps=np.array([[0,1],[-1,0]])
    Hu=W@hv@Ru.conj().T;Hd=W@(eps@hv.conjugate())@Rd.conj().T
    invariants=[(Hu,U,Yq,Q),(Hd,D,Yq,Q),(Hd,E,Yl,Ll),(Hu,N,Yl,Ll)];res=[]
    for qv,hv0 in [(1/6,.5),(.2,.7)]:
        GY=qv*C@C.conj().T-3*qv*L@L.conj().T-hv0*Ru@Ru.conj().T+hv0*Rd@Rd.conj().T
        for G in col+weak+[GY]:
            for A,B,Y,K in invariants:
                comm=lambda X:G@X-X@G
                variation=np.trace(comm(A)@B@Y@K+A@comm(B)@Y@K+A@B@Y@comm(K))
                res.append(abs(variation));assert abs(variation)<1e-10
    # SU2 pseudoreal conjugacy is checked separately on full finite transformations.
    for T in [W.conj().T@G@W for G in weak]:
        Ug=expm(.23j*T);assert np.linalg.norm(eps@Ug.conjugate()-Ug@eps)<1e-12
    Yu=Yf;Yd=rand((3,3));Ye=rand((3,3));Yn=rand((3,3));M=block_diag(np.kron(Yu,np.eye(3)),np.kron(Yd,np.eye(3)),Ye,Yn)
    DF=np.block([[np.zeros((24,24)),M],[M.conj().T,np.zeros((24,24))]]);gamma=np.diag([1]*24+[-1]*24)
    assert np.linalg.norm(DF@gamma+gamma@DF)<1e-12
    vals=np.linalg.eigvalsh(DF);assert np.linalg.norm(vals+vals[::-1])<1e-12
    return dict(status='PASS',actual_gauge_variation_maximum=float(max(res)),Yukawa_terms=['Tr(Hu U^c Yq Q)','Tr(Hd D^c Yq Q)','Tr(Hd E^c Yl L)','Tr(Hu N^c Yl L)'],
        Higgs_relation='Hd=W epsilon h* Rd^dag, Hu=W h Ru^dag; epsilon=[[0,1],[-1,0]]',
        bare_Majorana_constraint='2(3q-h)=0 => h=3q for a nonzero uncharged coefficient',conditional_hypercharges=sm,
        native_neutrino_characters=chars,H27_invariant_symmetric_Majorana_mask=mask.tolist(),Majorana_singular_values=np.linalg.svd(mask,compute_uv=False).tolist(),Majorana_rank=2,
        neutral_Higgs_Dirac_block=cm(M),selfadjoint_Dirac_spectrum=vals.tolist(),balanced_internal_grading_index=0,
        scope='Actual native rectangular Yukawa compositions are gauge invariant. Nonzero bare singlet Majorana mass removes the anomaly-free B-L charge ambiguity, an established mechanism; allowing a charged Majorana scalar changes this conclusion. Unbroken native character assignments allow only a rank2 Majorana pair. Lorentzian Weyl roles, source orientation, coefficients and mirror selection remain supplied; selfadjoint balanced internal grading does not select a chiral spacetime theory.')

def messenger_locality(c):
    lines=c['lines'];A=np.array([[int(i!=j and bool(set(x)&set(y)))for j,y in enumerate(lines)]for i,x in enumerate(lines)],dtype=int)
    assert np.all(A.sum(axis=1)==12)and np.array_equal(A@A,8*np.eye(40,dtype=int)-2*A+4*np.ones((40,40),int))
    indices=[0,int(np.where(A[:,0])[0][0]),int(np.where((A[:,0]==0)&(np.arange(40)!=0))[0][0])];dist=[0,1,2]
    P=np.zeros((40,40));P[0,0]=1;rows=[]
    for eta in (.01,.003,.001,.0003):
        K=np.eye(40)-eta*A;R=np.linalg.inv(K);lam=.1
        M=np.block([[K,lam*P],[np.zeros((40,40)),K]])
        inv=np.linalg.inv(M);cross=inv[:40,40:];expected=-lam*R@P@R
        assert np.linalg.norm(cross-expected)<1e-12
        masses=-np.diag(cross)[indices]/lam
        den=1+2*eta-8*eta**2;alpha=(1+2*eta)/den;beta=eta/den;gamma0=4*eta**2/(den*(1-12*eta))
        closed=np.array([alpha+gamma0,beta+gamma0,gamma0])**2
        assert np.linalg.norm(masses-closed)<1e-13
        assert np.all(masses>0)
        # Norm tail of Neumann expansion, native maxdegree12.
        tail=(12*eta)**7/(1-12*eta);partial=sum(eta**n*np.linalg.matrix_power(A,n)for n in range(7))
        residual=np.linalg.norm(R-partial,ord=2);assert residual<=tail+1e-14
        rows.append(dict(eta=eta,endpoint_mass_over_lambda=masses.tolist(),near_over_eta2=float(masses[1]/eta**2),far_over_16eta4=float(masses[2]/(16*eta**4)),Neumann_tail_norm=float(residual),Neumann_tail_bound=tail,floating_replay_allowance=1e-14,heavy_min_singular_value=float(np.linalg.svd(M,compute_uv=False)[-1])))
    assert abs(rows[-1]['near_over_eta2']-1)<.02 and abs(rows[-1]['far_over_16eta4']-1)<.02
    # Link spurions transform independently at every vertex of both networks.
    phases=np.random.default_rng(11410).uniform(-np.pi,np.pi,(2,40));Ga,Gb=[np.diag(np.exp(1j*t))for t in phases]
    K=np.eye(40)-.01*A;T=.1*P;Ra=np.linalg.inv(Ga@K@Ga.conj().T);Rb=np.linalg.inv(Gb@K@Gb.conj().T)
    transformed=Ra@(Ga@T@Gb.conj().T)@Rb;plain=np.linalg.inv(K)@T@np.linalg.inv(K)
    assert np.linalg.norm(transformed-Ga@plain@Gb.conj().T)<1e-12
    shortcut=P.copy();shortcut[indices[2],indices[2]]=.001
    shortcut_mass=float((np.linalg.inv(K)@shortcut@np.linalg.inv(K))[indices[2],indices[2]])
    # Crucial compatibility audit: three attachments to ONE pin do not make
    # three independent masses. The complete family matrix is rankone.
    R=np.linalg.inv(K);common_family=.1*np.outer(R[indices,0],R[0,indices]);assert np.linalg.matrix_rank(common_family,tol=1e-12)==1
    # A repair needs three explicitly separated messenger copies, each with
    # its family attached at its own selected distance. This is an input.
    M=np.block([[K,.1*P],[np.zeros((40,40)),K]])
    replicated=np.kron(np.eye(3),M);attachment=np.zeros((240,3));extraction=np.zeros((3,240))
    for j,i in enumerate(indices):attachment[80*j+40+i,j]=1;extraction[j,80*j+i]=1
    corrected=-extraction@np.linalg.solve(replicated,attachment)
    assert np.linalg.norm(corrected-np.diag(.1*R[indices,0]**2))<1e-12 and np.linalg.matrix_rank(corrected,tol=1e-12)==3
    first=[]
    for i,d in zip(indices,dist):
        coefficients=[int(np.linalg.matrix_power(A,n)[i,0])for n in range(5)]
        assert next(n for n,x in enumerate(coefficients)if x)==d;first.append(coefficients)
    return dict(status='PASS',native_adjacency=A.tolist(),endpoint_indices=indices,distances=dist,path_count_coefficients=first,scans=rows,
        messenger_matrix='M=[[J I-h A,lambda Ppin],[0,J I-h A]]; two separate vectorlike networks; endpoint left/right couplings gL,gR are supplied',
        exact_zero_momentum_matching='m_i=gL*gR*lambda [(J I-h A)^(-1)]_(i,pin)^2 (up to removable sign)',
        spurion_group='independent U(1) at all40 sites of each network; link couplings transform bifundamentally, the sole cross-network source is at the pin',
        minimum_link_spurion_degrees=[0,2,4],leading_mass_coefficients=[1,1,16],
        exact_native_resolvent='R=alpha I+beta A+gamma 11^T; den=1+2eta-8eta^2, alpha=(1+2eta)/den, beta=eta/den, gamma=4eta^2/[den(1-12eta)]',
        analytic_protection='Every covariant analytic cross-network mass at endpointi requires a path to the pin in each network: at least2*distance link spurions. Loops and local counterterms respect this power if this enlarged microscopic selection rule and analyticity hold.',
        shortcut_far_mass_control=shortcut_mass,
        common_pin_three_family_mass=cm(common_family),common_pin_family_rank=1,
        repaired_three_copy_family_mass=cm(corrected),repaired_family_rank=3,repaired_heavy_Dirac_block_dimension=240,
        repair='Three independent copies of the two-network messenger architecture, with a conserved family-copy tag and one endpoint per copy; family replication is supplied.',
        joint_alignment_boundary='The separate competing-clock action is colour compatible, but its simultaneous microscopic compatibility with the family-copy/link spurion rules and all loops is not constructed.',
        scope='A constructed microscopic locality/spurion architecture protects hop degree conditionally. It is not protection from the previous residual H27 clock alone. One common pin gives rankone family matching, despite three hierarchical diagonal entries; three conserved messenger copies provide a full-rank conditional repair. Networks, charges, family-copy tags, pin bridges, messenger scale and attachments are supplied. Finite zero-momentum tree matching is computed, not a full gauge/scalar loop renormalization. A direct extra endpoint bridge defeats protection. Pin-local vacuum diagrams have no hopping requirement and vacuum energy remains unprotected.')

def coherent_seams(c):
    # Primitive FCC Cartesian lattice: x+y+z even. The native metric is
    # isometric to this prior11389 realization; orientation here is supplied.
    F=sp.Matrix([[1,1,0],[1,0,1],[0,1,1]])
    F_native=sp.Matrix(c['x']['line']['fcc_roots']).applyfunc(sp.Rational)
    basis_change=F_native.inv()*F
    assert abs(basis_change.det())==1 and all(x.q==1 for x in basis_change)
    R=sp.Matrix([[sp.Rational(3,5),-sp.Rational(4,5),0],[sp.Rational(4,5),sp.Rational(3,5),0],[0,0,1]])
    Rx=sp.Matrix([[1,0,0],[0,sp.Rational(3,5),-sp.Rational(4,5)],[0,sp.Rational(4,5),sp.Rational(3,5)]])
    rational_period=F.inv()*R.inv()*F
    # Exact congruence enumeration mod5 constructs the coincidence sublattice.
    generators=[]
    for t in product(range(5),repeat=3):
        v=sp.Matrix(t)
        if all(x.q==1 for x in rational_period*v):generators.append(v)
    generators.extend([5*sp.eye(3)[:,i]for i in range(3)])
    M=hermite_normal_form(sp.Matrix.hstack(*generators));N=rational_period*M
    assert abs(M.det())==5 and all(x.q==1 for x in N)and F*M==R*F*N
    # Common boundary plane z=0 has a rank2 coincidence lattice.
    plane=[F*M*sp.Matrix(t)for t in product(range(-2,3),repeat=3)if (F*M*sp.Matrix(t))[2]==0]
    assert sp.Matrix.hstack(*plane).rank()==2
    # The earlier pi/4 orientation has only the common z-axis, not a planar seam.
    R45=sp.Matrix([[sp.sqrt(2)/2,-sp.sqrt(2)/2,0],[sp.sqrt(2)/2,sp.sqrt(2)/2,0],[0,0,1]])
    # Irrationality forces x=y=0 if both x,y and their rotated coordinates are rational.
    mean=sp.Rational(3,7)*sp.eye(3)+sp.Rational(4,7)*sp.Matrix([[0,-1,0],[1,0,0],[0,0,1]])
    assert mean==R*sp.diag(sp.Rational(5,7),sp.Rational(5,7),1)
    V=sp.Matrix(c['x']['line']['fcc_harmonic_displacements']).applyfunc(sp.Rational)
    tensors=[]
    for x in V.tolist():
        v=sp.Matrix(x);S=v*v.T
        if S not in tensors:tensors.append(S)
    flat=lambda S:sp.Matrix([S[0,0],S[1,1],S[2,2],S[0,1],S[0,2],S[1,2]])
    rank=sp.Matrix.hstack(*[flat(Q*S*Q.T)for Q in [sp.eye(3),R,Rx]for S in tensors]).rank();assert rank==6
    # Nonlinear scalar-transport bracket with central difference along a native
    # primitive period. Prior11301 owns the general finite Leibniz obstruction.
    defects=[]
    for e in (.2,.1,.05,.025):
        k,l,p=e,2*e,3*e
        exact=-np.sin(p)*(np.sin(p+l)-np.sin(p+k))
        continuum_form=-np.sin(p)*(np.sin(l)-np.sin(k))
        defects.append(dict(epsilon=e,closure_defect=float(exact-continuum_form),defect_over_epsilon4=float((exact-continuum_form)/e**4)))
    assert abs(defects[-1]['defect_over_epsilon4']-27)<.2
    # Exact discrete plane-wave replay avoids interpreting a scalar norm as gravity.
    n=64;x=2*np.pi*np.arange(n)/n;D=(np.roll(np.eye(n),-1,axis=0)-np.roll(np.eye(n),1,axis=0))/2
    xi=np.exp(2j*x);eta=np.exp(3j*x);f=np.exp(4j*x)
    T=lambda a:np.diag(a)@D;lhs=(T(xi)@T(eta)-T(eta)@T(xi))@f;rhs=T(xi*(D@eta)-eta*(D@xi))@f
    closure=float(np.linalg.norm(lhs-rhs)/np.sqrt(n));assert closure>1e-5
    return dict(status='PASS',FCC_primitive=[[int(x)for x in row]for row in F.tolist()],native_to_second_FCC_basis=[[int(x)for x in row]for row in basis_change.tolist()],rational_rotation=[[str(x)for x in row]for row in R.tolist()],polar_mean_native_weights=dict(identity=3,quarter_turn=4),
        coincidence_left_periods=[[int(x)for x in row]for row in M.tolist()],coincidence_right_periods=[[int(x)for x in row]for row in N.tolist()],coincidence_index=5,planar_common_translation_rank=2,irrational_pi4_common_translation_rank=1,mixed_grain_strain_rank=rank,
        nonlinear_transport_scans=defects,discrete_plane_wave_bracket_defect=closure,
        scope='Exact common enlarged-period seam maps join two supplied FCC grain realizations, retaining all six strain directions across three orientations. This supplies translation matching, not a canonical microscopic coupling of mismatched boundary vertices. Local interpolation, seam energy and nonlinear gravitational action remain inputs. Central-difference transport fails the undeformed nonlinear bracket at finite spacing; prior11301 owns the abstract no-go, this is an explicit native-period refinement witness. O(epsilon4) closure defect does not prove nonlinear Einstein constraints.')

def invariant_star_frame(Hp,Hl,e):
    basis=[e.astype(complex)];i=0
    while i<len(basis):
        for H in [Hp,Hl]:
            v=H@basis[i]
            for _ in range(2):
                for b in basis:v-=b*np.vdot(b,v)
            if np.linalg.norm(v)>1e-10:basis.append(v/np.linalg.norm(v))
        i+=1
    Q=np.column_stack(basis)
    assert max(np.linalg.norm(H@Q-Q@(Q.conj().T@H@Q))for H in [Hp,Hl])<1e-10
    return Q

def unitary_calibration(S,axis,angle):
    T,V=schur(S,output='complex');assert np.linalg.norm(T-np.diag(np.diag(T)))<1e-10
    angles=-np.angle(np.diag(T));overlaps=abs(V.conj().T@axis)**2;j=int(np.argmax(overlaps))
    assert overlaps[j]>.95 and abs(angles[j])>1e-10
    ticks=int(round(angle/angles[j]));assert ticks>0
    alpha=np.sqrt(max(0.,1-overlaps[j]));bound=4*alpha+abs(ticks*angles[j]-angle)
    return angles,V,ticks,float(bound),float(1-overlaps[j])

def finite_ticks(c,p):
    D=np.array(c['D'],float);Hp=D[:40].T@D[:40];Hl=D[40:].T@D[40:]
    Q=read(p['local_gates']['logical_basis']);a,b=p['local_gates']['gate_pair'];ideal=read(p['local_gates']['ideal_qubit_gate']);phase=p['local_gates']['five_pulse_phase_angles']
    proj=np.array(json.loads((ROOT/'data/w33_pass11395_11402_joint_atlas_spin_portal.json').read_text())['axis_bridge']['cycle_projector_160_numerator'],dtype=float)/160
    strengths=[.025,.0125];rows=[];phases=[]
    for dt in (.04,.02,.01):
        Up=np.eye(160)+(np.exp(-2j*dt)-1)*Hp/4
        Ul=np.eye(160)+(np.exp(-4j*dt)-1)*Hl/4;clock=Up@Ul@Up
        assert np.linalg.norm(clock@proj-proj)<1e-12
        for u in strengths:
            state=Q.copy();tickcounts=[];bound=0;rounding=[]
            for angle,e in zip(phase,[a,b,a,b,a]):
                amp=np.copysign(u,angle);kick=np.ones(160,complex);kick[e]=np.exp(-.5j*dt*amp)
                S=kick[:,None]*clock*kick[None,:];axis=proj[:,e]/np.sqrt(proj[e,e])
                ang,V,N,B,dress=unitary_calibration(S,axis,angle)
                state=V@(np.exp(-1j*N*ang)[:,None]*(V.conj().T@state));tickcounts.append(N);bound+=B
                rounding.append(dict(ticks=N,dressed_phase_per_tick=float(ang[int(np.argmax(abs(V.conj().T@axis)**2))]),dressing_probability=dress))
            err=float(np.linalg.norm(state-Q@ideal,ord=2));leak=float(np.linalg.norm((np.eye(160)-proj)@state,ord=2))
            assert err<bound+1e-9
            rows.append(dict(dt=dt,strength=u,ticks_per_pulse=tickcounts,total_ticks=sum(tickcounts),total_encoded_error=err,massive_leakage=leak,dressing_and_rounding_bound=bound,calibration=rounding))
            # The pi/8 gate is diag(1,exp(i*pi/4)) up to global phase.
            kick=np.ones(160,complex);kick[a]=np.exp(-.5j*dt*u);S=kick[:,None]*clock*kick[None,:]
            axis=proj[:,a]/np.sqrt(proj[a,a]);ang,V,N,bound,dress=unitary_calibration(S,axis,np.pi/4)
            state=V@(np.exp(-1j*N*ang)[:,None]*(V.conj().T@Q));target=Q@np.diag([np.exp(-1j*np.pi/4),1])
            error=float(np.linalg.norm(state-target,ord=2));assert error<=bound+1e-9
            phases.append(dict(dt=dt,strength=u,ticks=N,total_encoded_error=error,dressing_and_rounding_bound=bound))
    # The two-cell discrete walk has a small EXACT joint star invariant frame.
    edge=np.eye(160)[:,a];B=invariant_star_frame(Hp,Hl,edge);Kp=B.conj().T@Hp@B;Kl=B.conj().T@Hl@B;f=B.conj().T@edge
    v=B.conj().T@(proj@edge/np.sqrt(proj[a,a]));delta=np.kron(f,f);logical=np.kron(v,v);ent=[]
    for dt in (.04,.02,.01):
        clock=expm(-.5j*dt*Kp)@expm(-1j*dt*Kl)@expm(-.5j*dt*Kp)
        U0=np.kron(clock,clock)
        for g in strengths:
            kick=np.eye(len(delta))+(np.exp(-.5j*dt*g)-1)*np.outer(delta,delta.conj())
            S=kick@U0@kick;ang,V,N,bound,dress=unitary_calibration(S,logical,np.pi)
            state=V@(np.exp(-1j*N*ang)*(V.conj().T@logical));amp=np.vdot(logical,state);err=float(np.linalg.norm(state+logical));leak=float(np.sqrt(max(0.,1-abs(amp)**2)))
            assert err<bound+1e-9
            ent.append(dict(dt=dt,coupling=g,ticks=N,total_encoded_error=err,massive_leakage=leak,logical_00_amplitude=[float(amp.real),float(amp.imag)],dressing_and_rounding_bound=bound,logical_plus_plus_survival=float((abs(amp)**2+3)/4),conditional_concurrence=float(2*abs(amp-1)/(abs(amp)**2+3))))
    return dict(status='PASS',star_split='S_e(dt)=edge_half_kick * point_half_star * line_full_star * point_half_star * edge_half_kick; star phases2dt,4dt,2dt',star_layers_per_tick=3,edge_kick_layers_per_tick=2,Hadamard_scans=rows,pi_over_8_gate_scans=phases,entangler_joint_star_dimension=B.shape[1],two_cell_reduced_dimension=B.shape[1]**2,entangler_star_frame=cm(B),entangler_scans=ent,
        exact_error_bound='Per pulse4 sqrt(1-|<dressed Floquet eigenvector|protected axis>|²)+|N*theta-target_angle|; telescope over pulses. Exact unitary powers, no long-time Trotter accumulation estimate used.',
        adversarial_noise_bound='Additional operator error <= total_ticks * epsilon_per_tick; no fault tolerance or noise threshold is implied.',
        scope='Full counted-tick native star primitives realize finite-penalty Hadamard and entangling phases with exact calibrated Floquet shifts and explicit leakage/rounding bounds. Tick phases, strength, addressing, state preparation and intercell density kicks are supplied. Star-frame reduction is exactly invariant; physical hardware and fault tolerance remain open.')

def produce():
    p=prior();c=load();parts={}
    for name,fn in [('competing_clocks',lambda:competing_clocks(p)),('yukawa_arrows',lambda:yukawa_arrows(p)),('messenger_locality',lambda:messenger_locality(c)),('coherent_seams',lambda:coherent_seams(c)),('finite_ticks',lambda:finite_ticks(c,p))]:
        parts[name]=fn();print(name,'PASS',flush=True)
    source='data/w33_pass11403_11407_native_bimodule_collective_gates.json';hashval=hashlib.sha256(json.dumps(p,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    out=dict(status='PASS',passes='11408-11412',reservation='383a35691',source_sha256={source:hashval},source_hash_convention='canonical sorted compact JSON',**parts,
        scope='Five concrete conditional follow-ups. No measured masses/mixing, selected physical chirality, emergent nonlinear Einstein theory, protected vacuum energy or fault-tolerant quantum machine.')
    OUT.write_text(json.dumps(out,indent=2)+'\n');return out
if __name__=='__main__':produce()
