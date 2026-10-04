#!/usr/bin/env python3
"""Five concrete follow-ups: native bimodule, spectral action, collective geometry,
conditional radiative hierarchy and physically local finite-penalty gates.
Prior owners BT1041/1042/1058, Pass4091, BT4048, Pass11316 and11389-11402.
"""
from __future__ import annotations
import json
import hashlib
from pathlib import Path
from itertools import product,combinations
import numpy as np
import sympy as sp
from scipy.linalg import expm
from scipy.optimize import least_squares
from w33_pass11390_11394_context_matter_clock import load
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/w33_pass11403_11407_native_bimodule_collective_gates.json'

def cm(A):return dict(real=np.asarray(A).real.tolist(),imag=np.asarray(A).imag.tolist())
def exact(A):return [[str(v)for v in A.row(i)]for i in range(A.rows)]
def orth(P,rank):
    cols=[]
    for v in P.T:
        v=v.copy().astype(complex)
        for u in cols:v-=u*np.vdot(u,v)
        if np.linalg.norm(v)>1e-9:
            v/=np.linalg.norm(v);i=np.where(abs(v)>1e-9)[0][0];v*=np.exp(-1j*np.angle(v[i]));cols.append(v)
        if len(cols)==rank:break
    F=np.array(cols).T
    assert F.shape[1]==rank and np.linalg.norm(F.conj().T@F-np.eye(rank))<1e-10
    return F

def extended(c):
    prior=json.loads((ROOT/'data/w33_pass11395_11402_joint_atlas_spin_portal.json').read_text())['chiral']['spread_extension']
    spreads=prior['spreads'];B=np.zeros((40,36))
    for j,S in enumerate(spreads):B[S,j]=1
    C=np.zeros((80,80))
    for a,b in json.loads((ROOT/'data/w33_pass11390_11394_context_matter_clock.json').read_text())['portal']['nonedge_orbit']:C[a,b]=C[b,a]=1
    H=np.zeros((116,116));H[:80,:80]=c['A']+.1*C;H[40:80,80:]=B;H[80:,40:80]=B.T
    si={tuple(S):i for i,S in enumerate(spreads)};gs=[]
    for g in c['x']['internal_h27']['group_permutations']:
        gp=list(g)+[80+si[tuple(sorted(g[40+l]-40 for l in S))]for S in spreads]
        gs.append(np.eye(116)[np.argsort(gp)])
    PN=sum(gs)/27;w,v=np.linalg.eigh(np.eye(116)-PN);T=v[:,w>.5]
    w,v=np.linalg.eigh(T.T@H@T);K=T@v[:,abs(w)<1e-9];assert K.shape[1]==34
    return H,gs,K@K.T

def gellmann():
    T=[]
    for i,j in combinations(range(3),2):
        X=np.zeros((3,3),complex);X[i,j]=X[j,i]=.5;T.append(X)
        X=np.zeros((3,3),complex);X[i,j]=-.5j;X[j,i]=.5j;T.append(X)
    T.extend([np.diag([1,-1,0])/2,np.diag([1,1,-2])/(2*np.sqrt(3))]);return T

def fermion_bimodule(c):
    H,gs,K=extended(c);omega=np.exp(2j*np.pi/3)
    def linear(u,v):return K@sum(w*omega**(-(u*a+v*b))for w,(a,b,z)in zip(gs,product(range(3),repeat=3)))/27
    Pc=K@(np.eye(116)+omega.conjugate()*gs[1]+omega*gs[2])/3
    seed=Pc@(np.eye(116)+gs[3]+gs[3]@gs[3])/3
    F0=orth(seed,3);Fc=[np.linalg.matrix_power(gs[9],j)@F0 for j in range(3)];colour=np.hstack(Fc)
    assert np.linalg.norm(colour.conj().T@colour-np.eye(9))<1e-10
    weak=orth(linear(1,0),2)
    lepton=np.hstack([orth(linear(*uv),2)[:,:1]for uv in [(1,1),(1,2),(2,0)]])
    r1=orth(linear(0,1),2)[:,:1];r2=orth(linear(2,2),2)[:,:1]
    allslots=np.hstack([colour,weak,lepton,r1,r2]);assert np.linalg.norm(allslots.conj().T@allslots-np.eye(16))<1e-10
    colour_g=[sum(F@T@F.conj().T for F in Fc)for T in gellmann()]
    pauli=[np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1])]
    weak_g=[weak@T@weak.conj().T/2 for T in pauli]
    # The two Higgs arrows have conjugate H27 characters, allowing one
    # weak doublet and its pseudoreal conjugate at representation level.
    weak_char=(1,0);up_char=(0,1);down_char=(2,2)
    assert tuple((a+b)%3 for a,b in zip(up_char,down_char))==tuple(2*a%3 for a in weak_char)
    right_pairs=[(u,tuple((2*w-x)%3 for w,x in zip(weak_char,u)))for u in product(range(3),repeat=2)if u not in [(0,0),weak_char,tuple(2*x%3 for x in weak_char)]]
    assert len(right_pairs)==6
    maxcomm=max(np.linalg.norm(g@T-T@g)for g in gs for T in colour_g+weak_g);assert maxcomm<1e-9
    family=np.column_stack([A[:,0]for A in Fc]);nativeX=family.conj().T@gs[9]@family;nativeZ=family.conj().T@gs[3]@family
    X=sp.Matrix([[0,0,1],[1,0,0],[0,1,0]]);zeta=(-1+sp.I*sp.sqrt(3))/2;Z=sp.diag(1,zeta,zeta**2)
    assert np.linalg.norm(nativeX-np.array(X,dtype=complex))<1e-10
    assert min(np.linalg.norm(nativeZ-np.array(Z,dtype=complex)),np.linalg.norm(nativeZ-np.array(Z.conjugate(),dtype=complex)))<1e-10
    eq=sp.Matrix.vstack(sp.kronecker_product(sp.eye(3),X)-sp.kronecker_product(X.T,sp.eye(3)),sp.kronecker_product(sp.eye(3),Z)-sp.kronecker_product(Z.T,sp.eye(3)))
    assert eq.rank()==8
    q,h=sp.symbols('q h');charges={'Q':q,'u_conjugate':-q-h,'d_conjugate':-q+h,'L':-3*q,'e_conjugate':3*q+h,'nu_conjugate':3*q-h}
    fields={'Q':(colour,weak),'u_conjugate':(r1,colour),'d_conjugate':(r2,colour),'L':(lepton,weak),'e_conjugate':(r2,lepton),'nu_conjugate':(r1,lepton)}
    dims={name:A.shape[1]*B.shape[1]for name,(A,B)in fields.items()};assert dims=={'Q':18,'u_conjugate':9,'d_conjugate':9,'L':6,'e_conjugate':3,'nu_conjugate':3}
    Q,U,D,L,E,N=[charges[name]for name in fields]
    anomalies=[3*(2*Q+U+D),3*(3*Q+L),sum(dims[k]*charges[k]for k in fields),sum(dims[k]*charges[k]**3 for k in fields)]
    assert all(sp.expand(x)==0 for x in anomalies)
    # Gauge action on actual native rectangular maps, not a dimension match.
    Y=colour@colour.conj().T/6-lepton@lepton.conj().T/2-r1@r1.conj().T/2+r2@r2.conj().T/2
    for name,(A,B)in fields.items():
        X=A[:,0:1]@B[:,0:1].conj().T;y=float(charges[name].subs({q:sp.Rational(1,6),h:sp.Rational(1,2)}))
        assert np.linalg.norm(Y@X-X@Y-y*X)<1e-10
    # Casimir on Q: colour fundamental and weak antifundamental (equivalent doublet).
    Qcas=sum(np.kron((colour.conj().T@T@colour),np.eye(2))@np.kron((colour.conj().T@T@colour),np.eye(2))for T in colour_g)
    Wcas=sum(np.kron(np.eye(9),-(weak.conj().T@T@weak).T)@np.kron(np.eye(9),-(weak.conj().T@T@weak).T)for T in weak_g)
    assert np.linalg.norm(Qcas-4*np.eye(18)/3)<1e-9 and np.linalg.norm(Wcas-3*np.eye(18)/4)<1e-9
    return dict(native_slot_frames={'colour_qutrit_seeds':[cm(F)for F in Fc],'weak':cm(weak),'lepton':cm(lepton),'right_up':cm(r1),'right_down':cm(r2)},
        field_rectangles={k:[A.shape[1],B.shape[1]]for k,(A,B)in fields.items()},field_dimensions=dims,
        family_clock_matrices={'X':cm(nativeX),'Z':cm(nativeZ)},exact_family_commutant_dimension=9-eq.rank(),left_Weyl_bimodule_dimension=48,adjoint_doubling_dimension=96,declared_families=3,
        exact_anomalies=[str(sp.expand(x))for x in anomalies],witten_doublets=12,
        charge_family={k:str(v)for k,v in charges.items()},native_gauge_H27_commutator_residual=float(maxcomm),
        single_Higgs_character_pairing={'weak':weak_char,'right_up':up_char,'right_down':down_char,'allowed_ordered_right_pairs':right_pairs},
        unbroken_H27_one_character_Higgs_quark_singular_values='three equal by Schur; a linear character twist of an irreducible qutrit admits only a scalar times unitary intertwiner',
        gauge_slots='M3 colour, M2 weak, C lepton, C right-up, C right-down embedded in the34D protected native kernel',
        scope='Actual rectangular operators in End(native kernel) realize three declared SM+nu representation copies and their adjoints. The bare34 states cannot carry a commuting irreducible(3,2); the operator bimodule can. Slot marks, right rankone choices, two charge parameters, Weyl assignment and selection against mirrors are inputs. No physical chirality, unique hypercharge, Yukawa fit or gauge dynamics derived. BT1041/1042 own generic left/right bimodules; Pass4091 and11398 own generic anomaly family.')

def spectral_selector(c):
    prior=json.loads((ROOT/'data/w33_pass11390_11394_context_matter_clock.json').read_text())['portal'];C=np.zeros((80,80))
    for a,b in prior['nonedge_orbit']:C[a,b]=C[b,a]=1
    H=c['A']+C;g=c['x']['internal_h27']['group_permutations'];Z=np.eye(80)[np.argsort(g[1])];Pc=(2*np.eye(80)-Z-Z@Z)/3
    xs=np.sort(np.roots([1,-15,54,-27]).real);assert np.all(xs>0)
    Ns=[]
    for j,x in enumerate(xs):
        f=np.eye(80)
        for k,y in enumerate(xs):
            if k!=j:f=f@(H@H-y*np.eye(80))/(x-y)
        N=Pc@H@f/np.sqrt(x);assert np.linalg.norm(N-N.T)<1e-10;Ns.append(N)
    for i,N in enumerate(Ns):
        assert np.linalg.norm(N@N@N-N)<1e-9
        for j,M in enumerate(Ns):
            if i!=j:assert np.linalg.norm(N@M)<1e-9
    blocks=[c['plus'].conj().T@c['S'].conj().T@N@c['S']@c['minus']for N in Ns]
    for B in blocks:assert np.linalg.matrix_rank(B,tol=1e-9)==1 and abs(np.linalg.norm(B)-1)<1e-10
    old=json.loads((ROOT/'data/w33_pass11395_11402_joint_atlas_spin_portal.json').read_text())['selector_portal']['three_channel_compiler']
    # Old SVD order is descending; branch rule here is algebraic ascending.
    oldblocks=[]
    for orb in old['orbit_pairs']:
        M=np.zeros((80,80))
        for a,b in orb:M[a,b]=M[b,a]=1
        oldblocks.append(c['plus'].conj().T@c['S'].conj().T@M@c['S']@c['minus'])
    residuals=[]
    for j,weights in enumerate(old['rank1_channel_weights']):
        B=sum(w*M for w,M in zip(weights,oldblocks));residuals.append(float(np.linalg.norm(B-blocks[2-j])))
    assert max(residuals)<1e-10
    return dict(reference_squared_roots=xs.tolist(),reference_polynomial='x^3-15x^2+54x-27',
        operator_rule='N_j=Pcentral Href product_(k!=j)(Href^2-x_k I)/(x_j-x_k)/sqrt(x_j)',
        operators=[cm(N)for N in Ns],central_blocks=[cm(B)for B in blocks],old_compiler_residuals=residuals,
        supplied_joint_action='Hsel tensor I + lambda sum_j D_j tensor N_(2-j)',
        shell_endpoint_weights=['pin projector','sum_neighbour_projectors/12','sum_disjoint_projectors/27'],
        reference_input='Href=actual native incidence + actual prior H27 nonedge orbit, both coefficient1',
        scope='Spectral functional calculus fixes the three channel orientations and amplitudes from the named native reference operator; no SVD frame or88 independent orbital weights are inputs. The marked orbit, pin, common coupling and conditional endpoint architecture remain supplied. This is an explicit Hamiltonian ansatz, not a variational proof that nature chooses it.')

def colour_safe_bridge(c,matter,spectral):
    read=lambda d:np.array(d['real'])+1j*np.array(d['imag'])
    Fc=[read(F)for F in matter['native_slot_frames']['colour_qutrit_seeds']];F=np.hstack(Fc)
    Gs=[sum(A@T@A.conj().T for A in Fc)for T in gellmann()]
    Ps=[A@A.conj().T for A in Fc];unsafe=[]
    for d in spectral['operators']:
        N=np.zeros((116,116),complex);N[:80,:80]=read(d)
        unsafe.append(float(max(np.linalg.norm(N@G-G@N)for G in Gs)))
    assert min(unsafe)>.1
    assert max(np.linalg.norm(P@G-G@P)for P in Ps for G in Gs)<1e-10
    H,group,K=extended(c);Z=group[3];X=group[9]
    assert max(np.linalg.norm(P@Z-Z@P)for P in Ps)<1e-10
    assert np.linalg.norm(Ps[0]@X-X@Ps[0])>1
    # Family basis is the native qutrit clock; colour remains a spectator.
    family=np.column_stack([A[:,0]for A in Fc]);xf=family.conj().T@X@family
    ew,V=np.linalg.eig(xf);order=np.argsort(np.angle(ew));V=V[:,order]
    for j in range(3):V[:,j]/=np.linalg.norm(V[:,j])
    assert np.linalg.norm(V.conj().T@V-np.eye(3))<1e-10
    assert np.max(abs(abs(V)**2-1/3))<1e-10
    jarl=float(np.imag(V[0,0]*V[1,1]*V[0,1].conjugate()*V[1,0].conjugate()))
    assert abs(abs(jarl)-1/(6*np.sqrt(3)))<1e-10
    old=json.loads((ROOT/'data/w33_pass11395_11402_joint_atlas_spin_portal.json').read_text())['selector_portal']['three_channel_compiler'];scans=[]
    for row in old['scans']:
        phi=sum(w*P for w,P in zip(row['probabilities'],Ps));actual=np.linalg.eigvalsh(F.conj().T@phi@F)
        expected=np.sort(np.repeat(row['probabilities'],3));assert np.linalg.norm(actual-expected)<1e-12
        scans.append(dict(eta=row['eta'],family_singular_values=row['probabilities'],colour_triplet_eigenvalues=actual.tolist()))
    return dict(old_spectral_branch_colour_commutator_norms=unsafe,
        colour_preserving_native_family_projectors=[cm(P)for P in Ps],scans=scans,
        corrected_joint_action='Hsel tensor I + lambda sum_j D_shell_j tensor E_family_j, with E_family_j native Z-clock rank3 projector on colour9',
        residual_family_symmetry='central character and one cyclic qutrit clock; full H27 is broken by the nondegenerate family operator',
        same_clock_up_down_mixing='identity, J=0',other_native_clock_mixing=cm(V),
        other_clock_squared_moduli='all1/3',other_clock_absolute_J='1/(6*sqrt(3))',
        scope='The earlier80-state spectral channels fail colour commutation under the new116-state colour assignment, so they cannot be reused as quark masses unchanged. Retargeting selector probabilities to native family-clock projectors preserves colour and splits families. Native two-clock choices give either no mixing or a mutually-unbiased mixing limit, not observed CKM. Pin, clock and coupling choices remain inputs. Generic Weyl-clock Fourier mixing is prior art.')

def collective_geometry(c):
    # Rotations drawn from the previously certified native Cartesian holonomy group.
    Rx=sp.Matrix([[1,0,0],[0,0,-1],[0,1,0]]);Ry=sp.Matrix([[0,0,1],[0,1,0],[-1,0,0]]);Rz=sp.Matrix([[0,-1,0],[1,0,0],[0,0,1]])
    prior=json.loads((ROOT/'data/w33_pass11390_11394_context_matter_clock.json').read_text())['transport'];F=sp.Matrix(c['x']['line']['fcc_roots'])*sp.Matrix(c['x']['line']['primitive_fcc_change']).inv()
    group=[F*sp.Matrix(M).applyfunc(sp.Rational).inv().T*F.inv()for M in prior['holonomy_matrices']]
    assert Rx in group and Ry in group and Rz in group
    def avg(axis,n):
        a=(n-1)/np.sqrt((n-1)**2+1);b=1/np.sqrt((n-1)**2+1)
        if axis=='x':return np.array([[1,0,0],[0,a,-b],[0,b,a]])
        return np.array([[a,0,b],[0,1,0],[-b,0,a]])
    Lx=np.array([[0,0,0],[0,0,-1],[0,1,0]]);Ly=np.array([[0,0,1],[0,0,0],[-1,0,0]]);limit=Lx@Ly-Ly@Lx;scans=[]
    for n in (20,40,80,160):
        X,Y=avg('x',n),avg('y',n);W=X@Y@X.T@Y.T;a=1/n
        M=(1-1/n)*np.eye(3)+np.array(Rx,dtype=float)/n;u,_,v=np.linalg.svd(M)
        assert np.linalg.norm(u@v-X)<1e-12
        scans.append(dict(block_size=n,collective_angle=float(np.arctan(1/(n-1))),curvature_limit_error=float(np.linalg.norm((W-np.eye(3))/a**2-limit))))
    assert scans[-1]['curvature_limit_error']<scans[0]['curvature_limit_error']/5
    E=sp.Matrix(c['x']['line']['fcc_harmonic_displacements']).applyfunc(sp.Rational)
    base=[]
    for row in E.tolist():
        v=sp.Matrix(row);S=v*v.T
        if S not in base:base.append(S)
    flatten=lambda S:sp.Matrix([S[0,0],S[1,1],S[2,2],S[0,1],S[0,2],S[1,2]])
    rank=lambda ss:sp.Matrix.hstack(*[flatten(S)for S in ss]).rank()
    assert rank(base)==4
    def half(R):
        # Polar average of I and a native quarter-turn is the pi/4 rotation.
        if R==Rx:return sp.Matrix([[1,0,0],[0,sp.sqrt(2)/2,-sp.sqrt(2)/2],[0,sp.sqrt(2)/2,sp.sqrt(2)/2]])
        return sp.Matrix([[sp.sqrt(2)/2,-sp.sqrt(2)/2,0],[sp.sqrt(2)/2,sp.sqrt(2)/2,0],[0,0,1]])
    expanded=base+[U*S*U.T for U in (half(Rx),half(Rz))for S in base];assert rank(expanded)==6
    # Conditional native-momentum linearized Einstein symbol: exact Ward/null quotient.
    eta=sp.diag(-1,1,1,1);p=sp.Matrix([1,0,0,1]);pl=eta*p;p2=(p.T*eta*p)[0]
    pairs=list(combinations(range(4),2));coords=[(i,j)for i in range(4)for j in range(i,4)]
    basis=[]
    for i,j in coords:
        h=sp.zeros(4);h[i,j]=h[j,i]=1;basis.append(h)
    def einstein(h):
        tr=sp.trace(eta*h);hp=h*p;php=(p.T*h*p)[0]
        return (p2*h-pl*hp.T-hp*pl.T+pl*pl.T*tr+eta*(php-p2*tr))/2
    K=sp.Matrix.hstack(*[sp.Matrix([einstein(h)[i,j]for i,j in coords])for h in basis])
    gauges=[]
    for i in range(4):
        xi=sp.eye(4)[:,i];h=pl*xi.T+xi*pl.T;gauges.append(sp.Matrix([h[a,b]for a,b in coords]))
    G=sp.Matrix.hstack(*gauges);assert K*G==sp.zeros(10,4)and K.rank()==4 and G.rank()==4
    return dict(collective_polar_rule='polar((1-1/N)I+Rquarter/N)',scans=scans,
        commutator_curvature_limit=limit.tolist(),native_metric_strain_rank=4,collective_rotated_strain_rank=6,
        metric_generators=[exact(S)for S in expanded],conditional_null_symbol=exact(K),conditional_gauge_map=exact(G),
        conditional_massless_polarizations=10-K.rank()-G.rank(),native_low_momentum_speed_squared='27/3200',
        native_momentum_Ward_samples=fcc_ward_samples(c),
        scope='Collective polar averages enlarge the finite holonomy group into continuous SO3; individual products remain finite. A common aligned block frame and collective order parameter are inputs. This construction supplies small-curvature plaquettes and all six tensor response directions in a mixed-grain ansatz with orientations relative to a common host; no lattice seam map or dynamical coarse-graining law is constructed. A supplied Einstein quadratic action with native FCC derivative has exact Ward identity and two null polarizations (generic owner11316); neither its action nor nonlinear gravity emerges from averaging alone.')

def fcc_ward_samples(c):
    F=np.array(sp.Matrix(c['x']['line']['fcc_roots']),dtype=float);eta=np.diag([-1.,1,1,1]);coords=[(i,j)for i in range(4)for j in range(i,4)];out=[]
    for momentum in ([.05,.09,.13],[.2,.03,-.08],[.17,-.11,.06]):
        d=np.sqrt(27/3200)*np.linalg.solve(F.T,np.sin(F.T@momentum));p=np.r_[np.linalg.norm(d),d];pl=eta@p;p2=p@pl
        def E(h):
            tr=np.trace(eta@h);hp=h@p
            return (p2*h-np.outer(pl,hp)-np.outer(hp,pl)+np.outer(pl,pl)*tr+eta*(p@h@p-p2*tr))/2
        basis=[]
        for i,j in coords:
            h=np.zeros((4,4));h[i,j]=h[j,i]=1;basis.append(h)
        K=np.column_stack([np.array([E(h)[i,j]for i,j in coords])for h in basis]);G=[]
        for xi in np.eye(4):
            h=np.outer(pl,xi)+np.outer(xi,pl);G.append([h[i,j]for i,j in coords])
        G=np.array(G).T;assert np.linalg.norm(K@G)<1e-13
        assert np.linalg.matrix_rank(K,tol=1e-10)==4 and np.linalg.matrix_rank(G,tol=1e-10)==4
        out.append(dict(momentum=momentum,lattice_derivative=d.tolist(),on_shell_frequency=float(p[0]),Ward_residual=float(np.linalg.norm(K@G)),physical_polarizations=2))
    return out

def selector_radiative():
    quotient=np.array([[0,np.sqrt(12),0],[np.sqrt(12),2,6],[0,6,8]])
    D=[np.diag([1,0,0]),np.diag([0,1/12,0]),np.diag([0,0,1/27])];lam=.1;scans=[]
    def cw(m,mu):return 0. if m==0 else -3*m**4*(np.log(m*m/(mu*mu))-1.5)/(16*np.pi**2)
    # Fermion-only declared inventory: three Dirac spectators for each active branch.
    localV=np.array([cw(lam,lam),cw(lam/12,lam),cw(lam/27,lam)])
    gaps=1+localV[1:]-localV[0];assert min(gaps)>.99
    for eta in (.01,.003,.001):
        H=np.diag([0,1,1])-eta*quotient;E0=np.linalg.eigvalsh(H)[0];masses=[];shifts=[]
        for d in D:
            ep=np.linalg.eigvalsh(H+lam*d)[0];em=np.linalg.eigvalsh(H-lam*d)[0]
            masses.append((ep-em)/2);shifts.append((ep+em)/2-E0)
        Hr=np.diag([0,gaps[0],gaps[1]])-eta*quotient;_,v=np.linalg.eigh(Hr);ps=v[:,0]**2/np.array([1,12,27])
        scans.append(dict(eta=eta,finite_coupling_odd_masses=masses,even_self_energies=shifts,
            near_mass_over_lambda_eta2=float(masses[1]/(lam*eta**2)),far_mass_over_16lambda_eta4=float(masses[2]/(16*lam*eta**4)),
            loop_shifted_selector_probabilities=ps.tolist(),low_energy_CW_density=float(sum(cw(m,lam)for m in masses))))
    assert abs(scans[-1]['near_mass_over_lambda_eta2']-1/(1-(lam/12)**2))<.01
    assert abs(scans[-1]['far_mass_over_16lambda_eta4']-1/(1-(lam/27)**2))<.03
    # Exact common commutant of the three spectral branches permits direct masses.
    return dict(lambda_over_J=lam,cell_volume_in_cutoff_units=1,renormalization_scale_over_J=lam,
        CW_inventory='three supplied Dirac spectators per branch; fermion contribution only',
        local_selector_CW_energies=localV.tolist(),loop_shifted_gaps=gaps.tolist(),scans=scans,
        leading_probabilities_with_gaps=['1','h^2/Jnear^2','16h^4/(Jnear^2 Jfar^2)'],
        exact_channel_separation='N_i N_j=0 for i!=j; full selector-matter action decomposes into Hsel+lambda D_j and Hsel-lambda D_j',
        vacuum_density_eta0_over_lambda4='9/(32*pi^2) in the specified mu=lambda scheme, plus an independent vacuum counterterm',
        allowed_H27_counterterms='independent direct m_j*N_j are symmetry-allowed and need an additional microscopic locality/spurion protection argument',
        colour_safe_residual_counterterms='direct m_j*E_family_j commute with colour and the retained clock; residual symmetry alone does not protect small masses',
        scope='Exact finite-selector virtual corrections preserve hierarchy powers while gaps stay open. A specified fermion-only CW contribution shifts shell gaps and vacuum energy; cell volume converts density to selector energy. No full gauge/scalar loop inventory, UV naturalness or cosmological-constant cancellation. This separates conditional endpoint protection from H27 symmetry protection.')

def modular_entangler_clock(c,P):
    # One excitation per native cell; a supplied adjacent-cell density coupling.
    energies=np.array([0,4-np.sqrt(6),4,4+np.sqrt(6),8]);weights=np.array([81,24,30,24,1])/160
    f=np.sqrt(weights);delta=np.kron(f,f);base=np.add.outer(energies,energies).ravel();e0=np.eye(25)[:,0];scans=[]
    for coupling in (.1,.05,.025,.0125):
        H=np.diag(base)+coupling*np.outer(delta,delta);w,V=np.linalg.eigh(H);shift=w[0];t=np.pi/shift
        state=V@(np.exp(-1j*w*t)*(V.T@e0));amp=state[0];leak=np.linalg.norm(state[1:])
        err=np.linalg.norm(state+e0);survival=(abs(amp)**2+3)/4;concurrence=2*abs(amp-1)/(abs(amp)**2+3)
        scans.append(dict(coupling=coupling,phase_time=float(t),logical_00_amplitude=[float(amp.real),float(amp.imag)],massive_leakage=float(leak),total_state_error=float(err),
            plus_plus_logical_survival=float(survival),conditional_plus_plus_concurrence=float(concurrence)))
    assert scans[-1]['total_state_error']<.015 and scans[-1]['conditional_plus_plus_concurrence']>.999
    D=np.array(c['D'],dtype=float);Hp=D[:40].T@D[:40];Hl=D[40:].T@D[40:]
    assert np.linalg.norm((Hp+Hl)@P)<1e-12
    clockscans=[]
    for phi in (.1,np.pi):
        Up=np.eye(160)+(np.exp(-1j*phi)-1)*Hp/4;Ul=np.eye(160)+(np.exp(-1j*phi)-1)*Hl/4;W=Ul@Up
        assert np.linalg.norm(W.conj().T@W-np.eye(160))<1e-12 and np.linalg.norm(W@P-P)<1e-12
        angles=np.angle(np.linalg.eigvals(W));zero=sum(abs(angles)<1e-9);gap=min(abs(a)for a in angles if abs(a)>1e-9)
        clockscans.append(dict(star_phase=float(phi),unit_phase_eigenspace_dimension=int(zero),nonzero_phase_gap=float(gap)))
    assert clockscans[0]['unit_phase_eigenspace_dimension']==81 and clockscans[1]['unit_phase_eigenspace_dimension']==82
    # The extra Grover zero-phase mode is the uniform bright flag amplitude.
    uniform=np.ones(160)/np.sqrt(160);assert np.linalg.norm(P@uniform)<1e-12
    return dict(two_cell_bright_Krylov_dimension=25,two_cell_penalty='Delta(Ledge tensor I+I tensor Ledge)',
        interaction='g n_e(cellA) n_e(cellB), with one excitation per cell',
        projected_interaction='g*(81/160)^2 |00><00| on the chosen native two-flag qubits',scans=scans,
        uniform_module_gap_over_Delta='4-sqrt(6), for a tensor sum of independent cell constraints at every volume',
        star_clock_scans=clockscans,star_clock='W_phi=product_line_stars exp(-i phi Pstar) product_point_stars exp(-i phi Pstar)',
        star_gate_support=4,star_layers_per_tick=2,grover_resonance='phi=pi leaves an extra uniform bright eigenphase0; phi=.1 isolates exactly the81D cycle eigenspace',
        scope='A supplied modular architecture replaces the single connected-cover penalty with independent native cell constraints, a prepared one-excitation-per-cell invariant sector and neighbouring-cell density interactions. It has a volume-independent gap and an actual entangling phase, completing conditional local qubit universality with the Hadamard/phase controls. Two layers of local four-port stars provide a generalized native walk clock. Extra physical couplings, phase choice, initialization and cell boundaries are inputs; no hardware, noise threshold or fault-tolerant scalable implementation is proved.')

def local_gates(c):
    D=np.array(c['D'],dtype=float);L=D.T@D
    p=json.loads((ROOT/'data/w33_pass11395_11402_joint_atlas_spin_portal.json').read_text());P=np.array(p['axis_bridge']['cycle_projector_160_numerator'],dtype=float)/160
    a,b=p['controls']['edge_control_pair'];u=P[:,a]/np.sqrt(P[a,a]);v=P[:,b]/np.sqrt(P[b,b]);Q=np.column_stack([u,(v+u/3)/np.sqrt(8/9)])
    assert np.linalg.norm(Q.T@Q-np.eye(2))<1e-12
    Pe=np.diag([1.,0]);Pf=np.outer(np.array([-1/3,np.sqrt(8)/3]),np.array([-1/3,np.sqrt(8)/3]));axes=[Pe,Pf,Pe,Pf,Pe]
    target=1j*np.array([[1,1],[1,-1]])/np.sqrt(2)
    def SU(th):
        U=np.eye(2,dtype=complex)
        for angle,A in zip(th,axes):U=expm(-1j*angle*(A-np.eye(2)/2))@U
        return U
    def residual(th):
        M=SU(th)-target;return np.r_[M.real.ravel(),M.imag.ravel()]
    rng=np.random.default_rng(11407);solution=None
    for _ in range(12):
        sol=least_squares(residual,rng.uniform(-np.pi,np.pi,5),gtol=1e-13,ftol=1e-13,xtol=1e-13,max_nfev=1200)
        if np.linalg.norm(residual(sol.x))<1e-10:solution=sol.x;break
    assert solution is not None
    ideal=np.eye(2,dtype=complex)
    for angle,A in zip(solution,axes):ideal=expm(-1j*angle*A)@ideal
    scans=[]
    for strength in (.1,.05,.025,.0125):
        state=Q.astype(complex);times=[];bound=0.
        for angle,e in zip(solution,[a,b,a,b,a]):
            amp=np.copysign(strength,angle);H=L.copy();H[e,e]+=amp;w,V=np.linalg.eigh(H)
            low=w[(abs(w)>1e-9)&(abs(w)<.5)];assert len(low)==1
            shift=low[0];axis=P[:,e]/np.sqrt(P[e,e]);vlow=V[:,np.argmin(abs(w-shift))];bound+=4*np.sqrt(max(0.,1-abs(np.vdot(vlow,axis))**2));t=angle/shift;assert t>0;times.append(float(t))
            state=V@(np.exp(-1j*w*t)[:,None]*(V.T@state))
        error=np.linalg.norm(Q.T@state-ideal,ord=2);leak=np.linalg.norm((np.eye(160)-P)@state,ord=2)
        total_error=np.linalg.norm(state-Q@ideal,ord=2);assert total_error<=bound+1e-10
        scans.append(dict(local_strength=strength,total_encoded_state_error=float(total_error),exact_dressing_telescoping_bound=float(bound),gate_error=float(error),massive_leakage=float(leak),total_time=float(sum(times))))
    assert scans[-1]['gate_error']<.02 and scans[-1]['massive_leakage']<.05
    assert scans[-1]['gate_error']<scans[0]['gate_error']
    # Every nonzero penalty term shares a native vertex, hence is physically one-hop.
    for i,j in combinations(range(160),2):
        if L[i,j]!=0:assert set(c['edges'][i])&set(c['edges'][j])
    from w33_pass11389_parabolic_spatial_cover import bloch
    voltage=sp.Matrix(c['x']['line']['fcc_harmonic_displacements']).applyfunc(sp.Rational)
    cover_gap_scans=[]
    for k in (.08,.04,.02,.01):
        gap=float(np.linalg.eigvalsh(bloch(c['edges'],voltage,[k,0,0]))[0])
        assert abs(gap/k**2-27/3200)<1e-5
        cover_gap_scans.append(dict(momentum=k,positive_penalty_gap_over_Delta=gap,ratio_to_k_squared=gap/k**2))
    return dict(modular_escape=modular_entangler_clock(c,P),cover_gap_scans=cover_gap_scans,scalable_penalty_gap='closes as (27/3200)|k|^2 on the connected spatial cover; fixed cell gap is not a volume-independent gap',local_penalty='Delta D^T D on native160 edge sites; finite support on edges sharing a vertex',
        gap_over_Delta='4-sqrt(6)',local_control='u(t)|e><e|',protected_effective_control='u Pi|e><e|Pi',
        gate_pair=[a,b],logical_basis=cm(Q),five_pulse_phase_angles=solution.tolist(),ideal_qubit_gate=cm(ideal),scans=scans,
        continuous_phase_and_Hadamard_qubit_universality=True,
        scope='Actual local finite-penalty Hamiltonians approximate a compiled five-pulse Hadamard on the native two-flag qubit; time is calibrated to the exact dressed low eigenvalue. This realizes the formerly projected BT4048 control in a large-gap limit, an established Zeno/degenerate-perturbation mechanism. Independent address control, preparation, penalty scale, calibration and pulse timing are inputs. The spatial-cover acoustic branch closes the penalty gap at growing volume, so fixed Delta does not provide uniform protection; finite-volume isolation, extra constraints or Delta scaling are required. No scalable entangler, fault tolerance, full four-qutrit compilation or laboratory fidelity claimed.')

def produce():
    c=load();parts={}
    for name,fn in [('fermion_bimodule',lambda:fermion_bimodule(c)),('spectral_selector',lambda:spectral_selector(c)),('collective_geometry',lambda:collective_geometry(c)),('selector_radiative',selector_radiative),('local_gates',lambda:local_gates(c))]:
        parts[name]=fn();print(name,'PASS',flush=True)
    parts['colour_safe_bridge']=colour_safe_bridge(c,parts['fermion_bimodule'],parts['spectral_selector']);print('colour_safe_bridge PASS',flush=True)
    sources=['data/w33_pass11389_parabolic_spatial_cover.json','data/w33_pass11390_11394_context_matter_clock.json','data/w33_pass11395_11402_joint_atlas_spin_portal.json','analysis/w33_pass11389_parabolic_spatial_cover.py','analysis/w33_pass11390_11394_context_matter_clock.py']
    def canonical_source(f):
        p=ROOT/f
        return json.dumps(json.loads(p.read_text()),sort_keys=True,separators=(',',':')).encode()if p.suffix=='.json'else p.read_bytes().replace(b'\r\n',b'\n')
    hashes={f:hashlib.sha256(canonical_source(f)).hexdigest()for f in sources}
    out=dict(source_hash_convention='canonical sorted compact JSON; LF-normalized script bytes',source_sha256=hashes,status='PASS',reservation='06cf9b68a',passes='11403-11407',**parts,scope='Five executed mathematical/conditional-action fronts; no completed TOE or physical Standard Model/gravity prediction.')
    OUT.write_text(json.dumps(out,indent=2)+'\n');return out
if __name__=='__main__':produce()
