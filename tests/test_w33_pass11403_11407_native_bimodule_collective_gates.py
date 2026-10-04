"""Independent native matrix, representation and two-cell Hamiltonian replays."""
import json
import hashlib
import sys
from pathlib import Path
from itertools import combinations,product
import numpy as np
import sympy as sp
from scipy.linalg import block_diag,expm
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
from w33_pass11390_11394_context_matter_clock import load,clock_matrices

def packet():return json.loads((ROOT/'data/w33_pass11403_11407_native_bimodule_collective_gates.json').read_text())
def read(d):return np.array(d['real'])+1j*np.array(d['imag'])
def frames():
    p=packet()['fermion_bimodule']['native_slot_frames'];Fc=[read(F)for F in p['colour_qutrit_seeds']]
    return Fc,np.hstack(Fc),read(p['weak']),read(p['lepton']),read(p['right_up']),read(p['right_down'])
def generators():
    Fc,C,W,L,U,D=frames();T8=np.diag([1,1,-2])/(2*np.sqrt(3));Tc=sum(F@T8@F.conj().T for F in Fc);Tw=W@np.diag([1,-1])@W.conj().T/2
    return Tc,Tw

def test_native_bimodule_rectangles_and_casimirs():
    Fc,C,W,L,U,D=frames();F=np.hstack([C,W,L,U,D]);assert F.shape==(116,16)and np.linalg.norm(F.conj().T@F-np.eye(16))<1e-10
    Tc,Tw=generators();Q=C[:,0:1]@W[:,0:1].conj().T
    # Both gauge factors act nontrivially on the same actual native operator.
    assert np.linalg.norm(Tc@Q-Q@Tc)>.1 and np.linalg.norm(Tw@Q-Q@Tw)>.1
    assert packet()['fermion_bimodule']['left_Weyl_bimodule_dimension']==48
    clocks=packet()['fermion_bimodule']['family_clock_matrices'];X,Z=read(clocks['X']),read(clocks['Z']);eq=np.vstack([np.kron(np.eye(3),X)-np.kron(X.T,np.eye(3)),np.kron(np.eye(3),Z)-np.kron(Z.T,np.eye(3))]);assert np.linalg.matrix_rank(eq,tol=1e-9)==8
    for name,digest in packet()['source_sha256'].items():
        p=ROOT/name;raw=json.dumps(json.loads(p.read_text()),sort_keys=True,separators=(',',':')).encode()if p.suffix=='.json'else p.read_bytes().replace(b'\r\n',b'\n')
        assert hashlib.sha256(raw).hexdigest()==digest

def test_anomalies_on_actual_48_operator_representation():
    Fc,C,W,L,U,D=frames();fields=[(C,W),(U,C),(D,C),(L,W),(D,L),(U,L)];Tc,Tw=generators()
    def representation(G):
        return block_diag(*[np.kron(A.conj().T@G@A,np.eye(B.shape[1]))-np.kron(np.eye(A.shape[1]),(B.conj().T@G@B).T)for A,B in fields])
    gc,gw=representation(Tc),representation(Tw)
    for q,h in [(1/6,.5),(.27,.43)]:
        Y=q*C@C.conj().T-3*q*L@L.conj().T-h*U@U.conj().T+h*D@D.conj().T;gy=representation(Y)
        coeff=[np.trace(gc@gc@gc),np.trace(gc@gc@gy),np.trace(gw@gw@gy),np.trace(gy),np.trace(gy@gy@gy)]
        assert max(abs(z)for z in coeff)<1e-10
    assert packet()['fermion_bimodule']['witten_doublets']==12

def test_higgs_conjugate_character_pairing():
    p=packet()['fermion_bimodule']['single_Higgs_character_pairing'];w=tuple(p['weak']);pairs=[]
    for u,v in product(product(range(3),repeat=2),repeat=2):
        if u!=(0,0)and v!=(0,0)and u!=w and v!=w and tuple((a+b)%3 for a,b in zip(u,v))==tuple(2*a%3 for a in w):pairs.append((u,v))
    assert len(pairs)==6 and (tuple(p['right_up']),tuple(p['right_down']))in pairs

def test_colour_safe_family_split_and_unsafe_old_channels():
    p=packet();Fc,C,W,L,U,D=frames();Tc,_=generators();Ps=[F@F.conj().T for F in Fc]
    assert all(np.linalg.norm(P@Tc-Tc@P)<1e-10 for P in Ps)
    for row in p['colour_safe_bridge']['scans']:
        phi=sum(w*P for w,P in zip(row['family_singular_values'],Ps));actual=np.linalg.eigvalsh(C.conj().T@phi@C)
        assert np.linalg.norm(actual-np.sort(np.repeat(row['family_singular_values'],3)))<1e-12
    norms=[]
    for N in p['spectral_selector']['operators']:
        H=np.zeros((116,116),complex);H[:80,:80]=read(N);norms.append(np.linalg.norm(H@Tc-Tc@H))
    assert max(norms)>.1
    V=read(p['colour_safe_bridge']['other_native_clock_mixing']);assert np.max(abs(abs(V)**2-1/3))<1e-10
    J=np.imag(V[0,0]*V[1,1]*V[0,1].conjugate()*V[1,0].conjugate());assert abs(abs(J)-1/(6*np.sqrt(3)))<1e-10

def test_native_spectral_operator_formula_in_another_basis():
    c=load();p=packet()['spectral_selector'];C=np.zeros((80,80))
    for a,b in json.loads((ROOT/'data/w33_pass11390_11394_context_matter_clock.json').read_text())['portal']['nonedge_orbit']:C[a,b]=C[b,a]=1
    H=c['A']+C;g=c['x']['internal_h27']['group_permutations'];Z=np.eye(80)[np.argsort(g[1])];P=(2*np.eye(80)-Z-Z@Z)/3
    reorder=np.random.default_rng(11404).permutation(80);H=H[np.ix_(reorder,reorder)];P=P[np.ix_(reorder,reorder)]
    xs=np.sort(np.roots([1,-15,54,-27]).real)
    for j,x in enumerate(xs):
        f=np.eye(80)
        for k,y in enumerate(xs):
            if k!=j:f=f@(H@H-y*np.eye(80))/(x-y)
        N=P@H@f/np.sqrt(x);stored=read(p['operators'][j])[np.ix_(reorder,reorder)]
        assert np.linalg.norm(N-stored)<1e-10

def test_collective_mean_and_full_metric_response():
    p=packet()['collective_geometry'];Rx=np.array([[1,0,0],[0,0,-1],[0,1,0]])
    for row in p['scans']:
        n=row['block_size'];M=(1-1/n)*np.eye(3)+Rx/n;u,_,v=np.linalg.svd(M);R=u@v
        assert abs(np.arctan2(R[2,1],R[1,1])-row['collective_angle'])<1e-12
    mats=[sp.Matrix(S).applyfunc(sp.sympify)for S in p['metric_generators']]
    columns=[sp.Matrix([S[0,0],S[1,1],S[2,2],S[0,1],S[0,2],S[1,2]])for S in mats]
    assert sp.Matrix.hstack(*columns).rank()==6 and p['native_metric_strain_rank']==4

def test_conditional_spin_two_null_quotient():
    p=packet()['collective_geometry'];K=sp.Matrix(p['conditional_null_symbol']).applyfunc(sp.sympify);G=sp.Matrix(p['conditional_gauge_map']).applyfunc(sp.sympify)
    assert K*G==sp.zeros(10,4)and 10-K.rank()-G.rank()==2
    # Named transverse plus/cross representatives survive after removing gauge.
    plus=sp.Matrix([0,0,0,0,1,0,0,-1,0,0]);cross=sp.Matrix([0,0,0,0,0,1,0,0,0,0])
    assert K*plus==sp.zeros(10,1)and K*cross==sp.zeros(10,1)
    assert G.row_join(plus).row_join(cross).rank()==6
    assert all(r['physical_polarizations']==2 and r['Ward_residual']<1e-13 for r in p['native_momentum_Ward_samples'])

def test_selector_virtual_and_CW_corrections():
    p=packet()['selector_radiative'];lam=p['lambda_over_J'];mu=p['renormalization_scale_over_J']
    cw=lambda m:-3*m**4*(np.log(m*m/(mu*mu))-1.5)/(16*np.pi**2)
    values=np.array([cw(lam),cw(lam/12),cw(lam/27)]);assert np.linalg.norm(values-p['local_selector_CW_energies'])<1e-15
    assert np.linalg.norm(1+values[1:]-values[0]-p['loop_shifted_gaps'])<1e-15
    Q=np.array([[0,np.sqrt(12),0],[np.sqrt(12),2,6],[0,6,8]])
    for row in p['scans']:
        H=np.diag([0,1,1])-row['eta']*Q;mass=[]
        for d in [np.diag([1,0,0]),np.diag([0,1/12,0]),np.diag([0,0,1/27])]:mass.append((np.linalg.eigvalsh(H+lam*d)[0]-np.linalg.eigvalsh(H-lam*d)[0])/2)
        assert np.linalg.norm(np.array(mass)-row['finite_coupling_odd_masses'])<1e-14
    assert abs(cw(lam)/lam**4-9/(32*np.pi**2))<1e-12

def test_full_native_local_hadamard():
    c=load();p=packet()['local_gates'];L=np.array(c['D'],dtype=float).T@np.array(c['D'],dtype=float);Q=read(p['logical_basis']);angles=p['five_pulse_phase_angles'];a,b=p['gate_pair']
    state=Q.copy();strength=p['scans'][-1]['local_strength']
    for angle,e in zip(angles,[a,b,a,b,a]):
        H=L.copy();H[e,e]+=np.copysign(strength,angle);w,v=np.linalg.eigh(H);shift=w[(abs(w)>1e-9)&(abs(w)<.5)][0];t=angle/shift
        state=v@(np.exp(-1j*w*t)[:,None]*(v.T@state))
    error=np.linalg.norm(state-Q@read(p['ideal_qubit_gate']),ord=2)
    assert abs(error-p['scans'][-1]['total_encoded_state_error'])<1e-10 and error<.009
    assert error<p['scans'][-1]['exact_dressing_telescoping_bound']
    assert all(set(c['edges'][i])&set(c['edges'][j])for i,j in combinations(range(160),2)if L[i,j]!=0)

def test_two_cell_native_interaction_reduction_and_entanglement():
    c=load();p=packet()['local_gates'];a=p['gate_pair'][0];L=np.array(c['D'],dtype=float).T@np.array(c['D'],dtype=float);w,v=np.linalg.eigh(L);energies=np.array([0,4-np.sqrt(6),4,4+np.sqrt(6),8]);F=[];weights=[]
    for energy in energies:
        T=v[:,abs(w-energy)<1e-8];f=T@T[a,:];weights.append(np.vdot(f,f).real);F.append(f/np.linalg.norm(f))
    F=np.column_stack(F);assert np.max(abs(np.array(weights)-np.array([81,24,30,24,1])/160))<1e-12
    coupling=p['modular_escape']['scans'][-1]['coupling'];delta=np.kron(np.sqrt(weights),np.sqrt(weights));H=np.diag(np.add.outer(energies,energies).ravel())+coupling*np.outer(delta,delta)
    # Check the exact reduction on a25600-component native two-particle wavefunction.
    C=np.random.default_rng(11407).normal(size=(5,5));psi=F@C@F.T;ee=np.zeros((160,160));ee[a,a]=1
    native=L@psi+psi@L.T+coupling*ee@psi@ee;reduced=F@(H@C.ravel()).reshape(5,5)@F.T
    assert np.linalg.norm(native-reduced)<1e-10
    w,v=np.linalg.eigh(H);state=v@(np.exp(-1j*w*np.pi/w[0])*v[0,:]);amp=state[0]
    assert np.linalg.norm(state+np.eye(25)[:,0])<.004
    concurrence=2*abs(amp-1)/(abs(amp)**2+3);assert concurrence>.999

def test_native_clock_resonance_and_generic_phase():
    c=load();D=np.array(c['D'],dtype=float);A=D[:40].T@D[:40]/4;B=D[40:].T@D[40:]/4
    grover=(np.eye(160)-2*B)@(np.eye(160)-2*A);U,_,_=clock_matrices(c,[0,0,0]);even=list(range(0,320,2))
    assert np.linalg.norm(grover-(U@U)[np.ix_(even,even)])<1e-12
    for phi,expected in [(.1,81),(np.pi,82)]:
        W=(np.eye(160)+(np.exp(-1j*phi)-1)*B)@(np.eye(160)+(np.exp(-1j*phi)-1)*A)
        assert np.linalg.norm(W.conj().T@W-np.eye(160))<1e-12
        assert sum(abs(np.angle(z))<1e-9 for z in np.linalg.eigvals(W))==expected

def test_connected_cover_gap_collapse():
    from w33_pass11389_parabolic_spatial_cover import bloch
    c=load();p=packet()['local_gates'];E=sp.Matrix(c['x']['line']['fcc_harmonic_displacements']).applyfunc(sp.sympify)
    for row in p['cover_gap_scans']:
        k=row['momentum'];gap=np.linalg.eigvalsh(bloch(c['edges'],E,[k,0,0]))[0]
        assert abs(gap-row['positive_penalty_gap_over_Delta'])<1e-12
        assert abs(gap/k**2-27/3200)<1e-5
    assert p['cover_gap_scans'][-1]['positive_penalty_gap_over_Delta']<p['cover_gap_scans'][0]['positive_penalty_gap_over_Delta']/60
