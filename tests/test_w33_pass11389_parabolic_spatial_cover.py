"""Independent checks of the stored period map and its propagation operator."""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
SPEC = importlib.util.spec_from_file_location("spatial_cover",ROOT/"analysis/w33_pass11389_parabolic_spatial_cover.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def payload():
    return json.loads((ROOT/"data/w33_pass11389_parabolic_spatial_cover.json").read_text())


def matrix(x):
    return sp.Matrix(x).applyfunc(sp.Rational)


def test_period_image_and_harmonic_representative():
    p = payload()["line"]
    V,H = matrix(p["integer_voltage"]),matrix(p["harmonic_voltage"])
    D = sp.zeros(80,160)
    for i,(a,b) in enumerate(p["edges"]):
        D[a,i],D[b,i] = -1,1
    assert D*H == sp.zeros(80,3)
    # Difference must be exact, rather than a different cohomology class.
    L = D*D.T
    phi = (L+sp.ones(80)/80).inv()*D*(H-V)
    assert D.T*phi == H-V
    stored = matrix(p["cell_positions_in_primitive_periods"])
    assert D.T*stored == H-V
    classes = {tuple(sp.frac(x) for x in stored.row(i)) for i in range(80)}
    assert len(classes)==14
    assert all(x.q == 1 for x in V)
    assert hermite_normal_form(V.T.applyfunc(sp.Integer)) == sp.eye(3)


def test_finite_covers_connected_and_four_regular():
    p = payload()["line"]
    V = np.array(matrix(p["integer_voltage"]),dtype=int)
    for n in (2,3):
        def label(v,t):
            x,y,z = (int(a)%n for a in t)
            return v+80*(x+n*y+n*n*z)
        adjacency = [[] for _ in range(80*n**3)]
        for x in range(n):
            for y in range(n):
                for z in range(n):
                    t = np.array([x,y,z])
                    for (a,b),voltage in zip(p["edges"],V):
                        aa,bb = label(a,t),label(b,t+voltage)
                        adjacency[aa].append(bb)
                        adjacency[bb].append(aa)
        assert all(len(row)==4 for row in adjacency)
        seen,queue = {0},[0]
        for a in queue:
            for b in adjacency[a]:
                if b not in seen:
                    seen.add(b)
                    queue.append(b)
        assert len(seen)==80*n**3


def test_primitive_fcc_metric_and_actual_levi_group():
    p = payload()["line"]
    K = matrix(p["diffusion_tensor"])
    U,F = matrix(p["primitive_fcc_change"]),matrix(p["fcc_roots"])
    assert abs(U.det())==1
    assert U.T*K.inv()*U == sp.Rational(3200,27)*F.T*F
    reps = [matrix(r) for r in p["levi_cohomology_matrices"]]
    keys = {tuple(r) for r in reps}
    assert len(keys)==24
    assert all(tuple(a*b) in keys for a in reps for b in reps)
    assert all(r.T*K*r==K for r in reps)
    assert sum(sp.trace(r)**2 for r in reps)==24
    assert payload()["point"]["center_order"]==3
    assert payload()["line"]["center_order"]==27


def test_bloch_gauge_equivalence_and_quartic_dispersion():
    p = payload()["line"]
    V,H,E = (np.array(matrix(p[k]),dtype=float) for k in
             ("integer_voltage","harmonic_voltage","fcc_harmonic_displacements"))
    for k in (np.array([.12,-.07,.19]),np.array([.31,.11,-.17])):
        integer_spectrum = np.linalg.eigvalsh(MODULE.bloch(p["edges"],V,k))
        harmonic_spectrum = np.linalg.eigvalsh(MODULE.bloch(p["edges"],H,k))
        assert np.max(abs(integer_spectrum-harmonic_spectrum))<2e-13
        assert integer_spectrum[0]>=-1e-13
    a,b = float(sp.Rational(p["cartesian_quartic_axis"])),float(sp.Rational(p["cartesian_quartic_mixed"]))
    for k in (np.array([.012,-.009,.007]),np.array([.015,.004,-.006])):
        actual = np.linalg.eigvalsh(MODULE.bloch(p["edges"],E,k))[0]
        predicted = 27/3200*(k@k)+a*sum(k**4)+b*sum(k[i]**2*k[j]**2 for i in range(3) for j in range(i+1,3))
        assert abs(actual-predicted)<2e-12


def test_only_four_metric_tangents_at_the_symmetric_cell():
    E = matrix(payload()["line"]["fcc_harmonic_displacements"])
    response = sp.Matrix([[x*x,y*y,z*z,x*y,x*z,y*z] for x,y,z in E.tolist()])
    assert response.rank()==4
    for tensor in (sp.Matrix([-1,1,0,0,0,0]),sp.Matrix([-1,0,1,0,0,0])):
        assert response*tensor==sp.zeros(160,1)


def internal_witness():
    p=payload()
    info=p["internal_h27"]
    PN=np.zeros((80,80))
    for orbit in info["vertex_orbits"]:
        PN[np.ix_(orbit,orbit)]=1/len(orbit)
    Q=np.eye(80)-PN
    S=np.array(info["qutrit_seed_real"])+1j*np.array(info["qutrit_seed_imag"])
    X=np.zeros((80,80),dtype=complex)
    X[np.array(info["group_permutations"][9]),np.arange(80)]=1
    T=np.hstack([S,X@S,X@X@S])
    return p,PN,Q,T


def test_exact_forest_witness_and_actual_qutrit_factor():
    p,PN,Q,T=internal_witness()
    point,info=p["point"],p["internal_h27"]
    H=np.array(matrix(point["harmonic_voltage"]),dtype=float)
    psi=np.array(matrix(info["forest_phase_potential"]),dtype=float)
    L0=MODULE.bloch(point["edges"],H,[0,0,0])
    assert np.linalg.norm(T.conj().T@T-np.eye(21))<1e-12
    assert abs(np.trace(PN)-14)<1e-12
    for k in ([.17,-.09,.31],[1.1,-2.3,.87],[-2.2,.21,1.9]):
        L=MODULE.bloch(point["edges"],H,k)
        U=np.diag(np.exp(-1j*psi@k))
        assert np.linalg.norm((L-U@L0@U.conj().T)@Q)<2e-12
        A=T.conj().T@L@T
        assert np.linalg.norm(A-np.kron(np.eye(3),A[:7,:7]))<2e-12
        e=np.linalg.eigvalsh(A[:7,:7])
        assert np.max(abs(e-np.array([4-np.sqrt(6)]*2+[4]*3+[4+np.sqrt(6)]*2)))<2e-12


def test_additional_hopping_carries_the_spectator_qutrit():
    p,PN,Q,T=internal_witness()
    point,info=p["point"],p["internal_h27"]
    E=np.array(matrix(point["fcc_harmonic_displacements"]),dtype=float)
    roots=np.array(info["completion"]["fcc_positive_root_directions"])
    for k in (np.array([.011,-.009,.007]),np.array([.3,-.2,.17])):
        L=MODULE.bloch(point["edges"],E,k)
        ell=sum(2-2*np.cos(roots@k))
        shift=27/12800*ell
        completed=L+shift*Q
        assert np.linalg.eigvalsh(completed)[0]>=-1e-12
        assert np.linalg.norm((completed-L)@PN)<1e-12
        assert np.linalg.norm(T.conj().T@(completed-L)@T-shift*np.eye(21))<1e-12
        reduced=T.conj().T@completed@T
        gate=np.kron(np.diag([1,np.exp(1j*np.sqrt(2)),np.exp(1j*np.sqrt(3))]),np.eye(7))
        assert np.linalg.norm(gate@reduced-reduced@gate)<2e-12
        full_gate=np.eye(80)+T@(gate-np.eye(21))@T.conj().T
        assert np.linalg.norm(full_gate.conj().T@full_gate-np.eye(80))<2e-12
        assert np.linalg.norm(full_gate@completed-completed@full_gate)<3e-12
    tiny=np.array([.001,.002,-.001])
    ell=sum(2-2*np.cos(roots@tiny))
    assert abs((27/12800*ell)/(tiny@tiny)-27/3200)<1e-8


def test_small_prime_period_family_independent_of_the_q3_builder():
    witnesses=[MODULE.prime_period_witness(q) for q in (2,3,5)]
    assert [w["harmonic_rank"] for w in witnesses]==[2,3,5]
    assert [w["root_lattice"] for w in witnesses]==["A2","A3","A5"]
    assert all(w["gram_law"] and w["period_law"] for w in witnesses)
    assert witnesses[1]["raw_period_lattice"]==payload()["line"]["raw_period_lattice"]


def test_positive_gapless_action_and_index_preserving_rank_lift():
    p,PN,Q,T=internal_witness()
    point,info=p["point"],p["internal_h27"]
    E=np.array(matrix(point["fcc_harmonic_displacements"]),dtype=float)
    roots=np.array(info["completion"]["fcc_positive_root_directions"])
    Gamma=np.diag([1]*40+[-1]*40)
    for k in (np.array([0.,0.,0.]),np.array([.03,-.02,.01])):
        L=MODULE.bloch(point["edges"],E,k)
        A=4*np.eye(80)-L
        assert np.linalg.norm(Gamma@A+A@Gamma)<1e-12
        zero=Q-A@A@Q/6
        assert np.linalg.norm(zero@zero-zero)<2e-12
        assert abs(np.trace(zero)-26)<1e-12
        center_zero=T.conj().T@zero@T
        assert abs(np.trace(center_zero)-9)<1e-12
        assert abs(np.trace(T.conj().T@Gamma@zero@T)+3)<1e-12
        shift=27/12800*sum(2-2*np.cos(roots@k))
        gapless=PN@L@PN+Q@A@A@Q+shift*Q
        assert np.linalg.eigvalsh(gapless)[0]>=-2e-12
        assert np.linalg.norm(T.conj().T@gapless@T-(T.conj().T@A@T)@(T.conj().T@A@T)-shift*np.eye(21))<3e-12
    S=T[:,:7]
    ew,ev=np.linalg.eigh(S.conj().T@Gamma@S)
    plus,minus=ev[:,ew>.5],ev[:,ew<-.5]
    A0=4*np.eye(80)-MODULE.bloch(point["edges"],E,[0,0,0])
    B=plus.conj().T@(S.conj().T@A0@S)@minus
    u,d,v=np.linalg.svd(B,full_matrices=True)
    lifted=B+.1*np.outer(u[:,-1],v[2,:])
    chiral=np.block([[np.zeros((3,3)),lifted],[lifted.conj().T,np.zeros((4,4))]])
    square=np.linalg.eigvalsh(chiral@chiral)
    assert np.max(abs(square-np.array([0,.01,.01,6,6,6,6])))<2e-12
