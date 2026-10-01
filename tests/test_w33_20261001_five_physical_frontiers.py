"""Exact block replay, gauge normalization, independent lattice control and scale minors."""
from pathlib import Path
import sys,json
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_20261001_exact_cartan_mirror_completion as C
import w33_20261001_complete_gauge_scalar_fermion_loop as L
import w33_20261001_wilson_gravity_refinement as W
import w33_20261001_anomalous_scale_completion as D
import w33_20261001_global_e6_cartan_covariants as I

def stored(name):return json.loads((ROOT/'data'/('w33_20261001_'+name+'.json')).read_text())

def test_exact_plane_universal_blocks_and_three_null_completion():
    p=C.payload();assert C.compare_certificate(p,stored('exact_cartan_mirror_completion'))
    assert p['finite_completion_control']['rank_before_after']==[78,81]
    # The exact block proof must reject a sign corruption, even if dimensions agree.
    mats,_=C.operators();bad=[m.copy() for m in mats]
    i,j=np.argwhere(bad[0]!=0)[0];bad[0][i,j]*=-1;bad[0][j,i]*=-1
    rejected=False
    try:C.exact_blocks(bad)
    except AssertionError:rejected=True
    assert rejected

def test_exact_full_loop_gauge_identity_and_radial_scale_checks():
    p=L.payload();assert L.compare_certificate(p,stored('complete_gauge_scalar_fermion_loop'))
    assert p['gauge_mass_identity']['exact_quadratic_coefficient_matrices_checked']==9
    witnesses=p['full_one_loop']['interval_sign_proof']['sign_witnesses']
    assert witnesses[0]['Rayleigh_bounds'][1]<0<witnesses[1]['Rayleigh_bounds'][0]
    assert p['scale_dynamics']['heat_volume_over_Planck_fourth']=='144*pi**2/(kappa*theta)'

def test_independent_curved_lattice_matrix_and_exact_mass_grading():
    p=stored('wilson_gravity_refinement')
    assert W.independent_control()['error']<2e-5
    assert p['Richardson_absolute_error']<.05
    assert abs(p['doubling_control']['Richardson_ratio_to_single_continuum_species']-16)<.2
    assert p['doubling_control']['Wilson0_refinement_rows'][0]['ratio_to_single_continuum_species']<0
    assert W.finite_fiber(.2)['finite_Hilbert_dimension']==162
    # Independent brute-force momentum list versus compressed shell enumeration.
    N=8;h=2*np.pi/N;mom=2*np.pi*np.fft.fftfreq(N);si=np.sin(mom)/h;wi=(1-np.cos(mom))/h
    triples=np.array(list(__import__('itertools').product(range(N),repeat=3)))
    q=np.sum(si[triples]**2,axis=1);w=np.sum(wi[triples],axis=1)
    brute=W.block_perturbation(N,.2,q,w,np.ones(N**3))
    assert abs(brute-W.response(N,.2))<1e-9
    # Product-square identity on a concrete external matrix and graded finite block.
    A=np.array([[0,1j],[2,0]],complex);DF=np.block([[np.zeros_like(A),A],[A.conj().T,np.zeros_like(A)]])
    gamma=np.diag([1,1,-1,-1]);ext=np.array([[.3,1],[1,-.7]])
    total=np.kron(ext,gamma)+np.kron(np.eye(2),DF)
    assert np.max(abs(total@total-np.kron(ext@ext,np.eye(4))-np.kron(np.eye(2),DF@DF)))<1e-12

def test_exact_anomalous_scale_stability_and_nonzero_vacuum_boundary():
    p=D.payload();assert p==stored('anomalous_scale_completion')
    assert p['exact_determinant']=='32*B*lambda*mu**4*v**2'
    assert min(p['dimensionless_control_eigenvalues'])>0
    assert p['vacuum_energy']=='-B mu^4/2 before an independent constant'

def test_global_cubic_invariants_and_analytic_covectors():
    p=I.payload();assert I.compare_certificate(p,stored('global_e6_cartan_covariants'))
    assert p['exact_generator_invariance']['E6_lower_and_dual_generator_checks']==156
    assert p['numerical_controls']['T_normal_slice_gradient_error']<1e-10
    # Use the actual global differential circuit in the completed mass operator.
    E=C.plane()/np.sqrt(3);z=np.exp(2j*np.pi/9);q=np.array([1,z,z.conjugate()])/np.sqrt(3);phi=E@q
    _,_,g6,g12=I.evaluate(phi);_,op=C.operators();M=op(phi)
    lift=.001*np.outer(phi.conj(),phi.conj())+.002/16*np.outer(g6,g6)+.003/100*np.outer(g12,g12)
    assert np.max(abs(M.conj().T@lift))<1e-11
    vals=np.linalg.svd(M+lift,compute_uv=False)
    assert np.max(abs(np.sort(vals[78:])-[.001,.002,.003]))<1e-11
