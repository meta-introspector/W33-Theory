import sys,json,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11448_11455_subgroup_toron_auxiliary_control as N


def packet():return json.loads(N.OUT.read_text())


def test_source_binding():
 for n,h in packet()['source_sha256'].items():assert hashlib.sha256(json.dumps(json.loads((ROOT/n).read_text()),sort_keys=True,separators=(',',':')).encode()).hexdigest()==h


def test_su3_and_casimir_basis_invariance():
 K,B=N.P.gauge_fixed_space();p=packet()['subgroup'];R=np.linalg.qr(np.random.default_rng(1).normal(size=(8,8)))[0];Q=np.einsum('ab,bij->aij',R,K)
 assert np.linalg.norm(sum(k@k for k in K)-sum(k@k for k in Q))<1e-11
 assert p['E6_centralizer_dimension']==16
 assert p['center_dimension']==0 and p['family_component_norm']<1e-10 and p['fixed_dimension']==27
 assert sum(n for _,n in p['representation_Casimir_blocks'])==81


def test_mixed_soft_samples_are_scoped():
 p=packet()['soft'];assert len(p['rows'])==24
 assert all(np.isfinite(r['over_h4']) for r in p['rows'])
 assert 'do not certify' in p['scope']


def test_actual_toron_endpoint_covariance():
 E,_,G=N.flat_weyl(2,0);F,_,_=N.flat_weyl(2,2*np.pi)
 assert np.linalg.norm(F@F.conj().T-G[:,None]*(E@E.conj().T)*G.conj()[None,:])<1e-11
 for r in packet()['torons']['rows']:
  assert r['Wilson_gap']>.9 and r['endpoint_projector_error']<1e-10 and r['minimum_overlap_singular_value']>.1


def test_stationary_matter_force_translation_identity():
 X=np.array(N.previous()['lorentzian']['coordinates']);f=np.arange(5)/5;w=np.ones(5)/5
 K=N.P.scalar_refinement(X,w,.2)[0];S=N.P.scalar_refinement(X+np.array([.13,-.17,.31,.25]),w,.2)[0]
 assert np.linalg.norm(K-S)<1e-11
 assert max(r['step_error'] for r in packet()['stress']['rows'])<1e-7


def test_auxiliary_elimination_does_not_remove_determinant():
 p=packet()['auxiliary'];assert abs(p['factorization_error'])<1e-10 and p['derivative_error']<1e-6
 assert abs(p['heavy_logdet_phi_derivative'])>1
 # A normalized Gaussian is a different integral, already in scalar dimension.
 K=3.;J=.2;L=2.;full=np.array([[K,J],[J,L]])
 assert abs(np.linalg.det(full)-K*(L-J*J/K))<1e-12
 assert abs(np.linalg.det(full)-(L-J*J/K))>1


def test_joint_slice_is_not_physical_prediction():
 p=packet()['joint_slice'];assert .9<p['Majorana_scale']<1.1 and p['radial_curvature']>0 and abs(p['radial_force'])<1e-3
 assert len(p['pin_angles'])==2 and 'fixed inputs' in p['scope']


def test_reset_environment_minimum_from_actual_kraus_rank():
 # Keep logical span{e0,e1}; reset leakage e2..e11 to e0.
 P=np.diag([1.,1.]+[0.]*10);ks=[P]
 for j in range(2,12):
  K=np.zeros((12,12));K[0,j]=1;ks.append(K)
 assert np.linalg.norm(sum(K.T@K for K in ks)-np.eye(12))<1e-12
 assert np.linalg.matrix_rank(np.stack([K.ravel() for K in ks]))==11
 assert packet()['instrument']['reset_required_environment_dimension']==11
 for r in packet()['instrument']['rows']:assert r['TP_error']<1e-11 and r['logical_accept_error']<1e-11


def test_relay_faults_have_real_counterexamples():
 p=packet()['relay'];assert p['fault_count']==7*15 and p['relay_fault_count']>0
 assert p['three_site_fault_count']>0
 # A fault on relay after routing cannot be assumed absent on the next operation.
 assert any(r['output_support']==[1] for r in p['rows'])


def test_complete_two_cell_reset_target():
 p=N.reset_target();n=p['native_cell_dimension'];U=np.eye(n*n)[:,p['permutation']];ks=[U.reshape(n,n,n,n)[:,e,:,0] for e in range(n)]
 psi=np.zeros(n,complex);psi[:2]=[1/np.sqrt(2),1j/np.sqrt(2)];rho=np.outer(psi,psi.conj())
 assert np.linalg.norm(sum(K@rho@K.T for K in ks)-rho)<1e-12
 for j in range(2,n):
  v=np.eye(n)[:,j];out=sum(np.outer(K@v,K@v) for K in ks)
  assert np.linalg.norm(out-np.diag([1.]+[0.]*(n-1)))<1e-12
 assert p['Kraus_rank']==11 and p['unitarity_error']==0
