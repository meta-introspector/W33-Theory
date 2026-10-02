"""Independent tensor, symmetry, displacement and interval regression controls."""
from pathlib import Path
import json,sys
from fractions import Fraction
import numpy as np
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11270_full_condensate_moduli as C
import w33_pass11273_wilson_uniform_tail as W

def certificate(name):return json.loads((ROOT/'data'/name).read_text())

def test_full_moduli_compact_flavor_invariance_and_radial_kahler():
    p,N=C.reference();assert N.shape==(81,11)
    a=np.array([[.2,.1+.2j,0],[.1-.2j,-.1,.2j],[0,-.2j,-.1]])
    U=expm(1j*a);rot=(p.reshape(27,3)@U.T).ravel()
    for k in (0,.01,-.01):assert abs(C.potential(rot,k)-C.potential(p,k))<1e-10
    data=certificate('w33_pass11270_full_condensate_moduli.json')
    for row in data['rows']:
        ev=np.array(row['extrapolated_eigenvalues']);assert sum(abs(ev)<.002)==8 and sum(ev>.1)==14

def test_full_hessian_predicts_independent_displacements():
    data=certificate('w33_pass11270_full_condensate_moduli.json');p,N=C.reference();dirs=np.column_stack([N,1j*N])/np.sqrt(2)
    rng=np.random.default_rng(11270)
    for row in data['rows']:
        q=row['stationary_radius']*p;k=row['kappa'];steps=row['steps'];H=(4*np.array(steps[1]['Hessian'])-np.array(steps[0]['Hessian']))/3
        for _ in range(3):
            v=rng.normal(size=22);v/=np.linalg.norm(v);d=dirs@v;h=.0002
            direct=(C.potential(q+h*d,k)+C.potential(q-h*d,k)-2*C.potential(q,k))/h**2
            assert abs(direct-v@H@v)<.002

def test_canonical_SM_mass_and_second_basis_covariance():
    d=certificate('w33_pass11271_canonical_sm_higgs_bridge.json');_,tensor,_=C.D.I.tensors();M=tensor[:,:,0].astype(complex)
    Y=np.diag([float(Fraction(x)) for x in d['hypercharge_diagonal']]);B=np.load(ROOT/'artifacts/e6_27rep_basis_export/E6_basis_78.npy')
    assert np.linalg.matrix_rank(M)==10 and np.max(abs(Y.T@M+M@Y))==0
    heavy=d['exotic_mass_support'];assert heavy==list(range(17,27))
    assert np.all(M[:17,:]==0) and np.linalg.matrix_rank(np.kron(M,np.diag([1,2,3])))==30
    rng=np.random.default_rng(11271);U=np.linalg.qr(rng.normal(size=(27,27))+1j*rng.normal(size=(27,27)))[0]
    Mp=U.conj()@M@U.conj().T
    for j in d['unbroken_root_basis_indices']:
        X=U@B[j]@U.conj().T
        assert np.max(abs(X.T@Mp+Mp@X))<1e-12

def test_symmetric_family_operator_and_spectator_positive_mass():
    rng=np.random.default_rng(112711);_,d,_=C.D.I.tensors();psi=rng.normal(size=(27,3))+1j*rng.normal(size=(27,3));h=rng.normal(size=(27,3,3));h=h+h.transpose(0,2,1)
    U=np.diag(np.exp([.2,-.05,-.15]));inv=np.linalg.inv(U)
    val=lambda p,H:np.einsum('abc,ai,bj,cij->',d,p,p,H)
    transformed=np.einsum('ij,ajk,kl->ail',inv.T,h,inv)
    assert abs(val(psi@U.T,transformed)-val(psi,h))<1e-8
    c=certificate('w33_pass11271_chiral_symmetric_yukawa_completion.json');mass=np.array(c['spectator_mass_matrix'],float)
    assert min(np.linalg.eigvalsh(mass))>0

def test_soft_IR_is_recorded_without_claiming_poles():
    c=certificate('w33_pass11272_goldstone_ir_and_pole_audit.json');a=np.array(c['projected_linear_soft_slopes'])
    assert len(a)==79 and min(a)<0<max(a)
    assert abs(sum(a*a)/(32*np.pi**2)-c['log_abs_angular_displacement_curvature_coefficient'])<1e-12
    assert c['physical_poles_computed'] is False
    counter=c['same_static_curvature_counterexample'];M=counter['M_squared']
    poles=[M/counter['Z1'],M/counter['Z2']]
    assert poles[0]!=poles[1] and all(z*p-M==0 for z,p in zip([counter['Z1'],counter['Z2']],poles))

def test_computed_hard_tadpole_resums_global_goldstones():
    c=certificate('w33_pass11272_resummed_condensate_radial_control.json');m=c['parameters']['m'];r=c['shifted_radius']
    dV=m*m*(-14*81/(49*r**15)+6*27/(7*r**7))
    assert abs(dV/(2*r)-c['tree_Goldstone_mass_squared_at_shifted_radius'])<1e-13
    assert abs(c['tree_Goldstone_mass_squared_at_shifted_radius']+c['computed_hard_zero_momentum_shift'])<1e-12
    assert c['computed_hard_zero_momentum_shift']<0 and c['shifted_radius']>1
    assert c['nonzero_physical_poles_computed'] is False

def test_uniform_Wilson_bound_including_Brillouin_corners():
    for n in (8,32,128):
        h=2*np.pi/n
        for p in [np.array([0,0,0,0]),np.array([-n//2]*4),np.array([1,2,-1,0])]:
            x=h*p;E=(np.sum(np.sin(x)**2)+np.sum(1-np.cos(x))**2)/h**2
            assert .4*(p@p)-1e-12<=E<=11*(p@p)+1e-12
    t32,t64=W.uniform_tail(32),W.uniform_tail(64)
    assert Fraction(t64['upper_exact'])<Fraction(1,10**24)<Fraction(t32['upper_exact'])
    assert certificate('w33_pass11273_wilson_uniform_tail.json')['uniform_tail']==t64

def test_mass_supertrace_does_not_cancel_vacuum_running():
    c=certificate('w33_pass11274_scale_and_sequestering.json');b=[Fraction(162,49)]*10+[Fraction(486,49)]*2+[Fraction(648,7),Fraction(486,7)];f=[Fraction(81,49)]*8+[Fraction(324,49)]*2+[Fraction(81)]
    assert sum(b)-2*sum(f)==0
    assert sum(x*x for x in b)-2*sum(x*x for x in f)==Fraction(c['canonical_supertrace_m_fourth'])>0
    rho=[Fraction(2),Fraction(-5),Fraction(8)];weights=[Fraction(1,7),Fraction(2,7),Fraction(4,7)];shift=Fraction(13,3)
    avg=sum(w*x for w,x in zip(weights,rho));avg2=sum(w*(x+shift) for w,x in zip(weights,rho))
    assert [x-avg for x in rho]==[x+shift-avg2 for x in rho]
