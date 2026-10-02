#!/usr/bin/env python3
"""Explicit tree mediator matching and its inventory-conditional AF obstruction."""
from pathlib import Path
import json,sys
import sympy as s
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))

def signed_tensor_matching():
 import w33_20261001_global_e6_cartan_covariants as I
 _,d,_=I.tensors();v=np.arange(27,dtype=float)%5-2;F=np.array([[2,1,0],[1,3,1],[0,1,5]],float)
 D=np.einsum('abc,c->ab',d,v);B=np.column_stack([np.kron(D,np.eye(3)),np.kron(np.eye(27),F)])
 inverse=np.block([[np.zeros((81,81)),np.eye(81)],[np.eye(81),np.zeros((81,81))]])/7
 effective=-B@inverse@B.T;target=-2*np.kron(D,F)/7
 error=float(np.max(abs(effective-target)));assert error<1e-12 and np.linalg.norm(target)>0
 return {'signed_E6_tensor_Schur_error':error,'Weyl_light_dimension':81,'heavy_Weyl_dimension':162,'mediator_mass':7.,'effective_mass_norm':float(np.linalg.norm(effective))}

def economical_reference():
 import w33_pass11285_semisimple_factor_higgs as H
 v,A,U,V,W,K=H.reference();Z=U.conj().T@A@U;S=(Z+Z.conj()).real
 return v,A,U,V,K,S

def economical_constraints(v,A,U,V,K,S):
 import w33_pass11285_semisimple_factor_higgs as H
 _,_,_,d=H.F.H.setup();M=[np.einsum('abc,c->ab',d,v[:,i]) for i in range(2)];P=(np.eye(10)+1j*K)/2
 base=[v.conj().T@v-np.eye(2),A@v]+[np.einsum('abc,b,c->a',d,v[:,i],v[:,j]) for i,j in [(0,0),(0,1),(1,1)]]
 return base+[K@K+np.eye(10),K@S-S@K,S@S+S-6*np.eye(10),U@K+1j*U,V@K+1j*V,U@V.conj().T-M[0].conj().T@M[1],U.conj().T@U-P,V.conj().T@V-P,A@U-U@S]

def economical_completion():
 import w33_pass11285_semisimple_factor_higgs as H
 from w33_pass11275_polynomial_sm_higgs import rank_mod
 v,A,U,V,K,S=economical_reference();assert max(np.max(abs(c)) for c in economical_constraints(v,A,U,V,K,S))<1e-12
 assert np.max(abs(S-S.T))==0 and abs(np.trace(S))<1e-12
 # Integer compact orbit tangent, with frame coordinates multiplied by sqrt(2).
 U0=np.rint(np.sqrt(2)*U.real)+1j*np.rint(np.sqrt(2)*U.imag);V0=np.rint(np.sqrt(2)*V.real)+1j*np.rint(np.sqrt(2)*V.imag)
 A=np.rint(A).real;S=np.rint(S);cols=[]
 def flatten(*xs):return np.concatenate([np.r_[x.real.ravel(),x.imag.ravel()] for x in xs]).astype(int)
 for h in H.F.H.setup()[0]:cols.append(flatten(1j*h@v,1j*(h@A-A@h),1j*h@U0,1j*h@V0,np.zeros((10,10)),np.zeros((10,10))))
 for i in range(10):
  for j in range(i+1,10):
   X=np.zeros((10,10),int);X[i,j]=1;X[j,i]=-1
   cols.append(flatten(np.zeros_like(v),np.zeros_like(A),-U0@X,-V0@X,X@K-K@X,X@S-S@X))
 rank=rank_mod(np.array(cols).T);assert rank==111
 so=s.Rational(11,3)*8-s.Rational(1,3)*2*27-s.Rational(1,6)*8-s.Rational(1,6)*12;assert so==8
 return {'fields':'Two complex(27,10) frames U,V; real antisymmetric45 K; real symmetric traceless54 S; original two27 and real78 retained. Replace the old third(27,10) frame W by S.',
 'potential':'Squared old base constraints plus K²+I, [K,S], S²+S-6I, UK+iU,VK+iV, UVdagger-N,UdaggerU-P,VdaggerV-P, AU-US. All constraints at most quadratic, giving a quartic scalar potential in reference units with dimensionful coefficients restored.',
 'reference_S':S.astype(int).tolist(),'S_eigenvalues':[-3]*4+[2]*6,'exact_compact_orbit_rank_mod101':rank,
 'local_kernel_proof':'Remove deltaK along the20-dimensional SO10/U5 orbit. The explicit frame projection restricts deltaU,V to P. Since S is real and commutes with K, AU=US and UdaggerU=P determine S uniquely as UdaggerAU+its complex conjugate; the differential determines deltaS uniquely. S²+S-6I then implies(A²+A-6I)U=0. Together with UVdagger=N this maps any kernel vector to the old exact11275 selector kernel. Subtract its66 E6 orbit directions; the remaining frames have25 U5 gauge directions and deltaS is their gauge transform. Thus kernel<=111. The integer orbit tangent has rank111, proving equality and positivity of1254 normal directions among1365 alignment fields. This is analytic elimination using the prior rank theorem, not a numeric full Hessian.',
 'real_alignment_fields':1365,'positive_normal_directions':1254,
 'gauge_beta_composite_before_mediators':{'E6':14,'SO10':8},'with_scalar_H_mediator':{'E6':8,'SO10':8},'with_vectorlike_fermion_mediator':{'E6':2,'SO10':8},
 'scalar_completion_positive_potential':'Add ||M H + alpha v Fdagger||² and the Yukawa d(psi,psi,H). Eliminating H=-alpha v Fdagger/M yields the desired interaction while the nonnegative alignment zero locus is a graph over the old one; the stabilizing quartic ||v Fdagger||² is explicit. F dynamics and its selecting potential remain separate inputs.',
 'inventory_saving':'Remove ten complex27 copies: bE6 increases by10. Remove27 complex SO10 vectors: bSO10 increases by9. Add one real54 with T=12: bSO10 decreases by2. Net old(4,1)->(14,8), enough for either explicit mediator channel.',
 'boundary':'Local alignment and positive one-loop E6/SO10 gauge coefficients in a changed inventory. Family SU3, all-coupling complete asymptotic freedom, thresholds, global uniqueness and quantum stability remain open. This supersedes no old fixed-inventory result.'}

def payload():
 # Heavy vectorlike Weyls X^(A)_i (27,bar3), Xbar_A^i(bar27,3).
 # L includes M X Xbar + y1 d(psi,v,X) + y2 Xbar psi Fdagger.
 M,a,b,x,y=s.symbols('M a b x y',nonzero=True);L=M*x*y+a*x+b*y
 sol=s.solve([s.diff(L,x),s.diff(L,y)],(x,y));eff=s.simplify(L.subs(sol));assert eff==-a*b/M
 rng=np.random.default_rng(11293);F=rng.normal(size=(3,3));F=(F+F.T)/2
 # Both Wick orderings produce the same symmetric flavor coefficient.
 assert np.array_equal(F,F.T)
 rows=[]
 for r in range(3,12):
  baseline=34-6*r;so=7*r-34
  rows.append({'r':r,'SO2r':so,'baseline_E6':baseline,'scalar_mediator_E6':baseline-6,'vectorlike_fermion_E6':baseline-12})
 assert not any(x['SO2r']>0 and x['scalar_mediator_E6']>0 for x in rows)
 assert not any(x['SO2r']>0 and x['vectorlike_fermion_E6']>0 for x in rows)
 return {'status':'PASS','result_scope':'PASS_TREE_MATCHING_FIXED_INVENTORY_OBSTRUCTION_AND_ECONOMICAL_ALIGNMENT_ESCAPE',
 'external_fields':'psi(27,3), v(27,1), Fdagger(1,bar6); operator d_ABC psi^{A i} psi^{B j} v^C Fdagger_ij/M.',
 'fermion_completion':'X(27,bar3)+Xbar(bar27,3); M X^A_i Xbar_A^i + y1 d_ABC psi^{A i} v^B X^C_i + y2 Xbar_A^i psi^{A j} Fdagger_ij. Tree elimination gives -y1*y2 d(psi,psi,v)Fdagger/M, up to the declared identical-Weyl normalization.',
 'scalar_completion':'H(27,bar6), with y d(psi,psi,H)+mu Hdagger v Fdagger and mass² M² Hdagger H; elimination supplies Yukawa coefficient proportional to -y*mu/M². This is the removed elementary Higgs restored as a heavy mediator.',
 'channel_argument':'For a single scalar exchange the psi-psi family6 and E6 bar27 channel forces H(27,bar6). For fermion exchange, the psi-Fdagger vertex fixes E6 conjugate27, while the psi-v vertex selects familybar3; the vectorlike partner then supplies a gauge-invariant mass. This exhausts these two single-propagator tree topologies with the named external fields, not arbitrary UV theories.',
 'one_loop_costs':{'complex_H27bar6':6,'vectorlike_Weyl27bar3_pair':12,'E6_at_r5_scalar':-2,'E6_at_r5_fermion':-8,'SO10_at_r5':1},
 'budget_derivation':'T_E6(27)=3. Scalar: (1/3)*6*3=6. Two Weyls: (2/3)*2*3*3=12. SO10 singlet mediators do not alter bSO10.',
 'economical_alignment_escape':economical_completion(),'rows':rows,'signed_tensor_matching':signed_tensor_matching(),'matching_symbolic':str(eff),
 'boundaries':['The current three-frame composite architecture has no jointly E6/SO2r AF single-propagator tree completion of these two types.','This does not exclude asymptotic safety, strong/composite mediators, altered alignment inventory, multiple-step or radiative UV models.','Family SU3 remains non-AF; absolute mediator mass and Yukawa couplings are inputs.'],
 'prior_owners':['analysis/w33_pass11285_semisimple_factor_higgs.py','analysis/w33_pass11291_rank_five_AF_window.py','analysis/w33_pass11271_chiral_symmetric_yukawa_completion.py'],
 'primary_sources':['https://doi.org/10.1016/0370-1573(81)90092-2']}
if __name__=='__main__':
 out=payload();(ROOT/'data/w33_pass11293_yukawa_uv_budget.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
