#!/usr/bin/env python3
"""Rank-five fundamental factor lift: exact local kernel and audited gauge cost."""
from pathlib import Path
import sys,json
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11275_polynomial_sm_higgs as H

def reference():
    B,v,A,d=H.setup();N=d[:,:,0].T@d[:,:,1];rows=np.where(np.any(N,axis=1))[0]
    U=np.eye(27,dtype=complex)[:,rows];V=N.T@U;W=A@U
    return v,A,U,V,W

def constraints(v,A,U,V,W,Lambda=1.):
    _,_,_,d=H.setup();M0=np.einsum('abc,c->ab',d,v[:,0]);M1=np.einsum('abc,c->ab',d,v[:,1]);N=M0.conj().T@M1
    return U@V.conj().T-N,U.conj().T@U-np.eye(5),V.conj().T@V-np.eye(5),Lambda*W-A@U,(A+np.eye(27))@W-6*U/Lambda

def payload():
    v,A,U,V,W=reference();assert all(not np.any(z) for z in constraints(v,A,U,V,W))
    old=json.loads((ROOT/'data/w33_pass11275_polynomial_sm_higgs.json').read_text());assert old['constraint_rank_mod101']==120
    N=U@V.conj().T;assert not np.any(N.imag) and H.rank_mod(N.real)==5
    e6=s.Integer(28)-s.Rational(1,3)*15*3;aux=s.Rational(11,3)*5-s.Rational(1,3)*3*27*s.Rational(1,2)
    assert e6==13 and aux==s.Rational(29,6)
    return {'status':'PASS','result_scope':'PASS_EXACT_LOCAL_FUNDAMENTAL_HIGGS_FACTOR_LIFT',
      'fields':'Three complex27x5 fields U,V,W; E6 acts on the left and an added U(5) gauge group acts on their common right frame. All are family-neutral. Keep the original two27 Higgs vectors and real78 adjoint.',
      'potential':'Retain ||vdagger v-I2||²+sum||d(v_i,v_j)||²+||Av||². Add squared norms of UVdagger-N, Udagger U-I5, Vdagger V-I5, Lambda W-AU, (A+aI)W+(b/Lambda)U; N=d(v0)dagger d(v1), a=1,b=-6 in reference units.',
      'power_counting':'All constraints are at most quadratic in elementary scalar fields; the potential has degree at most4, with appropriate dimensionful constants.',
      'reference_factor_image_indices':[17,18,19,20,22],'reference_U':U.real.astype(int).tolist(),'reference_V':V.real.astype(int).tolist(),
      'E6_beta_combined':str(e6),'aux_SU5_beta':str(aux),'aux_U1_beta_charge1':'-135',
      'inventory':'The factor fields are15 complex27s, total E6 scalar index45 instead of324 for two End27 matrices. The auxiliary SU5 sees81 fundamental complex scalars, index81/2. Its central U1 is normalized with scalar charge1 and remains non-asymptotically free.',
      'real_scalar_dimensions':996,'gauge_orbit_dimension':91,'positive_normal_directions':905,
      'kernel_proof':'The last two equations imply(A²+aA+bI)U=0; multiplying by Vdagger and linearizing UVdagger=N recovers the previous degree8 selector and its exact66-dimensional E6 orbit kernel. Subtract that gauge variation. With delta(v,A)=0, linearizing UVdagger=N and the two orthonormality constraints gives deltaU=UK, deltaV=VK, Kdagger=-K; deltaW=WK follows. These25 directions are exactly the additional U5 gauge frame. Thus the kernel is66+25=91, and the sum-of-squares Hessian is positive in905 normal directions. No numerical996-dimensional Hessian is asserted.',
      'stabilizer':'The surviving12-dimensional SM algebra is diagonal between E6 and the appropriate U5 action on the two five-frames; matter neutral under the auxiliary group retains its canonical E6 hypercharges.',
      'coupling_boundary':'Diagonal gauge breaking changes low-energy coupling matching. The auxiliary coupling and all vev scales remain inputs; the U1 running remains a UV problem.',
      'rank_lower_bound':'At least5 factor columns are required within this architecture because rank(UVdagger)<=number of columns and the canonical N has rank5. No global minimality over all possible models is claimed.',
      'global_boundary':'The added partial-isometry conditions can alter distant zero orbits. Local isolation is proved; global equivalence to the old potential is not.',
      'prior_owners':['analysis/w33_pass11275_polynomial_sm_higgs.py','analysis/w33_pass11271_canonical_sm_higgs_bridge.py']}
if __name__=='__main__':
    out=payload();(ROOT/'data/w33_pass11280_factor_higgs_mediator.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
