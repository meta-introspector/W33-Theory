#!/usr/bin/env python3
"""Polynomial local alignment of the canonical SM orbit; no global uniqueness."""
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_20261001_global_e6_cartan_covariants as I

def rank_mod(a,p=101):
    a=np.array(a,dtype=np.int64)%p;m,n=a.shape;r=0
    for c in range(n):
        ix=next((j for j in range(r,m) if a[j,c]),None)
        if ix is None:continue
        a[[r,ix]]=a[[ix,r]];a[r]=a[r]*pow(int(a[r,c]),-1,p)%p
        for j in range(r+1,m):
            if a[j,c]:a[j]=(a[j]-a[j,c]*a[r])%p
        r+=1
        if r==m:break
    return r

def setup():
    B=np.load(ROOT/'artifacts/e6_27rep_basis_export/E6_basis_78.npy').real.astype(int)
    H=[];seen=set()
    for j,b in enumerate(B):
        if j in seen:continue
        if np.array_equal(b,b.T):H.append(b.astype(complex));seen.add(j)
        else:
            k=next(k for k,c in enumerate(B) if np.array_equal(c,b.T) or np.array_equal(c,-b.T))
            seen|={j,k};H.extend([b+b.T,1j*(b-b.T)])
    H=np.array(H);assert len(H)==78
    ids=[j for j,b in enumerate(B) if np.array_equal(b,np.diag(np.diag(b)))]
    y=sum(c*np.diag(B[j]) for c,j in zip([2,4,6,0,3,0],ids))
    return H,np.eye(27,dtype=complex)[:,:2],np.diag(y),I.tensors()[1]

def jacobians():
    H,v,Y,d=setup();M0=d[:,:,0];M1=d[:,:,1];N=M0.T@M1;R=Y@Y+Y-6*np.eye(27)
    assert not np.any(R@N)
    def differential(dv,a):
        gram=v.conj().T@dv+dv.conj().T@v
        cubic=np.concatenate([np.einsum('abc,b,c->a',d,dv[:,i],v[:,j])+np.einsum('abc,b,c->a',d,v[:,i],dv[:,j]) for i,j in [(0,0),(0,1),(1,1)]])
        dm0=np.einsum('abc,c->ab',d,dv[:,0]);dm1=np.einsum('abc,c->ab',d,dv[:,1])
        dn=dm0.conj().T@M1+M0.T@dm1
        select=(Y@a+a@Y+a)@N+R@dn
        r=np.r_[gram.ravel(),cubic,(Y@dv+a@v).ravel(),select.ravel()]
        return np.r_[r.real,r.imag].astype(int)
    cols=[]
    for c in range(108):
        dv=np.zeros((27,2),complex);dv.ravel()[c%54]=1 if c<54 else 1j
        cols.append(differential(dv,np.zeros((27,27))))
    for a in H:cols.append(differential(np.zeros((27,2)),a))
    J=np.array(cols).T;J=J[np.any(J,axis=1)]
    # Coordinates of actual compact orbit tangent: dv=iXv, dA=i[X,Y].
    Hb=np.column_stack([np.r_[h.real.ravel(),h.imag.ravel()] for h in H])
    tang=[]
    for h in H:
        dv=1j*h@v;da=1j*(h@Y-Y@h)
        coords=np.linalg.lstsq(Hb,np.r_[da.real.ravel(),da.imag.ravel()],rcond=None)[0]
        assert np.max(abs(coords-np.rint(coords)))<1e-10
        assert np.array_equal(Hb@np.rint(coords),np.r_[da.real.ravel(),da.imag.ravel()])
        tang.append(np.r_[dv.real.ravel(),dv.imag.ravel(),np.rint(coords)].astype(int))
    T=np.array(tang).T;assert not np.any(J@T)
    return J,T

def potential(v,A):
    _,_,Y,d=setup();M0=np.einsum('abc,c->ab',d,v[:,0]);M1=np.einsum('abc,c->ab',d,v[:,1])
    select=(A@A+A-6*np.eye(27))@M0.conj().T@M1
    cubic=[np.einsum('abc,b,c->a',d,v[:,i],v[:,j]) for i,j in [(0,0),(0,1),(1,1)]]
    return float(np.linalg.norm(v.conj().T@v-np.eye(2))**2+sum(np.linalg.norm(x)**2 for x in cubic)+np.linalg.norm(A@v)**2+np.linalg.norm(select)**2)

def lifted_constraints(v,A,X,Z,Lambda=1.):
    # X,Z are complex End(27) scalars, each1+78+650 under E6.
    _,_,_,d=setup();M0=np.einsum('abc,c->ab',d,v[:,0]);M1=np.einsum('abc,c->ab',d,v[:,1])
    N=M0.conj().T@M1
    return Lambda*X-N,Lambda*Z-A@X,(A+np.eye(27))@Z-6*X/Lambda

def payload():
    H,v,Y,d=setup();J,T=jacobians();rj=rank_mod(J);rt=rank_mod(T)
    assert (rj,rt)==(120,66) and rj+rt==186
    assert potential(v,Y)==0
    N=d[:,:,0].T@d[:,:,1];image=np.where(np.any(N,axis=1))[0].tolist()
    assert rank_mod(N)==5 and image==[17,18,19,20,22]
    assert set(Y.diagonal()[image])=={2,-3}
    assert all(not np.any(z) for z in lifted_constraints(v,Y,N,Y@N))
    return {'status':'PASS','result_scope':'PASS_EXACT_LOCAL_POLYNOMIAL_SM_ALIGNMENT',
      'potential':'V=||Vdagger V-I2||²+sum_i<=j||d(v_i,v_j)||²+||AV||²+||(A²+A-6I)N||_F²; N=d(v0)dagger d(v1), A Hermitian in compact e6. Positive dimensionful coefficients/cutoff factors are understood.',
      'cubic_selected_block':{'image_indices':[17,18,19,20,22],'image_dimension':5,'adjoint_eigenvalues_on_image':[2,-3],'covariance':'N transforms to U N Udagger under compact E6, so the last term is invariant. The reference quadratic is (A-2I)(A+3I).'},
      'reference':'V=(e0,e1), A=6Y with Y from Pass11271. Fundamental fields and adjoint are family-neutral extra Higgs fields.',
      'real_field_dimension':186,'constraint_Jacobian_shape':list(J.shape),'constraint_rank_mod101':rj,'compact_orbit_rank_mod101':rt,'exact_integer_J_times_T_zero':True,
      'exact_rank_argument':'A nonzero120-minor modulo101 gives rank_Q J>=120; a nonzero66-minor gives rank_Q T>=66; integer JT=0 bounds rank_Q J<=120 and rank_Q T<=66. Hence both ranks exact. The complete quadratic selector differential includes variations of N, not merely fixed-block variations.',
      'local_theorem':'The real Hessian of the nonnegative sum of squares has kernel exactly the66-dimensional E6 gauge orbit and is positive on its120-dimensional normal space. The compact slice theorem and a full-rank constraint minor isolate this orbit locally among zero-energy minima.',
      'polynomial_degree':8,'SM_stabilizer_dimension':12,
      'renormalizable_lift':{'new_scalars':'two complex End(27)=1+78+650 fields X,Z, family neutral',
        'replaces_last_term':'||Lambda X-N||²+||Lambda Z-AX||²+||(A+a I)Z+(b/Lambda)X||², with a=1,b=-6 in reference units. Each constraint is at most quadratic in fields, hence potential degree at most4.',
        'zero_locus':'X=N/Lambda, Z=AN/Lambda²; the remaining equation is (A²+aA+bI)N=0. This is exact for zeros, not an equality of off-shell effective potentials.',
        'real_field_dimension':3102,'gauge_kernel_dimension':66,'positive_normal_directions':3036,
        'local_rank_proof':'The first two linearized constraints solve deltaX and deltaZ uniquely in terms of deltaV,deltaA. The third reduces to the degree8 selector differential. Thus the previous exact kernel theorem extends; no numerical3102-dimensional Hessian is claimed.',
        'E6_Dynkin_index_each_End27':162,'E6_beta_combined_inventory':'-80',
        'UV_cost':'Two complex End27 scalars contribute-108 to bE6; the real78 Higgs inventory had bE6=28. The renormalizable scalar construction therefore loses E6 asymptotic freedom and is not a viable UV completion by itself.'},
      'boundary':['Local isolation only; remote disconnected zero orbits are not excluded.','Positive coefficients and spectral roots are imposed. No W33-derived scale or naturalness claim. The degree4 lifted scalar model is renormalizable by power counting, but its huge inventory worsens UV running.','The last term is dimension8 in fields and needs appropriate EFT cutoff factors; changing positive coefficients preserves the local rank theorem.'],
      'prior_owners':['analysis/w33_pass11271_canonical_sm_higgs_bridge.py'],
      'primary_sources':['https://arxiv.org/abs/1504.00904']}
if __name__=='__main__':
    out=payload();(ROOT/'data/w33_pass11275_polynomial_sm_higgs.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'],out['constraint_rank_mod101'])
