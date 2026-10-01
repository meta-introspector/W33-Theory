#!/usr/bin/env python3
"""Explicit global81-field I6/I12 polynomial circuits and analytic differentials.

Uses the canonical signed E6 cubic and its exact dual, not a gauge-fixing oracle.
Restriction to the new sparse plane is proved symbolically; all E6 generators
preserve both cubic tensors exactly. Evaluators use complex floating arithmetic.
"""
from pathlib import Path
from functools import lru_cache
import sys,json,argparse,itertools
import numpy as np
import sympy as s
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_20261001_exact_cartan_mirror_completion as C
import w33_pass11255_11259_g26_common as G
from w33_20261001_isovolume_dirac_gravity_response import compare_certificate
OUT=ROOT/'data/w33_20261001_global_e6_cartan_covariants.json'

@lru_cache(None)
def tensors():
    mod=C.P._load(ROOT/'tools/toe_e8_z3graded_bracket_jacobi.py','global_cubic_source')
    tri=mod._load_signed_cubic_triads();d=np.zeros((27,27,27),int);eps=np.zeros((3,3,3),int)
    for i,j,k,sg in tri:
        for ix in itertools.permutations((i,j,k)):d[ix]=sg
    for p in itertools.permutations(range(3)):
        eps[p]=(-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
    return tri,d,eps

@lru_cache(None)
def aronhold_terms():
    _,_,eps=tensors();perms=list(itertools.permutations(range(3)));indices=[];signs=[]
    for p1,p2,p3,p4 in itertools.product(perms,repeat=4):
        a,d,g=p1;b,e,j=p2;c,h,k=p3;f,i,l=p4
        indices.append(((a,b,c),(d,e,f),(g,h,i),(j,k,l)))
        signs.append(eps[p1]*eps[p2]*eps[p3]*eps[p4])
    return np.array(indices),np.array(signs)

def evaluate(phi):
    """Return I6,I12 and holomorphic covector derivatives on all81 coordinates."""
    V=np.asarray(phi,dtype=complex).reshape(27,3);tri,d,eps=tensors()
    X=np.einsum('ijk,ja,kb->iab',d,V,V,optimize=True)
    raw6=0j;GX=np.zeros_like(X)
    mixed=lambda a,b,c:np.einsum('ij,kl,mn,ikm,jln->',a,b,c,eps,eps,optimize=True)
    adj=lambda b,c:np.einsum('ikm,jln,kl,mn->ij',eps,eps,b,c,optimize=True)
    for i,j,k,sg in tri:
        raw6+=6*sg*mixed(X[i],X[j],X[k])
        GX[i]+=6*sg*adj(X[j],X[k]);GX[j]+=6*sg*adj(X[i],X[k]);GX[k]+=6*sg*adj(X[i],X[j])
    g6raw=2*np.einsum('ijk,iab,kb->ja',d,GX,V,optimize=True)
    F=np.einsum('ijk,ia,jb,kc->abc',d,V,V,V,optimize=True)
    idx,signs=aronhold_terms();vals=np.array([F[tuple(idx[:,j].T)] for j in range(4)])
    raw12=np.sum(signs*np.prod(vals,axis=0));GF=np.zeros_like(F)
    for j in range(4):np.add.at(GF,tuple(idx[:,j].T),signs*np.prod(np.delete(vals,j,axis=0),axis=0))
    g12raw=(np.einsum('abc,ijk,jb,kc->ia',GF,d,V,V,optimize=True)
      +np.einsum('abc,ijk,ia,kc->jb',GF,d,V,V,optimize=True)
      +np.einsum('abc,ijk,ia,jb->kc',GF,d,V,V,optimize=True))
    I6=9*raw6/4;g6=9*g6raw/4
    I12=(152*I6**2-3645*raw12)/32
    g12=(304*I6*g6-3645*g12raw)/32
    return I6,I12,g6.ravel(),g12.ravel()

def generator_certificate():
    _,d,eps=tensors();B=np.load(ROOT/'artifacts/e6_27rep_basis_export/E6_basis_78.npy')
    assert not np.any(B.imag) and np.max(abs(B.real-np.rint(B.real)))==0
    B=B.real.astype(int)
    def act(t,b):
        return (np.einsum('ai,ajk->ijk',b,t)+np.einsum('aj,iak->ijk',b,t)+np.einsum('ak,ija->ijk',b,t))
    for b in B:assert not np.any(act(d,b)) and not np.any(act(d,b.T))
    sl=[]
    for i in range(3):
        for j in range(3):
            if i!=j:
                b=np.zeros((3,3),int);b[i,j]=1;sl.append(b)
    sl+=[np.diag([1,-1,0]),np.diag([0,1,-1])]
    for b in sl:assert not np.any(act(eps,b))
    return {'E6_lower_and_dual_generator_checks':156,'SL3_volume_form_checks':8,
      'invariance_proof':'I6 contracts three dual-valued cross products with the contravariant E6 cubic and two SL3 volume forms. The Aronhold scalar contracts four E6-invariant ternary cubics with four volume forms. All canonical tensor generators annihilate the defining tensors exactly.'}

def restriction_certificate():
    tri,d,eps=tensors();q=G.VARS;V=s.Matrix(27,3,list(s.Matrix(C.plane())*s.Matrix(q)))
    X=[s.zeros(3) for _ in range(27)];F={ix:0 for ix in itertools.product(range(3),repeat=3)}
    for i,j,k,sg in tri:
        for a,b,c in itertools.permutations((i,j,k)):
            for u,v in itertools.product(range(3),repeat=2):X[a][u,v]+=sg*V[b,u]*V[c,v]
            for u,v,w in F:F[u,v,w]+=sg*V[a,u]*V[b,v]*V[c,w]
    raw6=0
    for i,j,k,sg in tri:
        for p1,p2 in itertools.product(itertools.permutations(range(3)),repeat=2):
            a,b,c=p1;e,f,g=p2
            raw6+=6*sg*int(eps[p1]*eps[p2])*X[i][a,e]*X[j][b,f]*X[k][c,g]
    idx,signs=aronhold_terms();raw12=sum(int(sg)*s.prod(F[tuple(ix)] for ix in row) for row,sg in zip(idx,signs))
    u6,u12,_=G.invariants()
    assert s.expand(raw6-12*u6)==0
    assert s.expand(raw12-(152*u6*u6-32*u12)/5)==0
    return {'unnormalized_Q_restriction':['raw6=12u6','AronholdS=(152u6^2-32u12)/5'],
      'canonical_E_Q_over_sqrt3_normalizations':['I6=(9/4)raw6','I12=(152I6^2-3645AronholdS)/32'],
      'exact_restrictions':'I6(Eq)=u6(q), I12(Eq)=u12(q), E=Q/sqrt(3)',
      'coefficient_circuit':'45 signed E6 triads, canonical3D epsilon;1296 signed four-F Aronhold terms. No expanded monomial table is needed to evaluate the global polynomials or gradients.'}

def payload():
    gc=generator_certificate();rc=restriction_certificate();rng=np.random.default_rng(11270)
    phi=.08*(rng.normal(size=81)+1j*rng.normal(size=81));a,b,g6,g12=evaluate(phi)
    derivative_error=0.
    for _ in range(4):
        v=rng.normal(size=81)+1j*rng.normal(size=81);v/=np.linalg.norm(v);h=1e-5
        plus=evaluate(phi+h*v);minus=evaluate(phi-h*v)
        derivative_error=max(derivative_error,abs((plus[0]-minus[0])/(2*h)-g6@v),abs((plus[1]-minus[1])/(2*h)-g12@v))
    assert derivative_error<1e-7
    B=np.load(ROOT/'artifacts/e6_27rep_basis_export/E6_basis_78.npy').real
    L=expm(.04*B[0]);R=np.diag(np.exp([.08,-.03,-.05]));U=np.kron(L,R.T)
    aa,bb,gg6,gg12=evaluate(U@phi)
    inv=np.linalg.inv(U);invariance=max(abs(aa-a),abs(bb-b))
    covariance=max(float(np.max(abs(gg6-inv.T@g6))),float(np.max(abs(gg12-inv.T@g12))))
    assert max(invariance,covariance)<1e-9
    E=C.plane()/np.sqrt(3);z=np.exp(2j*np.pi/9);q=np.array([1,z,z.conjugate()])/np.sqrt(3)
    vals=evaluate(E@q);local,_=C.gradient_data()
    grad_error=max(float(np.max(abs(vals[j+2]-E@local[j]))) for j in range(2))
    assert grad_error<1e-10 and abs(vals[0])<1e-10 and abs(vals[1])<1e-10
    return {'schema':'w33.20261001.global-e6-cartan-covariants.v1',
      'status':'PASS_EXPLICIT_GLOBAL_I6_I12_CIRCUITS_EXACT_TENSOR_INVARIANCE_AND_RESTRICTION',
      'defining_tensors':{'signed_E6_triads':[list(t) for t in tensors()[0]],'epsilon':'epsilon012=+1, fully alternating',
        'cross_product':'X_a(alpha,beta)=d_aij Phi_i(alpha)Phi_j(beta)',
        'raw6':'d^abc epsilon(alpha,gamma,eta)epsilon(beta,delta,theta)X_a(alpha,beta)X_b(gamma,delta)X_c(eta,theta)',
        'F':'F(alpha,beta,gamma)=d_ijk Phi_i(alpha)Phi_j(beta)Phi_k(gamma)',
        'AronholdS':'Fabc Fdef Fghi Fjkl epsilon(adg)epsilon(bej)epsilon(chk)epsilon(fil)'},
      'exact_generator_invariance':gc,'exact_restriction':rc,
      'global_evaluator':{'domain':'every complex81-vector, including non-Dflat and nonregular fields',
        'differentials':'analytic polynomial chain rule; no finite-difference derivative implementation',
        'coefficient_scope':'integer/rational contraction circuit; exact generator and symbolic restriction proofs, floating evaluation controls'},
      'numerical_controls':{'nonplane_complex_derivative_error':float(derivative_error),'complexified_gauge_invariance_error':float(invariance),'global_covector_transport_error':covariance,'T_normal_slice_gradient_error':grad_error},
      'mass_completion_boundary':'Global operators are now explicit and gauge invariant. Exact heavy/null separation and three specified light masses concern the D-flat T plane and its compact unitary orbit, not every background.',
      'full_classical_selector':{'potential':'lambda(||Phi||^2-r0^2)^2+alpha|I6(Phi)|^2+beta|I12(Phi)|^2+kappa sum moment_a(Phi)^2, all coefficients positive',
        'global_proof':'All terms are nonnegative; EqT at radiusr0 has value0. Every zero has moment0 and is semisimple. Conjugacy of Cartan subspaces and uniqueness of a moment-zero compact orbit representative reduce it to the prior exact72-ray Cartan intersection.',
        'gauge_boundary':'The72 Cartan rays are little-Weyl-related representatives of gauge orbits, not72 distinct physical vacua or prepared magic states. A common phase remains; quantum stability and state preparation are additional problems.',
        'inputs':'Uses prior exact72-ray intersection and theta/closed-orbit theorems. Radius/couplings are external; chosen zero energy is not a cosmological-constant prediction.'},
      'prior_owners':['analysis/w33_20261001_exact_cartan_mirror_completion.py','analysis/w33_pass11269_g26_qutrit_dictionary.py','tools/toe_e8_z3graded_bracket_jacobi.py','analysis/w33_pfaffian_doily_e6_cubic_bridge.py'],
      'physical_boundary':'Coefficients, radial scales, particle assignments and full higher-operator loops remain inputs/open. Aronhold contractions are classical invariant theory; the executable signed81-field identification and differential circuit are the result here.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args();out=payload()
    if args.check:assert compare_certificate(out,json.loads(OUT.read_text()))
    else:OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(out['status'])
