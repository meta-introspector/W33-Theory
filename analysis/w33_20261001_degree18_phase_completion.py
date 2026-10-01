#!/usr/bin/env python3
"""Global degree18 circuit via the Hessian of the E6 ternary cubic.

Classical Aronhold/Hessian construction; new signed81-field normalization.
Exact coefficient restrictions, analytic derivatives, no fitted gauge oracle.
"""
from pathlib import Path
from functools import lru_cache
import sys,itertools,json
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_20261001_global_e6_cartan_covariants as I

@lru_cache(None)
def cubic_circuit():
    keys=list(itertools.combinations_with_replacement(range(3),3));z=s.symbols('f0:10');a=s.symbols('a:3')
    F={ix:z[keys.index(tuple(sorted(ix)))] for ix in itertools.product(range(3),repeat=3)}
    f=sum(F[ix]*s.prod(a[j] for j in ix) for ix in F)
    h=s.Poly(s.det(s.hessian(f,a)),a)
    H={k:h.coeff_monomial(s.prod(a[j] for j in k))/ (s.factorial(3)/s.prod(s.factorial(k.count(j)) for j in range(3))) for k in keys}
    idx,sg=I.aronhold_terms();S=s.expand(sum(int(t)*s.prod(F[tuple(k)] for k in row) for row,t in zip(idx,sg)))
    T=s.expand(sum(s.diff(S,v)*H[k] for v,k in zip(z,keys)))
    calc=s.lambdify(z,[T]+[s.diff(T,v) for v in z],'numpy',cse=True)
    return keys,z,S,T,calc

def evaluate(phi):
    V=np.asarray(phi,complex).reshape(27,3);_,d,_=I.tensors();keys,z,S,T,calc=cubic_circuit()
    F=np.einsum('ijk,ia,jb,kc->abc',d,V,V,V,optimize=True)
    ans=calc(*[F[k] for k in keys]);raw=ans[0];GF=np.zeros((3,3,3),complex)
    for k,v in zip(keys,ans[1:]):
        ps=list(set(itertools.permutations(k)))
        for p in ps:GF[p]=v/len(ps)
    graw=3*np.einsum('abc,ijk,jb,kc->ia',GF,d,V,V,optimize=True).ravel()
    u6,u12,g6,g12=I.evaluate(phi)
    val=(19683*raw-s.Rational(48384,5)*u6**3+s.Rational(13824,5)*u6*u12)/1492992
    grad=(19683*graw-float(s.Rational(145152,5))*u6*u6*g6+float(s.Rational(13824,5))*(u12*g6+u6*g12))/1492992
    return complex(val),grad

def exact_restriction():
    keys,z,S,T,_=cubic_circuit();q=I.G.VARS;V=s.Matrix(27,3,list(s.Matrix(I.C.plane())*s.Matrix(q)))
    F={k:0 for k in keys}
    for i,j,k,sg in I.tensors()[0]:
        for a,b,c in itertools.permutations((i,j,k)):
            for u,v,w in F:F[u,v,w]+=sg*V[a,u]*V[b,v]*V[c,w]
    raw=s.expand(T.subs(dict(zip(z,[F[k] for k in keys])),simultaneous=True))
    u6,u12,u18=I.G.invariants()
    assert s.expand(raw-s.Rational(48384,5)*u6**3+s.Rational(13824,5)*u6*u12-1492992*u18)==0
    return {'raw_T_on_Q':'48384*u6^3/5-13824*u6*u12/5+1492992*u18',
      'I18':'(19683*T-48384*I6^3/5+13824*I6*I12/5)/1492992',
      'canonical_restriction':'I18(Qq/sqrt3)=u18(q)','classical_T_coefficient_terms':len(s.Poly(T,z).terms()),
      'construction':'T=dS_F[H_F], H_F is the coefficient tensor of det Hessian(Fabc a_a a_b a_c)'}

def payload():
    exact=exact_restriction();I.generator_certificate();rng=np.random.default_rng(1801)
    p=.15*(rng.normal(size=81)+1j*rng.normal(size=81));v=rng.normal(size=81)+1j*rng.normal(size=81);v/=np.linalg.norm(v)
    t,g=evaluate(p);h=1e-5;err=abs((evaluate(p+h*v)[0]-evaluate(p-h*v)[0])/(2*h)-g@v);assert err<1e-7
    E=I.C.plane()/np.sqrt(3);z=np.exp(2j*np.pi/9);q=np.array([1,z,z.conjugate()])/np.sqrt(3)
    value,grad=evaluate(E@q);assert abs(value+1/729)<1e-12
    assert abs(grad@(1j*E@q)-18j*value)<1e-12
    q2=np.array([.2+.1j,-.4+.2j,.6-.3j]);u=complex(I.G.invariants()[2].subs(dict(zip(I.G.VARS,q2))))
    assert abs(evaluate(E@q2)[0]-u)<1e-9
    b=np.load(ROOT/'artifacts/e6_27rep_basis_export/E6_basis_78.npy').real[0];from scipy.linalg import expm
    U=np.kron(expm(.04*b),np.diag(np.exp([.1,-.04,-.06])))
    t2,g2=evaluate(U@p);ge=abs(t2-t);ce=float(np.max(abs(g2-np.linalg.solve(U.T,g))));assert max(ge,ce)<1e-8
    return {'status':'PASS_GLOBAL_I18_AND_PHASE_LIFT','exact':exact,'numerical_controls':{'complex_derivative_error':float(err),'gauge_invariance_error':float(ge),'gradient_transport_error':ce},
      'phase_selector':{'term':'gamma*|I18(Phi)-tau|^2, tau=-r0^18/729 real, gamma>0',
       'vacua':'At the prior norm/I6/I12/moment-zero minima, exp(18i theta)=1. All resulting representatives share the invariant triple and lie in one G26, hence compact gauge, orbit; no eighteen physical vacua are claimed.',
       'canonical_phase_mass_squared':'324*gamma*r0^34/729^2 (kinetic phase coordinate sqrt2*r0*theta)',
       'CP_boundary':'T and its conjugate have the same real invariant triple (0,0,-r0^18/729); the prior finite G26 quotient identifies their gauge orbits. This real phase selector does not establish physical spontaneous CP violation.',
       'inputs':'tau is an imposed dimensionful EFT coefficient, not a W33-selected clock or measured CP phase.'},
      'prior_owners':['analysis/w33_20261001_global_e6_cartan_covariants.py','analysis/w33_20261001_balanced_g26_tmagic_vacuum.py','analysis/w33_pass11255_11259_g26_common.py'],
      'external_source':'https://people.dimai.unifi.it/ottaviani/tesi/tesidr_maurizio.pdf'}
if __name__=='__main__':
    out=payload();(ROOT/'data/w33_20261001_degree18_phase_completion.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
