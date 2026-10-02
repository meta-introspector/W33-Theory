#!/usr/bin/env python3
"""A one-sextet renormalizable obstruction, explicit higher-degree selector and bilinear repair."""
from pathlib import Path
import sys,json
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
from w33_pass11275_polynomial_sm_higgs import rank_mod

def potential(F):
    S=F@F.conj().T;c=[14,98,794]
    return float(sum(abs(np.trace(np.linalg.matrix_power(S,k))-v)**2 for k,v in enumerate(c,1))+abs(np.linalg.det(F)-6)**2)

def payload():
    x,y,z,m,lam,k,mu=s.symbols('x y z m lambda kappa mu',real=True)
    p=-m*(x*x+y*y+z*z)+lam*(x*x+y*y+z*z)**2+k*(x**4+y**4+z**4)-2*mu*x*y*z
    eq=[s.diff(p,a)/(2*a) for a in (x,y,z)]
    assert s.factor(eq[0]-eq[1])==(x-y)*(x+y)*(2*k*x*y+mu*z)/(x*y)
    # For three distinct positive singular values, subtracting all pairs forces
    # 2k+mu*z/(xy)=2k+mu*x/(yz)=0 =>mu*(z²-x²)=0 =>mu=k=0.
    assert s.simplify((eq[0]-eq[1])/(x*x-y*y)-(eq[1]-eq[2])/(y*y-z*z)-mu*(z-x)*(z+x)/(x*y*z))==0
    F=np.diag([1.,2.,3.]).astype(complex);dirs=[]
    for i in range(3):
        for j in range(i,3):
            a=np.zeros((3,3),complex);a[i,j]=a[j,i]=1;dirs.extend([a,1j*a])
    # Exact differential of the three spectral moments plus complex determinant.
    S=F@F.conj().T;cols=[]
    for a in dirs:
        dS=a@F.conj().T+F@a.conj().T
        det=6*np.trace(np.linalg.solve(F,a))
        col=[int(round((t*np.trace(np.linalg.matrix_power(S,t-1)@dS)).real)) for t in (1,2,3)]+[int(round(det.real)),int(round(det.imag))]
        cols.append(col)
    J=np.array(cols).T
    compact=[]
    for i,j in [(0,1),(0,2),(1,2)]:
        a=np.zeros((3,3),complex);a[i,j]=1;a[j,i]=-1;compact.append(a)
        a=np.zeros((3,3),complex);a[i,j]=a[j,i]=1j;compact.append(a)
    compact.extend([1j*np.diag([1,-1,0]),1j*np.diag([0,1,-1])])
    T=[]
    for a in compact:
        df=a@F+F@a.T;T.append([int(df[i,j].real) if phase==0 else int(df[i,j].imag) for i in range(3) for j in range(i,3) for phase in (0,1)])
    T=np.array(T).T
    assert rank_mod(J)==4 and rank_mod(T)==8 and not np.any(J@T) and potential(F)<1e-20
    # Mixed E6 singlet avoids the previous d(Hdagger)^3 zero.
    v=np.eye(27)[:,0];H=np.zeros((27,3,3),complex);H[0]=F.conj();composite=np.einsum('a,aij->ij',v,H.conj());assert np.array_equal(composite,F)
    return {'status':'PASS','result_scope':'PASS_EXACT_FAMILY_SELECTOR_AND_RENORMALIZABLE_SINGLE_SEXTET_OBSTRUCTION',
      'renormalizable_single_sextet_potential':'-m² Tr(FFdagger)+lambda[Tr(FFdagger)]²+kappa Tr((FFdagger)²)-[mu detF+h.c.], the SU3 invariant terms through degree4.',
      'stationarity_obstruction':'At a full-rank stationary point the determinant phase aligns or anti-aligns. For positive Takagi singular values, pair subtraction gives(si²-sj²)(2kappa+mu_eff sk/(si sj))=0. Three distinct singular values force mu_eff=kappa=0; then all fixed-norm shapes are flat. Thus the single-sextet tree-level renormalizable polynomial potential cannot isolate diag(1,2,3).',
      'polynomial_selector':'sum_k=1..3|Tr[(FFdagger)^k]-c_k|²+gamma|detF-6|², c=(14,98,794), gamma>0. Maximal field degree12.',
      'exact_rank_argument':'Integer J T=0; rank_mod101(J)=4 and rank_mod101(T)=8 in12 real fields give matching rational upper/lower rank bounds.',
      'exact_normal_rank':4,'real_sextet_dimension':12,'compact_SU3_orbit_dimension':8,
      'vacuum_scope':'The spectral moments fix eigenvalues1,4,9; detF fixes its phase. Takagi factorization identifies one compact SU3 orbit of diag(1,2,3), permutations included. The normal Hessian is positive; the8 orbit directions are gauge or global Goldstones according to the declared theory.',
      'bilinear_composite':'F_eff^{ij}=v^A Hdagger_A^{ij}, using the existing family-neutral27 Higgs v and H(27,6bar). This is an E6 singlet in family6 and is nonzero at the canonical SM point v=e0, H=e0 tensor diag(1,2,3).',
      'composite_spectator_operator':'chi chi [Sym(F_eff^3)]/M_UV^5, dimension9; an alternative to adding an elementary sextet. Its rank-ten mass follows from11276, but its extra suppression and UV matching are inputs.',
      'renormalizable_matching_if_F_added':'||Lambda F-v Hdagger||² is a degree4 invariant scalar coupling. It fixes F to the mixed bilinear at a zero without using the vanishing E6 cubic of three SM-singlet Hdaggers.',
      'Yukawa_boundary':'If Yu=a Fdagger and Yd=b Fdagger, their left Grams commute exactly: mixing angles vanish. Three distinct singular values alone do not produce CKM mixing or physical CP; an additional misaligned spurion/field and dynamics are required.',
      'inputs':'The moment targets, determinant target and dimensionful scales are imposed. This selects a tested vacuum in a declared EFT but does not predict observed family mass ratios.',
      'prior_owners':['analysis/w33_pass11276_sextet_composite_spectator.py','analysis/w33_pass11271_chiral_symmetric_yukawa_completion.py']}
if __name__=='__main__':
    out=payload();(ROOT/'data/w33_pass11281_family_sextet_selection.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
