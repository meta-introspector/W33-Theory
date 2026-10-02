#!/usr/bin/env python3
"""Explicit symmetric-family Yukawa and spectator mass completion, not fitted SM data."""
from pathlib import Path
import sys,json,itertools,math
import sympy as s
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_20261001_global_e6_cartan_covariants as I

def payload():
    I.generator_certificate();_,d,_=I.tensors()
    # Symmetric Weyl flavor and symmetric d require a bar6 Higgs, not epsilon3.
    Y=s.diag(1,2,3);sl=[]
    for i,j in itertools.permutations(range(3),2):
        x=s.zeros(3);x[i,j]=1;sl.append(x)
    sl.extend([s.diag(1,-1,0),s.diag(0,1,-1)])
    # Compact stabilizer: real antisymmetric plus imaginary symmetric/traceless.
    compact=[]
    for i,j in itertools.combinations(range(3),2):
        a=s.zeros(3);a[i,j]=1;a[j,i]=-1;compact.append(a)
        b=s.zeros(3);b[i,j]=s.I;b[j,i]=s.I;compact.append(b)
    compact.extend([s.I*s.diag(1,-1,0),s.I*s.diag(0,1,-1)])
    # Coefficients of compact generators are REAL. A complex rank would instead
    # find the3-dimensional complex SO3 stabilizer and answer a different question.
    columns=[s.Matrix(list(x.T*Y+Y*x)) for x in compact]
    A=s.Matrix.hstack(*[c.applyfunc(s.re).col_join(c.applyfunc(s.im)) for c in columns]);assert A.rank()==8
    # The full signed E6 Yukawa matrix is symmetric for arbitrary symmetric Y.
    v=np.zeros(27);v[0]=1;M=np.kron(np.einsum('ijk,k->ij',d,v),np.array(Y,int));assert np.array_equal(M,M.T)
    weights=list(itertools.combinations_with_replacement(range(3),3))
    def wick(ix):
        counts=[ix.count(i) for i in range(3)]
        if any(n%2 for n in counts):return 0
        return math.prod(math.prod(range(1,n,2)) for n in counts)
    spectator=s.Matrix([[wick(a+b) for b in weights] for a in weights]);assert spectator.rank()==10
    # SU5 adjoint centralizer: multiplicity 3+2, with traceless condition.
    adj=s.diag(2,2,2,-3,-3);basis=[]
    for i,j in itertools.permutations(range(5),2):
        x=s.zeros(5);x[i,j]=1;basis.append(x)
    for i in range(4):
        x=s.zeros(5);x[i,i]=1;x[4,4]=-1;basis.append(x)
    C=s.Matrix.hstack(*[s.Matrix(list(x*adj-adj*x)) for x in basis]);assert 24-C.rank()==12
    Yu=s.diag(1,2,4);Yd=s.Matrix([[2,1,0],[1,3,1],[0,1,5]])
    comm=Yu*Yu.T*Yd*Yd.T-Yd*Yd.T*Yu*Yu.T;assert comm!=s.zeros(3)
    return {'status':'PASS','result_scope':'PASS_EXPLICIT_SYMMETRIC_CHIRAL_YUKAWA_AND_SPECTATOR_MASS_INTERFACE',
      'fermions':'one Weyl(27,3) plus one Weyl(1,10bar); anomaly27-27=0 as in prior candidate',
      'new_scalars':'(27,6bar) for Yukawa/exotic mass; (1,28) for10bar spectator mass; further Higgs fields required for the named E6-to-SM chain',
      'Yukawa_operator':'d_ABC psi_i^A psi_j^B H^{C,ij}, H^{ij}=H^{ji} in family6bar. Symmetric in the combined fermion indices; does not vanish like the old d times epsilon map.',
      'family_compact_stabilizer_dimension_at_diag123':0,
      'spectator_operator':'chi_ijk chi_lmn Sigma^{ijklmn}, Sigma totally symmetric in6 upper family indices (28). The stored Wick tensor vev gives a full-rank10 mass matrix.',
      'spectator_mass_matrix':[[str(x) for x in row] for row in spectator.tolist()],'spectator_mass_determinant':str(spectator.det()),
      'SM_preserving_exotic_mass':'Under E6->SO10, a scalar1 in27 couples fermion10 times10; its vev masses the vectorlike5+5bar and leaves the chiral16. This is standard branching, not a stabilizer of the regular T vacuum.',
      'named_breaking_chain':'E6->SO10 using27 singlet; SO10->SU5 using16 singlet; SU5->SU3xSU2xU1 using adjoint diag(2,2,2,-3,-3). Alignment and a potential selecting all vevs are not supplied.',
      'explicit_SU5_adjoint_centralizer_dimension':12,'hypercharge_generator':'-diag(2,2,2,-3,-3)/6',
      'Yukawa_test_matrices':{'Yu':[[int(x) for x in row] for row in Yu.tolist()],'Yd':[[int(x) for x in row] for row in Yd.tolist()],'noncommuting_left_Gram':True},
      'updated_beta_subset':{'E6':'32','SU3_family':'-93/2','inventory':'one Weyl81+Weyl10bar; complex scalar(27,6bar)+(1,28), excluding further GUT-breaking scalars. T6=5/2,T28=63.'},
      'boundary':['A constructed anomaly-free operator/Higgs interface, not a dynamically selected chiral Standard Model vacuum.','SM branching is imported from prior E6 work; this script verifies the SU5 centralizer but does not map its generators into the canonical signed27 basis.','Yukawa matrices are independent coefficients; no observed masses, CKM angles or CP phase are predicted.','Extra scalar representations change beta functions and worsen family UV running; the old beta certificate is not inherited.'],
      'prior_owners':['analysis/w33_e6_27_standard_model.py','analysis/w33_20261001_chiral_decuplet_and_condensate_bridge.py','analysis/w33_20261001_global_e6_cartan_covariants.py'],'primary_source':'https://doi.org/10.1016/0370-1573(81)90092-2'}
if __name__=='__main__':
    out=payload();(ROOT/'data/w33_pass11271_chiral_symmetric_yukawa_completion.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
