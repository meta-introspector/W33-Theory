#!/usr/bin/env python3
"""An explicit SM-preserving Higgs orbit in the repository's signed E6 basis.

Standard E6->SO10->SU5->SM branching is prior art. This names its actual
matrix generators and exotic mass block, without a fitted family spectrum.
"""
from pathlib import Path
import sys,json
from collections import Counter
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_20261001_global_e6_cartan_covariants as I

def payload():
    B=np.load(ROOT/'artifacts/e6_27rep_basis_export/E6_basis_78.npy').real.astype(int)
    ids=[j for j,b in enumerate(B) if np.array_equal(b,np.diag(np.diag(b)))];H=[s.Matrix(B[j]) for j in ids]
    weights=s.Matrix([[h[i,i] for h in H] for i in range(27)])
    coeff=s.Matrix([s.Rational(1,3),s.Rational(2,3),1,0,s.Rational(1,2),0]);Y=sum((c*h for c,h in zip(coeff,H)),s.zeros(27))
    assert Y[0,0]==Y[1,1]==0
    roots=[];root_ids=[];su5_roots=[]
    for j,b in enumerate(B):
        if j in ids:continue
        if np.any(b[[0,1],:]) or np.any(b[:,[0,1]]):continue
        su5_roots.append(s.Matrix(b))
        if s.Matrix(b)*Y!=Y*s.Matrix(b):continue
        roots.append(s.Matrix(b));root_ids.append(j)
    assert len(roots)==8
    cs=weights[[0,1],:].nullspace();assert len(cs)==4
    cartan=[sum((c*h for c,h in zip(v,H)),s.zeros(27)) for v in cs]
    assert len(su5_roots)==20
    su5_span=s.SparseMatrix.hstack(*[flat for flat in [s.Matrix(list(a)) for a in su5_roots+cartan]])
    assert su5_span.rank()==24
    generators=roots+cartan
    flat=lambda a:s.Matrix(list(a))
    span=s.SparseMatrix.hstack(*[flat(a) for a in generators]);assert span.rank()==12
    # Closed under Hermitian conjugation, so this complex reductive algebra
    # is the complexification of an actual12-dimensional compact stabilizer.
    assert all(span.row_join(flat(a.T)).rank()==12 for a in generators)
    assert all(not any(a[:,0]) and not any(a[:,1]) and a*Y==Y*a for a in generators)
    brackets=[flat(a*b-b*a) for a in generators for b in generators]
    assert s.SparseMatrix.hstack(*brackets).rank()==11
    # The8 roots in the4D Cartan form A2+A1, plus a1D center.
    rootweights=[]
    for a in roots:
        i,j=next((i,j) for i in range(27) for j in range(27) if a[i,j])
        rootweights.append(tuple(h[i,i]-h[j,j] for h in cartan))
    assert s.Matrix(rootweights).rank()==3
    charge_sizes=[];seen=set()
    for start in range(27):
        if start in seen:continue
        component={start};front=[start]
        while front:
            i=front.pop()
            for a in roots:
                for j in range(27):
                    if (a[i,j] or a[j,i]) and j not in component:component.add(j);front.append(j)
        seen|=component
        assert len({Y[i,i] for i in component})==1
        charge_sizes.append({'indices':sorted(component),'dimension':len(component),'Y':str(Y[start,start])})
    expected=Counter({'1/6':6,'-2/3':3,'1/3':6,'-1/2':4,'1':1,'0':2,'-1/3':3,'1/2':2})
    assert Counter(str(Y[i,i]) for i in range(27))==expected
    _,d,_=I.tensors();mass=s.Matrix(d[:,:,0]);assert mass.rank()==10
    assert Y.T*mass+mass*Y==s.zeros(27)
    assert all(a.T*mass+mass*a==s.zeros(27) for a in generators)
    heavy=sorted({i for i in range(27) if any(mass[i,:])});assert len(heavy)==10
    return {'status':'PASS','result_scope':'PASS_EXACT_CANONICAL_SM_STABILIZER_AND_EXOTIC_MASS',
      'Higgs_reference':'two fundamental27 Higgs vectors e0,e1 plus an adjoint Higgs equal to Y; these are added scalar fields, not the regular Cartan T vacuum',
      'Cartan_basis_export_indices':ids,'hypercharge_Cartan_coefficients':list(map(str,coeff)),
      'hypercharge_diagonal':[str(Y[i,i]) for i in range(27)],'unbroken_root_basis_indices':root_ids,'unbroken_Cartan_coefficients':[list(map(str,v)) for v in cs],
      'unbroken_compact_dimension':12,'derived_algebra_dimension':11,'root_span_rank':3,'center_dimension':1,
      'e0_e1_SU5_stabilizer_dimension':24,'e0_e1_SU5_root_count':20,
      'SM_identification':'Inside the standard SU5 stabilizer of e0,e1, the8 roots form A2+A1; Cartan dimension4 and center1 give su3+su2+u1. Hypercharge branching is verified in the canonical signed27 basis.',
      'SM_components':charge_sizes,'singlet_exotic_mass_rank':10,'exotic_mass_support':heavy,
      'three_family_symmetric_Y_mass_rank':30,'remaining_81_Weyl_components':51,
      'mass_identity':'M_AB=d_AB0; each stored unbroken generator satisfies X^T M+M X=0 exactly. Symmetric family Y=diag(1,2,3) gives M tensorY, rank30, leaving three16s and threeE6 singlets.',
      'named_selecting_potential':'V_H=min_{U in compact E6}(||phi1-Ue0||²+||phi2-Ue1||²+||A-UYUdagger||_F²). Nonnegative and zero precisely on the prescribed compact Higgs orbit; a smooth squared-distance potential in its tubular neighborhood.',
      'boundary':['The orbit-distance potential is an explicitly imposed nonpolynomial Higgs EFT, not a W33-selected renormalizable potential or a naturalness result.','Extra Higgs scalars change UV running; supersymmetric D-flatness would require a separate completion.','This proves a consistent SM-preserving matrix/Higgs and mass interface. It does not predict vev scales, observed Yukawa matrices, mixing or vacuum energy.'],
      'prior_owners':['analysis/w33_e6_27_standard_model.py','analysis/w33_20261001_chiral_mass_assignment_audit.py','analysis/w33_pass11271_chiral_symmetric_yukawa_completion.py'],'primary_source':'https://doi.org/10.1016/0370-1573(81)90092-2'}
if __name__=='__main__':
    out=payload();(ROOT/'data/w33_pass11271_canonical_sm_higgs_bridge.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
