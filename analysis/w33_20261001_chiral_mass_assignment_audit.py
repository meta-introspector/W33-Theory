#!/usr/bin/env python3
"""Exact little-group action on the three completed cubic null modes.

Audits the declared vectorlike four-Weyl model, not every E6 construction.
The earlier standard27 branching remains prior art, not a chirality proof here.
"""
from pathlib import Path
import sys,json,itertools
from collections import Counter
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_20261001_exact_cartan_mirror_completion as C

def stabilizer():
    B=np.load(ROOT/'artifacts/e6_27rep_basis_export/E6_basis_78.npy').real.astype(int);sl=[]
    for i,j in itertools.permutations(range(3),2):
        a=np.zeros((3,3),int);a[i,j]=1;sl.append(a)
    sl +=[np.diag([1,-1,0]),np.diag([0,1,-1])]
    Rs=[s.SparseMatrix(np.kron(b,np.eye(3,dtype=int))) for b in B]+[s.SparseMatrix(np.kron(np.eye(27,dtype=int),a.T)) for a in sl]
    Q=s.Matrix(C.plane());p=Q*s.Matrix([1,2,3]);A=s.SparseMatrix.hstack(*[r*p for r in Rs]);ns=A.nullspace()
    assert len(ns)==8
    acts=[]
    for c in ns:
        r=sum((v*b for v,b in zip(c,Rs)),s.zeros(81));assert r*Q==s.zeros(81,3)
        acts.append({'coefficients':[str(v) for v in c],'restriction_to_three_null_modes':'zero3x3'})
    _,op=C.operators();assert s.Matrix(op(np.array(p,dtype=int).ravel())).rank()==78
    return {'reference_point':[1,2,3],'exact_gauge_orbit_rank':78,'stabilizer_dimension':8,'generators':acts,
      'T_extension':'The eight exact generators annihilate the entire Q plane. At regular T the prior exact gauge-mass identity and mirror nonvanishing give rank78 again, so they exhaust its unbroken Lie algebra.',
      'null_mode_charges':'All three normal directions are singlets of the unbroken connected gauge group.',
      'SM_obstruction':'Unbroken SU3xSU2xU1 has dimension12; this regular vacuum has only8 unbroken generators. The three light modes cannot be three charged Standard Model families in this vacuum.'}

def branching():
    # (SU3 rep, SU2 dimension, rational hypercharge); standard 16+10+1.
    rows=[('3',2,s.Rational(1,6)),('3bar',1,-s.Rational(2,3)),('3bar',1,s.Rational(1,3)),('1',2,-s.Rational(1,2)),('1',1,s.Integer(1)),('1',1,s.Integer(0)),
          ('3',1,-s.Rational(1,3)),('3bar',1,s.Rational(1,3)),('1',2,s.Rational(1,2)),('1',2,-s.Rational(1,2)),('1',1,s.Integer(0))]
    conj=lambda r:({'3':'3bar','3bar':'3','1':'1'}[r[0]],r[1],-r[2])
    dims={'3':3,'3bar':3,'1':1};assert sum(dims[a]*b for a,b,c in rows)==27
    full=Counter(rows*6+[conj(r) for r in rows]*6) # two species times three copies, plus conjugates
    assert all(full[r]==full[conj(r)] for r in full)
    original=Counter(rows*3)
    index={str(r):original[r]-original[conj(r)] for r in original if original[r]!=original[conj(r)]}
    return {'declared_model':'two(27,3) Weyl multiplets plus two conjugates','Weyl_components':324,
      'chiral_index_after_any_subgroup_restriction':0,'paired_multiplicities':{str(k):v for k,v in full.items()},
      'alternative_single_81_chiral_indices':index,
      'boundary':'A single(27,3) has SU3-family cubic gauge anomaly27 in fundamental units; dropping conjugates requires a new anomaly-canceling completion. Making the threefold factor global instead of gauged changes the theory. Gauge-preserving mass pairings cannot turn a zero representation index into three chiral generations.'}
if __name__=='__main__':
    out={'status':'PASS_EXACT_LIGHT_SINGLET_AND_DECLARED_CHIRAL_INDEX_OBSTRUCTION','little_group':stabilizer(),'standard_branching_audit':branching(),
      'mixing_boundary':'Singular values determine no CKM/PMNS matrix without named charged-current generators and left/right particle embeddings.',
      'prior_owners':['analysis/w33_e6_27_standard_model.py','analysis/w33_20261001_complete_gauge_scalar_fermion_loop.py','analysis/w33_20261001_exact_cartan_mirror_completion.py']}
    (ROOT/'data/w33_20261001_chiral_mass_assignment_audit.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
