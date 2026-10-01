#!/usr/bin/env python3
"""A separate chiral anomaly repair and E6-only enhancement/degree18 bridge.

Classical anomaly coefficients and published near-SUSY E6 dynamics are cited.
No replacement of the vectorlike EFT, no inferred non-SUSY condensate minimum.
"""
from pathlib import Path
import sys,json,itertools
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_20261001_exact_cartan_mirror_completion as C

def symmetric_cube_weights():
    weights=[(1,1),(-1,1),(0,-2)]
    sums=[tuple(sum(weights[i][a] for i in ix) for a in range(2)) for ix in itertools.combinations_with_replacement(range(3),3)]
    assert len(sums)==10
    tr2=lambda rows:sum(b*b for a,b in rows)
    tr3=lambda rows:sum(b**3 for a,b in rows)
    assert s.Rational(tr2(sums),tr2(weights))==15
    assert s.Rational(tr3(sums),tr3(weights))==27
    return {'fundamental_weights':weights,'decuplet_weights':sums,'T10':'15/2','A10':27,
      'chiral_completion':'one Weyl(27,3) plus one Weyl(1,10bar); one complex scalar(27,3)',
      'perturbative_family_anomaly':'27*A3+A10bar=27-27=0','E6_anomaly':'E6 has no adjoint cubic Casimir, so the chiral27 is perturbatively anomaly free.',
      'one_loop_b0':{'E6':'35','SU3_family':'-15/2'},'total_Weyl_components':91,
      'mass_interface_boundary':'The old antisymmetric cubic matrix requires two distinct Weyl species; a single species bilinear would vanish by spinor symmetry. The single81 chiral candidate therefore does not inherit that matrix. Two species plus two10bar spectators give six27 generations and bE6=29,bSU3=-43/2, not three observed families.',
      'boundary':'A distinct anomaly-free chiral candidate, not the earlier paired model. The spectator is Standard-Model neutral under a chosen E6 embedding, but its mass and family gauge breaking require new operators. This candidate does not repair the regular T vacuum\'s insufficient unbroken gauge algebra.'}

def enhancement():
    B=np.load(ROOT/'artifacts/e6_27rep_basis_export/E6_basis_78.npy').real.astype(int);Q=C.plane();rows=[]
    for q,expected in [((1,2,3),8),((0,1,2),8),((1,1,2),14),((1,0,0),28),((1,1,1),28)]:
        V=(Q@np.array(q)).reshape(27,3);A=s.SparseMatrix(np.array([(b@V).ravel() for b in B]).T)
        unbroken=78-A.rank();assert unbroken==expected
        u=[str(x.subs(dict(zip(C.G.VARS,q)))) for x in C.G.invariants()]
        rows.append({'q':q,'unbroken_E6_dimension':unbroken,'u6_u12_u18':u})
    # Entire generic SIC wall, over Q(z), not a single numerical drawing.
    z=s.Symbol('z');V=s.Matrix(27,3,list(s.Matrix(Q)*s.Matrix([1,1,z])))
    A=s.SparseMatrix.hstack(*[s.Matrix(list(s.Matrix(b)*V)) for b in B]);null=A.nullspace()
    assert len(null)==14
    assert all(A*v==s.zeros(81,1) for v in null)
    return {'exact_points':rows,'generic_SIC_wall':'q=(1,1,z) over Q(z); exact E6 stabilizer dimension14',
      'MUB_only_control':'q=(0,1,2) has a coordinate MUB mirror but u18!=0 and E6 stabilizer remains8',
      'denominator_identification':'Conditional on the cited condensate denominator being a degree18 E6xSL3 invariant that vanishes on enhanced E6 loci: the9 SIC walls form one G26 orbit with order2 reflections. Invariance forces even vanishing order on each wall. Degree18 then forces a scalar multiple of their squared product u18, hence of the compiled global I18.',
      'nonperturbative_interface':'W=k*(Lambda^27/I18)^(1/3) on a chosen nonzero-I18 branch; degree3, homogeneous degree-6 in fields. k and branch normalization must be matched to the published tensor conventions.',
      'regime':'Published theory: E6 gauge, family SU3 global, supersymmetric chiral multiplets plus small AMSB. Weakly gauging family SU3 with an anomaly spectator changes the dynamics and is not proved equivalent.',
      'T_class_connection':'I6=I12=0,I18!=0 is the class with zero degree6 singlet and zero degree12 Aronhold invariant, but nonzero degree18 Aronhold/Hessian invariant. It corresponds algebraically to the published fully broken global-family class after normalization, not to an observed Standard Model vacuum.',
      'boundary':'No normalization of the published effective superpotential, Kähler metric, AMSB soft terms, quantum minimum, mass ratios or cosmological constant is asserted here.'}
if __name__=='__main__':
    out={'status':'PASS_CHIRAL_DECUplet_ANOMALY_AND_EXACT_E6_SIC_ENHANCEMENT','anomaly_repair':symmetric_cube_weights(),'condensate_bridge':enhancement(),
      'prior_owners':['analysis/w33_anomaly_cancellation.py','analysis/w33_e6_27_standard_model.py','analysis/w33_pass11269_g26_qutrit_dictionary.py','analysis/w33_20261001_degree18_phase_completion.py'],
      'primary_source':'https://arxiv.org/html/2505.07931v1'}
    (ROOT/'data/w33_20261001_chiral_decuplet_and_condensate_bridge.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
