#!/usr/bin/env python3
"""Canonical-kinetic D-flatness and resource extraction for the actual Cartan T ray.
Prior: Pass11260/61 signed embedding, balanced G26 vacuum and framed T port.
The kinetic norm here is a declared canonical 27x3 norm, not a mass prediction.
"""
from pathlib import Path
import sys,json,argparse
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11261_g26_cartan_coordinate_map as P
OUT=ROOT/'data/w33_20261001_cartan_dflat_resource_bridge.json'

def payload():
    K,_,_=P.slice_matrices();Q=s.Matrix(K).T;gram=Q.T*Q
    V=s.Matrix(27,3,K[2]);G=V.T*V;norm=s.trace(G)
    assert G==s.Matrix([[10,0,0],[0,13,5],[0,5,7]])
    mu=G-norm*s.eye(3)/3
    assert s.trace(mu*mu)==68 and G.det()==660
    rho=G/norm;purity=s.trace(rho*rho)
    assert purity==s.Rational(92,225)
    assert s.trace((rho-s.eye(3)/3)**2)==s.Rational(17,225)
    # Explicit SL(3,C) right-column balancing, not a compact gauge rotation.
    vn=np.array(V,float);gn=np.array(G,float);ev,U=np.linalg.eigh(gn)
    scale=float(G.det())**(1/3)
    h=U@np.diag(np.sqrt(scale/ev))@U.T
    balanced=vn@h;Gb=balanced.T@balanced
    assert np.max(abs(Gb-scale*np.eye(3)))<1e-12
    assert abs(np.linalg.det(h)-1)<1e-12
    # Reuse actual signed records to check covariance of the tangent operator.
    old=P._load(ROOT/'analysis/w33_pass11218_e8_cubic_cartan_kernel.py','dflat_old')
    parent=old.load_parent();records=parent.build_records()
    def op(v):
        M=np.zeros((81,81))
        for u,x,o,c in records:M[o,x]+=c*v[u]
        return M
    J=op(vn.ravel());Jb=op(balanced.ravel());U81=np.kron(np.eye(27),h.T)
    Ui=np.linalg.inv(U81)
    err=float(np.max(abs(Jb-Ui.T@J@Ui)))
    assert err<1e-12
    assert np.linalg.matrix_rank(J,tol=1e-9)==np.linalg.matrix_rank(Jb,tol=1e-9)==78
    # Optimal one-sided pure target extraction; the measurement is an input resource.
    z=np.exp(2j*np.pi/9);T=np.array([1,z,z.conjugate()])/np.sqrt(3)
    gi=np.linalg.inv(gn);den=float(np.real(T.conj()@gi@T))
    eta=vn@gi@T.conj()/np.sqrt(den)
    amp=eta.conj()@vn/np.sqrt(float(norm));prob=float(np.vdot(amp,amp).real)
    assert abs(np.vdot(eta,eta)-1)<1e-12
    assert np.max(abs(amp-np.sqrt(prob)*T))<1e-12
    return {
      'schema':'w33.20261001.cartan-dflat-resource-bridge.v1',
      'status':'PASS_EXACT_SL3_D_OBSTRUCTION_AND_EXPLICIT_BALANCING_EXTRACTION_MAPS',
      'assumption':'canonical positive kinetic norm Tr(V^dag V) on the 27x3 field',
      'exact':{'ambient_Cartan_Gram':gram.tolist(),'T_column_Gram':G.tolist(),
        'T_norm_squared':30,'Gram_determinant':660,
        'Gram_eigenvalues':['10','10-sqrt(34)','10+sqrt(34)'],
        'SL3_moment_norm_squared':68,
        'normalized_SL3_moment_norm_squared':'17/225',
        'normalized_color_purity':'92/225','color_Schmidt_rank':3,
        'grade_pair_Gram_relation':'the companion reconstructed grade-pair Gram equals 10 times this exact ambient Gram'},
      'SL3_balancing':{'map':'V -> V h, h=det(G)^(1/6) G^(-1/2)',
        'det_h':float(np.linalg.det(h)),'balanced_column_Gram':'660^(1/3) I3',
        'norm_squared_after':float(np.trace(Gb)),
        'norm_reduction':float(norm)-float(np.trace(Gb)),
        'Gram_error':float(np.max(abs(Gb-scale*np.eye(3)))),
        'cubic_operator_congruence_error':err,'cubic_operator_rank_before_after':[78,78],
        'scope':'noncompact SL3 balancing only; the E6 moment map and full vacuum equations are not solved'},
      'one_sided_T_extraction':{'formula':'p_max=1/(Tr(G) * T^dag G^(-1) T)',
        'exact_denominator':'T^dag G^(-1) T = 1/30 + 10/99 - (5/99)*cos(4*pi/9)',
        'probability':prob,'measurement_vector_real':eta.real.tolist(),
        'measurement_vector_imag':eta.imag.tolist(),
        'balanced_Dflat_probability':'1/3 for any target with an arbitrary chosen rank-one measurement',
        'resource_boundary':'the measurement contains target-dependent ninth-root phases; its implementation/Clifford cost is not supplied or assumed free'},
      'interpretation':[
        'The abstract Cartan qutrit coordinate and the literal external color qutrit are different objects.',
        'Interpreting normalized V as a 27 tensor 3 quantum amplitude gives the displayed mixed color marginal, not a pure T state.',
        'Canonical SU3 D-flatness forces V^dag V proportional to I3 and hence a maximally mixed color marginal under that quantum interpretation.',
        'A classical scalar field configuration is not automatically a quantum resource state.',
        'The global Cartan selector remains exact in its declared unitary metric; its compatibility with full-field kinetics and D-flatness is separate.',
        'No physical mass, compact E6 vacuum, magic distillation protocol or state-preparation probability is predicted.'
      ],
      'prior_owners':['analysis/w33_pass11260_exact_semisimple_cartan.py','analysis/w33_pass11261_g26_cartan_coordinate_map.py','analysis/w33_20261001_g26_global_polynomial_vacuum.py','analysis/w33_20261001_tmagic_orbit_factory_contract.py','analysis/PASS10962_HIDDEN_SECTOR_DFLAT_TIERS.md'],
      'checks':{'exact_SL3_obstruction':True,'exact_mixed_color_marginal':True,'SL3_balancing_control':True,'operator_covariance_control':True,'target_extraction_control':True}
    }

def main():
    a=argparse.ArgumentParser();a.add_argument('--check',action='store_true');args=a.parse_args()
    p=payload();txt=json.dumps(p,indent=2,sort_keys=True,default=str)+'\n'
    if args.check:assert OUT.read_text()==txt
    else:OUT.write_text(txt)
    print(p['status'])
if __name__=='__main__':main()
