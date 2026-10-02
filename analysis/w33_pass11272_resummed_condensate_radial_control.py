#!/usr/bin/env python3
"""Hard1-loop radial shift and Ward-resummed global Goldstones in the11-modulus EFT.

Distinct from the paired nonsupersymmetric EFT. Heavy gauge multiplets and
threshold/Kähler matching are excluded; nonzero physical poles are not inferred.
"""
from pathlib import Path
import sys,json
import numpy as np
import sympy as s
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_20261001_completed_eft_one_loop as L

def matrices():
    q,_=L.tangent();fs=L.I.G.invariants();vs=L.I.G.VARS;at=dict(zip(vs,q));i=complex(fs[2].subs(at))
    g=np.array([complex(s.diff(fs[2],x).subs(at)) for x in vs])
    H=np.array([[complex(s.diff(fs[2],x,y).subs(at)) for y in vs] for x in vs])
    wh=(3/14)*(4*np.outer(g,g)/(9*i*i)-H/(3*i))
    gs,_=L.I.C.gradient_data();b=np.array([v.conj()/np.linalg.norm(v) for v in gs]).T
    dirs=np.column_stack([b[:,0],1j*b[:,0],b[:,1],1j*b[:,1],q,1j*q])
    soft=-9*np.real(dirs.T@wh@dirs)
    total=np.array([[float(s.Rational(v)) for v in row] for row in json.loads((ROOT/'data/w33_20261001_condensate_cartan_potential.json').read_text())['exact_slice_Hessian_matrix_over_m_squared']])
    # C.gradient_data phases agree with the invariant-gradient slice convention.
    assert np.max(abs(np.linalg.eigvalsh(soft[:4,:4])-np.array([-162/7]*2+[162/7]*2)))<1e-9
    return total-soft,soft

def payload(m=.01):
    F,A=matrices();mu=m
    ferm=np.array([81/49]*8+[324/49]*2+[81.])*m*m
    def masses(r):
        six=np.linalg.eigvalsh(m*m*(F*r**-16+A*r**-8))
        partners=np.full(8,m*m*((729/49)*r**-16-(81/7)*r**-8))
        return np.r_[six,partners],ferm*r**-16
    def f(x):
        assert min(x)>0
        return np.sum(x*x*(np.log(x/(mu*mu))-1.5))
    def hard(r):
        bos,fer=masses(r);return float((f(bos)-2*f(fer))/(64*np.pi**2))
    def tree(r):return m*m*(81/(49*r**14)-27/(7*r**6))
    def d1(fun,r,h=1e-5):return (fun(r-2*h)-8*fun(r-h)+8*fun(r+h)-fun(r+2*h))/(12*h)
    def dtree(r):return m*m*(-1134/(49*r**15)+162/(7*r**7))
    total=lambda r:tree(r)+hard(r)
    r=brentq(lambda x:dtree(x)+d1(hard,x),.99,1.01,xtol=1e-13)
    G=m*m*(81/7)*(r**-8-r**-16);Delta=d1(hard,r)/(2*r)
    assert abs(G+Delta)<1e-12
    curv=(total(r+.0001)+total(r-.0001)-2*total(r))/(2*.0001**2)
    assert curv>0
    # Independent Ward identity: tree derivative determines all8 compact flavor
    # Goldstone masses, whose soft self-energy is the hard tadpole divided by2r.
    assert abs(G-dtree(r)/(2*r))<1e-14
    return {'status':'PASS','result_scope':'PASS_HARD_MODULUS_LOOP_RADIAL_SHIFT_AND_WARD_RESUMMED_GOLDSTONES',
      'parameters':{'m':m,'mu':mu,'tree_radius':1.,'Lambda_ninth':'3m/14'},'hard_scalar_modes':14,'Weyl_modes':11,'soft_global_Goldstones':8,
      'shifted_radius':r,'radius_shift':r-1,'hard_one_loop_at_tree_radius':hard(1),'hard_one_loop_at_shifted_radius':hard(r),
      'tree_Goldstone_mass_squared_at_shifted_radius':G,'computed_hard_zero_momentum_shift':Delta,'resummed_Goldstone_mass_squared':G+Delta,
      'radial_static_curvature':curv,'tree_radial_static_curvature':(648/7)*m*m,
      'scalar_homogeneity':'M_scalar²(r)=m²[F*r^-16+A*r^-8]; F is obtained by subtracting the explicit AMSB holomorphic Hessian from the previously certified6-direction Hessian. Eight noncompact partners have(729/49)r^-16-(81/7)r^-8. Eight compact Goldstones have(81/7)(r^-8-r^-16).',
      'hard_determinant':'[sum14 mS^4(log(mS²/mu²)-3/2)-2sum11 mF^4(log(mF²/mu²)-3/2)]/(64pi²); soft8 Goldstones treated by the Ward-resummed kernel.',
      'Ward_resummation':'Delta_hard=V1hard_prime/(2r). At the hard-shifted stationary point, Gtree+Delta_hard=0. The8 f(Gbar) soft terms and their first derivatives vanish there, avoiding insertion of an arbitrary positive IR regulator. Double-counting subtraction is needed when adding explicit higher-loop terms.',
      'protected_physical_poles':'If the globalSU3-breaking vacuum persists, its8 Goldstone poles are exactly massless by the continuous symmetry. This protection is distinct from the unphysical gauge Goldstones of the paired EFT.',
      'nonzero_physical_poles_computed':False,
      'boundary':['One-loop canonical11-modulus effective theory with specifiedm/mu and tree radius, not the full E6 ultraviolet theory.','Heavy gauge multiplet thresholds, Kähler running/matching, higher-loop terms and gauge dynamics are omitted.','Radial curvature is a zero-momentum static quantity; nonzero pole masses require momentum-dependent self-energies.','The paired EFT hardGoldstone1PI problem in the companion audit is not solved by substituting this different theory.'],
      'prior_owners':['analysis/w33_20261001_condensate_cartan_potential.py','analysis/w33_pass11270_full_condensate_moduli.py','analysis/w33_pass11272_goldstone_ir_and_pole_audit.py'],'primary_sources':['https://arxiv.org/html/2505.07931v1','https://arxiv.org/abs/1406.2355']}
if __name__=='__main__':
    out=payload();(ROOT/'data/w33_pass11272_resummed_condensate_radial_control.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
