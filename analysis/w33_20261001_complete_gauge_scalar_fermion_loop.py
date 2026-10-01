#!/usr/bin/env python3
"""Exact Cartan gauge-mass identity and full declared renormalizable one-loop test.

Landau gauge, MSbar; one canonical complex81 scalar, moment-square potential,
E6 x SU3 gauge fields, two81 Weyl multiplets and their two conjugates.
No claim of a Standard Model identification or a full higher-dimensional EFT loop.
"""
from pathlib import Path
import sys,json,argparse
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_20261001_exact_cartan_mirror_completion as C
from w33_20261001_dflat_cartan_mass_loop import loop_saddle
from w33_20261001_isovolume_dirac_gravity_response import compare_certificate
OUT=ROOT/'data/w33_20261001_complete_gauge_scalar_fermion_loop.json'

def rational_hermitian_basis():
    B=np.load(ROOT/'artifacts/e6_27rep_basis_export/E6_basis_78.npy').real.astype(int)
    unique={}
    for b in B:
        for kind,a in ((0,b+b.T),(1,b-b.T)):
            nz=np.flatnonzero(a)
            if not len(nz):continue
            gcd=int(np.gcd.reduce(abs(a.ravel()[nz])));a=a//gcd
            if a.ravel()[nz[0]]<0:a=-a
            unique[(kind,tuple(a.ravel()))]=s.SparseMatrix(a)*(s.I if kind else 1)
    H=list(unique.values());assert len(H)==78
    F=[]
    for i in range(3):
        for j in range(i+1,3):
            a=s.zeros(3);a[i,j]=a[j,i]=1;F.append(s.SparseMatrix(a))
            a=s.zeros(3);a[i,j]=-s.I;a[j,i]=s.I;F.append(s.SparseMatrix(a))
    F+=[s.SparseMatrix(s.diag(1,-1,0)),s.SparseMatrix(s.diag(1,1,-2))]
    gram=lambda hs:s.SparseMatrix([[s.trace(a*b) for b in hs] for a in hs])
    GE,GS=gram(H),gram(F)
    assert GE.rank()==78 and GS.rank()==8
    return H,F,s.SparseMatrix(GE.inv()),s.SparseMatrix(GS.inv())

def exact_gauge_identity():
    Q=C.plane();Ms,_=C.operators();H,F,GEi,GSi=rational_hermitian_basis()
    Vs=[s.SparseMatrix(Q[:,j].reshape(27,3)) for j in range(3)]
    def flatten(A):return s.Matrix(list(A))
    TE=[s.SparseMatrix.hstack(*[flatten(h*v) for h in H]) for v in Vs]
    TS=[s.SparseMatrix.hstack(*[flatten(v*f) for f in F]) for v in Vs]
    verified=0
    for j in range(3):
        for k in range(3):
            lhs=2*TE[j]*GEi*TE[k].conjugate().T+s.Rational(1,3)*TS[j]*GSi*TS[k].conjugate().T
            rhs=s.SparseMatrix(Ms[k].T@Ms[j])/3
            assert lhs==rhs;verified+=1
    return {'exact_quadratic_coefficient_matrices_checked':verified,
      'relative_gauge_couplings':'g_SU3=g_E6/sqrt(6) in trace-normalized Hermitian generators',
      'identity':'2 sum_E(T Phi)(T Phi)^dag + (1/3) sum_SU3(T Phi)(T Phi)^dag = M(Phi)^dag M(Phi)/3 on the exact Cartan plane',
      'method':'rational Hermitian bases, exact inverse trace Grams; no rounding or floating eigenspectrum proof',
      'gauge_spectrum':'78 nonzero squared gauge masses equal g_E6^2/3 times the cubic squared singular masses; eight gauge zeros',
      'scalar_spectrum':'V_D=(1/2)sum_a g_a^2 mu_a^2 gives the same78 massive real scalar eigenvalues,78 Goldstone zeros and6 Cartan zeros',
      'normalization_boundary':'This is a specified relative gauge normalization, not a prediction of measured gauge couplings.'}

def scale_checks():
    r,mu,A,B=s.symbols('r mu A B',positive=True)
    V=r**4*(A+B*s.log(r*r/(mu*mu)))
    stationary_log=-A/B-s.Rational(1,2)
    assert s.simplify(s.diff(V,r).subs(s.log(r*r/(mu*mu)),stationary_log))==0
    radial=s.simplify(s.diff(V,r,2).subs(s.log(r*r/(mu*mu)),stationary_log)/2)
    assert radial==4*B*r*r
    g,y=s.symbols('g y',real=True)
    coefficient=4*g**4/9-8*y**4
    Bloop=s.factor(8*coefficient/(64*s.pi**2))
    x,alpha,beta,kappa,N,theta=s.symbols('x alpha beta kappa N theta',positive=True)
    heat_volume=kappa*alpha**4*theta/(4*s.pi**2)
    planck_squared=kappa*alpha**2*theta/(24*s.pi**2)
    ratio=s.factor(heat_volume/planck_squared**2)
    assert ratio==144*s.pi**2/(kappa*theta)
    return {'CW_radial_potential':'r^4[A(q)+B log(r^2/mu^2)]',
      'B_in_this_model':str(Bloop),'stationary_radius':'r0=mu exp[-A/(2B)-1/4], stable radially iff B>0',
      'canonical_radial_mass_squared':str(radial),'vacuum_value':'-B r0^4/2 before an independent vacuum counterterm',
      'angular_boundary':'For B>0 the T saddle remains; radial transmutation does not supply a joint stable T vacuum.',
      'dilaton_model':'r=x*s, Lambda_heat=alpha*s, Lambda_EFT=beta*s, V=s^4 U(x), M_Pl^2=s^2 K(x)',
      'dilaton_stationarity':'dV/dr=0 => Uprime(x)=0; dV/ds=0 then requires U(x)=0. A nonzero stationary s remains flat; a classical homogeneous potential cannot fix the absolute scale.',
      'positive_single_heat_action':'Leading volume and two-derivative EH terms: M_Pl^2=kappa Lambda_heat^2 Theta_F/(24 pi^2); V_volume=kappa Lambda_heat^4 Theta_F/(4 pi^2)',
      'heat_volume_over_Planck_fourth':str(ratio),
      'scale_boundary':'Cutoff promotion supplies a named dynamical variable, not an absolute scale or tiny cosmological constant. An anomaly/transmutation scale or independent dimensionful datum remains required; loop corrections/counterterms can change the displayed single-heat-action relation.'}

def payload():
    identity=exact_gauge_identity();saddle=loop_saddle()
    return {'schema':'w33.20261001.complete-gauge-scalar-fermion-loop.v1',
      'status':'PASS_EXACT_GAUGE_MASS_IDENTITY_AND_DECLARED_FULL_LOOP_T_SADDLE',
      'field_content':{'scalar':'one complex Phi in(27,3) with canonical kinetic norm',
        'potential':'V_D=(1/2)sum g_a^2 mu_a^2; no polynomial T selector added',
        'fermions':'chi,psi in(27,3) and chi_tilde,psi_tilde in(27*,3*) with conjugate Yukawa tensors and equal |y|',
        'anomaly_check':'Representations paired with their conjugates cancel perturbative gauge anomalies; no Standard Model assignment. A Z2 distinguishing original and conjugate fermions excludes direct cross-sector bare bilinears in the declared model.',
        'higher_operators':'Dimension-five and gradient-completion coefficients set to zero for this renormalizable loop test; their loop corrections are not included.'},
      'gauge_mass_identity':identity,
      'full_one_loop':{'scheme':'MSbar Landau gauge on D-flat constant backgrounds; 3 massive-vector,1 real-scalar and8 paired-fermion degrees per matched cubic singular value',
        'formula':'64 pi^2 V1 = 3 sum mV^4(log(mV^2/mu^2)-5/6)+sum mS^4(log(mS^2/mu^2)-3/2)-8 sum mF^4(log(mF^2/mu^2)-3/2)',
        'angular_reduction':'64 pi^2 V1_angular=(4g^4/9-8|y|^4) F(q), F=sum sigma^4 log(sigma^2); all subtraction/scale terms are radial by the exact mirror design law',
        'result':'T is a saddle whenever the coefficient is nonzero; at its cancellation it is angularly flat at one loop, not an isolated minimum.',
        'cancellation':'g^4=18|y|^4','interval_sign_proof':saddle,
        'boundary':'Full one-loop field determinant for this declared model on the D-flat plane. No full higher-order EFT, two-loop calculation, off-plane global minimization or measured vacuum is claimed.'},
      'scale_dynamics':scale_checks(),
      'prior_owners':['analysis/w33_20261001_exact_cartan_mirror_completion.py','analysis/w33_20261001_dflat_cartan_mass_loop.py','analysis/BT1130_ricci_flat_seed_paradox_resolution.md'],
      'external_sources':['https://arxiv.org/abs/hep-ph/0111209','https://arxiv.org/abs/hep-th/0512169','https://arxiv.org/abs/hep-th/9606001']}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args();out=payload()
    if args.check:assert compare_certificate(out,json.loads(OUT.read_text()))
    else:OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(out['status'])
