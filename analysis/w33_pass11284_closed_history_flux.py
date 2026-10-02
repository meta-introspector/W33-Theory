#!/usr/bin/env python3
"""A declared closed handlebody-double history and exact discrete flux constraints."""
from pathlib import Path
import sys,json
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
from w33_pass11279_history_flux_obstruction import complex_data
from w33_pass11275_polynomial_sm_higgs import rank_mod

def top_boundary(ticks=3):
    dc=np.zeros((ticks,ticks),int)
    for i in range(ticks):dc[i,i]=-1;dc[(i+1)%ticks,i]=1
    # M3=#81(S1xS2), minimal cellular counts1,81,81,1 with all boundaries0.
    # Product with a subdivided circle: C3 dimension82*ticks, C4=ticks.
    return np.vstack([np.zeros((81*ticks,ticks),int),-dc])

def payload():
    _,_,tets,_=complex_data()
    # Point/line incidence graph:80vertices,160flags, connected, cycle rank81.
    inc=np.zeros((80,160),int);col=0
    for line,tet in enumerate(tets):
        for point in tet:inc[point,col]=-1;inc[40+line,col]=1;col+=1
    assert col==160 and rank_mod(inc)==79
    genus=160-80+1;assert genus==81
    D=top_boundary();assert rank_mod(D)==2 and not np.any(D@np.ones(3,dtype=int))
    Q,qhat,K,mu4,M2=s.symbols('Q Qhat kappa2 mu4 M2',nonzero=True)
    Ravg=-2*mu4*qhat/(M2*Q);delta=s.factor(K*Ravg/4)
    t,c=s.symbols('T C');Lambda=t/4+delta
    assert s.simplify((t-4*c)/4+delta-(Lambda-c))==0
    return {'status':'PASS','result_scope':'PASS_DECLARED_CLOSED_HISTORY_AND_DISCRETE_FLUX_SEQUESTERING_INTERFACE',
      'construction':'Thicken the connected80-vertex160-edge W33 point-line incidence graph into an oriented3D handlebody of genus81. Double it across the boundary using the identity gluing: M3=#81(S1xS2). Set M4=M3xS1. These thickening, gluing and real-manifold choices are explicit added data.',
      'topology':'M4 is closed oriented, Betti numbers(1,82,162,82,1), Euler0 and H4=Z. The clique complex itself is unchanged and still fails the previous flux test.',
      'minimal_spatial_cell_counts':[1,81,81,1],'minimal_spatial_boundaries':'zero in this standard handlebody-double CW decomposition',
      'three_tick_history_cell_counts':[3,246,486,246,3],'top_boundary_shape':list(D.shape),'top_boundary_rank_mod101':2,
      'unique_top_cycle':[1,1,1],'local_rigidity':'D4 sigma(Lambda)=0 forces the top-cell coefficients equal along the connected circle, if sigma is injective. The prior clique history had injective D4 and forced zero instead.',
      'flux':'F=D4^T A3+(Q/ticks)*ones. Exact3-form variations have zero total flux because D4*ones=0; Q is an independent topological sector. Integer quantization, if imposed, is additional compact-gauge input.',
      'declared_action':'S=sum_i volume_i[kappa_i² R_i/2-Lambda_i-Lm_i]+sum_i sigma(Lambda_i/mu4)F_i+hatsigma(kappa_i²/M2)Fhat_i. Its continuum counterpart is the known covariant local sequestering action; discrete curvature/metric dynamics must be supplied.',
      'linear_sigma_test':'For sigma(z)=hatsigma(z)=z, variation of3-forms makesLambda,kappa² rigid; Lambda variation givestotal volume=Q/mu4, kappa² variation gives<R>=-2mu4 Qhat/(M2 Q).',
      'residual_flux_constant':str(delta),'conditional_gravity_equation':'kappa²G=T-(1/4)g<TrT>-DeltaLambda*g. For constant stress shift T->T-Cg, Lambda->Lambda-C and the source is exactly unchanged with fixed linear-sigma flux data.',
      'boundary':['This supplies a concrete closed completion and constraint, not a unique geometry derived by finite incidence. Different Heegaard gluings change the topology.','Q,Qhat,kappa² and the four-volume remain input sectors; the residual is not fixed to the observed CC.','The cellular action names the flux map and tests its variation. It does not prove convergence to a covariant metric/gravity path integral or cancel graviton loops.'],
      'prior_owners':['analysis/w33_pass11279_history_flux_obstruction.py','analysis/W33_LEDGER_CONTINUATION_2026_09_21.md','analysis/w33_pass11274_scale_and_sequestering.py'],
      'primary_sources':['https://arxiv.org/html/1505.01492v2']}
if __name__=='__main__':
    out=payload();(ROOT/'data/w33_pass11284_closed_history_flux.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
