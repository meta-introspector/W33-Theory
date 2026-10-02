#!/usr/bin/env python3
"""Leading complex pole blocks in an explicitly matched global-family quotient EFT."""
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11287_quotient_self_energy_matrix as Q
import w33_pass11290_gauged_family_decay as V

def inputs():
 d=json.loads((ROOT/'data/w33_pass11287_quotient_self_energy_matrix.json').read_text());H=np.array(d['Hessian']);C=np.array(d['cubic_tensor']);M=np.array(d['Weyl_mass_real'])+1j*np.array(d['Weyl_mass_imag']);Y=np.array(d['Yukawa_real'])+1j*np.array(d['Yukawa_imag']);return Q.cut_groups(H,C,M,Y)

def ward_matching(groups,sm,s0=-.0001,projector=None):
 P=np.diag((sm<1e-10).astype(float)) if projector is None else projector;SP=np.zeros((22,22));skipped=0.
 for g in groups:
  th=(np.sqrt(g['a'])+np.sqrt(g['b']))**2
  if th<1e-12:
   coupling=max(np.max(abs(g['R']@P)),np.max(abs(g['I']@P)));assert coupling<1e-12;skipped=max(skipped,float(coupling));continue
  SP+=Q.dispersion(g,0.,s0)@P
 C=-(SP+SP.T-P@SP);C=(C+C.T)/2
 assert np.max(abs((SP+C@P)))<1e-10
 return C,SP,P,skipped

def sigma(groups,t,C,s0=-.0001):
 return C+sum((Q.dispersion(g,t,s0)+1j*Q.rho(g,t) for g in groups),start=np.zeros((22,22),complex))

def payload():
 sm,fm,O,groups=inputs();C,SP,P,skip=ward_matching(groups,sm);rows=[]
 for t in sorted(set(round(float(x),12) for x in sm if x>1e-10)):
  ix=np.where(abs(sm-t)<1e-10)[0];S=sigma(groups,t,C);block=S[np.ix_(ix,ix)];shift=np.linalg.eigvals(block);poles=t-shift
  assert max(poles.imag)<1e-10
  rows.append({'tree_mass_squared':t,'multiplicity':len(ix),'pole_squared_real':poles.real.tolist(),'pole_squared_imag':poles.imag.tolist(),'leading_width_over_mass':(-poles.imag/t).tolist(),'nondegenerate_mixing_norm':float(np.linalg.norm(S[np.ix_(ix,np.where(abs(sm-t)>1e-10)[0])]))})
 # Global family remains global. Do not transplant gauged-family spectrum.
 ev,_=V.vector_masses(.5,0.);heavy=ev[ev>1e-10];assert len(heavy)==70
 minthreshold=4*min(heavy);assert minthreshold>max(sm)
 p,N,H,Ri,A=Q.geometry();acts=np.einsum('aij,j->ai',H,p);metric=(acts.conj()@acts.T).real
 w,U=np.linalg.eigh(metric);null=U[:,w<1e-10];unbroken=np.einsum('ab,aij->bij',null,H)
 neutrality=float(np.max(abs(np.einsum('aij,jk->aik',unbroken,N))));assert neutrality<1e-10
 assert sum(x['multiplicity'] for x in rows)==14
 return {'status':'PASS','result_scope':'PASS_LEADING_MATCHED_EFT_POLE_BLOCKS_WITH_LOCAL_TWO_POINT_WARD_CONDITIONS',
 'rows':rows,'Goldstone_poles':[0.]*8,'Ward_residual':float(np.max(abs(SP+C@P))),'massless_channel_Goldstone_projection_error':skip,'counterterm':C.tolist(),
 'matching':'At spacelike s0=-0.0001, normal-normal analytic mass and kinetic matching coefficients are set to zero as an input scheme. C=-(Sigma0 P+P Sigma0-P Sigma0 P) enforces (Sigma0+C)P=0. Only finite Goldstone-column limits are taken; divergent normal-normal massless scalar logs at s=0 are not evaluated.',
 'pole_method':'Inverse propagator s-Mtree²+Sigma(s). At one-loop order, diagonalize Sigma(Mtree²) inside each degenerate tree block; s_pole=Mtree²-eigenvalue. Off-block mixing and iterating momentum dependence change poles at higher loop order. Absorptive matrix PSD implies nonpositive pole imaginary parts.',
 'vector_control':{'gE_trace_normalized':.5,'global_family_gF':0.,'massive_E6_vectors':len(heavy),'smallest_two_vector_threshold_squared':float(minthreshold),'largest_scalar_tree_mass_squared':float(max(sm)),'unbroken_generator_slice_action_error':neutrality},
 'UV_boundary':'All massive E6 vector pair cuts are closed at the declared gE. The eight unbroken E6 generators act trivially on the quotient slice. Heavy vector analytic matching and gaugino/UV thresholds are not computed: changing the normal matching coefficients shifts poles at this same order.',
 'Ward_boundary':'These are local two-point Ward conditions with a stationary-vacuum matching prescription, not a derivation of an invariant counterterm functional or full tadpole equations. Therefore the poles are conditional matched-EFT outputs, not completed physical predictions of the UV theory.',
 'prior_owners':['analysis/w33_pass11287_quotient_self_energy_matrix.py','analysis/w33_pass11282_all_modulus_radial_self_energy.py','analysis/w33_pass11290_gauged_family_decay.py'],
 'primary_sources':['https://arxiv.org/abs/1609.06977','https://arxiv.org/abs/1910.02094']}
if __name__=='__main__':
 out=payload();(ROOT/'data/w33_pass11295_ward_matched_modulus_poles.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'],out['Ward_residual'])
