#!/usr/bin/env python3
"""Does gauging the eight global family Goldstones kinematically close their decay?"""
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11270_full_condensate_moduli as Q

def vector_masses(gE,gF):
 p,_=Q.reference();T=Q.L.generators();gs=np.r_[np.full(78,gE),np.full(8,gF)]
 acts=np.einsum('aij,j->ai',T,p)*gs[:,None]
 M=2*np.real(acts.conj()@acts.T);w=np.linalg.eigvalsh(M);w[abs(w)<1e-10]=0.;assert min(w)>-1e-10
 return w,M

def width(m2,M2,r0=1.):
 x=m2/M2
 return 0. if x>=.25 or x<=0 else M2**1.5/(64*np.pi*r0*r0)*(1-4*x+12*x*x)*np.sqrt(1-4*x)

def payload():
 M2=648/7*.01**2;rows=[]
 for gF in (.001,.01,.05,.1,.3):
  ev,M=vector_masses(.5,gF);positive=ev[ev>1e-10];assert len(positive)==78
  rows.append({'gE_trace_normalized':.5,'gF_trace_normalized':gF,'vector_mass_squared':positive.tolist(),'open_vector_channels':int(sum(positive<M2/4)),'radial_vector_width':sum(width(v,M2) for v in positive)})
 # As gF->0, the8 additional family longitudinal channels approach8 Goldstones.
 ev,_=vector_masses(.5,1e-4);light=ev[ev>1e-12][:8]
 ratio=sum(width(v,M2) for v in light)/(M2**1.5/(8*np.pi));assert abs(ratio-1)<1e-5
 assert rows[0]['open_vector_channels']==8 and rows[-1]['open_vector_channels']==0
 return {'status':'PASS','result_scope':'PASS_GAUGED_FLAVOR_VECTOR_THRESHOLD_AND_EQUIVALENCE_THEOREM_CONTROL',
 'rows':rows,'small_gF_Goldstone_width_ratio':ratio,
 'model':'Gauge the former global family SU3, using actual E6xSU3 action on the condensate. Full86x86 vector mass matrix includes cross mixing; generic gF>0 has78 massive vectors and8 unbroken generators.',
 'vertex':'Every massive eigenvalue scales asr² along the radial orbit, so L has h V V vertex sqrt2*mV²/r0. For each identical real vector, Gamma=Mh³/(64pi r0²)(1-4x+12x²)sqrt(1-4x).',
 'insight':'The8 longitudinal channels reproduce the earlier8-Goldstone width asgF->0, then close when their physical vector thresholds exceedMh/2. Gauging does not automatically remove the decay: it moves it into longitudinal vectors.',
 'boundaries':['Kinematics of a separately declared gauged-family completion; no gauge coupling is predicted.','Family anomalies require additional chiral matter, such as the previously constructed10bar spectator. Its masses, gaugino mixing and new thresholds must be counted in the full completion.','The old global-family scalar/fermion pole is not transferred unchanged to this gauged theory.','The stored trace-normalized generators fix the coupling convention; numbers must not be read as measured gauge couplings.'],
 'prior_owners':['analysis/w33_pass11282_all_modulus_radial_self_energy.py','analysis/w33_pass11271_chiral_symmetric_yukawa_completion.py','analysis/w33_pass11270_full_condensate_moduli.py'],
 'primary_source':'https://arxiv.org/html/2505.07931v1'}
if __name__=='__main__':
 out=payload();(ROOT/'data/w33_pass11290_gauged_family_decay.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
