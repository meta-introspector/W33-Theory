"""Gildener-Weinberg scale in the isolated54 restriction; no full inventory claim."""
from pathlib import Path
import sys,json
import numpy as np
from scipy.integrate import solve_ivp
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11313_symmetric_adjoint_portal_closure as P

def payload():
 r=73/90;c=.1;g=.2;b=8/3;R=365/123;wave=120-2*b-8*R
 S=np.diag([9.]+[-1.]*9)/np.sqrt(90);K=np.zeros((10,10));H=P.hessians(S,K)
 lam=g*g*np.array([r*c,-c]);ms=np.linalg.eigvalsh(lam[0]*H[0,:54,:54]+lam[1]*H[1,:54,:54]);ms[abs(ms)<1e-12]=0
 mv=g*g*np.linalg.eigvalsh(P.vector_gram(S,K));mv[abs(mv)<1e-12]=0
 mf=g*g*R*np.repeat(np.diag(S)**2,2)
 def contribution(w,multiplicity,constant):
  w=w[w>1e-12];return multiplicity*np.sum(w*w*(np.log(w)-constant)),multiplicity*np.sum(w*w)
 terms=[contribution(ms,1,1.5),contribution(mv,3,5/6),contribution(mf,-2,1.5)]
 A=sum(x[0] for x in terms)/(64*np.pi**2);B=sum(x[1] for x in terms)/(64*np.pi**2)
 def flow(t,v):
  x,y=v;return [124*x*x+44.8*x*y+3.18*y*y-wave*x+18,24*x*y+25.4*y*y-wave*y+60-8*R*R]
 beta_flat=np.array(flow(0,[r*c,-c]))@np.array([1,r]);err=abs(B-beta_flat*g**4/(128*np.pi**2));assert err<1e-12 and B>0
 sol=solve_ivp(flow,[0,1.5],[r*c,-c],rtol=1e-11,atol=1e-12,dense_output=True)
 assert sol.success;track=sol.sol(np.linspace(0,1.5,151));bounded=np.minimum(track[0]+track[1]/10,track[0]+r*track[1]);assert min(bounded)>-1e-10
 v=np.exp(-.25-A/(2*B));m2=8*B*v*v
 assert sum(ms>0)==44 and sum(mv>0)==9
 return {'status':'PASS','result_scope':'PASS_CONDITIONAL_AF_CONNECTED_RADIATIVE_SCALE_IN_ISOLATED_54_MODEL','inputs':{'g':g,'lambda2_over_g2':-c,'mu_GW':1.,'y2_over_g2':R},'orbit':'SO10 -> SO9, S=diag(9,-1^9)/sqrt90; not the prior6+4 alignment','scalar_squared_masses':ms.tolist(),'vector_squared_masses':mv.tolist(),'Weyl_squared_masses':mf.tolist(),'A':float(A),'B':float(B),'vacuum_energy_over_v_fourth':float(-B/2),'mass_supertrace_beta_error':float(err),'v_over_mu_GW':float(v),'scalon_mass_squared_over_mu_GW_squared':float(m2),'scalon_mass_over_v':float(np.sqrt(8*B)),'UV_trajectory_terminal_ratios':track[:,-1].tolist(),'UV_trajectory_min_bounded_margin':float(min(bounded)),'boundary':['Scale is selected relative to an input RG crossing scale, not an absolute observed mass.','The induced vacuum energy is -B v4/2, not a prediction or cancellation of the observed cosmological constant.','Omitted45/frame/E6/family portals do not vanish under RG; this is an isolated restriction, not the full TOE vacuum.','Positive scalon curvature is leading one loop in MS-bar/Landau gauge; no higher-loop or full pole assertion.','44 positive scalar normals and9 vector modes; angular orbit differs from required6+4 frame alignment.'],'prior_owners':['analysis/w33_pass11303_scalar_source_and_uv_barrier.py','analysis/PASS20261001_FIVE_PHYSICAL_FRONTIERS.md'],'primary_sources':['https://doi.org/10.1103/PhysRevD.13.3333']}
if __name__=='__main__':
 d=payload();(ROOT/'data/w33_pass11318_radiative_scale_on_af_trajectory.json').write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['B'],d['v_over_mu_GW'],d['UV_trajectory_terminal_ratios'])
