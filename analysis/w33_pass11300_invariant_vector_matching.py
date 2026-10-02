#!/usr/bin/env python3
"""An invariant vector-CW local counterterm and a conditional heavy-vector radial cut."""
from pathlib import Path
import sys,json
import numpy as np
import mpmath as mp
from scipy.integrate import quad
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11290_gauged_family_decay as V

def hermite(nodes,scale=.04,mu=.1):
 mp.mp.dps=70;n=2*len(nodes)+1;A=mp.matrix(n);b=mp.matrix(n,1);A[0,0]=1;b[0]=0
 g=lambda t:mp.mpf(str(scale))**2*t*t*(mp.log(mp.mpf(str(scale))*t/mp.mpf(str(mu))**2)-mp.mpf(5)/6)
 for i,x in enumerate(nodes):
  t=mp.mpf(str(x))/mp.mpf(str(scale))
  for k in range(n):A[2*i+1,k]=t**k;A[2*i+2,k]=k*t**(k-1) if k else 0
  b[2*i+1]=mp.diff(g,t);b[2*i+2]=mp.diff(g,t,2)
 c=mp.lu_solve(A,b);res=max(abs(a) for a in A*c-b)
 return [mp.mpf(0)]+[c[k]/(k+1) for k in range(n)],float(res)

def vector_dispersion(s,m2):
 threshold=4*m2;assert 0<s<threshold
 def spectral(t):
  x=m2/t;return t*t/(64*np.pi)*(1-4*x+12*x*x)*np.sqrt(max(0,1-4*x))
 return s**3/np.pi*quad(lambda t:spectral(t)/(t**3*(t-s)),threshold,np.inf,epsabs=1e-11,limit=200)[0]

def payload():
 w,_=V.vector_masses(.5,0.);positive=w[w>1e-10];nodes=[];multiplicities=[]
 for x in positive:
  if not nodes or abs(x-nodes[-1])>1e-9:nodes.append(float(x));multiplicities.append(1)
  else:multiplicities[-1]+=1
 c,error=hermite(nodes);scale=.04;mu=.1;f=lambda t:mp.mpf(str(scale))**2*t*t*(mp.log(mp.mpf(str(scale))*t/mp.mpf(str(mu))**2)-mp.mpf(5)/6)
 p=lambda t:sum(a*t**i for i,a in enumerate(c));deriv=0.
 for x in nodes:
  t=mp.mpf(str(x))/mp.mpf(str(scale));deriv=max(deriv,float(abs(mp.diff(f,t)-mp.diff(p,t))),float(abs(mp.diff(f,t,2)-mp.diff(p,t,2))))
 assert deriv<1e-40
 s=648/7*.01**2;shift=sum(n*vector_dispersion(s,x) for n,x in zip(multiplicities,nodes));assert shift>0
 remainder=3/(64*np.pi**2)*sum(n*float(f(mp.mpf(str(x))/scale)-p(mp.mpf(str(x))/scale)) for n,x in zip(multiplicities,nodes))
 return {'status':'PASS','result_scope':'PASS_INVARIANT_HEAVY_VECTOR_LOCAL_JET_MATCHING_AND_RADIAL_DISPERSION',
 'vector_groups':[{'mass_squared':x,'multiplicity':n} for x,n in zip(nodes,multiplicities)],'polynomial_coefficients':[str(a) for a in c],'polynomial_degree':len(c)-1,'Hermite_residual':error,'first_second_derivative_error':deriv,
 'invariant_counterterm':'Let X(q) be the actual78x78 E6 vector mass-squared matrix. V_CW=3Tr f(X/scale)/(64pi²), f(t)=scale²t²[log(scale*t/mu²)-5/6]. Add Vct=-3Tr P(X/scale)/(64pi²). Pprime(0)=0 and Pprime,Psecond match fprime,fsecond at every positive vacuum eigenvalue. A matrix spectral function trace is gauge/family invariant.',
 'local_jet_proof':'First derivatives depend on fprime at eigenvalues; second derivatives use fsecond and divided differences of fprime. Hermite interpolation matches both, hence cancels the full vector-CW tadpole/Hessian jet on the regular quotient slice. The vector matrix is a positive Gram matrix: its first variation restricted to the kernel is zero and masses emerging from the kernel are O(q²). Since f(x)=O(x² log x), those modes contribute no first or second field jet. Pprime(0)=0 matches their first-derivative term; no fsecond(0) value is assumed.',
 'vacuum_energy_remainder':remainder,'CC_boundary':'Matching tadpoles and curvature does not fix the residual constant or observed vacuum energy; an independent vacuum-energy matching coefficient remains.',
 'radial_vector_self_energy_real':float(shift),'leading_squared_pole_shift':float(-shift),'vector_added_width':0.,
 'dispersion_scheme':'Physical radial two-vector cut rho(t)=t²(1-4mV²/t+12mV4/t²)sqrt(1-4mV²/t)/(64pi r0²). Three subtractions at zero make its t² growth integrable. Local momentum coefficients through s² are declared zero; all pair cuts are closed at the scalar tree mass.',
 'boundaries':['Uses the earlier condensate EFT of11290; the enlarged Higgs/SO10 vector spectrum of11293 is not included.','An explicitly invariant higher-degree EFT counterterm functional in a declared Landau/MS vector-CW convention, not a renormalizable UV counterterm basis.','Scalar/Weyl tadpoles, full external vector/ghost/Goldstone momentum matrices and UV matching are not completed here.','The radial dispersion is a conditional matched EFT contribution; its subtraction freedom is not a physical mass prediction.'],
 'prior_owners':['analysis/w33_pass11295_ward_matched_modulus_poles.py','analysis/w33_pass11290_gauged_family_decay.py'],
 'primary_sources':['https://arxiv.org/abs/1910.02094','https://arxiv.org/abs/1609.06977']}
if __name__=='__main__':
 out=payload();(ROOT/'data/w33_pass11300_invariant_vector_matching.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'],out['leading_squared_pole_shift'])
