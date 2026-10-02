"""Constructed local throat and independently derived Einstein tensor."""
from pathlib import Path
import json
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
def payload():
 t,x,th,ph,a=s.symbols('t x theta phi a',real=True);coords=[t,x,th,ph];r2=x*x+a*a;g=s.diag(-1,1,r2,r2*s.sin(th)**2);gi=g.inv();n=4
 G=[[[s.simplify(sum(gi[i,l]*(s.diff(g[l,k],coords[j])+s.diff(g[l,j],coords[k])-s.diff(g[j,k],coords[l])) for l in range(n))/2) for k in range(n)] for j in range(n)] for i in range(n)]
 Ric=s.zeros(n)
 for i in range(n):
  for j in range(n):
   Ric[i,j]=s.simplify(sum(s.diff(G[k][i][j],coords[k])-s.diff(G[k][i][k],coords[j])+sum(G[k][i][j]*G[l][k][l]-G[l][i][k]*G[k][j][l] for l in range(n)) for k in range(n)))
 R=s.simplify(s.trace(gi*Ric));Ein=s.simplify(Ric-g*R/2);rho=Ein[0,0];pr=Ein[1,1];pt=s.simplify(Ein[2,2]/r2);null=s.simplify(rho+pr)
 assert s.simplify(R+2*a*a/r2**2)==0 and s.simplify(null+2*a*a/r2**2)==0
 return {'status':'PASS','result_scope':'PASS_LOCAL_SMOOTH_NECK_REQUIRES_RADIAL_NULL_ENERGY_VIOLATION','metric':'-dt²+dx²+(x²+a²)dOmega2², a>0','scalar_curvature':str(R),'stress_times_8piG':{'rho':str(rho),'radial_pressure':str(pr),'angular_pressure':str(pt),'radial_null_contraction':str(null)},'scope':['A nonsingular local Lorentzian Ellis throat, not the global81-handle metric or a Euclidean membrane instanton.','The required negative radial null contraction cannot be supplied by a cosmological constant or positive-energy canonical scalar/Maxwell matter alone.','Quantum stress, higher curvature, wall topology and non-static necks are untested alternatives; no complete gravity model is asserted.'],'prior_owners':['analysis/w33_pass11317_history_einstein_obstruction.py'],'primary_sources':['https://doi.org/10.1063/1.1666161']}
if __name__=='__main__':
 d=payload();(ROOT/'data/w33_pass11319_smooth_neck_stress.json').write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['stress_times_8piG'])
