#!/usr/bin/env python3
"""Concrete native Levi graph matter-bracket audit before a gravity completion."""
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
from w33_pass11289_cycle_gram_gluing_flux import cycle_basis

def bracket(phi,p,N,M,edges):
 def gradient(f):
  out=f*(phi+phi**3)
  for i,j in edges:
   a=.5*(f[i]+f[j])*(phi[i]-phi[j]);out[i]+=a;out[j]-=a
  return out
 lhs=float(gradient(N)@(M*p)-(N*p)@gradient(M))
 rhs=sum(.5*(N[i]*M[j]-M[i]*N[j])*(p[i]+p[j])*(phi[j]-phi[i]) for i,j in edges)
 return lhs,float(rhs)

def cycle_error(n):
 h=2*np.pi/n;x=np.arange(n)*h;phi=np.sin(x)+.2*np.cos(2*x);p=1+.2*np.cos(x);N=1+.3*np.sin(x);M=1+.2*np.sin(x)+.1*np.cos(2*x)
 # Correct canonical cell momentum is h*p. Link gradient weight is1/h.
 edges=[(i,(i+1)%n) for i in range(n)]
 b=sum(.5/h*(N[i]*M[j]-M[i]*N[j])*(p[i]+p[j])*(phi[j]-phi[i]) for i,j in edges)
 xx=np.arange(4096)*2*np.pi/4096;expected=np.mean((1+.2*np.cos(xx))*(np.cos(xx)-.4*np.sin(2*xx))*((1+.3*np.sin(xx))*(.2*np.cos(xx)-.2*np.sin(2*xx))-(1+.2*np.sin(xx)+.1*np.cos(2*xx))*.3*np.cos(xx)))*2*np.pi
 return abs(float(b-expected))

def payload():
 inc,_=cycle_basis();edges=[tuple(np.where(inc[:,j])[0]) for j in range(160)];rng=np.random.default_rng(11301);phi,p,N,M=rng.normal(size=(4,80));a,b=bracket(phi,p,N,M,edges);assert abs(a-b)<1e-10
 errs=[cycle_error(n) for n in (16,32,64,128)];assert all(errs[i+1]<errs[i]/3 for i in range(3))
 return {'status':'PASS','result_scope':'PASS_NATIVE_GRAPH_ADM_MATTER_BRACKET_OBSTRUCTION_AND_REFINEMENT_CONTROL',
 'native_graph':{'vertices':80,'edges':160},'Hamiltonian':'H[N]=sum_i N_i(p_i²/2+V(phi_i))+sum_edges(N_i+N_j)(phi_i-phi_j)²/4. Canonical{phi_i,p_j}=deltaij; V=phi²/2+phi4/4.',
 'exact_bracket':'{H[N],H[M]}=sum_edges(N_i M_j-M_i N_j)D_ij, D_ij=(p_i+p_j)(phi_j-phi_i)/2. Potential cancels exactly. This creates edge-current generators absent from the site Hamiltonian ansatz.',
 'closure_boundary':'The site Hamiltonians alone are not closed under field-independent lapse brackets: their quadratic parts contain pp and phiphi, while the nonzero bracket contains mixed pphi. Adding edge constraints and gravitational variables is a further construction, not an inherited continuum ADM proof.',
 'finite_derivation_proof':'For a pointwise finite vertex algebra, e_i²=e_i implies D(e_i)=2e_i D(e_i). Outside componenti it forces zero; componenti gives x=2x. Hence every exact Leibniz derivation is zero. A nonzero graph difference cannot be the undeformed continuum derivative on all vertex fields.',
 'bracket_control':a,'edge_control':b,'error':abs(a-b),'cycle_refinement_errors':errs,
 'scope':'A native W33-Levi regulator audit and a periodic-cycle O(h²) continuum control; not a completed discretized gravity theory. The fixed Levi graph has no supplied real3D refinement or metric embedding.',
 'prior_owners':['analysis/w33_pass11296_covariant_matter_constraints.py','analysis/w33_pass11289_cycle_gram_gluing_flux.py'],
 'primary_sources':['https://arxiv.org/abs/0810.2360']}
if __name__=='__main__':
 out=payload();(ROOT/'data/w33_pass11301_graph_adm_regulator.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
