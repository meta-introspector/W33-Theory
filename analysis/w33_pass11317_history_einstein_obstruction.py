"""Exact topological obstruction for smooth closed Riemannian Einstein histories."""
from pathlib import Path
import json
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
def payload():
 g=81;betti=[1,g+1,2*g,g+1,1];chi=sum((-1)**i*b for i,b in enumerate(betti));Q=s.zeros(2*g)
 for i in range(g):Q[i,g+i]=Q[g+i,i]=1
 assert Q*Q==s.eye(2*g) and chi==0
 return {'status':'PASS','result_scope':'PASS_NO_SMOOTH_CLOSED_RIEMANNIAN_EINSTEIN_METRIC_ON_DECLARED_HISTORY','history':'(#81(S1 x S2)) x S1','betti':betti,'Euler_characteristic':chi,'intersection_form':'81 hyperbolic planes','signature':0,'fundamental_group':'F81 x Z','free_group_reduced_word_count_length3':2*g*(2*g-1)**2,'proof':['For an Einstein4-metric the Chern-Gauss-Bonnet identity is8pi²chi=integral(|W|²+R²/24).','chi=0 forces W=0 and R=0, hence the metric is flat.','A closed flat4-manifold is finitely covered byT4 and has virtually abelian fundamental group.','A flat4-manifold also has b1<=4, immediately contradicting b1=82.','F81 x Z has exponential reduced-word growth and cannot be virtually abelian. Thus no such metric exists for any Einstein constant.'],'membrane_boundary':['A degree-one collapse toS4 transfers the top-dimensional flux pairing, not a nonsingular round Einstein metric.','The alternative11289 nonidentity Heegaard gluing changes the rational Betti numbers; this proof is for the declared identity-gluing81-handle history.','This excludes the smooth closed source-free Einstein saddle on this particular history. It does not exclude piecewise Einstein saddles with thin membranes, anisotropic matter, Lorentzian solutions or different topology.','No full wall saddle or its determinant has been constructed. The prior spherical two-cap reduced radial result remains scoped to its supplied geometry.'],'prior_owners':['analysis/w33_pass11307_quantized_membrane_radial_mode.py','analysis/PASS11298_11302_RUNNING_FLAVOR_GRAVITY_AND_JUNCTIONS.md'],'primary_sources':['https://www.math.stonybrook.edu/~claude/cambridge.pdf']}
if __name__=='__main__':
 d=payload();(ROOT/'data/w33_pass11317_history_einstein_obstruction.json').write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['betti'])
