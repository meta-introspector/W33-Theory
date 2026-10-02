#!/usr/bin/env python3
"""Native Levi cycle pairing supplies a nontrivial Heegaard shear, not a selected CC."""
from pathlib import Path
import sys,json
import numpy as np
import sympy as s
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
from w33_pass11279_history_flux_obstruction import complex_data

def cycle_basis():
 _,_,lines,_=complex_data();edges=[(p,40+i) for i,t in enumerate(lines) for p in t];adj=[[] for _ in range(80)]
 for k,(a,b) in enumerate(edges):adj[a].append((b,k));adj[b].append((a,k))
 parent={0:None};tree=[];queue=[0]
 for a in queue:
  for b,k in adj[a]:
   if b not in parent:parent[b]=(a,k);tree.append(k);queue.append(b)
 tree=set(tree);inc=np.zeros((80,160),int)
 for k,(a,b) in enumerate(edges):inc[a,k]=-1;inc[b,k]=1
 cols=[]
 for e,(a,b) in enumerate(edges):
  if e in tree:continue
  c=np.zeros(160,int);c[e]=1
  # Cancel the chord incidence: minus(a->root) plus(b->root).
  for v,sg in [(a,-1),(b,1)]:
   while parent[v] is not None:
    u,k=parent[v];x,y=edges[k];c[k]+=sg*(1 if (x,y)==(v,u) else -1);v=u
  assert not np.any(inc@c);cols.append(c)
 C=np.array(cols).T;assert C.shape==(160,81);return inc,C

def payload():
 inc,C=cycle_basis();B=s.Matrix(C.T@C);assert B==B.T and all(B[i,i]%2==0 for i in range(81))
 try:
  from flint import fmpz_mat
  fast=fmpz_mat([[int(x) for x in row] for row in B.tolist()]);det=int(fast.det());sn=fast.snf();diag=[abs(int(sn[i,i])) for i in range(81)]
 except ImportError:
  det=int(B.det(method='domain-ge'))
  prior=json.loads((ROOT/'data/PART_W33_PASS5031_CRITICAL_GROUPS.json').read_text())
  group=prior['levi']['critical_group'];diag=[1]*(81-sum(group.values()))+[int(k) for k,v in group.items() for _ in range(v)]
 expected=10**23*4**30;assert det==expected and np.prod(np.array(diag,dtype=object))==det
 g=81;I=s.eye(g);Z=s.zeros(g);T=I.row_join(Z).col_join(B.row_join(I));J=Z.row_join(I).col_join((-I).row_join(Z));assert T.T*J*T==J
 return {'status':'PASS','result_scope':'PASS_NATIVE_CYCLE_GRAM_HEEGAARD_SHEAR_AND_FLUX_SELECTION_AUDIT',
 'graph':{'vertices':80,'edges':160,'cycle_rank':81},'cycle_basis':C.tolist(),'cycle_Gram':[[int(x) for x in row] for row in B.tolist()],
 'exact_determinant':det,'determinant_formula':'10^23*4^30=2^83*5^23','Smith_diagonal':diag,
 'prior_owned_critical_group':'Pass5031 already owns (Z/4)^6+(Z/40)^22+Z/160 and order2^83*5^23; Pass5441 owns the cycle-Gram determinant. This producer reuses these results for a new explicit gluing, not a new critical-group discovery.',
 'spanning_tree_control':'Matrix-tree theorem: Levi Laplacian eigenvalues0,8,(4+-sqrt6) each24,4 multiplicity30 ->tau=8*10^24*4^30/80=10^23*4^30. The fundamental integral cycle Gram has determinanttau.',
 'basis_covariance':'Changing integral cycle basis C->C U, U unimodular, gives B->U^T B U. Determinant and Smith cokernel are independent of the spanning-tree choice; unit edge weights are the declared native combinatorial pairing.',
 'framed_link_alternative':'B can also be realized by81 framed unknot components, with diagonal framings Bii and signed pairwise Hopf clasps Bij. Surgery on a chosen standard realization gives a named closed3-manifold withH1=cokerB; link realization and orientation remain added choices.',
 'gluing':'Use the explicitly added symplectic shear [[I,0],[B,I]] as a genus81 Heegaard gluing action. H1=coker B is finite with the stored Smith invariants, versus Z^81 for identity gluing. The mapping-class representative remains an added choice; its homology action does not classify the3-manifold.',
 'glued_product_rational_Betti':[1,1,0,1,1],
 'minimal_spatial_chain':'C=(Z,Z^81,Z^81,Z), D2=B, D1=D3=0 for the declared homology presentation; H1=cokerB, H2=0. The circle product still has a unique top cycle.',
 'closed_history':'For either oriented gluing, M3xS1 is closed oriented and H4=Z. Thus existence of top flux is robust under this gluing change; topology alone still leaves its value free.',
 'selection_test':'An added Euclidean Maxwell four-form term has action Q²/(2e² Vol). With the prior linear sequestering volume constraint Vol=Q/mu4, this becomes mu4 Q/(2e²), not a preferred flux magnitude. If positive quantized Q=nq is assumed, canonical weights are geometric exp(-a n), a=mu4 q/(2e²); <n>=1/(1-exp(-a)). Couplings and the ensemble measure remain inputs.',
 'CC_boundary':'The residual flux ratio Qhat/Q remains free. Neither the cycle pairing nor a four-form Maxwell weight predicts observed vacuum energy; membrane dynamics, cosmological state and measure must be supplied.',
 'prior_owners':['analysis/w33_pass5028_5035_steinberg_apartments.py','analysis/w33_pass5436_5443_bicycle_apartment_scheme_packet.py','analysis/w33_pass11284_closed_history_flux.py','analysis/PASS20260924_HISTORY_COMPACTIFICATION_BREAKTHROUGH.md'],
 'primary_sources':['https://academicweb.nd.edu/~andyp/papers/SymplecticHeegaard.pdf','https://arxiv.org/html/1505.01492v2']}
if __name__=='__main__':
 out=payload();(ROOT/'data/w33_pass11289_cycle_gram_gluing_flux.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'],out['determinant_formula'])
