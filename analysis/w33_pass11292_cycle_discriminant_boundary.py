#!/usr/bin/env python3
"""A prior-owned cycle lattice defines an added Abelian CS model with a boundary obstruction."""
from pathlib import Path
import json,math
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
def payload():
 d=json.loads((ROOT/'data/w33_pass11289_cycle_gram_gluing_flux.json').read_text());B=s.Matrix(d['cycle_Gram']);n=d['exact_determinant'];assert n==2**83*5**23 and math.isqrt(n)**2!=n
 assert B==B.T and all(B[i,i]%2==0 for i in range(81))
 inv=B.inv(method='DM');assert B*inv==s.eye(81)
 sample=[str(s.Mod(inv[i,i]/2,1)) for i in range(5)]
 return {'status':'PASS','result_scope':'PASS_ADDED_CYCLE_LATTICE_CS_BOUNDARY_OBSTRUCTION_AND_DOUBLING',
 'declared_action':'S=(1/4pi) integral B_IJ a^I da^J with81 compact U1 connections and K=B. This is an added2+1D bosonic Abelian Chern-Simons action, not an automatic W33 spacetime or particle model.',
 'anyon_group':'Z^81/B Z^81; invariant factors are the prior-owned Levi critical group(Z/4)^6+(Z/40)^22+Z/160.',
 'anyon_count':n,'linking_pairing':'b([x],[y])=x^T B^-1 y mod1','topological_spin':'q([x])=x^T B^-1 x/2 mod1; even B makes this well-defined.','sample_generator_spins':sample,
 'signature':81,'chiral_central_charge_mod8':1,
 'gapped_boundary_obstruction':'A Lagrangian subgroupL in a nondegenerate finite pairing satisfies|L|²=|A|. Here|A|=2^83*5^23 is not a square, so no Lagrangian subgroup and no fully gapped topological boundary of this model exists. Positive signature provides a separate chiral obstruction.',
 'doubled_repair':'For K=B direct_sum(-B), the diagonal subgroup{([x],[x])} has order|A|, trivial pairing and spin, hence is Lagrangian. The doubled model has zero signature and admits the standard diagonal boundary.',
 'computational_boundary':'Abelian braiding supplies phase data and finite topological sectors; it is not universal non-Abelian quantum computation or a derivation of gravitational dynamics.',
 'prior_owners':['analysis/w33_pass5028_5035_steinberg_apartments.py','analysis/w33_pass5436_5443_bicycle_apartment_scheme_packet.py','analysis/w33_pass11289_cycle_gram_gluing_flux.py'],
 'primary_source':'https://arxiv.org/abs/1008.0654'}
if __name__=='__main__':
 out=payload();(ROOT/'data/w33_pass11292_cycle_discriminant_boundary.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
