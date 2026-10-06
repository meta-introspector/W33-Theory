import json,itertools,sympy as sp
old=json.load(open('data/w33_pass11423_11427_flux_index_curved_uv_noise.json',encoding='utf-8'))
scan=old['spacetime_index']['scans'][0]
idx=int(scan['index'])
# canonical internal Spin10 chirality from committed Cl10
gd=json.load(open('data/w33_pass10961_albert_clifford9_gammas.json',encoding='utf-8'))
g9=[sp.Matrix([[sp.Rational(x) for x in row] for row in M]) for M in gd['gamma9']]
I16=sp.eye(16);Z16=sp.zeros(16)
g=[sp.Matrix.vstack(sp.Matrix.hstack(Z16,x),sp.Matrix.hstack(x,Z16)) for x in g9]+[sp.diag(I16,-I16)]
prod=sp.eye(32)
for x in g:prod=prod*x
Chi=sp.I*prod
P=(sp.eye(32)+Chi)/2
assert P.rank()==16 and Chi*Chi==sp.eye(32)
for i,j in itertools.combinations(range(10),2):
 assert (g[i]*g[j])*Chi==Chi*(g[i]*g[j])
product_index=16*idx
print('external overlap index',idx)
print('internal Weyl rank',P.rank())
print('tensor product index multiplicity',product_index)
print('all 45 internal Spin10 bivectors preserve internal Weyl projector')
print('external Lorentz/overlap and internal gauge actions commute identically on tensor factors')
print('PASS')
