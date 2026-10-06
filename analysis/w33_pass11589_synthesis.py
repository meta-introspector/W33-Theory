# synthesis assertions for the executable post-11579 chain
dims=dict(spin10=45,so6=15,su2L=3,su2R=3,ps=21,sm=12)
assert dims['ps']==dims['so6']+dims['su2L']+dims['su2R']
assert dims['sm']==8+3+1
charges={'Q':1,'L':-3,'u^c':-4,'d^c':2,'e^c':6,'nu^c':0}
assert sum([6*charges['Q'],2*charges['L'],3*charges['u^c'],3*charges['d^c'],charges['e^c'],charges['nu^c']])==0
print('internal chain: Spin10 -> Spin6 x Spin4 -> SU3 x SU2L x SU2R x U1_B-L -> [SU3 x SU2L x U1Y]/Z6')
print('external chain: 4D overlap/Cl4 provides spacetime chirality independently')
print('single-family matter: one Spin10 Weyl16 with nu^c')
print('mass: unique 10_H channel for one family; multi-family Yukawa texture still open')
print('gravity: coarse variable-frame heat trace is not solely EH; refinement theorem still open')
print('clock: canonical internal chirality is not invariant under one/two tick dynamics')
print('PASS')
