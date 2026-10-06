import json
p50=json.load(open('data/w33_pass10950_clock_albert_lorentz_spinor.json',encoding='utf-8'))
tri=json.load(open('data/w33_e6_trinification_schlafli.json',encoding='utf-8'))
assert p50['lorentz']['spinor_module_dimension']==16
decomp=tri['trinification']['decomposition']
assert '9+9+9' in decomp or '9 + 9 + 9' in decomp
one27=27;spinor=16;vector=10;singlet=1
assert spinor+vector+singlet==one27
nonets=[9,9,9]
assert sum(nonets)==27 and all(x!=16 for x in nonets)
three_families=3*27
print('same E6 27: Spin10=',[16,10,1],'trinification=',nonets)
print('three independent E6 families require',three_families,'states before symmetry breaking, not 27')
print('CORRECTION: trinification S3 permutes gauge factors/nonets inside one 27; it is not by itself a 3-generation family symmetry')
print('PASS')
