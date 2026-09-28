import re, glob, os, sys
G = os.path.expanduser('~/orb/build7/Geometry/')
FAM = {  # family: (geometry, scan dir, old twist line regex, new twist line)
 'Z12-I':   ('Geometry_Z12-I_E6.txt',        'z12',   r'^\s*0\s+1/12\s+-5/12\s+1/3\s*$',  '    0   1/12  -5/12  -2/3'),
 'Z2xZ6-I': ('Geometry_Z2xZ6-I_SO4xG2xG2.txt','z2z6',  r'^\s*0\s+0\s+1/6\s+-1/6\s*$',      '    0     0   7/6  -1/6'),
 'Z3xZ6':   ('Geometry_Z3xZ6_SU3xG2xG2.txt', 'z3z6',  r'^\s*0\s+0\s+1/6\s+-1/6\s*$',      '    0     0   7/6  -1/6'),
 'Z6xZ6':   ('Geometry_Z6xZ6_G2^3.txt',      'z6z6',  r'^\s*0\s+1/6\s+0\s+-1/6\s*$',      '    0   7/6     0  -1/6'),
}
out = os.path.expanduser('~/orb/p1109x/fam')
os.makedirs(out, exist_ok=True)
for fam, (geo, d, old, new) in FAM.items():
    L = open(G + geo).read().splitlines()
    hits = [i for i, l in enumerate(L) if re.match(old, l)]
    # only inside the twist block
    tb = L.index('begin twist'); te = L.index('end twist')
    hits = [i for i in hits if tb < i < te]
    assert len(hits) == 1, (fam, hits)
    L[hits[0]] = new
    ngeo = geo.replace('.txt', '_NS.txt')
    open(G + ngeo, 'w').write('\n'.join(L) + '\n')
    models = []
    for f in sorted(glob.glob(os.path.expanduser(f'~/orb/scan/{d}/sm/*.txt'))):
        t = open(f).read()
        models += ['begin model' + b for b in t.split('begin model')[1:]]
    tw = [m.replace('Geometry/' + geo, 'Geometry/' + ngeo) for m in models]
    su = models
    open(f'{out}/{fam}_twins.txt', 'w').write('\n'.join(tw))
    open(f'{out}/{fam}_susy.txt', 'w').write('\n'.join(su))
    print(fam, len(models), 'models ->', ngeo)
