"""Pass 11100: build SO(16)xSO(16) base models (Z2W x ...) for the W(3,3) shifts of Z3xZ3, Z6-I and Z12-I.
For each base (its shifts V1[, V2]) find Witten shifts V0 = (a; b), a, b in (1/2)(norm-4 E8 vectors), with V0.Vk in Z for all
shifts; try candidates (standard-like first) until the non-SUSY orbifolder loads the model.  Writes base files + a table.
run in WSL: python3 so16_bases.py FAMILY"""
import itertools
import json
import os
import re
import subprocess
import sys
from fractions import Fraction as F

os.environ['LD_LIBRARY_PATH'] = os.path.expanduser('~/orb/sysroot/usr/lib/x86_64-linux-gnu')
OUT = os.path.expanduser('~/orb/p1109x/so16fam')
os.makedirs(OUT, exist_ok=True)


def e8_norm(n2):
    out = []
    for half in (False, True):
        rng = [F(k, 2) for k in range(-5, 6, 2)] if half else [F(k) for k in range(-2, 3)]
        for v in itertools.product(rng, repeat=8):
            if sum(x * x for x in v) == n2 and sum(v) % 2 == 0:
                out.append(v)
    return out


HALF = [tuple(x / 2 for x in v) for v in e8_norm(4)]
kind = lambda v: (0 if any(abs(x) == 1 for x in v) else 1 if any(x == 0 for x in v) else 2)
HALF.sort(key=kind)
dot = lambda u, v: sum(a * b for a, b in zip(u, v))


def candidates(shifts, limit=40):
    """admissible (a, b), standard-like first"""
    ra = {}
    for a in HALF:
        ra.setdefault(tuple(dot(a, V[:8]) % 1 for V in shifts), []).append(a)
    out = []
    for b in HALF:
        need = tuple((-dot(b, V[8:])) % 1 for V in shifts)
        for a in ra.get(need, [])[:3]:
            out.append(a + b)
            if len(out) >= 3 * limit:
                break
        if len(out) >= 3 * limit:
            break
    out.sort(key=lambda v: kind(v[:8]) + kind(v[8:]))
    return out[:limit]


def fmt(v):
    return ' '.join(str(x) for x in v)


def model_text(label, geom, rows):
    return f"begin model\nLabel:{label}\nSpaceGroup:Geometry/{geom}\nLattice:E8xE8\nShifts and Wilsonlines:\n" + \
        '\n'.join(fmt(r) for r in rows) + "\nend model\n"


def loads(path):
    r = subprocess.run(['./nsodump', path], cwd=os.path.expanduser('~/orb/nsobuild'), capture_output=True, text=True, timeout=600)
    g = [l for l in r.stdout.splitlines() if l.startswith('GAUGE')]
    t = [l for l in r.stdout.splitlines() if l.startswith('T_NFIELDS')]
    return (g[0] if g and 'LOADFAIL' not in r.stdout and g[0].strip() != 'GAUGE  NU1 0' else None), (t[0] if t else None)


def read_bases(family):
    Z = [F(0)] * 16
    if family == 'Z3xZ3':
        files = sorted(os.listdir(os.path.expanduser('~/orb/scan/z3z3/in')))
        out = []
        for f in files:
            t = open(os.path.expanduser('~/orb/scan/z3z3/in/' + f)).read()
            rows = t.split('Shifts and Wilsonlines:')[1].split('end model')[0].strip().splitlines()
            V1 = [F(x.strip(',')) for x in rows[0].replace(',', ' ').split()]
            V2 = [F(x.strip(',')) for x in rows[1].replace(',', ' ').split()]
            out.append((f[:-4], [V1, V2]))
        return out, 'Geometry_Z3xZ3_1_1.txt'
    if family == 'Z6-I':
        out = []
        for i in range(58):
            t = open(os.path.expanduser(f'~/orb/scan/Z6I_{i:02d}.txt')).read()
            rows = t.split('Shifts and Wilsonlines:')[1].split('end model')[0].strip().splitlines()
            out.append((f'Z6I_{i:02d}', [[F(x) for x in rows[0].split()]]))
        return out, 'Geometry_Z6-I_1_1.txt'
    if family == 'Z12-I':
        vs = []
        for f in sorted(os.listdir(os.path.expanduser('~/orb/scan/z12/sm'))):
            t = open(os.path.expanduser('~/orb/scan/z12/sm/' + f)).read()
            for blk in t.split('begin model')[1:]:
                rows = blk.split('Shifts and Wilsonlines:')[1].split('end model')[0].strip().splitlines()
                V = tuple(F(x.strip(',')) for x in rows[0].replace(',', ' ').split())
                if V not in vs:
                    vs.append(V)
        return [(f'Z12I_sm{i:02d}', [list(V)]) for i, V in enumerate(vs)], 'Geometry_Z12-I_1_1.txt'


def main(family):
    bases, geom = read_bases(family)
    table = {}
    for label, shifts in bases:
        found = None
        for V0 in candidates(shifts):
            rows = [list(V0)] + [s for s in shifts] + [[F(0)] * 16] * (2 - len(shifts)) + [[F(0)] * 16] * 6
            path = f'{OUT}/{family}_{label}.txt'
            open(path, 'w').write(model_text(f'{family}_{label}', geom, rows))
            g, t = loads(path)
            if g:
                found = dict(V0=fmt(V0), gauge=g, tachyons=t, file=path)
                break
        table[label] = found
        print(family, label, found['gauge'] if found else 'NO LOADABLE V0', found['tachyons'] if found else '', flush=True)
    json.dump(table, open(f'{OUT}/{family}_bases.json', 'w'), indent=1)


if __name__ == '__main__':
    main(sys.argv[1])
