"""Pass 11094: Standard-Model quantum numbers of the theta-sector tachyons of the 87 Z6-I twins.
inputs: smbasis dump (SM configuration: factors, U(1)s, labelled fields) and the mass-level -1/6 dump (tachyon fields in the
standard configuration, with its U(1) basis).  Exact rational arithmetic throughout."""
import json
import re
import sys
from collections import Counter
from fractions import Fraction as F

import sympy as sp

SMY = {"q": F(1, 6), "bu": F(-2, 3), "bd": F(1, 3), "l": F(-1, 2), "be": F(1), "bl": F(1, 2), "d": F(-1, 3),
       "u": F(2, 3), "bq": F(-1, 6), "e": F(-1)}


def ex(x):
    f = F(float(x)).limit_denominator(10 ** 5)
    assert abs(float(f) - float(x)) < 1e-9, x
    assert 1296 % f.denominator == 0, (x, f)
    return f


def dimof(t):
    return int(re.match(r'-?(\d+)', t).group(1))


def parse_blocks(path):
    blocks, cur = {}, None
    for ln in open(path, errors='ignore'):
        m = re.match(r'MODEL (\d+) (\S+)', ln)
        if m:
            cur = int(m.group(1))
            blocks[cur] = dict(label=m.group(2), factors=[], u1=[], fields=[])
            continue
        if cur is None:
            continue
        b = blocks[cur]
        m = re.match(r'FACTOR (\d+) (\S+) \|(.*)', ln)
        if m:
            roots = [[ex(x) for x in r.split()] for r in m.group(3).split('|')]
            b['factors'].append((m.group(2), roots))
            continue
        m = re.match(r'U1(?:STD)? (\d+) \|(.*)', ln)
        if m:
            b['u1'].append([ex(x) for x in m.group(2).split()])
            continue
        m = re.match(r'L (\S+) dim=(\S+) q=(\S+)', ln)
        if m:
            b['fields'].append(dict(lab=m.group(1), dim=m.group(2).split(','), q=[ex(x) for x in m.group(3).split(',')]))
            continue
        m = re.match(r'S k=(\d+) l=(\d+) m=(\d+) dim=(\S+) q=(\S+)', ln)
        if m:
            b['fields'].append(dict(k=int(m.group(1)), l=int(m.group(2)), m=int(m.group(3)), dim=m.group(4).split(','),
                                    q=[ex(x) for x in m.group(5).split(',')]))
    return blocks


def R(x):
    return sp.Rational(x.numerator, x.denominator)


def main(smb, lev, out):
    SM = parse_blocks(smb)
    LV = parse_blocks(lev)
    res = {}
    for i in sorted(SM):
        s, t = SM[i], LV[i]
        assert [a for a, _ in s['factors']] == [a for a, _ in t['factors']]
        assert all(ra == rb for (_, ra), (_, rb) in zip(s['factors'], t['factors']))   # same factors, same order
        lab = [f for f in s['fields'] if f['lab'] in SMY]
        c3 = next(j for j, x in enumerate(next(f for f in lab if f['lab'] == 'q')['dim']) if dimof(x) == 3)
        w2 = next(j for j, x in enumerate(next(f for f in lab if f['lab'] == 'l')['dim']) if dimof(x) == 2)
        A = sp.Matrix([[R(x) for x in f['q']] for f in lab])
        bY = sp.Matrix([R(SMY[f['lab']]) for f in lab])
        y, params = A.gauss_jordan_solve(bY)
        y = y.subs({p: 0 for p in params})
        assert A * y == bY
        TSM = sp.Matrix([[R(x) for x in u] for u in s['u1']])        # nSM x 16
        tY = TSM.T * y                                                # 16
        TST = sp.Matrix([[R(x) for x in u] for u in t['u1']])        # nstd x 16
        c, params = TST.T.gauss_jordan_solve(tY)
        c = c.subs({p: 0 for p in params})
        assert TST.T * c == tY
        # orthogonality check: tY is orthogonal to every non-Abelian root
        for _, roots in s['factors']:
            for r in roots:
                assert sum(R(a) * b for a, b in zip(r, tY)) == 0
        tach = [f for f in t['fields'] if (f['k'], f['l']) == (1, 0)]
        rows = []
        for f in tach:
            Y = sum(R(q) * c[j] for j, q in enumerate(f['q']))
            col, w = dimof(f['dim'][c3]), dimof(f['dim'][w2])
            hidden = [f['dim'][j] for j in range(len(f['dim'])) if j not in (c3, w2)]
            charges = [Y + sp.Rational(tt, 2) for tt in range(-(w - 1), w, 2)]
            rows.append(dict(colour=f['dim'][c3], su2=f['dim'][w2], Y=str(Y), hidden=hidden,
                             mult=(dimof(f['dim'][c3]) * w * int(sp.prod([dimof(h) for h in hidden]))),
                             sm_singlet=(col == 1 and w == 1 and Y == 0),
                             has_neutral_colourless_component=(col == 1 and any(q == 0 for q in charges)),
                             fractional_charge=(col == 1 and any(sp.Rational(q).q != 1 for q in charges))))
        res[i] = dict(label=s['label'], tY=[str(x) for x in tY], colour_roots=[[str(x) for x in r] for r in s['factors'][c3][1]],
                      su2_roots=[[str(x) for x in r] for r in s['factors'][w2][1]], tachyon_fields=len(rows), tachyon_states=sum(r['mult'] for r in rows), fields=rows,
                      all_sm_singlet=all(r['sm_singlet'] for r in rows),
                      coloured=sum(1 for r in rows if r['colour'] not in ('1',)),
                      sm_charged=sum(1 for r in rows if not r['sm_singlet']),
                      without_neutral_component=sum(1 for r in rows if not r['has_neutral_colourless_component']),
                      fractional=sum(1 for r in rows if r['fractional_charge']))
    json.dump(res, open(out, 'w'), indent=1, default=str)
    v = list(res.values())
    print('models', len(v))
    print('all tachyons SM singlets:', sum(x['all_sm_singlet'] for x in v))
    print('models with a coloured tachyon:', sum(1 for x in v if x['coloured']))
    print('models with an SM-charged tachyon:', sum(1 for x in v if x['sm_charged']))
    print('models with a tachyon lacking a neutral colourless component:', sum(1 for x in v if x['without_neutral_component']))
    print('models with a fractionally charged tachyon:', sum(1 for x in v if x['fractional']))
    reps = Counter((r['colour'], r['su2'], r['Y']) for x in v for r in x['fields'])
    print('tachyon SM reps (colour, SU2, Y):', reps.most_common(12))


if __name__ == '__main__':
    main(*sys.argv[1:4])
