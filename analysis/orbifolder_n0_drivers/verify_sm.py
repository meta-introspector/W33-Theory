"""Independent check of a non-SUSY-orbifolder SM-like model: take ITS hypercharge (16D, from its labelled fields) and apply
it to OUR fermion spectrum (patched orbifolder 1.2.1, levdump at mu = 0).  Factors are matched by root span.
usage: verify_sm.py nso_sm.dump our_lev.dump [model_index]"""
import json
import re
import sys
from collections import Counter
from fractions import Fraction as F

import sympy as sp

sys.path.insert(0, sys.path[0])
from tach_charges import parse_blocks, R, dimof  # noqa: E402

SMY = {"q": F(1, 6), "bu": F(-2, 3), "bd": F(1, 3), "l": F(-1, 2), "be": F(1), "bl": F(1, 2), "d": F(-1, 3),
       "u": F(2, 3), "bq": F(-1, 6), "e": F(-1), "bn": F(0), "n": F(0)}


def span_equal(A, B):
    MA, MB = sp.Matrix([[R(x) for x in r] for r in A]), sp.Matrix([[R(x) for x in r] for r in B])
    return MA.rank() == MB.rank() == sp.Matrix.vstack(MA, MB).rank()


def verify(th, us):
    lab = [f for f in th['fields'] if f['lab'] in SMY and f['lab'] not in ('n', 'bn')]
    c3 = next(j for j, x in enumerate(next(f for f in lab if f['lab'] == 'q')['dim']) if dimof(x) == 3)
    w2 = next(j for j, x in enumerate(next(f for f in lab if f['lab'] == 'l')['dim']) if dimof(x) == 2)
    A = sp.Matrix([[R(x) for x in f['q']] for f in lab])
    bY = sp.Matrix([R(SMY[f['lab']]) for f in lab])
    y, params = A.gauss_jordan_solve(bY)
    y = y.subs({p: 0 for p in params})
    assert A * y == bY
    tY = sp.Matrix([[R(x) for x in u] for u in th['u1']]).T * y
    YY = (tY.T * tY)[0]
    # our factor matching colour / SU(2)_L
    oc = [j for j, (a, roots) in enumerate(us['factors']) if span_equal(roots, th['factors'][c3][1])]
    ow = [j for j, (a, roots) in enumerate(us['factors']) if span_equal(roots, th['factors'][w2][1])]
    assert len(oc) == 1 and len(ow) == 1, (oc, ow)
    oc, ow = oc[0], ow[0]
    TST = sp.Matrix([[R(x) for x in u] for u in us['u1']])
    c, params = TST.T.gauss_jordan_solve(tY)
    c = c.subs({p: 0 for p in params})
    assert TST.T * c == tY
    cnt, scal = Counter(), Counter()
    frac_chiral = Counter()
    for f in us['fields']:
        if (f['k'], f['l']) == (0, 0) and f['m'] in (4, 5):
            continue
        Y = sum(R(q) * c[j] for j, q in enumerate(f['q']))
        key = (f['dim'][oc], dimof(f['dim'][ow]), Y)
        hid = 1
        for j, x in enumerate(f['dim']):
            if j not in (oc, ow):
                hid *= dimof(x)
        if f['m'] == 2:
            cnt[key] += hid
        elif f['m'] == 6:
            scal[key] += hid
    # net chirality: pair each rep with its conjugate (colour token negated, Y negated)
    conj = lambda k: (k[0][1:] if k[0].startswith('-') else ('-' + k[0] if k[0] not in ('1',) and not k[0].endswith('adj') else k[0]), k[1], -k[2])
    net = {}
    for k in set(cnt) | {conj(k) for k in cnt}:
        n = cnt.get(k, 0) - cnt.get(conj(k), 0)
        if n > 0:
            net[k] = n
    # electric charges of the net chiral states
    frac = {str(k): v for k, v in net.items() if dimof(k[0]) == 1 and any((k[2] + sp.Rational(t, 2)).q != 1 for t in range(-(k[1] - 1), k[1], 2))}
    higgs = sum(v for k, v in scal.items() if dimof(k[0]) == 1 and k[1] == 2 and abs(k[2]) == sp.Rational(1, 2))
    isfrac = lambda k: dimof(k[0]) == 1 and any((k[2] + sp.Rational(t, 2)).q != 1 for t in range(-(k[1] - 1), k[1], 2))
    frac_total = sum(v for k, v in cnt.items() if isfrac(k))
    frac_scal = sum(v for k, v in scal.items() if isfrac(k))
    return dict(YY=str(YY), net_chiral={f"({k[0]},{k[1]})_{k[2]}": v for k, v in sorted(net.items(), key=str)},
                net_fractional_colourless=frac, higgs_doublet_scalars=higgs,
                fractional_colourless_fermion_states=frac_total, fractional_colourless_scalar_states=frac_scal,
                sm_scalars_like_matter={f"({k[0]},{k[1]})_{k[2]}": v for k, v in scal.items() if dimof(k[0]) == 3})


if __name__ == '__main__':
    TH = parse_blocks(sys.argv[1])
    US = parse_blocks(sys.argv[2])
    out = {}
    for i in sorted(TH):
        if not TH[i]['factors']:
            continue
        r = verify(TH[i], US[i])
        out[i] = r
        print(i, TH[i]['label'], 'Y.Y =', r['YY'], '| net chiral:', r['net_chiral'], '| net fractional colourless:', r['net_fractional_colourless'], '| Higgs scalars:', r['higgs_doublet_scalars'])
    if len(sys.argv) > 3:
        json.dump(out, open(sys.argv[3], 'w'), indent=1, default=str)
