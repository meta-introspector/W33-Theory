"""Pass 11131 scan (needs the orbifolder dumps a8_our_v2/a8_theirs and all387): couplings unlocked by the lowest SM-neutral tachyon in all 21 neutral-exit survivors, hidden-gauge check enforced; frozen as data/w33_pass11131_unlock_scan_21.json"""
import sys, json, re
from fractions import Fraction as F
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, r'C:\Repos\Theory of Everything\analysis')
import w33_so16_selection_rules as S
import w33_pass11097_fractional_fermion_masses as P97
import w33_pass11102_higher_order_quark_hierarchy as H
import w33_pass11107_winding_tachyons as T7
import w33_pass11098_yukawa_textures as Y98
import w33_so16_one_loop as E
R = r'C:\Repos\Theory of Everything\data'
D1 = (r'C:/Users/wiljd/AppData/Local/Temp/a8_our_v2.dump', r'C:/Users/wiljd/AppData/Local/Temp/a8_theirs.dump')
D2 = (r'C:/Users/wiljd/AppData/Local/Temp/rescan2/all387_our.dump', r'C:/Users/wiljd/AppData/Local/Temp/rescan2/all387_theirs.dump')
_P = {}
def parsed(src):
    if src not in _P: _P[src] = (S.parse_ours(src[0]), S.parse_theirs(src[1]))
    return _P[src]
def index(src):
    return {l: int(i) for i, l in (re.match(r'MODEL (\d+) (\S+)', x).groups() for x in open(src[0]) if x.startswith('MODEL'))}
def rows_of(lab):
    r = json.load(open(R + r'\w33_pass11108_rescan_levels.json'))
    if lab in r: return r[lab]['rows'], D2
    s = json.load(open(R + r'\w33_pass11106_survivor_shifts.json'))['models']
    return next(v['rows'] for v in s.values() if v['label'] == lab), D1
def run(a):
    lab, t = a
    rows, src = rows_of(lab)
    US, TH = parsed(src); idx = index(src)[lab]; us, th = US[idx], TH[idx]
    algs, hid, nu1 = P97.prepare(us, th)
    c, oc, ow, qcol = S.sm_data(th, us); qb = P97.conj(qcol, 'A2')
    ferm = [f for f in us['fields'] if f['m'] == 2]; scal = [f for f in us['fields'] if f['m'] == 6]
    pick = lambda L, col, w, Yv: [f for f in L if f['col'] == col and f['w'] == w and f['Y'] == Yv]
    Q, U, Dd = pick(ferm, qcol, 2, F(1, 6)), pick(ferm, qb, 1, F(-2, 3)), pick(ferm, qb, 1, F(1, 3))
    L, Ee = pick(ferm, '1', 2, F(-1, 2)), pick(ferm, '1', 1, F(1))
    Nn = [f for f in ferm if S.dimof(f['col']) == 1 and f['w'] == 1 and f['Y'] == 0 and f['hidden_singlet']]
    Hu, Hd = pick(scal, '1', 2, F(1, 2)), pick(scal, '1', 2, F(-1, 2))
    Sv = [s for s in scal if S.dimof(s['col']) == 1 and s['w'] == 1 and s['Y'] == 0 and s['hidden_singlet']]
    Tw = [1j, 3j] if t == 0 else [3j, 1j]
    Uf = np.array([[float(x) for x in u] for u in us['u1']])
    fac = [np.array([[float(x) for x in r] for r in roots]) for _, roots in us['factors']]
    neu = []
    for r in T7.tachyons('x', Tw, T7.RHO, rows=rows):
        l = np.array(r['l']); Yv = float(sum(float(ci) * qi for ci, qi in zip(c, Uf @ l)))
        if abs(Yv) < 1e-9 and np.all(np.abs(fac[ow] @ l) < 1e-9) and np.all(np.abs(fac[oc] @ l) < 1e-9):
            neu.append(r)
    st = min(neu, key=lambda r: r['Delta'])
    l = [F(x).limit_denominator(36) for x in st['l']]
    assert max(abs(float(a) - b) for a, b in zip(l, st['l'])) < 1e-9
    qT = [sum(F(ui) * li for ui, li in zip(u, l)) for u in us['u1']]
    V0_, V_, W_ = E.parse_model(rows)
    wlt = [i for i in (0, 2, 4) if np.any(W_[i] != 0)]
    N = st['N']; a_ = wlt[t] // 2
    assert N[1 - t] % 3 == 0 and N[t] % 3 != 0, N
    cls = [F(0)] * 3; cls[a_] = F(N[t] % 3, 3)
    ccl = [F(0)] * 3; ccl[a_] = F((-N[t]) % 3, 3)
    tvec = qT + [F(1, 2), F(0)] + cls + [F(0)] * 3
    tcon = [-x for x in qT] + [F(1, 2), F(0)] + ccl + [F(0)] * 3
    base = [H.mod1(S.charge_vector(s, True), nu1) for s in Sv]
    E0 = S.FastOrders(base, nu1, [1] * 8); E1 = S.FastOrders(base + [tvec, tcon], nu1, [1] * 8)
    w0 = [F(0)] * (nu1 + 8)
    order = lambda EE, fs: EE.coupling([H.mod1(S.charge_vector(x, True), nu1) for x in fs], w0, min_total=0)
    out = {}
    for name, combos in (('up', [(a, b, h) for a in Q for b in U for h in Hu]), ('down', [(a, b, h) for a in Q for b in Dd for h in Hd]),
                         ('lepton', [(a, b, h) for a in L for b in Ee for h in Hd]), ('nu_dirac', [(a, b, h) for a in L for b in Nn[:12] for h in Hu]),
                         ('mu_HuHd', [(a, b) for a in Hu for b in Hd])):
        new, new_hid, orders, tw = 0, 0, {}, 0
        for fs in combos:
            o0, o1 = order(E0, fs), order(E1, fs)
            if o0 is None and isinstance(o1, int):
                new += 1
                if Y98.hidden_ok(list(fs), hid, algs):
                    new_hid += 1; orders[o1] = orders.get(o1, 0) + 1; tw += any(f['l'] != 0 for f in fs)
        out[name] = dict(combos=len(combos), newly_allowed=new, newly_allowed_hidden_ok=new_hid, orders=orders, with_twisted=tw)
    ctl = (E0.order([6 * x for x in qT] + [F(0)] * 8), E1.order([6 * x for x in qT] + [F(0)] * 8))
    return lab, dict(torus=t, N=list(N), n_neutral_states=len(neu), Delta=st['Delta'], control=ctl, results=out)
if __name__ == '__main__':
    tor = json.load(open(R + r'\w33_pass11119_torus_resolved_tachyons_491.json'))
    labs = json.load(open(R + r'\w33_pass11119_neutral_tori_across_491.json'))['passing_gauntlet']
    jobs = [(lab, 0 if 'neutral' in tor[lab]['torus1']['kinds'] else 1) for lab in labs]
    with Pool(4) as p:
        res = dict(p.map(run, jobs, chunksize=1))
    json.dump(res, open('unlock21.json', 'w'), indent=1, default=str)
    for k, v in res.items():
        print(k, v['torus'], v['N'], v['control'], {n: (r['newly_allowed'], r['newly_allowed_hidden_ok']) for n, r in v['results'].items()}, flush=True)
