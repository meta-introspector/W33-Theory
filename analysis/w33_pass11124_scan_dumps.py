"""Pass 11124 scan (heavy; needs the orbifolder field dumps all387_our/theirs.dump of Pass 11108, not in the repo):
couplings unlocked by the neutral condensate in the six neutral-exit survivors.  The tachyon T enters the FastOrders
engine as an extra VEV field with U(1) charges U.l (exact rationals), Witten sector k = 1, Z3 sector l = 0, R = 0 (NS
ground state, q_sh = 0) and its winding class N mod 3 in the space-group class coordinate of its own torus (winding mod
(1 - theta) Lambda is the fixed-point class).  Output frozen as data/w33_pass11124_unlock_spacegroup.json; see
analysis/w33_pass11124_condensate_unlocked_couplings.py for the reading."""
import sys, json, re
from fractions import Fraction as F
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, r'C:\Repos\Theory of Everything\analysis')
import w33_so16_selection_rules as S
import w33_pass11097_fractional_fermion_masses as P97
import w33_pass11102_higher_order_quark_hierarchy as H
import w33_pass11107_winding_tachyons as T7
import w33_so16_one_loop as E
INFO = {}
D = (r'C:/Users/wiljd/AppData/Local/Temp/rescan2/all387_our.dump', r'C:/Users/wiljd/AppData/Local/Temp/rescan2/all387_theirs.dump')
R = r'C:\Repos\Theory of Everything\data'
SIX = ['A8SM_20260942_36621_TF', 'A8SM_20260942_46043_TF', 'A8SM_20260951_24165_TF', 'A8SM_20260971_40521_TF',
       'A8SM_20260971_5904_TF', 'A8SM_20260972_17224_TF']
_P = {}
def parsed():
    if 'x' not in _P: _P['x'] = (S.parse_ours(D[0]), S.parse_theirs(D[1]))
    return _P['x']
def ext(v, w):     # append winding-class coordinate (already divided by its period 3)
    return list(v) + [F(w, 3)]
def run(lab):
    resc = json.load(open(R + r'\w33_pass11108_rescan_levels.json'))
    idx = {l: int(i) for i, l in (re.match(r'MODEL (\d+) (\S+)', x).groups() for x in open(D[0]) if x.startswith('MODEL'))}[lab]
    US, TH = parsed(); us, th = US[idx], TH[idx]
    algs, hid, nu1 = P97.prepare(us, th)
    c, oc, ow, qcol = S.sm_data(th, us); qb = P97.conj(qcol, 'A2')
    ferm = [f for f in us['fields'] if f['m'] == 2]; scal = [f for f in us['fields'] if f['m'] == 6]
    pick = lambda L, col, w, Yv: [f for f in L if f['col'] == col and f['w'] == w and f['Y'] == Yv]
    Q, U, Dd = pick(ferm, qcol, 2, F(1, 6)), pick(ferm, qb, 1, F(-2, 3)), pick(ferm, qb, 1, F(1, 3))
    L, Ee = pick(ferm, '1', 2, F(-1, 2)), pick(ferm, '1', 1, F(1))
    Nn = [f for f in ferm if S.dimof(f['col']) == 1 and f['w'] == 1 and f['Y'] == 0 and f['hidden_singlet']]
    Hu, Hd = pick(scal, '1', 2, F(1, 2)), pick(scal, '1', 2, F(-1, 2))
    Sv = [s for s in scal if S.dimof(s['col']) == 1 and s['w'] == 1 and s['Y'] == 0 and s['hidden_singlet']]
    # the neutral tachyon: charges U.l (exact rationals), k = 1, l = 0, bulk classes, R = 0, winding class 1 (conj: 2)
    rs = T7.tachyons('x', [1j, 3j], T7.RHO, rows=resc[lab]['rows']) or T7.tachyons('x', [3j, 1j], T7.RHO, rows=resc[lab]['rows'])
    l = np.array(rs[0]['l'])
    qT = [F(round(float(sum(F(ui) * F(li).limit_denominator(36) for ui, li in zip(u, l))) * 36)) / 36 for u in us['u1']]
    base = [H.mod1(S.charge_vector(s, True), nu1) for s in Sv]
    vecs0 = base
    V0_, V_, W_ = E.parse_model(resc[lab]['rows'])
    wlt = [i for i in (0, 2, 4) if np.any(W_[i] != 0)]
    N = rs[0]['N']
    tor = [t for t in (0, 1) if N[t] % 3]
    assert len(tor) == 1, N
    a = wlt[tor[0]] // 2
    cls = [F(0)] * 3; cls[a] = F(N[tor[0]] % 3, 3)
    ccl = [F(0)] * 3; ccl[a] = F((-N[tor[0]]) % 3, 3)
    tvec = qT + [F(1, 2), F(0)] + cls + [F(0)] * 3
    tcon = [-x for x in qT] + [F(1, 2), F(0)] + ccl + [F(0)] * 3
    INFO[lab] = dict(N=list(N), torus_index=a, n_tachyon_states=len(rs))
    E0 = S.FastOrders(vecs0, nu1, [1] * 8)
    E1 = S.FastOrders(vecs0 + [tvec, tcon], nu1, [1] * 8)
    w0 = [F(0)] * (nu1 + 8)
    def order(E, fields):
        vs = [H.mod1(S.charge_vector(x, True), nu1) for x in fields]
        return E.coupling(vs, w0, min_total=0)
    out = {}
    for name, combos in (('up', [(a, b, h) for a in Q for b in U for h in Hu]), ('down', [(a, b, h) for a in Q for b in Dd for h in Hd]),
                         ('lepton', [(a, b, h) for a in L for b in Ee for h in Hd]), ('nu_dirac', [(a, b, h) for a in L for b in Nn[:12] for h in Hu]),
                         ('mu_HuHd', [(a, b) for a in Hu for b in Hd])):
        new, lowered, orders, tw_new, higgs = 0, 0, {}, 0, set()
        for fs in combos:
            o0, o1 = order(E0, fs), order(E1, fs)
            if o0 is None and isinstance(o1, int):
                new += 1
                orders[o1] = orders.get(o1, 0) + 1
                tw_new += any(f['l'] != 0 for f in fs)
                higgs.add(id(fs[-1]))
            elif isinstance(o0, int) and isinstance(o1, int) and o1 < o0: lowered += 1
        out[name] = dict(combos=len(combos), newly_allowed=new, order_lowered=lowered, new_orders=orders,
                         new_with_theta_twisted_sm_field=tw_new, higgs_involved=len(higgs),
                         sm_fields_twisted=[[sum(1 for f in {id(x): x for x in col}.values() if f['l'] != 0), len({id(x) for x in col})] for col in zip(*combos)] if combos else [])
    ctl = (E0.order([6 * x for x in qT] + [F(0)] * 8), E1.order([6 * x for x in qT] + [F(0)] * 8), E1.order(list(qT) + [F(1, 2), F(0)] + cls + [F(0)] * 3))
    return lab, dict(info=INFO[lab], control=ctl, results=out)
if __name__ == '__main__':
    with Pool(6) as p:
        res = dict(p.map(run, SIX))
    json.dump(res, open('unlock6.json', 'w'), indent=1, default=str)
    for k, v in res.items(): print(k, v)
