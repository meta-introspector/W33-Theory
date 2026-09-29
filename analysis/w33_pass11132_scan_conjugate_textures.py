"""Pass 11132/11133 scan (needs the orbifolder dumps; unlock21.py = analysis/w33_pass11131_scan_unlock_21.py in the working
directory): textures of conj(top Higgs) in all 21 neutral-exit survivors, without and with the condensate, and with a
two-scale weight r = log eps_T / log eps_S; frozen as data/w33_pass11132_conjugate_textures_21.json"""
exec(open('unlock21.py').read().split("    base = [H.mod1(S.charge_vector(s, True), nu1) for s in Sv]")[0].replace("def run(a):", "def prep(a):")
     + "    return us, th, algs, hid, nu1, Q, U, Dd, L, Ee, Hu, Hd, Sv, tvec, tcon, st" + chr(10))
from scipy.optimize import milp, LinearConstraint
RS = (0.5, 1.0, 2.0, 3.0)
def worder(EE, target, wT, nT=2):
    """min sum w_i e_i with the last nT columns (T, Tbar) weighted wT; returns (weighted total, n_S, n_T) or None"""
    t = [F(x) * EE.den for x in target]
    if any(x.denominator != 1 for x in t): return None
    b = np.array([float(int(x)) for x in t])
    c = EE.cobj.copy(); c[EE.n - nT:EE.n] = wT
    r = milp(c, constraints=[LinearConstraint(EE.A, b, b), LinearConstraint(EE.cobj.reshape(1, -1), 1, np.inf)], integrality=np.ones(EE.n + EE.nd), bounds=EE.bounds,
             options={"time_limit": 30})
    if r.status != 0: return None
    e = [int(round(x)) for x in r.x[:EE.n]]
    nt = sum(e[EE.n - nT:]); ns = sum(e) - nt
    return ns + wT * nt, ns, nt
def cl(f): return [int(x) for x in S.classes(f['n'])]
def run(a):
    lab, t = a
    us, th, algs, hid, nu1, Q, U, Dd, L, Ee, Hu, Hd, Sv, tvec, tcon, st = prep(a)
    base = [H.mod1(S.charge_vector(s, True), nu1) for s in Sv]
    E0 = S.FastOrders(base, nu1, [1] * 8)
    E1 = S.FastOrders(base + [tvec, tcon], nu1, [1] * 8)
    assert E1.vecs[-2:] == [list(tvec), list(tcon)] or E1.n == E0.n + 2
    w0 = [F(0)] * (nu1 + 8)
    key = lambda f: tuple(str(x) for x in f['q']); neg = lambda f: tuple(str(-x) for x in f['q'])
    tree_u = [k for k, h in enumerate(Hu) if any(o == 0 for o in H.order_matrix(Q, U, h, Sv, E0, nu1, hid, algs)[0].values())]
    conj = {k: [j for j, h in enumerate(Hd) if neg(h) == key(Hu[k]) and all((x + y) % 3 == 0 for x, y in zip(cl(h), cl(Hu[k])))] for k in tree_u}
    rng = np.random.default_rng(11132)
    out = dict(tree_up=tree_u, conj=conj, rows={})
    for k, js in conj.items():
        if len(js) != 1: out['rows'][k] = dict(error=len(js)); continue
        h = Hd[js[0]]; rec = {}
        for name, A, B in (('down', Q, Dd), ('lepton', L, Ee)):
            ent = {}
            for i, x in enumerate(A):
                for j, y in enumerate(B):
                    fs = [x, y, h]
                    if not Y98.hidden_ok(fs, hid, algs): continue
                    vs = [S.charge_vector(z, True) for z in fs]
                    if P97.cubic_ok(vs, nu1):
                        ent[(i, j)] = ('cubic', sum(1 for tt in range(3) if not (cl(x)[tt] == cl(y)[tt] == cl(h)[tt])))
                        continue
                    target = [F(w) - sum(F(v[c]) for v in [H.mod1(z, nu1) for z in vs]) for c, w in enumerate(w0)]
                    o0 = E0.order(target, 1)
                    ws = {r: worder(E1, target, r) for r in RS}
                    ent[(i, j)] = ('hi', o0 if isinstance(o0, int) else None, ws)
            def expo(mode, r=1.0):
                n, geo = {}, {}
                for key_, v in ent.items():
                    if v[0] == 'cubic': n[key_] = 0; geo[key_] = v[1]
                    elif mode == '0' and v[1] is not None: n[key_] = v[1]
                    elif mode == 'T' and v[2][r] is not None: n[key_] = v[2][r][0]
                return H.exponents(n, geo, len(A), len(B), 1.0, rng) if n else None
            rec[name] = dict(exp0=expo('0'), expT={str(r): expo('T', r) for r in RS},
                             nT_used=sorted({v[2][1.0][2] for v in ent.values() if v[0] == 'hi' and v[2][1.0]}))
        out['rows'][k] = dict(conj_Hd=js[0], **rec)
    return lab, out
if __name__ == '__main__':
    tor = json.load(open(R + r'\w33_pass11119_torus_resolved_tachyons_491.json'))
    labs = json.load(open(R + r'\w33_pass11119_neutral_tori_across_491.json'))['passing_gauntlet']
    jobs = [(lab, 0 if 'neutral' in tor[lab]['torus1']['kinds'] else 1) for lab in labs]
    with Pool(6) as p:
        res = {}
        for lab, o in p.imap_unordered(run, jobs):
            res[lab] = o
            print(lab, o['tree_up'], {k: (v.get('conj_Hd'), v.get('down', {}).get('exp0'), v.get('down', {}).get('expT'), v.get('lepton', {}).get('exp0'), v.get('lepton', {}).get('expT')) for k, v in o['rows'].items()}, flush=True)
            json.dump(res, open('conj21.json', 'w'), indent=1, default=str)
