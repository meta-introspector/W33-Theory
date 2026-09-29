"""Pass 11138 scan (needs the orbifolder dumps; conj21.py / unlock21.py in the working directory): conjugate-sector textures with 12 random per-singlet VEV weights; frozen as data/w33_pass11138_singlet_vev_draws.json"""
exec(open('conj21.py').read().split(chr(10) + 'def run(a):')[0])
MODELS = ['A8SM_20260951_24165_TF', 'A8SM_20260971_40521_TF', 'A8SM_20260928_2233_TF']
NDRAW = 12
def wmin(EE, target, c):
    t = [F(x) * EE.den for x in target]
    if any(x.denominator != 1 for x in t): return None
    b = np.array([float(int(x)) for x in t])
    r = milp(c, constraints=[LinearConstraint(EE.A, b, b), LinearConstraint(EE.cobj.reshape(1, -1), 1, np.inf)],
             integrality=np.ones(EE.n + EE.nd), bounds=EE.bounds, options={"time_limit": 30})
    return None if r.status != 0 else float(c @ r.x)
def run(lab):
    tor = json.load(open(R + r'\w33_pass11119_torus_resolved_tachyons_491.json'))
    t = 0 if 'neutral' in tor[lab]['torus1']['kinds'] else 1
    us, th, algs, hid, nu1, Q, U, Dd, L, Ee, Hu, Hd, Sv, tvec, tcon, st = prep((lab, t))
    base = [H.mod1(S.charge_vector(s, True), nu1) for s in Sv]
    E0 = S.FastOrders(base, nu1, [1] * 8)
    E1 = S.FastOrders(base + [tvec, tcon], nu1, [1] * 8)
    key = lambda f: tuple(str(x) for x in f['q']); neg = lambda f: tuple(str(-x) for x in f['q'])
    tree_u = [k for k, h in enumerate(Hu) if any(o == 0 for o in H.order_matrix(Q, U, h, Sv, E0, nu1, hid, algs)[0].values())]
    conj = {k: [j for j, h in enumerate(Hd) if neg(h) == key(Hu[k]) and all((x + y) % 3 == 0 for x, y in zip(cl(h), cl(Hu[k])))][0] for k in tree_u}
    rng = np.random.default_rng(11138)
    draws = []
    for dr in range(NDRAW):
        c = np.concatenate([rng.uniform(1.0, 3.0, E1.n - 2), [1.0, 1.0], np.zeros(E1.nd)])
        draws.append(c)
    w0 = [F(0)] * (nu1 + 8)
    out = {}
    for k, j in conj.items():
        h = Hd[j]; res = {}
        for name, A, B in (('down', Q, Dd), ('lepton', L, Ee)):
            ent = {}
            for i, x in enumerate(A):
                for jj, y in enumerate(B):
                    fs = [x, y, h]
                    if not Y98.hidden_ok(fs, hid, algs): continue
                    vs = [S.charge_vector(z, True) for z in fs]
                    if P97.cubic_ok(vs, nu1):
                        ent[(i, jj)] = ('cubic', sum(1 for tt in range(3) if not (cl(x)[tt] == cl(y)[tt] == cl(h)[tt])))
                        continue
                    target = [F(w) - sum(F(v[cc]) for v in [H.mod1(z, nu1) for z in vs]) for cc, w in enumerate(w0)]
                    ent[(i, jj)] = ('hi', [wmin(E1, target, c) for c in draws])
            pats = []
            for dr in range(NDRAW):
                n, geo = {}, {}
                for kk, v in ent.items():
                    if v[0] == 'cubic': n[kk] = 0; geo[kk] = v[1]
                    elif v[1][dr] is not None: n[kk] = v[1][dr]
                e = H.exponents(n, geo, len(A), len(B), 1.0, np.random.default_rng(dr)) if n else None
                pats.append(e)
            res[name] = pats
        out[str(k)] = dict(conj_Hd=j, **res)
    return lab, out
if __name__ == '__main__':
    with Pool(3) as p:
        res = dict(p.map(run, MODELS))
    json.dump(res, open('vevs.json', 'w'), indent=1, default=str)
    for m, v in res.items():
        for k, r in v.items():
            print(m, k, 'down', r['down'], flush=True)
            print(m, k, 'lepton', r['lepton'], flush=True)
