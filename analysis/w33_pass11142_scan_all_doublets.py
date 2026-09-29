"""Pass 11142 scan (needs the orbifolder dumps; conj21.py = analysis/w33_pass11132_scan_conjugate_textures.py and unlock21.py in the working directory): up/down/lepton order matrices for all 9 Hu / 9 Hd with the condensate; frozen as data/w33_pass11142_all_doublet_textures.json"""
"""Pass 11142: all-doublet textures (up from Hu_k, down/lepton from conj(Hu_k)) with the condensate, for Higgs-mixing search"""
exec(open('conj21.py').read().split(chr(10) + 'def run(a):')[0])
MODELS = ['A8SM_20260942_36621_TF', 'A8SM_20260951_24165_TF', 'A8SM_20260971_40521_TF', 'A8SM_20260928_2233_TF']
def run(lab):
    tor = json.load(open(R + r'\w33_pass11119_torus_resolved_tachyons_491.json'))
    t = 0 if 'neutral' in tor[lab]['torus1']['kinds'] else 1
    us, th, algs, hid, nu1, Q, U, Dd, L, Ee, Hu, Hd, Sv, tvec, tcon, st = prep((lab, t))
    base = [H.mod1(S.charge_vector(s, True), nu1) for s in Sv]
    E1 = S.FastOrders(base + [tvec, tcon], nu1, [1] * 8)
    key = lambda f: tuple(str(x) for x in f['q']); neg = lambda f: tuple(str(-x) for x in f['q'])
    conj = {k: [j for j, h in enumerate(Hd) if neg(h) == key(Hu[k]) and all((x + y) % 3 == 0 for x, y in zip(cl(h), cl(Hu[k])))] for k in range(len(Hu))}
    enc = lambda n, g: dict(orders={f"{i},{j}": o for (i, j), o in n.items()}, geo={f"{i},{j}": v for (i, j), v in g.items()})
    out = dict(conj={k: v for k, v in conj.items()}, up={}, down={}, lepton={}, sizes=dict(Q=len(Q), U=len(U), D=len(Dd), L=len(L), E=len(Ee)))
    for k, h in enumerate(Hu):
        out['up'][str(k)] = enc(*H.order_matrix(Q, U, h, Sv, E1, nu1, hid, algs))
    for j, h in enumerate(Hd):
        out['down'][str(j)] = enc(*H.order_matrix(Q, Dd, h, Sv, E1, nu1, hid, algs))
        out['lepton'][str(j)] = enc(*H.order_matrix(L, Ee, h, Sv, E1, nu1, hid, algs))
    return lab, out
if __name__ == '__main__':
    with Pool(3) as p:
        res = dict(p.map(run, MODELS))
    json.dump(res, open('mix.json', 'w'), indent=1, default=str)
    print('done', flush=True)
