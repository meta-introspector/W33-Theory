"""Pass 11145 scan (needs the orbifolder dumps; conj21.py / unlock21.py in the working directory): neutrino Dirac textures L N^c h for the tree-level (top) doublets, with and without the condensate; post-processed by nuan (charged-lepton conjugate textures from Pass 11142) into data/w33_pass11145_neutrino_textures.json"""
exec(open('conj21.py').read().split(chr(10) + 'def run(a):')[0])
MODELS = ['A8SM_20260942_36621_TF', 'A8SM_20260951_24165_TF', 'A8SM_20260971_40521_TF']
def run(lab):
    tor = json.load(open(R + r'\w33_pass11119_torus_resolved_tachyons_491.json'))
    t = 0 if 'neutral' in tor[lab]['torus1']['kinds'] else 1
    us, th, algs, hid, nu1, Q, U, Dd, L, Ee, Hu, Hd, Sv, tvec, tcon, st = prep((lab, t))
    ferm = [f for f in us['fields'] if f['m'] == 2]
    Nn = [f for f in ferm if S.dimof(f['col']) == 1 and f['w'] == 1 and f['Y'] == 0 and f['hidden_singlet']]
    base = [H.mod1(S.charge_vector(s, True), nu1) for s in Sv]
    E0 = S.FastOrders(base, nu1, [1] * 8); E1 = S.FastOrders(base + [tvec, tcon], nu1, [1] * 8)
    tree_u = [k for k, h in enumerate(Hu) if any(o == 0 for o in H.order_matrix(Q, U, h, Sv, E0, nu1, hid, algs)[0].values())]
    enc = lambda n, g: dict(orders={f"{i},{j}": o for (i, j), o in n.items()}, geo={f"{i},{j}": v for (i, j), v in g.items()})
    out = dict(n_N=len(Nn), L_twisted=[i for i, f in enumerate(L) if f['l'] != 0], rows={})
    for k in tree_u:
        out['rows'][str(k)] = dict(nu0=enc(*H.order_matrix(L, Nn, Hu[k], Sv, E0, nu1, hid, algs)),
                                   nuT=enc(*H.order_matrix(L, Nn, Hu[k], Sv, E1, nu1, hid, algs)))
    return lab, out
if __name__ == '__main__':
    with Pool(3) as p:
        res = dict(p.map(run, MODELS))
    json.dump(res, open('nu.json', 'w'), indent=1, default=str)
    print('done', flush=True)
