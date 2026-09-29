"""Pass 11127/11130 scan (needs the Pass 11108 orbifolder dumps and the Pass 11124 scan file in the working directory as
unlock6.py): down and lepton textures of the six with and without the condensate, per Higgs doublet, using the Pass 11102
order_matrix (hidden-gauge check + exact cubic rules); frozen as data/w33_pass11127_textures_six.json"""
exec(open('unlock6.py').read().split("    INFO[lab] = dict(")[0].replace("def run(lab):", "def prep(lab):")
     + "    return us, th, algs, hid, nu1, Q, U, Dd, L, Ee, Hu, Hd, Sv, vecs0, tvec, tcon" + chr(10))
import w33_pass11098_yukawa_textures as Y98
def run(lab):
    us, th, algs, hid, nu1, Q, U, Dd, L, Ee, Hu, Hd, Sv, vecs0, tvec, tcon = prep(lab)
    E0 = S.FastOrders(vecs0, nu1, [1] * 8)
    E1 = S.FastOrders(vecs0 + [tvec, tcon], nu1, [1] * 8)
    rng = np.random.default_rng(11127)
    ltw = [i for i, f in enumerate(L) if f['l'] != 0]
    out = dict(L_twisted=ltw, n_L=len(L), n_E=len(Ee), n_Hd=len(Hd), higgs=[])
    for k, h in enumerate(Hd):
        rec = dict(k=k, h_twisted=h['l'] != 0, h_classes=S.classes(h['n']))
        for name, A, B in (('down', Q, Dd), ('lepton', L, Ee)):
            for tag, EE in (('0', E0), ('T', E1)):
                n, geo = H.order_matrix(A, B, h, Sv, EE, nu1, hid, algs)
                rec[f'{name}_{tag}'] = dict(orders={f"{i},{j}": o for (i, j), o in sorted(n.items())},
                                            geo={f"{i},{j}": v for (i, j), v in sorted(geo.items())},
                                            exp=H.exponents(n, geo, len(A), len(B), 1.0, rng) if n else None,
                                            min_order=min(n.values()) if n else None)
                if name == 'lepton' and n:
                    sub = {(ltw.index(i), j): o for (i, j), o in n.items() if i in ltw}
                    sg = {(ltw.index(i), j): geo.get((i, j), 0) for (i, j) in n if i in ltw}
                    rec[f'{name}_{tag}']['exp_twisted_block'] = H.exponents(sub, sg, len(ltw), len(B), 1.0, rng) if sub else None
        out['higgs'].append(rec)
    return lab, out
if __name__ == '__main__':
    with Pool(6) as p:
        res = dict(p.map(run, SIX))
    json.dump(res, open('lep.json', 'w'), indent=1, default=str)
    for lab, v in res.items():
        print(lab, 'Ltw', v['L_twisted'], v['n_L'], v['n_E'])
        for r in v['higgs']:
            print('   H', r['k'], 'down0', r['down_0']['min_order'], r['down_0']['exp'], 'downT', r['down_T']['min_order'], r['down_T']['exp'],
                  '| lep0', r['lepton_0']['min_order'], 'lepT', r['lepton_T']['min_order'], r['lepton_T']['exp'], r['lepton_T'].get('exp_twisted_block'))
