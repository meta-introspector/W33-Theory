"""Pass 11130 scan (needs the orbifolder dumps and unlock6.py = analysis/w33_pass11124_scan_dumps.py in the working
directory): the conjugation map between the Y = +1/2 and Y = -1/2 scalar doublets of the six (U(1) charges negated,
space-group classes negated) and the doublets with a tree-level up coupling; frozen as
data/w33_pass11130_higgs_conjugation_six.json"""
exec(open('unlock6.py').read().split("    INFO[lab] = dict(")[0].replace("def run(lab):", "def prep(lab):")
     + "    return us, th, algs, hid, nu1, Q, U, Dd, L, Ee, Hu, Hd, Sv, vecs0, tvec, tcon" + chr(10))
import json
out = {}
for lab in SIX:
    us, th, algs, hid, nu1, Q, U, Dd, L, Ee, Hu, Hd, Sv, vecs0, tvec, tcon = prep(lab)
    key = lambda f: tuple(str(x) for x in f['q'])
    neg = lambda f: tuple(str(-x) for x in f['q'])
    hu = {key(f): i for i, f in enumerate(Hu)}
    pairs = {k: hu.get(neg(h)) for k, h in enumerate(Hd)}
    # which Hu have a tree-level up coupling: reuse 11102 order_matrix
    E0 = S.FastOrders(vecs0, nu1, [1] * 8)
    tree_u = [k for k, h in enumerate(Hu) if any(o == 0 for o in H.order_matrix(Q, U, h, Sv, E0, nu1, hid, algs)[0].values())]
    cl = lambda f: [int(x) for x in S.classes(f['n'])]
    conj_of_tree = {k: [j for j, h in enumerate(Hd) if neg(h) == key(Hu[k]) and all((a + b) % 3 == 0 for a, b in zip(cl(h), cl(Hu[k])))] for k in tree_u}
    conj_of_tree_eq = {k: [j for j, h in enumerate(Hd) if neg(h) == key(Hu[k]) and cl(h) == cl(Hu[k])] for k in tree_u}
    out[lab] = dict(conj_of_tree_neg=conj_of_tree, conj_of_tree_eq=conj_of_tree_eq, Hu_cl=[cl(h) for h in Hu], Hd_cl=[cl(h) for h in Hd], n_Hu=len(Hu), n_Hd=len(Hd), Hd_to_conj_Hu=pairs, tree_up_Hu=tree_u,
                    Hd_l=[h['l'] for h in Hd], Hu_l=[h['l'] for h in Hu], Hd_k=[h['k'] for h in Hd], Hu_k=[h['k'] for h in Hu])
    print(lab, out[lab], flush=True)
json.dump(out, open('hconj2.json', 'w'), indent=1)
