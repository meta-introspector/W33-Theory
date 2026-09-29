"""Pass 11135 scan (needs the orbifolder dumps; unlock21.py = analysis/w33_pass11131_scan_unlock_21.py in the working directory): Tr Q_anom over left-handed fermions and the anomalous charge of the neutral tachyon in the 21 survivors; frozen as data/w33_pass11135_fi_data_21.json"""
exec(open('unlock21.py').read().split("    base = [H.mod1(S.charge_vector(s, True), nu1) for s in Sv]")[0].replace("def run(a):", "def prep(a):")
     + "    return us, th, algs, hid, nu1, Q, U, Dd, L, Ee, Hu, Hd, Sv, tvec, tcon, st, l, qT" + chr(10))
import math
def mult(f):
    m = 1
    for x in f['dim']: m *= S.dimof(x)
    return m
def run(a):
    lab, t = a
    us, th, algs, hid, nu1, Q, U, Dd, L, Ee, Hu, Hd, Sv, tvec, tcon, st, l, qT = prep(a)
    ferm = [f for f in us['fields'] if f['m'] == 2]
    nu = len(us['u1'])
    tr = [sum(F(f['q'][k]) * mult(f) for f in ferm) for k in range(nu)]
    anom = [k for k in range(nu) if tr[k] != 0]
    assert len(anom) == 1, anom
    A = anom[0]
    uA = [F(x) for x in us['u1'][A]]
    norm2 = sum(x * x for x in uA)
    qTA = qT[A]
    # singlet VEV fields: anomalous charges (both signs available?)
    sq = sorted({F(s['q'][A]) for s in Sv})
    ratio = float(tr[A]) / (192 * math.pi ** 2 * abs(float(qTA))) if qTA != 0 else None
    return lab, dict(anom_index=A, TrQA=str(tr[A]), uA_norm2=str(norm2), qT_A=str(qTA),
                     T_over_MP_g2half=math.sqrt(0.5 * ratio) if ratio else None,
                     T_over_MP_g1=math.sqrt(ratio) if ratio else None,
                     singlet_qA_min_max=[str(sq[0]), str(sq[-1])] if sq else None,
                     sign_opposite_available=(qTA != 0), Delta=st['Delta'])
if __name__ == '__main__':
    tor = json.load(open(R + r'\w33_pass11119_torus_resolved_tachyons_491.json'))
    labs = json.load(open(R + r'\w33_pass11119_neutral_tori_across_491.json'))['passing_gauntlet']
    jobs = [(lab, 0 if 'neutral' in tor[lab]['torus1']['kinds'] else 1) for lab in labs]
    with Pool(3) as p:
        res = dict(p.map(run, jobs, chunksize=1))
    json.dump(res, open('fi.json', 'w'), indent=1, default=str)
    for k, v in res.items(): print(k, v, flush=True)
