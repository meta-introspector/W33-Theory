#!/usr/bin/env python3
"""Pass 11104: the CKM texture and the charged leptons of the 12 surviving SO(16)xSO(16) A8 models.

For every pair (H_u, H_d) of localized light Higgs doublets with a tree-level top and bottom (Passes 11101-11102), the up
and down Yukawa matrices are built to all orders (hidden-singlet VEVs, all selection rules; engine of Pass 11102), and
    V = U_u^dagger U_d (left rotations, heaviest first),  |V_ij| ~ eps^(e_ij)  (eps = 1e-3, 1e-6; random O(1) coefficients;
    eps_geo = eps^g, g = 1, 4).
The charged leptons L ebar H_d use the SAME H_d; the neutrino Dirac couplings L N H_u use the H_u (N = SM-singlet
fermions).  Frozen in data/w33_pass11104_ckm_leptons_12.json (needs the field dumps to recompute).

Results (12 models, 9 (H_u, H_d) pairs each):
  * The heavy top and the heavy bottom lie in the SAME quark doublet only when H_u and H_d sit at the same fixed point
    of the family torus: 36 of 108 pairs.  In the other 72, |V_tb| ~ eps^2: excluded.  The vacuum must align H_u and
    H_d at one point.
  * The 36 aligned pairs give V_tb = O(1), |V_cb|, |V_ub|, |V_ts|, |V_td| ~ eps^2, and an O(1) Cabibbo block (24/36; eps
    in 12/36).  Hence V_cb ~ m_c/m_t (both eps^2 in the large-area regime, where m_c/m_t ~ eps^2): observed
    V_cb / (m_c/m_t) ~ 11 at M_Z -- an order-of-magnitude tension, not a contradiction, with O(1) coefficients;
    V_ub ~ V_cb where the data have V_ub/V_cb ~ 0.09.
  * Charged leptons: the same H_d gives a single heavy tau in 36/36 aligned pairs, with mu and e at the SAME order
    ((0,1,1) or (0,2,2)): the Delta(54) degeneracy of Passes 11102-11103 hits the leptons too; m_mu/m_e ~ 207 is not
    generated.
  * Neutrinos: 39 or 51 SM-singlet fermions and tree-level Dirac couplings L N H_u in all 12 models: light neutrinos
    need a seesaw (Majorana masses of N from the scalar VEVs; not computed here).
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
FROZEN = ROOT / "data" / "w33_pass11104_ckm_leptons_12.json"
OUT = ROOT / "data" / "w33_pass11104_ckm_and_leptons.json"


def mat(n, geo, n1, n2, g, eps, coef):
    M = np.zeros((n1, n2), dtype=complex)
    for (i, j), o in n.items():
        M[i, j] = coef[(i, j)] * eps ** o * (eps ** (g * geo.get((i, j), 0)) if o == 0 else 1.0)
    return M


def ckm_exponents(nu, gu, nd, gd, nQ, nU, nD, g, rng):
    cu = {k: rng.standard_normal() + 1j * rng.standard_normal() for k in nu}
    cd = {k: rng.standard_normal() + 1j * rng.standard_normal() for k in nd}
    Vs = []
    for eps in (1e-3, 1e-6):
        Uu = np.linalg.svd(mat(nu, gu, nQ, nU, g, eps, cu))[0]
        Ud = np.linalg.svd(mat(nd, gd, nQ, nD, g, eps, cd))[0]
        Vs.append(np.abs(Uu[:, :3].conj().T @ Ud[:, :3]))
    E = np.log(np.maximum(Vs[1], 1e-300) / np.maximum(Vs[0], 1e-300)) / np.log(1e-3)
    return np.round(E, 2).tolist()


def run(m, ours, theirs):
    import w33_pass11097_fractional_fermion_masses as P97
    import w33_pass11102_higher_order_quark_hierarchy as H
    import w33_so16_selection_rules as S
    US, TH = S.parse_ours(ours), S.parse_theirs(theirs)
    us, th = US[m], TH[m]
    rng = np.random.default_rng(11104 + m)
    algs, hid, nu1 = P97.prepare(us, th)
    c, oc, ow, qcol = S.sm_data(th, us)
    qb = P97.conj(qcol, 'A2')
    ferm = [f for f in us['fields'] if f['m'] == 2]
    scal = [f for f in us['fields'] if f['m'] == 6]
    pick = lambda L, col, w, Yv: [f for f in L if f['col'] == col and f['w'] == w and f['Y'] == Yv]
    Q, U, D = pick(ferm, qcol, 2, F(1, 6)), pick(ferm, qb, 1, F(-2, 3)), pick(ferm, qb, 1, F(1, 3))
    L, E = pick(ferm, '1', 2, F(-1, 2)), pick(ferm, '1', 1, F(1))
    Nn = [f for f in ferm if S.dimof(f['col']) == 1 and f['w'] == 1 and f['Y'] == 0 and f['hidden_singlet']]
    Hu, Hd = pick(scal, '1', 2, F(1, 2)), pick(scal, '1', 2, F(-1, 2))
    Sv = [s for s in scal if S.dimof(s['col']) == 1 and s['w'] == 1 and s['Y'] == 0 and s['hidden_singlet']]
    E0 = S.FastOrders([H.mod1(S.charge_vector(s, True), nu1) for s in Sv], nu1, [1] * 8)
    om = lambda A, B, h: H.order_matrix(A, B, h, Sv, E0, nu1, hid, algs)
    top = lambda n, geo: [k for k, o in n.items() if o == 0 and geo.get(k) == 0]
    ups = {k: om(Q, U, h) for k, h in enumerate(Hu)}
    downs = {k: om(Q, D, h) for k, h in enumerate(Hd)}
    leps = {k: om(L, E, h) for k, h in enumerate(Hd)}
    rec = dict(label=us['label'], nQ=len(Q), nU=len(U), nD=len(D), nL=len(L), nE=len(E), nN=len(Nn),
               nHu=len(Hu), nHd=len(Hd), pairs=[])
    for ku, (nu, gu) in ups.items():
        tu = top(nu, gu)
        if len(tu) != 1:
            continue
        for kd, (nd, gd) in downs.items():
            tb = top(nd, gd)
            if len(tb) != 1:
                continue
            nl, gl = leps[kd]
            row = dict(Hu=ku, Hd=kd, top_Q=tu[0][0], bottom_Q=tb[0][0], same_doublet=tu[0][0] == tb[0][0],
                       tau_heavy=len(top(nl, gl)) == 1,
                       lepton_orders={f"{i},{j}": o for (i, j), o in sorted(nl.items())},
                       lepton_geo={f"{i},{j}": o for (i, j), o in sorted(gl.items())},
                       ckm_exponents={str(g): ckm_exponents(nu, gu, nd, gd, len(Q), len(U), len(D), g, rng) for g in (1.0, 4.0)})
            if nl:
                row['lepton_exponents'] = {str(g): H.exponents(nl, gl, len(L), len(E), g, rng) for g in (1.0, 4.0)}
            rec['pairs'].append(row)
    rec['neutrino_dirac'] = {}
    for ku, (nu, gu) in ups.items():
        if len(top(nu, gu)) == 1:
            nn, _ = om(L, Nn, Hu[ku])
            rec['neutrino_dirac'][ku] = dict(n_entries=len(nn), min_order=min(nn.values()) if nn else None)
    return m, rec


def summarize():
    d = json.loads(FROZEN.read_text())
    pairs = Counter()
    ckm = {g: Counter() for g in ('1.0', '4.0')}
    lep = Counter()
    for r in d.values():
        for p in r['pairs']:
            pairs[f"same_doublet={p['same_doublet']}|tau_heavy={p['tau_heavy']}"] += 1
            if p['same_doublet']:
                for g in ckm:
                    ckm[g][json.dumps([[round(x) for x in row] for row in p['ckm_exponents'][g]])] += 1
                lep[json.dumps(p['lepton_exponents'])] += 1
    # the top (bottom) is the coupling at a common fixed point, so top_Q (bottom_Q) is the doublet at H_u's (H_d's) point:
    # 'same doublet' = 'same point'.  Each H_u has exactly one aligned H_d partner.
    one_partner = all(sum(1 for p in r['pairs'] if p['Hu'] == hu and p['same_doublet']) == 1
                      for r in d.values() for hu in {p['Hu'] for p in r['pairs']})
    res = dict(pass_id=11104, models=len(d), pairs=dict(pairs), ckm_exponents_aligned={g: dict(c) for g, c in ckm.items()},
               lepton_exponents_aligned=dict(lep), each_Hu_has_one_aligned_Hd=one_partner,
               singlet_fermions=dict(Counter(str(r['nN']) for r in d.values())),
               tree_level_neutrino_dirac=sum(1 for r in d.values() if min(v['min_order'] for v in r['neutrino_dirac'].values()) == 0),
               observed=dict(Vcb=0.041, Vub=0.0037, mc_over_mt_MZ=0.0036, Vcb_over_mc_mt=round(0.041 / 0.0036, 1)))
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    print(json.dumps(res, indent=1))
    return res


def main(ours, theirs):
    from multiprocessing import Pool
    import w33_pass11102_higher_order_quark_hierarchy as H
    with Pool(6) as p:
        out = dict(p.starmap(run, [(m, ours, theirs) for m in H.BEST]))
    FROZEN.write_text(json.dumps({str(k): v for k, v in out.items()}, indent=1, default=str))


if __name__ == "__main__":
    if len(sys.argv) > 2:
        main(*sys.argv[1:3])
    summarize()
