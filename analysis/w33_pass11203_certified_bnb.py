"""Pass 11203: a rounding-certified replay of the Pass 11187 branch and bound (4/sqrt15 on the torus-symmetric class).

Pass 11187 closed every box of the torus-symmetric class x = (s1, s2, t1, t2), c = 1 - sum(x) >= 0, with
    f(x) = A(s1,t1) + A(s2,t2) + sqrt(s1 s2) + A(t1,s1) + A(t2,s2) + sqrt(t1 t2),   A(p,q) = (sqrt(q^2 + 4cp) - q)/2,
using float64 with a 1e-9 margin (the prior-art sweep of Passes 11188-11192 flagged: not directed rounding).

Here the replay is rigorous: every closing decision is evaluated with outward rounding.  IEEE-754 +, -, *, / and sqrt
are correctly rounded (error <= 1/2 ulp), so nudging each result one ulp outward with nextafter gives a guaranteed
enclosure.  Only two closing rules are used:

  * bound: the upward-rounded monotone bound (A increasing in p and c, decreasing in q) is < 4/sqrt15 rounded down;
  * concavity core: the Pass 11187 interval Hessian test, evaluated with directed rounding on k^3 cells whose union
    provably covers the hull (endpoints pinned), after which Pass 11187/11167's slice argument applies unchanged.

Infeasible boxes are discarded only if sum(lo), rounded down, exceeds 1.  Box splitting shares the midpoint exactly.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "w33_pass11203_certified_bnb.json"
INF = np.inf


def up(x):
    return np.nextafter(x, INF)


def dn(x):
    return np.nextafter(x, -INF)


TARGET_LO = dn(4.0 / up(up(np.sqrt(15.0))))    # value strictly below 4/sqrt15


def A_up(p_hi, q_lo, c_hi):
    """upper bound of (sqrt(q^2 + 4 c p) - q)/2 for p <= p_hi, q >= q_lo, c <= c_hi (all >= 0)"""
    g = up(np.sqrt(up(up(q_lo * q_lo) + up(up(4.0 * c_hi) * p_hi))))
    return up(up(g - q_lo) / 2.0)


def monotone_up(lo, hi):
    s1l, s2l, t1l, t2l = lo.T
    s1h, s2h, t1h, t2h = hi.T
    c_hi = np.maximum(up(1.0 - dn(dn(dn(s1l + s2l) + t1l) + t2l)), 0.0)
    tot = up(A_up(s1h, t1l, c_hi) + A_up(s2h, t2l, c_hi))
    tot = up(tot + up(np.sqrt(up(s1h * s2h))))
    tot = up(tot + A_up(t1h, s1l, c_hi))
    tot = up(tot + A_up(t2h, s2l, c_hi))
    tot = up(tot + up(np.sqrt(up(t1h * t2h))))
    return tot


def hessian_test_rig(s_lo, s_hi, t_lo, t_hi, c_lo, c_hi):
    """directed-rounding version of the Pass 11187 test: H11 < 0, H22 < 0 and H12^2 < H11*H22 on the whole cell"""
    ok = (c_lo > 0) & (s_lo > 0) & (t_lo > 0)
    s_lo, t_lo, c_lo = np.where(ok, s_lo, 1.0), np.where(ok, t_lo, 1.0), np.where(ok, c_lo, 1.0)
    G1_lo = dn(np.sqrt(dn(dn(t_lo * t_lo) + dn(dn(4.0 * c_lo) * s_lo))))
    G1_hi = up(np.sqrt(up(up(t_hi * t_hi) + up(up(4.0 * c_hi) * s_hi))))
    G2_lo = dn(np.sqrt(dn(dn(s_lo * s_lo) + dn(dn(4.0 * c_lo) * t_lo))))
    G2_hi = up(np.sqrt(up(up(s_hi * s_hi) + up(up(4.0 * c_hi) * t_hi))))
    cube_up = lambda g: up(up(g * g) * g)
    cube_dn = lambda g: dn(dn(g * g) * g)
    # H11 <= -4 c_lo^2 / G1_hi^3 + 4 c_hi t_hi / G2_lo^3
    neg1 = dn(dn(dn(4.0 * c_lo) * c_lo) / cube_up(G1_hi))
    pos1 = up(up(up(4.0 * c_hi) * t_hi) / cube_dn(G2_lo))
    H11_hi = up(pos1 - neg1)
    neg2 = dn(dn(dn(4.0 * c_lo) * c_lo) / cube_up(G2_hi))
    pos2 = up(up(up(4.0 * c_hi) * s_hi) / cube_dn(G1_lo))
    H22_hi = up(pos2 - neg2)
    H12_abs = up(up(up(2.0 * c_hi) * t_hi) / cube_dn(G1_lo) + up(up(2.0 * c_hi) * s_hi) / cube_dn(G2_lo))
    H12_abs = up(H12_abs)
    # H11, H22 both negative: |H11| >= -H11_hi, |H22| >= -H22_hi, so H11*H22 >= H11_hi*H22_hi (product of negatives)
    prod_dn = dn(H11_hi * H22_hi)
    return ok & (H11_hi < 0) & (H22_hi < 0) & (up(H12_abs * H12_abs) < prod_dn)


def concave_core_rig(lo, hi, k=6):
    s_lo = np.minimum(lo[:, 0], lo[:, 1]); s_hi = np.maximum(hi[:, 0], hi[:, 1])
    t_lo = np.minimum(lo[:, 2], lo[:, 3]); t_hi = np.maximum(hi[:, 2], hi[:, 3])
    c_lo = np.maximum(dn(1.0 - up(up(up(hi[:, 0] + hi[:, 1]) + hi[:, 2]) + hi[:, 3])), 0.0)
    c_hi = np.maximum(up(1.0 - dn(dn(dn(lo[:, 0] + lo[:, 1]) + lo[:, 2]) + lo[:, 3])), 0.0)

    def cuts(a, b):
        pts = [a + (b - a) * (i / k) for i in range(k + 1)]
        pts[0], pts[-1] = a, b                     # pin endpoints: the cells provably cover [a, b]
        return pts
    S, T, C = cuts(s_lo, s_hi), cuts(t_lo, t_hi), cuts(c_lo, c_hi)
    ok = np.ones(len(lo), bool)
    for a in range(k):
        for b in range(k):
            for e in range(k):
                ok &= hessian_test_rig(np.minimum(S[a], S[a + 1]), np.maximum(S[a], S[a + 1]),
                                       np.minimum(T[b], T[b + 1]), np.maximum(T[b], T[b + 1]),
                                       np.minimum(C[e], C[e + 1]), np.maximum(C[e], C[e + 1]))
                if not ok.any():
                    return ok
    return ok


def run(max_frontier=6_000_000):
    t0 = time.time()
    lo, hi = np.zeros((1, 4)), np.ones((1, 4))
    processed = by_bound = by_core = empty = 0
    while len(lo):
        processed += len(lo)
        infeasible = dn(dn(dn(lo[:, 0] + lo[:, 1]) + lo[:, 2]) + lo[:, 3]) > 1.0
        empty += int(infeasible.sum())
        lo, hi = lo[~infeasible], hi[~infeasible]
        kill = monotone_up(lo, hi) < TARGET_LO
        core = ~kill & concave_core_rig(lo, hi)
        by_bound += int(kill.sum()); by_core += int(core.sum())
        keep = ~(kill | core)
        lo, hi = lo[keep], hi[keep]
        if not len(lo):
            break
        w = hi - lo
        kk = w.argmax(-1)
        idx = np.arange(len(lo))
        m = (lo[idx, kk] + hi[idx, kk]) / 2
        lo2, hi2 = lo.copy(), hi.copy()
        hi[idx, kk] = m
        lo2[idx, kk] = m
        lo, hi = np.concatenate([lo, lo2]), np.concatenate([hi, hi2])
        if len(lo) > max_frontier:
            return dict(complete=False, processed=processed, frontier=len(lo))
    return dict(complete=True, processed=processed, closed_by_bound=by_bound, closed_by_concavity=by_core,
                empty=empty, seconds=round(time.time() - t0, 1))


def controls():
    """the rounding helpers enclose: A_up >= exact value on random points (mpmath, 50 digits); the rigorous Hessian test
    still rejects a non-concave cell and accepts the optimum's neighbourhood"""
    import mpmath as mp
    mp.mp.dps = 50
    rng = np.random.default_rng(11203)
    ok = True
    for _ in range(300):
        p, q, c = rng.uniform(0, 0.5, 3)
        exact = (mp.sqrt(mp.mpf(q) ** 2 + 4 * mp.mpf(c) * mp.mpf(p)) - mp.mpf(q)) / 2
        ok &= mp.mpf(float(A_up(np.array([p]), np.array([q]), np.array([c]))[0])) >= exact
    one = lambda v: np.array([v])
    rejects = not hessian_test_rig(one(0.399), one(0.401), one(0.049), one(0.051), one(0.299), one(0.301))[0]
    accepts = bool(concave_core_rig(np.full((1, 4), 2 / 15 - 0.005), np.full((1, 4), 2 / 15 + 0.005))[0])
    target_ok = mp.mpf(float(TARGET_LO)) < 4 / mp.sqrt(15)
    return dict(A_up_encloses=bool(ok), negative_control_rejected=bool(rejects),
                optimum_neighbourhood_accepted=accepts, target_below_4_over_sqrt15=bool(target_ok))


def main():
    res = dict(pass_id=11203, **controls(), **run())
    OUT.write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
