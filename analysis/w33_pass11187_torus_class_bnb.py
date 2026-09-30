#!/usr/bin/env python3
"""Pass 11187: 4/sqrt15 is the maximum of N_AB + N_AC on the torus-symmetric class -- a certified branch and bound.

The class.  States fixed by the torus (D, conj D, conj D), D = diag(1, e^{i a}, e^{i b}) -- the maximal torus of the U(2)
symmetry of the optimum psi* (Pass 11162): psi = alpha|000> + sum_{a=1,2} (beta_a |aa0> + gamma_a |a0a>), phases removable.
With s_a = beta_a^2, t_a = gamma_a^2, c = alpha^2 = 1 - s1 - s2 - t1 - t2 the partial transposes split into 2 x 2 blocks:
    N_AB = sum_a (sqrt(t_a^2 + 4 c s_a) - t_a)/2 + sqrt(s1 s2),     N_AC = sum_a (sqrt(s_a^2 + 4 c t_a) - s_a)/2 + sqrt(t1 t2)
(checked against direct negativity).  Pass 11167 settled the 2-parameter U(2) class (s1 = s2, t1 = t2); this pass the
4-parameter torus class.
Proof structure (all bounds evaluated in floating point with an explicit safety margin 1e-9 over the target):
  (1) AM-GM: sqrt(s1 s2) <= (s1+s2)/2, sqrt(t1 t2) <= (t1+t2)/2, so f <= Psi = (1/2) sum_a phi(s_a, t_a; c),
      phi(s, t; c) = sqrt(t^2 + 4cs) + sqrt(s^2 + 4ct), with equality on the symmetric slice.
  (2) CORE: on a box where phi(., .; c) is concave in (s, t) for every c in the box's range (interval Hessian test:
      H11 < 0, H22 < 0, H12^2 < H11 H22 on the rectangle hull of both (s_a, t_a) ranges), Jensen gives
      Psi(s1,t1,s2,t2) <= Psi(sym) = f(sym) <= 4/sqrt15 (Pass 11167, exact on the symmetric slice).
  (3) Elsewhere a box is discarded when an upper bound for f is below 4/sqrt15: the minimum of the monotone corner bound
      (each term is monotone in each variable, c decreasing in all) and the first-order centred form
      f(mid) + sum_i max|df/dx_i| * halfwidth_i (interval gradient).
The branch and bound terminates with every box closed by (2) or (3): a computer-assisted proof that the torus-class
maximum is 4/sqrt15.  Scope: this is still a restricted class (the torus-invariant states), not the global problem.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "w33_pass11187_torus_class_bnb.json"
TARGET = 4 / np.sqrt(15)
MARGIN = 1e-9


def f_point(x):
    s1, s2, t1, t2 = x.T
    c = np.clip(1 - x.sum(-1), 0, None)
    A = lambda p, q: (np.sqrt(q * q + 4 * c * p) - q) / 2
    return A(s1, t1) + A(s2, t2) + np.sqrt(s1 * s2) + A(t1, s1) + A(t2, s2) + np.sqrt(t1 * t2)


def monotone_bound(lo, hi):
    """A(p, q; c) = (sqrt(q^2 + 4cp) - q)/2 is increasing in p and c, decreasing in q"""
    c_hi = np.clip(1 - lo.sum(-1), 0, None)
    A = lambda p, q: (np.sqrt(q * q + 4 * c_hi * p) - q) / 2
    s1l, s2l, t1l, t2l = lo.T
    s1h, s2h, t1h, t2h = hi.T
    return (A(s1h, t1l) + A(s2h, t2l) + np.sqrt(s1h * s2h) + A(t1h, s1l) + A(t2h, s2l) + np.sqrt(t1h * t2h))


def grad_bound(lo, hi):
    """interval bound on max |df/dx_i| over the box (inf where a sqrt(x) derivative blows up)"""
    c_lo = np.clip(1 - hi.sum(-1), 0, None)
    c_hi = np.clip(1 - lo.sum(-1), 0, None)

    def pieces(p_lo, p_hi, q_lo, q_hi):
        G_lo = np.sqrt(q_lo ** 2 + 4 * c_lo * p_lo)
        G_hi = np.sqrt(q_hi ** 2 + 4 * c_hi * p_hi)
        G_lo = np.maximum(G_lo, 1e-300)
        Ap = (c_lo / G_hi, c_hi / G_lo)                       # dA/dp = c/G
        Aq = ((q_lo / G_hi - 1) / 2, (q_hi / G_lo - 1) / 2)    # dA/dq = (q/G - 1)/2
        Ac = (p_lo / G_hi, p_hi / G_lo)                       # dA/dc = p/G
        return Ap, Aq, Ac
    s1l, s2l, t1l, t2l = lo.T
    s1h, s2h, t1h, t2h = hi.T
    P1 = pieces(s1l, s1h, t1l, t1h)       # A(s1, t1)
    P2 = pieces(s2l, s2h, t2l, t2h)       # A(s2, t2)
    Q1 = pieces(t1l, t1h, s1l, s1h)       # A(t1, s1)
    Q2 = pieces(t2l, t2h, s2l, s2h)       # A(t2, s2)
    Ac_lo = P1[2][0] + P2[2][0] + Q1[2][0] + Q2[2][0]
    Ac_hi = P1[2][1] + P2[2][1] + Q1[2][1] + Q2[2][1]
    with np.errstate(divide='ignore', invalid='ignore'):
        r = lambda a_lo, a_hi, b_lo, b_hi: (0.5 * np.sqrt(b_lo / np.maximum(a_hi, 1e-300)), 0.5 * np.sqrt(np.where(a_lo > 0, b_hi / np.where(a_lo > 0, a_lo, 1), np.inf)))
        rs1, rs2 = r(s1l, s1h, s2l, s2h), r(s2l, s2h, s1l, s1h)
        rt1, rt2 = r(t1l, t1h, t2l, t2h), r(t2l, t2h, t1l, t1h)
    # df/ds1 = dA(s1,t1)/dp + dA(t1,s1)/dq + (1/2) sqrt(s2/s1) - sum dA/dc
    d = []
    for Pp, Qq, rr in ((P1, Q1, rs1), (P2, Q2, rs2)):
        lo_ = Pp[0][0] + Qq[1][0] + rr[0] - Ac_hi
        hi_ = Pp[0][1] + Qq[1][1] + rr[1] - Ac_lo
        d.append(np.maximum(np.abs(lo_), np.abs(hi_)))
    for Qp, Pq, rr in ((Q1, P1, rt1), (Q2, P2, rt2)):
        lo_ = Qp[0][0] + Pq[1][0] + rr[0] - Ac_hi
        hi_ = Qp[0][1] + Pq[1][1] + rr[1] - Ac_lo
        d.append(np.maximum(np.abs(lo_), np.abs(hi_)))
    return np.stack([d[0], d[1], d[2], d[3]], -1)            # order s1, s2, t1, t2


def centred_bound(lo, hi):
    mid = (lo + hi) / 2
    hw = (hi - lo) / 2
    val = f_point(np.clip(mid, 0, None))
    with np.errstate(all='ignore'):
        g = grad_bound(lo, hi)
        ub = val + (np.where(np.isfinite(g), g, np.inf) * hw).sum(-1)
    # the centred form is used only on boxes lying entirely inside the simplex (c = 1 - sum > 0 on the whole box),
    # where f is smooth and the mean value theorem applies with the gradient formulas above
    inside = (hi.sum(-1) < 1) & (lo > 0).all(-1)
    return np.where(inside & np.isfinite(ub), ub, np.inf)


def _hessian_test(s_lo, s_hi, t_lo, t_hi, c_lo, c_hi):
    """interval test that phi(s,t;c) is concave on every cell [s]x[t] for every c in [c]"""
    return _hessian_test_raw(s_lo, s_hi, t_lo, t_hi, c_lo, c_hi)


def _hessian_test_raw(s_lo, s_hi, t_lo, t_hi, c_lo, c_hi):
    ok = (c_lo > 0) & (s_lo > 0) & (t_lo > 0)
    s_lo, t_lo, c_lo = np.where(ok, s_lo, 1.0), np.where(ok, t_lo, 1.0), np.where(ok, c_lo, 1.0)
    G1_lo = np.sqrt(t_lo ** 2 + 4 * c_lo * s_lo); G1_hi = np.sqrt(t_hi ** 2 + 4 * c_hi * s_hi)
    G2_lo = np.sqrt(s_lo ** 2 + 4 * c_lo * t_lo); G2_hi = np.sqrt(s_hi ** 2 + 4 * c_hi * t_hi)
    # H1 = [[-4c^2, -2ct], [-2ct, 4cs]]/G1^3,  H2 = [[4ct, -2cs], [-2cs, -4c^2]]/G2^3
    H11_hi = -4 * c_lo ** 2 / G1_hi ** 3 + 4 * c_hi * t_hi / G2_lo ** 3
    H22_hi = 4 * c_hi * s_hi / G1_lo ** 3 - 4 * c_lo ** 2 / G2_hi ** 3
    H12_abs = 2 * c_hi * t_hi / G1_lo ** 3 + 2 * c_hi * s_hi / G2_lo ** 3
    return ok & (H11_hi < 0) & (H22_hi < 0) & (H12_abs ** 2 < H11_hi * H22_hi)


def concave_core(lo, hi, k=6):
    """phi(., .; c) concave on the rectangle hull of both (s_a, t_a) boxes, for every c in the box's range; the hull is
    cut into k x k x k cells, each tested with interval bounds (all must pass)"""
    s_lo = np.minimum(lo[:, 0], lo[:, 1]); s_hi = np.maximum(hi[:, 0], hi[:, 1])
    t_lo = np.minimum(lo[:, 2], lo[:, 3]); t_hi = np.maximum(hi[:, 2], hi[:, 3])
    c_lo = np.clip(1 - hi.sum(-1), 0, None); c_hi = np.clip(1 - lo.sum(-1), 0, None)
    ok = np.ones(len(lo), bool)
    fr = np.linspace(0, 1, k + 1)
    for a in range(k):
        for b in range(k):
            for e in range(k):
                sl = s_lo + (s_hi - s_lo) * fr[a]; sh = s_lo + (s_hi - s_lo) * fr[a + 1]
                tl = t_lo + (t_hi - t_lo) * fr[b]; th = t_lo + (t_hi - t_lo) * fr[b + 1]
                cl = c_lo + (c_hi - c_lo) * fr[e]; ch = c_lo + (c_hi - c_lo) * fr[e + 1]
                ok &= _hessian_test(sl, sh, tl, th, cl, ch)
                if not ok.any():
                    return ok
    return ok


def run(max_boxes=40_000_000):
    t0 = time.time()
    lo = np.zeros((1, 4))
    hi = np.full((1, 4), 0.5)          # s_a, t_a <= 1/2 automatically (2 s_a <= s_a + ... <= 1) -- use the simplex below
    hi = np.full((1, 4), 1.0)
    processed = discarded_bound = discarded_core = 0
    empty = 0
    while len(lo):
        processed += len(lo)
        if processed > max_boxes:
            return dict(complete=False, processed=processed)
        feasible = lo.sum(-1) <= 1
        empty += int((~feasible).sum())
        lo, hi = lo[feasible], hi[feasible]
        ub = np.minimum(monotone_bound(lo, hi), centred_bound(lo, hi))
        kill = ub < TARGET - MARGIN
        core = ~kill & concave_core(lo, hi)
        discarded_bound += int(kill.sum())
        discarded_core += int(core.sum())
        keep = ~(kill | core)
        lo, hi = lo[keep], hi[keep]
        if not len(lo):
            break
        # split the longest edge
        w = hi - lo
        k = w.argmax(-1)
        m = (lo[np.arange(len(lo)), k] + hi[np.arange(len(lo)), k]) / 2
        lo2, hi2 = lo.copy(), hi.copy()
        hi[np.arange(len(lo)), k] = m
        lo2[np.arange(len(lo)), k] = m
        lo, hi = np.concatenate([lo, lo2]), np.concatenate([hi, hi2])
        if len(lo) > 4_000_000:
            return dict(complete=False, processed=processed, frontier=len(lo))
    return dict(complete=True, processed=processed, closed_by_bound=discarded_bound, closed_by_concavity=discarded_core,
                empty=empty, seconds=round(time.time() - t0, 1))


def hessian_checks(seed=5):
    """(a) the analytic Hessian of phi(s,t;c) matches finite differences; (b) the interval test is not vacuous: it
    rejects cells where phi is not concave (negative control) and accepts the optimum's neighbourhood"""
    rng = np.random.default_rng(seed)
    phi = lambda s, t, c: np.sqrt(t * t + 4 * c * s) + np.sqrt(s * s + 4 * c * t)

    def H_formula(s, t, c):
        G1, G2 = np.sqrt(t * t + 4 * c * s), np.sqrt(s * s + 4 * c * t)
        H1 = np.array([[-4 * c * c, -2 * c * t], [-2 * c * t, 4 * c * s]]) / G1 ** 3
        H2 = np.array([[4 * c * t, -2 * c * s], [-2 * c * s, -4 * c * c]]) / G2 ** 3
        return H1 + H2
    maxerr, h = 0.0, 1e-5
    for _ in range(50):
        s, t, c = rng.uniform(0.02, 0.3, 3)
        num = np.zeros((2, 2))
        num[0, 0] = (phi(s + h, t, c) - 2 * phi(s, t, c) + phi(s - h, t, c)) / h ** 2
        num[1, 1] = (phi(s, t + h, c) - 2 * phi(s, t, c) + phi(s, t - h, c)) / h ** 2
        num[0, 1] = num[1, 0] = (phi(s + h, t + h, c) - phi(s + h, t - h, c) - phi(s - h, t + h, c) + phi(s - h, t - h, c)) / (4 * h * h)
        maxerr = max(maxerr, float(np.abs(num - H_formula(s, t, c)).max()))
    one = lambda v: np.array([v])
    non_concave_point = np.linalg.eigvalsh(H_formula(0.4, 0.05, 0.3)).max() > 0
    rejects = not _hessian_test(one(0.399), one(0.401), one(0.049), one(0.051), one(0.299), one(0.301))[0]
    accepts = bool(concave_core(np.full((1, 4), 2 / 15 - 0.005), np.full((1, 4), 2 / 15 + 0.005))[0])
    return dict(hessian_fd_max_error=maxerr, negative_control_not_concave=bool(non_concave_point),
                negative_control_rejected=bool(rejects), optimum_neighbourhood_accepted=accepts)


def summarize():
    rng = np.random.default_rng(11187)
    # sanity: the closed form equals the direct negativity (checked in Pass 11167 for the symmetric slice; here random)
    import sys
    sys.path.insert(0, str(ROOT / "analysis"))
    import w33_pass11167_polygamy_global as Gp
    ok = True
    for _ in range(20):
        x = rng.dirichlet(np.ones(5))[:4]
        s1, s2, t1, t2 = x
        psi = np.zeros(27)
        psi[0] = np.sqrt(1 - x.sum())
        psi[1 * 9 + 1 * 3], psi[2 * 9 + 2 * 3], psi[1 * 9 + 1], psi[2 * 9 + 2] = np.sqrt([s1, s2, t1, t2])
        ok &= abs(sum(Gp.negs(psi)) - f_point(x[None])[0]) < 1e-12
    res = dict(pass_id=11187, closed_form_checked=bool(ok), target=TARGET, margin=MARGIN, **hessian_checks(), **run())
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
