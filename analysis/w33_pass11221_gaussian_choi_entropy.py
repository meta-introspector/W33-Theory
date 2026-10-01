#!/usr/bin/env python3
"""Pass 11221: the Gaussian arrow is an entanglement rate -- the log-divergence of the Choi state's mode entanglement.

Pass 11219 proved, over the reals, A(S) = n - c(S) = min over mode decompositions of sum_k rank S[P_k-perp <- P_k], the
number of quadrature channels by which modes must leak into each other.  A first attempt to read this as state
entanglement (the vacuum) failed.  The right entropic object is the OPERATOR entanglement of the Gaussian unitary U_S,
measured on its Choi state with finite squeezing r:

    rho_S(r) = (U_S (x) 1) |TMSV(r)><TMSV(r)|^{(x) n} (U_S (x) 1)^dagger,

cut as (system mode k + its reference) | (all other modes and references).  Claim, checked here:

    S_k(r) = 2 rank S[P_k-perp <- P_k] * r + O(1)      (nats),

i.e. each coupling channel carries one infinitely squeezed EPR link across the cut (a TMSV has entropy 2r + O(1)).
Hence
    A(S) = (1/2) min over mode decompositions of  lim_{r -> oo} (d/dr) sum_k S_k(r),

the arrow of a Gaussian dynamics is half its minimal total operator-entanglement rate per unit squeezing.  It is zero
exactly for normal-mode-decomposable dynamics (elliptic and hyperbolic), and two per loxodromic quartet.

Why: in the limit the reference of mode k purifies its input, and S_k counts the independent quadratures of mode k's
output that depend on other modes' inputs (rank S[P_k <- P_k-perp], equal to rank S[P_k-perp <- P_k] for symplectic S);
each is one divergent symplectic eigenvalue nu ~ e^{2r} of the reduced covariance, contributing 2r.

Checks: every (S, split) from Pass 11219's generic and non-semisimple families, with the constructed optimal split and
two random splits each: |slope_k - 2 rank_k| < 0.05 at r = 5 -> 6 for every mode (float64; larger r amplifies rounding by e^{4r}); random splits have slope 4 per mode.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11219_gaussian_arrow as G  # noqa: E402

OUT = ROOT / "data" / "w33_pass11221_gaussian_choi_entropy.json"


def choi_cov(S, r):
    """covariance (x_all, p_all) of the Choi state on 2n modes (system 0..n-1, references n..2n-1)"""
    n = len(S) // 2
    c, s = np.cosh(2 * r), np.sinh(2 * r)
    Z = np.diag([1.0] * n + [-1.0] * n)
    sig = 0.5 * np.block([[c * np.eye(2 * n), s * Z], [s * Z, c * np.eye(2 * n)]])
    T = np.block([[S, np.zeros((2 * n, 2 * n))], [np.zeros((2 * n, 2 * n)), np.eye(2 * n)]])
    sig = T @ sig @ T.T
    perm = list(range(n)) + list(range(2 * n, 3 * n)) + list(range(n, 2 * n)) + list(range(3 * n, 4 * n))
    return sig[np.ix_(perm, perm)]


def entropy(sig, modes, N):
    rows = list(modes) + [N + i for i in modes]
    red = sig[np.ix_(rows, rows)]
    m = len(modes)
    Om = np.block([[np.zeros((m, m)), np.eye(m)], [-np.eye(m), np.zeros((m, m))]])
    nu = np.sort(np.abs(np.linalg.eigvals(1j * Om @ red)))[::2]
    tot = 0.0
    for v in nu:
        v = max(v, 0.5 + 1e-15)
        a, b = v + 0.5, v - 0.5
        tot += a * np.log(a) - (b * np.log(b) if b > 1e-300 else 0.0)
    return tot


def slopes(Sp, r1=5.0, r2=6.0):
    n = len(Sp) // 2
    s1, s2 = choi_cov(Sp, r1), choi_cov(Sp, r2)
    return [(entropy(s2, [k, n + k], 2 * n) - entropy(s1, [k, n + k], 2 * n)) / (r2 - r1) for k in range(n)]


def in_split_basis(S, planes):
    n = len(S) // 2
    J = G.Jm(n)
    es, fs = [], []
    for P in planes:
        e, f = P[:, 0], P[:, 1]
        es.append(e)
        fs.append(f / (e @ J @ f))
    M = np.column_stack(es + fs)
    assert np.allclose(M.T @ J @ M, J, atol=1e-7)
    return np.linalg.inv(M) @ S @ M, M


def random_split(n, rng):
    M = G.rand_symplectic(n, rng, scale=0.8)
    return [np.stack([M[:, k], M[:, n + k]], 1) for k in range(n)]


def build(e, h, l, rng):
    n = e + h + 2 * l
    blocks = [G.elliptic(rng.uniform(0.2, 3.0)) for _ in range(e)] + \
             [G.hyperbolic(rng.uniform(1.2, 2.5) * rng.choice([-1, 1])) for _ in range(h)] + \
             [G.loxodromic(rng.uniform(1.1, 2.0), rng.uniform(0.3, 2.8)) for _ in range(l)]
    S0, idx = G.embed(blocks, n)
    planes0 = []
    for (ii, k), (Sb, _) in zip(idx, blocks):
        if k == 1:
            P = np.zeros((2 * n, 2))
            P[ii[0], 0] = 1
            P[ii[1], 1] = 1
            planes0.append(P)
        else:
            sub, _ = G.tz_split_quartet(Sb)
            for Pl in sub:
                P = np.zeros((2 * n, 2))
                P[ii, :] = Pl
                planes0.append(P)
    M = G.rand_symplectic(n, rng, scale=0.3)
    return M @ S0 @ np.linalg.inv(M), [M @ P for P in planes0]


def check(S, planes):
    ranks = G.coupling_ranks(S, planes)
    Sp, _ = in_split_basis(S, planes)
    sl = slopes(Sp)
    return ranks, sl, all(abs(s - 2 * r) < 0.05 for s, r in zip(sl, ranks))


def run(seed=11221):
    rng = np.random.default_rng(seed)
    rows = []
    for (e, h, l) in [(2, 0, 0), (0, 2, 0), (1, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1), (0, 0, 2), (1, 1, 1),
                      (0, 0, 3)]:
        for _ in range(2):
            S, planes = build(e, h, l, rng)
            r_opt, s_opt, ok_opt = check(S, planes)
            rand = [check(S, random_split(len(S) // 2, rng)) for _ in range(2)]
            rows.append(dict(e=e, h=h, l=l, n=len(S) // 2, A=2 * l, optimal_ranks=r_opt,
                             optimal_slopes=[round(x, 4) for x in s_opt], optimal_ok=ok_opt,
                             optimal_total_over_2=round(sum(s_opt) / 2, 3),
                             random_slopes=[[round(x, 3) for x in rs[1]] for rs in rand],
                             random_ok=all(rs[2] for rs in rand)))
    ns = []
    for name, v in G.nonsemisimple_cases(np.random.default_rng(seed)).items():
        ns.append(dict(case=name, total_rank=v["total_rank"]))
    # re-run the non-semisimple splits through the entropy check
    rng2 = np.random.default_rng(seed + 1)
    ns_rows = []
    J = G.Jm(2)
    Hhh = np.zeros((4, 4))
    Hhh[1, 2] = Hhh[2, 1] = 1.0
    Hhh[0, 3] = Hhh[3, 0] = -1.0
    Hhh[0, 0] = Hhh[1, 1] = 1.0
    from scipy.linalg import expm
    for name, S0 in {"elliptic Jordan pair (Krein collision)": expm(0.7 * J @ Hhh)}.items():
        M = G.rand_symplectic(2, rng2, scale=0.3)
        S = M @ S0 @ np.linalg.inv(M)
        planes, _ = G.find_bilagrangian(S, rng2)
        ranks, sl, ok = check(S, planes)
        ns_rows.append(dict(case=name, ranks=ranks, slopes=[round(x, 4) for x in sl], ok=ok))
    res = dict(pass_id=11221, cases=rows, nonsemisimple=ns_rows,
               slope_equals_twice_rank=all(r["optimal_ok"] and r["random_ok"] for r in rows)
               and all(r["ok"] for r in ns_rows),
               optimal_half_total_equals_A=all(abs(r["optimal_total_over_2"] - r["A"]) < 0.05 for r in rows))
    return res


def main():
    res = run()
    OUT.write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps({k: v for k, v in res.items() if k != "cases"}, indent=1, default=float))
    for r in res["cases"]:
        print(r["e"], r["h"], r["l"], "ranks", r["optimal_ranks"], "slopes", r["optimal_slopes"], "random",
              r["random_slopes"][0])


if __name__ == "__main__":
    main()
