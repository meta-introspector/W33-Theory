"""Pass 11418: is the degree-8 witness J_8 complete on PU(3)?  (Pass 11369: J_6 is not; conjecture t* = 4.)

Three tests, each with J_6 as the positive control (it MUST fail them, since Pass 11369 found its spurious family):
  (a) RATIO.  R_t(U) = J_2t(U) / dist(U)^2, minimised over {dist(U) >= delta} (penalty method, random Haar starts),
      dist = distance to the reversible set (Pass 11369).  The reported quantity that matters is J itself at the
      minimiser: an incomplete witness reaches J = 0 far from the reversible set.  CORRECTION (found here): the RATIO
      cannot be bounded below even for a complete witness -- near the identity every J_2t vanishes like dist^12 along
      pseudo-reflections (identity_asymptotics, --asym), so min R_t -> 0 as delta -> 0 regardless of completeness.
  (b) SEEDED.  J_8 descents started AT certified spurious zeros of J_6 (both populations): do they end on the
      reversible set (complete) or at a non-reversible zero of J_8?
  (c) PSEUDO-REFLECTIONS.  U = 1 + (e^{i b} - 1)|psi><psi| (the degenerate-spectrum locus that carries one J_6
      component, Pass 11419).  Test whether J_2t factorises as f(b) g(psi) and minimise g over CP^2: J_6's g vanishes
      on a real hypersurface of non-real-type rays; does J_8's?
"""

from __future__ import annotations

import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11357_jarlskog_degree_by_dimension as P7  # noqa: E402
import w33_pass11369_j6_completeness as M  # noqa: E402

OUT = ROOT / "data" / "w33_pass11418_j8_completeness.json"
_S = {}


def _init():
    _S["Cl"] = P7.clifford_group(3)
    _S["G"] = M.hermitian_basis(3)


def _ratio_job(args):
    t, delta, seed = args
    Cl, G = _S["Cl"], _S["G"]
    U0 = M.haar(3, np.random.default_rng(seed))

    def U(th):
        return U0 @ M.expi(np.einsum('a,aij->ij', th, G))

    def obj(th):
        V = U(th)
        d = M.rev_distance(V, Cl)
        pen = 1e3 * max(0.0, delta - d) ** 2
        return P7.J(V, t, Cl) / max(d, 1e-6) ** 2 + pen

    best = None
    for _ in range(3):                                    # restarts from the current point
        r = minimize(obj, np.zeros(8) if best is None else best.x, method='Nelder-Mead',
                     options=dict(maxiter=4000, xatol=1e-9, fatol=1e-14))
        best = r if best is None or r.fun < best.fun else best
    V = U(best.x)
    d = M.rev_distance(V, Cl)
    return dict(t=t, delta=delta, ratio=float(P7.J(V, t, Cl) / d ** 2), dist=float(d), J=float(P7.J(V, t, Cl)))


def ratio_test(ts=(3, 4), deltas=(0.1, 0.3), starts=200, nproc=8):
    out = {}
    with Pool(nproc, initializer=_init) as pool:
        for t in ts:
            for delta in deltas:
                rows = pool.map(_ratio_job, [(t, delta, 114180000 + 1000 * t + int(100 * delta) * 7 + i)
                                             for i in range(starts)])
                ok = [r for r in rows if r["dist"] >= 0.9 * delta]
                best = min(ok, key=lambda r: r["ratio"])
                out[f"J{2 * t}, dist>={delta}"] = dict(starts=starts, feasible=len(ok), min_ratio=best["ratio"],
                                                     at_dist=best["dist"], J_there=best["J"])
                print(f"J{2 * t} delta {delta}", out[f"J{2 * t}, dist>={delta}"], flush=True)
    return out


def _seeded_job(seed):
    Cl, G = _S["Cl"], _S["G"]
    U, J6 = M.descend(M.haar(3, np.random.default_rng(seed)), 3, Cl, G)
    d6 = M.rev_distance(U, Cl)
    if not (abs(J6) < 1e-12 and d6 > 0.05 and P7.J(U, 4, Cl) > 1e-6):
        return None                                          # not a certified J_6-spurious point
    gap = float(np.min(np.diff(np.sort(np.angle(np.linalg.eigvals(U / np.linalg.det(U) ** (1 / 3)))))))
    V, J8 = M.descend(U, 4, Cl, G)
    return dict(seed=seed, start_dist=d6, start_J8=float(P7.J(U, 4, Cl)), near_degenerate=gap < 0.02 * 2 * np.pi,
                end_J8=float(J8), end_dist=M.rev_distance(V, Cl), end_J10=float(P7.J(V, 5, Cl)))


def seeded(n=2000, nproc=8):
    with Pool(nproc, initializer=_init) as pool:
        rows = [r for r in pool.map(_seeded_job, range(113690000, 113690000 + n), chunksize=8) if r]
    zeros = [r for r in rows if r["end_J8"] < 1e-12]
    spurious = [r for r in zeros if r["end_dist"] > 1e-3 and r["end_J10"] > 1e-6]
    return dict(seeds=n, j6_spurious_starts=len(rows),
                near_degenerate_starts=sum(r["near_degenerate"] for r in rows),
                j8_zeros_reached=len(zeros), certified_j8_spurious=len(spurious),
                max_end_dist_at_zeros=max((r["end_dist"] for r in zeros), default=None),
                min_start_J8=min((r["start_J8"] for r in rows), default=None))


def pseudo_reflection(t, psi, b, Cl):
    psi = psi / np.linalg.norm(psi)
    return P7.J(np.eye(3) + (np.exp(1j * b) - 1) * np.outer(psi, psi.conj()), t, Cl)


def real_type_overlap(psi, Cl):
    """max_C |<conj psi | C psi>|: equals 1 iff conj(psi) is in the Clifford orbit of psi (a reversible ray)"""
    psi = psi / np.linalg.norm(psi)
    return float(np.abs(np.einsum('i,cij,j->c', psi, Cl, psi)).max())


def pseudo_reflection_test(starts=150, seed=11418):
    Cl = P7.clifford_group(3)
    rng = np.random.default_rng(seed)
    out = {}
    for t in (3, 4):
        psis = [rng.normal(size=3) + 1j * rng.normal(size=3) for _ in range(6)]
        bs = [1.2, 1.9, 2.6, 3.0]
        T = np.array([[pseudo_reflection(t, p, b, Cl) for b in bs] for p in psis])
        sv = np.linalg.svd(T, compute_uv=False)
        fact = float(sv[1] / sv[0])                       # rank 1 (second singular value ~ 0) <=> J = f(b) g(psi)
        par = lambda x: x[:3] + 1j * x[3:]
        zeros = []
        for _ in range(starts):
            r = minimize(lambda x: pseudo_reflection(t, par(x), 2.0, Cl), rng.normal(size=6), method='BFGS',
                         options=dict(gtol=1e-14))
            zeros.append((float(r.fun), real_type_overlap(par(r.x), Cl)))
        z = [o for f, o in zeros if f < 1e-12]
        out[f"J{2 * t}"] = dict(second_singular_value_relative=fact, minimisations=starts, zeros=len(z),
                                zeros_off_real_type=sum(o < 0.999 for o in z),
                                min_overlap_at_zeros=min(z) if z else None,
                                smallest_minimum=min(f for f, _ in zeros))
        print(f"J{2 * t}", out[f"J{2 * t}"], flush=True)
    return out


def identity_asymptotics(rays=3, betas=(1.0, 0.6, 0.4, 0.25), seed=4):
    """why the ratio test cannot bound J_8 below near the identity: on pseudo-reflections (Pass 11419)
    J_8(1 + lam P) = |lam|^12 (252 Delta_6 - 24 |lam|^2 Delta_7 + |lam|^4 Delta_8), so J_8 / (|lam|^12 Delta_6) -> 252,
    while the distance to the reversible set is only O(|lam|): every witness vanishes like dist^12 there"""
    import w33_pass11419_spurious_family as S
    Cl = P7.clifford_group(3)
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(rays):
        psi = rng.normal(size=3) + 1j * rng.normal(size=3)
        psi /= np.linalg.norm(psi)
        d6 = S.ray_witness(psi, Cl)
        rows.append([float(P7.J(np.eye(3) + (np.exp(1j * b) - 1) * np.outer(psi, psi.conj()), 4, Cl)
                           / (abs(np.exp(1j * b) - 1) ** 12 * d6)) for b in betas])
    return dict(betas=list(betas), J8_over_lam12_Delta6=rows, limit_predicted=252)


def run():
    res = dict(pass_id=11418)
    res["pseudo_reflections"] = pseudo_reflection_test()
    res["seeded_from_J6_spurious"] = seeded()
    print(res["seeded_from_J6_spurious"], flush=True)
    res["ratio_test"] = ratio_test()
    return res


def main():
    if "--asym" in sys.argv:
        res = json.load(open(OUT))
        res["identity_asymptotics"] = identity_asymptotics()
        print(res["identity_asymptotics"])
        json.dump(res, open(OUT, "w"), indent=1, default=str)
        return
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
