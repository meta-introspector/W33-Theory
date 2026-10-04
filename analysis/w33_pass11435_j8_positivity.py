"""Pass 11435: J_8 is positive on both components of J_6's blind spot (Passes 11369, 11419), with an exact law on the
pseudo-reflection component and high-precision values elsewhere.

CORRECTION of Passes 11418/11419 (displayed formula only).  For U = 1 + lam P, P = |psi><psi|, every trace is real:
tr(U^dag C U C^dag) = 3 - |lam|^2 (1 - x), x = |<psi|C psi>|^2 (and y with conj psi), so
        J_2t(U) = sum_k C(t2,k) 3^(t2-k) (-|lam|^2)^k E_k,   t2 = 2t,   E_k = avg (1-x)^k - avg (1-y)^k,
with E_k = sum_j C(k,j)(-1)^j Delta_j and Delta_j = avg x^j - avg y^j.  Since Delta_1..5 = 0 (Pass 11419):
        J_8(U) = |lam|^12 [252 Delta_6 - 24 |lam|^2 (7 Delta_6 - Delta_7) + |lam|^4 (28 Delta_6 - 8 Delta_7 + Delta_8)],
not |lam|^12 (252 Delta_6 - 24 |lam|^2 Delta_7 + ...) as displayed before; the limit 252 is unchanged.

THEOREM (exact, given the signs below).  On J_6's pseudo-reflection component {Delta_6 = 0}:
        J_8 = |lam|^14 [ 24 Delta_7 + |lam|^2 (Delta_8 - 8 Delta_7) ],   |lam|^2 in (0, 4],
which is linear in |lam|^2, hence positive for every beta != 0 iff Delta_7 > 0 and Delta_8 > 2 Delta_7.
COMPUTED: both hold at every sampled non-real-type ray of {h6 = 0}.

Elsewhere: J_8 at certified J_6-spurious points of the generic component, and a 50-digit evaluation (mpmath) of J_8
at the small-distance minimiser of Pass 11418's ratio test (whose double-precision value 3.6e-12 was at rounding).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11357_jarlskog_degree_by_dimension as P7  # noqa: E402
import w33_pass11369_j6_completeness as M  # noqa: E402
import w33_pass11418_j8_completeness as E  # noqa: E402

OUT = ROOT / "data" / "w33_pass11435_j8_positivity.json"


def deltas(psi, Cl, ks=(6, 7, 8)):
    psi = psi / np.linalg.norm(psi)
    x = np.abs(np.einsum('i,cij,j->c', psi.conj(), Cl, psi)) ** 2
    y = np.abs(np.einsum('i,cij,j->c', psi.conj(), Cl, psi.conj())) ** 2
    return [float((x ** k).mean() - (y ** k).mean()) for k in ks]


def pseudo_reflection_component(Cl, starts=300, seed=11435):
    rng = np.random.default_rng(seed)
    par = lambda v: v[:3] + 1j * v[3:]
    pts = []
    for _ in range(starts):
        r = minimize(lambda v: E.pseudo_reflection(3, par(v), 2.0, Cl), rng.normal(size=6), method='BFGS',
                     options=dict(gtol=1e-14))
        psi = par(r.x) / np.linalg.norm(par(r.x))
        ov = E.real_type_overlap(psi, Cl)
        if r.fun < 1e-12 and ov < 0.999:
            d6, d7, d8 = deltas(psi, Cl)
            errs = []
            for b in (2.0, 2.9):                                   # the law, where J_8 is far above rounding
                lam2 = abs(np.exp(1j * b) - 1) ** 2
                pred = lam2 ** 7 * (24 * d7 + lam2 * (d8 - 8 * d7))
                U = np.eye(3) + (np.exp(1j * b) - 1) * np.outer(psi, psi.conj())
                errs.append(abs(P7.J(U, 4, Cl) - pred) / pred)
            pts.append(dict(d6=d6, d7=d7, d8=d8, overlap=ov, law_rel_err=max(errs)))
    d7 = np.array([p["d7"] for p in pts])
    d8 = np.array([p["d8"] for p in pts])
    return dict(points=len(pts), max_abs_Delta6=max(abs(p["d6"]) for p in pts), min_Delta7=float(d7.min()),
                min_Delta8_over_Delta7=float((d8 / d7).min()), condition_holds_everywhere=bool((d7 > 0).all() and
                                                                                             (d8 > 2 * d7).all()),
                max_law_rel_err=max(p["law_rel_err"] for p in pts),
                min_J8_at_beta_pi=float(min(4 ** 7 * (24 * a + 4 * (b - 8 * a)) for a, b in zip(d7, d8))))


def generic_component(Cl, G, seeds=range(113690000, 113690300)):
    vals = []
    for s in seeds:
        U, J = M.descend(M.haar(3, np.random.default_rng(s)), 3, Cl, G)
        if abs(J) < 1e-12 and M.rev_distance(U, Cl) > 0.05:
            V = U / np.linalg.det(U) ** (1 / 3)
            ph = np.sort(np.angle(np.linalg.eigvals(V)))
            gap = float(np.min(np.diff(np.concatenate([ph, [ph[0] + 2 * np.pi]]))) / (2 * np.pi))
            if gap > 0.01:
                vals.append(P7.J(U, 4, Cl))
    return dict(points=len(vals), min_J8=float(min(vals)), median_J8=float(np.median(vals)))


def j8_mpmath(U, Cl, dps=50):
    import mpmath as mp
    mp.mp.dps = dps
    Um = mp.matrix(U.tolist())
    Ud = Um.H
    Ut = Um.T
    a = b = mp.mpf(0)
    for C in Cl:
        Cm = mp.matrix(C.tolist())
        X = Cm * Um * Cm.H
        Y = Cm * Ut * Cm.H
        a += abs(sum(Ud[i, j] * X[j, i] for i in range(3) for j in range(3))) ** 8
        b += abs(sum(Ud[i, j] * Y[j, i] for i in range(3) for j in range(3))) ** 8
    return (a - b) / len(Cl)


def small_distance_minimiser(Cl, delta=0.1, starts=40, seed=114350):
    """re-find a low-J8 point at dist >= delta (Pass 11418's ratio test) and evaluate J8 there in high precision.
    NOTE: U is stored as double; the mpmath value is of J8 at that double-precision U (exactly representable)."""
    G = M.hermitian_basis(3)
    best = None
    for i in range(starts):
        U0 = M.haar(3, np.random.default_rng(seed + i))

        def obj(th):
            V = U0 @ M.expi(np.einsum('a,aij->ij', th, G))
            d = M.rev_distance(V, Cl)
            return P7.J(V, 4, Cl) + 1e3 * max(0.0, delta - d) ** 2

        r = minimize(obj, np.zeros(8), method='Nelder-Mead', options=dict(maxiter=6000, xatol=1e-10, fatol=1e-16))
        V = U0 @ M.expi(np.einsum('a,aij->ij', r.x, G))
        d = M.rev_distance(V, Cl)
        if d >= 0.95 * delta and (best is None or P7.J(V, 4, Cl) < best[0]):
            best = (P7.J(V, 4, Cl), d, V)
    J8_hp = j8_mpmath(best[2], Cl)
    return dict(delta=delta, starts=starts, double_J8=float(best[0]), dist=float(best[1]), mpmath_J8=str(J8_hp),
                mpmath_positive=bool(J8_hp > 0))


def run():
    Cl = P7.clifford_group(3)
    G = M.hermitian_basis(3)
    res = dict(pass_id=11435)
    res["pseudo_reflection_component"] = pseudo_reflection_component(Cl)
    print(res["pseudo_reflection_component"], flush=True)
    res["generic_component"] = generic_component(Cl, G)
    print(res["generic_component"], flush=True)
    res["small_distance"] = [small_distance_minimiser(Cl, delta) for delta in (0.1, 0.2)]
    print(res["small_distance"], flush=True)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
