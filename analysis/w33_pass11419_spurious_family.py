"""Pass 11419: the anatomy of J_6's spurious zero set on PU(3) (Pass 11369).

Pass 11369 found non-reversible zeros of J_6 in two populations: a generic-spectrum one (refined point: smooth,
4-dimensional) and a near-degenerate-spectrum one (35% of descents).  Here:
  (a) PSEUDO-REFLECTIONS.  Up to phase, a unitary with a repeated eigenvalue is P_{psi,b} = 1 + (e^{ib} - 1)|psi><psi|,
      and P^T = P_{conj psi, b}: it is reversible iff conj(psi) is in the Clifford orbit of psi (a 2-dimensional set
      of 'real-type' rays in CP^2).  Computed: J_6(P_{psi,b}) = f(b) g(psi) EXACTLY (ratios over b independent of psi),
      and the zero set of g contains non-real-type rays and has codimension 1 in CP^2 (Hessian rank 1 at its zeros).
      So {P_{psi,b} : g(psi) = 0} is a 4-dimensional spurious COMPONENT (3 + 1), inside the 5-dimensional
      degenerate-spectrum locus.
  (b) The near-degenerate descents land on it: projecting them to the nearest pseudo-reflection (merge the two close
      eigenvalues) keeps J_6 at rounding level while the distance to the reversible set stays finite.
  (c) Generic population: distinct component (refined point of Pass 11369: generic spectrum, rank dF = 4).
  (d) Control: the natural cubic Clifford invariant S3(psi) = sum_p <psi|W(p)|psi>^3 is REAL for every ray (so
      tau-even), hence is not the odd invariant behind g; g's identification is open.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11252_exact_reversibility as R  # noqa: E402
import w33_pass11357_jarlskog_degree_by_dimension as P7  # noqa: E402
import w33_pass11369_j6_completeness as M  # noqa: E402
import w33_pass11418_j8_completeness as E  # noqa: E402

OUT = ROOT / "data" / "w33_pass11419_spurious_family.json"


def factorisation(Cl, rng, rays=12, betas=(1.2, 1.9, 2.6, 3.0)):
    """J_6(P_{psi,b}) = f(b) g(psi)  <=>  the matrix T[ray, b] has rank 1: report its normalised singular values (small
    b are avoided: there J_6 ~ b^12-small and rounding dominates any ratio)"""
    psis = [rng.normal(size=3) + 1j * rng.normal(size=3) for _ in range(rays)]
    T = np.array([[E.pseudo_reflection(3, psi, b, Cl) for b in betas] for psi in psis])
    sv = np.linalg.svd(T, compute_uv=False)
    prof = T[0] / T[0, -1]
    law = np.abs(1 - np.exp(1j * np.array(betas))) ** 12
    law = law / law[-1]
    return dict(rays=rays, betas=list(betas), singular_values_relative=[float(x) for x in sv / sv[0]],
                beta_profile_relative=[float(x) for x in prof],
                max_rel_dev_from_abs_1_minus_e_ib_pow12=float(np.max(np.abs(prof / law - 1))))


def ray_witness(psi, Cl):
    """g(psi) = J_6 of the projector |psi><psi| = avg_C |<psi|C psi>|^12 - avg_C |<psi|C conj(psi)>|^12"""
    psi = psi / np.linalg.norm(psi)
    a = np.abs(np.einsum('i,cij,j->c', psi.conj(), Cl, psi)) ** 12
    b = np.abs(np.einsum('i,cij,j->c', psi.conj(), Cl, psi.conj())) ** 12
    return float(a.mean() - b.mean())


def projector_identity(Cl, rng, rays=10, betas=(0.7, 1.6, 2.4, 3.1)):
    """J_6(1 + lam P) = |lam|^12 J_6(P) for P = |psi><psi|, lam = e^{ib} - 1 (and J_6(P) = ray_witness(psi))"""
    dev = dev_w = 0.0
    rel_raw = 0.0
    jps = []
    for _ in range(rays):
        psi = rng.normal(size=3) + 1j * rng.normal(size=3)
        psi /= np.linalg.norm(psi)
        P = np.outer(psi, psi.conj())
        jp = P7.J(P, 3, Cl)
        jps.append(jp)
        dev_w = max(dev_w, abs(jp - ray_witness(psi, Cl)) / abs(jp))
        for b in betas:
            lam = np.exp(1j * b) - 1
            e = abs(P7.J(np.eye(3) + lam * P, 3, Cl) / abs(lam) ** 12 - jp)
            dev = max(dev, e)
            rel_raw = max(rel_raw, e / abs(jp))
    return dict(rays=rays, betas=list(betas),
                max_abs_dev_of_J6U_over_abs_lam12_minus_J6P_relative_to_max_J6P=dev / max(jps),
                max_raw_relative_dev=rel_raw, note="raw relative deviation is rounding-dominated where J6(P) and "
                "|lam|^12 are both small (J6 is a difference of averages with ~1e-13 absolute error)",
                max_rel_dev_J6P_equals_ray_witness=dev_w)


def odd_ray_invariants(kmax=7):
    """Clifford-invariant polynomials of bidegree (k,k) on C^3 (= Clifford invariants in Sym^k V (x) Sym^k Vbar, Pass
    11372's Space) and the action of tau: f(psi) -> f(conj psi), which transposes the coefficient tensor.  Counts the
    tau-odd ones (tau-eigenvalue -1)."""
    import w33_pass11372_irrep_contraction_rates as IR
    reps = IR.coset_reps()
    out = {}
    for k in range(1, kmax + 1):
        B = IR.Space(k, k, reps).clifford_invariants()
        ev = np.linalg.eigvals(np.einsum('aij,bji->ab', B.conj(), B)).real
        out[str(k)] = dict(invariants=len(B), tau_odd=int((ev < -0.5).sum()), tau_even=int((ev > 0.5).sum()))
    return out


def moments(Cl, rng, rays=8, kmax=7):
    """Delta_k(psi) = avg_C |<psi|C psi>|^(2k) - avg_C |<psi|C conj psi>|^(2k): zero for k <= 5, nonzero from k = 6"""
    mx = {}
    for _ in range(rays):
        psi = rng.normal(size=3) + 1j * rng.normal(size=3)
        psi /= np.linalg.norm(psi)
        x = np.abs(np.einsum('i,cij,j->c', psi.conj(), Cl, psi)) ** 2
        y = np.abs(np.einsum('i,cij,j->c', psi.conj(), Cl, psi.conj())) ** 2
        for k in range(1, kmax + 1):
            mx[str(k)] = max(mx.get(str(k), 0.0), float(abs((x ** k).mean() - (y ** k).mean())))
    return mx


def g_zero_set(Cl, rng, starts=200):
    par = lambda x: x[:3] + 1j * x[3:]
    zeros = []
    for _ in range(starts):
        r = minimize(lambda x: E.pseudo_reflection(3, par(x), 2.0, Cl), rng.normal(size=6), method='BFGS',
                     options=dict(gtol=1e-14))
        psi = par(r.x) / np.linalg.norm(par(r.x))
        if r.fun < 1e-12:
            zeros.append((psi, E.real_type_overlap(psi, Cl)))
    off = [z for z in zeros if z[1] < 0.999]
    ranks = []
    for psi, _ in off[:5]:
        # Hessian of g on the 4-dim tangent space of CP^2 at psi (orthonormal real basis of psi-perp)
        Bc = np.linalg.svd(np.eye(3) - np.outer(psi, psi.conj()))[0][:, :2]
        B = [Bc[:, 0], 1j * Bc[:, 0], Bc[:, 1], 1j * Bc[:, 1]]
        h = 1e-4
        f = lambda d: E.pseudo_reflection(3, psi + d, 2.0, Cl)
        H = np.array([[(f(h * u + h * v) - f(h * u - h * v) - f(-h * u + h * v) + f(-h * u - h * v)) / (4 * h * h)
                       for v in B] for u in B])
        ev = np.sort(np.abs(np.linalg.eigvalsh((H + H.T) / 2)))[::-1]
        ranks.append(dict(eigs=[float(f"{x:.3e}") for x in ev], rank=int((ev > 1e-4 * ev[0]).sum())))
    return dict(minimisations=starts, zeros=len(zeros), zeros_off_real_type=len(off),
                min_overlap_off_real_type=min((z[1] for z in off), default=None), hessian_at_zeros=ranks)


def nearest_pseudo_reflection(U):
    """merge the two closest eigenvalues of U (mod phase): returns (psi, b) of the pseudo-reflection
    1 + (e^{ib} - 1)|psi><psi| (up to phase) whose distinct eigenvector psi is U's isolated one"""
    V = U / np.linalg.det(U) ** (1 / 3)
    w, Q = np.linalg.eig(V)
    Q, _ = np.linalg.qr(Q)                                       # eigenvectors orthonormal (U normal)
    ph = np.angle(w)
    pairs = [(abs(np.angle(np.exp(1j * (ph[i] - ph[j])))), i, j) for i in range(3) for j in range(i + 1, 3)]
    _, i, j = min(pairs)
    k = 3 - i - j
    m = np.angle(np.exp(1j * (ph[i] + np.angle(np.exp(1j * (ph[j] - ph[i]))) / 2)))
    return Q[:, k], float(np.angle(np.exp(1j * (ph[k] - m))))


def populations(Cl, G, seeds=range(113690000, 113690400)):
    gen, deg = [], []
    for s in seeds:
        U, J = M.descend(M.haar(3, np.random.default_rng(s)), 3, Cl, G)
        if not (abs(J) < 1e-12 and M.rev_distance(U, Cl) > 0.05 and P7.J(U, 4, Cl) > 1e-6):
            continue
        V = U / np.linalg.det(U) ** (1 / 3)
        ph = np.sort(np.angle(np.linalg.eigvals(V)))
        gap = float(np.min(np.diff(np.concatenate([ph, [ph[0] + 2 * np.pi]]))) / (2 * np.pi))
        if gap < 0.005:
            psi, b = nearest_pseudo_reflection(U)
            par = lambda x: x[:3] + 1j * x[3:]
            x0 = np.concatenate([psi.real, psi.imag])
            r = minimize(lambda x: E.pseudo_reflection(3, par(x), b, Cl), x0, method='BFGS', options=dict(gtol=1e-14))
            psi2 = par(r.x) / np.linalg.norm(par(r.x))
            move = float(1 - abs(np.vdot(psi, psi2)))
            deg.append(dict(gap=gap, J6_on_component=float(r.fun), ray_moved=move,
                            real_type_overlap=E.real_type_overlap(psi2, Cl)))
        else:
            gen.append(gap)
    return dict(certified_spurious=len(gen) + len(deg), generic_spectrum=len(gen), near_degenerate=len(deg),
                min_gap_generic=min(gen) if gen else None,
                snapped_to_pseudo_reflection_component=dict(
                    reached_g_zero=sum(d["J6_on_component"] < 1e-12 for d in deg),
                    max_ray_movement=max((d["ray_moved"] for d in deg), default=None),
                    off_real_type=sum(d["real_type_overlap"] < 0.999 for d in deg)))


def s3_control(rng, n=200):
    wl = R.Weyl(1)
    mx = 0.0
    for _ in range(n):
        psi = rng.normal(size=3) + 1j * rng.normal(size=3)
        psi /= np.linalg.norm(psi)
        chi = np.einsum('i,pij,j->p', psi.conj(), wl.W, psi)
        mx = max(mx, abs(np.sum(chi ** 3).imag))
    return dict(rays=n, max_abs_imag_S3=mx)


def run():
    Cl = P7.clifford_group(3)
    G = M.hermitian_basis(3)
    rng = np.random.default_rng(11419)
    res = dict(pass_id=11419)
    res["factorisation"] = factorisation(Cl, rng)
    print(res["factorisation"], flush=True)
    res["projector_identity"] = projector_identity(Cl, rng)
    print(res["projector_identity"], flush=True)
    res["g_zero_set"] = g_zero_set(Cl, rng)
    print(res["g_zero_set"], flush=True)
    res["populations"] = populations(Cl, G)
    print(res["populations"], flush=True)
    res["S3_control"] = s3_control(rng)
    print(res["S3_control"], flush=True)
    return res


def main():
    if "--fact" in sys.argv:
        res = json.load(open(OUT))
        res["factorisation"] = factorisation(P7.clifford_group(3), np.random.default_rng(11419))
        res["projector_identity"] = projector_identity(P7.clifford_group(3), np.random.default_rng(114191))
        res["odd_ray_invariants"] = odd_ray_invariants()
        res["moment_differences_max"] = moments(P7.clifford_group(3), np.random.default_rng(114192))
        print(res["factorisation"], res["projector_identity"], res["odd_ray_invariants"], res["moment_differences_max"])
        json.dump(res, open(OUT, "w"), indent=1, default=str)
        return
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
