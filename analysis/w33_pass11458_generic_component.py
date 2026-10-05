"""Pass 11458: the generic-spectrum component of J_6's blind spot -- refined points, local dimension, and the achiral
eigenvector.

Pass 11434 observed that 79% of generic-spectrum J_6-spurious points (descents converged to J_6 ~ 1e-14, i.e. about
1e-7 from the zero set) have an eigenvector with Delta_6 < 1e-8 (Haar: 13%).  Here several such points are REFINED by
Gauss-Newton on the full residual F(U) = A_3(U) - A_3(U^T) (729 x 729, Pass 11369) until ||F|| ~ 1e-14, then:
  * the rank of dF there (the local codimension of the zero set);
  * Delta_6 = 349920 h6^2 of each eigenvector (Pass 11434), now resolved far below the earlier 1e-8 ambiguity;
  * the distance to the reversible set and J_8 (certifying the point non-reversible).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11357_jarlskog_degree_by_dimension as P7  # noqa: E402
import w33_pass11369_j6_completeness as M  # noqa: E402
import w33_pass11419_spurious_family as S  # noqa: E402

OUT = ROOT / "data" / "w33_pass11458_generic_component.json"


def residual_factory(Cl):
    def pi3(V):
        a = np.kron(np.kron(V, V), V)
        return np.kron(a, a.conj())

    def F(U):
        A = sum(pi3(C @ U @ C.conj().T) for C in Cl) / len(Cl)
        B = sum(pi3(C @ U.T @ C.conj().T) for C in Cl) / len(Cl)
        return (A - B).ravel()
    return F


def refine(U, Cl, G, F, iters=3):
    norms, sv = [], None
    for _ in range(iters):
        f0 = F(U)
        h = 1e-6
        Jm = np.array([(F(U @ M.expi(h * g)) - F(U @ M.expi(-h * g))) / (2 * h) for g in G]).T
        A = np.concatenate([Jm.real, Jm.imag])
        sv = np.linalg.svd(A, compute_uv=False)
        norms.append(float(np.linalg.norm(f0)))
        step = np.linalg.lstsq(A, -np.concatenate([f0.real, f0.imag]), rcond=1e-9)[0]
        U = U @ M.expi(np.einsum('a,aij->ij', step, G))
    norms.append(float(np.linalg.norm(F(U))))
    return U, norms, sv


def eig_delta6(U, Cl):
    V = U / np.linalg.det(U) ** (1 / 3)
    w, Q = np.linalg.eig(V)
    Q, _ = np.linalg.qr(Q)
    return sorted(float(S.ray_witness(Q[:, j], Cl)) for j in range(3))


def continuous_symmetry_test(pairs=80, seed=1):
    """HIDDEN-SYMMETRY TEST.  If a group G' > Cl had the same degree-(3,3) commutant, A_3 would be G'-conjugation
    invariant and every U with U^T ~ g U g^-1 (g in G') would have J_6 = 0 -- a candidate explanation of the blind spot.
    A continuous G' needs H in su(3) with d/de f_Y(e^{ieH} U e^{-ieH}) = 0 for all the spanning invariants
    f_Y(U) = avg_C |tr(Y^dag C U C^dag)|^6.  Returns the singular values of that derivative map (8 su(3) directions);
    no zero singular value = no continuous hidden symmetry.  (A finite G' strictly containing the Hessian group of
    order 216 does not exist in PU(3): Blichfeldt's classification of finite primitive subgroups.)"""
    Cl = P7.clifford_group(3)
    G = M.hermitian_basis(3)
    rng = np.random.default_rng(seed)

    def f(Y, U):
        return np.mean(np.abs(np.einsum('ij,cjk,kl,cil->c', Y.conj().T, Cl, U, Cl.conj())) ** 6)

    rows = []
    for _ in range(pairs):
        Y, U = M.haar(3, rng), M.haar(3, rng)
        h = 1e-5
        rows.append([(f(Y, M.expi(h * H) @ U @ M.expi(-h * H)) - f(Y, M.expi(-h * H) @ U @ M.expi(h * H))) / (2 * h)
                     for H in G])
    sv = np.linalg.svd(np.array(rows), compute_uv=False)
    return [float(x) for x in sv / sv[0]]


def run(max_points=6):
    Cl = P7.clifford_group(3)
    G = M.hermitian_basis(3)
    F = residual_factory(Cl)
    rows = []
    for s in range(113690000, 113690400):
        U, J = M.descend(M.haar(3, np.random.default_rng(s)), 3, Cl, G)
        if not (abs(J) < 1e-12 and M.rev_distance(U, Cl) > 0.05 and P7.J(U, 4, Cl) > 1e-6):
            continue
        V = U / np.linalg.det(U) ** (1 / 3)
        ph = np.sort(np.angle(np.linalg.eigvals(V)))
        if np.min(np.diff(np.concatenate([ph, [ph[0] + 2 * np.pi]]))) / (2 * np.pi) < 0.01:
            continue                                                    # pseudo-reflection component (Pass 11419)
        before = eig_delta6(U, Cl)
        U2, norms, sv = refine(U, Cl, G, F)
        after = eig_delta6(U2, Cl)
        rows.append(dict(seed=s, residual_norms=norms, dF_rank=int((sv > 1e-6 * sv[0]).sum()),
                         singular_values_relative=[float(x) for x in sv / sv[0]],
                         eigvec_Delta6_before=before, eigvec_Delta6_after=after,
                         distance=M.rev_distance(U2, Cl), J8=float(P7.J(U2, 4, Cl))))
        print(rows[-1], flush=True)
        if len(rows) >= max_points:
            break
    return dict(pass_id=11458, points=rows, continuous_symmetry_singular_values=continuous_symmetry_test(),
                min_eigvec_Delta6_after=[r["eigvec_Delta6_after"][0] for r in rows],
                dF_ranks=[r["dF_rank"] for r in rows])


def main():
    if "--sym" in sys.argv:
        res = json.load(open(OUT))
        res["continuous_symmetry_singular_values"] = continuous_symmetry_test()
        json.dump(res, open(OUT, "w"), indent=1, default=str)
        print(res["continuous_symmetry_singular_values"])
        return
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
