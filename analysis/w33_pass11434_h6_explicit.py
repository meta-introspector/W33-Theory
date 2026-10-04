"""Pass 11434: the T-odd ray invariant h6 made explicit, and a first look at the generic component.

Pass 11419: the first T-odd Clifford-invariant polynomial of a qutrit ray has bidegree (6,6) and is unique there (h6);
J_6(1 + lam|psi><psi|) = |lam|^12 Delta_6(psi), Delta_6 = 2|h6|^2.

EXPLICIT FORM.  The 12 stabiliser states s (the four MUBs Z, X, XZ, XZ^2) determine every ray through
p_s(psi) = |<s|psi>|^2.  The extended Clifford group (216 unitaries C, 216 antiunitaries C K; sign +1 / -1) permutes
them.  For a monomial m = prod_i p_{s_i}, the odd Reynolds average
        R_m(psi) = (1/432) sum_g sgn(g) m(g psi)
is T-odd and Clifford-invariant of bidegree (6,6) when deg m = 6, so it is a multiple of h6.  All 74 orbits of degree-6
monomials are scanned; the simplest nonzero one is
        m = p_a^3 p_b^2 p_c,   a, b, c stabiliser states from three DIFFERENT MUBs (Z, X, XZ),
and  Delta_6 = 349920 * R_m^2  exactly (constant ratio on random rays).  So h6 is a signed sum over the Clifford images
of an ORIENTED triple of mutually unbiased stabiliser states with exponents 3, 2, 1: a chirality of the state relative
to the MUB geometry.

GENERIC COMPONENT (observation).  For J_6-spurious points with generic spectrum, the smallest Delta_6 among the three
eigenvectors is compared with that of Haar-random unitaries.
"""

from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11357_jarlskog_degree_by_dimension as P7  # noqa: E402
import w33_pass11369_j6_completeness as M  # noqa: E402
import w33_pass11419_spurious_family as S  # noqa: E402

OUT = ROOT / "data" / "w33_pass11434_h6_explicit.json"


def stabiliser_states():
    w = np.exp(2j * np.pi / 3)
    X = np.roll(np.eye(3), 1, axis=0)
    Z = np.diag([1, w, w * w])
    st, mub = [], []
    for m, A in enumerate((Z, X, X @ Z, X @ Z @ Z)):
        _, V = np.linalg.eig(A)
        for j in range(3):
            st.append(V[:, j] / np.linalg.norm(V[:, j]))
            mub.append(m)
    return np.array(st), mub


def group_permutations(Cl, st):
    def idx(v):
        o = np.abs(st.conj() @ v)
        k = int(np.argmax(o))
        assert abs(o[k] - 1) < 1e-8
        return k
    perm = [np.array([idx(C @ s) for s in st]) for C in Cl]
    conj = np.array([idx(s.conj()) for s in st])
    return [(p, 1) for p in perm] + [(p[conj], -1) for p in perm]


def odd_reynolds(m, psi, Cl, st):
    psi = psi / np.linalg.norm(psi)
    tot = 0.0
    for C in Cl:
        a = np.abs(st.conj() @ (C @ psi)) ** 2
        b = np.abs(st.conj() @ (C @ psi.conj())) ** 2
        tot += np.prod(a[list(m)]) - np.prod(b[list(m)])
    return tot / (2 * len(Cl))


def h6(psi, Cl=None, st=None):
    """h6 up to normalisation: the odd Reynolds average of p_a^3 p_b^2 p_c (a, b, c in the Z, X, XZ bases)"""
    Cl = P7.clifford_group(3) if Cl is None else Cl
    st = stabiliser_states()[0] if st is None else st
    return odd_reynolds((0, 0, 0, 3, 3, 6), psi, Cl, st)


def run():
    Cl = P7.clifford_group(3)
    st, mub = stabiliser_states()
    G = group_permutations(Cl, st)
    res = dict(pass_id=11434, extended_clifford_order=len(G), distinct_permutations=len({tuple(g[0]) for g in G}))
    seen, reps = set(), []
    for m in itertools.combinations_with_replacement(range(12), 6):
        if m in seen:
            continue
        orb = {tuple(sorted(p[list(m)])) for p, _ in G}
        seen |= orb
        reps.append((m, len(orb)))
    rng = np.random.default_rng(11434)
    psis = [rng.normal(size=3) + 1j * rng.normal(size=3) for _ in range(6)]
    d6 = [S.ray_witness(p, Cl) for p in psis]
    nonzero = []
    for m, size in reps:
        vals = [odd_reynolds(m, p, Cl, st) for p in psis]
        if max(abs(v) for v in vals) > 1e-10:
            ratios = [d6[i] / vals[i] ** 2 for i in range(len(psis))]
            nonzero.append(dict(monomial=m, mubs=[mub[i] for i in m], orbit=size,
                                ratio_mean=float(np.mean(ratios)), ratio_rel_spread=float(np.ptp(ratios) / np.mean(ratios))))
    res["monomial_orbits"] = len(reps)
    res["nonzero_odd_orbit_sums"] = nonzero
    res["all_nonzero_proportional_to_h6"] = all(n["ratio_rel_spread"] < 1e-6 for n in nonzero)
    simplest = nonzero[0]
    res["explicit_h6"] = dict(monomial="p_a^3 p_b^2 p_c", states=simplest["monomial"], mubs=simplest["mubs"],
                              Delta6_over_h6_squared=simplest["ratio_mean"])
    # generic component: smallest eigenvector Delta_6
    G8 = M.hermitian_basis(3)
    spur, rand = [], []
    for s in range(113690000, 113690300):
        U, J = M.descend(M.haar(3, np.random.default_rng(s)), 3, Cl, G8)
        if abs(J) < 1e-12 and M.rev_distance(U, Cl) > 0.05 and P7.J(U, 4, Cl) > 1e-6:
            V = U / np.linalg.det(U) ** (1 / 3)
            ph = np.sort(np.angle(np.linalg.eigvals(V)))
            if np.min(np.diff(np.concatenate([ph, [ph[0] + 2 * np.pi]]))) / (2 * np.pi) > 0.01:
                _, Q = np.linalg.eig(V)
                Q, _ = np.linalg.qr(Q)
                spur.append(min(S.ray_witness(Q[:, j], Cl) for j in range(3)))
    for _ in range(300):
        _, Q = np.linalg.eig(M.haar(3, rng))
        Q, _ = np.linalg.qr(Q)
        rand.append(min(S.ray_witness(Q[:, j], Cl) for j in range(3)))
    spur, rand = np.array(spur), np.array(rand)
    res["generic_component_eigenvectors"] = dict(
        spurious_points=len(spur), haar_controls=len(rand),
        fraction_min_Delta6_below_1e_8=dict(spurious=float((spur < 1e-8).mean()), haar=float((rand < 1e-8).mean())),
        median_min_Delta6=dict(spurious=float(np.median(spur)), haar=float(np.median(rand))))
    print(json.dumps(res, indent=1, default=str), flush=True)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
