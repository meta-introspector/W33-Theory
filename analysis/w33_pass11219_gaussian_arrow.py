#!/usr/bin/env python3
"""Pass 11219: the arrow law over the reals -- Gaussian (free bosonic) dynamics.

Pass 11207 proved A(S) = n - c(S) for Sp(2n, F_q), q odd, and Pass 11217 for q = 2.  Every step of the odd proof is
field-general for characteristic != 2: the characteristic-free lower bound, the Wall/Milnor orthogonal decomposition,
the bi-Lagrangian radical lemma, Taussky-Zassenhaus graphs, the Cayley transform, the eigenvector reduction and the
Hermitian real form.  So the law holds for S in Sp(2n, R): the symplectic matrices of Gaussian unitaries, i.e. of free
bosonic dynamics.

What A means over R.  For a mode decomposition V = P_1 + ... + P_n (symplectic planes = modes), dim(P cap SP) =
2 - rank S[P-perp <- P], the rank of the block coupling mode P to the rest.  So
    A(S) = min over mode decompositions of sum_k rank S[P_k-perp <- P_k]
counts the quadrature channels by which modes must leak into each other per tick, and A = n - c(S) with, for the real
Jordan form,
    c(S) = #elliptic size-1 pairs (e^{+-i theta}) + #hyperbolic size-1 pairs (lambda, 1/lambda real) + sum_{+-1} (m_1/2 + m_2).
For semisimple S without +-1 this is A = 2 x (number of loxodromic quartets r e^{+-i theta}, r^-1 e^{+-i theta}).
  * Stable free dynamics (positive Hamiltonian): every mode is elliptic, c = n, A = 0 -- the normal modes are the
    subsystems it keeps to itself.  Pure squeezing (hyperbolic) also has A = 0.
  * The arrow is carried by complex instability (loxodromic quartets, two channels each) and by Jordan defects.
  * At a Krein collision of two elliptic modes of OPPOSITE Krein signature (Hamiltonian-Hopf), the generic collision is
    non-semisimple and A jumps 0 -> 2, staying 2 on the loxodromic side; modes of the SAME signature pass through each
    other with A = 0 throughout (Krein's theorem, here as a statement about the arrow).

Checks (numerical, tolerance 1e-8 relative):
  1. random S = M S0 M^-1 (M random symplectic) of every type mix (e elliptic, h hyperbolic, l loxodromic, n <= 6):
     an explicit split with total coupling rank exactly 2l = n - c is built (invariant planes + a Taussky-Zassenhaus
     split per quartet) and verified;
  2. non-semisimple pieces (elliptic Jordan block at a Krein collision, V(4) at +1, a hyperbolic Jordan pair): a
     transverse bi-Lagrangian is found numerically and its split has coupling rank exactly n;
  3. Krein-collision families (modes w1 = 1, w2 = 1.2, coupling eps (x1 x2 - p1 p2)): opposite signature collides at
     eps = 0.1 -- non-semisimple there, loxodromic beyond with growth sqrt(eps^2 - 0.01) -- and A switches 0 -> 2
     exactly there; equal signature never collides and keeps A = 0, and an exact equal-signature degeneracy is
     semisimple (A = 0).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm, null_space, schur
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "w33_pass11219_gaussian_arrow.json"
TOL = 1e-8


def Jm(n):
    return np.block([[np.zeros((n, n)), np.eye(n)], [-np.eye(n), np.zeros((n, n))]])


def rand_symplectic(n, rng, scale=0.5):
    H = rng.normal(size=(2 * n, 2 * n))
    return expm(scale * Jm(n) @ (H + H.T) / 2)


def rot(t):
    return np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])


def num_rank(X, tol=TOL):
    s = np.linalg.svd(X, compute_uv=False)
    return int((s > tol * max(s[0], 1.0)).sum())


def coupling_ranks(S, planes):
    """rank of the block of S mapping mode P into the complement, for each plane (columns = basis of P)"""
    n = len(S) // 2
    J = Jm(n)
    out = []
    for k, P in enumerate(planes):
        others = np.hstack([Q for i, Q in enumerate(planes) if i != k]) if n > 1 else np.zeros((2 * n, 0))
        # coordinates of S P in the basis (P, others): block on 'others' is the coupling
        B = np.hstack([P, others])
        coeff = np.linalg.solve(B, S @ P)
        out.append(num_rank(coeff[2:]) if n > 1 else 0)
    return out


def is_split(planes, n):
    J = Jm(n)
    ok = True
    for i, P in enumerate(planes):
        ok &= abs(P[:, 0] @ J @ P[:, 1]) > 1e-6
        for Q in planes[i + 1:]:
            ok &= np.abs(P.T @ J @ Q).max() < 1e-8 * (1 + np.abs(P).max() * np.abs(Q).max())
    return bool(ok)


# ---------------------------------------------------------------- normal-form building blocks (in local coordinates)
def embed(blocks, n):
    """blocks: list of (S_block (2k x 2k) in (x..., p...) order, k).  Returns S0 on R^2n in (x1..xn, p1..pn) order and
    the index lists of each block"""
    S0 = np.zeros((2 * n, 2 * n))
    pos, idx = 0, []
    for Sb, k in blocks:
        ii = list(range(pos, pos + k)) + list(range(n + pos, n + pos + k))
        S0[np.ix_(ii, ii)] = Sb
        idx.append((ii, k))
        pos += k
    return S0, idx


def elliptic(theta):
    """one mode rotated by theta in its (x, p) plane"""
    return np.array([[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]]), 1


def hyperbolic(lam):
    return np.diag([lam, 1 / lam]), 1


def loxodromic(r, theta):
    A = r * rot(theta)
    return np.block([[A, np.zeros((2, 2))], [np.zeros((2, 2)), np.linalg.inv(A).T]]), 2


def tz_split_quartet(Sb):
    """Taussky-Zassenhaus graph split of a loxodromic quartet in local (x1, x2, p1, p2) coordinates"""
    A = Sb[:2, :2]
    basis = [np.array([[1., 0], [0, 0]]), np.array([[0., 1], [1, 0]]), np.array([[0., 0], [0, 1]])]
    K = null_space(np.array([(E @ A - A.T @ E).ravel() for E in basis]).T)
    Phi = sum(c * E for c, E in zip(K[:, 0], basis))
    L = np.vstack([np.eye(2), Phi])
    return planes_from_bilagrangian(Sb, L)


def planes_from_bilagrangian(S, L):
    """diagonalise b(x, y) = omega(x, Sy) on a transverse bi-Lagrangian L (columns) and return the planes <e, Se>"""
    n = len(S) // 2
    J = Jm(n)
    B = L.T @ J @ S @ L
    w, Q = np.linalg.eigh((B + B.T) / 2)
    E = L @ Q
    return [np.stack([E[:, i], S @ E[:, i]], 1) for i in range(E.shape[1])], dict(
        lagrangian=float(np.abs(L.T @ J @ L).max()), bi=float(np.abs(L.T @ J @ (S + np.linalg.inv(S)) @ L).max()),
        b_min=float(np.abs(w).min()))


def build_generic(e, h, l, rng):
    """random S of type (e elliptic, h hyperbolic, l loxodromic) and the constructed split; returns verification"""
    n = e + h + 2 * l
    blocks = [elliptic(rng.uniform(0.2, 3.0)) for _ in range(e)] + \
             [hyperbolic(rng.uniform(1.2, 2.5) * rng.choice([-1, 1])) for _ in range(h)] + \
             [loxodromic(rng.uniform(1.1, 2.0), rng.uniform(0.3, 2.8)) for _ in range(l)]
    S0, idx = embed(blocks, n)
    planes0 = []
    for (ii, k), (Sb, _) in zip(idx, blocks):
        if k == 1:
            P = np.zeros((2 * n, 2))
            P[ii[0], 0] = 1
            P[ii[1], 1] = 1
            planes0.append(P)
        else:
            sub, _ = tz_split_quartet(Sb)
            for Pl in sub:
                P = np.zeros((2 * n, 2))
                P[ii, :] = Pl
                planes0.append(P)
    M = rand_symplectic(n, rng)
    S = M @ S0 @ np.linalg.inv(M)
    planes = [M @ P for P in planes0]
    ranks = coupling_ranks(S, planes)
    return dict(e=e, h=h, l=l, n=n, c=e + h, split=is_split(planes, n), total_rank=int(sum(ranks)),
                law=int(sum(ranks)) == n - (e + h), ranks=ranks,
                symplectic_err=float(np.abs(S.T @ Jm(n) @ S - Jm(n)).max()))


# ---------------------------------------------------------------- non-semisimple pieces: numeric bi-Lagrangian
def find_bilagrangian(S, rng, tries=40):
    """L = graph of symmetric Z over the x-Lagrangian, solving omega(., T .) = 0 on L, transverse to SL"""
    n = len(S) // 2
    J = Jm(n)
    T = S + np.linalg.inv(S)
    iu = np.triu_indices(n)

    def L_of(z):
        Z = np.zeros((n, n))
        Z[iu] = z
        Z = Z + Z.T - np.diag(np.diag(Z))
        return np.vstack([np.eye(n), Z])

    def res(z):
        L = L_of(z)
        R = L.T @ J @ T @ L
        return R[np.triu_indices(n, 1)]
    for _ in range(tries):
        sol = least_squares(res, rng.normal(size=len(iu[0])), xtol=1e-15, ftol=1e-15, gtol=1e-15)
        L = L_of(sol.x)
        if np.abs(res(sol.x)).max() > 1e-10:
            continue
        sv = np.linalg.svd(np.hstack([L, S @ L]), compute_uv=False)
        if sv[-1] / sv[0] > 1e-6:                              # transverse
            planes, info = planes_from_bilagrangian(S, L)
            if info["b_min"] > 1e-6:
                return planes, info
    return None, None


def nonsemisimple_cases(rng):
    """the c = 0 pieces that are not semisimple, each conjugated by a random symplectic M"""
    J = Jm(2)
    # Hamiltonian-Hopf normal form at criticality (Meyer-Hall): H = w (x2 p1 - x1 p2) + (x1^2 + x2^2)/2, coordinates
    # (x1, x2, p1, p2): a double non-semisimple elliptic pair +-i w -- the Krein collision point
    w = 1.0
    Hhh = np.zeros((4, 4))
    Hhh[1, 2] = Hhh[2, 1] = w
    Hhh[0, 3] = Hhh[3, 0] = -w
    Hhh[0, 0] = Hhh[1, 1] = 1.0
    # V(4): H = p1 x2 + p2^2/2 is a regular nilpotent Hamiltonian matrix (p1 -> p2 -> x2 -> x1)
    Hv = np.zeros((4, 4))
    Hv[1, 2] = Hv[2, 1] = 1.0
    Hv[3, 3] = 1.0
    A = 1.7 * np.array([[1.0, 1.0], [0.0, 1.0]])
    hyp_jordan = np.block([[A, np.zeros((2, 2))], [np.zeros((2, 2)), np.linalg.inv(A).T]])
    cases = {"elliptic Jordan pair (Krein collision point)": expm(0.7 * J @ Hhh),
             "unipotent V(4)": expm(J @ Hv),
             "hyperbolic Jordan pair": hyp_jordan}
    out = {}
    for name, S0 in cases.items():
        M = rand_symplectic(2, rng)
        S = M @ S0 @ np.linalg.inv(M)
        ev = np.linalg.eigvals(S)
        jordan = max(len(S) - num_rank(S - z * np.eye(4), 1e-7) for z in ev)  # max geometric multiplicity
        planes, info = find_bilagrangian(S, rng)
        ranks = coupling_ranks(S, planes) if planes else None
        out[name] = dict(eigenvalues=sorted([[round(float(z.real), 6), round(float(z.imag), 6)] for z in ev]),
                         max_geometric_multiplicity=int(jordan), found=planes is not None,
                         split=bool(planes and is_split(planes, 2)),
                         total_rank=None if ranks is None else int(sum(ranks)), info=info)
    return out


def c_two_modes(S, tol=1e-6):
    """c for a 4x4 symplectic S without eigenvalues +-1: 2 if semisimple with unimodular or real spectrum (two invariant
    planes), else 0 (a loxodromic quartet or a non-semisimple pair is one 4-dimensional indecomposable)"""
    ev = np.linalg.eigvals(S)
    if np.any(np.abs(np.abs(ev) - 1) > 1e-7) and np.any(np.abs(ev.imag) > 1e-7):
        return 0, "loxodromic"
    for z in ev:
        alg = int(np.sum(np.abs(ev - z) < 1e-5))
        geo = len(S) - num_rank(S - z * np.eye(len(S)), tol)
        if geo < alg:
            return 0, "non-semisimple"
    return 2, "semisimple"


def krein_family(sign, eps_list, w1=1.0, w2=1.2, dt=0.6):
    """modes of frequency w1, w2 with Krein signatures + and sign, coupled by eps (x1 x2 - p1 p2).  For opposite
    signature the pair collides at eps = |w1 - w2|/2 and turns loxodromic; for equal signature it never collides"""
    rows = []
    for eps in eps_list:
        H = np.diag([w1, sign * w2, w1, sign * w2])
        H[0, 1] = H[1, 0] = eps
        H[2, 3] = H[3, 2] = -eps
        K = Jm(2) @ H
        S = expm(dt * K)
        c, kind = c_two_modes(S)
        rows.append(dict(eps=eps, kind=kind, growth=float(np.abs(np.linalg.eigvals(K).real).max()), c=c, A=2 - c))
    return rows


def collision_points(dt=0.6):
    """exact degeneracies: opposite signature at eps = |w1 - w2|/2 is non-semisimple (A = 2); equal signature with
    w1 = w2 is semisimple (A = 0) -- Krein's theorem read as a statement about the arrow"""
    out = {}
    H = np.diag([1.0, -1.2, 1.0, -1.2])
    H[0, 1] = H[1, 0] = 0.1
    H[2, 3] = H[3, 2] = -0.1
    c, kind = c_two_modes(expm(dt * Jm(2) @ H), tol=1e-6)
    out["opposite signature, eps = |w1-w2|/2"] = dict(kind=kind, c=c, A=2 - c)
    c, kind = c_two_modes(expm(dt * Jm(2) @ np.eye(4)), tol=1e-6)
    out["equal signature, w1 = w2"] = dict(kind=kind, c=c, A=2 - c)
    return out


def run(seed=11219):
    rng = np.random.default_rng(seed)
    res = dict(pass_id=11219)
    gen = []
    for (e, h, l) in [(2, 0, 0), (0, 2, 0), (1, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1), (2, 1, 1), (0, 0, 2),
                      (1, 1, 2), (0, 0, 3), (2, 2, 1)]:
        for _ in range(3):
            gen.append(build_generic(e, h, l, rng))
    res["generic"] = gen
    res["generic_law"] = all(g["law"] and g["split"] for g in gen)
    res["nonsemisimple"] = nonsemisimple_cases(rng)
    res["nonsemisimple_law"] = all(v["found"] and v["split"] and v["total_rank"] == 2 for v in res["nonsemisimple"].values())
    # Krein: opposite signature (sign = -1) vs same signature (sign = +1)
    eps = [0.0, 0.03, 0.06, 0.09, 0.099, 0.101, 0.11, 0.15, 0.2, 0.3, 0.5]
    res["krein_opposite"] = krein_family(-1, eps, w1=1.0, w2=1.2)
    res["krein_same"] = krein_family(+1, eps, w1=1.0, w2=1.2)
    res["krein_collision_points"] = collision_points()
    op = res["krein_opposite"]
    res["krein_arrow_switches_on_at_collision"] = all((r["A"] == 2) == (r["eps"] > 0.1) for r in op)
    res["krein_same_never"] = all(r["A"] == 0 for r in res["krein_same"])
    res["growth_matches_sqrt"] = all(abs(r["growth"] - np.sqrt(max(r["eps"] ** 2 - 0.01, 0))) < 1e-9 for r in op)
    return res


def main():
    res = run()
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True, default=float))
    print(json.dumps({k: v for k, v in res.items() if k != "generic"}, indent=1, default=float))
    print("generic law:", res["generic_law"], [(g["e"], g["h"], g["l"], g["total_rank"]) for g in res["generic"]])


if __name__ == "__main__":
    main()
