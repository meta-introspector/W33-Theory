"""Pass 11502: (A) the rank-one cancellation in Delta2 = Delta6/108 proved by an exact census; (B) the generator degrees of
the time-even ring and the time-odd module of qutrit-state Clifford invariants.

(A) Pass 11492 derived Delta2(psi x s) = (1/120 + (27/40) 3^-6) Delta6(psi) = Delta6/108 from the compressions
A_C = (1 x <0|) C (1 x |0>), but the vanishing of the rank-one terms' odd part was only numerical.  Here: every rank-one
compression is |x><y| (unit singular value) with x, y STABILISER states, and the 1,259,712 of them are spread EXACTLY
uniformly over the 144 ordered stabiliser pairs.  The odd part is then const * (sum_x p_x^6)(sum_y p_y^6 - sum_y p_conj(y)^6) = 0,
since complex conjugation permutes the 12 stabiliser states.  With Pass 11492's 1/120, 27/40 and uniformity, the identity
Delta2 = Delta6/108 is proved.

(B) The 12 stabiliser probabilities p_s = |<s|psi>|^2 span all Hermitian quadratic forms, so the Reynolds averages
R_m = sum_g sgn^eps(g) m(g psi) of p-monomials m span all even (eps = 0) and odd (eps = 1) invariants of each degree.  Evaluated
on NON-normalised random psi (degrees stay visible), they give bases of E_k and O_k (dimensions checked against Pass 11491's
Molien series).  New generators in degree k:
    even ring:   E_k - dim span{ e_i e_j : deg e_i + deg e_j = k, both >= 1 }
    odd module:  O_k - dim span{ e o : e in E_j (j >= 1), o in O_(k-j) }
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11357_jarlskog_degree_by_dimension as P7  # noqa: E402
import w33_pass11434_h6_explicit as H  # noqa: E402

OUT = ROOT / "data" / "w33_pass11502_time_odd_module.json"
_M = json.load(open(ROOT / "data" / "w33_pass11491_time_odd_molien.json"))["sequences"]
MOLIEN_E, MOLIEN_O = _M["even"], _M["odd"]                 # Pass 11491 (proved closed forms)


def rank_one_census():
    import w33_pass11330_orbit_census as O
    import w33_pass11350_linear_decider as L
    D = L.Decider(2)
    Ms = np.array(O.all_symplectic(D.wl)[0]) % 3
    st, _ = H.stabiliser_states()
    conj = np.array([int(np.argmax(np.abs(st.conj() @ s.conj()))) for s in st])
    Iso = np.kron(np.eye(3), np.array([[1.0], [0.0], [0.0]]))
    N = np.zeros((12, 12), np.int64)
    nonstab = 0
    scales = set()
    for M in Ms:
        A = np.einsum('ip,apq,qj->aij', Iso.T, np.einsum('apq,qr->apr', D.wl.W, D.weil(M)), Iso)
        sv = np.linalg.svd(A, compute_uv=False)
        for Ai, s_ in zip(A, sv):
            if s_[0] > 1e-9 and s_[1] < 1e-9:
                u, _, vh = np.linalg.svd(Ai)
                ox, oy = np.abs(st.conj() @ u[:, 0]), np.abs(st.conj() @ vh[0].conj())
                if abs(ox.max() - 1) > 1e-8 or abs(oy.max() - 1) > 1e-8:
                    nonstab += 1
                    continue
                N[int(np.argmax(ox)), int(np.argmax(oy))] += 1
                scales.add(round(float(s_[0]), 9))
    return dict(rank_one=int(N.sum()) + nonstab, non_stabiliser=nonstab, singular_values=sorted(scales),
                pair_counts_min=int(N.min()), pair_counts_max=int(N.max()), uniform_over_144_pairs=bool(N.min() == N.max()),
                conjugation_symmetric=bool((N == N[:, conj]).all()))


class Invariants:
    def __init__(self, n_points=160, seed=11502):
        Cl = P7.clifford_group(3)
        self.st, _ = H.stabiliser_states()
        self.perms = H.group_permutations(Cl, self.st)                        # 432 (permutation, sign)
        rng = np.random.default_rng(seed)
        psi = rng.normal(size=(n_points, 3)) + 1j * rng.normal(size=(n_points, 3))
        psi /= np.linalg.norm(psi, axis=1, keepdims=True)   # a bidegree-(k,k) form is fixed by its values on the sphere
        self.p = np.abs(psi @ self.st.conj().T) ** 2
        self.Pp = np.stack([self.p[:, perm] for perm, _ in self.perms])        # (432, points, 12)
        self.sign = np.array([s for _, s in self.perms], float)
        self.rng = rng

    def reynolds(self, m, odd):
        vals = np.prod(self.Pp[:, :, list(m)], axis=2)                        # (432, points)
        w = self.sign if odd else np.ones(len(self.sign))
        return (w[:, None] * vals).sum(axis=0) / len(w)

    def basis(self, k, odd, target, max_tries=30000):
        if target == 0:
            return np.zeros((self.p.shape[0], 0))
        cols = []
        B = np.zeros((self.p.shape[0], 0))
        for _ in range(max_tries):
            m = tuple(sorted(self.rng.integers(12, size=k)))
            v = self.reynolds(m, odd)
            if np.linalg.norm(v) < 1e-12:
                continue
            v = v / np.linalg.norm(v)
            C = np.concatenate([B, v[:, None]], 1)
            if rank(C) > B.shape[1]:
                B = C
                cols.append(m)
            if B.shape[1] == target:
                break
        return B


GAPS = []


def rank(A, rel=1e-9):
    if A.shape[1] == 0:
        return 0
    A = A / np.linalg.norm(A, axis=0, keepdims=True)
    s = np.linalg.svd(A, compute_uv=False)
    r = int((s > rel * s[0]).sum())
    GAPS.append((float(s[r - 1] / s[0]), float(s[r] / s[0]) if r < len(s) else 0.0))
    return r


def module_structure(kmax=13):
    I = Invariants()
    E = {k: I.basis(k, False, MOLIEN_E[k]) for k in range(0, kmax + 1)}
    O = {k: I.basis(k, True, MOLIEN_O[k]) for k in range(0, kmax + 1)}
    out = {}
    for k in range(1, kmax + 1):
        assert E[k].shape[1] == MOLIEN_E[k] and O[k].shape[1] == MOLIEN_O[k], (k, E[k].shape, O[k].shape)
        assert rank(E[k]) == MOLIEN_E[k] and rank(O[k]) == MOLIEN_O[k] if MOLIEN_O[k] else True
        prods_e = [E[j][:, a] * E[k - j][:, b] for j in range(1, k) for a in range(E[j].shape[1]) for b in range(E[k - j].shape[1])]
        prods_o = [E[j][:, a] * O[k - j][:, b] for j in range(1, k + 1) for a in range(E[j].shape[1])
                   for b in range(O[k - j].shape[1])]
        re = rank(np.stack(prods_e, 1)) if prods_e else 0
        ro = rank(np.stack(prods_o, 1)) if prods_o else 0
        g = GAPS[-2:]
        out[str(k)] = dict(E=MOLIEN_E[k], O=MOLIEN_O[k], new_even_generators=MOLIEN_E[k] - re, gaps=g,
                           new_odd_generators=MOLIEN_O[k] - ro, odd_products_rank=ro)
        print(k, out[str(k)], flush=True)
    return out


def main():
    res = dict(pass_id=11502)
    res["rank_one_census"] = rank_one_census()
    print(res["rank_one_census"], flush=True)
    res["module_structure"] = module_structure()
    res["even_generator_degrees"] = sum(([int(k)] * v["new_even_generators"] for k, v in res["module_structure"].items()), [])
    res["odd_generator_degrees"] = sum(([int(k)] * v["new_odd_generators"] for k, v in res["module_structure"].items()), [])
    print(res["even_generator_degrees"], res["odd_generator_degrees"], flush=True)
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
