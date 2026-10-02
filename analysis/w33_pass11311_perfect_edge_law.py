"""Pass 11311: the perfect-tick magic edge law, tested over ALL two-qutrit Clifford ticks (it mostly fails).

Pass 11268 (scoped by the audit to two fixed circulant representatives) found that a perfect tick dressed with a cubic
phase on exactly one of its legs is always substrate-time-reversible, and with magic on 2n - 1 legs never.  Here the
statement is tested over the whole two-qutrit Clifford group, split by Codex's compiler classes (Pass 11193):
    local (incl. swap)  1152 symplectic classes,   perfect  13824 (of which 64 F9-linear, the class of Pass 11170's
    circulants, kept separate as 'perfectF9'),   other entangling  36864,
with the perfectness of each symplectic class verified on both nontrivial 2|2 leg cuts.

For each class: uniform random Clifford ticks V = P_a V_M (Pauli coset uniform), dressings
    U = (T^{a_1} (x) T^{a_2}) V (T^{b_1} (x) T^{b_2}),  exactly one / exactly three of the four exponents magic (not 0 mod
    3), the others uniform in {0, 3, 6} (Paulis Z^k);
every U decided exactly (Pass 11252 criterion).  Plus the k = 0..4 magic-leg profile for perfect ticks.
"""

from __future__ import annotations

import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11227_two_qutrit_cp_mixing as P2  # noqa: E402
import w33_pass11252_exact_reversibility as R  # noqa: E402

OUT = ROOT / "data" / "w33_pass11311_perfect_edge_law.json"
MAGIC = [1, 2, 4, 5, 7, 8]
CLIF = [0, 3, 6]
_S = {}


def perfect(V):
    Tn = V.reshape(3, 3, 3, 3)
    for perm in ((0, 2, 1, 3), (0, 3, 1, 2)):
        M = np.transpose(Tn, perm).reshape(9, 9)
        G = M @ M.conj().T
        if not np.allclose(G, G[0, 0] * np.eye(9), atol=1e-8):
            return False
    return True


def is_local(V):
    """V = A (x) B or (A (x) B) SWAP: operator-Schmidt rank 1 in either leg pairing"""
    Tn = V.reshape(3, 3, 3, 3)
    for perm in ((0, 2, 1, 3), (0, 3, 1, 2)):
        M = np.transpose(Tn, perm).reshape(9, 9)
        if np.linalg.matrix_rank(M, tol=1e-8) == 1:
            return True
    return False


def symplectic_part(wl, C):
    M = np.zeros((4, 4), int)
    for j in range(4):
        e = np.zeros(4, int)
        e[j] = 1
        img = C @ wl.W[wl.index(e)] @ C.conj().T
        M[:, j] = wl.labels[int(np.argmax(np.abs(np.einsum('pij,ij->p', wl.W.conj(), img))))]
    return M


JB = np.array([[0, 2, 0, 0], [1, 0, 0, 0], [0, 0, 0, 2], [0, 0, 1, 0]])     # multiplication by i in F9 on each qutrit


def classes(reps):
    wl = R.Weyl(2)
    lab = np.empty(len(reps), dtype="<U10")
    for i, V in enumerate(reps):
        if is_local(V):
            lab[i] = "local"
        elif perfect(V):
            M = symplectic_part(wl, V)
            lab[i] = "perfectF9" if (((M @ JB - JB @ M) % 3) == 0).all() else "perfect"
        else:
            lab[i] = "other"
    return lab


def local_T(e1, e2):
    T = P2.T
    return np.kron(np.linalg.matrix_power(T, e1), np.linalg.matrix_power(T, e2))


def _init():
    R.WEYL[2] = R.Weyl(2)
    _S["reps"] = np.load(P2.CACHE)


def _job(args):
    seed, idx, k, count = args
    rng = np.random.default_rng(seed)
    reps = _S["reps"]
    v = 0
    for _ in range(count):
        V = P2.PA[rng.integers(81)] @ reps[idx[rng.integers(len(idx))]]
        legs = rng.permutation(4)[:k]
        e = [int(rng.choice(CLIF)) for _ in range(4)]
        for l in legs:
            e[l] = int(rng.choice(MAGIC))
        U = local_T(e[0], e[1]) @ V @ local_T(e[2], e[3])
        if R.decide(U, 2, rng)[0] is False:
            v += 1
    return v


def run(per_cell=6000):
    reps = np.load(P2.CACHE)
    lab = classes(reps)
    res = dict(pass_id=11311, class_sizes={c: int(np.sum(lab == c)) for c in ("local", "perfectF9", "perfect", "other")})
    print(res["class_sizes"], flush=True)
    cells = {}
    with Pool(11, initializer=_init) as pool:
        for c in ("local", "perfectF9", "perfect", "other"):
            idx = np.flatnonzero(lab == c)
            for k in (0, 1, 2, 3, 4):
                jobs = [(11311 + 1000 * k + 100 * ["local", "perfectF9", "perfect", "other"].index(c) + j, idx, k, per_cell // 12) for j in range(12)]
                v = sum(pool.map(_job, jobs))
                n = 12 * (per_cell // 12)
                cells[f"{c}:k={k}"] = dict(violating=v, samples=n, fraction=v / n)
                print(c, k, cells[f"{c}:k={k}"], flush=True)
    res["cells"] = cells
    res["perfect_one_magic_leg_never"] = cells["perfect:k=1"]["violating"] == 0 and cells["perfectF9:k=1"]["violating"] == 0
    res["perfect_three_magic_legs_always"] = cells["perfect:k=3"]["violating"] == cells["perfect:k=3"]["samples"]
    res["F9_perfect_three_magic_legs_always"] = cells["perfectF9:k=3"]["violating"] == cells["perfectF9:k=3"]["samples"]
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
