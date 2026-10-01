"""Pass 11268: fixed perfect (maximally scrambling) magic ticks -- which dressings break substrate time reversal?

Perfect ticks (perfect tensors / AME(2n,3) as gates) exist for 1, 2, 3 and 5 qutrits (Pass 11170).  Perfectness is
invariant under local unitaries on every leg, so dressing a perfect Clifford tick V with local cubic phases,
    U(a, b) = (T^{a_1} (x) ... (x) T^{a_n}) V (T^{b_1} (x) ... (x) T^{b_n}),   a_i, b_i in Z_9,
gives perfect MAGIC ticks (T^3 = Z is Clifford, so a_i mod 3 != 0 carries the magic).  Each U is decided exactly with
Pass 11252's criterion: all 9^4 = 6561 dressings of a perfect two-qutrit tick, and a sample of the 9^6 dressings of a
perfect three-qutrit tick.  The perfect Clifford ticks are built from Pass 11170's perfect F_9-linear circulants
(embedded symplectic matrix -> Weil twirl), and perfectness is verified numerically on every n|n cut of the 2n legs.
Baseline: Pass 11252's census of generic Clifford+T words.
"""

from __future__ import annotations

import itertools
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11170_perfect_gate_ladder as PL  # noqa: E402
import w33_pass11227_two_qutrit_cp_mixing as P2  # noqa: E402
import w33_pass11252_exact_reversibility as R  # noqa: E402

OUT = ROOT / "data" / "w33_pass11268_perfect_magic_ticks.json"
T = P2.T


def perfect_symplectic(m):
    hits = PL.circulant_count(m)
    for idx in hits:
        S = PL.embed(PL.circulant(PL.F9[list(idx)]))
        wl = R.WEYL[m]
        if ((S.T @ wl.Om @ S - wl.Om) % 3 == 0).all():
            return S, idx
    raise RuntimeError("no symplectic embedding in this convention")


def is_perfect(U, n):
    """every n|n cut of the 2n legs (n output + n input legs) gives a matrix proportional to a unitary"""
    D = 3 ** n
    Tn = U.reshape([3] * (2 * n))
    legs = list(range(2 * n))
    for X in itertools.combinations(legs, n):
        Y = [l for l in legs if l not in X]
        M = np.transpose(Tn, list(X) + Y).reshape(D, D)
        G = M @ M.conj().T
        if not np.allclose(G, G[0, 0] * np.eye(D), atol=1e-8):
            return False
    return True


def local_T(n, exps):
    M = np.array([[1.0 + 0j]])
    for e in exps:
        M = np.kron(M, np.linalg.matrix_power(T, int(e)))
    return M


_W = {}


def _init3(Vr, Vi):
    R.WEYL[3] = R.Weyl(3)
    _W["V"] = Vr + 1j * Vi


def _decide3(e):
    U = local_T(3, e[:3]) @ _W["V"] @ local_T(3, e[3:])
    return R.decide(U, 3, np.random.default_rng(sum(e)))[0] is False


def exhaustive_edges(V):
    """n = 3: every dressing with exactly one magic leg, and every dressing with exactly five magic legs"""
    from multiprocessing import Pool
    magic = [1, 2, 4, 5, 7, 8]
    clif = [0, 3, 6]
    one, five = [], []
    for leg in range(6):
        for m in magic:
            for rest in itertools.product(clif, repeat=5):
                e = list(rest)
                e.insert(leg, m)
                one.append(tuple(e))
        for c in clif:
            for rest in itertools.product(magic, repeat=5):
                e = list(rest)
                e.insert(leg, c)
                five.append(tuple(e))
    with Pool(11, initializer=_init3, initargs=(V.real.copy(), V.imag.copy())) as pool:
        v1 = sum(pool.map(_decide3, one, chunksize=200))
        v5 = sum(pool.map(_decide3, five, chunksize=500))
    return dict(one_magic_leg=[v1, len(one)], five_magic_legs=[v5, len(five)])


def run():
    res = dict(pass_id=11268)
    rng = np.random.default_rng(11268)
    out = {}
    for n, mode in ((2, "all"), (3, 1500)):
        R.WEYL[n] = R.Weyl(n)
        S, idx = perfect_symplectic(n)
        V = R.weil(R.WEYL[n], S, rng)
        perf = is_perfect(V, n)
        assert perf, "constructed tick is not perfect"
        if mode == "all":
            grid = list(itertools.product(range(9), repeat=2 * n))
        else:
            grid = [tuple(rng.integers(9, size=2 * n)) for _ in range(mode)]
        verdict = Counter()
        by_magic = Counter()
        tot_magic = Counter()
        perfect_checked = 0
        for e in grid:
            U = local_T(n, e[:n]) @ V @ local_T(n, e[n:])
            if perfect_checked < 20:
                assert is_perfect(U, n)
                perfect_checked += 1
            v = R.decide(U, n, rng)[0]
            verdict[v] += 1
            k = sum(1 for x in e if x % 3)                    # number of legs carrying magic
            tot_magic[k] += 1
            by_magic[k] += (v is False)
        out[str(n)] = dict(circulant=list(map(int, idx)), dressings=len(grid), exhaustive=(mode == "all"),
                           violating=verdict[False], reversible=verdict[True], undecided=verdict[None],
                           violating_by_magic_legs={str(k): [by_magic[k], tot_magic[k]] for k in sorted(tot_magic)},
                           perfectness_rechecked_on=perfect_checked)
        print(n, out[str(n)], flush=True)
        if n == 3:
            out["3"]["exhaustive_edges"] = exhaustive_edges(V)
            print(out["3"]["exhaustive_edges"], flush=True)
    res["perfect_magic"] = out
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
