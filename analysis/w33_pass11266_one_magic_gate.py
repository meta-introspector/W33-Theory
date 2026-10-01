"""Pass 11266: how often does a single magic (cubic) gate break substrate time reversal?  Exact counts.

Decisions use the exact criterion of Pass 11252 (c(Lp) = mu omega^{<b,Lp>} c(p), L anti-symplectic).

  * ONE QUTRIT, k = 1: all 216 Cliffords C (mod phase), U = C T.  Exact rule found: C T violates iff the symplectic part
    of C is a shear [[1,0],[c,1]] (C fixes Z: a diagonal Clifford times a Pauli) AND the Pauli part shifts X (a != 0).
    Count 3 x 6 = 18, probability exactly 1/12.  The parity-twisted shears [[2,0],[c,2]] never violate.
  * ONE QUTRIT, k = 2, 3: words C_k T ... C_1 T (C_0 is absorbed by conjugation), all 216^k.
  * TWO QUTRITS, k = 1: 200,000 Cliffords drawn exactly uniformly from the 51840 x 81 cosets (the exhaustive
    count needs ~37 CPU-hours at ~32 ms per decision and is left as a follow-up).
  * THREE QUTRITS, k = 1: uniform-ish random Cliffords (long random words), 2000 samples.  (Four qutrits: the
    decision is exact in principle, but the vectorised search needs ~1 GB per level at n = 4; not run here.)
"""

from __future__ import annotations

import itertools
import json
import sys
from fractions import Fraction
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11213_cubic_t_violation as P1  # noqa: E402
import w33_pass11227_two_qutrit_cp_mixing as P2  # noqa: E402
import w33_pass11252_exact_reversibility as R  # noqa: E402

OUT = ROOT / "data" / "w33_pass11266_one_magic_gate.json"


def _init(n):
    R.WEYL[n] = R.Weyl(n)


def _one_qutrit_word(idx):
    Cf = _one_qutrit_word.Cf
    U = np.eye(3, dtype=complex)
    for i in idx:
        U = Cf[i] @ P1.T1 @ U
    return R.decide(U, 1, _one_qutrit_word.rng)[0] is False


def _init1():
    R.WEYL[1] = R.Weyl(1)
    _one_qutrit_word.Cf = P1.clifford1()
    _one_qutrit_word.rng = np.random.default_rng(0)


def _two_qutrit_chunk(seed_count):
    seed, count = seed_count
    reps = _two_qutrit_chunk.reps
    T1 = np.kron(P2.T, P2.I3)
    rng = np.random.default_rng(seed)
    v = 0
    for _ in range(count):
        C = P2.PA[rng.integers(81)] @ reps[rng.integers(len(reps))]    # exactly uniform over the Clifford group mod phase
        if R.decide(C @ T1, 2, rng)[0] is False:
            v += 1
    return v


def _init2():
    R.WEYL[2] = R.Weyl(2)
    _two_qutrit_chunk.reps = np.load(P2.CACHE)


def one_qutrit_rule():
    wl = R.Weyl(1)
    R.WEYL[1] = wl
    rng = np.random.default_rng(0)
    Cf = P1.clifford1()
    out = []
    for C in Cf:
        M = np.zeros((2, 2), int)
        for j in range(2):
            e = np.zeros(2, int)
            e[j] = 1
            img = C @ wl.W[wl.index(e)] @ C.conj().T
            M[:, j] = wl.labels[int(np.argmax(np.abs(np.einsum('pij,ij->p', wl.W.conj(), img))))]
        viol = R.decide(C @ P1.T1, 1, rng)[0] is False
        shear = M[0, 0] == 1 and M[0, 1] == 0 and M[1, 1] == 1
        # X-shift of the Pauli part: C X^0 ... use the action on the computational basis: C|0> support
        xshift = int(np.argmax(np.abs(C[:, 0]))) != 0 if shear else None
        out.append((viol, bool(shear), xshift))
    n_v = sum(v for v, _, _ in out)
    rule = all(v == (s and bool(x)) for v, s, x in out)
    return dict(violating=n_v, total=len(out), probability=str(Fraction(n_v, len(out))), rule_holds=rule)


def run():
    res = dict(pass_id=11266)
    res["one_qutrit_k1"] = one_qutrit_rule()
    print(res["one_qutrit_k1"], flush=True)
    with Pool(11, initializer=_init1) as pool:
        for k in (2, 3):
            words = list(itertools.product(range(216), repeat=k))
            v = sum(pool.map(_one_qutrit_word, words, chunksize=2000))
            res[f"one_qutrit_k{k}"] = dict(violating=v, total=len(words), probability=v / len(words))
            print(k, res[f"one_qutrit_k{k}"], flush=True)
    with Pool(11, initializer=_init2) as pool:
        N2 = 200_000
        jobs = [(1000 + i, 500) for i in range(N2 // 500)]
        v = sum(pool.map(_two_qutrit_chunk, jobs))
        p2 = v / N2
        res["two_qutrit_k1"] = dict(violating=v, samples=N2, probability=p2, stderr=float(np.sqrt(p2 * (1 - p2) / N2)),
                                    sampling="exactly uniform over the 51840 x 81 Clifford cosets",
                                    exhaustive_count="not run: ~37 CPU-hours at ~32 ms per decision")
        print(res["two_qutrit_k1"], flush=True)
    for n, N in ((3, 2000),):
        R.WEYL[n] = R.Weyl(n)
        rng = np.random.default_rng(11266 + n)
        v = 0
        for _ in range(N):
            C = R.random_clifford(n, rng, length=40 * n)
            U = C @ R.local(n, 0, P2.T)
            if R.decide(U, n, rng)[0] is False:
                v += 1
        p = v / N
        res[f"n{n}_k1_sampled"] = dict(violating=v, samples=N, probability=p, stderr=float(np.sqrt(p * (1 - p) / N)))
        print(n, res[f"n{n}_k1_sampled"], flush=True)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
