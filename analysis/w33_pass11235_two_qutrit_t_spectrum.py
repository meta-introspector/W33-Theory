"""Pass 11235: random two-qutrit Clifford+cubic dynamics -- how often T breaks, and is 2pi/9 the smallest quantum?

Pass 11227 found, by exact search over the anti-unitary two-qutrit Clifford group, that (T(x)T) SUM breaks substrate
time-reversal symmetry, that single-sector placements do not, and that every violating circuit tested had best
time-reversal fidelity F_T = (1 + 2cos 2pi/9)/3 = 0.8440296..., the one-qutrit minimal value (Pass 11228).

Here: random words W = C_d T_{q_d} ... C_1 T_{q_1} C_0, with C_i uniformly random two-qutrit Cliffords (one of the 51840
symplectic representatives times a random Pauli) and the cubic phase T on a random qutrit, d = 1..4 cubic gates.  For
each word: the exact verdict and the exact best fidelity
    F_T(U) = max over 51840 x 81 anti-unitary Cliffords of |tr(V U^* V^dag U)| / 9,
computed fast as max_{r,a} |tr(B_r Y_a)| / 9 with B_r = C_r U^* C_r^dag and Y_a = P_a^dag U P_a.
Question: is any violation milder than the 2pi/9 quantum (F_T > 0.8440296 but < 1)?
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11227_two_qutrit_cp_mixing as P2  # noqa: E402

OUT = ROOT / "data" / "w33_pass11235_two_qutrit_t_spectrum.json"
QUANTUM = (1 + 2 * np.cos(2 * np.pi / 9)) / 3


def fidelity_fast(U, reps, chunk=4000):
    Pa = P2.PA
    Y = np.einsum('aji,jk,akl->ail', Pa.conj(), U, Pa)         # P_a^dag U P_a
    best = 0.0
    for s in range(0, len(reps), chunk):
        C = reps[s:s + chunk]
        B = np.einsum('rij,jk,rlk->ril', C, U.conj(), C.conj())
        tr = np.einsum('rjk,akj->ra', B, Y)
        best = max(best, float(np.abs(tr).max() / 9))
    return best


def random_word(rng, reps, d):
    T1 = np.kron(P2.T, P2.I3)
    T2 = np.kron(P2.I3, P2.T)

    def cliff():
        return P2.PA[rng.integers(81)] @ reps[rng.integers(len(reps))]
    U = cliff()
    for _ in range(d):
        U = cliff() @ (T1 if rng.integers(2) == 0 else T2) @ U
    return U


def run(per_depth=60, seed=11235):
    reps = np.load(P2.CACHE)
    rng = np.random.default_rng(seed)
    res = dict(pass_id=11235, quantum=QUANTUM)
    # controls against Pass 11227
    res["control_TT_SUM"] = fidelity_fast(np.kron(P2.T, P2.T) @ P2.SUM, reps)
    res["control_SUM"] = fidelity_fast(P2.SUM, reps)
    by = {}
    allf = []
    for d in range(1, 5):
        fs = [fidelity_fast(random_word(rng, reps, d), reps) for _ in range(per_depth)]
        viol = [f for f in fs if f < 1 - 1e-9]
        allf += viol
        by[str(d)] = dict(words=per_depth, violating=len(viol),
                          fidelities=dict(Counter(round(f, 9) for f in viol)))
        print(d, by[str(d)], flush=True)
    res["by_depth"] = by
    res["max_violator_fidelity"] = max(allf) if allf else None
    res["milder_than_quantum"] = sum(1 for f in allf if f > QUANTUM + 1e-9)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
