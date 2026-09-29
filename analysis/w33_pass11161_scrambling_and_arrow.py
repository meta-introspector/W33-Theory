#!/usr/bin/env python3
"""Pass 11161: perfect gates are the maximal scramblers, the determinant law forbids information from ever being
concentrated (I3 <= 0 for every two-qutrit Clifford gate), and a maximally entangling tick followed by forgetting the
partner is exactly the paper's arrow of time (Theorem 4.7: two trits).

For S in Sp(4,3) (blocks S_AA, S_AB, S_BA, S_BB), the Choi state's mutual informations (units of log 3) follow from the
stabilizer rank formula of Pass 11156:
    I(A_in : A_out) = 2 - rank S_BA,     I(A_in : B_out) = 2 - rank S_AA,     I(A_in : A_out B_out) = 2,
so the tripartite information is
    I3(A_in : A_out : B_out) = 2 - rank S_AA - rank S_BA.
With det S_AA + det S_BA = 1 (Pass 11156) the pair (rank S_AA, rank S_BA) is one of (2,0), (0,2), (2,1), (1,2), (2,2):
    I3 = 0  for 1152 gates (local and swap-type),  I3 = -1 for 36864,  I3 = -2 (the minimum) for the 13824 perfect ones.
So I3 <= 0 for EVERY two-qutrit Clifford gate -- information about A's past is never concentrated in A's or B's future
alone beyond what the pair holds -- and the perfect gates are exactly the maximal scramblers.  Checked against
brute-force entropies of actual Choi states.
THE ARROW: if S_BA is invertible (rank 2 -- every perfect gate), then with the partner B maximally mixed the channel on A,
    rho -> Tr_B[ U (rho (x) I/3) U^dagger ],
is completely depolarising (checked numerically): one tick with a partner, then forgetting the partner (a maximally mixed
qutrit together with its purification -- a 9-dimensional environment, two trits), reproduces exactly the two-trit
forgetting of the paper's Theorem 4.7 (Pauli twirl = complete depolarisation, Landauer cost 2 k_B T ln 3).  The arrow of
time appears as soon as a maximally entangling tick is followed by discarding what it entangled with.
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11156_perfect_spacetime_gates as G  # noqa: E402

OUT = ROOT / "data" / "w33_pass11161_scrambling_and_arrow.json"


def choi_mutual_infos(U):
    omega = np.zeros(9); omega[[0, 4, 8]] = 1
    om2 = np.kron(omega, omega).reshape(3, 3, 3, 3).transpose(0, 2, 1, 3).reshape(81) / 3
    psi = (np.kron(np.eye(9), U) @ om2).reshape(3, 3, 3, 3)          # Ai, Bi, Ao, Bo

    def S(keep):
        rest = [i for i in range(4) if i not in keep]
        M = psi.transpose(list(keep) + rest).reshape(3 ** len(keep), -1)
        sv = np.linalg.svd(M, compute_uv=False) ** 2
        sv = sv[sv > 1e-12]
        return float(-(sv * np.log(sv)).sum() / np.log(3))
    I_aa = S((0,)) + S((2,)) - S((0, 2))
    I_ab = S((0,)) + S((3,)) - S((0, 3))
    I_a_ab = S((0,)) + S((2, 3)) - S((0, 2, 3))
    return I_aa, I_ab, I_a_ab, I_aa + I_ab - I_a_ab


def channel_on_A(U, rho):
    big = U @ np.kron(rho, np.eye(3) / 3) @ U.conj().T
    return big.reshape(3, 3, 3, 3).trace(axis1=1, axis2=3)


def summarize():
    Gp = G.sp43()
    census = Counter(2 - G.rank3(S[:2, :2]) - G.rank3(S[2:, :2]) for S in Gp)
    checks, depol = [], []
    rng = np.random.default_rng(3)
    for U in G.two_qutrit_cliffords(14, seed=5):
        S = G.symplectic_of(U)
        Iaa, Iab, Iaab, I3 = choi_mutual_infos(U)
        rAA, rBA = G.rank3(S[:2, :2]), G.rank3(S[2:, :2])
        checks.append(abs(Iaa - (2 - rBA)) < 1e-6 and abs(Iab - (2 - rAA)) < 1e-6 and abs(Iaab - 2) < 1e-6
                      and abs(I3 - (2 - rAA - rBA)) < 1e-6)
        if rBA == 2:
            X = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)); rho = X @ X.conj().T; rho /= np.trace(rho)
            depol.append(float(np.abs(channel_on_A(U, rho) - np.eye(3) / 3).max()))
    res = dict(pass_id=11161, I3_census={str(k): v for k, v in sorted(census.items())}, I3_never_positive=max(census) <= 0,
               formulas_checked=all(checks), n_checked=len(checks), depolarising_cases=len(depol),
               depolarising_max_dev=max(depol) if depol else None)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
