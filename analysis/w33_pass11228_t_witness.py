"""Pass 11228: a complete measure of T-violation on the substrate, and why polynomial Pauli witnesses fail.

Pass 11213 found that the cubic layer breaks substrate time-reversal symmetry (no anti-unitary Clifford Theta = V K
inverts the tick).  This pass asks for a witness.

1. The natural cubic invariant J3(U) = Im sum_{P,Q} tr(U P U^dag Q)^3 is invariant under Cliffords on both sides,
   and satisfies J3(U^-1) = J3(U), J3(U^*) = -J3(U), so it vanishes on every reversible tick -- but it vanishes
   IDENTICALLY: pairing P with P^dagger makes the sum real.  Any T-witness built from the Pauli transfer matrix must be
   odd under exchanging its rows and columns (U <-> U^dagger).
2. Transpose-odd magnitude moments (sorted row vs column moments of |R(P,Q)|^4, ^6) are sound witnesses (zero on all
   reversible ticks) but incomplete: they certify only part of the violators and miss the minimal one T X.
3. A complete, operational measure: the best time-reversal fidelity
       F_T(U) = max_{V in Clifford} |tr(V U^* V^dag U)| / d,
   which equals 1 iff some substrate time reversal inverts U (V U^* V^dag U is unitary; |tr| = d iff it is scalar).
   1 - F_T quantifies T-violation; computed exactly for the 18 minimal violators and its distribution by depth.
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11213_cubic_t_violation as P  # noqa: E402

OUT = ROOT / "data" / "w33_pass11228_t_witness.json"
PAULIS = [np.linalg.matrix_power(P.X1, a) @ np.linalg.matrix_power(P.Z1, b) for a in range(3) for b in range(3)]


def J3(U):
    s = 0
    for Pm in PAULIS:
        M = U @ Pm @ U.conj().T
        for Q in PAULIS:
            s += np.trace(M @ Q) ** 3
    return float(s.imag)


def moment_witness(U):
    R = np.array([[np.trace(U @ Pm @ U.conj().T @ Q.conj().T) / 3 for Q in PAULIS] for Pm in PAULIS])
    A = np.abs(R) ** 2
    out = [np.abs(np.sort((A ** k).sum(1)) - np.sort((A ** k).sum(0))).max() for k in (2, 3)]
    return float(max(out))


def fidelity(U, C):
    M = np.einsum('vij,jk,vlk,lm->vim', C, U.conj(), C.conj(), U)
    return float(np.abs(np.einsum('vii->v', M)).max() / 3)


def run(samples=1500, seed=11228):
    C = P.clifford1()
    rng = np.random.default_rng(seed)
    res = dict(pass_id=11228)
    # 1. J3 identically zero; symmetries
    Us = [P.random_word(rng, C, d) for d in range(1, 7) for _ in range(40)]
    res["J3_max_abs_on_240_words"] = max(abs(J3(U)) for U in Us)
    U = Us[7]
    res["J3_symmetries"] = dict(clifford_invariant=abs(J3(C[3] @ U @ C[11]) - J3(U)) < 1e-9,
                                inverse=abs(J3(U.conj().T) - J3(U)) < 1e-9, conjugate=abs(J3(U.conj()) + J3(U)) < 1e-9)
    # 2 + 3. witnesses vs exact reversibility, by depth
    table, fid = Counter(), {}
    for d in range(0, 9):
        fs = []
        for _ in range(samples if d else 100):
            W = P.random_word(rng, C, d)
            inv = P.t_invertible(W, C)
            f = fidelity(W, C)
            assert inv == (f > 1 - 1e-9)
            mw = moment_witness(W) > 1e-9
            table[("reversible" if inv else "violating", "moment!=0" if mw else "moment=0")] += 1
            if not inv:
                fs.append(f)
        fid[str(d)] = dict(violating=len(fs), mean_fidelity=float(np.mean(fs)) if fs else None,
                           min_fidelity=float(np.min(fs)) if fs else None)
    res["moment_witness_vs_exact"] = {f"{k[0]}|{k[1]}": v for k, v in table.items()}
    res["fidelity_equals_one_iff_reversible"] = True
    res["violator_fidelity_by_depth"] = fid
    # minimal violators T X^a Z^b S^c, a != 0
    mins = Counter()
    for a in (1, 2):
        for b in range(3):
            for c in range(3):
                V = np.linalg.matrix_power(P.X1, a) @ np.linalg.matrix_power(P.Z1, b) @ np.linalg.matrix_power(P.S1, c)
                mins[round(fidelity(P.T1 @ V, C), 12)] += 1
    res["minimal_violator_fidelities"] = {str(k): v for k, v in mins.items()}
    res["T_X_fidelity"] = fidelity(P.T1 @ P.X1, C)
    res["T_X_moment_witness"] = moment_witness(P.T1 @ P.X1)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
