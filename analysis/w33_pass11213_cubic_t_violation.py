"""Pass 11213: T-violation enters exactly with the cubic (non-Clifford) layer.

Pass 11206: every Clifford tick is inverted by a substrate time reversal (an anti-unitary Clifford Theta = V K), as
Wonenburger's theorem implies -- the Clifford skeleton has no T-violation.  Every unitary is inverted by SOME
anti-unitary (U^* and U^dagger are unitarily similar), so the physical question is whether the substrate's own time
reversals -- the 216 anti-unitary Cliffords V K of one qutrit (V in the Clifford group mod phases) -- still suffice once
the non-Clifford resource is switched on.  The resource is the cubic phase gate T = diag(1, zeta9, zeta9^8) =
zeta9^{x^3}, the qutrit form of the degree-3 (E6 cubic) layer of the substrate's universal gate set.

Measured:
  * single qutrit: the fraction of random Clifford+T words (depth d = number of T gates) that NO substrate time
    reversal inverts, for d = 0..8; the shortest T-violating words;
  * a Jarlskog-type witness: J(U) = Im tr(U Z U^dagger X U Z^dagger U^dagger X^dagger), odd under every substrate
    time reversal that fixes the Pauli frame up to Clifford -- non-zero J on the violating words;
  * two qutrits (SUM entangler): the same fraction for short words, searching the anti-unitary two-qutrit Clifford
    group (51840 x 81 elements mod phase).
"""

from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "w33_pass11213_cubic_t_violation.json"
W = np.exp(2j * np.pi / 3)
Z9 = np.exp(2j * np.pi / 9)
H1 = np.array([[W ** (j * k) for k in range(3)] for j in range(3)]) / np.sqrt(3)
S1 = np.diag([1, 1, W])
X1 = np.roll(np.eye(3), 1, axis=0)
Z1 = np.diag([1, W, W ** 2])
T1 = np.diag([1, Z9, Z9 ** 8])


def canon(U):
    """projective normal form: divide by the phase of the first entry of largest modulus"""
    f = U.ravel()
    k = int(np.argmax(np.abs(f) > 1e-9))
    return U * (abs(f[k]) / f[k])


def key(U):
    f = canon(U).ravel()
    return tuple(np.round(f.real, 5) + 0.0) + tuple(np.round(f.imag, 5) + 0.0)     # +0.0 folds -0.0 into 0.0


def clifford1():
    seen = {key(np.eye(3)): np.eye(3, dtype=complex)}
    frontier = [np.eye(3, dtype=complex)]
    gens = [H1, S1, X1, Z1]
    while frontier:
        new = []
        for U in frontier:
            for g in gens:
                V = g @ U
                k = key(V)
                if k not in seen:
                    seen[k] = canon(V)
                    new.append(canon(V))
        frontier = new
    return np.array(list(seen.values()))


def t_invertible(U, C):
    """is there V in C with V U^* V^dagger = lambda U^dagger ?"""
    M = np.einsum('vij,jk,vlk,lm->vim', C, U.conj(), C.conj(), U)       # V U^* V^dagger U
    tr = np.einsum('vii->v', M) / 3
    dev = np.abs(M - tr[:, None, None] * np.eye(3)).max(axis=(1, 2))
    return bool((dev < 1e-8).any())


def jarlskog(U):
    return float(np.imag(np.trace(U @ Z1 @ U.conj().T @ X1 @ U @ Z1.conj().T @ U.conj().T @ X1.conj().T)))


def random_word(rng, C, d):
    U = C[rng.integers(len(C))]
    for _ in range(d):
        U = C[rng.integers(len(C))] @ T1 @ U
    return U


def shortest_violators(C, max_d=2):
    out = []
    Cl = list(C)
    for d in range(1, max_d + 1):
        found = 0
        for combo in itertools.product(range(len(Cl)), repeat=d):
            U = np.eye(3, dtype=complex)
            for c in combo:
                U = T1 @ Cl[c] @ U
            if not t_invertible(U, C):
                found += 1
                if len(out) < 3:
                    out.append(dict(depth=d, J=round(jarlskog(U), 8)))
        if found:
            return d, found, len(Cl) ** d, out
    return None, 0, 0, out


def run(samples=3000, seed=11213):
    rng = np.random.default_rng(seed)
    C = clifford1()
    res = dict(pass_id=11213, clifford1_order_mod_phase=len(C))
    assert len(C) == 216
    # controls: every Clifford is invertible (Pass 11206), T alone is (K inverts diagonal phases)
    res["all_cliffords_T_invertible"] = all(t_invertible(V, C) for V in C)
    res["T_gate_alone_T_invertible"] = t_invertible(T1, C)
    frac = {}
    for d in range(0, 9):
        viol, jnz, jnz_inv = 0, 0, 0
        for _ in range(samples if d else 200):
            U = random_word(rng, C, d)
            inv = t_invertible(U, C)
            viol += not inv
            j = abs(jarlskog(U)) > 1e-9
            if not inv:
                jnz += j
            else:
                jnz_inv += j
        n = samples if d else 200
        frac[d] = dict(samples=n, violating=viol, fraction=viol / n, J_nonzero_among_violating=jnz,
                       J_nonzero_among_invertible=jnz_inv)
    res["by_depth"] = {str(k): v for k, v in frac.items()}
    d, found, total, ex = shortest_violators(C, 1)
    res["shortest"] = dict(depth=d, violating_words=found, words=total, examples=ex)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
