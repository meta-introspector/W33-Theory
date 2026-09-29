#!/usr/bin/env python3
"""Pass 11166: the paper's spectral clock tick is a Clifford gate -- the free propagator of the light-cone quadric -- and
no translation-invariant clock can be a perfect (maximally scrambling) gate.

The null-history clock (analysis/w33_20260924_null_history_spectral_clock.py): 27 events F_3^3, edges along the 8 null
vectors of q(v) = v0 v2 - v1^2, Laplacian L, tick U = exp(-2 pi i L / 9), U^3 = I.  Treat the 27 events as three qutrits
(one per coordinate).  Results:
  * U is a Clifford gate: it sends every Pauli X^x Z^z to a single Pauli (checked on all 729), and its symplectic matrix is
        X_a -> X_a (all shifts conserved: the clock is translation invariant),
        Z0 -> Z0 X2,  Z1 -> Z1 X1,  Z2 -> Z2 X0     (up to phases),
    i.e. S = [[I, 0], [M, I]] in (x | z) form with M the Gram matrix of the quadric's polar form
    B(u, v) = q(u+v) - q(u) - q(v) = u0 v2 + u2 v0 + u1 v1 (M_02 = M_20 = M_11 = 1):
    the tick is the Gaussian propagator exp(-i q(p)) of the Lorentzian form, coupling the two light-cone coordinates and
    acting locally on the transverse one.
  * Party blocks: diagonal det 1 (rank 2), (0,2) and (2,0) rank 1, all others 0 -- the column law sum_i det S_ij = 1 holds,
    and the tick is NOT perfect: it exchanges one trit between the light-cone qutrits, none with the transverse one.
  * Theorem (general): any Clifford that commutes with all translations X_a has S = [[I, 0], [M, I]] with M symmetric, so
    every off-diagonal party block is [[0, 0], [m, 0]] (rank <= 1, det 0).  A translation-invariant (spectral, Cayley-graph)
    clock is never perfect: momentum conservation caps the scrambling at one trit per pair per tick.  Checked for all
    3^6 = 729 symmetric M.
  * Information flow (brute-force Choi entropies, cross-checked with the stabilizer formula of Pass 11165): the transverse
    qutrit's past stays in its own future (I = 2); each light-cone qutrit's past is split one trit in its own future, and
    the full two trits sit in the light-cone pair.  So the clock does not secret-share (contrast Pass 11165): forgetting a
    partner after a clock tick erases only one trit, not two.
"""
from __future__ import annotations

import itertools
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_20260924_null_history_spectral_clock as C  # noqa: E402
import w33_pass11165_three_qutrit_arrow as A3  # noqa: E402

OUT = ROOT / "data" / "w33_pass11166_clock_tick_clifford.json"
W = np.exp(2j * np.pi / 3)


def clock():
    A = C.adjacency()
    L = 8 * np.eye(27) - A
    ev, V = np.linalg.eigh(L)
    return V @ np.diag(np.exp(-1j * 2 * math.pi / 9 * ev)) @ V.conj().T


def pauli(x, z):
    """X^x Z^z on F_3^3 in the ordering of C.V (x, z in F_3^3)"""
    P = np.zeros((27, 27), complex)
    for s in C.V:
        t = tuple((s[i] + x[i]) % 3 for i in range(3))
        P[C.VID[t], C.VID[s]] = W ** (sum(z[i] * s[i] for i in range(3)) % 3)
    return P


PAULIS = {(x, z): pauli(x, z) for x in itertools.product(range(3), repeat=3) for z in itertools.product(range(3), repeat=3)}
KEYS = list(PAULIS)
STACK = np.array([PAULIS[k] for k in KEYS])


def image(U, x, z):
    """the unique Pauli (up to phase) equal to U P U^dagger, or None if U P U^dagger is not a single Pauli"""
    Q = U @ PAULIS[(x, z)] @ U.conj().T
    ov = np.abs(np.einsum('kij,ij->k', STACK.conj(), Q))
    hits = np.flatnonzero(np.abs(ov - 27) < 1e-6)
    return KEYS[hits[0]] if len(hits) == 1 else None


def symplectic(U):
    """columns = images of the generators in per-qutrit ordering (x0, z0, x1, z1, x2, z2)"""
    S = np.zeros((6, 6), int)
    for q in range(3):
        for t in range(2):
            x, z = [0] * 3, [0] * 3
            (x if t == 0 else z)[q] = 1
            ix, iz = image(U, tuple(x), tuple(z))
            col = [ix[0], iz[0], ix[1], iz[1], ix[2], iz[2]]
            S[:, 2 * q + t] = col
    return S % 3


def choi_infos(U):
    """brute force: I(in_a : Y) for Y subsets of outputs, from the Choi state of U (legs in0 in1 in2 out0 out1 out2)"""
    psi = np.zeros((27, 27), complex)
    for i in range(27):
        psi[i] = U[:, i]
    psi = psi.reshape([3] * 6) / math.sqrt(27)            # legs: in0 in1 in2 out0 out1 out2

    def S(keep):
        rest = [i for i in range(6) if i not in keep]
        M = psi.transpose(list(keep) + rest).reshape(3 ** len(keep), -1)
        sv = np.linalg.svd(M, compute_uv=False) ** 2
        sv = sv[sv > 1e-12]
        return float(-(sv * np.log(sv)).sum() / np.log(3))
    out = {}
    for a in range(3):
        for k in range(1, 4):
            for Y in itertools.combinations(range(3, 6), k):
                out[f"in{a}|out{tuple(y - 3 for y in Y)}"] = round(S((a,)) + S(Y) - S((a,) + Y), 6)
    return out


def blocks(S):
    det = [[int(S[2 * i, 2 * j] * S[2 * i + 1, 2 * j + 1] - S[2 * i, 2 * j + 1] * S[2 * i + 1, 2 * j]) % 3 for j in range(3)] for i in range(3)]
    rk = [[A3.G2.rank3(S[2 * i:2 * i + 2, 2 * j:2 * j + 2]) for j in range(3)] for i in range(3)]
    return det, rk


def translation_invariant_theorem():
    """all symmetric M over F3: S = [[I,0],[M,I]] (x|z form) -> per-party blocks; max off-diagonal rank, any perfect?"""
    max_off, perfect = 0, 0
    for m in itertools.product(range(3), repeat=6):
        M = np.array([[m[0], m[1], m[2]], [m[1], m[3], m[4]], [m[2], m[4], m[5]]])
        S = np.zeros((6, 6), int)
        for i in range(3):
            S[2 * i, 2 * i] = S[2 * i + 1, 2 * i + 1] = 1
            for j in range(3):
                S[2 * i + 1, 2 * j] = (S[2 * i + 1, 2 * j] + M[i, j]) % 3   # x_j -> z_i component (either convention: one entry)
        det, rk = blocks(S)
        max_off = max(max_off, max(rk[i][j] for i in range(3) for j in range(3) if i != j))
        perfect += all(d != 0 for row in det for d in row)
    return max_off, perfect


def summarize():
    U = clock()
    order3 = float(np.abs(np.linalg.matrix_power(U, 3) - np.eye(27)).max())
    clifford = all(image(U, x, z) is not None for (x, z) in PAULIS)
    S = symplectic(U)
    J = np.zeros((6, 6), int)
    for k in range(3):
        J[2 * k, 2 * k + 1], J[2 * k + 1, 2 * k] = 1, -1
    det, rk = blocks(S)
    shifts_conserved = all(image(U, x, (0, 0, 0)) == (x, (0, 0, 0)) for x in itertools.product(range(3), repeat=3))
    brute = choi_infos(U)
    stab = A3.infos(S)
    max_off, n_perfect = translation_invariant_theorem()
    res = dict(pass_id=11166, U_cubed_minus_I=order3, is_clifford=clifford, symplectic=S.tolist(),
               symplectic_ok=bool(np.array_equal((S.T @ J @ S) % 3, J % 3)), shifts_conserved=shifts_conserved,
               block_dets=det, block_ranks=rk, column_law=[sum(det[i][j] for i in range(3)) % 3 for j in range(3)],
               perfect=all(d != 0 for row in det for d in row),
               infos_bruteforce=brute, infos_match_stabilizer_formula=all(abs(brute[k] - stab[k]) < 1e-6 for k in brute),
               translation_invariant_max_offdiag_rank=max_off, translation_invariant_perfect=n_perfect)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    r = summarize()
    print(json.dumps({k: v for k, v in r.items() if k != 'infos_bruteforce'}, indent=1))
    print({k: v for k, v in r['infos_bruteforce'].items() if v})
