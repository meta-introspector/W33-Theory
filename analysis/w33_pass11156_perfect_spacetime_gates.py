#!/usr/bin/env python3
"""Pass 11156: perfect two-qutrit gates -- the Clifford gates of W(3,3) that are maximally entangling across time, across
space AND across the space-time diagonal at once (Choi state = AME(4,3), a perfect tensor).

A two-qutrit Clifford gate U is a symplectic S in Sp(4,3) acting on Pauli labels (a, b) in F_3^2 + F_3^2 (form w + w,
the spatial W(3,3) of Pass 11143), with blocks S = [[S_AA, S_AB], [S_BA, S_BB]].  Its Choi state lives on four qutrits
A_in, B_in, A_out, B_out and is a stabilizer state with Lagrangian L = {(v, S v)}; for a subsystem X,
S(X) = |X| - dim L_X (in units of log 3), L_X the elements of L supported on X.  The three 2|2 cuts are:
    time    {A_in B_in | A_out B_out}:  entropy 2 always (unitarity);
    space   {A_in A_out | B_in B_out}:  entropy rank(S_BA)  (the gate's operator entanglement);
    diagonal{A_in B_out | B_in A_out}:  entropy rank(S_AA)  (A's past with B's future).
So U is PERFECT (its Choi state is absolutely maximally entangled, AME(4,3)) iff S_AA and S_BA are both invertible.
Checked against brute-force von Neumann entropies of the actual Choi states of random Clifford unitaries.
Results (exact, all 51840 elements of Sp(4,3); the formula checked against brute-force Choi entropies of 10 random
Clifford unitaries):
  * COMPLEMENTARITY LAW: det S_AA + det S_BA = 1 (mod 3) for every symplectic S (it is the symplectic condition on the
    first block column, since w(Mu, Mu') = det(M) w(u, u') on F_3^2).  Hence the space entanglement rank(S_BA) and the
    diagonal entanglement rank(S_AA) can never both be deficient: the census is (0,2): 576 local gates, (2,0): 576
    swap-type, (1,2): 18432, (2,1): 18432, (2,2): 13824 -- no gate is weak in both;
  * PERFECT gates: exactly those with det S_AA = det S_BA = -1: 13824 = 24^3 of 51840 (4/15); all four blocks are then
    invertible, so the Choi Lagrangian is an MDS code of length 4 over F_3^2 (perfect = AME(4,3) = MDS);
  * UNIQUENESS: the 13824 perfect gates form ONE orbit under local Clifford operations on both sides (local stabilizer of
    order 24) -- the Clifford-level instance of the known uniqueness of the four-qutrit AME state up to local unitaries
    (I. Tan, arXiv:2601.19677);
  * THE TETRACODE: for anti-unitary (time-reversing) Clifford operations, similitudes of multiplier -1, the law reads
    det S_AA + det S_BA = -1, and again 13824 are perfect; the tetracode gate K (x) I, K = [[1,1],[1,-1]] (the tetracode's
    redundancy matrix acting on phase-space labels), is one of them; composed with local time reversal
    D = diag(1,-1,1,-1) it becomes a symplectic, perfect -- i.e. the unique -- unitary perfect gate.
    That AME(4,3) comes from the tetracode / orthogonal array OA(9,4,3,2) is known (Goyeneche-Zyczkowski, Helwig et al.);
    new here is the W(3,3) reading: the repo's "every point of W(3,3) carries a tetracode" has a dynamical face -- the
    tetracode is the perfect space-time gate, up to a local time reversal, and the complementarity law above.
  * orders of perfect gates: 3, 4, 5, 6, 8, 9, 10, 12, 18 (never 1 or 2).
"""
from __future__ import annotations

import itertools
import json
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "w33_pass11156_perfect_spacetime_gates.json"
J = np.zeros((4, 4), int)
J[0, 1], J[1, 0], J[2, 3], J[3, 2] = 1, -1, 1, -1


def rank3(M):
    M = [list(map(int, r)) for r in np.array(M) % 3]
    r = 0
    rows, cols = len(M), len(M[0])
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i][c] % 3), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = 1 if M[r][c] % 3 == 1 else 2
        M[r] = [(x * inv) % 3 for x in M[r]]
        for i in range(rows):
            if i != r and M[i][c] % 3:
                f = M[i][c]
                M[i] = [(a - f * b) % 3 for a, b in zip(M[i], M[r])]
        r += 1
    return r


def sp43():
    gens = [np.array(g) % 3 for g in (
        [[1, 1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], [[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]],
        [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 1], [0, 0, 0, 1]], [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]],
        [[1, 0, 0, 0], [0, 1, 1, 0], [0, 0, 1, 0], [1, 0, 0, 1]])]
    seen = {tuple(np.eye(4, dtype=int).ravel())}
    frontier = [np.eye(4, dtype=int)]
    while frontier:
        nxt = []
        for g in frontier:
            for h in gens:
                m = (h @ g) % 3
                k = tuple(m.ravel())
                if k not in seen:
                    seen.add(k); nxt.append(m)
        frontier = nxt
    out = [np.array(k).reshape(4, 4) for k in seen]
    assert all(np.array_equal((g.T @ J @ g) % 3, J % 3) for g in out[:500])
    return out


def order(S):
    M, k = S.copy(), 1
    I = np.eye(4, dtype=int)
    while not np.array_equal(M % 3, I):
        M = (M @ S) % 3; k += 1
    return k


# ---- brute-force check of the entropy formula on actual Choi states -------------------------------------------------
W3 = np.exp(2j * np.pi / 3)
Xq = np.roll(np.eye(3), 1, axis=0)
Zq = np.diag([1, W3, W3 * W3])


def two_qutrit_cliffords(n=12, seed=0):
    """random products of the generators F(x)I, I(x)F, S(x)I, I(x)S, CSUM"""
    F = np.array([[W3 ** (j * k) for k in range(3)] for j in range(3)]) / np.sqrt(3)
    Sg = np.diag([1, 1, W3])
    I3 = np.eye(3)
    CSUM = np.zeros((9, 9))
    for a in range(3):
        for b in range(3):
            CSUM[a * 3 + (a + b) % 3, a * 3 + b] = 1
    gens = [np.kron(F, I3), np.kron(I3, F), np.kron(Sg, I3), np.kron(I3, Sg), CSUM]
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        U = np.eye(9, dtype=complex)
        for _ in range(rng.integers(3, 25)):
            U = gens[rng.integers(len(gens))] @ U
        out.append(U)
    return out


def pauli2(v):
    P = lambda x, z: np.linalg.matrix_power(Xq, x % 3) @ np.linalg.matrix_power(Zq, z % 3)
    return np.kron(P(v[0], v[1]), P(v[2], v[3]))


def symplectic_of(U):
    cols = []
    for e in np.eye(4, dtype=int):
        Pe = U @ pauli2(e) @ U.conj().T
        hit = [v for v in itertools.product(range(3), repeat=4) if abs(abs(np.trace(pauli2(v).conj().T @ Pe)) - 9) < 1e-6]
        cols.append(hit[0])
    return np.array(cols).T % 3


def choi_entropies(U):
    """von Neumann entropies (log 3 units) of the Choi state across the three 2|2 cuts; parties order Ai, Bi, Ao, Bo"""
    omega = np.zeros(9); omega[[0, 4, 8]] = 1
    om2 = np.kron(omega, omega).reshape(3, 3, 3, 3).transpose(0, 2, 1, 3).reshape(81) / 3   # |Ai Bi Ao' Bo'>
    psi = (np.kron(np.eye(9), U) @ om2).reshape(3, 3, 3, 3)                                    # Ai, Bi, Ao, Bo
    def S(keep):
        rest = [i for i in range(4) if i not in keep]
        M = psi.transpose(list(keep) + rest).reshape(9, 9)
        sv = np.linalg.svd(M, compute_uv=False) ** 2
        sv = sv[sv > 1e-12]
        return float(-(sv * np.log(sv)).sum() / np.log(3))
    return dict(time=S((0, 1)), space=S((0, 2)), diagonal=S((0, 3)))


def det2(M):
    return int(round(M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0])) % 3


def local_orbits(perfect):
    """double cosets L S L' of the perfect gates under local Cliffords (block-diagonal SL2 x SL2 on both sides)"""
    SL2 = [np.array(m).reshape(2, 2) for m in itertools.product(range(3), repeat=4) if (m[0] * m[3] - m[1] * m[2]) % 3 == 1]
    gens = []
    for g in SL2[:6]:
        for side in (0, 1):
            Lm = np.eye(4, dtype=int); Lm[2 * side:2 * side + 2, 2 * side:2 * side + 2] = g
            gens.append(Lm)
    key = lambda M: tuple((M % 3).ravel())
    remaining = {key(S): S for S in perfect}
    orbits = []
    while remaining:
        k0, S0 = next(iter(remaining.items()))
        orb = {k0}; frontier = [S0]
        while frontier:
            nxt = []
            for S in frontier:
                for g in gens:
                    for M in ((g @ S) % 3, (S @ g) % 3):
                        k = key(M)
                        if k not in orb:
                            orb.add(k); nxt.append(M)
            frontier = nxt
        orbits.append(len(orb))
        for k in orb:
            remaining.pop(k, None)
    return sorted(orbits, reverse=True)


def summarize():
    G = sp43()
    perfect, ranks = [], Counter()
    for S in G:
        rs, rd = rank3(S[2:, :2]), rank3(S[:2, :2])
        ranks[(rs, rd)] += 1
        if rs == 2 and rd == 2:
            perfect.append(S)
    # brute-force validation of the formula
    checks = []
    for U in two_qutrit_cliffords(10):
        S = symplectic_of(U)
        e = choi_entropies(U)
        checks.append(abs(e['space'] - rank3(S[2:, :2])) < 1e-6 and abs(e['diagonal'] - rank3(S[:2, :2])) < 1e-6
                      and abs(e['time'] - 2) < 1e-6)
    orders = Counter(order(S) for S in perfect)
    det_law = all((det2(S[:2, :2]) + det2(S[2:, :2])) % 3 == 1 for S in G)
    perfect_dets = Counter((det2(S[:2, :2]), det2(S[2:, :2])) for S in perfect)
    # multiplier -1 similitudes (anti-unitary Clifford operations): D = diag(1,-1,1,-1) has D^T J D = -J
    D = np.diag([1, 2, 1, 2])
    anti = [(S @ D) % 3 for S in G]
    anti_law = all((det2(S[:2, :2]) + det2(S[2:, :2])) % 3 == 2 for S in anti)
    anti_perfect = [S for S in anti if rank3(S[:2, :2]) == 2 and rank3(S[2:, :2]) == 2]
    tet = (np.kron(np.array([[1, 1], [1, 2]]), np.eye(2, dtype=int))) % 3
    tet_is_anti = np.array_equal((tet.T @ J @ tet) % 3, (2 * J) % 3)
    tet_perfect = any(np.array_equal(tet, S) for S in anti_perfect)
    orbs = local_orbits(perfect)
    tetD = (tet @ D) % 3
    tetD_symplectic = np.array_equal((tetD.T @ J @ tetD) % 3, J % 3)
    tetD_perfect = rank3(tetD[:2, :2]) == 2 and rank3(tetD[2:, :2]) == 2
    local_stabilizer = 576 * 576 // orbs[0] if len(orbs) == 1 else None
    also_ABBB = sum(1 for S in perfect if rank3(S[:2, 2:]) == 2 and rank3(S[2:, 2:]) == 2)
    res = dict(pass_id=11156, sp43_order=len(G), formula_checked_on_random_cliffords=all(checks), n_checks=len(checks),
               rank_census={f"space{a}_diag{b}": c for (a, b), c in sorted(ranks.items())},
               perfect_count=len(perfect), perfect_fraction=len(perfect) / len(G),
               perfect_all_blocks_invertible=also_ABBB, perfect_orders=dict(sorted(orders.items())),
               det_law_detAA_plus_detBA_eq_1=det_law, perfect_block_dets={str(k): v for k, v in perfect_dets.items()},
               antiunitary_det_law_eq_minus1=anti_law, antiunitary_perfect_count=len(anti_perfect),
               tetracode_gate_is_antiunitary=tet_is_anti, tetracode_gate_is_perfect=tet_perfect,
               perfect_local_double_cosets=orbs, tetracode_times_local_time_reversal_is_symplectic=tetD_symplectic,
               tetracode_times_local_time_reversal_is_perfect=tetD_perfect, local_stabilizer_order=local_stabilizer)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
