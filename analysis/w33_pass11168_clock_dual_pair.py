#!/usr/bin/env python3
"""Pass 11168: the paper's clock, its Fourier-dual potential kick and the Lorentz group of the finite light cone never
make a perfect tick; a kick-tick-kick step with a Lorentz-breaking potential can -- and then it always time-reverses along
the light-cone reflection v0 <-> v2.

Objects (three qutrits = the coordinates (v0, v1, v2) of F_3^3, q(v) = v0 v2 - v1^2, Gram matrix of the polar form
M = [[0,0,1],[0,1,0],[1,0,0]] = M^{-1}):
  K      = the clock tick of Theorem 4.3 (Pass 11166): X_a -> X_a, Z_b -> Z_b X_{Mb}  (kinetic, e^{-i q(p)}),
  V(N)   = a potential kick e^{-i N(x)} for a symmetric N: X_a -> X_a Z_{Na}, Z_b -> Z_b;  V = V(M) is the Fourier dual of K,
  L      = the Lorentz transformations O(q) (48 linear maps preserving q), acting as X_a -> X_{La}, Z_b -> Z_{L^{-T} b}.
Results (party blocks S_ij; perfect = all nine invertible; pi as in Pass 11169):
  * <K, V> has order 24 (a copy of SL(2,3): K = I + M(x)E12, V = I + M(x)E21, M^2 = I); every element has blocks built from
    M, and M_01 = M_12 = 0 kills the four blocks between the transverse qutrit and the light-cone pair -- no perfect element.
  * With the Lorentz group, <K, V, O(q)> has order 576 (the Lorentz maps commute with K and V: L M L^T = M) and still NO
    perfect element.  The Lorentz-invariant potentials are exactly N = cM, so V(N) stays inside this group.
  * One clock tick is never enough: K V(N) has off-diagonal determinants -M_ij N_ij, zero at (0,1), for all 729 N.
    Tick-kick-tick K V(N) K is never perfect either: its (0,1) determinant is N_21 N_01 - N_21 N_01 = 0 identically.
  * Kick-tick-kick V(N) K V(N) is perfect for exactly 108 of the 729 potentials (all Lorentz-breaking), and EVERY one has
    pi = (v0 -> v2, v1 -> v1, v2 -> v0): the light-cone diagonal blocks have determinant (1 + N_02)^2, a square, so they can
    never reverse orientation, while the transverse block has determinant 1 + N_01 N_12, which perfectness forces to -1.
    So the perfect interacting tick time-reverses each light-cone coordinate into the OTHER one and the transverse one into
    itself -- pi is the light-cone reflection v0 <-> v2, the only nontrivial coordinate permutation preserving q.
  * Also: V(N1) K V(N2) is perfect for 23328 = 2^5 3^6 ordered pairs; K V(N) K V(N) for 30 potentials.
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11166_clock_tick_clifford as C  # noqa: E402
import w33_pass11169_perfect_gate_local_orbits as O  # noqa: E402

OUT = ROOT / "data" / "w33_pass11168_clock_dual_pair.json"
MQ = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]])


def xz_to_party(A, B, Cm, D):
    """build a per-qutrit (x0,z0,x1,z1,x2,z2) symplectic matrix from the (x|z) block form [[A, B], [Cm, D]]"""
    S = np.zeros((6, 6), int)
    for i in range(3):
        for j in range(3):
            S[2 * i, 2 * j], S[2 * i, 2 * j + 1] = A[i, j], B[i, j]
            S[2 * i + 1, 2 * j], S[2 * i + 1, 2 * j + 1] = Cm[i, j], D[i, j]
    return S % 3


I3, Z3 = np.eye(3, dtype=int), np.zeros((3, 3), int)


def kick(N):
    return xz_to_party(I3, Z3, N % 3, I3)               # x column gains z = N x


def lorentz_group():
    out = []
    for m in itertools.product(range(3), repeat=9):
        L = np.array(m).reshape(3, 3)
        if int(round(np.linalg.det(L))) % 3 == 0:
            continue
        if np.array_equal((L.T @ MQ @ L) % 3, MQ):
            out.append(L)
    return out


def inv3(L):
    d = int(round(np.linalg.det(L))) % 3
    adj = np.round(np.linalg.inv(L) * np.linalg.det(L)).astype(int)
    return (pow(d, -1, 3) * adj) % 3


def lorentz_gate(L):
    return xz_to_party(L, Z3, Z3, inv3(L).T)


def perfect(S):
    return all(int(S[2 * i, 2 * j] * S[2 * i + 1, 2 * j + 1] - S[2 * i, 2 * j + 1] * S[2 * i + 1, 2 * j]) % 3 != 0
               for i in range(3) for j in range(3))


def closure(gens):
    key = lambda S: S.tobytes()
    seen = {key(np.eye(6, dtype=np.int64)): np.eye(6, dtype=np.int64)}
    frontier = list(seen.values())
    gens = [g.astype(np.int64) for g in gens]
    while frontier:
        new = []
        for S in frontier:
            for g in gens:
                T = (g @ S) % 3
                k = key(T)
                if k not in seen:
                    seen[k] = T
                    new.append(T)
        frontier = new
    return list(seen.values())


def summarize():
    Kc = C.symplectic(C.clock())
    # check the clock is the kinetic form with Gram matrix MQ in one of the two conventions
    K_kin = xz_to_party(I3, MQ, Z3, I3)                   # z column gains x = M z
    V = kick(MQ)
    clock_is_kinetic = bool(np.array_equal(Kc, K_kin))
    G2 = closure([K_kin, V])
    Lg = lorentz_group()
    G3 = closure([K_kin, V] + [lorentz_gate(L) for L in Lg[:6]] + [lorentz_gate(L) for L in Lg])
    perfect_kv = sum(perfect(S) for S in G2)
    perfect_kvl = sum(perfect(S) for S in G3)
    # transverse decoupling inside <K,V>
    zero_blocks = all(not S[0:2, 2:4].any() and not S[2:4, 0:2].any() and not S[2:4, 4:6].any() and not S[4:6, 2:4].any() for S in G2)
    sym = lambda n: np.array([[n[0], n[1], n[2]], [n[1], n[3], n[4]], [n[2], n[4], n[5]]])
    allN = list(itertools.product(range(3), repeat=6))
    kicks = {n: kick(sym(n)) for n in allN}
    two_step = [n for n in allN if perfect((K_kin @ kicks[n]) % 3)]
    kvk = [n for n in allN if perfect((K_kin @ kicks[n] @ K_kin) % 3)]
    vkv = [n for n in allN if perfect((kicks[n] @ K_kin @ kicks[n]) % 3)]
    kvkv = [n for n in allN if perfect((K_kin @ kicks[n] @ K_kin @ kicks[n]) % 3)]
    v1kv2 = sum(perfect((kicks[a] @ K_kin @ kicks[b]) % 3) for a in allN for b in allN)
    pis = sorted({O.pattern((kicks[n] @ K_kin @ kicks[n]) % 3) for n in vkv})

    def bdet(S, i, j):
        return int(S[2 * i, 2 * j] * S[2 * i + 1, 2 * j + 1] - S[2 * i, 2 * j + 1] * S[2 * i + 1, 2 * j]) % 3
    formulas = all(bdet((kicks[n] @ K_kin @ kicks[n]) % 3, 0, 0) == (1 + n[2]) ** 2 % 3
                   and bdet((kicks[n] @ K_kin @ kicks[n]) % 3, 2, 2) == (1 + n[2]) ** 2 % 3
                   and bdet((kicks[n] @ K_kin @ kicks[n]) % 3, 1, 1) == (1 + n[1] * n[4]) % 3
                   and bdet((K_kin @ kicks[n] @ K_kin) % 3, 0, 1) == 0 for n in allN)
    lorentz_breaking = all(n not in [(0, 0, c, c, 0, 0) for c in range(3)] for n in vkv)
    lorentz_inv_potentials = [n for n in itertools.product(range(3), repeat=6)
                              if all(np.array_equal((L.T @ np.array([[n[0], n[1], n[2]], [n[1], n[3], n[4]], [n[2], n[4], n[5]]]) @ L) % 3,
                                                    np.array([[n[0], n[1], n[2]], [n[1], n[3], n[4]], [n[2], n[4], n[5]]])) for L in Lg)]
    res = dict(pass_id=11168, clock_is_kinetic_with_gram_M=clock_is_kinetic, order_K_V=len(G2), perfect_in_K_V=perfect_kv,
               transverse_blocks_vanish_in_K_V=zero_blocks, lorentz_order=len(Lg), order_K_V_Lorentz=len(G3),
               perfect_in_K_V_Lorentz=perfect_kvl, K_V_perfect=len(two_step), K_V_K_perfect=len(kvk),
               V_K_V_perfect=len(vkv), V_K_V_pis=[list(p) for p in pis], V_K_V_all_lorentz_breaking=lorentz_breaking,
               K_V_K_V_perfect=len(kvkv), V1_K_V2_perfect_pairs=v1kv2, determinant_formulas_hold=formulas,
               lorentz_invariant_potentials=[list(n) for n in lorentz_inv_potentials],
               V_K_V_example=list(vkv[0]) if vkv else None)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
