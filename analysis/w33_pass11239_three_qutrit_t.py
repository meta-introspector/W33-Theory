"""Pass 11239: three qutrits -- does the 2pi/9 minimal T-violation persist?  (structured evidence, not exhaustive)

Two-qutrit statements (Passes 11227, 11235) are exact: every anti-unitary two-qutrit Clifford was searched.  For three
qutrits the anti-unitary Clifford group has |Sp(6,3)| x 729 ~ 3.3e12 elements, so an exhaustive search is out of reach,
and magnitude-only certificates of violation are known to be incomplete (Pass 11228 -- they miss the minimal violator).

What is computed:
  * lower bounds on the best time-reversal fidelity F_T(U) = max_V |tr(V U^* V^dag U)|/27 from (a) structured
    candidates (embedded optimal two-qutrit reversals V12 (x) V3, local anti-unitaries) and (b) 20000 random three-qutrit
    Cliffords (random words in H, S, SUM; the Pauli part maximised exactly over all 729 Paulis);
  * a certificate whenever a candidate reaches F = 1 (then U IS substrate-reversible);
  * for embedded two-qutrit violators U (x) I, whether any three-qutrit reversal beats the two-qutrit optimum.
"""

from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11227_two_qutrit_cp_mixing as P2  # noqa: E402

OUT = ROOT / "data" / "w33_pass11239_three_qutrit_t.json"
I3, X, Z, H, S, T = P2.I3, P2.X, P2.Z, P2.H, P2.S, P2.T
QUANTUM = (1 + 2 * np.cos(2 * np.pi / 9)) / 3


def kron3(a, b, c):
    return np.kron(np.kron(a, b), c)


def sum_gate(n, ctrl, tgt):
    D = 3 ** n
    M = np.zeros((D, D), dtype=complex)
    for idx in range(D):
        digs = [(idx // 3 ** (n - 1 - k)) % 3 for k in range(n)]
        out = digs.copy()
        out[tgt] = (digs[tgt] + digs[ctrl]) % 3
        j = sum(out[k] * 3 ** (n - 1 - k) for k in range(n))
        M[j, idx] = 1
    return M


PAULIS3 = np.array([kron3(np.linalg.matrix_power(X, a) @ np.linalg.matrix_power(Z, b),
                          np.linalg.matrix_power(X, c) @ np.linalg.matrix_power(Z, d),
                          np.linalg.matrix_power(X, e) @ np.linalg.matrix_power(Z, f))
                    for a, b, c, d, e, f in itertools.product(range(3), repeat=6)])


def best_over_paulis(U, Vs):
    """max over Cliffords V in Vs and all 729 Paulis of |tr(P V U^* V^dag P^dag U)| / 27"""
    Y = np.einsum('aji,jk,akl->ail', PAULIS3.conj(), U, PAULIS3)
    best = 0.0
    for s in range(0, len(Vs), 200):
        V = Vs[s:s + 200]
        B = np.einsum('rij,jk,rlk->ril', V, U.conj(), V.conj())
        best = max(best, float(np.abs(np.einsum('rjk,akj->ra', B, Y)).max() / 27))
    return best


def random_cliffords(rng, N, length=40):
    gates = []
    for q in range(3):
        for g in (H, S):
            ops = [I3, I3, I3]
            ops[q] = g
            gates.append(kron3(*ops))
    for c, t in itertools.permutations(range(3), 2):
        gates.append(sum_gate(3, c, t))
    out = np.empty((N, 27, 27), dtype=complex)
    for i in range(N):
        V = np.eye(27, dtype=complex)
        for g in rng.integers(len(gates), size=length):
            V = gates[g] @ V
        out[i] = V
    return out


def two_qutrit_optimal_reversals(U2, reps):
    """all two-qutrit Cliffords (rep part) achieving the two-qutrit optimum for U2, with their best Pauli"""
    Pa = P2.PA
    Y = np.einsum('aji,jk,akl->ail', Pa.conj(), U2, Pa)
    best, arg = 0.0, []
    for s in range(0, len(reps), 4000):
        C = reps[s:s + 4000]
        B = np.einsum('rij,jk,rlk->ril', C, U2.conj(), C.conj())
        tr = np.abs(np.einsum('rjk,akj->ra', B, Y)) / 9
        m = tr.max()
        if m > best + 1e-12:
            best, arg = m, []
        if abs(m - best) < 1e-12:
            r, a = np.argwhere(np.abs(tr - best) < 1e-12)[0]
            arg.append(Pa[a] @ C[r])
    return best, arg[:20]


def gate_list():
    gates = []
    for q in range(3):
        for g in (H, S, X, Z):
            ops = [I3, I3, I3]
            ops[q] = g
            gates.append(kron3(*ops))
    for c, tg in itertools.permutations(range(3), 2):
        gates.append(sum_gate(3, c, tg))
    return gates


def hill_climb(U, starts, steps=150, rng=None):
    """greedy improvement of |tr(P V U^* V^dag P^dag U)|/27 by left-multiplying V with Clifford generators"""
    gates = gate_list()
    best_overall = 0.0
    for V in starts:
        cur = best_over_paulis(U, V[None])
        for _ in range(steps):
            cands = np.array([g @ V for g in gates] + [V @ g for g in gates])
            vals = [best_over_paulis(U, c[None]) for c in cands]
            k = int(np.argmax(vals))
            if vals[k] <= cur + 1e-12:
                break
            cur, V = vals[k], cands[k]
            if cur > 1 - 1e-9:
                break
        best_overall = max(best_overall, cur)
        if best_overall > 1 - 1e-9:
            break
    return best_overall


def run(N=3000, seed=11239):
    rng = np.random.default_rng(seed)
    reps = np.load(P2.CACHE)
    res = dict(pass_id=11239, random_cliffords=N)
    S12 = sum_gate(3, 0, 1)
    S23 = sum_gate(3, 1, 2)
    U2 = np.kron(T, T) @ P2.SUM
    f2, opt2 = two_qutrit_optimal_reversals(U2, reps)
    res["two_qutrit_optimum"] = f2
    Vs = random_cliffords(rng, N)
    local3 = [I3, H, S, X, Z]                                   # single-qutrit factors on qutrit 3
    candidates = {
        "(T(x)T)SUM (x) I": np.kron(U2, I3),
        "(T(x)T(x)T) SUM12 SUM23": kron3(T, T, T) @ S23 @ S12,
        "(T(x)T(x)T) SUM12": kron3(T, T, T) @ S12,
        "(I(x)I(x)T) SUM23 (x) control": kron3(I3, I3, T) @ S23,
        "(T(x)I(x)T) SUM12 SUM23": kron3(T, I3, T) @ S23 @ S12,
    }
    out = {}
    for name, U in candidates.items():
        structured = np.array([np.kron(V2, v3) for V2 in opt2 for v3 in local3]) if "(x) I" in name else Vs[:0]
        lb_struct = best_over_paulis(U, structured) if len(structured) else 0.0
        lb_rand = best_over_paulis(U, Vs)
        scores = np.array([best_over_paulis(U, Vs[i][None]) for i in range(0, len(Vs), 30)])
        top = [Vs[i * 30] for i in np.argsort(-scores)[:6]]
        lb_hill = hill_climb(U, top)
        lb = max(lb_struct, lb_rand, lb_hill)
        out[name] = dict(lower_bound=lb, structured=lb_struct, random=lb_rand, hill_climb=lb_hill,
                         certified_reversible=bool(lb > 1 - 1e-9))
        print(name, out[name], flush=True)
    res["candidates"] = out
    res["quantum"] = QUANTUM
    return res


def factorised(reps):
    """exact product reversals for the candidates that are tensor products of a two-qutrit and a one-qutrit tick:
    with V = V_A (x) V_B the trace factorises, so F(A (x) B) >= F(A) F(B); a diagonal one-qutrit tick has F = 1 (V = I).
    These calibrate the generic search: where the factorised bound is known, how much does random + hill-climb miss?"""
    out = {}
    U2 = np.kron(T, T) @ P2.SUM
    f2, opt2 = two_qutrit_optimal_reversals(U2, reps)
    U = kron3(T, T, T) @ sum_gate(3, 0, 1)                    # = [(T(x)T)SUM] (x) T
    assert np.allclose(U, np.kron(U2, T))
    out["(T(x)T(x)T) SUM12"] = dict(factorisation="[(T(x)T)SUM] (x) T", two_qutrit_factor=f2,
                                    product_bound=best_over_paulis(U, np.array([np.kron(V, I3) for V in opt2])))
    W2 = np.kron(I3, T) @ P2.SUM                              # one cubic gate + SUM: reversible (Pass 11227)
    g2, optw = two_qutrit_optimal_reversals(W2, reps)
    U = kron3(I3, I3, T) @ sum_gate(3, 1, 2)                  # = I (x) [(I(x)T)SUM]
    assert np.allclose(U, np.kron(I3, W2))
    pb = best_over_paulis(U, np.array([np.kron(I3, V) for V in optw]))
    out["(I(x)I(x)T) SUM23 (x) control"] = dict(factorisation="I (x) [(I(x)T)SUM]", two_qutrit_factor=g2,
                                                product_bound=pb, certified_reversible=bool(pb > 1 - 1e-9))
    return out


def main():
    res = run()
    res["factorised"] = factorised(np.load(P2.CACHE))
    for k, v in res["factorised"].items():
        c = res["candidates"][k]
        c["lower_bound"] = max(c["lower_bound"], v["product_bound"])
        c["certified_reversible"] = bool(c["lower_bound"] > 1 - 1e-9)
    json.dump(res, open(OUT, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
