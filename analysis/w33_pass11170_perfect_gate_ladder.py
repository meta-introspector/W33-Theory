#!/usr/bin/env python3
"""Pass 11170: perfect qutrit gates exist for 1, 2, 3 and 5 qutrits but NOT for 4 -- a perfect four-qutrit gate of any kind
would have an AME(8,3) Choi state, which the shadow inequalities exclude (Huber, Eltschka, Siewert, Guhne, arXiv:1708.06298; the
Huber-Wyderka table lists AME(8,3) as 'No (Shadow)').  We check the F9-linear slice independently and build an explicit
perfect five-qutrit Clifford gate.

F9-linear perfect gates on m qutrits = unitary m x m matrices over F9 = F3[i] with every entry and every 2x2 minor nonzero
(Pass 11164; larger minors follow by Jacobi's identity for unitary matrices).  Exhaustive clique search (rows = unit
vectors with all entries nonzero; edges = orthogonal with all 2x2 minors nonzero):
    m = 2: 32 row sets (= 64 / 2!),   m = 3: 2048 (= 12288 / 3!),   m = 4: 0.
The m = 4 zero is the F9 shadow of the AME(8,3) no-go (equivalently: no Hermitian self-dual [8,4,5]_9 code, since such a
code has generator [I | eta K] with K unitary and superregular).
Circulant search (K = sum_k c_k P^k, c in F9^m): m = 3: 24 (Pass 11164, e.g. P^2 + i*1); m = 4: 0; m = 5: 80.  The first
m = 5 example, with first row (i, i, 1-i, -1, 1-i), embeds as S in Sp(10,3) whose Choi state has entropy 5 (trits) on all
252 cuts 5|5 -- AME(10,3), a perfect five-qutrit gate.  (AME(10,3) is known: [[10,0,6]]_3 from the circulant Hermitian
self-dual [10,5,6]_9 code of Grassl-Gulliver, used by Grassl-Roetteler arXiv:1502.05267; Danielsen arXiv:1106.2428.)
Reading: a register of four qutrits cannot scramble perfectly in one tick; one, two, three and five can.
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11156_perfect_spacetime_gates as G2  # noqa: E402

OUT = ROOT / "data" / "w33_pass11170_perfect_gate_ladder.json"
F9 = np.array([(a, b) for a in range(3) for b in range(3)])


def mul(a, b):
    return np.stack([(a[..., 0] * b[..., 0] - a[..., 1] * b[..., 1]) % 3, (a[..., 0] * b[..., 1] + a[..., 1] * b[..., 0]) % 3], -1)


def conj(a):
    return np.stack([a[..., 0] % 3, (-a[..., 1]) % 3], -1)


def clique_count(m):
    nz = [tuple(x) for x in F9 if tuple(x) != (0, 0)]
    norm = {x: (x[0] ** 2 + x[1] ** 2) % 3 for x in nz}
    R = np.array([r for r in itertools.product(nz, repeat=m) if sum(norm[x] for x in r) % 3 == 1])
    a, b = R[..., 0], R[..., 1]
    re = (a[:, None, :] * a[None, :, :] + b[:, None, :] * b[None, :, :]).sum(-1) % 3
    im = (a[:, None, :] * b[None, :, :] - b[:, None, :] * a[None, :, :]).sum(-1) % 3
    ok = (re == 0) & (im == 0)
    for i, j in itertools.combinations(range(m), 2):
        p = mul(R[:, None, i], R[None, :, j])
        q = mul(R[:, None, j], R[None, :, i])
        ok &= ((p - q) % 3 != 0).any(-1)
    np.fill_diagonal(ok, False)
    adj = [set(np.flatnonzero(ok[i])) for i in range(len(R))]

    def extend(cl, cand):
        if len(cl) == m:
            return 1
        return sum(extend(cl + [v], cand & adj[v]) for v in cand if v > cl[-1])
    return len(R), sum(extend([i], adj[i]) for i in range(len(R)))


def circulant(c):
    m = len(c)
    return np.array([[c[(s - r) % m] for s in range(m)] for r in range(m)])


def is_perfect_unitary(K):
    m = K.shape[0]
    Gm = (mul(conj(K)[:, :, None, :], K[:, None, :, :]).sum(0)) % 3          # (K^dagger K)[a,b]
    if not all((Gm[a, b] == ([1, 0] if a == b else [0, 0])).all() for a in range(m) for b in range(m)):
        return False
    if (K == 0).all(-1).any():
        return False
    for (i, j) in itertools.combinations(range(m), 2):
        for (k, l) in itertools.combinations(range(m), 2):
            if ((mul(K[i, k], K[j, l]) - mul(K[i, l], K[j, k])) % 3 == 0).all():
                return False
    return True


def circulant_count(m):
    hits = []
    for idx in itertools.product(range(9), repeat=m):
        if 0 in idx:
            continue
        K = circulant(F9[list(idx)])
        if is_perfect_unitary(K):
            hits.append(idx)
    return hits


def embed(K):
    m = K.shape[0]
    S = np.zeros((2 * m, 2 * m), int)
    for i in range(m):
        for j in range(m):
            a, b = K[i, j]
            S[2 * i:2 * i + 2, 2 * j:2 * j + 2] = [[a, -b], [b, a]]
    return S % 3


def ame_cuts(S):
    n = S.shape[0] // 2
    Gm = np.vstack([np.eye(2 * n, dtype=int), S % 3])
    ents = set()
    for X in itertools.combinations(range(2 * n), n):
        comp = [p for p in range(2 * n) if p not in X]
        rows = [2 * p + t for p in comp for t in range(2)]
        ents.add(len(X) - 2 * n + G2.rank3(Gm[rows, :]))
    return sorted(ents), len(list(itertools.combinations(range(2 * n), n)))


def summarize():
    cl = {m: clique_count(m) for m in (2, 3, 4)}
    circ = {m: circulant_count(m) for m in (3, 4, 5)}
    K5 = circulant(F9[list(circ[5][0])])
    K3 = circulant(F9[list(circ[3][0])])
    J = lambda n: np.kron(np.eye(n, dtype=int), np.array([[0, 1], [-1, 0]]))
    S5, S3 = embed(K5), embed(K3)
    res = dict(pass_id=11170,
               f9_rows={str(m): v[0] for m, v in cl.items()}, f9_perfect_row_sets={str(m): v[1] for m, v in cl.items()},
               circulant_perfect={str(m): len(v) for m, v in circ.items()},
               five_qutrit_first_row=[[int(x) for x in F9[i]] for i in circ[5][0]],
               five_qutrit_symplectic=bool(np.array_equal((S5.T @ J(5) @ S5) % 3, J(5) % 3)),
               five_qutrit_cut_entropies=ame_cuts(S5)[0], five_qutrit_cuts=ame_cuts(S5)[1],
               three_qutrit_control_cut_entropies=ame_cuts(S3)[0],
               ame_8_3='No (Shadow): Huber et al. arXiv:1708.06298 (Huber-Wyderka table)',
               ladder={'1': True, '2': True, '3': True, '4': False, '5': True})
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
