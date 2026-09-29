#!/usr/bin/env python3
"""Pass 11165: in a perfect three-qutrit tick each qutrit's past is secret-shared among ALL three futures -- forgetting any
single one of them erases it completely, and the erased record is still two trits.

Stabilizer formula (validated against brute-force entropies for two qutrits in Passes 11156, 11161): for an n-qutrit
Clifford gate S, the Choi state's Lagrangian is spanned by the columns of G = [I; S]; for a set X of the 2n legs,
S(X) = |X| - 2n + rank(G restricted to the legs outside X)   (units of log 3).
For each input leg A_in and each nonempty subset Y of the three output legs we compute I(A_in : Y) = S(A_in) + S(Y) -
S(A_in Y).  Gates: the three-qutrit tetracode analogue K = P^2 + i*1 over F9 (Pass 11164), perfect gates from the Pass 11158
sample, and non-perfect controls.
Results: see summarize().  For every perfect gate: I(A_in : Y) = 0 for every proper subset Y of the outputs and 2 for
all three -- a (3,3)-threshold sharing of every input's past (the quantum secret-sharing property of AME(6,3)).  So the
arrow of Theorem 4.7 is triggered by forgetting ANY one output, not both partners; and the erased record about A is two
trits, exactly as for two qutrits.  Non-perfect gates leak: some proper subset of outputs keeps information.
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

OUT = ROOT / "data" / "w33_pass11165_three_qutrit_arrow.json"
SAMPLE = ROOT / "data" / "w33_pass11158_sp63_sample.json"


def entropy(Gm, X, n):
    comp = [p for p in range(2 * n) if p not in X]
    rows = [2 * p + t for p in comp for t in range(2)]
    r = G2.rank3(Gm[rows, :]) if rows else 0
    return len(X) - 2 * n + r


def infos(S):
    n = S.shape[0] // 2
    Gm = np.vstack([np.eye(2 * n, dtype=int), S % 3])
    out = {}
    for a in range(n):
        for k in range(1, n + 1):
            for Y in itertools.combinations(range(n, 2 * n), k):
                I = entropy(Gm, (a,), n) + entropy(Gm, Y, n) - entropy(Gm, (a,) + Y, n)
                out[f"in{a}|out{tuple(y - n for y in Y)}"] = I
    return out


def f9_block(x):
    a, b = x
    return np.array([[a, -b], [b, a]]) % 3


def tetracode3():
    """K = P^2 + i*ones over F9 (entries: i off the shift, 1+i on it); embedded as a 6x6 symplectic matrix"""
    K = [[(1, 1) if (j - i) % 3 == 2 else (0, 1) for j in range(3)] for i in range(3)]
    S = np.zeros((6, 6), int)
    for i in range(3):
        for j in range(3):
            S[2 * i:2 * i + 2, 2 * j:2 * j + 2] = f9_block(K[i][j])
    return S % 3


def random_symplectic(rng, n=3, steps=80):
    D = 2 * n
    J = np.zeros((D, D), int)
    for k in range(n):
        J[2 * k, 2 * k + 1], J[2 * k + 1, 2 * k] = 1, -1
    M = np.eye(D, dtype=int)
    for _ in range(steps):
        v = rng.integers(0, 3, D); c = rng.integers(1, 3)
        M = (M + c * np.outer(v, (J @ v) @ M)) % 3
    return M, J


def summarize():
    J = np.zeros((6, 6), int)
    for k in range(3):
        J[2 * k, 2 * k + 1], J[2 * k + 1, 2 * k] = 1, -1
    Tk = tetracode3()
    tk_symplectic = bool(np.array_equal((Tk.T @ J @ Tk) % 3, J % 3))
    ex = np.array(json.loads(SAMPLE.read_text())['example']) % 3
    rng = np.random.default_rng(11165)
    perfect, controls = [Tk, ex], []
    while len(perfect) < 12 or len(controls) < 12:
        M, _ = random_symplectic(rng)
        dets = [[int(round(M[2*i, 2*j] * M[2*i+1, 2*j+1] - M[2*i, 2*j+1] * M[2*i+1, 2*j])) % 3 for j in range(3)] for i in range(3)]
        if all(d != 0 for row in dets for d in row):
            if len(perfect) < 12:
                perfect.append(M)
        elif len(controls) < 12:
            controls.append(M)

    def shared(S):
        """(max information about a single input in any PROPER subset of outputs, min and max over all three outputs)"""
        I = infos(S)
        full = [v for k, v in I.items() if 'out(0, 1, 2)' in k]
        return max(v for k, v in I.items() if 'out(0, 1, 2)' not in k), min(full), max(full)
    perf = [shared(S) for S in perfect]
    ctrl = [shared(S) for S in controls]
    res = dict(pass_id=11165, tetracode3_symplectic=tk_symplectic, tetracode3_infos=infos(Tk),
               perfect_max_info_in_proper_subsets=max(p[0] for p in perf), perfect_full_info=sorted({p[1] for p in perf} | {p[2] for p in perf}),
               controls_with_leak=sum(1 for c in ctrl if c[0] > 0), n_perfect=len(perf), n_controls=len(ctrl))
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    r = summarize()
    print(json.dumps({k: v for k, v in r.items() if k != 'tetracode3_infos'}, indent=1))
