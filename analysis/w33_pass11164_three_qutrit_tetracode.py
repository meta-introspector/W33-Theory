#!/usr/bin/env python3
"""Pass 11164: the three-qutrit analogue of the tetracode gate lives over F9 -- K = P^2 + i*1 (a cyclic shift plus i times
the all-ones matrix) is unitary because 1^2 = 3*1 = 0 in characteristic 3, and it is perfect.

Why F9: over F3 no 3x3 matrix has all entries and all 2x2 minors nonzero (normalising each row to first entry 1, three
rows would need pairwise distinct ratios in {+1, -1}), so there is no scalar-block (F3-tetracode-like) perfect three-qutrit
gate -- equivalently, no [6,3,4] MDS code over F3.  Over F9 = F3[i] (i^2 = -1), with conjugation x -> x^3 and Hermitian
form <u, v> = sum conj(u_i) v_i, the imaginary part of <,> is the symplectic form of three qutrits; so a unitary K over F9
embeds (x = a + b i -> [[a, -b], [b, a]]) as a symplectic 6x6 matrix with F9-linear blocks, and it is perfect iff every
entry is nonzero (block det = norm).
Enumeration (analysis/w33_pass11164_scan_f9_unitary.py; frozen data/w33_pass11164_f9_unitary.json):
  * U(2, F9): order 96; 64 have all entries nonzero (F9-linear perfect two-qutrit gates; the tetracode gate times a
    primitive F9 scalar is among them);
  * U(3, F9): order 24192 (= |U(3,3)|); 12288 have all entries nonzero -- F9-linear perfect three-qutrit gates exist;
  * every column has norms (1, 1, 2) -- the (1, 1, -1) determinant pattern of Pass 11158;
  * 24 of them are circulant; one is
        K = P^2 + i * 1,     K^dagger K = (P - i 1)(P^2 + i 1) = I + 1^2 = I   (1^2 = 3 * 1 = 0 in characteristic 3),
    entries i and 1 + i, and all nine 2x2 minors nonzero, so its graph {(u, Ku)} is a [6,3,4]_9 MDS code; its Choi state has maximal entropy on every 3|3 cut (AME(6,3); its information pattern is
    checked by the stabilizer formula in Pass 11165).  The all-ones matrix is nilpotent exactly because the number of parties equals the characteristic.
Prior art (Pass 11170 correction): with eta = 1 + i (norm -1), [I | eta K] generates a Hermitian self-dual [6,3,4]_9
double-circulant code -- Grassl-Gulliver, Des. Codes Cryptogr. 52 (2009) 57-81, used for [[6,0,4]]_3 by Grassl-Roetteler
arXiv:1502.05267.  The closed form is an instance of that construction, not new.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "w33_pass11164_f9_unitary.json"
OUT = ROOT / "data" / "w33_pass11164_three_qutrit_tetracode.json"


def mul(x, y):
    return ((x[0] * y[0] - x[1] * y[1]) % 3, (x[0] * y[1] + x[1] * y[0]) % 3)


def conj(x):
    return (x[0] % 3, (-x[1]) % 3)


def summarize():
    d = json.loads(DATA.read_text())
    # K = P^2 + i*1: entry (r, c) = 1 + i if c - r = 2 mod 3 else i
    K = [[(1, 1) if (c - r) % 3 == 2 else (0, 1) for c in range(3)] for r in range(3)]
    KdK = [[(0, 0)] * 3 for _ in range(3)]
    for a in range(3):
        for b in range(3):
            s = (0, 0)
            for r in range(3):
                t = mul(conj(K[r][a]), K[r][b])
                s = ((s[0] + t[0]) % 3, (s[1] + t[1]) % 3)
            KdK[a][b] = s
    unitary = all(KdK[a][b] == ((1, 0) if a == b else (0, 0)) for a in range(3) for b in range(3))
    norms = sorted(((x[0] ** 2 + x[1] ** 2) % 3 for x in [K[r][0] for r in range(3)]))
    minors = []
    for a, b in ((0, 1), (0, 2), (1, 2)):
        for c, e in ((0, 1), (0, 2), (1, 2)):
            p, q = mul(K[a][c], K[b][e]), mul(K[a][e], K[b][c])
            minors.append(((p[0] - q[0]) % 3, (p[1] - q[1]) % 3))
    superregular = all(m != (0, 0) for m in minors)          # graph {(u, Ku)} is a [6,3,4]_9 MDS code
    res = dict(pass_id=11164, U2_order=d['2']['order'], U2_perfect=d['2']['all_entries_nonzero'],
               U3_order=d['3']['order'], U3_perfect=d['3']['all_entries_nonzero'], embedding_symplectic=d['3']['embedding_symplectic'],
               K_shift_plus_i_ones_unitary=unitary, K_column_norms=norms,
               K_superregular_graph_is_MDS_6_3_4_over_F9=superregular)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
