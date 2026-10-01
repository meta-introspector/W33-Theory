#!/usr/bin/env python3
"""Pass 11245: the counting lemma behind Pass 11244's Hesselink identity -- alternating forms are 2^-k of the
nondegenerate symmetric forms over F_2.

Pass 11244 found (on all 53 Jordan types with a part 2, m <= 6) that the fraction of unipotents of Sp(2m,2) with
chi_2 = 0 is 2^{-m_2} for m_2 even and 0 for m_2 odd.  chi_2 = 0 means the symmetric form b(v, w) = omega(v, (u-1)w)
induced on the m_2-dimensional multiplicity space of the size-2 blocks is alternating.  If that form is a uniformly
random nondegenerate symmetric form, the identity is exactly the statement
        #{nondegenerate alternating k x k over F_2} / #{nondegenerate symmetric k x k over F_2} = 2^{-k} (k even), 0 (odd).
This pass proves the counting statement (closed forms: |GL_k(2)|/|Sp_k(2)| and the MacWilliams count of symmetric
matrices) and checks it by brute force for k <= 5.  The uniformity of the induced form is NOT proved here -- that is
the remaining step for identity 2.  Identity 1 (Jordan types) is left as verified-only.
"""
from __future__ import annotations

import itertools
import json
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "w33_pass11245_symmetric_form_lemma.json"


def gl(k):
    out = 1
    for i in range(k):
        out *= 2 ** k - 2 ** i
    return out


def sp(k):
    m = k // 2
    out = 2 ** (m * m)
    for j in range(1, m + 1):
        out *= 4 ** j - 1
    return out


def nondeg_symmetric(k):
    """MacWilliams: number of nonsingular symmetric k x k matrices over F_q, q = 2"""
    q = 2
    out = q ** (k * (k + 1) // 2)
    for i in range(1, (k + 1) // 2 + 1):
        out = out * (q ** (2 * i - 1) - 1) // q ** (2 * i - 1)
    return out


def rank2(M):
    M = M.copy() % 2
    r = 0
    for c in range(M.shape[1]):
        nz = [i for i in range(r, M.shape[0]) if M[i, c]]
        if not nz:
            continue
        M[[r, nz[0]]] = M[[nz[0], r]]
        for i in range(M.shape[0]):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        r += 1
    return r


def brute(k):
    iu = [(i, j) for i in range(k) for j in range(i, k)]
    sym = alt = 0
    for bits in itertools.product(range(2), repeat=len(iu)):
        M = np.zeros((k, k), np.int64)
        for (i, j), b in zip(iu, bits):
            M[i, j] = M[j, i] = b
        if rank2(M) == k:
            sym += 1
            alt += not any(M[i, i] for i in range(k))
    return sym, alt


def run():
    rows = []
    for k in range(1, 13):
        s = nondeg_symmetric(k)
        a = gl(k) // sp(k) if k % 2 == 0 else 0
        rows.append(dict(k=k, nondeg_symmetric=s, nondeg_alternating=a, ratio=str(Fraction(a, s)),
                         equals_2_to_minus_k=Fraction(a, s) == (Fraction(1, 2 ** k) if k % 2 == 0 else 0)))
    br = []
    for k in range(1, 6):
        s, a = brute(k)
        br.append(dict(k=k, symmetric=s, alternating=a, matches_formula=(s == nondeg_symmetric(k) and
                                                                          a == (gl(k) // sp(k) if k % 2 == 0 else 0))))
    return dict(pass_id=11245, closed_forms=rows, all_k_le_12=all(r["equals_2_to_minus_k"] for r in rows), brute=br)


def main():
    res = run()
    OUT.write_text(json.dumps(res, indent=1))
    print(res["all_k_le_12"], [b["matches_formula"] for b in res["brute"]])


if __name__ == "__main__":
    main()
