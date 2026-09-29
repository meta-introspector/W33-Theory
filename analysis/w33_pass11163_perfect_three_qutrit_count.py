#!/usr/bin/env python3
"""Pass 11163: the exact number of perfect three-qutrit Clifford gates is 286 654 464 = 2^17 3^7, a fraction 128/4095 of
Sp(6,3) (two qutrits: 13824 = 2^9 3^3, fraction 4/15).

Method (analysis/w33_pass11163_scan_perfect_count.py): a symplectic S is an ordered symplectic basis
(s1, s2 | s3, s4 | s5, s6).  The determinant of party block i of a column pair (a, b) is the LOCAL symplectic form
w_i(a, b) (so the column law sum_i det = w(a, b) = 1 is automatic).  S is perfect iff all nine party blocks are
invertible.  (Proof of the 'empirical' step of Pass 11158: a 3|3 cut with k inputs needs the 2k x 2k submatrix of S on
k output and k input parties invertible; for k = 2, Jacobi's identity with det S = 1 and S^{-1} = -J S^T J gives
det S[{i,i'},{j,j'}] = +-det S_{i''j''}, the complementary block -- so the nine blocks suffice.)
Enumerate the 176904 symplectic pairs (a, b); 41472 have all three local forms nonzero.  For each, enumerate the pairs
(c, d) in the 4-dimensional complement with all local forms nonzero AND 1 - w_i(a,b) - w_i(c,d) != 0 (the third column's
block determinants, fixed by the row law); the third pair is then any symplectic basis of the remaining plane
(|SL(2,3)| = 24 choices, block determinants unchanged).  Total good (first, second) pairs: 11 943 936, so
        #perfect = 24 x 11 943 936 = 286 654 464 = 2^17 3^7,   fraction 128/4095 of |Sp(6,3)| = 9 170 703 360,
consistent with the Pass 11158 sample estimate 0.03108 +- 0.0004 (exact 0.0312576).
Known-positive control (recomputed here): the same method for two qutrits -- first pairs with local forms (-1, -1), times
24 for the second pair -- returns 13824 = the Pass 11156 enumeration.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "w33_pass11163_perfect_count.json"
OUT = ROOT / "data" / "w33_pass11163_perfect_three_qutrit_count.json"


def two_qutrit_control():
    import itertools
    import numpy as np
    V = np.array(list(itertools.product(range(3), repeat=4)))
    w = [(V[:, None, 2 * i] * V[None, :, 2 * i + 1] - V[:, None, 2 * i + 1] * V[None, :, 2 * i]) % 3 for i in range(2)]
    good = (w[0] == 2) & (w[1] == 2)          # both blocks invertible and summing to 1: (-1, -1); second pair: dets 1-(-1) = -1
    return 24 * int(good.sum())


def summarize():
    d = json.loads(DATA.read_text())
    frac = Fraction(d['perfect_count'], d['sp63'])
    res = dict(pass_id=11163, **d, fraction_exact=str(frac), count_is_24_times=d['perfect_count'] == 24 * d['good_first_two'],
               within_sample_estimate=abs(float(frac) - 0.03108) < 3 * 0.00039,
               two_qutrit_control=two_qutrit_control())
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
