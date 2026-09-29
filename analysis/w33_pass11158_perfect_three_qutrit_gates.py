#!/usr/bin/env python3
"""Pass 11158: perfect three-qutrit Clifford gates exist (Choi state = AME(6,3)), and the determinant law generalises --
every column and every row of blocks sums to 1.

GENERAL LAW (proof): for S in Sp(2n, 3) acting on n qutrits with 2x2 blocks S_ij (output party i, input party j), the
symplectic condition restricted to input party j reads sum_i w(S_ij u, S_ij u') = w(u, u'), and w(M u, M u') = det(M) w(u, u')
on F_3^2, so
        sum_i det S_ij = 1  (mod 3)  for every input party j,   and (S^T is symplectic) sum_j det S_ij = 1 for every i.
For n = 2 this is Pass 11156's det S_AA + det S_BA = 1.
Scan: analysis/w33_pass11158_scan_sp63.py -- 200000 uniform samples of Sp(6,3) (random products of 80 symplectic
transvections, vectorised; |Sp(6,3)| = 9 170 703 360 is too large to enumerate).  Frozen: data/w33_pass11158_sp63_sample.json.
Results:
  * the column and row laws hold in 100% of samples (as proved);
  * a three-qutrit Clifford gate is perfect (its Choi state maximal across all 20 cuts 3|3 -- AME(6,3)) iff every 1x1 and
    every 2x2 block minor is invertible; empirically the 2x2 minors are then automatic (Jacobi's complementary-minor
    identity for symplectic S): perfect <=> all nine blocks invertible;
  * perfect three-qutrit gates EXIST: fraction 0.03108 +- 0.0004 (6216 of 200000), i.e. about 2.85e8 gates;
  * every perfect gate has, in each column, determinants (1, 1, -1) in some order -- the only way three nonzero values
    in {1, 2} can sum to 1 mod 3.  The two-qutrit rule "both -1" becomes "each column has exactly one -1".
Scope: exact law; the fraction is a sampling estimate (exact count open).  The existence of AME(6,3) is known in the
literature; new here: that W(3,3)'s three-qutrit Clifford group realises it, and the determinant law.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "w33_pass11158_sp63_sample.json"
OUT = ROOT / "data" / "w33_pass11158_perfect_three_qutrit_gates.json"
SP63 = 3 ** 9 * (3 ** 2 - 1) * (3 ** 4 - 1) * (3 ** 6 - 1)


def summarize():
    d = json.loads(DATA.read_text())
    f = d['perfect_fraction']
    err = math.sqrt(f * (1 - f) / d['samples'])
    res = dict(pass_id=11158, sp63_order=SP63, samples=d['samples'], column_law=d['column_law'], row_law=d['row_law'],
               perfect_fraction=f, perfect_fraction_stderr=err, perfect_estimate=f * SP63,
               perfect_equals_single_blocks=abs(d['perfect_fraction'] - d['single_blocks_invertible']) < 1e-12,
               column_patterns=d['perfect_column_det_patterns'], example=d['example'])
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    r = summarize()
    print(json.dumps({k: v for k, v in r.items() if k != 'example'}, indent=1))
