"""Pass 11206: T-violation on the substrate -- every Clifford tick has a time-reversal symmetry, and for about half of
them it can only be anti-unitary.

A tick S is time-reversal symmetric if some time reversal T (anti-unitary Clifford = anti-symplectic map) inverts it:
T S T^-1 = S^-1 (projectively).  GAP (analysis/gap/w33_pass11206_time_reversal_invertible.g) tests every conjugacy class
of PSp(4,3) and PSp(6,3): g is T-invertible iff g^tau is PSp-conjugate to g^-1; it also records whether g is real
(inverted by a unitary, i.e. conjugate to g^-1 inside PSp).

Result: EVERY class is T-invertible (n = 2, 3) -- consistent with Wonenburger's theorem (every symplectic map is a
product of two skew-symplectic involutions, M. J. Wonenburger 1966), which this re-derives on the substrate.  The
physics is in the complement of the real classes: a tick that is not real can be run backwards only anti-unitarily.
Fractions of such ticks: 77/162 (two qutrits), 84013/157464 (three qutrits).  So the Clifford substrate has no
T-violation, and Wigner's anti-unitarity is forced for about half of its dynamics.

Independent check for two qutrits in Python: for each GAP class representative (Pass 11188 file) search the 25920
projective anti-symplectic maps for one that inverts it.
"""

from __future__ import annotations

import json
import re
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11156_perfect_spacetime_gates as G  # noqa: E402
import w33_pass11188_exact_arrow_by_class as E  # noqa: E402

GAPOUT = ROOT / "data" / "w33_pass11206_gap_time_reversal.txt"
OUT = ROOT / "data" / "w33_pass11206_t_violation.json"
TAU = np.diag([1, 2, 1, 2]).astype(np.int64)


def parse():
    text = re.sub(r"\s+", " ", GAPOUT.read_text().replace("\\\n", ""))
    rows = re.findall(r"CLASS n=(\d) order=(\d+) size=(\d+) T_invertible=(\w+) real=(\w+)", text)
    out = {2: [], 3: []}
    for n, o, s, t, r in rows:
        out[int(n)].append(dict(order=int(o), size=int(s), T_invertible=t == "true", real=r == "true"))
    return out


def python_check_two_qutrits():
    J = np.zeros((4, 4), np.int64)
    J[0, 1], J[1, 0], J[2, 3], J[3, 2] = 1, -1, 1, -1
    anti = []
    seen = set()
    for S in G.sp43():
        T = (TAU @ S) % 3
        k = min(tuple(T.ravel()), tuple(((-T) % 3).ravel()))
        if k not in seen:
            seen.add(k)
            anti.append(T)
    A = np.array(anti)
    Ainv = np.einsum('ij,nkj,kl->nil', J, A, J) % 3          # anti-symplectic: T^-1 = J T^t J
    ok_all, real_count = True, 0
    for c in E.parse_classes()[2]:
        S = c["S"]
        Sinv = (-J @ S.T @ J) % 3
        conj = np.einsum('nij,jk,nkl->nil', A, S, Ainv) % 3
        hit = ((conj == Sinv).all(axis=(1, 2)) | (conj == (-Sinv) % 3).all(axis=(1, 2))).any()
        ok_all &= bool(hit)
    return dict(every_class_inverted_by_an_anti_symplectic_map=bool(ok_all), anti_symplectic_projective=len(anti))


def run():
    cls = parse()
    order = {2: 25920, 3: 4585351680}
    res = dict(pass_id=11206)
    for n in (2, 3):
        rows = cls[n]
        assert sum(r["size"] for r in rows) == order[n]
        res[f"n{n}"] = dict(classes=len(rows), all_T_invertible=all(r["T_invertible"] for r in rows),
                           real_fraction=str(Fraction(sum(r["size"] for r in rows if r["real"]), order[n])),
                           only_antiunitary_fraction=str(Fraction(sum(r["size"] for r in rows if not r["real"]), order[n])),
                           non_real_orders=sorted({r["order"] for r in rows if not r["real"]}))
    res["python_check_n2"] = python_check_two_qutrits()
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
