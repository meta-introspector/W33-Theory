"""Pass 11202: where an E6 reflection is local -- a double-six and the generalized quadrangle GQ(2,2).

Pass 11192 found that among the time reversals of the two-qutrit substrate (the anti-symplectic half of W(E6)) exactly
the 36 reflections have the maximal arrow A = 4, although each is block-monomial (it swaps the two factors) in 15 of the
45 splits.  In the cubic-surface dictionary (Pass 11177 / memory: 45 splits = tritangent planes, 27 frames = lines, a
frame = 5 splits whose octets partition the 40 points; line in plane <=> split in frame) this pass identifies them:

  * each reflection fixes exactly 15 frames and moves 12, and the 12 moved frames form a DOUBLE-SIX (two sixes of
    pairwise skew lines, each line meeting all lines of the other six but its partner), the reflection swapping
    partners -- the classical reflection <-> double-six correspondence of W(E6);
  * the 15 fixed splits are exactly the tritangent planes made of three fixed lines, and fixed lines with fixed planes
    form the generalized quadrangle GQ(2,2) = W(2) (15 points, 15 lines, 3 per line and per point, GQ axiom checked);
  * the map reflection -> double-six is a bijection onto the 36 double-sixes.
"""

from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11156_perfect_spacetime_gates as G  # noqa: E402
import w33_pass11180_mereology as M  # noqa: E402

OUT = ROOT / "data" / "w33_pass11202_reflection_double_six.json"
TAU = np.diag([1, 2, 1, 2]).astype(np.int64)


def norm(v):
    i = next(k for k in range(len(v)) if v[k] % 3)
    return tuple(int(x) for x in (v * pow(int(v[i]), -1, 3)) % 3)


def plane_pts(u, v):
    return frozenset(norm((a * u + b * v) % 3) for a in range(3) for b in range(3) if ((a * u + b * v) % 3).any())


def splits():
    out = []
    for B in M.factorisations(2):
        out.append(frozenset([plane_pts(B[:, 0], B[:, 1]), plane_pts(B[:, 2], B[:, 3])]))
    return out


def frames(S):
    octet = [p1 | p2 for p1, p2 in (tuple(s) for s in S)]
    n = len(S)
    disjoint = [[not (octet[i] & octet[j]) for j in range(n)] for i in range(n)]
    out = []
    for c in itertools.combinations(range(n), 5):
        if all(disjoint[a][b] for a, b in itertools.combinations(c, 2)):
            out.append(frozenset(c))
    return out


def act_split(T, s):
    return frozenset(frozenset(norm((T @ np.array(p)) % 3) for p in plane) for plane in s)


def is_gq(points, lines):
    if any(len(l) != 3 for l in lines):
        return False
    if any(sum(p in l for l in lines) != 3 for p in points):
        return False
    for l in lines:
        for p in points:
            if p in l:
                continue
            meets = [q for q in l if any(p in m and q in m for m in lines)]
            if len(meets) != 1:
                return False
    return True


def run():
    S = splits()
    idx = {s: i for i, s in enumerate(S)}
    Fr = frames(S)
    assert len(S) == 45 and len(Fr) == 27
    fidx = {f: i for i, f in enumerate(Fr)}
    meet = lambda a, b: a != b and bool(Fr[a] & Fr[b])          # lines meet iff they share a tritangent plane
    J = M.form(2)
    # reflections: anti-symplectic involutions (projective) moving exactly 12 frames -- identified below as class 36
    refl, seen = [], set()
    for Ssp in G.sp43():
        T = (TAU @ Ssp) % 3
        key = min(tuple(T.ravel()), tuple(((-T) % 3).ravel()))
        if key in seen:
            continue
        seen.add(key)
        T2 = (T @ T) % 3
        if not (np.array_equal(T2, np.eye(4, dtype=np.int64)) or np.array_equal(T2, (2 * np.eye(4, dtype=np.int64)) % 3)):
            continue
        sp = [idx[act_split(T, s)] for s in S]
        fp = [fidx[frozenset(sp[i] for i in f)] for f in Fr]
        refl.append((T, sp, fp))
    by_moved = {}
    for T, sp, fp in refl:
        by_moved.setdefault(sum(fp[i] != i for i in range(27)), []).append((T, sp, fp))
    res = dict(pass_id=11202, outer_involutions={int(k): len(v) for k, v in sorted(by_moved.items())})
    R = by_moved[12]
    assert len(R) == 36
    ok_ds = ok_planes = ok_gq = ok_local = True
    doublesixes = set()
    for T, sp, fp in R:
        moved = [i for i in range(27) if fp[i] != i]
        fixed = [i for i in range(27) if fp[i] == i]
        # double-six: moved lines split into two sixes of pairwise skew lines, partner = image under T
        colour = {moved[0]: 0}
        stack = [moved[0]]
        while stack:                                   # 2-colour the meeting graph on the moved lines
            a = stack.pop()
            for b in moved:
                if meet(a, b) and b not in colour:
                    colour[b] = 1 - colour[a]
                    stack.append(b)
        six1 = [m for m in moved if colour.get(m) == 0]
        six2 = [m for m in moved if colour.get(m) == 1]
        ok = (len(six1) == 6 and all(not meet(a, b) for a, b in itertools.combinations(six1, 2))
              and all(not meet(a, b) for a, b in itertools.combinations(six2, 2))
              and all(fp[a] in six2 for a in six1)
              and all(meet(a, b) == (fp[a] != b) for a in six1 for b in six2))
        ok_ds &= ok
        doublesixes.add(frozenset([frozenset(six1), frozenset(six2)]))
        fixed_splits = [i for i in range(45) if sp[i] == i]
        planes_of_fixed = [i for i in range(45) if all(i in Fr[l] for l in [l for l in fixed if i in Fr[l]])
                           and sum(i in Fr[l] for l in fixed) == 3]
        ok_planes &= (len(fixed_splits) == 15 and sorted(fixed_splits) == sorted(planes_of_fixed))
        pts = fixed
        lines = [frozenset(l for l in fixed if s in Fr[l]) for s in fixed_splits]
        ok_gq &= is_gq(pts, lines)
        # block-monomial (swap) in the fixed splits
        loc = 0
        for s in fixed_splits:
            B = M.factorisations(2)[s]
            Binv = (-J @ B.T @ J) % 3
            TF = (Binv @ T @ B) % 3
            loc += int(not TF[:2, :2].any() and not TF[2:, 2:].any())
        ok_local &= loc == 15
    res.update(reflections=len(R), double_six_each=bool(ok_ds), distinct_double_sixes=len(doublesixes),
               fixed_splits_are_planes_of_fixed_lines=bool(ok_planes), fixed_geometry_is_GQ22=bool(ok_gq),
               swaps_factors_in_all_15=bool(ok_local))
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1)
    print(res)


if __name__ == "__main__":
    main()
