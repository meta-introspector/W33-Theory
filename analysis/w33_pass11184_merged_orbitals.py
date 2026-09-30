#!/usr/bin/env python3
"""Pass 11184: separating the three-qutrit relations that the Choi signature merges (Pass 11178).

Independent recomputation (no GAP): the stabiliser of the standard split F0 -- local SL(2,3) on each qutrit and the
relabellings of the three qutrits -- acts on the 110565 splits; union-find over its generators gives the orbits, which
must be the 20 suborbits of Pass 11178 (1, 36, 96, 256, 256, 768, 864, 2304 x3, 3456, 4608, 6912 x6, 10368, 41472).
For a split F with adapted basis B_F, the block B_F[P_i, Q_j] is the projection of F's plane Q_j onto the standard plane
P_i.  Beyond its rank/orientation code (the signature), a rank-one block has an IMAGE point in P_i and a KERNEL point in
Q_j.  Invariants of the pair (F0, F) under the stabiliser: for each row, the number of distinct image points among its
rank-one blocks; for each column, the number of distinct kernel points; and for each rank-one block, whether its kernel
point is mapped by the other blocks of its column onto image points of rank-one blocks (recorded as incidence counts).
RESULTS.
  * The 20 orbits are recomputed independently (union-find over the stabiliser's generators) with exactly GAP's sizes.
  * The fine invariant is constant on orbits and separates more: the class {6912, 6912} splits, and the five-orbital
    class splits into {256, 256}, {2304}, {6912, 6912}.
  * 14 relations are symmetric (the reversed relation (F, F0) lies in the same orbital); the other 6 form three pairs of
    mutually REVERSED relations -- a gate and its inverse -- of sizes 256, 2304, 6912.  The two ties left by the fine
    invariant are exactly two of these reversal pairs (the third, 2304, is already split by the signature).
So the Choi-entanglement data (block ranks, orientations, image/kernel incidences) classify the relation between two
three-qutrit subsystem splits completely, UP TO THE DIRECTION OF TIME: for two pairs of relations, a gate and its inverse
look identical to every one of these measures, and only the geometry of the relation knows which way it runs.
"""
from __future__ import annotations

import itertools
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11182_paper_ticks_mereology as P  # noqa: E402

OUT = ROOT / "data" / "w33_pass11184_merged_orbitals.json"
SUBDEGREES = sorted([1, 36, 96, 256, 256, 768, 864, 2304, 2304, 2304, 3456, 4608, 6912, 6912, 6912, 6912, 6912, 6912, 10368, 41472])


def norm(v):
    v = np.asarray(v) % 3
    i = next(k for k in range(len(v)) if v[k])
    return tuple((v * pow(int(v[i]), -1, 3)) % 3)


def plane_key(u, v):
    return frozenset(norm((a * u + b * v) % 3) for a in range(3) for b in range(3) if ((a * u + b * v) % 3).any())


def fac_key(B):
    return frozenset(plane_key(B[:, 2 * k], B[:, 2 * k + 1]) for k in range(3))


def generators():
    gens = []
    for q in range(3):
        for m in ([[1, 1], [0, 1]], [[1, 0], [1, 1]]):
            g = np.eye(6, dtype=np.int64)
            g[2 * q:2 * q + 2, 2 * q:2 * q + 2] = m
            gens.append(g)
    for perm in ([1, 0, 2], [1, 2, 0]):
        g = np.zeros((6, 6), np.int64)
        for new, old in enumerate(perm):
            g[2 * new:2 * new + 2, 2 * old:2 * old + 2] = np.eye(2, dtype=np.int64)
        gens.append(g)
    return gens


def orbits(Bs, return_index=False):
    keys = [fac_key(B) for B in Bs]
    index = {k: i for i, k in enumerate(keys)}
    parent = list(range(len(Bs)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for g in generators():
        GB = np.einsum('ij,fjk->fik', g, Bs) % 3
        for i in range(len(Bs)):
            j = index[fac_key(GB[i])]
            a, b = find(i), find(j)
            if a != b:
                parent[a] = b
    groups = defaultdict(list)
    for i in range(len(Bs)):
        groups[find(i)].append(i)
    if return_index:
        return list(groups.values()), index
    return list(groups.values())


def reverse_orbit(B, index, orb_of):
    """orbit of the reversed relation: g = B^-1 maps F to F0, so (F, F0) ~ (F0, g F0); g F0 has adapted basis B^-1"""
    J = np.zeros((6, 6), np.int64)
    for k in range(3):
        J[2 * k, 2 * k + 1], J[2 * k + 1, 2 * k] = 1, -1
    Binv = (-J @ B.T @ J) % 3
    return orb_of[index[fac_key(Binv)]]


def invariants(B):
    blocks = B.reshape(3, 2, 3, 2).transpose(0, 2, 1, 3)
    det = (blocks[..., 0, 0] * blocks[..., 1, 1] - blocks[..., 0, 1] * blocks[..., 1, 0]) % 3
    zero = (blocks == 0).all(axis=(-1, -2))
    code = np.where(zero, 0, np.where(det == 0, 1, np.where(det == 1, 2, 3)))
    img, ker = {}, {}
    for i in range(3):
        for j in range(3):
            if code[i, j] == 1:
                Mb = blocks[i, j]
                col = Mb[:, 0] if Mb[:, 0].any() else Mb[:, 1]
                img[(i, j)] = norm(col)
                ker[(i, j)] = next(norm(c) for c in ([1, 0], [0, 1], [1, 1], [1, 2]) if not ((Mb @ np.array(c)) % 3).any())
    row_distinct = sorted(len({img[(i, j)] for j in range(3) if (i, j) in img}) for i in range(3))
    col_distinct = sorted(len({ker[(i, j)] for i in range(3) if (i, j) in ker}) for j in range(3))
    # incidence: for a rank-one block (i,j), apply the other blocks (i',j), i' != i, to its kernel vector, and ask whether
    # the result is an image point of a rank-one block in row i'
    inc = 0
    for (i, j), kv in ker.items():
        for i2 in range(3):
            if i2 != i and (i2, j) in img:
                w = (blocks[i2, j] @ np.array(kv)) % 3
                if w.any() and norm(w) == img[(i2, j)]:
                    inc += 1
    canon = min(tuple(code[list(r)][:, list(c)].flatten()) for r in itertools.permutations(range(3))
                for c in itertools.permutations(range(3)))
    return canon, (tuple(row_distinct), tuple(col_distinct), inc)


def summarize():
    Bs = P.load_facs()
    orbs, index = orbits(Bs, return_index=True)
    orb_of = {}
    for k, o in enumerate(orbs):
        for i in o:
            orb_of[i] = k
    paired = {k: reverse_orbit(Bs[o[0]], index, orb_of) for k, o in enumerate(orbs)}
    self_paired = sum(1 for k, v in paired.items() if v == k)
    sizes = sorted(len(o) for o in orbs)
    table = []
    for k, o in enumerate(orbs):
        invs = Counter(invariants(Bs[i]) for i in o)
        table.append(dict(orbit=k, size=len(o), reverse=paired[k],
                          invariants=[[list(kk[0]), [list(kk[1][0]), list(kk[1][1]), kk[1][2]], v] for kk, v in invs.items()]))
    by_sig = defaultdict(list)
    for row in table:
        sigs = {tuple(int(x) for x in s[0]) for s in row['invariants']}
        fines = {json.dumps(x[1]) for x in row['invariants']}
        by_sig[tuple(sorted(sigs))].append((row['size'], sorted(fines), row['orbit'], row['reverse']))
    merged = {str(k): [(a, b) for a, b, _, _ in v] for k, v in by_sig.items() if len(v) > 1}
    separated = all(len({json.dumps(f) for _, f, _, _ in v}) == len(v) for v in by_sig.values())
    # remaining ties: are they exactly the pairs {O, O^T} of mutually reversed (non-self-paired) orbitals?
    ties = []
    for v in by_sig.values():
        groups = defaultdict(list)
        for size, f, o, rev in v:
            groups[json.dumps(f)].append((o, rev))
        for g in groups.values():
            if len(g) > 1:
                ties.append(g)
    ties_are_reversal_pairs = all(len(g) == 2 and g[0][1] == g[1][0] and g[1][1] == g[0][0] and g[0][0] != g[1][0]
                                  for g in ties)
    constant = all(len({json.dumps(x[1]) for x in row['invariants']}) == 1 for row in table)
    res = dict(pass_id=11184, orbits=len(orbs), orbit_sizes=sizes, matches_gap_subdegrees=sizes == SUBDEGREES,
               fine_invariant_constant_on_orbits=constant, fine_invariant_separates_merged=separated,
               merged_classes=merged, self_paired_orbitals=self_paired, remaining_ties=len(ties),
               remaining_ties_are_reversal_pairs=ties_are_reversal_pairs,
               non_self_paired_sizes=sorted(len(orbs[k]) for k, v in paired.items() if v != k))
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
