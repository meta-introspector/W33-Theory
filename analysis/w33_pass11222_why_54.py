#!/usr/bin/env python3
"""Pass 11222: why chi = 54 Q -- the Maslov chirality as a signed count of orbits of the relation's symmetry group.

Pass 11212 found chi = 54 Q on every member of the 256 pair (chi = Pass 11209's profile-oriented Maslov count, Q =
master Pass 11189's twist).  This pass explains the factor 27 and isolates what remains.

For a relation (F0, F) its symmetry group Gamma = Stab(F0) n Stab(F) (local Cliffords and qutrit relabellings fixing
both splits) acts on oriented triples (L1, L2, M) of product Lagrangians (L1, L2 of F0 with profile(L1) > profile(L2),
M of F), preserving the Kashiwara-Maslov class.  So chi = sum over Gamma-orbits of odd-class triples of
(+1 for class 1, -1 for class 3) x orbit size.

Results:
  * 256 pair: |Gamma| = 324 = 2^2 3^4.  The odd triples fall into 409 orbits of sizes 162, 81, 54, 27 -- every size a
    multiple of 27 (no odd-class triple is fixed by an element of order 3 beyond a group of order 3), hence 27 | chi.
    Classes 1 and 3 balance orbit-size by orbit-size except for three orbits, and chi = -162 + 81 + 27 = -54 on the
    Q = -1 orbital (mirrored on Q = +1).  So 54 = 2 * 27: the 27 is the 3-part of the odd triples' orbits, the 2 is
    the net (-6 + 3 + 1) of the unbalanced orbits in units of 27.
  * Only four (profile, intersection) cells carry net chirality, all with L2 and M of the 'unique partner' profile
    (45, 18, 1, 0); their nets are +-162, -+54, -+243, +-81.
  * 6912 pair (control): |Gamma| = 12, so no 27-divisibility is forced; chi = +-18 there.
A conceptual reason for the net 2 (equivalently: why the unbalanced orbits are one each of sizes 162, 81 and 27) is
not found here.
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
import w33_pass11184_merged_orbitals as O  # noqa: E402
import w33_pass11189_oriented_holonomy as H  # noqa: E402
import w33_pass11209_maslov_chirality as M  # noqa: E402

REPS = ROOT / "data" / "w33_pass11212_one_chirality.json"
OUT = ROOT / "data" / "w33_pass11222_why_54.json"
SL = [np.array([[a, b], [c, d]]) for a in range(3) for b in range(3) for c in range(3) for d in range(3)
      if (a * d - b * c) % 3 == 1]


def relation_stabilizer(B):
    key = O.fac_key(B)
    out = []
    for perm in itertools.permutations(range(3)):
        P = np.zeros((6, 6), np.int64)
        for new, old in enumerate(perm):
            P[2 * new:2 * new + 2, 2 * old:2 * old + 2] = np.eye(2, dtype=np.int64)
        for a in SL:
            for b in SL:
                for c in SL:
                    g = np.zeros((6, 6), np.int64)
                    g[0:2, 0:2], g[2:4, 2:4], g[4:6, 4:6] = a, b, c
                    g = (g @ P) % 3
                    if O.fac_key((g @ B) % 3) == key:
                        out.append(g)
    return out


def lkey(L):
    return frozenset(tuple((L @ np.array(c)) % 3) for c in itertools.product(range(3), repeat=3))


def analyse(B):
    stab = relation_stabilizer(B)
    P0 = M.product_lagrangians(np.eye(6, dtype=np.int64))
    PF = M.product_lagrangians(B)
    k0 = {lkey(L): i for i, L in enumerate(P0)}
    kF = {lkey(L): i for i, L in enumerate(PF)}
    act0 = [[k0[lkey((g @ L) % 3)] for L in P0] for g in stab]
    actF = [[kF[lkey((g @ L) % 3)] for L in PF] for g in stab]
    i, m = (a.ravel() for a in np.meshgrid(np.arange(64), np.arange(64), indexing="ij"))
    dim = (6 - M.rank_batch(np.concatenate([P0[i], PF[m]], 2))).reshape(64, 64)
    prof0 = [tuple(sorted(r)) for r in dim.tolist()]
    profF = [tuple(sorted(c)) for c in dim.T.tolist()]
    pairs = [(a, b) for a in range(64) for b in range(64) if prof0[a] > prof0[b]]
    i1 = np.repeat([a for a, _ in pairs], 64)
    i2 = np.repeat([b for _, b in pairs], 64)
    mm = np.tile(np.arange(64), len(pairs))
    cls = M.kashiwara(P0[i1], P0[i2], PF[mm])
    cl = {(a, b, c): int(k) for a, b, c, k in zip(i1.tolist(), i2.tolist(), mm.tolist(), cls.tolist())}
    seen, orbits = set(), Counter()
    for t, c in cl.items():
        if t in seen or c % 2 == 0:
            continue
        orb = {(act0[g][t[0]], act0[g][t[1]], actF[g][t[2]]) for g in range(len(stab))}
        seen |= orb
        assert len({cl[x] for x in orb}) == 1, "Maslov class must be constant on orbits"
        orbits[(len(orb), c)] += 1
    chi = sum(n * size * (1 if c == 1 else -1) for (size, c), n in orbits.items())

    def comp(p):
        cc = Counter(p)
        return tuple(cc.get(k, 0) for k in range(4))
    net = defaultdict(int)
    for (a, b, c), k in cl.items():
        if k % 2:
            net[(comp(prof0[a]), comp(prof0[b]), comp(profF[c]), int(dim[a, c]), int(dim[b, c]))] += 1 if k == 1 else -1
    unbalanced = {}
    for size in sorted({s for s, _ in orbits}):
        d = orbits.get((size, 1), 0) - orbits.get((size, 3), 0)
        if d:
            unbalanced[str(size)] = d
    return dict(stabilizer_order=len(stab), chi=chi, twist=H.twist(B),
                orbits={f"{s}:{c}": n for (s, c), n in sorted(orbits.items())},
                orbit_sizes=sorted({s for s, _ in orbits}),
                all_sizes_divisible_by_27=all(s % 27 == 0 for s, _ in orbits),
                unbalanced_by_size=unbalanced,
                net_cells=[dict(L1=list(k[0]), L2=list(k[1]), M=list(k[2]), dims=[k[3], k[4]], net=v)
                           for k, v in sorted(net.items()) if v])


def run():
    reps = {int(k): v for k, v in json.loads(REPS.read_text())["representatives"].items()}
    res = dict(pass_id=11222, orbitals={})
    for o, v in sorted(reps.items()):
        if v["tau_image"] == o:
            continue
        a = analyse(np.array(v["B"], np.int64))
        a["size"] = v["size"]
        res["orbitals"][str(o)] = a
        print(o, v["size"], {k: a[k] for k in ("stabilizer_order", "chi", "twist", "orbit_sizes",
                                                "all_sizes_divisible_by_27", "unbalanced_by_size")}, flush=True)
    p256 = [a for a in res["orbitals"].values() if a["size"] == 256]
    res["chi_is_54Q_on_256"] = all(a["chi"] == 54 * (1 if a["twist"] == 1 else -1) for a in p256)
    res["divisible_by_27_on_256"] = all(a["all_sizes_divisible_by_27"] for a in p256)
    return res


def main():
    res = run()
    OUT.write_text(json.dumps(res, indent=1, default=int))
    print(json.dumps({k: v for k, v in res.items() if k != "orbitals"}, indent=1))


if __name__ == "__main__":
    main()
