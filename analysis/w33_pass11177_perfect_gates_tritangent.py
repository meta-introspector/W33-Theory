#!/usr/bin/env python3
"""Pass 11177: a two-qutrit Clifford gate is perfect EXACTLY when it moves the tensor factorisation to a collinear point
of GQ(4,2) -- a factorisation sharing a complete frame, i.e. a tritangent plane meeting the original in a line of the
cubic surface.  The three space-time entanglement classes of Pass 11161 are the three relations of this geometry.

Known ingredients (Holotrade track, cited in the memory dictionary): a tensor factorisation C^9 = C^3 (x) C^3 is a pair of
orthogonal nondegenerate planes {L, L-perp} of F_3^4, i.e. an 'octet' of 8 points of W(3,3); there are 45; a complete
factorisation frame is 5 factorisations whose octets partition the 40 points; there are 27, the lines of the cubic
surface (GQ(2,4)), and the 45 factorisations are its tritangent planes.
New here (exact, all 51840 elements of Sp(4,3)):
  * the octet-disjointness graph on the 45 factorisations is SRG(45,12,3,3), the collinearity graph of GQ(4,2); every
    disjoint pair lies in exactly one of the 27 frames; each factorisation lies in 3 frames;
  * for S in Sp(4,3) with blocks S_AA, S_BA, the relation between F0 and S(F0) is determined by the Choi entanglement:
        rank (S_AA, S_BA) = (2, 2)  [perfect, I3 = -2]      <=>  S(F0) collinear with F0 (octets disjoint)  13824 gates
        rank (1, 2) or (2, 1)       [I3 = -1]               <=>  S(F0) non-collinear (octets meet)          36864 gates
        rank (2, 0) or (0, 2)       [local or swap, I3 = 0] <=>  S(F0) = F0                                  1152 gates
  * the stabiliser H of F0 (order 1152) is transitive on the 12 collinear factorisations (point stabiliser order 96), so
    the perfect gates form ONE double coset H g H;
  * GAP (analysis/gap/w33_pass11177_factorisation_orbitals.g, frozen in data/w33_pass11177_gap_orbitals.txt): Sp(4,3) on
    the 45 factorisations has rank 3, subdegrees 1, 12, 32 -- the rank-3 action of U4(2) on the 45 tritangent planes
    (the index-45 subgroup is unique up to conjugacy) -- and Sp(6,3) on the 110565 three-qutrit factorisations has rank
    20, the perfect relation being the single suborbit of size 3456 (= 128/4095 of the whole, Pass 11163).
Reading: which way the register splits into 'A' and 'B' is a tritangent plane; a perfectly scrambling tick moves the
split to a neighbouring tritangent plane, one that shares a line with it.  The E6 cubic has one monomial per tritangent
plane, so a perfect tick links two monomials that share a variable.
"""
from __future__ import annotations

import itertools
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11156_perfect_spacetime_gates as G  # noqa: E402

OUT = ROOT / "data" / "w33_pass11177_perfect_gates_tritangent.json"
GAPOUT = ROOT / "data" / "w33_pass11177_gap_orbitals.txt"
J = np.array([[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]])
V = [np.array(v) for v in itertools.product(range(3), repeat=4) if any(v)]


def norm(v):
    v = np.array(v) % 3
    i = next(k for k in range(4) if v[k])
    return tuple((v * pow(int(v[i]), -1, 3)) % 3)


PTS = sorted({norm(v) for v in V})
PID = {p: i for i, p in enumerate(PTS)}


def w(u, v):
    return int(u @ J @ v) % 3


def span_pts(B):
    s = set()
    for a, b in itertools.product(range(3), repeat=2):
        v = (a * B[0] + b * B[1]) % 3
        if v.any():
            s.add(PID[norm(v)])
    return frozenset(s)


def perp(P):
    vs = [np.array(PTS[i]) for i in P]
    return frozenset(i for i, p in enumerate(PTS) if all(w(np.array(p), q) == 0 for q in vs))


def factorisations():
    planes = {span_pts(np.array([u, v])) for u in V for v in V if w(u, v)}
    facs = list({frozenset([P, perp(P)]) for P in planes})
    return len(planes), facs


def summarize():
    n_planes, facs = factorisations()
    octs = [frozenset().union(*f) for f in facs]
    n = len(facs)
    A = np.array([[1 if i != j and not (octs[i] & octs[j]) else 0 for j in range(n)] for i in range(n)])
    A2 = A @ A
    lam = sorted({int(A2[i, j]) for i in range(n) for j in range(n) if A[i, j]})
    mu = sorted({int(A2[i, j]) for i in range(n) for j in range(n) if i != j and not A[i, j]})
    frames = [c for c in itertools.combinations(range(n), 5) if all(A[a, b] for a, b in itertools.combinations(c, 2))]
    edge_frames = Counter(p for c in frames for p in itertools.combinations(c, 2))
    per_fac = Counter(i for c in frames for i in c)
    e1 = np.array([[1, 0, 0, 0], [0, 1, 0, 0]])
    e2 = np.array([[0, 0, 1, 0], [0, 0, 0, 1]])
    F0 = frozenset([span_pts(e1), span_pts(e2)])
    i0 = facs.index(F0)
    census = Counter()
    H, perfect_images = [], Counter()
    for S in G.sp43():
        img = frozenset([span_pts((S @ e1.T).T % 3), span_pts((S @ e2.T).T % 3)])
        j = facs.index(img)
        rel = 'equal' if j == i0 else ('collinear' if A[i0, j] else 'noncollinear')
        census[f"{G.rank3(S[:2, :2])}{G.rank3(S[2:, :2])}|{rel}"] += 1
        if j == i0:
            H.append(S)
    nbrs = [j for j in range(n) if A[i0, j]]
    target = nbrs[0]

    # H-orbit of one collinear neighbour (acting on factorisations through their octets)
    def act(S, f):
        out = []
        for P in f:
            pts = [np.array(PTS[p]) for p in P]
            out.append(frozenset(PID[norm((S @ q) % 3)] for q in pts))
        return frozenset(out)
    orbit = {act(S, facs[target]) for S in H}
    stab = sum(1 for S in H if act(S, facs[target]) == facs[target])
    gap = GAPOUT.read_text() if GAPOUT.exists() else ''
    res = dict(pass_id=11177, nondegenerate_planes=n_planes, factorisations=n, degree=sorted(set(A.sum(1).tolist())),
               srg_lambda=lam, srg_mu=mu, frames=len(frames), frames_partition_40=all(len(frozenset().union(*[octs[i] for i in c])) == 40 for c in frames),
               disjoint_pairs=int(A.sum() // 2), pairs_per_frame_count=dict(Counter(edge_frames.values())),
               frames_per_factorisation=sorted(set(per_fac.values())), census=dict(sorted(census.items())),
               H_order=len(H), H_orbit_on_collinear=len(orbit), H_stabiliser_of_neighbour=stab,
               gap_rank3_45='subdegrees=[ 1, 12, 32 ]' in gap, gap_rank20_110565='rank=20' in gap and '3456' in gap)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
