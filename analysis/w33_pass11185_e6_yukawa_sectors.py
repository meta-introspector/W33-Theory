#!/usr/bin/env python3
"""Pass 11185: the E6 Yukawa reading of two-qutrit splits and gate classes -- a counting dictionary, no dynamics claimed.

Known (Holotrade dictionary, cited in memory): the 27 complete factorisation frames are the 27 lines of the cubic surface
= the 27 weights of E6; relative to one frame Phi they split 1 + 10 + 16 (frames sharing all / one / no factorisation),
the SO(10) decomposition of one E6 generation.  The 45 factorisations are the tritangent planes = the 45 monomials of
the E6 cubic invariant.
Checked here (exact):
  * each factorisation lies in exactly 3 frames, pairwise sharing it (a triangle of the frame graph), and every triangle
    of the frame graph is a factorisation -- the 45 cubic monomials;
  * relative to Phi: the 5 factorisations through Phi have their other two frames in the '10' -- the 5 monomials
    1.10.10; each of the other 40 has exactly one frame in the '10' and two in the '16' -- the 40 monomials 16.16.10;
  * sector transitions of the three gate classes (Passes 11161, 11177), for a split F in sector S1 (1.10.10) or S16
    (16.16.10):
        from S1:  collinear (perfect) -> 4 in S1, 8 in S16;   non-collinear (I3 = -1) -> 0 in S1, 32 in S16;
        from S16: collinear (perfect) -> 1 in S1, 11 in S16;  non-collinear (I3 = -1) -> 4 in S1, 28 in S16;
    at the gate level (base split F0 in S1): the 13824 perfect gates send F0 into S1 for 4608 and into S16 for 9216;
    the 36864 partially entangling gates never keep it in S1 (all into S16); local gates fix it;
  * the stabiliser of the frame Phi in Sp(4,3) has order 1920 (W(D5), the SO(10) Weyl group; 960 in PSp(4,3), as in the
    Holotrade dictionary).  Relative to a split F0 in the singlet sector, EVERY element is either local (384) or perfect
    (1536) -- never partially entangling -- because it maps F0 to a split of the same frame, and those are pairwise
    collinear.  768 of its elements fix no split at all (intrinsically entangling ticks inside W(D5)).
Reading, stated at the level the computation supports: a merely entangling two-qutrit tick always moves a singlet-sector
split (a 1.10.10 monomial) into the matter sector (16.16.10); only a perfect tick can move it within the singlet sector.
No claim is made about fields, charges or Yukawa values.
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
import w33_pass11177_perfect_gates_tritangent as T  # noqa: E402

OUT = ROOT / "data" / "w33_pass11185_e6_yukawa_sectors.json"


def summarize():
    npl, facs = T.factorisations()
    octs = [frozenset().union(*f) for f in facs]
    n = len(facs)
    A = np.array([[1 if i != j and not (octs[i] & octs[j]) else 0 for j in range(n)] for i in range(n)])
    frames = [c for c in itertools.combinations(range(n), 5) if all(A[a, b] for a, b in itertools.combinations(c, 2))]
    fr_of = {i: [k for k, c in enumerate(frames) if i in c] for i in range(n)}
    # frame graph: frames adjacent iff they share a factorisation
    FG = np.zeros((27, 27), int)
    for i in range(n):
        for a, b in itertools.combinations(fr_of[i], 2):
            FG[a, b] = FG[b, a] = 1
    triangles = {frozenset(t) for t in itertools.combinations(range(27), 3) if FG[t[0], t[1]] and FG[t[1], t[2]] and FG[t[0], t[2]]}
    fac_triples = {frozenset(fr_of[i]) for i in range(n)}
    phi = 0
    ten = [k for k in range(27) if FG[phi, k]]
    sixteen = [k for k in range(27) if k != phi and not FG[phi, k]]
    sector = {}
    ok_split = True
    for i in range(n):
        fs = fr_of[i]
        if phi in fs:
            sector[i] = 'S1'
            ok_split &= all(k in ten for k in fs if k != phi)
        else:
            sector[i] = 'S16'
            ok_split &= sum(k in ten for k in fs) == 1 and sum(k in sixteen for k in fs) == 2
    trans = Counter()
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            rel = 'collinear' if A[i, j] else 'noncollinear'
            trans[f"{sector[i]}->{sector[j]}|{rel}"] += 1
    per_split = {k: v // sum(1 for i in range(n) if sector[i] == k.split('->')[0]) for k, v in trans.items()}
    # gate level: base split F0 in S1
    F0 = next(i for i in range(n) if sector[i] == 'S1')
    fid = {f: i for i, f in enumerate(facs)}

    def image(S, i):
        pp = [T.PID[T.norm((S @ np.array(p)) % 3)] for p in T.PTS]
        return fid[frozenset(frozenset(pp[q] for q in P) for P in facs[i])]
    gate = Counter()
    frame_stab, frame_stab_prof = 0, Counter()
    phi_set = frozenset(frames[phi])
    for S in G.sp43():
        pp = [T.PID[T.norm((S @ np.array(p)) % 3)] for p in T.PTS]
        img = lambda i: fid[frozenset(frozenset(pp[q] for q in P) for P in facs[i])]
        j = img(F0)
        rel = 'equal' if j == F0 else ('collinear' if A[F0, j] else 'noncollinear')
        gate[f"{rel}->{sector[j]}"] += 1
        if frozenset(img(i) for i in frames[phi]) == phi_set:
            frame_stab += 1
            fixed = sum(img(i) == i for i in range(n))
            frame_stab_prof[f"fixed{fixed}|F0:{rel}"] += 1
    res = dict(pass_id=11185, factorisations=n, frames=len(frames), frames_per_factorisation=sorted({len(v) for v in fr_of.values()}),
               factorisations_are_frame_triangles=fac_triples == triangles, triangles=len(triangles),
               frame_graph_degree=sorted(set(FG.sum(1).tolist())), decomposition=[1, len(ten), len(sixteen)],
               sectors=dict(Counter(sector.values())), sector_rule_holds=bool(ok_split), transitions_per_split=per_split,
               gate_level_from_S1=dict(gate), frame_stabiliser_order=frame_stab,
               frame_stabiliser_profile=dict(frame_stab_prof))
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
