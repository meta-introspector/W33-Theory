"""Pass 11199: the local geometry of AME(10,3) is a regular spread of PG(3,3) whose reguli are the Steiner blocks.

Pass 11190 showed that in every AME(10,3) stabilizer state the uniform 4-sets of parties form a Steiner system
S(3,4,10) -- the inversive plane of order 3 -- and asked (open) whether the ten parties can be identified with the ten
lines of a regular spread of PG(3,3) (Thas 1997: the reguli of a regular spread form the Miquelian inversive plane).

Test, for every 3-set A of parties.  W_A = stabilisers trivial on A is 4-dim, so PG(W_A) = PG(3,3).  For each of the
seven parties x outside A, the kernel of restriction to x is a line L_x; the seven lines are pairwise skew (Pass 11190).
  (1) the seven lines extend to a spread (10 pairwise skew lines covering the 40 points) in how many ways?
  (2) is that spread regular (every regulus through 3 of its lines lies in it)?
  (3) is there a labelling of the 3 remaining lines by the 3 parties of A such that for EVERY Steiner block the
      four corresponding lines form a regulus -- i.e. the local spread IS the inversive plane of the state?
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
import w33_pass11175_scan_five_qutrit as F  # noqa: E402
import w33_pass11190_ame10_sign_structure as P  # noqa: E402

OUT = ROOT / "data" / "w33_pass11199_ame10_local_spreads.json"
NQ = 10
V4 = [np.array(v) for v in itertools.product(range(3), repeat=4) if any(v)]


def norm(v):
    i = next(k for k in range(len(v)) if v[k] % 3)
    return tuple(int(x) for x in (v * pow(int(v[i]), -1, 3)) % 3)


POINTS = sorted({norm(v) for v in V4})
PIDX = {p: i for i, p in enumerate(POINTS)}


def line_of(u, v):
    return frozenset(PIDX[norm((a * u + b * v) % 3)] for a in range(3) for b in range(3) if ((a * u + b * v) % 3).any())


LINES = sorted({line_of(np.array(POINTS[i]), np.array(POINTS[j])) for i in range(40) for j in range(i + 1, 40)},
               key=sorted)
assert len(POINTS) == 40 and len(LINES) == 130


def regulus(l1, l2, l3):
    trans = [m for m in LINES if all(len(m & l) == 1 for l in (l1, l2, l3))]
    return frozenset(m for m in LINES if all(len(m & t) == 1 for t in trans)) if len(trans) == 4 else None


def local_lines(L, A):
    P7 = [q for q in range(NQ) if q not in A]
    R = list(A) + P7[:2]
    cols = [2 * q + t for q in R for t in range(2)]
    Ri = F.inv_mod3(L[:, cols])
    W = Ri[-4:, :] @ L % 3                                   # 4 x 20: basis of W_A
    assert not W[:, [2 * q + t for q in A for t in range(2)]].any()
    lines = {}
    for x in P7:
        Wx = W[:, [2 * x, 2 * x + 1]]
        ker = [PIDX[norm(v)] for v in V4 if not ((v @ Wx) % 3).any()]
        lines[x] = frozenset(ker)
        assert len(lines[x]) == 4
    return lines


def analyse(G, blocks):
    L = P.stabiliser(G)
    stats = Counter()
    for A in itertools.combinations(range(NQ), 3):
        lines = local_lines(L, A)
        Ls = list(lines.values())
        assert all(not (a & b) for a, b in itertools.combinations(Ls, 2))
        free = set(range(40)) - set().union(*Ls)
        cand = [m for m in LINES if m <= free]
        spreads = [c for c in itertools.combinations(cand, 3) if len(set().union(*c)) == 12]
        stats[f"spreads_{len(spreads)}"] += 1
        if len(spreads) != 1:
            continue
        spread = Ls + list(spreads[0])
        regular = all(regulus(*t) <= frozenset(spread) for t in itertools.combinations(spread, 3))
        stats[f"regular_{regular}"] += 1
        # labellings of the three extra lines by the parties of A making every block a regulus
        good = 0
        for perm in itertools.permutations(spreads[0]):
            lab = dict(lines)
            lab.update(dict(zip(A, perm)))
            if all(lab[K[3]] in regulus(lab[K[0]], lab[K[1]], lab[K[2]]) for K in blocks):
                good += 1
        stats[f"block_labellings_{good}"] += 1
        # reguli among the seven kernel lines versus Steiner blocks and the sign split of the complementary 6-set
        P7 = sorted(lines)
        for K in itertools.combinations(P7, 4):
            reg = lines[K[3]] in regulus(lines[K[0]], lines[K[1]], lines[K[2]])
            sg = P.signs(L, K)
            U = sorted(sg)
            if len(set(sg.values())) == 1:
                kind = "block"
            else:
                cls = {q for q in U if sg[q] == sg[U[0]]}
                kind = f"split_A{min(len(cls & set(A)), 3 - len(cls & set(A)))}"
            stats[f"regulus_{reg}|{kind}"] += 1
    return dict(stats)


def run(n_graphs=71):
    D = json.load(open(P.TABU))
    out = Counter()
    for w in D["graphs"][:n_graphs]:
        G = F.to_mat(np.array(w))
        L = P.stabiliser(G)
        blocks = []
        for K in itertools.combinations(range(NQ), 4):
            c = P.signs(L, K)
            if len(set(c.values())) == 1:
                blocks.append(K)
        assert len(blocks) == 30
        out.update(analyse(G, blocks))
    return dict(pass_id=11199, graphs=n_graphs, three_sets=120 * n_graphs, counts=dict(out))


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1)
    print(res)


if __name__ == "__main__":
    main()
