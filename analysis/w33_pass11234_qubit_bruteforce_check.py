"""Pass 11234: an independent brute-force check of the qubit protection numbers (Pass 11220) and of A = n - c.

The parallel session's Pass 11217 proved the arrow law A(S) = n - c(S) for qubits. Its Pass 11220 computed the mean
number of protected qubits E_n[c] from GAP's conjugacy classes and the Jordan closed form
c = m1(1)/2 + chi2 m2(1) + m1(x^2+x+1), with E_2[c] = 29/40 and E_3[c] = 5981/8960; its brute-force check at n = 6 had to
be stopped.  Following the repo rule that a cross-check must be independent, this pass uses neither the closed form nor
the class list:
  * every element of Sp(4,2) (720) and Sp(6,2) (1 451 520) is enumerated by breadth-first search from transvections;
  * c(S) is computed directly: the S-invariant nondegenerate planes of F2^{2n} (span(v, Sv) with S^2 v in it, or two
    fixed vectors), then the largest mutually orthogonal family;
  * A(S) is computed directly as the minimum over ALL qubit factorisations of sum_q rank S_F[others, q]
    (10 for two qubits, 1120 for three), on every element of Sp(4,2) and a random sample of Sp(6,2).
"""

from __future__ import annotations

import itertools
import json
import random
from collections import Counter
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "w33_pass11234_qubit_bruteforce_check.json"


def omega(x, y, n):
    s = 0
    for k in range(n):
        x0, x1 = (x >> (2 * k)) & 1, (x >> (2 * k + 1)) & 1
        y0, y1 = (y >> (2 * k)) & 1, (y >> (2 * k + 1)) & 1
        s ^= (x0 & y1) ^ (x1 & y0)
    return s


def apply(M, v):
    """M = tuple of column images (ints); v an int vector"""
    out, k = 0, 0
    while v:
        if v & 1:
            out ^= M[k]
        v >>= 1
        k += 1
    return out


def compose(A, B):
    return tuple(apply(A, b) for b in B)


def transvection(v, n):
    D = 2 * n
    return tuple((1 << j) ^ (v if omega(1 << j, v, n) else 0) for j in range(D))


def enumerate_sp(n):
    D = 2 * n
    gens = [transvection(1 << j, n) for j in range(D)] + [transvection((1 << j) | (1 << (j + 2)), n) for j in range(D - 2)]
    ident = tuple(1 << j for j in range(D))
    seen = {ident}
    frontier = [ident]
    while frontier:
        new = []
        for M in frontier:
            for g in gens:
                N = compose(g, M)
                if N not in seen:
                    seen.add(N)
                    new.append(N)
        frontier = new
    return list(seen)


def protected(M, n):
    D = 2 * n
    vecs = range(1, 1 << D)
    planes = set()
    fixed = []
    for v in vecs:
        w = apply(M, v)
        if w == v:
            fixed.append(v)
            continue
        w2 = apply(M, w)
        if w2 in (0, v, w, v ^ w) and omega(v, w, n):
            planes.add(frozenset((v, w, v ^ w)))
    for a, b in itertools.combinations(fixed, 2):
        if omega(a, b, n):
            planes.add(frozenset((a, b, a ^ b)))
    P = list(planes)
    orth = [[all(omega(x, y, n) == 0 for x in P[i] for y in P[j]) for j in range(len(P))] for i in range(len(P))]
    best = [0]

    def rec(k, cand):
        if k > best[0]:
            best[0] = k
        if best[0] == n or k + len(cand) <= best[0]:
            return
        for idx, c in enumerate(cand):
            rec(k + 1, [d for d in cand[idx + 1:] if orth[c][d]])
    rec(0, list(range(len(P))))
    return best[0]


def factorisations(n):
    """all ordered-free splits into n mutually orthogonal hyperbolic planes, each given by a basis pair (u, v)"""
    D = 2 * n
    planes = {}
    for u in range(1, 1 << D):
        for v in range(u + 1, 1 << D):
            if omega(u, v, n):
                planes.setdefault(frozenset((u, v, u ^ v)), (u, v))
    plist = list(planes.items())
    out = set()

    def rec(chosen, start):
        if len(chosen) == n:
            out.add(frozenset(chosen))
            return
        for i in range(start, len(plist)):
            S, _ = plist[i]
            if all(all(omega(x, y, n) == 0 for x in S for y in T) for T in chosen):
                rec(chosen + [S], i + 1)
    rec([], 0)
    return [[planes[S] for S in F] for F in out]


def rank2(vs):
    """rank over F2 of a list of ints"""
    basis = []
    for v in vs:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
    return len(basis)


def arrow(M, facs, n):
    best = 10 ** 9
    for F in facs:
        tot = 0
        for q, (u, v) in enumerate(F):
            others = [p for i, p in enumerate(F) if i != q]
            # component of S u, S v along the other planes: express in the adapted basis via omega pairings
            imgs = [apply(M, u), apply(M, v)]
            comps = []
            for x in imgs:
                c = 0
                for k, (a, b) in enumerate(others):
                    # coordinates of x on plane (a, b): x = alpha a + beta b + ..., alpha = omega(x, b), beta = omega(a, x)
                    c |= (omega(x, b, n) << (2 * k)) | (omega(a, x, n) << (2 * k + 1))
                comps.append(c)
            tot += rank2(comps)
            if tot >= best:
                break
        best = min(best, tot)
    return best


def run(sample=1500, seed=11234):
    res = dict(pass_id=11234)
    rng = random.Random(seed)
    for n, expect_order, ref in ((2, 720, Fraction(29, 40)), (3, 1451520, Fraction(5981, 8960))):
        G = enumerate_sp(n)
        assert len(G) == expect_order, len(G)
        dist = Counter(protected(M, n) for M in G)
        E = Fraction(sum(c * k for c, k in dist.items()), len(G))
        facs = factorisations(n)
        test = G if n == 2 else rng.sample(G, sample)
        law = Counter()
        for M in test:
            law[arrow(M, facs, n) == n - protected(M, n)] += 1
        res[f"n{n}"] = dict(order=len(G), c_distribution={str(k): v for k, v in sorted(dist.items())},
                            E_c=str(E), matches_pass11220=(E == ref), factorisations=len(facs),
                            arrow_law_checked=sum(law.values()), arrow_law_violations=law[False])
        print(n, res[f"n{n}"], flush=True)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
