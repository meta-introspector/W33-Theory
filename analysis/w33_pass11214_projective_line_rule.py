"""Pass 11214: the sign structures of AME(6,3) and AME(10,3) are cross-ratio geometry on PG(1,5) and PG(1,9).

Passes 11190 and 11200 found that the sign structures of AME(6,3) and AME(10,3) stabilizer states are unique with
symmetry PGL(2,5) and PGL(2,9).  This pass writes them down explicitly on the projective lines.

AME(6,3) (six qutrits): label the parties by PG(1,5).  Every 4-set U of parties splits 2|2 (Pass 11200); the split
is the HARMONIC pairing -- the unique pairing {a,b}|{c,d} of U with cross-ratio (a,b;c,d) = -1.  Exactly |PGL(2,5)| =
120 of the 720 labellings do this.

AME(10,3) (ten qutrits): label the parties by PG(1,9) = F9 u {inf}, F9 = F3[i], i^2 = -1.
  * the uniform 4-sets (the Steiner system S(3,4,10), Pass 11190) are exactly the Baer sublines PG(1,3): the harmonic
    quadruples (cross-ratio in F3);
  * for a 4-set K that is not a subline, its three pairings have cross-ratios in the three classes {+-i}, {1+i, -1+i},
    {1-i, -1-i} of F9 \\ F3; two points x, y outside K lie in the same sign class of the 6-set K^c iff {x, y} is
    harmonic to exactly one pair {a,b} of K and the pairing {a,b}|{c,d} of K is NOT of class {1+i, -1+i}.
    The rule distinguishes 1+i from its Frobenius conjugate 1-i: the sign structure is chiral, which is why its
    symmetry is PGL(2,9) and not PGamma L(2,9); the second structure on the same Steiner system (Pass 11190: exactly two)
    is the Frobenius image.
Checked on every AME(6,3) sample and all 71 AME(10,3) graph states.
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
import w33_pass11190_ame10_sign_structure as P10  # noqa: E402
import w33_pass11200_sign_method_ame6_ame8 as P6  # noqa: E402

OUT = ROOT / "data" / "w33_pass11214_projective_line_rule.json"


class Field:
    """F_p (p prime) or F9 = F3[i]; elements are ints; INF = size"""

    def __init__(self, p, ext=False):
        self.p, self.ext = p, ext
        self.size = p * p if ext else p
        self.INF = self.size
        self.ONE = self.enc(1, 0)

    def enc(self, a, b=0):
        return (a % self.p) + self.p * (b % self.p) if self.ext else a % self.p

    def dec(self, x):
        return (x % self.p, x // self.p) if self.ext else (x, 0)

    def add(self, x, y):
        (a, b), (c, d) = self.dec(x), self.dec(y)
        return self.enc(a + c, b + d)

    def neg(self, x):
        a, b = self.dec(x)
        return self.enc(-a, -b)

    def mul(self, x, y):
        (a, b), (c, d) = self.dec(x), self.dec(y)
        return self.enc(a * c - b * d, a * d + b * c) if self.ext else self.enc(a * c)

    def inv(self, x):
        return next(y for y in range(1, self.size) if self.mul(x, y) == self.ONE)

    def cr(self, a, b, c, d):
        """cross-ratio (a,b;c,d) = (c-a)(d-b) / ((c-b)(d-a)), homogeneous handling of infinity"""
        vals = []
        for x, y in ((c, a), (d, b), (c, b), (d, a)):
            vals.append(None if self.INF in (x, y) else self.add(x, self.neg(y)))
        n1, n2, d1, d2 = [self.ONE if v is None else v for v in vals]
        return self.mul(self.mul(n1, n2), self.inv(self.mul(d1, d2)))


def ame6_check(samples=12, seed=11214):
    K5 = Field(5)
    pts = list(range(6))                     # 0..4 and INF = 5
    rng = np.random.default_rng(seed)
    counts = []
    for _ in range(samples):
        st = P6.structure(P6.random_ame(6, rng), 3)
        hits = 0
        for perm in itertools.permutations(pts):
            ok = True
            for U, cls in st.items():
                img = [perm[x] for x in U]
                a, b, c, d = img
                harm = [pr for pr in (((a, b), (c, d)), ((a, c), (b, d)), ((a, d), (b, c)))
                        if K5.cr(pr[0][0], pr[0][1], pr[1][0], pr[1][1]) == K5.neg(K5.ONE)]
                pairing = {frozenset(perm[x] for x in cls), frozenset(img) - frozenset(perm[x] for x in cls)}
                if len(harm) != 1 or {frozenset(harm[0][0]), frozenset(harm[0][1])} != pairing:
                    ok = False
                    break
            hits += ok
        counts.append(hits)
    return counts


def ame10_rule(G, K9):
    L = P10.stabiliser(G)
    signs = {K: P10.signs(L, K) for K in itertools.combinations(range(10), 4)}
    blocks = [K for K, c in signs.items() if len(set(c.values())) == 1]
    m1 = K9.neg(K9.ONE)

    def harmonic(Q):
        a, b, c, d = Q
        return any(K9.cr(*p) == m1 for p in ((a, b, c, d), (a, c, b, d), (a, d, b, c)))
    labellings = []
    others = [x for x in range(10) if x not in (0, K9.ONE, K9.INF)]
    for perm in itertools.permutations(others):
        s = {0: 0, 1: K9.ONE, 2: K9.INF}
        s.update(zip(range(3, 10), perm))
        if all(harmonic([s[x] for x in B]) for B in blocks):
            labellings.append(s)
    results = []
    for s in labellings:
        results.append(tuple(rule_holds(s, signs, blocks, K9, bad) for bad in
                             ({K9.enc(1, 1), K9.enc(-1, 1)}, {K9.enc(1, -1), K9.enc(-1, -1)})))
    return len(labellings), results


def rule_holds(s, signs, blocks, K9, bad_class):
    m1 = K9.neg(K9.ONE)
    if True:
        ok = True
        for K, c in signs.items():
            if K in blocks:
                continue
            Kp = [s[k] for k in K]
            U = [x for x in range(10) if x not in K]
            for x, y in itertools.combinations(U, 2):
                X, Y = s[x], s[y]
                hits = [(a, b) for a, b in itertools.combinations(range(4), 2) if K9.cr(Kp[a], Kp[b], X, Y) == m1]
                if len(hits) == 1:
                    a, b = hits[0]
                    cc, dd = [i for i in range(4) if i not in (a, b)]
                    lam = K9.cr(Kp[a], Kp[b], Kp[cc], Kp[dd])
                    pred = lam not in bad_class and K9.inv(lam) not in bad_class
                else:
                    pred = False
                if pred != (c[x] == c[y]):
                    ok = False
                    break
            if not ok:
                break
        return ok


def run():
    res = dict(pass_id=11214)
    res["ame6_harmonic_labellings_per_sample"] = ame6_check()
    K9 = Field(3, ext=True)
    D = json.load(open(P10.TABU))
    per = Counter()
    for w in D["graphs"]:
        n_lab, ok = ame10_rule(F.to_mat(np.array(w)), K9)
        per[(n_lab, tuple(ok))] += 1
    res["ame10_labellings_and_rule"] = {f"{k[0]} labellings, rule holds {k[1]}": v for k, v in per.items()}
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
