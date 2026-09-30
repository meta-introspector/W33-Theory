#!/usr/bin/env python3
"""Pass 11181: which qutrit dynamics entangle in EVERY subsystem decomposition -- local dynamics only knows the primes 2
and 3.

GAP (analysis/gap/w33_pass11180_mereology_classes.g, frozen in data/w33_pass11180_gap_classes.txt): for every conjugacy
class of PSp(2n,3) acting on the stabilizer-compatible tensor factorisations (45 for n = 2, 110565 for n = 3), the number
of factorisations fixed by the class -- the subsystem splits in which the tick is LOCAL (a product of single-qutrit
gates up to relabelling).
Exact fractions of ticks local in 0 / exactly 1 / at least 2 splits:
    two qutrits:    19/45      19/48     131/720        (20 classes)
    three qutrits:  7922/12285  7/48     41143/196560   (74 classes)
So an INTRINSICALLY ENTANGLING tick -- one that no choice of subsystems makes non-entangling -- is the typical case
already for three qutrits (64.5%), and a tick that selects a unique emergent split becomes rare (19/48 -> 7/48).
WHY (the prime rule).  A tick is local in F iff it lies in the stabiliser of F, (SL(2,3)^n : S_n)/<-1>, whose order is
|PSp|/#factorisations = 25920/45 = 2^6 3^2 (n = 2) or 4585351680/110565 = 2^9 3^4 (n = 3) -- a {2,3}-group; its element
orders are exactly the orders of the classes with fixed points (recorded).  Hence:
  * every tick whose period has a prime factor >= 5 (orders 5, 7, 10, 13, 14, 15, 20, 30) is intrinsically entangling;
  * for two qutrits the stabiliser's Sylow 3-subgroup is C3 x C3 (exponent 3), so every tick of order 9 is intrinsically
    entangling as well: the fixed-point-free orders are exactly {5, 9};
  * for three qutrits the Sylow 3-subgroup is C3 wr C3 (exponent 9), so order-9 ticks split: 6 of the 8 order-9 classes
    are fixed-point-free, 2 fix 3 splits; the fixed-point-free orders are {5, 7, 9, 10, 13, 14, 15, 18, 20, 30, 36}.
Reading: which subsystems exist is decided by the dynamics (quantum mereology; Zanardi 2001, Zanardi-Lidar-Lloyd 2004,
Carroll-Singh 2021), and on the W(3,3) substrate the answer is arithmetic -- a tick has a local description in some
subsystem split only if its period is built from the primes 2 and 3 (with a cap on the power of 3 for two qutrits).
"""
from __future__ import annotations

import json
import re
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CLASSES = ROOT / "data" / "w33_pass11180_gap_classes.txt"
OUT = ROOT / "data" / "w33_pass11181_intrinsic_entanglement.json"


def parse():
    t = CLASSES.read_text()
    out = {}
    for part in re.split(r'(?=n=\d \|PSp)', t):
        m = re.search(r'n=(\d) \|PSp\|=(\d+) factorisations=(\d+) classes=(\d+)', part)
        if not m:
            continue
        n, G, nf, ncl = map(int, m.groups())
        cl = [tuple(map(int, x)) for x in re.findall(r'CLASS order=(\d+) size=(\d+) fixed=(\d+)', part)]
        out[n] = dict(order=G, factorisations=nf, classes=cl)
    return out


def summarize():
    data = parse()
    res = dict(pass_id=11181)
    for n, d in sorted(data.items()):
        G = d['order']
        cl = d['classes']
        assert sum(s for _, s, _ in cl) == G
        free = sum(s for o, s, f in cl if f == 0)
        one = sum(s for o, s, f in cl if f == 1)
        orders_free = sorted({o for o, s, f in cl if f == 0})
        orders_nonfree = sorted({o for o, s, f in cl if f > 0})
        prime_rule = all(f == 0 for o, s, f in cl if any(o % p == 0 for p in (5, 7, 11, 13)))
        nonfree_orders_are_23 = all(all(p in (2, 3) for p in _primes(o)) for o in orders_nonfree)
        # Burnside: average number of fixed factorisations = number of orbits = 1 (transitive action)
        burnside = Fraction(sum(s * f for o, s, f in cl), G)
        stab = G // d['factorisations']
        assert stab * d['factorisations'] == G
        res[f'n{n}'] = dict(psp_order=G, factorisations=d['factorisations'], classes=len(cl),
                           intrinsically_entangling=str(Fraction(free, G)), unique_split=str(Fraction(one, G)),
                           at_least_two=str(Fraction(G - free - one, G)), orders_fixed_point_free=orders_free,
                           orders_with_local_splits=orders_nonfree, prime_rule_holds=prime_rule,
                           local_orders_are_2_3_numbers=nonfree_orders_are_23, burnside_average_fixed=str(burnside),
                           stabiliser_order=stab, stabiliser_primes=sorted(_primes(stab)),
                           stabiliser_has_order_9=9 in orders_nonfree)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


def _primes(m):
    ps, p = set(), 2
    while m > 1:
        while m % p == 0:
            ps.add(p)
            m //= p
        p += 1
    return ps


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
