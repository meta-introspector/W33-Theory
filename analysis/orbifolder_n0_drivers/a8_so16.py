"""SO(16)xSO(16) (Witten Z2) completions of the W(3,3) A8 Kac pair of T6/Z3.
V1 = (1/6^7, 5/6; 0^7, 2/3).  V0 = (a; b) with a, b in (1/2)*{norm-4 vectors of E8} (one Weyl orbit, 2160 each, so every
choice gives SO16 x SO16 in 10D).  Condition V0.V1 in Z (non-SUSY orbifolder, eq. (1): V0.Vi = 0 mod 1).
For every admissible V0: the 4D gauge group = E8xE8 roots with p.V0 in Z and p.V1 in Z (counted, with rank-by-factor via
connected components of the root system)."""
import itertools
import json
from collections import Counter
from fractions import Fraction as F


def e8_norm(n2):
    out = []
    for half in (False, True):
        rng = [F(k, 2) for k in range(-5, 6, 2)] if half else [F(k) for k in range(-2, 3)]
        for v in itertools.product(rng, repeat=8):
            if sum(x * x for x in v) == n2 and sum(v) % 2 == 0:
                out.append(v)
    return out


ROOTS8 = e8_norm(2)
NORM4 = e8_norm(4)
assert len(ROOTS8) == 240 and len(NORM4) == 2160
V1a = (F(1, 6),) * 7 + (F(5, 6),)
V1b = (F(0),) * 7 + (F(2, 3),)
dot = lambda u, v: sum(a * b for a, b in zip(u, v))
HALF = [tuple(x / 2 for x in v) for v in NORM4]


def kept(a, V):
    return [r for r in ROOTS8 if dot(r, a).denominator == 1 and dot(r, V).denominator == 1]


def components(roots):
    """simple-ish classification: number of roots and rank of the root system, by connected components of the
    (non-orthogonality) graph"""
    roots = list(roots)
    idx = {r: i for i, r in enumerate(roots)}
    seen, comps = set(), []
    for r in roots:
        if r in seen:
            continue
        stack, comp = [r], []
        seen.add(r)
        while stack:
            x = stack.pop()
            comp.append(x)
            for y in roots:
                if y not in seen and dot(x, y) != 0:
                    seen.add(y)
                    stack.append(y)
        comps.append(len(comp))
    return sorted(comps, reverse=True)


NAMES = {2: 'A1', 6: 'A2', 12: 'A3', 20: 'A4', 30: 'A5', 42: 'A6', 56: 'A7', 72: 'A8', 24: 'D4', 40: 'D5', 60: 'D6', 84: 'D7',
         112: 'D8', 72.1: 'E6', 126: 'E7', 240: 'E8'}


def name(comps):
    return ' x '.join(NAMES.get(c, f'[{c}]') for c in comps) or '-'


oka = [a for a in HALF if dot(a, V1a).denominator in (1, 2, 3, 6)]
res = Counter()
examples = {}
fa = Counter(dot(a, V1a) % 1 for a in HALF)
fb = Counter(dot(b, V1b) % 1 for b in HALF)
print('a.V1a mod 1 distribution', {str(k): v for k, v in fa.items()})
print('b.V1b mod 1 distribution', {str(k): v for k, v in fb.items()})
admissible = [(a, b) for a in HALF for b in HALF if (dot(a, V1a) + dot(b, V1b)).denominator == 1]
print('admissible V0 = (a;b):', len(admissible), 'of', len(HALF) ** 2)
groups = Counter()
ga, gb = {}, {}
for a in {a for a, _ in admissible}:
    ga[a] = name(components(kept(a, V1a)))
for b in {b for _, b in admissible}:
    gb[b] = name(components(kept(b, V1b)))
for a, b in admissible:
    g = (ga[a], gb[b])
    groups[g] += 1
    examples.setdefault(g, [str(x) for x in a + b])
for g, c in groups.most_common():
    print(c, g, 'e.g. V0 =', examples[g])
json.dump(dict(admissible=len(admissible), groups={' | '.join(g): c for g, c in groups.items()},
               examples={' | '.join(g): e for g, e in examples.items()}),
          open(r'C:\Users\wiljd\AppData\Local\Temp\claude\c--Repos-Theory-of-Everything\ad42d157-70aa-40e7-80cd-3fad6695948b\scratchpad\p1109x\a8_so16.json', 'w'), indent=1)
