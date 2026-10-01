#!/usr/bin/env python3
"""Pass 11231: a dynamical area law -- entanglement in circuits of the substrate's ticks equals the spacetime min-cut.

Pass 11226 found the static Ryu-Takayanagi law S(A) = min-cut exactly in networks of the perfect tick p.  A dynamical
space-time needs time evolution.  Here a chain of N qutrits starts in the product state |0...0> and evolves by a
brick-wall circuit of a fixed tick.  The circuit is itself a tensor network in spacetime, so for every region A of the
output, S(A) <= min-cut: the least number of worldline segments that must be cut to separate A's output legs from the
rest (initial product states are free ends).  We compute S(A) exactly (stabilizer ranks, trits) for every contiguous
interval at every depth and compare.

Ticks: the perfect tick p (Pass 11193; AME(4,3) Choi state), SUM (not perfect), random two-qutrit Cliffords, and the
three-qutrit perfect interacting tick VKV of Pass 11182 (the tick that carries an arrow, A = 2, Pass 11223) in a
three-site brick pattern.
Outcome per tick: the fraction of (interval, depth) pairs with S(A) = min-cut, and the entanglement velocity of a
half-chain cut.  Results (N = 12, depth 10, 570 pairs): no tick saturates the min-cut everywhere.  Random perfect gates
(p dressed by random local Cliffords) from solvable Bell-pair states come closest, 548/570, with the half-chain entropy
growing at the maximal rate to saturation; the arrow-carrying three-qutrit tick VKV from product states gives 534-539;
the bare perfect tick p is fine-tuned (CSS-like: |0...0> stays a product, 15/570); SUM is worst.  Deficits are not only
at the open boundary.  So circuits of the geometry's ticks obey an approximate, not exact, dynamical area law.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import networkx as nx
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11226_holographic_area_law as H  # noqa: E402

OUT = ROOT / "data" / "w33_pass11231_dynamical_area_law.json"


def apply_gate(L, S, sites):
    cols = [2 * s + t for s in sites for t in range(2)]
    L = L.copy()
    L[:, cols] = (L[:, cols] @ S.T) % 3
    return L


def brick_layers(N, k, depth):
    """k-site bricks; layer t starts at offset t mod k"""
    layers = []
    for t in range(depth):
        off = t % k
        layers.append([list(range(s, s + k)) for s in range(off, N - k + 1, k)])
    return layers


def circuit_mincut(N, layers, A, bell_pairs=()):
    G = nx.Graph()
    last = {q: ("init", q) for q in range(N)}
    for a, b in bell_pairs:                          # an initial Bell pair is one tensor joining two worldlines
        last[a] = last[b] = ("bell", a)
    for t, layer in enumerate(layers):
        for gi, sites in enumerate(layer):
            node = ("g", t, gi)
            for q in sites:
                u = last[q]
                if G.has_edge(u, node):
                    G[u][node]["capacity"] += 1
                else:
                    G.add_edge(u, node, capacity=1)
                last[q] = node
    for q in range(N):
        tgt = "SRC" if q in A else "SNK"
        u = last[q]
        if G.has_edge(u, tgt):
            G[u][tgt]["capacity"] += 1
        else:
            G.add_edge(u, tgt, capacity=1)
    return int(nx.minimum_cut(G, "SRC", "SNK")[0])


DIRS = [(1, 0), (0, 1), (1, 1), (1, 2)]               # the four single-qutrit stabilizer directions X, Z, XZ, XZ^2


def run_tick(name, gates, k, N=12, depth=10, seed=0, init="xz"):
    """init: 'z' (|0...0>, a Z eigenstate), 'xz' (every qutrit an XZ eigenstate) or 'random' (random directions)"""
    rng = np.random.default_rng(seed)
    L = np.zeros((N, 2 * N), np.int64)
    bell = ()
    if init == "bell":
        # Bell pairs on (1,2), (3,4), ...: offset from the first brick layer (0,1), (2,3), ...
        bell = tuple((q, q + 1) for q in range(1, N - 1, 2))
        r = 0
        for a, b in bell:
            L[r, 2 * a], L[r, 2 * b] = 1, 1                       # X_a X_b
            L[r + 1, 2 * a + 1], L[r + 1, 2 * b + 1] = 1, 2       # Z_a Z_b^-1
            r += 2
        for q in (0, N - 1):
            L[r, 2 * q + 1] = 1
            r += 1
    else:
        for q in range(N):
            d = (0, 1) if init == "z" else ((1, 1) if init == "xz" else DIRS[rng.integers(4)])
            L[q, 2 * q], L[q, 2 * q + 1] = d
    layers = brick_layers(N, k, depth)
    equal = total = 0
    below = 0
    deficits = []
    half = []
    for t in range(depth):
        for sites in layers[t]:
            S = gates[rng.integers(len(gates))] if len(gates) > 1 else gates[0]
            L = apply_gate(L, S, sites)
        for a in range(N):
            for b in range(a + 1, N + 1):
                if b - a > N // 2:
                    continue
                A = list(range(a, b))
                s = H.entropy(L, A)
                c = circuit_mincut(N, layers[:t + 1], set(A), bell)
                total += 1
                equal += s == c
                below += s < c
                if s < c:
                    deficits.append((a, b, t + 1, int(s), int(c)))
                assert s <= c, "entropy cannot exceed the min-cut"
        half.append(H.entropy(L, list(range(N // 2))))
    return dict(tick=name, init=init, N=N, depth=depth, brick=k, pairs=total, rt_exact=equal, entropy_below_cut=below,
                fraction_exact=equal / total, half_chain_entropy=half,
                deficits_touching_boundary=sum(1 for a, b, *_ in deficits if a == 0 or b == N),
                deficits_interior=sum(1 for a, b, *_ in deficits if a > 0 and b < N),
                deficit_depths=sorted({d[2] for d in deficits}))


def random_cliffords(count, seed=1):
    sys.path.insert(0, str(ROOT / "analysis"))
    import w33_pass11223_named_tick_arrows as N
    G, _ = N.sp43()
    rng = np.random.default_rng(seed)
    return [G[i] for i in rng.choice(len(G), count, replace=False)]


def dressed_perfect(p, count, seed=7):
    """local (x) p (x) local with random single-qutrit Cliffords: still perfect (AME(4,3) Choi state)"""
    SL = [np.array([[a, b], [c, d]]) for a in range(3) for b in range(3) for c in range(3) for d in range(3)
          if (a * d - b * c) % 3 == 1]
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(count):
        l = np.zeros((4, 4), np.int64)
        r = np.zeros((4, 4), np.int64)
        l[:2, :2], l[2:, 2:] = SL[rng.integers(24)], SL[rng.integers(24)]
        r[:2, :2], r[2:, 2:] = SL[rng.integers(24)], SL[rng.integers(24)]
        out.append((l @ p @ r) % 3)
    return out


def run():
    import w33_pass11193_optimal_perfect_gate_compiler as C
    import w33_pass11182_paper_ticks_mereology as PT
    p = C.fixed_perfect_gate() % 3
    SUM = np.array([[1, 0, 0, 0], [0, 1, 0, 2], [1, 0, 1, 0], [0, 0, 0, 1]], np.int64).T
    vkv = PT.ticks()["perfect_tick_VKV"] % 3
    runs = []
    for init in ("z", "xz", "random", "bell"):
        runs += [run_tick("perfect tick p", [p], 2, init=init, seed=5),
                 run_tick("SUM", [SUM], 2, init=init, seed=5),
                 run_tick("random two-qutrit Cliffords", random_cliffords(200), 2, seed=3, init=init),
                 run_tick("random perfect gates (dressed p)", dressed_perfect(p, 200), 2, seed=4, init=init),
                 run_tick("three-qutrit perfect tick VKV", [vkv], 3, init=init, seed=5)]
    res = dict(pass_id=11231, runs=runs)
    return res


def main():
    res = run()
    OUT.write_text(json.dumps(res, indent=1))
    for r in res["runs"]:
        print(r["init"], r["tick"], r["rt_exact"], "/", r["pairs"], "below:", r["entropy_below_cut"], "half-chain S(t):",
              r["half_chain_entropy"])


if __name__ == "__main__":
    main()
