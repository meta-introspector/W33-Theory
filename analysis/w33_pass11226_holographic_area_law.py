#!/usr/bin/env python3
"""Pass 11226: an exact area law -- Ryu-Takayanagi in qutrit holographic codes built from perfect tensors.

Perfect tensors (AME states: every half of the legs maximally entangled with the other half) are the building blocks of
the holographic codes of Pastawski, Yoshida, Harlow and Preskill (HaPPY, 2015), in which the entanglement entropy of a
boundary region equals the size of a minimal cut through the network: a discrete Ryu-Takayanagi formula, S = |gamma|
(in units of log q), the information-theoretic face of 'area = entropy'.

This pass builds such codes from qutrit perfect tensors of the kind the corpus already has (AME(6,3); the AME(4,3)
perfect gate of Pass 11193; the AME(10,3) Glynn state of Pass 11224) and checks RT EXACTLY:
  * stabilizer contraction over F_3 (projecting contracted leg pairs onto Bell pairs) gives the boundary state;
  * S(A) = rank of the stabilizer Lagrangian restricted to A minus |A| (trits), for every contiguous boundary interval;
  * the minimal cut separating A from the rest is computed by max-flow on the network graph (bonds of capacity 1,
    bulk legs fixed in a product state so they are never cut);
  * reported: the number of intervals with S(A) = min-cut, and any violations.
Results: RT holds exactly for every contiguous boundary interval on a HaPPY {5,4} patch of AME(6,3) tensors (bulk legs
in a product state), on a flat {4,4} patch of the same tensors, and on a {4,5} patch built ONLY from the substrate's own
perfect two-qutrit tick p (Pass 11193; its Choi state is AME(4,3)), which has no bulk legs.
What this does NOT show: flat grids of the same perfect tick (2x2 ... 5x5, no interior boundary legs) satisfy RT on
every interval too, so the area law does not select negative curvature and the sign of the cosmological constant is not
decided here (a first plan to contrast {5,4} with {4,4} was dropped when both satisfied RT); and a discrete area law is
not Einstein's equations -- those follow from the entanglement first law only for holographic CFTs.
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import networkx as nx
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "w33_pass11226_holographic_area_law.json"


# ---------------------------------------------------------------- F_3 linear algebra
def rref3(M):
    M = np.array(M, np.int64) % 3
    if M.ndim == 1:
        M = M[None]
    r = 0
    for c in range(M.shape[1]):
        piv = [i for i in range(r, M.shape[0]) if M[i, c]]
        if not piv:
            continue
        M[[r, piv[0]]] = M[[piv[0], r]]
        M[r] = (M[r] * (1 if M[r, c] == 1 else 2)) % 3
        for i in range(M.shape[0]):
            if i != r and M[i, c]:
                M[i] = (M[i] - M[i, c] * M[r]) % 3
        r += 1
        if r == M.shape[0]:
            break
    return M[:r]


def rank3(M):
    return len(rref3(M)) if len(M) else 0


def nullspace3(M):
    M = np.array(M, np.int64) % 3
    R = rref3(M)
    piv = [int(np.argmax(row != 0)) for row in R]
    free = [c for c in range(M.shape[1]) if c not in piv]
    out = []
    for f in free:
        x = np.zeros(M.shape[1], np.int64)
        x[f] = 1
        for row, p in zip(R, piv):
            x[p] = (-row[f]) % 3
        out.append(x)
    return np.array(out, np.int64).reshape(len(out), M.shape[1])


# ---------------------------------------------------------------- graph-state perfect tensors
def graph_lagrangian(G):
    n = len(G)
    L = np.zeros((n, 2 * n), np.int64)
    for v in range(n):
        L[v, 2 * v] = 1
        for u in range(n):
            L[v, 2 * u + 1] = G[v, u] % 3
    return L


def entropy(L, A):
    """trits of entanglement of a pure stabilizer state (Lagrangian L, rows) across region A (list of legs)"""
    cols = [2 * a + t for a in A for t in range(2)]
    return rank3(L[:, cols]) - len(A)


def is_perfect(L, n):
    return all(entropy(L, list(A)) == n // 2 for A in itertools.combinations(range(n), n // 2))


def find_ame6(seed=1):
    rng = np.random.default_rng(seed)
    pairs = list(itertools.combinations(range(6), 2))
    for _ in range(200000):
        w = rng.integers(0, 3, len(pairs))
        G = np.zeros((6, 6), np.int64)
        for (i, j), x in zip(pairs, w):
            G[i, j] = G[j, i] = x
        if is_perfect(graph_lagrangian(G), 6):
            return G
    raise RuntimeError("no AME(6,3) graph found")


# ---------------------------------------------------------------- tensor networks of stabilizer states
def contract(L, legs, a, b):
    """project legs a, b onto the qutrit Bell pair (stabilised by X_a X_b and Z_a Z_b^-1) and drop them.
    L: Lagrangian rows over columns (x_l, z_l) for l in legs order."""
    ia, ib = legs.index(a), legs.index(b)
    cols_ab = [2 * ia, 2 * ia + 1, 2 * ib, 2 * ib + 1]
    # stabilizer elements whose (a, b)-part lies in the Bell group span{(1,0,1,0), (0,1,0,2)}
    bell = np.array([[1, 0, 1, 0], [0, 1, 0, 2]], np.int64)
    # v = c L ; condition: v[cols_ab] in rowspace(bell)  <=>  v[cols_ab] . bell_perp^T = 0 with bell_perp spanning
    # the orthogonal complement (standard dot product) of rowspace(bell)
    perp = nullspace3(bell)
    C = (L[:, cols_ab] @ perp.T) % 3                       # (rows of L) x 2
    coeffs = nullspace3(C.T)                               # c with c C = 0
    V = (coeffs @ L) % 3
    keep = [i for i in range(len(legs)) if i not in (ia, ib)]
    cols = [2 * i + t for i in keep for t in range(2)]
    W = rref3(V[:, cols])
    return W, [legs[i] for i in keep]


def build_network(tensors, bonds, tensor_L):
    """tensors: dict name -> list of leg labels (length 6); bonds: list of (leg, leg) pairs.  Returns the Lagrangian on
    the remaining (dangling) legs."""
    legs, blocks = [], []
    for name, tl in tensors.items():
        legs += tl
        blocks.append(tensor_L)
    n = len(legs)
    L = np.zeros((n, 2 * n), np.int64)
    r = c = 0
    for B in blocks:
        k = B.shape[0]
        L[r:r + k, 2 * c:2 * c + 2 * k] = B
        r += k
        c += k
    for a, b in bonds:
        L, legs = contract(L, legs, a, b)
    return L, legs


def fix_bulk(L, legs, bulk):
    """put every bulk leg into |0> (add Z on it) and drop it: boundary state stays pure"""
    for bl in bulk:
        i = legs.index(bl)
        # stabilizer elements with x-part 0 on leg i (commuting with Z_i); keep and restrict
        col = 2 * i
        coeffs = nullspace3(L[:, [col]].T)
        V = (coeffs @ L) % 3
        keep = [j for j in range(len(legs)) if j != i]
        cols = [2 * j + t for j in keep for t in range(2)]
        L = rref3(V[:, cols])
        legs = [legs[j] for j in keep]
    return L, legs


def min_cut(net_edges, boundary_legs, A, leg_owner):
    """minimal number of bonds separating boundary legs A from the other boundary legs"""
    Gf = nx.Graph()
    for u, v in net_edges:
        if Gf.has_edge(u, v):
            Gf[u][v]["capacity"] += 1
        else:
            Gf.add_edge(u, v, capacity=1)
    src, snk = "SRC", "SNK"
    for leg in boundary_legs:
        t = leg_owner[leg]
        target = src if leg in A else snk
        # each boundary leg is a bond of capacity 1 between its tensor and src/snk
        if Gf.has_edge(t, target):
            Gf[t][target]["capacity"] += 1
        else:
            Gf.add_edge(t, target, capacity=1)
    value, _ = nx.minimum_cut(Gf, src, snk)
    return int(value)


# ---------------------------------------------------------------- tilings
def pentagon_patch(layers=1):
    """HaPPY-like {5,4} patch: a central pentagon tensor, a ring of 5 edge-neighbours and 5 corner tensors.  Every tensor
    has 6 legs: one bulk leg and five edge legs; edge legs shared between tensors become bonds, the rest boundary."""
    T = {}
    edges_by_tensor = {}
    # central tensor C with edges e0..e4; neighbour N_k shares edge e_k with C; corner K_k sits between N_k and N_{k+1}
    T["C"] = ["C_bulk"] + [f"C_e{k}" for k in range(5)]
    for k in range(5):
        T[f"N{k}"] = [f"N{k}_bulk"] + [f"N{k}_e{j}" for j in range(5)]
        T[f"K{k}"] = [f"K{k}_bulk"] + [f"K{k}_e{j}" for j in range(5)]
    bonds, net = [], []
    for k in range(5):
        bonds.append((f"C_e{k}", f"N{k}_e0"))
        net.append(("C", f"N{k}"))
        # N_k edge 1 meets corner K_k; N_{k+1} edge 4 meets corner K_k  ({5,4}: four pentagons around a vertex)
        bonds.append((f"N{k}_e1", f"K{k}_e0"))
        net.append((f"N{k}", f"K{k}"))
        bonds.append((f"N{(k + 1) % 5}_e4", f"K{k}_e1"))
        net.append((f"N{(k + 1) % 5}", f"K{k}"))
    return T, bonds, net


def square_patch():
    """flat {4,4} patch with the same tensors: a 3x3 grid of 6-leg tensors (4 edge legs used in-plane, plus bulk and one
    extra boundary leg each)"""
    T, bonds, net = {}, [], []
    for i in range(3):
        for j in range(3):
            T[f"S{i}{j}"] = [f"S{i}{j}_bulk"] + [f"S{i}{j}_{d}" for d in ("n", "e", "s", "w", "x")]
    for i in range(3):
        for j in range(3):
            if j < 2:
                bonds.append((f"S{i}{j}_e", f"S{i}{j + 1}_w"))
                net.append((f"S{i}{j}", f"S{i}{j + 1}"))
            if i < 2:
                bonds.append((f"S{i}{j}_s", f"S{i + 1}{j}_n"))
                net.append((f"S{i}{j}", f"S{i + 1}{j}"))
    return T, bonds, net


def boundary_order(T, bonds, boundary):
    """order boundary legs around the patch (by tensor angle, then leg index) for contiguous intervals"""
    return boundary


def rt_check(T, bonds, net, tensor_L, order_fn):
    L, legs = build_network(T, bonds, tensor_L)
    bulk = [l for l in legs if l.endswith("_bulk")]
    L, legs = fix_bulk(L, legs, bulk)
    owner = {leg: name for name, tl in T.items() for leg in tl}
    boundary = order_fn(legs)
    n = len(boundary)
    assert rank3(L) == n, "boundary state must be pure"
    rows = []
    for length in range(1, n // 2 + 1):
        for start in range(n):
            A = [boundary[(start + t) % n] for t in range(length)]
            idx = [legs.index(a) for a in A]
            S = entropy(L, idx)
            cut = min_cut(net, boundary, set(A), owner)
            rows.append((length, start, S, cut))
    agree = sum(1 for r in rows if r[2] == r[3])
    below = sum(1 for r in rows if r[2] < r[3])
    return dict(boundary_legs=n, intervals=len(rows), rt_exact=agree, entropy_below_cut=below,
                entropy_above_cut=len(rows) - agree - below,
                by_length={str(l): [sum(1 for r in rows if r[0] == l and r[2] == r[3]),
                                    sum(1 for r in rows if r[0] == l)] for l in range(1, n // 2 + 1)})


def pentagon_order(legs):
    def key(leg):
        name, rest = leg.split("_", 1)
        if name.startswith("N"):
            k = int(name[1:])
            return (2 * k, int(rest[1:]))
        if name.startswith("K"):
            k = int(name[1:])
            return (2 * k + 1, int(rest[1:]))
        return (99, 0)
    return sorted(legs, key=key)


def square_order(legs):
    # walk the perimeter of the 3x3 grid clockwise; at each tensor list its free legs
    perim = [(0, 0), (0, 1), (0, 2), (1, 2), (2, 2), (2, 1), (2, 0), (1, 0), (1, 1)]
    out = []
    for i, j in perim:
        out += sorted(l for l in legs if l.startswith(f"S{i}{j}_"))
    return out


def choi_lagrangian(p):
    rows = []
    for k in range(4):
        e = np.zeros(4, np.int64)
        e[k] = 1
        inp = e.copy()
        inp[1::2] = (-inp[1::2]) % 3
        rows.append(np.concatenate([inp, (p @ e) % 3]))
    return np.array(rows, np.int64)


def perfect_tick_choi():
    """Lagrangian of the Choi state of the perfect tick p on legs (in1, in2, out1, out2): rows (e_k^*, p e_k)"""
    sys.path.insert(0, str(ROOT / "analysis"))
    import w33_pass11193_optimal_perfect_gate_compiler as C
    p = C.fixed_perfect_gate() % 3
    rows = []
    for k in range(4):
        e = np.zeros(4, np.int64)
        e[k] = 1
        inp = e.copy()
        inp[1::2] = (-inp[1::2]) % 3
        rows.append(np.concatenate([inp, (p @ e) % 3]))
    return np.array(rows, np.int64), p


def square_hyperbolic_patch():
    """{4,5} patch of 4-leg tensors (no bulk): central square C, edge-neighbours N_k, and at each corner of C the two
    further squares K_k, K'_k that complete the five squares around that vertex"""
    T, bonds, net = {}, [], []
    T["C"] = [f"C_e{k}" for k in range(4)]
    for k in range(4):
        T[f"N{k}"] = [f"N{k}_e{j}" for j in range(4)]
        T[f"K{k}"] = [f"K{k}_e{j}" for j in range(4)]
        T[f"Q{k}"] = [f"Q{k}_e{j}" for j in range(4)]
    for k in range(4):
        bonds.append((f"C_e{k}", f"N{k}_e0"))
        net.append(("C", f"N{k}"))
        bonds.append((f"N{k}_e1", f"K{k}_e0"))
        net.append((f"N{k}", f"K{k}"))
        bonds.append((f"N{(k + 1) % 4}_e3", f"Q{k}_e0"))
        net.append((f"N{(k + 1) % 4}", f"Q{k}"))
        bonds.append((f"K{k}_e1", f"Q{k}_e3"))
        net.append((f"K{k}", f"Q{k}"))
    return T, bonds, net


def hyperbolic_square_order(legs):
    def key(leg):
        name, rest = leg.split("_", 1)
        k = int(name[1:]) if name[0] in "NKQ" else 0
        rank = {"N": 0, "K": 1, "Q": 2}.get(name[0], 9)
        return (k, rank, int(rest[1:]))
    return sorted(legs, key=key)


def flat_grid(nr, nc):
    """flat {4,4} grid of 4-leg tensors with no boundary legs in the interior: legs (W, N, E, S)"""
    T, bonds, net = {}, [], []
    for i in range(nr):
        for j in range(nc):
            T[f"G{i}_{j}"] = [f"G{i}_{j}_W", f"G{i}_{j}_N", f"G{i}_{j}_E", f"G{i}_{j}_S"]
    for i in range(nr):
        for j in range(nc):
            if j < nc - 1:
                bonds.append((f"G{i}_{j}_E", f"G{i}_{j + 1}_W"))
                net.append((f"G{i}_{j}", f"G{i}_{j + 1}"))
            if i < nr - 1:
                bonds.append((f"G{i}_{j}_S", f"G{i + 1}_{j}_N"))
                net.append((f"G{i}_{j}", f"G{i + 1}_{j}"))
    return T, bonds, net


def grid_order(nr, nc):
    def order(legs):
        out = ["G0_0_W"] + [f"G0_{j}_N" for j in range(nc)] + [f"G{i}_{nc - 1}_E" for i in range(nr)] + \
              [f"G{nr - 1}_{j}_S" for j in reversed(range(nc))] + [f"G{i}_0_W" for i in reversed(range(1, nr))]
        assert sorted(out) == sorted(legs)
        return out
    return order


def run():
    G6 = find_ame6()
    L6 = graph_lagrangian(G6)
    res = dict(pass_id=11226, ame6_graph=G6.tolist(), ame6_perfect=is_perfect(L6, 6))
    T, bonds, net = pentagon_patch()
    res["pentagon_{5,4}"] = rt_check(T, bonds, net, L6, pentagon_order)
    T, bonds, net = square_patch()
    res["square_{4,4}"] = rt_check(T, bonds, net, L6, square_order)
    Lp, p = perfect_tick_choi()
    res["perfect_tick_choi_is_AME43"] = is_perfect(Lp, 4)
    T, bonds, net = square_hyperbolic_patch()
    res["perfect_tick_{4,5}"] = rt_check(T, bonds, net, Lp, hyperbolic_square_order)
    SUM = np.array([[1, 0, 0, 0], [0, 1, 0, 2], [1, 0, 1, 0], [0, 0, 0, 1]], np.int64).T
    Ls = choi_lagrangian(SUM)
    res["control_SUM_choi_is_AME43"] = is_perfect(Ls, 4)
    res["control_SUM_{4,5}"] = rt_check(T, bonds, net, Ls, hyperbolic_square_order)
    res["perfect_tick_flat_grids"] = {}
    for nr in (2, 3, 4, 5):
        T, bonds, net = flat_grid(nr, nr)
        res["perfect_tick_flat_grids"][f"{nr}x{nr}"] = rt_check(T, bonds, net, Lp, grid_order(nr, nr))
    res["flat_grids_rt_exact"] = all(v["rt_exact"] == v["intervals"] for v in res["perfect_tick_flat_grids"].values())
    res["rt_exact_everywhere"] = all(res[k]["rt_exact"] == res[k]["intervals"]
                                     for k in ("pentagon_{5,4}", "square_{4,4}", "perfect_tick_{4,5}"))
    return res


def main():
    res = run()
    OUT.write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
