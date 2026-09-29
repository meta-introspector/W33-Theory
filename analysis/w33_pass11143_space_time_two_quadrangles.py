#!/usr/bin/env python3
"""Pass 11143: space and time are two generalized quadrangles on the same 40 Pauli classes -- they share exactly the 16
product lines, and partial transposition exchanges their 24 entangled lines (two-qutrit Bell frames <-> one-qutrit
Clifford channels).

Setting.  P_(x,z) = X^x Z^z on one qutrit; a two-qutrit (or two-time) Pauli label is (a, b) in F_3^2 + F_3^2, and the
40 projective points of F_3^4 are the 40 points of W(3,3).  Two symplectic forms live on the same space:
    spatial   w_s((a,b),(a',b')) = [a,a'] + [b,b']   (commutation of P_a (x) P_b: Section 2 of the paper),
    temporal  w_t((a,b),(a',b')) = [a,a'] - [b,b']   (left/right action on one qutrit's operators: Theorem 3.1).
Each makes the 40 points a generalized quadrangle GQ(3,3); its lines are the 40 Lagrangian planes of that form.

Results (all exact, brute force over F_3):
  * both structures are GQ(3,3) with collinearity SRG(40,12,2,4); they share EXACTLY 16 lines -- the product planes
    l1 + l2 (4 x 4) -- and differ in 24 each;
  * the 24 spatial-only lines are the graphs {(a, N a)} of the 24 maps with det N = -1 (anti-symplectic); their
    stabilizer states are the 216 maximally entangled two-qutrit stabilizer states (Schmidt coefficients 1/3,1/3,1/3),
    the 16 common lines carry the 144 product states (360 = 9 x 40, cf. BT821);
  * the 24 temporal-only lines are the graphs of SL(2,3) = the one-qutrit Clifford group modulo Paulis (generated here
    from the Fourier and phase gates: 216 unitaries mod phase, 24 symplectic classes); the two-time pseudo-density
    operator R_U = J(U)/3 of the channel rho -> U rho U^dagger (maximally mixed input; J = Jamiolkowski matrix = partial
    transpose of the Choi state) is supported exactly on the temporal line graph(-M_U), has spectrum (+1/3)^6 (-1/3)^3
    and temporal negativity 1 -- equal to the spatial negativity of a maximally entangled pair;
  * partial transposition (x,z) -> (x,-z) on the second leg maps the 24 spatial lines bijectively onto the 24
    temporal lines; on states, rho^(T_B) of every maximally entangled stabilizer state IS, operator for operator, the
    two-time pseudo-density operator R_U of a one-qutrit Clifford gate -- a BIJECTION between the 216 maximally
    entangled two-qutrit stabilizer states and the 216 Clifford unitaries (mod phase);
  * the collineations preserving BOTH structures (exhaustive search over all 103680 similitudes of the spatial form)
    are exactly the block-diagonal maps with equal multipliers and the factor swap: order 2304 (1152 projectively) --
    local operations and exchange of the parties.  Every entangling Clifford operation distinguishes space from time.
Reading: a maximally entangled pair of qutrits and ONE qutrit persisting through a Clifford gate are the same point set
seen through the two forms; the entangled lines are exactly where space and time disagree.  This is the W(3,3) form of
the "partial transpose as a space-time swap" (Fullwood-Li, arXiv:2508.12256; Marcovitch-Reznik, arXiv:1107.2186) --
those works are general; the two-quadrangle statement, the 16/24 split and the common symmetry group are specific to
W(3,3).  Nothing here moves information backwards in time: temporal correlations may signal, spatial ones may not.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "w33_pass11143_space_time_two_quadrangles.json"
W = np.exp(2j * np.pi / 3)
X = np.roll(np.eye(3), 1, axis=0)
Z = np.diag([1, W, W * W])


def pauli(u):
    return np.linalg.matrix_power(X, u[0] % 3) @ np.linalg.matrix_power(Z, u[1] % 3)


def br(a, b):
    return (a[0] * b[1] - a[1] * b[0]) % 3


def form(sign):
    return lambda p, q: (br(p[:2], q[:2]) + sign * br(p[2:], q[2:])) % 3


VECS = [v for v in itertools.product(range(3), repeat=4) if any(v)]


def normal(v):
    v = np.array(v) % 3
    i = next(k for k in range(4) if v[k])
    return tuple((v * pow(int(v[i]), -1, 3)) % 3)


POINTS = sorted({normal(v) for v in VECS})


def lagrangians(f):
    lines = set()
    for p, q in itertools.combinations(POINTS, 2):
        if f(p, q) == 0:
            span = {normal((a * np.array(p) + b * np.array(q)) % 3) for a in range(3) for b in range(3) if a or b}
            lines.add(frozenset(span))
    return lines


def srg(f):
    A = np.array([[1 if p != q and f(p, q) == 0 else 0 for q in POINTS] for p in POINTS])
    k = A.sum(1)[0]
    lam = {int(A[i] @ A[j]) for i in range(40) for j in range(40) if i != j and A[i, j]}
    mu = {int(A[i] @ A[j]) for i in range(40) for j in range(40) if i != j and not A[i, j]}
    return int(k), sorted(lam), sorted(mu)


def graph_line(N):
    return frozenset(normal((a[0], a[1], *(np.array(N) @ np.array(a) % 3))) for a in itertools.product(range(3), repeat=2) if any(a))


def mats2(det):
    return [np.array(m).reshape(2, 2) for m in itertools.product(range(3), repeat=4)
            if (m[0] * m[3] - m[1] * m[2]) % 3 == det]


def stab_states(line):
    ops = [np.kron(pauli(p[:2]), pauli(p[2:])) for p in line]
    rng = np.random.default_rng(0)
    H = sum(rng.normal() * (o + o.conj().T) + 1j * rng.normal() * (o - o.conj().T) for o in ops)
    _, V = np.linalg.eigh(H)
    return [V[:, i] for i in range(9)]


def schmidt(psi):
    return np.round(np.linalg.svd(psi.reshape(3, 3), compute_uv=False) ** 2, 9)


def cliffords():
    F = np.array([[W ** (j * k) for k in range(3)] for j in range(3)]) / np.sqrt(3)
    S = np.diag([1, 1, W])
    def key(U):
        i = np.argmax(np.abs(U.ravel()) > 1e-9)
        return tuple(np.round((U / (U.ravel()[i] / abs(U.ravel()[i]))).ravel(), 6))
    seen, frontier = {key(np.eye(3)): np.eye(3)}, [np.eye(3)]
    while frontier:
        nxt = []
        for U in frontier:
            for G in (F, S, X, Z):
                V = G @ U
                k = key(V)
                if k not in seen:
                    seen[k] = V; nxt.append(V)
        frontier = nxt
    return list(seen.values())


def symplectic_class(U):
    cols = []
    for e in ((1, 0), (0, 1)):
        Pe = U @ pauli(e) @ U.conj().T
        hit = [u for u in itertools.product(range(3), repeat=2) if abs(abs(np.trace(pauli(u).conj().T @ Pe)) - 3) < 1e-6]
        cols.append(hit[0])
    return np.array(cols).T % 3


def jamiolkowski(U):
    J = np.zeros((9, 9), complex)
    for i in range(3):
        for j in range(3):
            Eji = np.zeros((3, 3)); Eji[j, i] = 1
            Eij = np.zeros((3, 3)); Eij[i, j] = 1
            J += np.kron(Eij, U @ Eji @ U.conj().T)
    return J


def weyl_support(R):
    return frozenset(normal((*a, *b)) for a in itertools.product(range(3), repeat=2) for b in itertools.product(range(3), repeat=2)
                     if (any(a) or any(b)) and abs(np.trace(np.kron(pauli(a), pauli(b)).conj().T @ R)) > 1e-9)


def partial_transpose(rho):
    return rho.reshape(3, 3, 3, 3).transpose(0, 3, 2, 1).reshape(9, 9)


def gsp_spatial():
    Js = np.zeros((4, 4), int)
    Js[0, 1], Js[1, 0], Js[2, 3], Js[3, 2] = 1, -1, 1, -1
    gens = [np.array(g) % 3 for g in (
        [[1, 1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], [[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]],
        [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 1], [0, 0, 0, 1]], [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]],
        [[1, 0, 0, 0], [0, 1, 1, 0], [0, 0, 1, 0], [1, 0, 0, 1]], [[1, 0, 0, 0], [0, 2, 0, 0], [0, 0, 1, 0], [0, 0, 0, 2]])]
    seen = {tuple(np.eye(4, dtype=int).ravel())}
    frontier = [np.eye(4, dtype=int)]
    while frontier:
        nxt = []
        for g in frontier:
            for h in gens:
                m = (h @ g) % 3
                k = tuple(m.ravel())
                if k not in seen:
                    seen.add(k); nxt.append(m)
        frontier = nxt
    return [np.array(k).reshape(4, 4) for k in seen], Js


def summarize():
    ws, wt = form(1), form(-1)
    Ls, Lt = lagrangians(ws), lagrangians(wt)
    common = Ls & Lt
    prod = {frozenset(normal((*(a * np.array(l1)), *(b * np.array(l2)))) for a in range(3) for b in range(3) if a or b)
            for l1 in [(1, 0), (0, 1), (1, 1), (1, 2)] for l2 in [(1, 0), (0, 1), (1, 1), (1, 2)]}
    sp_graphs = {graph_line(N) for N in mats2(2)}
    tm_graphs = {graph_line(N) for N in mats2(1)}
    # stabilizer states
    sch = {'entangled': set(), 'product': set()}
    for line in Ls:
        kind = 'product' if line in common else 'entangled'
        for psi in stab_states(line):
            sch[kind].add(tuple(schmidt(psi)))
    # Clifford channels
    Us = cliffords()
    classes = {}
    for U in Us:
        classes.setdefault(tuple(symplectic_class(U).ravel()), U)
    pdm = []
    for k, U in classes.items():
        R = jamiolkowski(U) / 3
        ev = np.round(np.linalg.eigvalsh(R), 9)
        M = np.array(k).reshape(2, 2)
        pdm.append(dict(support_is_graph_minus_M=weyl_support(R) == graph_line((-M) % 3),
                        spectrum=sorted(set(ev.tolist())), neg=float((np.abs(ev).sum() - 1) / 2)))
    # partial transpose of entangled stabilizer states -> PDM of a Clifford channel
    by_support = {}
    for U in Us:
        by_support.setdefault(weyl_support(jamiolkowski(U) / 3), []).append(jamiolkowski(U) / 3)
    pt_ok, pt_lines = 0, set()
    for line in Ls - common:
        for psi in stab_states(line):
            rho = np.outer(psi, psi.conj()); ptr = partial_transpose(rho)
            sup = weyl_support(ptr); pt_lines.add(sup)
            pt_ok += any(np.allclose(ptr, R, atol=1e-9) for R in by_support.get(sup, []))
    Rall = [jamiolkowski(U) / 3 for U in Us]
    matched = set()
    for line in Ls - common:
        for psi in stab_states(line):
            ptr = partial_transpose(np.outer(psi, psi.conj()))
            m = [i for i, R in enumerate(Rall) if np.allclose(ptr, R, atol=1e-9)]
            if len(m) == 1: matched.add(m[0])
    # common automorphisms: similitudes of the spatial form that are also similitudes of the temporal form
    Sp, Js = gsp_spatial()
    Jt = Js.copy(); Jt[2, 3], Jt[3, 2] = -1, 1
    both = 0
    for g in Sp:
        a = (g.T @ Jt @ g) % 3
        if np.array_equal(a, (1 * Jt) % 3) or np.array_equal(a, (2 * Jt) % 3):
            both += 1
    res = dict(pass_id=11143, n_points=len(POINTS), spatial_lines=len(Ls), temporal_lines=len(Lt),
               spatial_srg=srg(ws), temporal_srg=srg(wt), common_lines=len(common), common_are_products=common == prod,
               spatial_only_are_detminus1_graphs=(Ls - common) == sp_graphs,
               temporal_only_are_SL23_graphs=(Lt - common) == tm_graphs,
               schmidt_entangled=sorted(sch['entangled']), schmidt_product=sorted(sch['product']),
               clifford_unitaries_mod_phase=len(Us), clifford_symplectic_classes=len(classes),
               pdm_all_supported_on_temporal_lines=all(p['support_is_graph_minus_M'] for p in pdm),
               pdm_spectra=sorted({tuple(p['spectrum']) for p in pdm}), pdm_negativity=sorted({round(p['neg'], 9) for p in pdm}),
               partial_transpose_lines_are_temporal=pt_lines == (Lt - common), pt_states_equal_to_a_clifford_pdm=pt_ok, bijection_states_to_cliffords=len(matched) == 216 == len(Us),
               gsp_spatial_order=len(Sp), common_similitudes=both)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True, default=str))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1, default=str))
