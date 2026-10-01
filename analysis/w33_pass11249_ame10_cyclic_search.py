#!/usr/bin/env python3
"""Pass 11249: an algebraic search for AME(10,3) stabilizer states -- every code with a twisted cyclic symmetry.

Passes 11224/11229 found only F9-linear AME(10,3) states (Glynn's arc) by random search.  Here the search is exhaustive
inside a structured family: all ten-qutrit stabilizer states (Lagrangian subspaces C of F_3^20 with the symplectic form)
invariant under a twisted cyclic shift
        sigma = (cyclic shift of the 10 qutrits) with a local twist T in SL(2,3) on the wrap-around,
for the semisimple twists T in {I, -I, J} (J^2 = -1; the other classes of SL(2,3) have order divisible by 3).  Every
automorphism of a state that permutes the qutrits as a 10-cycle is LC-conjugate to such a sigma (move all local parts
onto one coordinate), and Glynn's state has 10-cycles in its permutation group PGL(2,9) -- a positive control.

Method (exact over F_3): sigma is semisimple, F_3^20 = sum of V_f = ker f(sigma) over irreducible f | x^N - 1; V_f is a
vector space over K_f = F_3[x]/f; omega pairs V_f with V_{f*} (f* reciprocal).  An invariant Lagrangian is a direct sum
of K-subspaces: for f = f*, a Lagrangian K-subspace of V_f; for f != f*, a K-subspace of V_f together with its
annihilator in V_{f*}.  All are enumerated; each is tested for minimum distance 6 (AME) and for F9-linearity: the local
endomorphism algebra {A = diag(A_1..A_10) : A C <= C} is computed by linear algebra and searched for an element with
A_i^2 = -1 for all i.
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "w33_pass11249_ame10_cyclic_search.json"
p = 3
n = 10
OMEGA = np.zeros((2 * n, 2 * n), np.int64)
for i in range(n):
    OMEGA[2 * i, 2 * i + 1], OMEGA[2 * i + 1, 2 * i] = 1, -1


# ----------------------------------------------------------------------------- mod-3 linear algebra
def rref(M):
    M = np.array(M, np.int64) % p
    r, piv = 0, []
    for c in range(M.shape[1]):
        nz = [i for i in range(r, M.shape[0]) if M[i, c]]
        if not nz:
            continue
        M[[r, nz[0]]] = M[[nz[0], r]]
        M[r] = (M[r] * pow(int(M[r, c]), -1, p)) % p
        for i in range(M.shape[0]):
            if i != r and M[i, c]:
                M[i] = (M[i] - M[i, c] * M[r]) % p
        piv.append(c)
        r += 1
        if r == M.shape[0]:
            break
    return M[:r], piv


def nullspace(M):
    R, piv = rref(M)
    free = [c for c in range(M.shape[1]) if c not in piv]
    out = []
    for f in free:
        v = np.zeros(M.shape[1], np.int64)
        v[f] = 1
        for i, c in enumerate(piv):
            v[c] = (-R[i, f]) % p
        out.append(v)
    return np.array(out, np.int64).reshape(len(out), M.shape[1])


def key(rows):
    R, _ = rref(rows)
    return tuple(R.ravel())


# ----------------------------------------------------------------------------- the twisted shift
TWISTS = {"I": np.eye(2, dtype=np.int64), "-I": (2 * np.eye(2, dtype=np.int64)) % 3,
          "J": np.array([[0, 2], [1, 0]], np.int64)}


def sigma(T):
    S = np.zeros((2 * n, 2 * n), np.int64)
    for i in range(n):
        j = (i + 1) % n
        blk = T if j == 0 else np.eye(2, dtype=np.int64)
        S[2 * j:2 * j + 2, 2 * i:2 * i + 2] = blk
    assert not ((S.T @ OMEGA @ S - OMEGA) % p).any()
    return S


def order(S):
    M, k = S.copy(), 1
    while not np.array_equal(M % p, np.eye(2 * n, dtype=np.int64)):
        M = (M @ S) % p
        k += 1
    return k


def poly_at(coeffs, S):
    """coeffs high->low"""
    M = np.zeros_like(S)
    for c in coeffs:
        M = (M @ S + c * np.eye(len(S), dtype=np.int64)) % p
    return M


def components(S):
    N = order(S)
    x = sp.symbols("x")
    facs = sp.factor_list(sp.Poly(x ** N - 1, x, modulus=p))[1]
    comps = []
    for f, e in facs:
        assert e == 1
        co = [int(c) % p for c in f.all_coeffs()]
        V = nullspace(poly_at(co, S))
        if len(V):
            rec = [int(c) % p for c in reversed(co)]
            lead = rec[0]
            rec = [(c * pow(lead, -1, p)) % p for c in rec]
            comps.append(dict(f=tuple(co), fstar=tuple(rec), deg=len(co) - 1, V=V))
    assert sum(len(c["V"]) for c in comps) == 2 * n
    return N, comps


def k_subspaces(V, S, d):
    """all K-subspaces of V (K = F_3[sigma]|V of degree d): returns list of basis arrays"""
    kd = len(V) // d
    subs = [np.zeros((0, 2 * n), np.int64), V]
    if kd == 2:
        seen = set()
        for c in itertools.product(range(p), repeat=len(V)):
            if not any(c):
                continue
            v = (np.array(c) @ V) % p
            L = [v]
            for _ in range(d - 1):
                L.append((S @ L[-1]) % p)
            L = np.array(L)
            R, _ = rref(L)
            k = tuple(R.ravel())
            if k not in seen:
                seen.add(k)
                subs.append(R)
    elif kd > 2:
        raise ValueError("K-dimension > 2 not handled")
    return subs


def isotropic(B):
    return not ((B @ OMEGA @ B.T) % p).any()


def annihilator(B, V):
    """vectors of span V orthogonal (omega) to all rows of B"""
    if len(B) == 0:
        return V
    coef = nullspace((B @ OMEGA @ V.T) % p)
    return (coef @ V) % p if len(coef) else np.zeros((0, 2 * n), np.int64)


# ----------------------------------------------------------------------------- tests
ALL = np.array(list(itertools.product(range(p), repeat=n)), np.int64)


def min_weight(B):
    W = (ALL @ B) % p
    w = ((W[:, 0::2] != 0) | (W[:, 1::2] != 0)).sum(1)
    return int(w[1:].min()) if len(w) > 1 else 0


def f9_linear(B):
    """the local algebra {diag(A_i): A C <= C}; is there an element with A_i^2 = -1 for all i?"""
    rows = []
    # unknowns: a_i = (A_i)_{rs}, 4 per qutrit; condition omega(b', A b) = 0 for all basis b, b'
    for b in B:
        for b2 in B:
            row = np.zeros(4 * n, np.int64)
            w = (b2 @ OMEGA) % p
            for i in range(n):
                for r in range(2):
                    for s in range(2):
                        row[4 * i + 2 * r + s] = (w[2 * i + r] * b[2 * i + s]) % p
            rows.append(row)
    Nsp = nullspace(np.array(rows))
    dim = len(Nsp)
    if dim > 8:
        return dim, None
    for c in itertools.product(range(p), repeat=dim):
        a = (np.array(c) @ Nsp) % p if dim else np.zeros(4 * n, np.int64)
        ok = True
        for i in range(n):
            A = a[4 * i:4 * i + 4].reshape(2, 2)
            if not np.array_equal((A @ A) % p, (2 * np.eye(2, dtype=np.int64)) % p):
                ok = False
                break
        if ok:
            return dim, True
    return dim, False


def search(tname):
    S = sigma(TWISTS[tname])
    N, comps = components(S)
    rec = dict(twist=tname, order=N, components=[dict(f=c["f"], deg=c["deg"], kdim=len(c["V"]) // c["deg"],
                                                       self_reciprocal=c["f"] == c["fstar"]) for c in comps])
    by = {c["f"]: c for c in comps}
    choice_lists, done = [], set()
    for c in comps:
        if c["f"] in done:
            continue
        if c["f"] == c["fstar"]:
            opts = [B for B in k_subspaces(c["V"], S, c["deg"]) if 2 * len(B) == len(c["V"]) and isotropic(B)]
            done.add(c["f"])
        else:
            d = by[c["fstar"]]
            opts = []
            for B in k_subspaces(c["V"], S, c["deg"]):
                A = annihilator(B, d["V"])
                opts.append(np.vstack([B, A]) if len(A) else B)
            done |= {c["f"], c["fstar"]}
        choice_lists.append(opts)
    rec["choices_per_block"] = [len(o) for o in choice_lists]
    codes = ame = f9 = 0
    found, keys = [], set()
    for pick in itertools.product(*choice_lists):
        B = np.vstack([b for b in pick if len(b)])
        if len(B) != n:
            continue
        assert isotropic(B)
        codes += 1
        if min_weight(B) < 6:
            continue
        k = key(B)
        if k in keys:
            continue
        keys.add(k)
        ame += 1
        dim, lin = f9_linear(B)
        f9 += bool(lin)
        found.append(dict(local_algebra_dim=dim, f9_linear=lin, basis=B.tolist()))
    rec.update(lagrangians=codes, ame_states=ame, f9_linear=f9, non_f9=[x for x in found if not x["f9_linear"]],
               examples=found[:2])
    return rec


def run():
    return dict(pass_id=11249, searches=[search(t) for t in ("I", "-I", "J")])


def main():
    res = run()
    OUT.write_text(json.dumps(res, indent=1, default=str))
    for s in res["searches"]:
        print(s["twist"], s["order"], s["choices_per_block"], "lagrangians", s["lagrangians"], "AME", s["ame_states"],
              "F9-linear", s["f9_linear"], "non-F9", len(s["non_f9"]))


if __name__ == "__main__":
    main()
