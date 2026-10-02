"""Pass 11310: the U81 cocycle orientation is NOT the orientation of the cubic gate (an identification ruled out).

Codex (Passes 11067-11069; w33_20261001_u81_antiunitary_orientation.py) put the chamber group on H27 x F3 with the law
(h,d)*_s(h',D) = (hh', d+D+s kappa(h,h')), kappa = A^2 b - 2Ac; s = +-1 give U81 (centre 3, derived 9, class 3, 36
elements of order 9), qutrit complex conjugation exchanges s = +1 and -1, and 'an energetic coupling selecting one
sign remains open'.  The corpus already identifies U81 as the chamber Sylow-3 subgroup of Sp(4,3) (Passes 11061-11075).

Tempting identification tested: complex conjugation also exchanges the cubic gate T with T^-1 = T^*, so is the abstract
cocycle orientation the choice of cube root of Z made by the magic gate?  The one-qutrit 'magic Sylow' group
<X, S, T> modulo scalars is computed exactly (elements X^b diag(zeta^{f(x)}), f mod constants):
  * order 81, centre 3 (= <Z>), derived subgroup 9, class 3 -- the same coarse invariants as U81 --
  * but 18 elements of order 9, against 36 for U81 (counted here from Codex's law for s = 1, 2).
So the two groups are NOT isomorphic, and no extension of e0 -> g0, e1 -> X is a homomorphism for either s (searched).
GAP (script in the certificate): <X,S,T>/scalars = SmallGroup(81,9); U81 (both s) = SmallGroup(81,7) = C3 wr C3 =
Sylow-3 of Sp(4,3) and of PSp(4,3).  The cocycle orientation therefore lives at the Clifford (symplectic) level of the
substrate, not in the magic layer; selecting it is not the same as choosing T over T^dagger.
"""

from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11068_ternary_cocycle_deformation as D68  # noqa: E402

OUT = ROOT / "data" / "w33_pass11310_magic_orientation_cocycle.json"


def norm(h):
    return tuple((v - h[0]) % 9 for v in h)


def mul(p, q):
    """(b, f): |x> -> zeta^{f[x]} |x + b>, modulo scalars"""
    b, f = p
    B, F = q
    return ((b + B) % 3, norm([(F[x] + f[(x + B) % 3]) % 9 for x in range(3)]))


ID = (0, (0, 0, 0))
X = (1, (0, 0, 0))
S = (0, norm([0, 0, 3]))


def T_gate(sign):
    return (0, norm([0, 1 * sign % 9, 8 * sign % 9]))


def generate(gens):
    G = {ID}
    fr = [ID]
    while fr:
        nf = []
        for M in fr:
            for g in gens:
                N = mul(g, M)
                if N not in G:
                    G.add(N)
                    nf.append(N)
        fr = nf
    return sorted(G)


def invariants(G):
    ix = {g: i for i, g in enumerate(G)}
    n = len(G)
    M = [[ix[mul(a, b)] for b in G] for a in G]
    e = ix[ID]
    inv = [next(j for j in range(n) if M[i][j] == e) for i in range(n)]

    def comm(a, b):
        return M[M[M[a][b]][inv[a]]][inv[b]]

    def clos(S0):
        S0 = set(S0) | {e}
        while True:
            new = {M[a][b] for a in S0 for b in S0}
            if new <= S0:
                return S0
            S0 |= new

    cent = [i for i in range(n) if all(M[i][j] == M[j][i] for j in range(n))]
    D = clos({comm(a, b) for a in range(n) for b in range(n)})
    G3 = clos({comm(a, b) for a in range(n) for b in D})

    def order(a):
        k, x = 1, a
        while x != e:
            x = M[x][a]
            k += 1
        return k
    return dict(order=n, center=len(cent), derived=len(D), gamma3=len(G3),
                order9=sum(order(a) == 9 for a in range(n)), center_elements=[G[i] for i in cent])


def u81_word_map(g0, g1, s):
    """try to extend e0 -> g0, e1 -> g1 to a homomorphism from the s-law onto the magic group; returns the map or None"""
    E = D68.ELTS
    e0, e1 = (1, 0, 0, 0), (0, 1, 0, 0)
    img = {D68.ID if hasattr(D68, "ID") else (0, 0, 0, 0): ID, e0: g0, e1: g1}
    fr = [e0, e1]
    while fr:
        nf = []
        for x in fr:
            for gen, gimg in ((e0, g0), (e1, g1)):
                y = D68.mul(gen, x, s)
                yi = mul(gimg, img[x])
                if y in img:
                    if img[y] != yi:
                        return None
                else:
                    img[y] = yi
                    nf.append(y)
        fr = nf
    if len(img) != len(E) or len(set(img.values())) != len(E):
        return None
    for x, y in itertools.product(E, repeat=2):                     # full homomorphism check
        if img[D68.mul(x, y, s)] != mul(img[x], img[y]):
            return None
    return img


def search(G, s, g1=X):
    """all g0 such that e0 -> g0, e1 -> X is an isomorphism for the s-law"""
    hits = []
    for g0 in G:
        if u81_word_map(g0, g1, s) is not None:
            hits.append(g0)
    return hits


def classify(g, sign):
    """is g in the coset T^k * (Clifford diagonal * Pauli) with k = 1 (T) or 2 (T^2)?  magic content = f mod 3-multiples"""
    b, f = g
    cubic = [(f[x] - 0) % 9 for x in range(3)]
    # the cubic part of f: f(x) = c x^3 + (multiples of 3) -> c = f(1) mod 3 (since x^3 = x mod 3 and f(0) = 0)
    return int(f[1] % 3)


def u81_order9(s):
    E = D68.ELTS

    def order(x):
        k, y = 1, x
        while y != (0, 0, 0, 0):
            y = D68.mul(y, x, s)
            k += 1
        return k
    return sum(order(x) == 9 for x in E)


GAP_RESULTS = {"magic_sylow_X_S_T_mod_scalars": [81, 9], "U81_s1": [81, 7], "U81_s2": [81, 7], "C3_wr_C3": [81, 7],
               "Sylow3_Sp43": [81, 7], "Sylow3_PSp43": [81, 7],
               "order9_counts_SmallGroup_81_7_8_9_10": [36, 54, 18, 72]}


def run():
    res = dict(pass_id=11310)
    out = {}
    for sign, name in ((1, "T"), (-1, "T_inverse")):
        G = generate([X, S, T_gate(sign)])
        inv = invariants(G)
        hits = {s: search(G, s) for s in (1, 2)}
        out[name] = dict(invariants={k: v for k, v in inv.items() if k != "center_elements"},
                         homomorphisms_with_e1_to_X={str(s): len(hits[s]) for s in (1, 2)})
        print(name, out[name], flush=True)
    res["magic_groups"] = out
    res["T_and_T_inverse_generate_same_group_mod_scalars"] = generate([X, S, T_gate(1)]) == generate([X, S, T_gate(-1)])
    res["U81_order9"] = {str(s): u81_order9(s) for s in (1, 2)}
    res["gap"] = GAP_RESULTS
    res["non_isomorphic"] = out["T"]["invariants"]["order9"] != res["U81_order9"]["1"]
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
