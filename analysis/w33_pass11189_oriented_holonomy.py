"""Pass 11189: the arrow that entanglement cannot see -- oriented holonomy on three-qutrit split relations.

Pass 11184 classified how two three-qutrit tensor factorisations F0, F can relate (20 orbitals of Sp(6,3)) by block
ranks, orientations and image/kernel incidences, up to two ties, each a relation and its time reverse (orbital sizes
256 and 6912).  Those invariants cannot see the direction of time: for T = B_F0^-1 B_F the reverse relation has
T^-1 with blocks (T^-1)_ji = adj(T_ij), and every trace of a loop of blocks is unchanged by that transposition.

The way out is that SL(2,3) does not fuse a unipotent u with u^-1, nor a nilpotent N with -N (they are GL(2,3)- but
not SL(2,3)-conjugate), while reversing a loop replaces its holonomy P by adj(P) = P^-1 (or -N).  So the
SL(2,3)-conjugacy classes of the holonomies

    P = T_{r1 c1} adj(T_{r2 c1}) T_{r2 c2} adj(T_{r3 c2}) ... adj(T_{r1 ck})

of closed alternating walks (rows = factors of F0, columns = factors of F) are local-invariant, oriented data.
This pass computes, for every orbital, (i) the multiset of rooted oriented walk holonomy classes (length <= 2k),
(ii) the reverse orbital, (iii) the image under the anti-unitary time reversal tau = (x,z) -> (x,-z), which fixes F0.
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
import w33_pass11182_paper_ticks_mereology as PT  # noqa: E402
import w33_pass11184_merged_orbitals as MO  # noqa: E402

OUT = ROOT / "data" / "w33_pass11189_oriented_holonomy.json"
CACHE = Path.home() / "AppData" / "Local" / "Temp" / "orbit_labels3.npy"
TAU = np.diag([1, 2, 1, 2, 1, 2]).astype(np.int64)
J6 = np.zeros((6, 6), np.int64)
for _k in range(3):
    J6[2 * _k, 2 * _k + 1], J6[2 * _k + 1, 2 * _k] = 1, -1

SL = [np.array([[a, b], [c, d]]) for a in range(3) for b in range(3) for c in range(3) for d in range(3)
      if (a * d - b * c) % 3 == 1]
SLI = [np.array([[g[1, 1], -g[0, 1]], [-g[1, 0], g[0, 0]]]) % 3 for g in SL]


def slclass(P):
    return min(tuple(int(x) for x in ((g @ P @ gi) % 3).ravel()) for g, gi in zip(SL, SLI))


def adj(M):
    return np.array([[M[1, 1], -M[0, 1]], [-M[1, 0], M[0, 0]]]) % 3


def blocks(B):
    return B.reshape(3, 2, 3, 2).transpose(0, 2, 1, 3) % 3


def codes(T):
    det = (T[..., 0, 0] * T[..., 1, 1] - T[..., 0, 1] * T[..., 1, 0]) % 3
    zero = (T == 0).all(axis=(-1, -2))
    return np.where(zero, 0, np.where(det == 0, 1, np.where(det == 1, 2, 3)))


def holonomies(B, kmax=4):
    """multiset over rooted oriented closed alternating walks of length 2k (k <= kmax) of
    (k, block codes along the walk, SL(2,3)-class of the holonomy)"""
    T = blocks(B)
    c = codes(T)
    out = Counter()
    for k in range(2, kmax + 1):
        for rows in itertools.product(range(3), repeat=k):
            for cols in itertools.product(range(3), repeat=k):
                if any(rows[t] == rows[(t + 1) % k] for t in range(k)):
                    continue                        # skip backtracking steps T adj(T) = det * I
                P = np.eye(2, dtype=np.int64)
                cs = []
                for t in range(k):
                    r, cc, rn = rows[t], cols[t], rows[(t + 1) % k]
                    P = P @ T[r, cc] @ adj(T[rn, cc]) % 3
                    cs += [int(c[r, cc]), int(c[rn, cc])]
                out[(k, tuple(cs), slclass(P))] += 1
    return out


J2 = np.array([[0, 1], [-1, 0]])


def twist(B):
    """tau-odd invariant for relations whose rank-1 blocks share images along rows and kernels along columns
    (the 256 pair, where every holonomy vanishes).  With R_i the invertible block of row i, u_i the common image of
    row i and E_{i' sigma(i)} = u_i' (x) phi the rank-1 blocks of column sigma(i), the functional
    psi_i = omega(u_i, R_i .) is proportional to phi: psi_i = rho(i,i') phi.  Q = rho(0,1) rho(1,2) rho(2,0) is
    invariant under all rescalings (each u_i enters twice) and local SL(2,3) changes, equals the product for the
    opposite cyclic order, and changes sign under tau (one omega per factor, three factors).  None if undefined."""
    T = blocks(B)
    c = codes(T)
    inv = [np.flatnonzero(c[i] >= 2) for i in range(3)]
    if any(len(x) != 1 for x in inv) or int((c == 1).sum()) != 6:
        return None
    sig = {i: int(inv[i][0]) for i in range(3)}
    u = {}
    for i in range(3):
        M = T[i, [j for j in range(3) if j != sig[i]][0]]
        u[i] = (M[:, 0] if M[:, 0].any() else M[:, 1]) % 3

    def rho(i, ip):
        E = T[ip, sig[i]]
        k = int(np.flatnonzero(u[ip])[0])
        phi = (E[k] * pow(int(u[ip][k]), -1, 3)) % 3
        if not np.array_equal(np.outer(u[ip], phi) % 3, E % 3):
            return None
        psi = (u[i] @ J2 @ T[i, sig[i]]) % 3
        return next((r for r in (1, 2) if np.array_equal((r * phi) % 3, psi)), None)
    qp = [rho(0, 1), rho(1, 2), rho(2, 0)]
    qm = [rho(0, 2), rho(2, 1), rho(1, 0)]
    if None in qp + qm:
        return None
    a, b = int(np.prod(qp)) % 3, int(np.prod(qm)) % 3
    return a if a == b else None


def canon_multiset(ms):
    """relabel-invariant summary: the multiset itself is already invariant under row/column relabelling because it
    ranges over all rooted walks; return a hashable form"""
    return tuple(sorted(ms.items()))


def run(sample=12, kmax=4):
    Bs = PT.load_facs()
    lab = np.load(CACHE)
    keys = [MO.fac_key(B) for B in Bs]
    index = {k: i for i, k in enumerate(keys)}
    orbitals = sorted(set(lab.tolist()))
    size = Counter(lab.tolist())
    rep = {o: int(np.flatnonzero(lab == o)[0]) for o in orbitals}
    res = []
    hol = {}
    for o in orbitals:
        i = rep[o]
        B = Bs[i]
        Binv = (-J6 @ B.T @ J6) % 3
        rev = int(lab[index[MO.fac_key(Binv)]])
        tau_img = int(lab[index[MO.fac_key((TAU @ B) % 3)]])
        members = np.flatnonzero(lab == o)[:sample]
        ms = {canon_multiset(holonomies(Bs[m], kmax)) for m in members}
        assert len(ms) == 1, "holonomy multiset not constant on the orbital"
        hol[o] = dict(next(iter(ms)))
        res.append(dict(orbital=o, size=size[o], reverse=rev, tau_image=tau_img, self_paired=(rev == o)))
    # separating features for each reversal pair
    pairs = sorted({tuple(sorted((r["orbital"], r["reverse"]))) for r in res if not r["self_paired"]})
    sep = {}
    for a, b in pairs:
        diff = {k: (hol[a].get(k, 0), hol[b].get(k, 0)) for k in set(hol[a]) | set(hol[b])
                if hol[a].get(k, 0) != hol[b].get(k, 0)}
        kmin = min((k[0] for k in diff), default=None)
        sep[f"{a}-{b}"] = dict(size=size[a], separated=bool(diff), shortest_walk_half_length=kmin,
                               n_features=len(diff),
                               examples=[[list(map(int, k[1])), list(k[2]), v] for k, v in sorted(diff.items())
                                         if k[0] == kmin][:6])
    all_distinct = len({canon_multiset(Counter(h)) for h in hol.values()}) == len(orbitals)
    twists = {}
    for o in orbitals:
        vals = Counter(twist(Bs[m]) for m in np.flatnonzero(lab == o))
        twists[o] = {str(k): v for k, v in vals.items()}
    for a, b in pairs:
        if not sep[f"{a}-{b}"]["separated"]:
            ta, tb = twists[a], twists[b]
            sep[f"{a}-{b}"]["twist"] = dict(orbital_a=ta, orbital_b=tb,
                                            separated=(len(ta) == 1 and len(tb) == 1 and ta.keys() != tb.keys()
                                                       and "None" not in ta))
    # complete classification: holonomy multiset plus twist separates all 20
    full = {(canon_multiset(Counter(hol[o])), tuple(sorted(twists[o]))) for o in orbitals}
    return dict(pass_id=11189, orbitals=res, reversal_pairs=sep, twist_by_orbital=twists,
                holonomy_plus_twist_classifies_all_20=(len(full) == len(orbitals)),
                holonomy_classifies_all_20=all_distinct,
                tau_exchanges_every_reversal_pair=all(r["tau_image"] == r["reverse"] for r in res),
                tau_fixes_every_orbital=all(r["tau_image"] == r["orbital"] for r in res))


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1)
    for k in ("holonomy_classifies_all_20", "tau_exchanges_every_reversal_pair", "tau_fixes_every_orbital"):
        print(k, res[k])
    for r in res["orbitals"]:
        print(r)
    for k, v in res["reversal_pairs"].items():
        print(k, v)


if __name__ == "__main__":
    main()
