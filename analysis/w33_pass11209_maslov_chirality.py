#!/usr/bin/env python3
"""Pass 11209: the direction of time between two subsystem splits is a chirality, and a Berry phase records it.

Pass 11184 left two ties: for the reversal pairs {256, 256} and {6912, 6912} of three-qutrit split relations, a gate
and its inverse look identical to every entanglement measure used there.  This pass explains why and separates them.

1. TIME REVERSAL = MIRROR.  Complex conjugation of qutrit states acts on the symplectic data as the anti-symplectic map
   tau: (x, z) -> (x, -z) on every qutrit (the unique outer automorphism of PSp(6,3); it fixes the standard split F0).
   Computed on the 20 orbitals of Stab(F0) on the 110565 splits: tau FIXES 16 of them and SWAPS exactly the two tied
   reversal pairs, i.e. for those O^T = tau(O) -- the reversed relation is the mirror image of the relation.  The third
   reversal pair (2304) is tau-fixed, which is why the signature already separates it.  Consequently no invariant that
   is blind to the sign of the symplectic form (ranks, points, incidences, determinants, holonomy traces, ...) can
   separate the tied pairs: they are enantiomers.
2. THE CHIRAL INVARIANT.  For Lagrangians L1, L2, L3 the Kashiwara-Maslov index tau(L1, L2, L3) is the Witt class in
   W(F_3) = Z/4 of Q(x1,x2,x3) = w(x1,x2) + w(x2,x3) + w(x3,x1) on L1 + L2 + L3 (r - 2[disc = 2] mod 4); it is cyclic,
   alternating, and odd under tau.  For a pair of splits (F0, F) take the 64 product Lagrangians of each, give every
   F0-Lagrangian its PROFILE (sorted intersection dimensions with the 64 F-Lagrangians), and count ordered triples
   (L1, L2, M) with profile(L1) > profile(L2) -- a gauge-invariant orientation that alternation cannot cancel.  The
   spectrum (profiles, intersection dimensions, class) is constant on orbitals; its NET CHIRALITY
        chi(F0, F) = #{class 1} - #{class 3}
   is 0 on all 16 tau-fixed orbitals (as tau-invariance forces) and -54 / +54 on the 256 pair, +18 / -18 on the 6912
   pair.  With Pass 11184's signature and fine invariant, (signature, fine invariant, chi) classifies all 20 relations.
3. IT IS A BERRY PHASE.  For the zero-shift stabilizer states |L> (projector (1/27) sum_{v in L} D(v), D the
   symmetric Weyl operators), the Bargmann invariant satisfies
        tr(P_L1 P_L2 P_L3) = |tr(...)| * (-i)^{tau(L1, L2, L3)}           (checked on random triples, all agree),
   the finite-field case of the classical link between the Maslov index and the phase of the Weil representation.  So
   chi is minus the net imaginary part of the Pancharatnam phases of stabilizer-state triangles spanned by the two
   splits: the direction of time between two subsystem descriptions is recorded in geometric phase, and only there.
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
import w33_pass11182_paper_ticks_mereology as P  # noqa: E402
import w33_pass11184_merged_orbitals as O  # noqa: E402

OUT = ROOT / "data" / "w33_pass11209_maslov_chirality.json"
J6 = np.zeros((6, 6), np.int64)
for _k in range(3):
    J6[2 * _k, 2 * _k + 1], J6[2 * _k + 1, 2 * _k] = 1, -1
TAU = np.diag([1, 2, 1, 2, 1, 2]).astype(np.int64)
INV3 = np.array([0, 1, 2])


def rank_batch(X):
    """rank over F_3 of a batch of matrices (N, r, c)"""
    X = X.copy() % 3
    N, r, c = X.shape
    ar = np.arange(N)
    rank = np.zeros(N, np.int64)
    row = np.zeros(N, np.int64)
    for col in range(c):
        cand = (X[:, :, col] != 0) & (np.arange(r)[None, :] >= row[:, None])
        has = cand.any(1)
        p = np.argmax(cand, 1)
        rr = np.minimum(row, r - 1)
        tmp = X[ar, rr].copy()
        X[ar, rr] = np.where(has[:, None], X[ar, p], X[ar, rr])
        X[ar, p] = np.where(has[:, None], tmp, X[ar, p])
        piv = X[ar, rr, col]
        f = (X[:, :, col] * INV3[piv][:, None]) % 3
        f[ar, rr] = 0
        f = np.where(has[:, None], f, 0)
        X = (X - f[:, :, None] * X[ar, rr][:, None, :]) % 3
        rank += has
        row += has
    return rank


def witt_batch(A):
    """Witt class in W(F_3) = Z/4 of symmetric forms A (N, d, d): rank - 2 [disc = 2]  (mod 4)"""
    A = A.copy() % 3
    N, d, _ = A.shape
    ar = np.arange(N)
    r = np.zeros(N, np.int64)
    disc = np.ones(N, np.int64)
    for k in range(d):
        idx = np.arange(k, d)
        dg = A[:, idx, idx]
        sub = A[:, k:, k:] != 0
        hasoff = sub.any((1, 2)) & ~(dg != 0).any(1)
        if hasoff.any():                                   # no diagonal pivot: row/col p += row/col q
            flat = np.argmax(sub.reshape(N, -1), 1)
            p, q = flat // (d - k) + k, flat % (d - k) + k
            s = np.flatnonzero(hasoff)
            A[s, p[s], :] = (A[s, p[s], :] + A[s, q[s], :]) % 3
            A[s, :, p[s]] = (A[s, :, p[s]] + A[s, :, q[s]]) % 3
            dg = A[:, idx, idx]
        j = np.argmax(dg != 0, 1) + k
        ok = (dg != 0).any(1)
        perm = np.tile(np.arange(d), (N, 1))
        perm[ar, k] = j
        perm[ar, j] = k
        A = A[ar[:, None, None], perm[:, :, None], perm[:, None, :]]
        piv = A[:, k, k]
        f = (A[:, :, k] * INV3[piv][:, None]) % 3
        f[:, :k + 1] = 0
        A = (A - f[:, :, None] * A[:, k][:, None, :]) % 3
        A = (A - f[:, None, :] * A[:, :, k][:, :, None]) % 3
        r += ok
        disc = np.where(ok, (disc * np.where(piv == 0, 1, piv)) % 3, disc)
    return (r - 2 * (disc == 2)) % 4


def kashiwara(L1, L2, L3):
    """Kashiwara-Maslov index in Z/4 of batches of Lagrangians (N, 6, 3)"""
    W = [np.einsum('nai,ab,nbj->nij', X, J6, Y) % 3 for X, Y in ((L1, L2), (L2, L3), (L3, L1))]
    K = np.zeros((len(L1), 9, 9), np.int64)
    K[:, 0:3, 3:6], K[:, 3:6, 6:9], K[:, 6:9, 0:3] = W
    return witt_batch((2 * (K + K.transpose(0, 2, 1))) % 3)


def product_lagrangians(B):
    pts = [np.array(p) for p in ([1, 0], [0, 1], [1, 1], [1, 2])]
    out = []
    for c in itertools.product(range(4), repeat=3):
        L = np.zeros((6, 3), np.int64)
        for i in range(3):
            L[2 * i:2 * i + 2, i] = pts[c[i]]
        out.append((B @ L) % 3)
    return np.array(out)


def spectrum(B):
    """profile-oriented Maslov spectrum of the pair (F0, F = B F0): Counter over ordered triples (L1, L2, M)"""
    P0, PF = product_lagrangians(np.eye(6, dtype=np.int64)), product_lagrangians(B)
    i, m = (a.ravel() for a in np.meshgrid(np.arange(64), np.arange(64), indexing='ij'))
    dim = (6 - rank_batch(np.concatenate([P0[i], PF[m]], 2))).reshape(64, 64)
    prof0 = [tuple(sorted(r)) for r in dim.tolist()]
    profF = [tuple(sorted(c)) for c in dim.T.tolist()]
    pairs = [(a, b) for a in range(64) for b in range(64) if prof0[a] > prof0[b]]
    i1 = np.repeat(np.array([a for a, _ in pairs], np.int64), 64)
    i2 = np.repeat(np.array([b for _, b in pairs], np.int64), 64)
    mm = np.tile(np.arange(64), len(pairs))
    if len(pairs) == 0:
        return Counter()
    cls = kashiwara(P0[i1], P0[i2], PF[mm])
    return Counter((prof0[a], prof0[b], profF[c], int(dim[a, c]), int(dim[b, c]), int(k))
                   for a, b, c, k in zip(i1.tolist(), i2.tolist(), mm.tolist(), cls.tolist()))


def chirality(mu):
    return sum(v for k, v in mu.items() if k[-1] == 1) - sum(v for k, v in mu.items() if k[-1] == 3)


def negated(mu):
    return Counter({k[:-1] + ((-k[-1]) % 4,): v for k, v in mu.items()})


# ---------- Berry phase: zero-shift stabilizer states on three qutrits ----------
_W = np.exp(2j * np.pi / 3)
_X1 = np.roll(np.eye(3), 1, 0)
_Z1 = np.diag([1, _W, _W * _W])


def weyl(v):
    """symmetric Weyl operator D(x, z) = w^{2 x.z} X^x Z^z on (C^3)^3, v = (x1, z1, x2, z2, x3, z3)"""
    op = np.array([[1.0 + 0j]])
    ph = 0
    for q in range(3):
        x, z = int(v[2 * q]) % 3, int(v[2 * q + 1]) % 3
        op = np.kron(op, np.linalg.matrix_power(_X1, x) @ np.linalg.matrix_power(_Z1, z))
        ph += 2 * x * z
    return _W ** (ph % 3) * op


def stabilizer_projector(L):
    vecs = {tuple((L @ np.array(c)) % 3) for c in itertools.product(range(3), repeat=3)}
    return sum(weyl(np.array(v)) for v in vecs) / 27


def bargmann_check(n=200, seed=11189):
    """arg tr(P1 P2 P3) = -(pi/2) * kashiwara(L1, L2, L3), on random mixed product Lagrangians of random split pairs"""
    Bs = P.load_facs()
    rng = np.random.default_rng(seed)
    P0 = product_lagrangians(np.eye(6, dtype=np.int64))
    agree, classes = 0, Counter()
    for _ in range(n):
        PF = product_lagrangians(Bs[rng.integers(len(Bs))])
        L = [P0[rng.integers(64)], P0[rng.integers(64)], PF[rng.integers(64)]]
        rng.shuffle(L)
        c = int(kashiwara(L[0][None], L[1][None], L[2][None])[0])
        Pr = [stabilizer_projector(x) for x in L]
        t = np.trace(Pr[0] @ Pr[1] @ Pr[2])
        agree += bool(abs(t) > 1e-9 and np.isclose(t / abs(t), (-1j) ** c))
        classes[c] += 1
    return agree, n, {str(k): v for k, v in sorted(classes.items())}


def summarize(members=4, seed=11189):
    Bs = P.load_facs()
    orbs, index = O.orbits(Bs, return_index=True)
    orb_of = np.zeros(len(Bs), np.int64)
    for k, o in enumerate(orbs):
        orb_of[o] = k
    rng = np.random.default_rng(seed)
    rows = []
    for k, o in enumerate(orbs):
        rev = int(O.reverse_orbit(Bs[o[0]], index, orb_of))
        tau_img = sorted({int(orb_of[index[O.fac_key((TAU @ Bs[i]) % 3)]]) for i in o[:64]})
        picks = [o[0]] + [o[int(t)] for t in rng.integers(0, len(o), members - 1)]
        mus = [spectrum(Bs[i]) for i in picks]
        mu = mus[0]
        sig, fine = O.invariants(Bs[o[0]])
        rows.append(dict(orbital=k, size=len(o), reverse=rev, tau_image=tau_img, rep=Bs[o[0]].tolist(),
                         signature=[int(x) for x in sig], fine=[list(fine[0]), list(fine[1]), int(fine[2])],
                         spectrum_constant=all(m == mu for m in mus), members_checked=len(picks),
                         negation_symmetric=mu == negated(mu), chirality=chirality(mu),
                         triples=sum(mu.values()),
                         class_counts={str(c): sum(v for kk, v in mu.items() if kk[-1] == c) for c in range(4)}))
    by = {r['orbital']: r for r in rows}
    tau_swapped = sorted(tuple(sorted((r['orbital'], r['tau_image'][0]))) for r in rows
                         if r['tau_image'] != [r['orbital']])
    tau_swapped = sorted(set(tau_swapped))
    rev_pairs = sorted({tuple(sorted((r['orbital'], r['reverse']))) for r in rows if r['reverse'] != r['orbital']})
    keys = {(tuple(r['signature']), json.dumps(r['fine']), r['chirality']) for r in rows}
    agree, n, cls = bargmann_check()
    res = dict(pass_id=11209, orbitals=len(rows), orbit_sizes=sorted(r['size'] for r in rows),
               reversal_pairs=[list(p) for p in rev_pairs], tau_swapped_pairs=[list(p) for p in tau_swapped],
               tau_fixed=sum(1 for r in rows if r['tau_image'] == [r['orbital']]),
               tied_pairs_are_mirror_pairs=all(by[a]['reverse'] == b and by[a]['tau_image'] == [b]
                                               for a, b in tau_swapped),
               spectrum_constant_on_orbitals=all(r['spectrum_constant'] for r in rows),
               achiral_iff_tau_fixed=all((r['chirality'] == 0 and r['negation_symmetric'])
                                         == (r['tau_image'] == [r['orbital']]) for r in rows),
               chirality_by_size={f"{r['size']}#{r['orbital']}": r['chirality'] for r in rows if r['chirality']},
               complete_classification=len(keys) == len(rows),
               bargmann=dict(agree=agree, total=n, classes=cls), rows=rows)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    r = summarize()
    print({k: v for k, v in r.items() if k != 'rows'})
    for row in r['rows']:
        print(row['orbital'], row['size'], 'rev', row['reverse'], 'tau', row['tau_image'], 'chi', row['chirality'],
              'const', row['spectrum_constant'], 'negsym', row['negation_symmetric'])
