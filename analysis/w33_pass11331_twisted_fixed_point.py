"""Pass 11331: one magic gate as a twisted fixed-point problem on the stabilizer of the magic axis -- an independent
second derivation of the census.

DERIVATION (exact).  Let U = C T1 (C Clifford, T1 = T (x) I) and suppose an anti-unitary Clifford reversal D K exists:
D U^* D^dag = lambda U^-1.  With U^* = Cbar T1^-1 this reads D Cbar T1^-1 D^dag = lambda T1^-1 C^-1.  Putting
E = D Cbar (a Clifford) and multiplying out:
        T1 E T1^-1  =  lambda C^-1 E C^T.                                  (*)
Conversely a Clifford E solving (*) gives the reversal D = E Cbar^-1.  So U is reversible iff (*) has a Clifford
solution E.  The right side is Clifford, so E must lie in N_T = {E Clifford : T1 E T1^-1 Clifford}.
FACT (computed): the symplectic parts of N_T are exactly the stabilizer of the magic axis z1 in Sp(2n,3) (two qutrits:
648 = |Stab(z1)|; one qutrit: the 3 unit shears), and every Pauli frame is allowed, so |N_T| = 81 * 648 (mod phase).
METHOD.  For each of the 4110 orbit representatives (M, a) of Pass 11330, all 52,488 candidates E = W(e) V_Q (Q z1 = z1)
are screened by the symplectic part of (*) (vectorised over F_3), and the survivors are checked as unitaries up to
phase.  This never calls the Weyl-coefficient criterion of Pass 11252, so agreement with Pass 11330 orbit by orbit is
an independent confirmation of the census 223/2430.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11252_exact_reversibility as R  # noqa: E402
import w33_pass11330_orbit_census as O  # noqa: E402

OUT = ROOT / "data" / "w33_pass11331_twisted_fixed_point.json"
T1 = O.T1
JJ = np.diag([1, 2, 1, 2])                                     # complex conjugation on labels: z -> -z


def sympl(wl, U):
    M = np.zeros((4, 4), int)
    for j in range(4):
        e = np.zeros(4, int)
        e[j] = 1
        Q = U @ wl.W[wl.index(e)] @ U.conj().T
        ov = np.abs(np.einsum('pij,ij->p', wl.W.conj(), Q)) / 9
        k = int(np.argmax(ov))
        if not np.isclose(ov[k], 1):
            return None
        M[:, j] = wl.labels[k]
    return M


def run():
    wl = R.Weyl(2)
    R.WEYL[2] = wl
    rng = np.random.default_rng(11331)
    Ms, keys = O.all_symplectic(wl)
    z = np.array([0, 1, 0, 0])
    stab = [i for i in range(len(Ms)) if tuple(Ms[i] @ z % 3) == tuple(z)]
    res = dict(pass_id=11331, stab_z1=len(stab))
    # N_T check: symplectic parts of N_T = Stab(z1)
    nt = []
    VQ = {}
    for i in range(len(Ms)):
        V = R.weil(wl, Ms[i], rng)
        if sympl(wl, T1 @ V @ T1.conj().T) is not None:
            nt.append(i)
            VQ[i] = V
    res["N_T_symplectic_classes"] = len(nt)
    res["N_T_equals_stab_z1"] = sorted(nt) == sorted(stab)
    print(res, flush=True)
    # candidates E = W(e) V_Q with LHS phi(E) = T1 E T1^-1 precomputed
    cand_E, cand_L, cand_Lsym = [], [], []
    for i in nt:
        for e in range(81):
            E = wl.W[e] @ VQ[i]
            L = T1 @ E @ T1.conj().T
            cand_E.append(E)
            cand_L.append(L)
            cand_Lsym.append(sympl(wl, L))
    cand_E = np.array(cand_E)
    cand_L = np.array(cand_L)
    cand_Lsym = np.array(cand_Lsym)
    cand_Q = np.array([Ms[i] for i in nt for _ in range(81)])
    # orbit representatives of Pass 11330
    import w33_pass11330_orbit_census as OC  # noqa: F401
    cert = json.load(open(ROOT / "data" / "w33_pass11330_orbit_census.json"))
    counts = np.load(ROOT / "data" / "w33_pass11330_class_counts.npy")
    A = np.array([(a // 27, a // 9 % 3, a // 3 % 3, a % 3) for a in range(81)])
    # rebuild the orbit labels exactly as Pass 11330 does, then test one representative per orbit
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    pow3 = np.array([27, 9, 3, 1])
    nM = len(Ms)
    rows, cols = [], []
    ids = np.arange(nM * 81).reshape(nM, 81)
    for N, d in O.k_generators(wl):
        Ninv = R._inv_mod3(N)
        imgM = np.array([keys[tuple(((N @ M @ Ninv) % 3).ravel())] for M in Ms])
        for i in range(nM):
            Mp = Ms[imgM[i]]
            rows.append(ids[i])
            cols.append(imgM[i] * 81 + ((A @ N.T + (d - Mp @ d)) % 3) @ pow3)
    G = coo_matrix((np.ones(sum(len(r) for r in rows), dtype=np.int8), (np.concatenate(rows), np.concatenate(cols))),
                   shape=(nM * 81, nM * 81)).tocsr()
    ncomp, lab = connected_components(G, directed=True, connection="weak")
    assert ncomp == cert["orbits"]
    first = np.full(ncomp, -1, dtype=np.int64)
    first[lab[np.argsort(lab, kind="stable")]] = np.argsort(lab, kind="stable")
    verdict = np.zeros(ncomp, bool)
    for o, el in enumerate(first):
        M = Ms[el // 81]
        a = el % 81
        C = wl.W[a] @ R.weil(wl, M, rng)
        Cinv = C.conj().T
        CT = C.T
        Minv = R._inv_mod3(M)
        target = (Minv[None] @ cand_Q @ (JJ @ Minv @ JJ)[None]) % 3      # symplectic part of C^-1 E C^T
        ok = np.flatnonzero((target == cand_Lsym).all(axis=(1, 2)))
        reversible = False
        for k in ok:
            Rm = Cinv @ cand_E[k] @ CT
            r = np.trace(Rm.conj().T @ cand_L[k]) / 9
            if abs(abs(r) - 1) < 1e-9 and np.allclose(cand_L[k], r * Rm, atol=1e-8):
                reversible = True
                break
        verdict[o] = not reversible
    counts2 = verdict[lab].reshape(nM, 81).sum(axis=1)
    res["orbits_tested"] = int(ncomp)
    res["violating_cliffords_by_twisted_method"] = int(counts2.sum())
    res["agrees_with_pass_11330_total"] = int(counts2.sum()) == cert["violating_cliffords_exact"]
    res["agrees_with_pass_11330_every_class"] = bool((counts2 == counts).all())
    print(res, flush=True)
    return res


def solvability():
    """split the classes by whether (*) has a solution at the SYMPLECTIC level (ignoring Pauli frames and phases)"""
    from collections import Counter
    wl = R.Weyl(2)
    R.WEYL[2] = wl
    rng = np.random.default_rng(1)
    Ms, keys = O.all_symplectic(wl)
    z = np.array([0, 1, 0, 0])
    H = [i for i in range(len(Ms)) if tuple(Ms[i] @ z % 3) == tuple(z)]
    Ls, Qs = [], []
    for i in H:
        V = R.weil(wl, Ms[i], rng)
        for e in range(81):
            Ls.append(sympl(wl, T1 @ wl.W[e] @ V @ T1.conj().T))
            Qs.append(Ms[i])
    pairs = np.unique(np.concatenate([np.array(Qs).reshape(-1, 16), np.array(Ls).reshape(-1, 16)], axis=1), axis=0)
    PQ, PL = pairs[:, :16].reshape(-1, 4, 4), pairs[:, 16:].reshape(-1, 4, 4)
    counts = np.load(ROOT / "data" / "w33_pass11330_class_counts.npy")
    solv = np.zeros(len(Ms), bool)
    for m, M in enumerate(Ms):
        Mi = R._inv_mod3(M)
        solv[m] = ((Mi[None] @ PQ @ (JJ @ Mi @ JJ)[None]) % 3 == PL).all(axis=(1, 2)).any()
    split = Counter(zip(solv.tolist(), counts.tolist()))
    return dict(distinct_symplectic_pairs=int(len(pairs)),
                split_solvable_by_frame_count={f"{'solvable' if s else 'unsolvable'}:{c}": int(n)
                                               for (s, c), n in sorted(split.items())},
                all_frames_violate_iff_symplectically_unsolvable=bool(((~solv) == (counts == 81)).all()))


def geometry():
    """the bad set in W(3,3) terms: where M sends the magic axis z1 (point of W(3,3)) and its preimage"""
    from collections import Counter, defaultdict
    wl = R.Weyl(2)
    Ms, keys = O.all_symplectic(wl)
    counts = np.load(ROOT / "data" / "w33_pass11330_class_counts.npy")
    z = np.array([0, 1, 0, 0])

    def om(u, v):
        return int(u @ wl.Om @ v) % 3

    tab = defaultdict(Counter)
    for i, M in enumerate(Ms):
        a = M @ z % 3
        b = R._inv_mod3(M) @ z % 3
        a2 = M @ a % 3
        if (a == z).all():
            key = "Mz1 = z1"
        elif ((a + z) % 3 == 0).all():
            key = "Mz1 = -z1"
        elif om(z, a) != 0:
            key = "Mz1 not collinear with z1"
        elif om(a, b) != 0:
            key = "Mz1, M^-1 z1 collinear with z1 on different lines"
        elif (a2 == z).all():
            key = "same line, M^2 z1 = z1"
        elif ((a2 + z) % 3 == 0).all():
            key = "same line, M^2 z1 = -z1"
        else:
            key = "same line, other"
        tab[key]["bad" if counts[i] > 0 else "good"] += 1
        tab[key]["all81"] += int(counts[i] == 81)
    return {k: dict(v) for k, v in sorted(tab.items())}


def main():
    if "--geometry-only" in sys.argv:
        res = json.load(open(OUT))
        res["w33_geometry_of_bad_set"] = geometry()
        json.dump(res, open(OUT, "w"), indent=1, default=str)
        print(json.dumps(res["w33_geometry_of_bad_set"], indent=1))
        return
    if "--solvability-only" in sys.argv:
        res = json.load(open(OUT))
        res["symplectic_level"] = solvability()
        json.dump(res, open(OUT, "w"), indent=1, default=str)
        print(json.dumps(res["symplectic_level"], indent=1))
        return
    res = run()
    res["symplectic_level"] = solvability()
    res["w33_geometry_of_bad_set"] = geometry()
    json.dump(res, open(OUT, "w"), indent=1, default=str)
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
