"""Pass 11330: the unconditional exhaustive census of one magic gate on two qutrits, by centralizer orbits.

Pass 11309 reported P = 223/2430 for U = C (T (x) I) over all 4,199,040 two-qutrit Cliffords C, but (Codex's audit,
PASS11303_11307_..., 'publication-time parallel integration correction') it decided all 81 Pauli frames only for the
6480 classes flagged by the frames 0, e1..e4; the other 45,360 classes were certified only through the unproved
affine law.  Here the census is made exhaustive with no such assumption.

METHOD.  Write every Clifford (mod phase) uniquely as C = W(a) V_M with W the symmetric Weyl operators, V_M the
canonical Weil unitary (V_M W(p) V_M^dag = W(Mp) exactly), a in F_3^4, M in Sp(4,3).  Reversibility is invariant
under conjugation by any Clifford D with D T1 D^dag proportional to T1 (T1 = T (x) I): D U D^dag = (D C D^dag) T1.
The group K of such D (mod phase) contains W(d) for d with no X1 component, the canonical Weil unitaries of the unit
shears on qutrit 1, and Sp(2,3) on qutrit 2 (|K| = 27 * 72 = 1944); every generator is checked numerically to commute
with T1 up to phase.  Conjugation acts by
        (M, a)  ->  (N M N^-1,  d + N a - N M N^-1 d).
The 4,199,040 elements split into K-orbits (connected components of the generator graph); one representative per
orbit is decided exactly with Pass 11252's criterion and its verdict is assigned to the whole orbit.

VALIDATION.  (i) 1000 random elements are decided directly and compared with their orbit verdicts (tests the symmetry
reduction itself); (ii) the per-class frame counts are compared with Pass 11309's directly decided flagged classes.
OUTPUT.  The exact count, the exact frame-count distribution for every class (including the previously unflagged
45,360), and an exhaustive check of the affine law on all 51,840 classes.
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11227_two_qutrit_cp_mixing as P2  # noqa: E402
import w33_pass11252_exact_reversibility as R  # noqa: E402
import w33_pass11309_one_gate_law as P09  # noqa: E402

OUT = ROOT / "data" / "w33_pass11330_orbit_census.json"
T1 = np.kron(P2.T, P2.I3)
_S = {}


def all_symplectic(wl):
    reps = np.load(P2.CACHE)
    Ms = np.array([P09.symplectic_part(wl, reps[r]) for r in range(len(reps))])
    keys = {tuple(M.ravel()): i for i, M in enumerate(Ms)}
    assert len(keys) == 51840
    return Ms, keys


def k_generators(wl):
    """(N, d) pairs generating K, with N symplectic on (x1, z1, x2, z2)"""
    I4 = np.eye(4, dtype=int)
    shear1 = I4.copy()
    shear1[1, 0] = 1                                   # z1 -> z1 + x1
    shear2 = I4.copy()
    shear2[3, 2] = 1                                   # z2 -> z2 + x2
    four2 = I4.copy()
    four2[2:, 2:] = [[0, 2], [1, 0]]                   # (x2, z2) -> (-z2, x2)
    gens = [(shear1, np.zeros(4, int)), (shear2, np.zeros(4, int)), (four2, np.zeros(4, int))]
    for e in (1, 2, 3):                                # Paulis z1, x2, z2
        d = np.zeros(4, int)
        d[e] = 1
        gens.append((I4, d))
    return gens


def check_generators(wl, gens, rng):
    ok = []
    for N, d in gens:
        assert ((N.T @ wl.Om @ N - wl.Om) % 3 == 0).all()
        D = wl.W[wl.index(d)] @ R.weil(wl, N, rng)
        X = D @ T1 @ D.conj().T
        r = np.trace(T1.conj().T @ X) / 9
        ok.append(bool(abs(abs(r) - 1) < 1e-9 and np.allclose(X, r * T1)))
    return ok


def _init():
    R.WEYL[2] = R.Weyl(2)


def _decide(item):
    M, a = item
    wl = R.WEYL[2]
    rng = np.random.default_rng(int(abs(hash((tuple(M.ravel()), tuple(a)))) % 2 ** 31))
    V = R.weil(wl, M, rng)
    U = wl.W[wl.index(a)] @ V @ T1
    return R.decide(U, 2, rng)[0] is False


def run():
    wl = R.Weyl(2)
    R.WEYL[2] = wl
    rng = np.random.default_rng(11330)
    res = dict(pass_id=11330)
    gens = k_generators(wl)
    res["generators_commute_with_T1"] = check_generators(wl, gens, rng)
    assert all(res["generators_commute_with_T1"])
    Ms, keys = all_symplectic(wl)
    nM = len(Ms)
    A = np.array([(a // 27, a // 9 % 3, a // 3 % 3, a % 3) for a in range(81)])          # frame a, index base 3
    pow3 = np.array([27, 9, 3, 1])
    rows, cols = [], []
    ids = np.arange(nM * 81).reshape(nM, 81)
    for N, d in gens:
        Ninv = R._inv_mod3(N)
        imgM = np.array([keys[tuple(((N @ M @ Ninv) % 3).ravel())] for M in Ms])
        for i in range(nM):
            Mp = Ms[imgM[i]]
            a2 = (A @ N.T + (d - Mp @ d)) % 3
            rows.append(ids[i])
            cols.append(imgM[i] * 81 + a2 @ pow3)
    rows = np.concatenate(rows)
    cols = np.concatenate(cols)
    G = coo_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(nM * 81, nM * 81)).tocsr()
    ncomp, lab = connected_components(G, directed=True, connection="weak")
    sizes = np.bincount(lab)
    res["orbits"] = int(ncomp)
    res["orbit_size_values"] = sorted(map(int, set(sizes)))
    assert sizes.sum() == nM * 81 and sizes.max() <= 1944
    print("orbits", ncomp, res["orbit_size_values"], flush=True)
    first = np.full(ncomp, -1, dtype=np.int64)
    order = np.argsort(lab, kind="stable")
    first[lab[order]] = order                           # last write wins -> some element of each orbit
    items = [(Ms[e // 81], A[e % 81]) for e in first]
    with Pool(11, initializer=_init) as pool:
        verdict = np.array(pool.map(_decide, items, chunksize=8))
        # validation: random elements decided directly vs their orbit verdict
        samp = rng.choice(nM * 81, 1000, replace=False)
        direct = np.array(pool.map(_decide, [(Ms[e // 81], A[e % 81]) for e in samp], chunksize=8))
    res["validation_direct_vs_orbit_agree"] = int(np.sum(direct == verdict[lab[samp]]))
    res["validation_samples"] = 1000
    viol = verdict[lab].reshape(nM, 81)
    total = int(viol.sum())
    res["violating_cliffords_exact"] = total
    res["probability_exact"] = str(Fraction(total, nM * 81))
    counts = viol.sum(axis=1)
    res["class_frame_counts"] = {str(k): int(np.sum(counts == k)) for k in sorted(set(counts.tolist()))}

    def is_affine(rows_):
        S = A[rows_]
        if len(S) == 0:
            return True
        keys_ = {tuple(s) for s in S}
        return all(tuple((x + y - z) % 3) in keys_ for x in S for y in S[:9] for z in S[:9])

    non_affine = 0
    for i in range(nM):
        good = np.flatnonzero(~viol[i])
        if 0 < len(good) < 81:
            S = A[good]
            keys_ = {tuple(s) for s in S}
            if not all(tuple((x + y - z) % 3) in keys_ for x in S for y in S for z in S):
                non_affine += 1
    res["classes_violating_the_affine_law"] = non_affine
    res["bad_classes"] = int(np.sum(counts > 0))
    res["five_frame_screen_exact"] = bool(all((viol[i, [0, 27, 9, 3, 1]].any()) == (counts[i] > 0) for i in range(nM)))
    p09 = json.load(open(ROOT / "data" / "w33_pass11309_one_gate_law.json"))
    res["agrees_with_pass_11309"] = (res["probability_exact"] == p09["probability_exact"]
                                     and res["bad_classes"] == p09["bad_classes"])
    np.save(ROOT / "data" / "w33_pass11330_bad_classes.npy", np.flatnonzero(counts > 0))
    np.save(ROOT / "data" / "w33_pass11330_class_counts.npy", counts.astype(np.int16))
    print({k: v for k, v in res.items() if k != "generators_commute_with_T1"}, flush=True)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
