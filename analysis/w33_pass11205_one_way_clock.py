"""Pass 11205: the one-way clock -- the spectrum of the time-oriented relations between three-qutrit splits.

The 110565 three-qutrit tensor factorisations carry a rank-20 coherent configuration under Sp(6,3) (Passes 11178,
11184).  Six of its relations are directed: R and its time reverse R^T differ (orbitals 3<->4 of size 256, 7<->9 of
size 2304, 12<->13 of size 6912; Pass 11189).  The adjacency matrix A_R of a directed relation is not symmetric, and
a walk that repeatedly steps along R never steps back.

The intersection numbers p^j_{ik} = #{z : (x,z) in R_i, (z,y) in R_k} for fixed (x,y) in R_j give the matrix
(L_i)_{jk} = p^j_{ik} of multiplication by A_i on the Bose-Mesner span; its eigenvalues are eigenvalues of A_i.  This
pass computes L_i for every directed relation (and one symmetric control), their spectra, and whether they commute.
Complex eigenvalues are the frequencies of a one-way rotation: the walk along R circulates, the walk along R^T turns
the other way (conjugate spectrum).
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11182_paper_ticks_mereology as PT  # noqa: E402
import w33_pass11184_merged_orbitals as MO  # noqa: E402

OUT = ROOT / "data" / "w33_pass11205_one_way_clock.json"
CACHE = Path.home() / "AppData" / "Local" / "Temp" / "orbit_labels3.npy"
J6 = np.zeros((6, 6), np.int64)
for _k in range(3):
    J6[2 * _k, 2 * _k + 1], J6[2 * _k + 1, 2 * _k] = 1, -1


def run(relations=(3, 4, 7, 9, 12, 13, 1)):
    Bs = PT.load_facs()
    lab = np.load(CACHE)
    index = {MO.fac_key(B): i for i, B in enumerate(Bs)}
    orbitals = sorted(set(lab.tolist()))
    size = Counter(lab.tolist())
    rep = {o: int(np.flatnonzero(lab == o)[0]) for o in orbitals}
    L = {}
    for i in relations:
        Mi = np.zeros((20, 20), np.int64)
        for z in np.flatnonzero(lab == i):
            Binv = (-J6 @ Bs[z].T @ J6) % 3
            for j in orbitals:
                k = int(lab[index[MO.fac_key((Binv @ Bs[rep[j]]) % 3)]])
                Mi[j, k] += 1
        assert (Mi.sum(1) == size[i]).all()          # every row counts all z in R_i(x)
        L[i] = Mi
    res = dict(pass_id=11205, sizes={int(k): int(v) for k, v in size.items()})
    spec = {}
    for i, Mi in L.items():
        ev = np.linalg.eigvals(Mi.astype(float))
        ev = ev[np.argsort(-np.abs(ev))]
        spec[i] = [[round(float(e.real), 6), round(float(e.imag), 6)] for e in ev]
    res["spectra"] = {str(k): v for k, v in spec.items()}
    res["reverse_is_conjugate"] = {
        f"{a}-{b}": bool(np.allclose(sorted(np.linalg.eigvals(L[a].astype(float)), key=lambda z: (round(z.real, 6), round(z.imag, 6))),
                                     sorted(np.conj(np.linalg.eigvals(L[b].astype(float))), key=lambda z: (round(z.real, 6), round(z.imag, 6))), atol=1e-6))
        for a, b in ((3, 4), (7, 9), (12, 13))}
    res["directed_have_complex_eigenvalues"] = {str(i): bool(np.abs(np.array(spec[i])[:, 1]).max() > 1e-6) for i in L}
    res["commute_3_12"] = bool(np.array_equal(L[3] @ L[12], L[12] @ L[3]))
    res["commute_3_4"] = bool(np.array_equal(L[3] @ L[4], L[4] @ L[3]))
    res["L"] = {str(k): v.tolist() for k, v in L.items()}
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1)
    for k in ("spectra", "reverse_is_conjugate", "directed_have_complex_eigenvalues", "commute_3_12", "commute_3_4"):
        print(k, res[k])


if __name__ == "__main__":
    main()
