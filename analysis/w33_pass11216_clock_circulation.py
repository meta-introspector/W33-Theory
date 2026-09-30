"""Pass 11216: the one-way clock circulates -- current spectrum and oriented triangles of the time-oriented relations.

Pass 11205 computed the intersection matrices L_R of the rank-20 configuration of three-qutrit splits and found non-real
spectra for the six directed relations (reverse = conjugate).  For a directed relation R with reverse R^T:
  * the current operator J_R = A_R - A_{R^T} is real skew-symmetric; A_R and A_{R^T} commute (Pass 11205), so on
    their joint eigenspaces J_R = lambda - conj(lambda) = 2 i Im(lambda): the circulation frequencies;
  * oriented triangles x -> z -> y -> x along R through a point number |R| p^{R^T}_{R R} = |R| (L_R)[R^T, R], and
    closed walks of length 3 along R minus those along R^T measure the net circulation (tr A_R^3 = tr A_{R^T}^3 for a
    real matrix, so the one-way signal is in the mixed cycles);
  * compared with the Kashiwara-Maslov chirality chi of Pass 11209 (merged from the parallel session): 256 pair
    chi = -54/+54, 6912 pair +18/-18.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "w33_pass11216_clock_circulation.json"
SRC = ROOT / "data" / "w33_pass11205_one_way_clock.json"


def run():
    d = json.load(open(SRC))
    L = {int(k): np.array(v, dtype=np.int64) for k, v in d["L"].items()}
    size = {int(k): v for k, v in d["sizes"].items()}
    res = dict(pass_id=11216)
    for a, b in ((3, 4), (12, 13), (7, 9)):
        La, Lb = L[a], L[b]
        Jm = (La - Lb).astype(float)
        ev = np.linalg.eigvals(Jm)
        freqs = sorted({round(abs(e.imag) / np.sqrt(3), 6) for e in ev if abs(e.imag) > 1e-8})
        tri_forward = int(size[a] * La[b, a])          # x->z->y->x all along R: (x,y) in R^T
        tri_mixed = int(size[a] * La[a, a])            # x->z->y along R with (x,y) in R: transitive triangles
        res[f"{a}-{b}"] = dict(size=size[a],
                               current_eigenvalues_over_i_sqrt3=freqs,
                               oriented_3_cycles_per_point=tri_forward,
                               transitive_triangles_per_point=tri_mixed,
                               commute=bool(np.array_equal(La @ Lb, Lb @ La)))
    try:
        m = json.load(open(ROOT / "data" / "w33_pass11209_maslov_chirality.json"))
        res["maslov_file_keys"] = list(m.keys())[:12]
    except FileNotFoundError:
        res["maslov_file_keys"] = None
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
