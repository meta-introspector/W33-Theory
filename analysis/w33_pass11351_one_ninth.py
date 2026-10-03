"""Pass 11351: the 1/9 on the non-collinear cell = (2/9) x (1/2): a twisted-centralizer count times the parity -I.

Pass 11331 split the 6480 bad two-qutrit classes over the W(3,3) geometry of the magic axis z1; on the cell 'M z1 not
collinear with z1' (34,992 classes) exactly 3888 = 1/9 are bad, uniformly across the invariants tried.  With the
F3-linear decider of Pass 11350 the reason is visible:
  * on this cell only k = 0 admits symplectic solutions Q of  M Q (J M J) = Q, Q z1 = z1;
  * the number of symplectic solutions is 1, 3, 4, 6 or 24 -- and the class is bad ONLY when it is exactly 3;
  * the 3 solutions always differ by a transvection (rank(Q1 Q0^-1 - 1) = 1): a coset of a transvection subgroup;
  * of these 2/9 of the cell, exactly half are bad, and the halves are exchanged by the central involution -I
    (M bad <=> -M good), while badness is invariant under M -> M^-1 and M -> J M J.
So 1/9 = (2/9)(1/2).  The -I exchange also maps 'z1 fixed' (all bad) onto 'z1 reversed' (all good) and preserves the
single-line cells ((-M)^2 = M^2); one qutrit shows the same exchange (unit shears bad, parity-twisted shears good).
"""

from __future__ import annotations

import itertools
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11252_exact_reversibility as R  # noqa: E402
import w33_pass11330_orbit_census as O  # noqa: E402
import w33_pass11350_linear_decider as L  # noqa: E402

OUT = ROOT / "data" / "w33_pass11351_one_ninth.json"


def run():
    D = L.Decider(2)
    wl = D.wl
    Ms, keys = O.all_symplectic(wl)
    counts = np.load(ROOT / "data" / "w33_pass11330_class_counts.npy")
    bad = counts > 0
    z = D.z1
    Om = wl.Om
    inv = R._inv_mod3
    J = np.diag([1, 2, 1, 2])

    def om(u, v):
        return int(u @ Om @ v) % 3

    def idx(X):
        return keys[tuple((X % 3).ravel())]

    def cell(M):
        a = M @ z % 3
        b = inv(M) @ z % 3
        if (a == z).all():
            return "fix"
        if ((a + z) % 3 == 0).all():
            return "rev"
        if om(z, a):
            return "non-collinear"
        if om(a, b):
            return "collinear, different lines"
        a2 = M @ a % 3
        if (a2 == z).all():
            return "same line, M^2 z1 = z1"
        if ((a2 + z) % 3 == 0).all():
            return "same line, M^2 z1 = -z1"
        return "same line, other"

    sol_table = defaultdict(Counter)
    k_nonzero = 0
    transvection_coset = 0
    three = []
    for i, M in enumerate(Ms):
        if cell(M) != "non-collinear":
            continue
        for k in (1, 2):
            if D.symplectic_solutions(M, k) is not None:
                k_nonzero += 1
        x0, b = D.symplectic_solutions(M, 0)
        sols = []
        for cs in itertools.product(range(3), repeat=len(b)):
            q = x0.copy()
            for c, v in zip(cs, b):
                q = (q + c * v) % 3
            Q = q.reshape(4, 4)
            if D.is_symplectic(Q):
                sols.append(Q)
        sol_table[len(sols)]["bad" if bad[i] else "good"] += 1
        if len(sols) == 3:
            three.append(i)
            g = (sols[1] @ inv(sols[0])) % 3
            transvection_coset += len(L.rref3((g - np.eye(4, dtype=np.int64)) % 3)[1]) == 1
    pair = defaultdict(Counter)
    for i, M in enumerate(Ms):
        pair[cell(M)][f"{int(bad[i])}{int(bad[idx(-M)])}"] += 1
    sym = Counter()
    for i in three:
        M = Ms[i]
        sym["-M complementary"] += bad[i] != bad[idx(-M)]
        sym["M^-1 same"] += bad[i] == bad[idx(inv(M))]
        sym["JMJ same"] += bad[i] == bad[idx(J @ M @ J)]
    res = dict(pass_id=11351,
               noncollinear_classes=int(sum(sum(v.values()) for v in sol_table.values())),
               classes_with_k_nonzero_solutions=k_nonzero,
               symplectic_solution_count_vs_verdict={str(k): dict(v) for k, v in sorted(sol_table.items())},
               three_solution_classes=len(three), three_solutions_differ_by_transvection=int(transvection_coset),
               three_solution_symmetries={k: int(v) for k, v in sym.items()},
               minus_identity_pairing_by_cell={k: dict(v) for k, v in pair.items()})
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
