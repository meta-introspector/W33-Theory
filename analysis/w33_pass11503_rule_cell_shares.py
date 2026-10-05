"""Pass 11503: the exact shares of the four one-gate rule cells, for every n.

With N = 3^(2n) - 1 (the nonzero vectors of F3^(2n)), M uniform in Sp(2n, 3) and z1 fixed:
    Pr(M^2 z1 = z1, v = z1 + M z1 != 0) = 2/N,   Pr(M z1 = -z1) = 1/N,   Pr(M^2 z1 = -z1) = 3/N,   Pr(M z1 = z1) = 1/N.
DERIVATION.  M z1 is uniform on the N nonzero vectors, and Sp(2n, 3) is transitive on ordered pairs (u, w) with a given
omega(u, w) and u, w independent (Witt).  For M^2 z1 = +-z1 with y = M z1:
  * y = +-z1: one choice each (M z1 = z1 lies in the M^2 z1 = z1 cell with v = 2 z1 != 0; M z1 = -z1 gives v = 0, and
    M^2 z1 = z1, so it is its own cell).
  * y independent of z1: omega(z1, y) = omega(M z1, M y) = omega(y, +-z1) forces omega(z1, y) = 0 for the + sign and leaves
    it free for the - sign.  For each admissible y, Pr(M z1 = y, M y = +-z1) = 1/#(images of the pair (z1, y)); summing over
    y gives 1/N for the isotropic pairs and 2/N for the two non-isotropic values of omega.
So M^2 z1 = z1 with v != 0 collects 1/N (y = z1) + 1/N (isotropic y); M^2 z1 = -z1 collects 1/N + 2/N.
Checked: n = 2 by counting all 51,840 elements; n = 3 by Pass 11373's exact orbit sizes over all 9,170,703,360 elements.
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11350_linear_decider as L  # noqa: E402

OUT = ROOT / "data" / "w33_pass11503_rule_cell_shares.json"


def cells(M, z):
    Mz = M @ z % 3
    M2 = M @ Mz % 3
    out = []
    if ((Mz - z) % 3 == 0).all():
        out.append("M z1 = z1")
    if ((M2 - z) % 3 == 0).all() and ((z + Mz) % 3).any():
        out.append("M^2 z1 = z1, v != 0")
    if ((Mz + z) % 3 == 0).all():
        out.append("M z1 = -z1")
    if ((M2 + z) % 3 == 0).all():
        out.append("M^2 z1 = -z1")
    return out


def shares(items, z, n):
    tot = 0
    acc = {k: 0 for k in ("M z1 = z1", "M^2 z1 = z1, v != 0", "M z1 = -z1", "M^2 z1 = -z1")}
    for M, size in items:
        tot += size
        for c in cells(np.asarray(M) % 3, z):
            acc[c] += size
    N = 3 ** (2 * n) - 1
    pred = {"M z1 = z1": Fraction(1, N), "M^2 z1 = z1, v != 0": Fraction(2, N), "M z1 = -z1": Fraction(1, N),
            "M^2 z1 = -z1": Fraction(3, N)}
    return dict(group_order=tot, shares={k: str(Fraction(v, tot)) for k, v in acc.items()},
                predicted={k: str(v) for k, v in pred.items()},
                all_match=all(Fraction(acc[k], tot) == pred[k] for k in acc))


def run():
    import w33_pass11330_orbit_census as O
    import w33_pass11373_three_qutrit_exact_fraction as X
    D2 = L.Decider(2)
    res = dict(pass_id=11503)
    res["n2"] = shares([(M, 1) for M in O.all_symplectic(D2.wl)[0]], D2.z1, 2)
    D3 = L.Decider(3)
    res["n3"] = shares(X.read_orbits(3), D3.z1, 3)
    print(json.dumps(res, indent=1), flush=True)
    return res


def main():
    res = run()
    assert res["n2"]["all_match"] and res["n3"]["all_match"]
    json.dump(res, open(OUT, "w"), indent=1)


if __name__ == "__main__":
    main()
