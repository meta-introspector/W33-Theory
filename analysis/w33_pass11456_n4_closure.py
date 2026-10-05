"""Pass 11456: the magic-axis law at four qutrits -- what is now settled, and the same-line cell.

FIXED-AXIS CELL (= Stab(z1), 549 conjugacy classes; Pass 11433): combine
  * the sparse decider (232 classes, 97.60% of H by class size, 0 failures), and
  * the magnitude theorem of Pass 11457 (proves F' for every class with Im(M - I) nondegenerate, and for frames orthogonal
    to the radical otherwise),
and report what remains open (classes neither decided nor fully covered by the theorem).

SAME-LINE CELL (M z1 = y on a line through z1, M^2 z1 = z1): its H-orbits are not enumerated here (that needs a GAP
orbit census over Sp(8,3)); instead 200 targeted classes (an exact n = 2 class of the cell embedded and conjugated by a
random element of Stab(z1), Pass 11421) are decided and checked against F' and against the theorem's coverage.
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11252_exact_reversibility as R  # noqa: E402
import w33_pass11350_linear_decider as L  # noqa: E402
import w33_pass11352_three_qutrit_geometry as GEO  # noqa: E402
import w33_pass11421_four_qutrits as F4  # noqa: E402
import w33_pass11433_magic_axis_law as A  # noqa: E402
import w33_pass11457_magic_axis_theorem as TH  # noqa: E402

OUT = ROOT / "data" / "w33_pass11456_n4_closure.json"


def fixed_cell_accounting():
    import ast
    lines = A.CLASSES.read_text().splitlines()
    D = F4.SparseDecider(4)
    L.CAP = 3 ** 11
    z, Om, lab = D.z1, D.wl.Om, D.wl.labels.astype(np.int64)
    tot = 0
    mass = Counter()
    cnt = Counter()
    for line in lines:
        m, size = line.rsplit(";", 1)
        M = np.array(ast.literal_eval(m), dtype=np.int64).T % 3
        size = int(size)
        tot += size
        dims = [len(s[1]) for s in (D.symplectic_solutions(M, k) for k in range(3)) if s is not None]
        reached = not dims or max(3 ** d for d in dims) <= L.CAP
        tf = TH.theorem_frames(M, z, Om, lab)
        wv = (lab @ Om @ ((z + M @ z) % 3)) % 3 != 0
        covered = tf is not None and bool(tf[0][wv].all())
        key = "decided by the decider" if reached else ("proved by the theorem" if covered else "open")
        cnt[key] += 1
        mass[key] += size
    return dict(classes=dict(cnt), mass_fraction={k: v / tot for k, v in mass.items()}, group_order=tot)


def _line_job(seed):
    F4._S.setdefault("D", None)
    D, gens = F4._S["D"], F4._S["gens"]
    rng = np.random.default_rng(seed)
    z, Om = D.z1, D.wl.Om
    lab = D.wl.labels.astype(np.int64)
    reps = F4._S["line_reps"]["M^2 z1 = z1"]
    Am = np.eye(8, dtype=np.int64)
    Am[:4, :4] = reps[int(rng.integers(len(reps)))]
    Am[4:, 4:] = GEO.random_symplectic(rng, GEO.gen_mats(2), length=100)
    h = GEO.random_symplectic(rng, gens, length=300)
    N = (GEO.map_to((h @ z) % 3, z, Om, rng) @ h) % 3
    M = (N @ Am @ R._inv_mod3(N)) % 3
    assert GEO.cell(M, z, Om) == "same line, M^2 z1 = z1"
    wv = (lab @ Om @ ((z + M @ z) % 3)) % 3 != 0
    tf = TH.theorem_frames(M, z, Om, lab)
    th = "theorem: fully" if tf is not None and tf[0][wv].all() else ("theorem: partly" if tf is not None else
                                                                       "theorem: silent")
    g = D.good_frames(M)
    if g is None:
        return ("undecided", th)
    bad = ~g
    return ("F' equality" if (bad == wv).all() else ("F' containment" if (bad >= wv).all() else "FAILS"), th)


def run():
    res = dict(pass_id=11456)
    res["fixed_axis_cell"] = fixed_cell_accounting()
    print(res["fixed_axis_cell"], flush=True)
    with Pool(6, initializer=F4._init) as pool:
        rows = pool.map(_line_job, range(114560000, 114560200), chunksize=4)
    res["same_line_cell_n4"] = dict(verdicts=dict(Counter(r[0] for r in rows)),
                                    theorem=dict(Counter(r[1] for r in rows)),
                                    joint=dict(Counter(f"{a} | {b}" for a, b in rows)))
    print(res["same_line_cell_n4"], flush=True)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
