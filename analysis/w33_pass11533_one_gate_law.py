"""Pass 11533: the one-gate law -- which frames of a one-magic-gate tick violate, by linear algebra alone.

THEOREM (sufficiency, any n; the argument of Pass 11498).  Let U = W(a) V_M T1.  Call v "universally anti-fixed" if
M v = v, omega(v, z1) = 0, and Q J v = -v for every Q in every solution space of (S) (k = 0, 1, 2).  Then every frame with
omega(v, a) != 0 violates.  (omega(v, z1) = 0 makes the shear fix v and kills the T-twist; r_Q = 0; pairing (F) with
omega(v, .) leaves 2 omega(v, a) = 0.)  If (S) has no solution at all, every frame violates.
The universally anti-fixed vectors form a SUBSPACE W, cut out by linear equations in v (M v = v, omega(v, z1) = 0, and for
each solution space x0 + span(b_i): x0 J v = -v, b_i J v = 0).  So the predicted violating set is the complement of W^perp,
with 3^(2n) - 3^(2n - dim W) frames, and the predicted class verdict is: bad iff W != 0 or (S) has no solution.

CONJECTURE (necessity): these are ALL the violating frames, except for a small exceptional family.
Checked here against the exact deciders: every class at n = 2, every orbit of Pass 11373 at n = 3 (with exact weights).
Then the law is used as a FAST estimator of the bad-class fraction P_n and of the cell fractions for n up to 8.
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11350_linear_decider as L  # noqa: E402
import w33_pass11352_three_qutrit_geometry as GEO  # noqa: E402

OUT = ROOT / "data" / "w33_pass11533_one_gate_law.json"


def W_dim(D, M):
    """dim of the space of universally anti-fixed vectors, or None if (S) has no solution for any k"""
    N2 = D.N2
    z = D.z1
    Om = D.wl.Om
    spaces = [s for s in (D.symplectic_solutions(M, k) for k in range(3)) if s is not None]
    if not spaces:
        return None
    rows = [(M - np.eye(N2, dtype=np.int64)) % 3, (z @ Om)[None, :] % 3]
    rhs = [np.zeros(N2, np.int64), np.zeros(1, np.int64)]
    for x0, B in spaces:
        X0 = x0.reshape(N2, N2)
        rows.append((X0 @ D.J + np.eye(N2, dtype=np.int64)) % 3)
        rhs.append(np.zeros(N2, np.int64))
        for b in B:
            rows.append((np.asarray(b).reshape(N2, N2) @ D.J) % 3)
            rhs.append(np.zeros(N2, np.int64))
    sol = L.solve_affine(np.concatenate(rows) % 3, np.concatenate(rhs))
    return len(sol[1]) if sol is not None else 0


def validate_n2():
    import w33_pass11330_orbit_census as O
    L.CAP = 3 ** 11
    D = L.Decider(2)
    lab = D.wl.labels.astype(np.int64)
    st = Counter()
    for M in np.array(O.all_symplectic(D.wl)[0]) % 3:
        g = D.good_frames(M)
        if g is None:
            st["undecided"] += 1
            continue
        w = W_dim(D, M)
        pred_bad_frames = 81 if w is None else 81 - 3 ** (4 - w)
        actual = int((~g).sum())
        cell = GEO.cell(M, D.z1, D.wl.Om)
        st["class verdict agrees"] += (actual > 0) == (pred_bad_frames > 0)
        st["frame count agrees"] += actual == pred_bad_frames
        if actual != pred_bad_frames:
            st[f"frame-count exceptions in {cell} (law {pred_bad_frames}, actual {actual})"] += 1
        st["classes"] += 1
    return dict(st)


def validate_n3():
    import w33_pass11373_three_qutrit_exact_fraction as X
    L.CAP = 3 ** 11
    D = L.Decider(3)
    st = Counter()
    bad_mass = Fraction(0)
    tot = 0
    law_bad_mass = Fraction(0)
    for M, size in X.read_orbits(3):
        M = M % 3
        w = W_dim(D, M)
        law_bad = (w is None) or w > 0
        tot += size
        law_bad_mass += size * law_bad
        g = D.good_frames(M)
        if g is None:
            st["undecided by the decider"] += 1
            continue
        actual_bad = bool((~g).any())
        st["orbit verdict agrees"] += actual_bad == law_bad
        st["orbits decided"] += 1
        pred = 729 if w is None else 729 - 3 ** (6 - w)
        st["frame count agrees"] += int((~g).sum()) == pred
    return dict(counts=dict(st), law_bad_class_fraction=str(law_bad_mass / tot),
                exact_pass11373="437/3276", law_equals_exact=(law_bad_mass / tot) == Fraction(437, 3276))


def estimate(n, samples, seed):
    """fast law-based estimate of P_n and of the bad fraction in each cell (uniform M in Sp(2n, 3))"""
    D = L.Decider(n) if n <= 3 else __import__("w33_pass11421_four_qutrits").SparseDecider(n)
    rng = np.random.default_rng(seed)
    gens = GEO.gen_mats(n)
    st = Counter()
    for _ in range(samples):
        M = GEO.random_symplectic(rng, gens) % 3
        w = W_dim(D, M)
        bad = (w is None) or w > 0
        cell = GEO.cell(M, D.z1, D.wl.Om)
        st[(cell, bad)] += 1
        st[("all", bad)] += 1
    out = {}
    for cell in sorted({c for c, _ in st}):
        b, g = st[(cell, True)], st[(cell, False)]
        p = b / (b + g)
        out[cell] = dict(bad=b, good=g, fraction=p, stderr=(p * (1 - p) / (b + g)) ** 0.5)
    return out


def _est_job(args):
    n, samples, seed = args
    raw = estimate(n, samples, seed)
    return {c: (v["bad"], v["good"]) for c, v in raw.items()}


def estimate_parallel(n, samples, nproc=6, chunks=60):
    from multiprocessing import Pool
    per = samples // chunks
    tot = Counter()
    with Pool(nproc) as pool:
        for part in pool.imap_unordered(_est_job, [(n, per, 1153300 + 1000 * n + i) for i in range(chunks)]):
            for c, (b, g) in part.items():
                tot[(c, True)] += b
                tot[(c, False)] += g
    out = {}
    for cell in sorted({c for c, _ in tot}):
        b, g = tot[(cell, True)], tot[(cell, False)]
        p = b / (b + g)
        out[cell] = dict(bad=b, good=g, fraction=p, stderr=(p * (1 - p) / (b + g)) ** 0.5)
    if "non-collinear" in out:
        pred = (1 - 9.0 ** (1 - n)) / 8
        nc = out["non-collinear"]
        out["non-collinear"]["conjecture_(1-9^(1-n))/8"] = pred
        out["non-collinear"]["z_vs_conjecture"] = (nc["fraction"] - pred) / nc["stderr"]
    return out


def main():
    res = json.load(open(OUT)) if OUT.exists() else dict(pass_id=11533)
    stage = sys.argv[1] if len(sys.argv) > 1 else "validate"
    if stage == "validate":
        n2 = validate_n2()
        print(n2, flush=True)
        n3 = validate_n3()
        print(n3, flush=True)
        res = json.load(open(OUT)) if OUT.exists() else dict(pass_id=11533)       # merge with whatever is there now
        res["n2"], res["n3"] = n2, n3
        json.dump(res, open(OUT, "w"), indent=1, default=str)
    else:
        n, samples = int(sys.argv[2]), int(sys.argv[3])
        est = estimate_parallel(n, samples)
        print(n, est, flush=True)
        res = json.load(open(OUT)) if OUT.exists() else dict(pass_id=11533)       # merge with whatever is there now
        res.setdefault("estimates", {})[str(n)] = est
        json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
