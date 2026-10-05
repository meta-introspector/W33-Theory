"""Pass 11490: re-audit of the fixed numerical cuts in the sampled passes, after Pass 11459's threshold artefact.

Pass 11459 found that a cut which separates cleanly at shallow depth (overlap > 2.999) fails at depth, because
Clifford+T words are dense and violators come arbitrarily close to the reversible value.  The same mechanism can bite
every SAMPLED decision with a hard cut.  Audited here:

  A. Pass 11369 (deep words): "J6 is complete on deep words" -- a violator counts as a J6 zero if J6 < 1e-9.  Its committed
     minimum violating J6 fell from 1.5e-4 (k = 5) to 9.1e-9 (k = 16), within one decade of the cut.
  B. The deep-word decider itself (R.decide, Pass 11252's Weyl criterion, used by 11312 and 11369) against the
     INDEPENDENT overlap decider of Theorem 1 (Pass 11355) at the tight cut 3 - 1e-7 (Pass 11459).
  C. Pass 11353 (level statistics): an eigenphase spacing counts as degenerate if s < 1e-9.

For each cut we report the two one-sided extremes (largest value on the 'zero' side, smallest on the 'positive' side) at
depths beyond those used, and whether a clean gap of several decades still surrounds the cut.
"""

from __future__ import annotations

import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11357_jarlskog_degree_by_dimension as P7  # noqa: E402

OUT = ROOT / "data" / "w33_pass11490_cut_audit.json"
TIGHT = 3 - 1e-7


def _job(args):
    k, seed, n = args
    import w33_pass11213_cubic_t_violation as P1
    import w33_pass11252_exact_reversibility as R
    import w33_pass11422_depth5_exact as E
    R.WEYL[1] = R.Weyl(1)
    E._init()
    CF = np.array(P1.clifford1())
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(n):
        U = np.eye(3, dtype=complex)
        for c in rng.integers(216, size=k):
            U = CF[c] @ P1.T1 @ U
        try:
            v = R.decide(U, 1, rng)[0]
        except AssertionError:
            v = "ambiguous"
        ov = float(E.overlaps(U[None])[0])
        rows.append((v if isinstance(v, str) else (None if v is None else bool(v)), 3 - ov, P7.J(U, 3, CF)))
    return k, rows


def deep_audit(ks=(12, 16, 20, 24, 30), per_k=24000, nproc=6, seed=114900000):
    jobs = [(k, seed + 1000 * k + i, per_k // nproc) for k in ks for i in range(nproc)]
    acc = {k: [] for k in ks}
    with Pool(nproc) as pool:
        for k, rows in pool.imap_unordered(_job, jobs):
            acc[k] += rows
    out = {}
    for k in ks:
        rows = acc[k]
        rd = np.array([r[0] for r in rows], dtype=object)
        deficit = np.array([r[1] for r in rows])
        j6 = np.array([r[2] for r in rows])
        ov_rev = deficit < 1e-7
        decided = np.array([x in (True, False) for x in rd])
        agree = int(sum(1 for x, o in zip(rd, ov_rev) if x in (True, False) and x == o))
        viol, rev = ~ov_rev, ov_rev
        out[str(k)] = dict(
            words=len(rows),
            reversible_by_overlap=int(rev.sum()),
            R_decide_decided=int(decided.sum()),
            R_decide_ambiguous=int(sum(1 for x in rd if x == "ambiguous")),
            R_decide_undecided=int(sum(1 for x in rd if x is None)),
            deciders_agree_on_decided=agree,
            deciders_disagree=int(decided.sum()) - agree,
            overlap_gap=dict(max_deficit_reversible=float(deficit[rev].max()) if rev.any() else None,
                             min_deficit_violator=float(deficit[viol].min())),
            j6_gap=dict(max_abs_j6_reversible=float(np.abs(j6[rev]).max()) if rev.any() else None,
                        min_j6_violator=float(j6[viol].min()),
                        violators_below_1e9=int((j6[viol] < 1e-9).sum()),
                        violators_below_1e7=int((j6[viol] < 1e-7).sum()),
                        violators_below_1e5=int((j6[viol] < 1e-5).sum())),
            j6_negative_violators=int((j6[viol] < 0).sum()),
        )
        print(k, out[str(k)], flush=True)
    return out


def level_audit(trials=4000, seed=11490):
    """Pass 11353's degenerate-spacing cut on its own two-qutrit magic circuits"""
    import w33_pass11353_level_statistics as LS
    rng = np.random.default_rng(seed)
    smin = []
    for layers in (4, 8, 16, 32):
        for _ in range(trials // 4):
            U = LS.magic_circuit(2, rng, layers)
            ph = np.sort(np.angle(np.linalg.eigvals(U)))
            s = np.diff(np.concatenate([ph, [ph[0] + 2 * np.pi]]))
            smin.append(float(s.min()))
    smin = np.array(smin)
    return dict(trials=len(smin), below_1e9=int((smin < 1e-9).sum()),
                between_1e12_and_1e6=int(((smin > 1e-12) & (smin < 1e-6)).sum()),
                max_below_1e9=float(smin[smin < 1e-9].max()) if (smin < 1e-9).any() else None,
                min_above_1e9=float(smin[smin >= 1e-9].min()))


def run():
    res = dict(pass_id=11490, tight_cut=TIGHT)
    res["deep_words"] = deep_audit()
    res["level_statistics"] = level_audit()
    print(res["level_statistics"], flush=True)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
