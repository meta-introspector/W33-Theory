"""Pass 11333: one magic gate on three qutrits -- is the bad fraction of symplectic classes still 1/8?

One qutrit (Pass 11266): bad classes 3/24 = 1/8.  Two qutrits (Pass 11330, exhaustive): bad classes 6480/51840 = 1/8.
Here Sp(6,3) (|Sp(6,3)| = 9,170,703,360) is sampled.  For a class M the canonical Weil unitary V_M is built by the twirl;
U = W(a) V_M (T (x) I (x) I) is decided exactly (Pass 11252) for:
  * all 729 frames a, for a subset of classes -- this tests the affine law at n = 3 directly (no five-frame shortcut);
  * the 7 frames 0, e_1..e_6 for every sampled class -- a class is flagged bad if one of them violates.  By the affine
    law (checked on the subset; proved exhaustively at n = 2 by Pass 11330 only) a class is bad iff flagged; the
    flagged fraction is reported with that condition stated.
Sampling: M from long random words in the generators of Sp(6,3) (H, S on each qutrit, SUM on each ordered pair), whose
distribution approaches the uniform one on Sp(6,3); the chi-square-style control is the fraction of M fixing z1,
expected 1/(3^6 - 1) = 1/728 under uniformity.
"""

from __future__ import annotations

import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11227_two_qutrit_cp_mixing as P2  # noqa: E402
import w33_pass11252_exact_reversibility as R  # noqa: E402

OUT = ROOT / "data" / "w33_pass11333_three_qutrit_one_gate.json"
_S = {}


def symplectic_part(wl, C):
    n2 = 2 * wl.n
    M = np.zeros((n2, n2), int)
    for j in range(n2):
        e = np.zeros(n2, int)
        e[j] = 1
        img = C @ wl.W[wl.index(e)] @ C.conj().T
        M[:, j] = wl.labels[int(np.argmax(np.abs(np.einsum('pij,ij->p', wl.W.conj(), img))))]
    return M


def _init():
    R.WEYL[3] = R.Weyl(3)
    _S["T1"] = R.local(3, 0, P2.T)


def _class(args):
    seed, full = args
    wl = R.WEYL[3]
    rng = np.random.default_rng(seed)
    C = R.random_clifford(3, rng, length=300)
    M = symplectic_part(wl, C)
    V = R.weil(wl, M, rng)
    frames = range(729) if full else [0] + [3 ** (5 - j) for j in range(6)]
    verdict = []
    for a in frames:
        U = wl.W[a] @ V @ _S["T1"]
        verdict.append(R.decide(U, 3, rng)[0] is False)
    fixes_z1 = bool(tuple(M[:, 1] % 3) == (0, 1, 0, 0, 0, 0))
    return dict(full=full, verdict=verdict, fixes_z1=fixes_z1)


def is_affine(points):
    keys = {tuple(p) for p in points}
    P = np.array(points)
    return all(tuple((x + y - z) % 3) in keys for x in P for y in P for z in P)


def _frame(args):
    seed, a = args
    wl = R.WEYL[3]
    rng = np.random.default_rng(seed)
    C = R.random_clifford(3, rng, length=300)
    M = symplectic_part(wl, C)
    V = R.weil(wl, M, np.random.default_rng(0))
    return R.decide(wl.W[a] @ V @ _S["T1"], 3, np.random.default_rng(a))[0] is False


def run(n_flag=1100, n_full_bad=2):
    res = dict(pass_id=11333)
    labels = np.array([[(a // 3 ** (5 - j)) % 3 for j in range(6)] for a in range(729)])
    seeds = [50000 + i for i in range(n_flag)]
    with Pool(11, initializer=_init) as pool:
        flag = pool.map(_class, [(s, False) for s in seeds], chunksize=10)
        bad_seeds = [s for s, r in zip(seeds, flag) if any(r["verdict"])]
        full = []
        for s in bad_seeds[:n_full_bad]:
            v = np.array(pool.map(_frame, [(s, a) for a in range(729)], chunksize=8))
            good = labels[~v]
            full.append(dict(seed=s, violating_frames=int(v.sum()),
                             nonviolating_affine=bool(len(good) == 0 or is_affine(good))))
    bad = len(bad_seeds)
    p = bad / n_flag
    res["flag_classes"] = n_flag
    res["flagged_bad"] = bad
    res["bad_fraction"] = p
    res["bad_fraction_stderr"] = float(np.sqrt(p * (1 - p) / n_flag))
    res["one_eighth_z_score"] = float((p - 0.125) / np.sqrt(0.125 * 0.875 / n_flag))
    res["fixes_z1_fraction"] = sum(r["fixes_z1"] for r in flag) / n_flag
    res["fixes_z1_expected_uniform"] = 1 / 728
    res["full_frame_checks_of_bad_classes"] = full
    res["condition"] = ("bad = flagged by frames 0, e_1..e_6; exact if the affine law holds at n = 3 (checked on all 729 "
                        "frames only for the classes listed in full_frame_checks_of_bad_classes)")
    print(res, flush=True)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
