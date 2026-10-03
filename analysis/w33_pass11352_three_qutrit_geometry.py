"""Pass 11352: one magic gate on three qutrits with the exact linear decider -- the bad fraction and the W(3,3)-type
geometry of the magic axis.

Pass 11333 sampled 1100 classes with the slow Weyl criterion (0.111 +- 0.009).  Pass 11350's F3-linear decider is exact
(validated on all 51,840 two-qutrit classes and on the two exhaustive three-qutrit classes of Pass 11333) and needs one
Weil unitary per candidate symplectic solution instead of 729 Weyl decisions.  Here:
  (a) a uniform sample of Sp(6,3) classes (long random words in the H, S, SUM generators) gives the bad fraction and its
      split over the geometric cells of the magic axis z1 (fixed / reversed / collinear on different lines / same line /
      non-collinear), as at n = 2 (Pass 11331);
  (b) targeted samples of the rare cells 'Mz1 = z1' and 'Mz1 = -z1' (M composed with a random symplectic g carrying
      M z1 back to +-z1) test the n = 2 rules 'fixed => bad' and 'reversed => good'.
Undecided classes (solution space of (S) beyond the cap) are reported, never guessed.
"""

from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11252_exact_reversibility as R  # noqa: E402
import w33_pass11350_linear_decider as L  # noqa: E402

OUT = ROOT / "data" / "w33_pass11352_three_qutrit_geometry.json"
_S = {}


def gen_mats(n):
    """symplectic matrices (on (x1,z1,...,xn,zn) column vectors) of H, S on each qutrit and SUM on each ordered pair"""
    N2 = 2 * n
    out = []
    for q in range(n):
        H = np.eye(N2, dtype=np.int64)
        H[2 * q:2 * q + 2, 2 * q:2 * q + 2] = [[0, 2], [1, 0]]
        S = np.eye(N2, dtype=np.int64)
        S[2 * q + 1, 2 * q] = 1
        out += [H, S]
    for c in range(n):
        for t in range(n):
            if c != t:
                M = np.eye(N2, dtype=np.int64)
                M[2 * t, 2 * c] = 1          # x_t += x_c
                M[2 * c + 1, 2 * t + 1] = 2  # z_c -= z_t
                out.append(M)
    return out


def random_symplectic(rng, gens, length=200):
    M = np.eye(gens[0].shape[0], dtype=np.int64)
    for g in rng.integers(len(gens), size=length):
        M = (gens[g] @ M) % 3
    return M


def transvection(u, Om):
    """tau_u(x) = x + omega(x, u) u  (symplectic)"""
    N2 = len(u)
    return (np.eye(N2, dtype=np.int64) + np.outer(u, u @ Om.T)) % 3


def map_to(v, w, Om, rng):
    """a product of transvections sending v to w (both nonzero) -- random, then composed with a random element of
    Stab(w) so the result is spread over the coset"""
    def om(a, b):
        return int(a @ Om @ b) % 3
    N2 = len(v)
    if (v == w).all():
        return np.eye(N2, dtype=np.int64)
    if om(v, w) == 1:
        for c in (1, 2):
            h = transvection((c * (w - v)) % 3, Om)
            if ((h @ v - w) % 3 == 0).all():
                return h
    while True:
        y = rng.integers(3, size=N2)
        if om(v, y) == 1 and om(y, w) == 1:
            h1 = map_to(v, y, Om, rng)
            h2 = map_to(y, w, Om, rng)
            return (h2 @ h1) % 3


def cell(M, z, Om):
    def om(u, v):
        return int(u @ Om @ v) % 3
    a = M @ z % 3
    b = R._inv_mod3(M) @ z % 3
    if (a == z).all():
        return "Mz1 = z1"
    if ((a + z) % 3 == 0).all():
        return "Mz1 = -z1"
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


def _init():
    L.CAP = 3 ** 11
    _S["D"] = L.Decider(3)
    _S["gens"] = gen_mats(3)


def _job(args):
    seed, target = args
    D = _S["D"]
    rng = np.random.default_rng(seed)
    gens = _S["gens"]
    M = random_symplectic(rng, gens)
    z = D.z1
    if target:
        want = (z if target == 1 else (2 * z)) % 3
        h = map_to(M @ z % 3, want, D.wl.Om, rng)
        M = (h @ M) % 3
    assert ((M.T @ D.wl.Om @ M - D.wl.Om) % 3 == 0).all()
    g = D.good_frames(M)
    c = cell(M, z, D.wl.Om)
    if g is None:
        return dict(cell=c, verdict="undecided")
    nb = int((~g).sum())
    return dict(cell=c, verdict="bad" if nb else "good", violating_frames=nb)


def summarize(rows):
    tab = defaultdict(Counter)
    for r in rows:
        if r is None:
            continue
        tab[r["cell"]][r["verdict"]] += 1
        if r["verdict"] == "bad":
            tab[r["cell"]][f"frames={r['violating_frames']}"] += 1
    return {k: dict(v) for k, v in sorted(tab.items())}


def run(n_uniform=8000, n_target=150):
    res = dict(pass_id=11352)
    with Pool(10, initializer=_init) as pool:
        uni = pool.map(_job, [(500000 + i, 0) for i in range(n_uniform)], chunksize=20)
        fix = pool.map(_job, [(700000 + i, 1) for i in range(n_target)], chunksize=5)
        rev = pool.map(_job, [(800000 + i, 2) for i in range(n_target)], chunksize=5)
    dec = [r for r in uni if r and r["verdict"] != "undecided"]
    bad = sum(r["verdict"] == "bad" for r in dec)
    p = bad / len(dec)
    res["uniform"] = dict(samples=n_uniform, decided=len(dec), undecided=n_uniform - len(dec), bad=bad, bad_fraction=p,
                          stderr=float(np.sqrt(p * (1 - p) / len(dec))),
                          z_vs_one_eighth=float((p - 1 / 8) / np.sqrt(p * (1 - p) / len(dec))),
                          z_vs_one_ninth=float((p - 1 / 9) / np.sqrt(p * (1 - p) / len(dec))),
                          cells=summarize(uni))
    res["targeted_fixed"] = summarize(fix)
    res["targeted_reversed"] = summarize(rev)
    print(json.dumps(res, indent=1), flush=True)
    return res


def _validate(seed):
    """independent check of the linear decider at n = 3: Weyl-criterion verdicts (Pass 11252) on the 7 screening frames
    and on 4 random frames, against the decider's good-frame set"""
    D = _S["D"]
    rng = np.random.default_rng(seed)
    M = random_symplectic(rng, _S["gens"])
    g = D.good_frames(M)
    if g is None:
        return None
    V = D.weil(M)
    frames = [0, 243, 81, 27, 9, 3, 1] + [int(x) for x in rng.integers(729, size=4)]
    agree = 0
    for a in frames:
        U = D.wl.W[a] @ V @ D.T1
        weyl_rev = R.decide(U, 3, np.random.default_rng(a))[0] is True
        agree += weyl_rev == bool(g[a])
    return dict(bad=bool((~g).any()), agree=agree, frames=len(frames))


def validate(n=60):
    with Pool(10, initializer=_init) as pool:
        rows = pool.map(_validate, [500000 + i for i in range(n)], chunksize=1)
    rows = [r for r in rows if r]
    return dict(classes=len(rows), bad_classes=sum(r["bad"] for r in rows),
                frame_checks=sum(r["frames"] for r in rows), agree=sum(r["agree"] for r in rows))


def extend(n_more=24000):
    """a second, independent uniform batch (seeds 600000+) merged with the first"""
    with Pool(10, initializer=_init) as pool:
        more = pool.map(_job, [(600000 + i, 0) for i in range(n_more)], chunksize=40)
    dec = [r for r in more if r and r["verdict"] != "undecided"]
    return dict(samples=n_more, decided=len(dec), bad=sum(r["verdict"] == "bad" for r in dec), cells=summarize(more))


def main():
    if "--extend-only" in sys.argv:
        res = json.load(open(OUT))
        ext = extend()
        res["uniform_batch2"] = ext
        b = res["uniform"]["bad"] + ext["bad"]
        n = res["uniform"]["decided"] + ext["decided"]
        p = b / n
        se = float(np.sqrt(p * (1 - p) / n))
        res["uniform_combined"] = dict(decided=n, bad=b, bad_fraction=p, stderr=se,
                                       z_vs_one_eighth=(p - 1 / 8) / se, z_vs_one_ninth=(p - 1 / 9) / se)
        json.dump(res, open(OUT, "w"), indent=1, default=str)
        print(res["uniform_combined"], ext["cells"])
        return
    if "--validate-only" in sys.argv:
        res = json.load(open(OUT))
        res["weyl_cross_check"] = validate()
        json.dump(res, open(OUT, "w"), indent=1, default=str)
        print(res["weyl_cross_check"])
        return
    res = run()
    res["weyl_cross_check"] = validate()
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
