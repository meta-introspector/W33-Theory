"""Pass 11421: four qutrits, one magic gate -- do the universal magic-axis rules of Pass 11373 hold at n = 4, and where
does the bad-class fraction go (1/8, 1/8, 437/3276 for n = 1, 2, 3)?

The F3-linear decider of Pass 11350 is exact for any n, but its Weyl tables (3^(2n) operators of size 3^n) cost 690 MB
at n = 4.  SparseDecider replaces every use of the table by the closed form
        W(p)|j> = omega^(2 x.z + z.j) |j + x>        (p = (x1, z1, ..., xn, zn), Pass 11252's convention),
so V_M = sum_p W(Mp)|j><0|W(p)^dag is assembled from 3^(2n) phased rank-one entries and frame_of reads x from one
column and z from n phase ratios.  It is validated class by class against the table-based decider at n = 2 and n = 3
(identical good-frame sets).

Sampling at n = 4 (Sp(8,3) has ~6.3e17 elements; no orbit census):
  * uniform classes (long random words in H, S, SUM): bad-class fraction and per-cell fractions;
  * targeted classes in the four rule cells (fixed / reversed / M^2 z1 = +-z1 on a line through z1), built with
    Pass 11352's transvection map;
  * classes whose (S) solution space exceeds the cap are reported as undecided, never guessed.
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
import w33_pass11227_two_qutrit_cp_mixing as P2  # noqa: E402
import w33_pass11252_exact_reversibility as R  # noqa: E402
import w33_pass11350_linear_decider as L  # noqa: E402
import w33_pass11352_three_qutrit_geometry as GEO  # noqa: E402

OUT = ROOT / "data" / "w33_pass11421_four_qutrits.json"
OM = np.exp(2j * np.pi / 3)


class LightWeyl:
    """labels, index and forms of R.Weyl without the operator table"""

    def __init__(self, n):
        import itertools
        self.n, self.D = n, 3 ** n
        self.labels = np.array(list(itertools.product(range(3), repeat=2 * n)))
        self.pow3 = 3 ** np.arange(2 * n - 1, -1, -1)
        om = np.zeros((2 * n, 2 * n), dtype=int)
        for i in range(n):
            om[2 * i, 2 * i + 1] = 1
            om[2 * i + 1, 2 * i] = -1
        self.Om = om
        self.J = np.diag([1, -1] * n) % 3
        self.digits = np.array(list(itertools.product(range(3), repeat=n)))      # basis |j>, qutrit 1 most significant
        self.dpow = 3 ** np.arange(n - 1, -1, -1)

    def index(self, v):
        return int((np.asarray(v) % 3) @ self.pow3)

    def act(self, p, j):
        """W(p)|j> = phase |k>: returns (phase, k) for label array p (..., 2n) and digit array j (..., n)"""
        x, z = p[..., 0::2], p[..., 1::2]
        ph = (2 * np.sum(x * z, axis=-1) + np.sum(z * j, axis=-1)) % 3
        k = ((j + x) % 3) @ self.dpow
        return OM ** ph, k

    def matrix(self, p):
        p = np.asarray(p) % 3
        ph, k = self.act(np.broadcast_to(p, (self.D, 2 * self.n)), self.digits)
        Wm = np.zeros((self.D, self.D), complex)
        Wm[k, np.arange(self.D)] = ph
        return Wm


class SparseDecider(L.Decider):
    def __init__(self, n, rng=None):                       # noqa: D401 -- deliberately not calling the table builder
        self.n = n
        self.wl = LightWeyl(n)
        self.rng = rng or np.random.default_rng(0)
        self.N2 = 2 * n
        self.T1 = R.local(n, 0, P2.T)
        self.J = np.diag([1, 2] * n)
        self.z1 = np.zeros(self.N2, dtype=np.int64)
        self.z1[1] = 1
        s = np.eye(self.N2, dtype=np.int64)
        s[1, 0] = 1
        self.s_pow = [np.linalg.matrix_power(s, k) % 3 for k in range(3)]
        self.sinv_pow = [self.s_pow[(-k) % 3] for k in range(3)]
        self.g = {}
        self._weil_cache = {}

    def weil(self, M):
        key = tuple(M.ravel())
        if key in self._weil_cache:
            return self._weil_cache[key]
        wl = self.wl
        P = wl.labels.astype(np.int64)
        MP = (P @ M.T) % 3
        zero = np.zeros((len(P), self.n), dtype=np.int64)
        ph0, x0 = wl.act(P, zero)                               # W(p)|0> = ph0 |x0>
        V = None
        for jidx in range(wl.D):
            j = np.broadcast_to(wl.digits[jidx], (len(P), self.n))
            ph, k = wl.act(MP, j)                               # W(Mp)|j> = ph |k>
            V0 = np.zeros((wl.D, wl.D), complex)
            np.add.at(V0, (k, x0), ph * np.conj(ph0))
            nrm = np.linalg.norm(V0)
            if nrm > 1e-6:
                V = V0 / nrm * np.sqrt(wl.D)
                break
        for i in range(self.N2):
            e = np.zeros(self.N2, dtype=np.int64)
            e[i] = 1
            assert np.allclose(V @ wl.matrix(e) @ V.conj().T, wl.matrix(M @ e % 3), atol=1e-8)
        self._weil_cache[key] = V
        return V

    def frame_of(self, U):
        wl = self.wl
        col = U[:, 0]
        k0 = int(np.argmax(np.abs(col)))
        x = wl.digits[k0]
        z = np.zeros(self.n, dtype=np.int64)
        for i in range(self.n):
            j = np.zeros(self.n, dtype=np.int64)
            j[i] = 1
            src = int(j @ wl.dpow)
            dst = int(((j + x) % 3) @ wl.dpow)
            r = U[dst, src] / col[k0]                           # = omega^{z_i} (the 2x.z phase cancels)
            z[i] = int(np.rint(np.angle(r) / (2 * np.pi / 3))) % 3
        lab = np.empty(self.N2, dtype=np.int64)
        lab[0::2], lab[1::2] = x, z
        Wq = wl.matrix(lab)
        ov = abs(np.vdot(Wq, U)) / wl.D
        assert np.isclose(ov, 1), "not a Weyl operator"
        return lab

    def gframe(self, e):
        key = tuple(e)
        if key not in self.g:
            k = int(e[0])
            X = self.T1 @ self.wl.matrix(e) @ self.T1.conj().T @ self.weil(self.sinv_pow[k]).conj().T
            self.g[key] = self.frame_of(X)
        return self.g[key]


def validate(n, classes=40, seed=0):
    rng = np.random.default_rng(seed)
    L.CAP = 3 ** 11
    A, B = L.Decider(n), SparseDecider(n)
    gens = GEO.gen_mats(n)
    agree = checked = 0
    for _ in range(classes):
        M = GEO.random_symplectic(rng, gens)
        ga, gb = A.good_frames(M), B.good_frames(M)
        if ga is None or gb is None:
            continue
        checked += 1
        agree += bool((ga == gb).all())
    return dict(n=n, classes_checked=checked, identical_good_frame_sets=agree)


_S = {}


def _init():
    import ast
    L.CAP = 3 ** 11
    _S["D"] = SparseDecider(4)
    _S["gens"] = GEO.gen_mats(4)
    d2 = L.Decider(2)
    reps = {"M^2 z1 = z1": [], "M^2 z1 = -z1": []}
    for line in (ROOT / "data" / "w33_pass11373_orbits_n2.txt").read_text().splitlines():
        m, _ = line.rsplit(";", 1)
        M = np.array(ast.literal_eval(m), dtype=np.int64).T % 3
        c = GEO.cell(M, d2.z1, d2.wl.Om)
        for key in reps:
            if c == "same line, " + key:
                reps[key].append(M)
    _S["line_reps"] = reps


def _target(M, target, rng):
    """compose M with a symplectic map so that the image of z1 is target(z1) (the rule cells)"""
    D = _S["D"]
    z = D.z1
    w = (target @ z) % 3 if target is not None else None
    return (GEO.map_to((M @ z) % 3, w, D.wl.Om, rng) @ M) % 3


def _job(args):
    seed, kind = args
    D, gens = _S["D"], _S["gens"]
    rng = np.random.default_rng(seed)
    M = GEO.random_symplectic(rng, gens, length=300)
    z, Om = D.z1, D.wl.Om
    if kind == "fixed":
        M = (GEO.map_to((M @ z) % 3, z, Om, rng) @ M) % 3
    elif kind == "reversed":
        M = (GEO.map_to((M @ z) % 3, (-z) % 3, Om, rng) @ M) % 3
    elif kind in ("line, M^2z1=z1", "line, M^2z1=-z1"):
        # an exact n = 2 class of that cell (Pass 11373 orbit list) embedded as A2 (+) B, B random in Sp(4,3) on
        # qutrits 3-4, then conjugated by a random N in Stab(z1) (preserves the cell, mixes all qutrits)
        reps = _S["line_reps"]["M^2 z1 = z1" if kind.endswith("=z1") else "M^2 z1 = -z1"]
        A = np.eye(D.N2, dtype=np.int64)
        A[:4, :4] = reps[int(rng.integers(len(reps)))]
        A[4:, 4:] = GEO.random_symplectic(rng, GEO.gen_mats(2), length=100)
        h = GEO.random_symplectic(rng, gens, length=300)
        N = (GEO.map_to((h @ z) % 3, z, Om, rng) @ h) % 3
        M = (N @ A @ R._inv_mod3(N)) % 3
    cell = GEO.cell(M, z, Om)
    g = D.good_frames(M)
    if g is None:
        return dict(kind=kind, cell=cell, undecided=True)
    return dict(kind=kind, cell=cell, bad=bool((~g).any()), bad_frames=int((~g).sum()))


def run(uniform=3000, targeted=60, nproc=4):
    res = dict(pass_id=11421)
    res["validation"] = [validate(2, 60), validate(3, 25)]
    print(res["validation"], flush=True)
    assert all(v["classes_checked"] == v["identical_good_frame_sets"] for v in res["validation"])
    jobs = [(114210000 + i, "uniform") for i in range(uniform)]
    for j, kind in enumerate(("fixed", "reversed", "line, M^2z1=z1", "line, M^2z1=-z1")):
        jobs += [(114215000 + 100 * j + i, kind) for i in range(targeted)]
    with Pool(nproc, initializer=_init) as pool:
        rows = pool.map(_job, jobs, chunksize=1)
    uni = [r for r in rows if r["kind"] == "uniform"]
    dec = [r for r in uni if "bad" in r]
    p = sum(r["bad"] for r in dec) / max(1, len(dec))
    cells = defaultdict(Counter)
    for r in dec:
        cells[r["cell"]]["bad" if r["bad"] else "good"] += 1
    res["uniform"] = dict(sampled=len(uni), decided=len(dec), undecided=len(uni) - len(dec), bad_class_fraction=p,
                          stderr=float(np.sqrt(p * (1 - p) / max(1, len(dec)))), cells={k: dict(v) for k, v in cells.items()},
                          bad_frames_histogram=dict(Counter(r["bad_frames"] for r in dec)))
    tg = defaultdict(Counter)
    for r in rows:
        if r["kind"] != "uniform":
            tg[r["kind"]]["skipped" if r.get("skipped") else ("undecided" if r.get("undecided") else
                                                               ("bad" if r["bad"] else "good"))] += 1
            if "cell" in r:
                tg[r["kind"]]["cell=" + r["cell"]] += 1
    res["targeted_rule_cells"] = {k: dict(v) for k, v in tg.items()}
    print(json.dumps({k: res[k] for k in ("uniform", "targeted_rule_cells")}, indent=1), flush=True)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
