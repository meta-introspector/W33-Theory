"""Pass 11531: why the cubed-phase test is exact on Clifford+T words but not on PU(3).

A relation (L, f) for U is: c(Lp) = mu omega^f(p) c(p) on all p (L anti-symplectic on F3^2, c the Weyl coefficients).
Cubing forgets whether f is affine; U is reversible iff some relation has f AFFINE on the support (Pass 11252).

(A) THE COUNTEREXAMPLE FAMILY.  For L = -swap and the quadratic f(a,b) = delta(b,2) + 2 delta(a,2) (and its double), the
    unitaries in the 6-dimensional eigenspace form a smooth family of real dimension 4 (3 after the global phase); every
    sampled point is non-reversible (Pass 11512 found the first one).
(B) THE WORDS.  Every one-qutrit Clifford+T operator of depth <= DMAX (distinct operators, Pass 11488's deduplication) is
    scanned for ALL relations (24 L's): which hold, whether f is affine, and -- for non-affine ones -- whether the operator is
    a Clifford or has T-count 1 (C1 T^(+-1) C2).  Result: non-affine relations occur ONLY for Cliffords and T-count-1
    operators, at full support; both kinds are reversible.  Operators of T-count >= 2 never satisfy a non-affine relation,
    so the discarded condition never matters on them -- that is why the test decides every deep word.
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11350_linear_decider as L  # noqa: E402
import w33_pass11422_depth5_exact as E  # noqa: E402
import w33_pass11436_depth_reduction as RD  # noqa: E402
import w33_pass11500_exact_depth9 as P9  # noqa: E402
import w33_pass11512_cubed_phase_exhaustive as C  # noqa: E402

OUT = ROOT / "data" / "w33_pass11531_cubed_test_words_vs_pu3.json"
DMAX = 7


def family(samples=300, seed=11531):
    E._init()
    anti = C.mats(2)
    w = np.exp(2j * np.pi / 3)
    rng = np.random.default_rng(seed)
    out = {}
    for f in ([0, 0, 1, 0, 0, 1, 2, 2, 0], [0, 0, 2, 0, 0, 2, 1, 1, 0]):
        Lp = C.perm_of(anti[0])
        Pm = np.zeros((9, 9), complex)
        for p in range(9):
            Pm[Lp[p], p] = w ** f[p]
        vals, vecs = np.linalg.eig(Pm)
        for mu in np.unique(np.round(vals, 8)):
            Q, _ = np.linalg.qr(vecs[:, np.abs(vals - mu) < 1e-6])
            d = Q.shape[1]
            if d != 6:
                continue

            def res(x):
                U = np.einsum('p,pij->ij', Q @ (x[:d] + 1j * x[d:]), P9.WP)
                R = U @ U.conj().T - np.eye(3)
                return np.concatenate([R.real.ravel(), R.imag.ravel()])

            dims, rev, ovs = Counter(), Counter(), []
            for _ in range(samples):
                r = least_squares(res, rng.normal(size=2 * d), method='lm')
                if np.abs(r.fun).max() > 1e-11:
                    continue
                s = np.linalg.svd(r.jac, compute_uv=False)
                dims[2 * d - int((s > 1e-8 * s[0]).sum())] += 1
                U = np.einsum('p,pij->ij', Q @ (r.x[:d] + 1j * r.x[d:]), P9.WP)
                ov = float(P9.X.overlaps(U[None])[0])
                ovs.append(ov)
                rev["reversible" if ov > 3 - 1e-7 else "not reversible"] += 1
            out[str(f)] = dict(L="-swap (index 0)", eigenspace_dim=d, mu=[float(mu.real), float(mu.imag)],
                               local_dimension_counts=dict(dims), reversibility=dict(rev),
                               overlap_min=min(ovs), overlap_max=max(ovs))
            print(f, out[str(f)], flush=True)
    return out


def words():
    E._init()
    RD._init()
    CT, RT, CF = RD._S["CT"], RD._S["RT"], E._S["CF"]
    labs = np.array([(a, b) for a in range(3) for b in range(3)])
    Aff = np.concatenate([np.ones((9, 1), np.int64), labs], 1) % 3
    T = E._S["T"]
    Tinv = [np.linalg.inv(np.linalg.matrix_power(T, j)) for j in (1, 2)]
    U = CT.copy()
    w = np.ones(len(U), np.int64)
    out = {}
    for k in range(1, DMAX + 1):
        if k > 1:
            U = np.einsum('rij,wjk->wrik', RT, U).reshape(-1, 3, 3)
            w = np.repeat(w, 8)
            U, w = P9.dedupe(U, w)
        st = Counter()
        for i in range(0, len(U), 50000):
            Uc = U[i:i + 50000]
            c = np.einsum('pij,wji->wp', P9.WP.conj(), Uc) / 3
            mod = np.abs(c)
            okm = (np.abs(mod[:, P9.PERM] - mod[:, None, :]) < 1e-9).all(axis=2)
            wi, li = np.nonzero(okm)
            if len(wi) == 0:
                continue
            cc = c[wi]
            supp = mod[wi] > 1e-9
            ratio = cc[np.arange(len(wi))[:, None], P9.PERM[li]] / np.where(supp, cc, 1)
            rel = ratio / ratio[np.arange(len(wi)), np.argmax(supp, axis=1)][:, None]
            ph = np.angle(rel) * 3 / (2 * np.pi)
            ok = np.where(supp, np.abs(ph - np.rint(ph)), 0).max(axis=1) < 1e-6
            for j in np.flatnonzero(ok):
                f = np.rint(ph[j]).astype(int) % 3
                sp_ = supp[j]
                if L.solve_affine(Aff[sp_], f[sp_]) is not None:
                    st["affine relations"] += 1
                    continue
                Uw = Uc[wi[j]]
                if np.abs(np.abs(np.einsum('cij,ij->c', CF.conj(), Uw)) - 3).min() < 1e-8:
                    st["non-affine: Clifford"] += 1
                    continue
                tc1 = any(np.abs(np.einsum('dij,cij->cd', CF.conj(),
                                           np.einsum('ij,cjk,kl->cil', Uw, CF.conj().transpose(0, 2, 1), Ti))).max() > 3 - 1e-8
                          for Ti in Tinv)
                st["non-affine: T-count 1" if tc1 else "non-affine: T-count >= 2 (would break the mechanism)"] += 1
                st[f"non-affine support {int(sp_.sum())}"] += 1
        out[str(k)] = dict(distinct_operators=len(U), **st)
        print(k, out[str(k)], flush=True)
    return out


def main():
    res = dict(pass_id=11531)
    res["counterexample_family"] = family()
    res["word_relations"] = words()
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
