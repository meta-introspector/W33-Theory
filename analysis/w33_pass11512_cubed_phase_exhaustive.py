"""Pass 11512: the cubed-phase test on PU(3), exhaustively up to symmetry.

Pass 11504 searched 4000 random (L, f) and found only 4 unitaries on non-affine phase relations, all reversible.  Here
EVERY pair is covered.  Two symmetries reduce the 24 x 3^9 pairs (L anti-symplectic on F3^2, f : F3^2 -> Z3):
  * f matters modulo affine functions (an affine part is absorbed into mu and the Weyl translation b);
  * Clifford conjugation U -> V_S U V_S^dag sends c(p) -> c(S^-1 p), hence (L, f) -> (S L S^-1, f o S^-1), S in SL(2,3).
f is canonicalised by subtracting the affine function through its values at 0, e1, e2.  For one representative of each
orbit with f non-affine, every eigenspace of the monomial matrix P_(L,f) is searched for unitaries
U = sum_p c(p) W(p) (BFGS on ||U U^dag - 1||^2 from many starts), and each unitary found is tested with Theorem 1.
"""

from __future__ import annotations

import itertools
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11422_depth5_exact as E  # noqa: E402
import w33_pass11500_exact_depth9 as P9  # noqa: E402

OUT = ROOT / "data" / "w33_pass11512_cubed_phase_exhaustive.json"
LABS = [(a, b) for a in range(3) for b in range(3)]
IDX = {p: i for i, p in enumerate(LABS)}


def mats(det):
    return [np.array([[a, b], [c, d]]) for a, b, c, d in itertools.product(range(3), repeat=4) if (a * d - b * c) % 3 == det]


def perm_of(S):
    return tuple(IDX[tuple(int(x) for x in (S @ np.array(p)) % 3)] for p in LABS)


def canon_f(f):
    """f minus the affine function agreeing with it at 0, e1 = (1,0), e2 = (0,1)"""
    f0, f1, f2 = f[IDX[(0, 0)]], f[IDX[(1, 0)]], f[IDX[(0, 1)]]
    return tuple(int((f[i] - (f0 + (f1 - f0) * a + (f2 - f0) * b)) % 3) for i, (a, b) in enumerate(LABS))


def orbit_reps():
    SL = mats(1)
    anti = mats(2)
    Lp = [perm_of(A) for A in anti]
    Sp = [perm_of(S) for S in SL]
    Sinvp = [perm_of(np.round(np.linalg.inv(S) * 1).astype(int) % 3 if False else
                     np.array([[S[1, 1], -S[0, 1]], [-S[1, 0], S[0, 0]]]) % 3) for S in SL]
    seen = set()
    reps = []
    for li, L in enumerate(anti):
        for vals in itertools.product(range(3), repeat=9):
            f = canon_f(vals)
            if not any(f):
                continue                                                         # affine
            key = (Lp[li], f)
            if key in seen:
                continue
            reps.append((li, f))
            for S, Sinv in zip(SL, Sinvp):                                       # mark the whole orbit
                L2 = perm_of((S @ L @ np.array([[S[1, 1], -S[0, 1]], [-S[1, 0], S[0, 0]]])) % 3)
                f2 = canon_f([f[Sinv[i]] for i in range(9)])                      # f o S^-1
                seen.add((L2, f2))
    return anti, reps


def search(anti, li, f, rng, starts=12):
    w = np.exp(2j * np.pi / 3)
    Lperm = perm_of(anti[li])
    Pm = np.zeros((9, 9), complex)
    for p in range(9):
        Pm[Lperm[p], p] = w ** f[p]
    vals, vecs = np.linalg.eig(Pm)
    found = []
    nspaces = 0
    for mu in np.unique(np.round(vals, 8)):
        Q, _ = np.linalg.qr(vecs[:, np.abs(vals - mu) < 1e-6])
        d = Q.shape[1]
        nspaces += 1

        def loss(x):
            U = np.einsum('p,pij->ij', Q @ (x[:d] + 1j * x[d:]), P9.WP)
            return np.linalg.norm(U @ U.conj().T - np.eye(3)) ** 2

        for _ in range(starts):
            r = minimize(loss, rng.normal(size=2 * d), method='BFGS')
            if r.fun < 1e-20:
                U = np.einsum('p,pij->ij', Q @ (r.x[:d] + 1j * r.x[d:]), P9.WP)
                found.append(float(P9.X.overlaps(U[None])[0]))
    return nspaces, found


def strong(starts=200, seed=115120):
    """second, stronger pass over the same 816 orbits: dimension-1 eigenspaces decided EXACTLY (U is fixed up to scale, so
    it is unitary iff U U^dag is a multiple of 1); higher-dimensional ones by Levenberg-Marquardt from `starts` starts,
    recording how many starts hit each unitary-containing eigenspace (the miss rate)"""
    from scipy.optimize import least_squares
    E._init()
    anti, reps = orbit_reps()
    rng = np.random.default_rng(seed)
    w = np.exp(2j * np.pi / 3)
    st = Counter()
    hit_rates = []
    found = []
    for li, f in reps:
        Lperm = perm_of(anti[li])
        Pm = np.zeros((9, 9), complex)
        for p in range(9):
            Pm[Lperm[p], p] = w ** f[p]
        vals, vecs = np.linalg.eig(Pm)
        for mu in np.unique(np.round(vals, 8)):
            Q, _ = np.linalg.qr(vecs[:, np.abs(vals - mu) < 1e-6])
            d = Q.shape[1]
            st[f"eigenspaces of dim {d}"] += 1
            Us = []
            if d == 1:
                U0 = np.einsum('p,pij->ij', Q[:, 0], P9.WP)
                G = U0 @ U0.conj().T
                lam = np.trace(G).real / 3
                if np.abs(G - lam * np.eye(3)).max() < 1e-10 * max(lam, 1e-300):
                    Us.append(U0 / np.sqrt(lam))
            else:
                def res(x):
                    U = np.einsum('p,pij->ij', Q @ (x[:d] + 1j * x[d:]), P9.WP)
                    R = U @ U.conj().T - np.eye(3)
                    return np.concatenate([R.real.ravel(), R.imag.ravel()])
                hits = 0
                for _ in range(starts):
                    r = least_squares(res, rng.normal(size=2 * d), method='lm')
                    if np.abs(r.fun).max() < 1e-11:
                        hits += 1
                        if not Us:
                            Us.append(np.einsum('p,pij->ij', Q @ (r.x[:d] + 1j * r.x[d:]), P9.WP))
                if hits:
                    hit_rates.append(hits / starts)
            for U in Us:
                ov = float(P9.X.overlaps(U[None])[0])
                st["unitary-containing eigenspaces"] += 1
                st["... reversible" if ov > 3 - 1e-7 else "... NOT reversible"] += 1
                found.append(dict(L=li, f=list(f), dim=d, overlap=ov))
    out = dict(starts=starts, counts=dict(st), hit_rates=hit_rates,
               min_hit_rate=min(hit_rates) if hit_rates else None, unitaries=found)
    print(json.dumps(out, indent=1), flush=True)
    return out


def counterexample(seed=115121):
    """an explicit non-reversible unitary passing the cubed-phase test: rebuild from the strong pass's first non-reversible
    hit, verify unitarity, the prefilter, Theorem 1's overlap < 3 and Pass 11252's independent decider"""
    from scipy.optimize import least_squares
    import w33_pass11252_exact_reversibility as R
    E._init()
    R.WEYL[1] = R.Weyl(1)
    d = json.load(open(OUT))["strong"]
    bad = [u for u in d["unitaries"] if u["overlap"] <= 3 - 1e-7]
    anti = mats(2)
    w = np.exp(2j * np.pi / 3)
    rng = np.random.default_rng(seed)
    out = []
    for b in bad:
        Lperm = perm_of(anti[b["L"]])
        Pm = np.zeros((9, 9), complex)
        for p in range(9):
            Pm[Lperm[p], p] = w ** b["f"][p]
        vals, vecs = np.linalg.eig(Pm)
        for mu in np.unique(np.round(vals, 8)):
            Q, _ = np.linalg.qr(vecs[:, np.abs(vals - mu) < 1e-6])
            dd = Q.shape[1]
            if dd != b["dim"]:
                continue

            def res(x):
                U = np.einsum('p,pij->ij', Q @ (x[:dd] + 1j * x[dd:]), P9.WP)
                Rm = U @ U.conj().T - np.eye(3)
                return np.concatenate([Rm.real.ravel(), Rm.imag.ravel()])

            for _ in range(100):
                r = least_squares(res, rng.normal(size=2 * dd), method='lm')
                if np.abs(r.fun).max() < 1e-11:
                    U = np.einsum('p,pij->ij', Q @ (r.x[:dd] + 1j * r.x[dd:]), P9.WP)
                    ov = float(P9.X.overlaps(U[None])[0])
                    if ov > 3 - 1e-7:
                        continue
                    try:
                        verdict = R.decide(U, 1, rng)[0]
                    except AssertionError:
                        verdict = "ambiguous"
                    c = np.einsum('pij,ji->p', P9.WP.conj(), U) / 3
                    out.append(dict(L=b["L"], f=b["f"], eigenspace_dim=dd, mu=[float(mu.real), float(mu.imag)],
                                    U_real=U.real.round(15).tolist(), U_imag=U.imag.round(15).tolist(),
                                    unitarity_error=float(np.abs(U @ U.conj().T - np.eye(3)).max()),
                                    prefilter_passes=bool(P9.prefilter(U[None])[0]), max_overlap=ov,
                                    pass11252_decider=verdict, weyl_moduli=np.abs(c).round(12).tolist()))
                    break
            break
        if out:
            break
    print(json.dumps(out, indent=1), flush=True)
    return out


def run(starts=12, seed=11512):
    E._init()
    anti, reps = orbit_reps()
    rng = np.random.default_rng(seed)
    st = Counter()
    st["orbit representatives (L, non-affine f mod affine)"] = len(reps)
    ce = []
    for li, f in reps:
        n, found = search(anti, li, f, rng, starts)
        st["eigenspaces"] += n
        st["unitaries found (with multiplicity over starts)"] += len(found)
        st["orbits containing a unitary"] += bool(found)
        bad = [ov for ov in found if ov <= 3 - 1e-7]
        st["non-reversible unitaries"] += len(bad)
        if bad and len(ce) < 5:
            ce.append(dict(L=li, f=list(f), overlaps=bad[:3]))
    res = dict(pass_id=11512, starts_per_eigenspace=starts, counts=dict(st), counterexamples=ce)
    print(json.dumps(res, indent=1), flush=True)
    return res


def main():
    if "counter" in sys.argv:
        res = json.load(open(OUT))
        res["counterexample"] = counterexample()
        json.dump(res, open(OUT, "w"), indent=1)
        return
    if "strong" in sys.argv:
        res = json.load(open(OUT))
        res["strong"] = strong()
        json.dump(res, open(OUT, "w"), indent=1)
        return
    res = run()
    json.dump(res, open(OUT, "w"), indent=1)


if __name__ == "__main__":
    main()
