"""Pass 11487: the magic-axis law when the magic axis lies inside Im(M - I) -- the case Pass 11457's theorem could not touch.

Setting as in Pass 11457: U = W(a) V_M T1, z1 the magic axis, M^2 z1 = z1, v = z1 + M z1 != 0, R_M = Im(M - I),
Rad = R_M cap R_M^perp.  If z1 is in R_M then v is in Rad, so Pass 11457's condition (a orthogonal to Rad) excludes every
F' frame and the theorem says nothing.

LEMMA (level hyperplanes).  If z1 is in R_M, m_G = |c_G| for G = V_M T1 is supported on R_M and equals
|R_M|^(-1/2) |sum_j tau_j omega^(j s(p) + q j^2)| with s affine on R_M.  For q != 0 the three values are distinct, so the
level sets of m_G are the three cosets of a hyperplane K of R_M (K = ker of the linear part of s).  [checked per class:
support = R_M and the level sets are three cosets of one subspace of index 3]

THEOREM (extension).  In that situation every frame a with omega(a, v) != 0 and omega(a, Rad cap K) = 0 is VIOLATING.
PROOF.  m_U(p) = m_G(p - a).  If U is reversible, Pass 11457 Lemma 2 gives an anti-symplectic L preserving m_U; the three
level sets carry distinct values, so L fixes each of a + X_k (X_k the cosets of K).  Then L(a + R_M) = a + R_M, so
L R_M = R_M and t := L a - a is in R_M; and L(a + r) - a = t + L r lies in the coset of r, so L r - r is in K (put r = 0:
t in K).  L anti-symplectic preserves R_M, hence R_M^perp and Rad; v in Rad, so r0 := L v - v is in Rad cap K.  Then
    -omega(a, v) = omega(L a, L v) = omega(a + t, v + r0) = omega(a, v) + omega(a, r0)        (t in R_M; v, r0 in R_M^perp)
and omega(a, r0) = 0 by hypothesis, so 2 omega(a, v) = 0: contradiction.  QED.
In particular, if Rad cap K = 0 (e.g. Rad = span(v) and v not in K) the class satisfies F' completely.
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11350_linear_decider as L  # noqa: E402
import w33_pass11352_three_qutrit_geometry as GEO  # noqa: E402
import w33_pass11421_four_qutrits as F4  # noqa: E402
import w33_pass11457_magic_axis_theorem as TH  # noqa: E402

OUT = ROOT / "data" / "w33_pass11487_magic_axis_extension.json"


def weyl_moduli(wl, G):
    """m_G(p) = |tr(W(p)^dag G)| / D for every label p, from the monomial action of W(p) (no operator table)"""
    P = wl.labels.astype(np.int64)
    ph, k = wl.act(P[:, None, :], wl.digits[None, :, :])                   # (N, D) phases and targets
    tr = np.einsum('pj,pj->p', np.conj(ph), G[k, np.arange(wl.D)[None, :]])
    return np.abs(tr) / wl.D


def level_hyperplane(m, lab, R, tol=1e-7):
    """if supp(m) = span(R) and the level sets are the three cosets of one index-3 subspace K of it, return a row basis of
    K; else None"""
    supp = np.flatnonzero(m > tol)
    if len(supp) != 3 ** len(R):
        return None
    vals = np.round(m[supp], 6)
    levels = [supp[vals == u] for u in np.unique(vals)]
    if len(levels) != 3 or len({len(x) for x in levels}) != 1:
        return None
    base = lab[levels[0][0]]
    K = TH.span_basis(((lab[levels[0]] - base) % 3).T) if len(levels[0]) > 1 else np.zeros((0, lab.shape[1]), np.int64)
    if 3 ** len(K) != len(levels[0]) or len(K) != len(R) - 1:
        return None
    for lev in levels:                                                       # each level set is a coset of K
        d = (lab[lev] - lab[lev[0]]) % 3
        if not all(TH.in_span(x, K) for x in d):
            return None
    return K


def intersect(A, B, N2):
    """row basis of span(A) cap span(B) over F3"""
    if len(A) == 0 or len(B) == 0:
        return np.zeros((0, N2), np.int64)
    sol = L.solve_affine(np.concatenate([A, -B % 3]).T % 3, np.zeros(N2, np.int64))
    if not sol[1]:
        return np.zeros((0, N2), np.int64)
    C = np.array(sol[1], dtype=np.int64)[:, :len(A)]
    X = (C @ A) % 3
    return TH.span_basis(X.T) if X.any() else np.zeros((0, N2), np.int64)


def extension_frames(D, M, G=None):
    """(mask of frames the extension proves violating, info) or (None, reason)"""
    wl, z, Om = D.wl, D.z1, D.wl.Om
    lab = wl.labels.astype(np.int64)
    v = (z + M @ z) % 3
    if not v.any() or ((M @ M @ z - z) % 3).any():
        return None, "not F' setting"
    R, Rad = TH.radical(M, Om)
    if not TH.in_span(z, R):
        return None, "z1 not in R_M (Pass 11457 case)"
    if G is None:
        G = D.weil(M) @ D.T1
    lw = F4.LightWeyl(D.n)
    assert (lw.labels == lab).all()
    K = level_hyperplane(weyl_moduli(lw, G), lab, R)
    if K is None:
        return None, "no level hyperplane (q = 0: magnitude-blind)"
    RK = intersect(Rad, K, len(z))
    wv = (lab @ Om @ v) % 3 != 0
    orth = ((lab @ Om @ RK.T) % 3 == 0).all(axis=1) if len(RK) else np.ones(len(lab), bool)
    return wv & orth, dict(rad=len(Rad), rad_cap_K=len(RK), v_in_K=bool(TH.in_span(v, K)))


def n2_check():
    """soundness against the exhaustive n = 2 verdicts and agreement with Pass 11457's exhaustive magnitude reach"""
    import w33_pass11252_exact_reversibility as RV
    import w33_pass11330_orbit_census as O
    L.CAP = 3 ** 11
    D = L.Decider(2)
    Ms, _ = O.all_symplectic(D.wl)
    z, Om, lab = D.z1, D.wl.Om, D.wl.labels.astype(np.int64)
    rng = np.random.default_rng(0)
    st = Counter()
    for M in Ms:
        M = M % 3
        if ((M @ M @ z - z) % 3).any() or not ((z + M @ z) % 3).any() or not TH.in_span(z, TH.radical(M, Om)[0]):
            continue
        st["z1 in R_M classes"] += 1
        mask, info = extension_frames(D, M)
        if mask is None:
            st[info] += 1
            continue
        g = D.good_frames(M)
        if g is None:
            V = D.weil(M)
            g = np.array([RV.decide(D.wl.W[a] @ V @ D.T1, 2, rng)[0] is True for a in range(len(lab))])
        wv = (lab @ Om @ ((z + M @ z) % 3)) % 3 != 0
        st["hyperplane classes"] += 1
        st["sound"] += bool((~g)[mask].all())
        st["F' true on class (all omega(a,v)!=0 frames violate)"] += bool((~g)[wv].all())
        st["fully certified" if mask[wv].all() else ("none certified" if not mask.any() else "partly certified")] += 1
        st[f"dim(Rad cap K) = {info['rad_cap_K']}"] += 1
    return dict(st)


def coverage(reps, D, label):
    lab = D.wl.labels.astype(np.int64)
    st, mass = Counter(), Counter()
    for M, size in reps:
        z = D.z1
        old = TH.theorem_frames(M, z, D.wl.Om, lab)
        wv = (lab @ D.wl.Om @ ((z + M @ z) % 3)) % 3 != 0
        if old is not None and old[0][wv].all():
            key = "fully covered by Pass 11457"
        else:
            mask, info = extension_frames(D, M)
            if mask is None:
                key = info if old is None else "partial (Pass 11457), extension n/a"
            else:
                key = "fully covered by extension" if mask[wv].all() else f"extension partial (dim Rad cap K = {info['rad_cap_K']})"
        st[key] += 1
        mass[key] += size
    tot = sum(mass.values())
    out = dict(set=label, classes=dict(st), mass_fraction={k: v / tot for k, v in mass.items()})
    print(out, flush=True)
    return out


def n3_soundness(D3):
    """n = 3 orbits with z1 in R_M: soundness wherever the table decider reaches"""
    st = Counter()
    lab = D3.wl.labels.astype(np.int64)
    for M, size in TH.n3_reps():
        mask, info = extension_frames(D3, M)
        if mask is None:
            continue
        g = D3.good_frames(M)
        if g is None:
            st["undecided"] += 1
            continue
        st["checked"] += 1
        st["sound"] += bool((~g)[mask].all())
        wv = (lab @ D3.wl.Om @ ((D3.z1 + M @ D3.z1) % 3)) % 3 != 0
        st["F' true on class"] += bool((~g)[wv].all())
    return dict(st)


def run(stages=("n2", "n3", "n4")):
    res = json.load(open(OUT)) if OUT.exists() else dict(pass_id=11487)
    if "n2" in stages:
        res["n2"] = n2_check()
        print(res["n2"], flush=True)
        json.dump(res, open(OUT, "w"), indent=1, default=str)
    if "n3" in stages:
        D3 = L.Decider(3)
        res["n3_soundness"] = n3_soundness(D3)
        print(res["n3_soundness"], flush=True)
        res["n3_orbits"] = coverage(TH.n3_reps(), D3, "n=3 orbits in the two bad cells (Pass 11373)")
        json.dump(res, open(OUT, "w"), indent=1, default=str)
    if "n4" in stages:
        D4 = F4.SparseDecider(4)
        slow = []
        for M, size in TH.n4_classes():
            dims = [len(sol[1]) for sol in (D4.symplectic_solutions(M, k) for k in range(3)) if sol is not None]
            if dims and max(3 ** d for d in dims) > 3 ** 11:
                slow.append((M, size))
        res["n4_unreached"] = coverage(slow, D4, "n=4: the 317 Stab(z1) classes the decider could not reach (Pass 11433)")
        res["n4_unreached"]["count"] = len(slow)
        json.dump(res, open(OUT, "w"), indent=1, default=str)
    return res


def main():
    stages = [a for a in sys.argv[1:] if a in ("n2", "n3", "n4")] or ["n2", "n3", "n4"]
    run(stages)


if __name__ == "__main__":
    main()
