"""Pass 11457: a proof of the magic-axis law F' whenever Im(M - I) is nondegenerate -- by Weyl-coefficient magnitudes.

Setting: U = W(a) V_M T1 (n qutrits), z1 the magic axis, v = z1 + M z1.  F' (Pass 11420): if M^2 z1 = z1 then every frame
with omega(a, v) != 0 violates.

LEMMA 1 (support; standard, checked).  The Weyl coefficients c_{V_M}(p) = tr(W(p)^dag V_M)/D vanish off R_M = Im(M - I)
and have constant modulus |R_M|^(-1/2) on it.
LEMMA 2 (magnitudes are a reversal invariant).  If C U C^dag ~ U^T for a Clifford C = W(b) V_S, then
|c_U(S^-1 p)| = |c_U(-J p)| for all p, i.e. m_U = |c_U| is invariant under the ANTI-symplectic linear map L = -J S.
LEMMA 3 (the three levels).  T1 = sum_j tau_j W(j z1) with tau_j = (1 + zeta^(1-3j) + zeta^(8-6j))/3, whose moduli
0.844, 0.449, 0.293 are distinct.  Hence c_U(p) = phase * sum_j tau_j c_{V_M}(p - a - j z1) and, if z1 is not in R_M, m_U
takes the three distinct nonzero values |tau_j| |R_M|^(-1/2) exactly on the three cosets R_M + j z1 + a.

THEOREM.  Let M^2 z1 = z1 (so v = z1 + M z1 is fixed by M, v in ker(M - I) = R_M^perp, and v in 2 z1 + R_M), v != 0,
and let Rad = R_M cap R_M^perp.  Then every frame a with omega(a, v) != 0 AND omega(a, Rad) = 0 is VIOLATING.
In particular, if R_M is nondegenerate (Rad = 0), F' holds for M -- every frame with omega(a, v) != 0 violates.
PROOF.  If z1 were in R_M then v in R_M cap R_M^perp = Rad; Rad = 0 is the main case, and in general the argument
below needs only z1 not in R_M, which holds when v not in Rad (checked per class).  Suppose U reversible; take L from
Lemma 2.  L permutes the level sets of m_U and they carry distinct values, so L fixes each coset R_M + j z1 + a.  Hence
L(R_M) = R_M, La - a in R_M, Lz1 - z1 in R_M.  L anti-symplectic and L(R_M) = R_M give L(R_M^perp) = R_M^perp, so
Lv in R_M^perp; and Lv in 2 Lz1 + R_M = v + R_M.  So r := Lv - v lies in Rad.  Then
    -omega(a, v) = omega(La, Lv) = omega(a + r2, v + r) = omega(a, v) + omega(a, r)        (r2 in R_M, v and r in R_M^perp)
and omega(a, r) = 0 because r in Rad and a is orthogonal to Rad.  So 2 omega(a, v) = 0, i.e. omega(a, v) = 0 mod 3: a
contradiction.  QED.

Computed here: Lemma 1 on random classes (n = 2, 3); soundness of the theorem against the exhaustive n = 2 verdicts;
coverage of F' by the theorem at n = 2 (all classes), n = 3 (all orbits of Pass 11373 in the two cells) and n = 4 (the
549 conjugacy classes of Stab(z1), Pass 11433), in particular on the classes the decider could not reach.
"""

from __future__ import annotations

import ast
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11350_linear_decider as L  # noqa: E402
import w33_pass11352_three_qutrit_geometry as GEO  # noqa: E402

OUT = ROOT / "data" / "w33_pass11457_magic_axis_theorem.json"


def taus():
    z = np.exp(2j * np.pi / 9)
    return [complex((1 + z ** (1 - 3 * j) + z ** (8 - 6 * j)) / 3) for j in range(3)]


def span_basis(B):
    """row basis of the column space of B over F3 (as rows)"""
    R, piv = L.rref3(B.T.copy() % 3)
    return R % 3


def perp_basis(rows, Om, N2):
    """basis of {x : omega(r, x) = 0 for all rows r}"""
    if len(rows) == 0:
        return np.eye(N2, dtype=np.int64)
    A = (rows @ Om) % 3
    sol = L.solve_affine(A, np.zeros(len(rows), dtype=np.int64))
    return np.array(sol[1], dtype=np.int64) if sol[1] else np.zeros((0, N2), dtype=np.int64)


def in_span(x, rows):
    if len(rows) == 0:
        return not (x % 3).any()
    return L.solve_affine(rows.T % 3, x % 3) is not None


def radical(M, Om):
    N2 = M.shape[0]
    R = span_basis((M - np.eye(N2, dtype=np.int64)) % 3)            # R_M
    P = perp_basis(R, Om, N2)                                         # R_M^perp
    # Rad = R cap R^perp: vectors of R orthogonal to R
    if len(R) == 0:
        return R, R
    G = (R @ Om @ R.T) % 3                                            # Gram of R
    sol = L.solve_affine(G.T, np.zeros(len(R), dtype=np.int64))
    coeffs = np.array(sol[1], dtype=np.int64) if sol[1] else np.zeros((0, len(R)), dtype=np.int64)
    Rad = (coeffs @ R) % 3 if len(coeffs) else np.zeros((0, N2), dtype=np.int64)
    return R, Rad


def theorem_frames(M, z1, Om, labels):
    """boolean mask of frames the theorem proves violating (or None if its hypotheses fail)"""
    v = (z1 + M @ z1) % 3
    if not v.any():
        return None
    R, Rad = radical(M, Om)
    if in_span(z1, R):                                                # level sets would merge
        return None
    wv = (labels @ Om @ v) % 3 != 0
    if len(Rad):
        orth = ((labels @ Om @ Rad.T) % 3 == 0).all(axis=1)
    else:
        orth = np.ones(len(labels), dtype=bool)
    return wv & orth, int(len(Rad))


def support_check(n, trials, seed=0):
    D = L.Decider(n)
    rng = np.random.default_rng(seed)
    gens = GEO.gen_mats(n)
    W, lab = D.wl.W, D.wl.labels.astype(np.int64)
    ok = 0
    for _ in range(trials):
        M = GEO.random_symplectic(rng, gens)
        c = np.abs(np.einsum('pji,ji->p', W.conj(), D.weil(M))) / D.wl.D
        supp = set(np.flatnonzero(c > 1e-9).tolist())
        B = (M - np.eye(2 * n, dtype=np.int64)) % 3
        img = {int(((B @ x) % 3) @ D.wl.pow3) for x in lab}
        ok += supp == img and np.ptp(c[list(supp)]) < 1e-9
    return dict(n=n, classes=trials, support_is_Im_M_minus_I_uniform=int(ok))


def n2_soundness_and_coverage():
    import w33_pass11330_orbit_census as O
    import w33_pass11252_exact_reversibility as R
    L.CAP = 3 ** 11
    D = L.Decider(2)
    Ms, _ = O.all_symplectic(D.wl)
    z, Om, lab = D.z1, D.wl.Om, D.wl.labels.astype(np.int64)
    rng = np.random.default_rng(0)
    st = Counter()
    for M in Ms:
        if not ((M @ M @ z - z) % 3 == 0).all() or not ((z + M @ z) % 3).any():
            continue
        tf = theorem_frames(M, z, Om, lab)
        g = D.good_frames(M)
        if g is None:
            V = D.weil(M)
            g = np.array([R.decide(D.wl.W[a] @ V @ D.T1, 2, rng)[0] is True for a in range(len(lab))])
        if tf is None:
            st["hypotheses fail"] += 1
            continue
        mask, rad = tf
        st["classes"] += 1
        st["sound (all theorem frames violating)"] += bool((~g)[mask].all())
        wv = (lab @ Om @ ((z + M @ z) % 3)) % 3 != 0
        st["F' fully covered (Rad orthogonal to every omega(a,v)!=0 frame)"] += bool(mask[wv].all())
        st[f"Rad dim {rad}"] += 1
    return dict(st)


def coverage_from_reps(reps, n):
    D = L.Decider(n)
    z, Om, lab = D.z1, D.wl.Om, D.wl.labels.astype(np.int64)
    st = Counter()
    mass = Counter()
    for M, size in reps:
        tf = theorem_frames(M, z, Om, lab)
        key = "hypotheses fail" if tf is None else (
            "fully covered" if tf[0][((lab @ Om @ ((z + M @ z) % 3)) % 3 != 0)].all() else f"partial (Rad dim {tf[1]})")
        st[key] += 1
        mass[key] += size
    tot = sum(mass.values())
    return dict(classes=dict(st), mass_fraction={k: v / tot for k, v in mass.items()})


def n3_reps():
    import w33_pass11373_three_qutrit_exact_fraction as X
    D = L.Decider(3)
    out = []
    for M, size in X.read_orbits(3):
        if GEO.cell(M, D.z1, D.wl.Om) in ("Mz1 = z1", "same line, M^2 z1 = z1"):
            out.append((M, size))
    return out


def n4_classes():
    lines = (ROOT / "data" / "w33_pass11433_stab_z1_classes_n4.txt").read_text().splitlines()
    out = []
    for line in lines:
        m, size = line.rsplit(";", 1)
        out.append((np.array(ast.literal_eval(m), dtype=np.int64).T % 3, int(size)))
    return out


def run():
    res = dict(pass_id=11457, tau_moduli=[abs(t) for t in taus()])
    res["support_lemma"] = [support_check(2, 400), support_check(3, 60)]
    print(res["support_lemma"], flush=True)
    res["n2"] = n2_soundness_and_coverage()
    print(res["n2"], flush=True)
    res["n3_orbits"] = coverage_from_reps(n3_reps(), 3)
    print(res["n3_orbits"], flush=True)
    res["n4_stab_z1_classes"] = coverage_from_reps(n4_classes(), 4)
    print(res["n4_stab_z1_classes"], flush=True)
    # the classes the decider could not reach in Pass 11433 (same cap test as the decider)
    import w33_pass11421_four_qutrits as F4
    D4 = F4.SparseDecider(4)
    slow = []
    for M, size in n4_classes():
        dims = [len(sol[1]) for sol in (D4.symplectic_solutions(M, k) for k in range(3)) if sol is not None]
        if dims and max(3 ** d for d in dims) > 3 ** 11:
            slow.append((M, size))
    res["n4_unreached_classes"] = dict(count=len(slow), coverage=coverage_from_reps(slow, 4))
    print(res["n4_unreached_classes"], flush=True)
    return res


def magnitude_obstruction_n2(sample=60, seed=11457):
    """the reach of ANY magnitude-based proof at n = 2.  For U = W(a) G, m_U(p) = m_G(p - a); an anti-symplectic L
    preserves m_U iff the AFFINE map q -> L q + (L a - a) preserves m_G.  So per class: find all affine automorphisms
    (L, t) of m_G with L anti-symplectic (L must preserve R_M, the direction space of the support), then a frame a is
    'magnitude-reversible' iff (L - I) a = t for one of them; otherwise magnitudes alone certify the violation.
    Exhaustive over all 51840 anti-symplectic maps, on a random sample of classes from the two bad cells."""
    import w33_pass11330_orbit_census as O
    D = L.Decider(2)
    Ms, _ = O.all_symplectic(D.wl)
    Ms = np.array(Ms) % 3
    z, Om, W = D.z1, D.wl.Om, D.wl.W
    lab = D.wl.labels.astype(np.int64)
    pow3 = D.wl.pow3
    anti = (D.J[None] @ Ms) % 3
    img = ((np.einsum('aij,pj->api', anti, lab) % 3) @ pow3)          # (51840, 81)
    add = (((lab[:, None, :] + lab[None, :, :]) % 3) @ pow3)            # (81, 81): index of x + t
    bad_cells = [i for i, M in enumerate(Ms) if GEO.cell(M, z, Om) in ("Mz1 = z1", "same line, M^2 z1 = z1")]
    rng = np.random.default_rng(seed)
    st = Counter()
    for i in rng.choice(bad_cells, size=min(sample, len(bad_cells)), replace=False):
        M = Ms[i]
        G = D.weil(M) @ D.T1
        mG = np.round(np.abs(np.einsum('pji,ji->p', W.conj(), G)) / D.wl.D, 9)
        R, Rad = radical(M, Om)
        autos = []
        for c0 in range(0, len(anti), 4096):                           # vectorised over maps and translations
            X = img[c0:c0 + 4096]                                        # (B, 81) index of L q
            vals = mG[add[X]]                                            # (B, 81 q, 81 t): m_G(L q + t)
            okt = (vals == mG[None, :, None]).all(axis=1)               # (B, 81 t)
            for bi, ti in zip(*np.nonzero(okt)):
                autos.append((anti[c0 + bi], lab[ti]))
        v = (z + M @ z) % 3
        frames = np.flatnonzero((lab @ Om @ v) % 3 != 0)
        certified = 0
        for a in frames:
            av = lab[a]
            if not any((((Lm - np.eye(4, dtype=np.int64)) @ av - tt) % 3 == 0).all() for Lm, tt in autos):
                certified += 1
        kind = "Rad = 0" if len(Rad) == 0 else ("z1 in R_M" if in_span(z, R) else "Rad != 0, z1 not in R_M")
        st[(kind, "all F' frames certified" if certified == len(frames) else
            ("none certified" if certified == 0 else "some certified"))] += 1
    return {f"{k[0]} | {k[1]}": v for k, v in sorted(st.items())}


def twisted_moduli():
    """extension route: |sum_j tau_j omega^(j s + q j^2)| for s = 0,1,2 and q = 0,1,2 (uniform for q = 0, three distinct
    values otherwise)"""
    w = np.exp(2j * np.pi / 3)
    t = taus()
    return {str(q): [abs(sum(t[j] * w ** (j * s + q * j * j) for j in range(3))) for s in range(3)] for q in range(3)}


def main():
    if "--magnitude" in sys.argv:
        res = json.load(open(OUT))
        res["magnitude_obstruction_n2"] = magnitude_obstruction_n2()
        json.dump(res, open(OUT, "w"), indent=1, default=str)
        print(json.dumps(res["magnitude_obstruction_n2"], indent=1))
        return
    if "--moduli" in sys.argv:
        res = json.load(open(OUT))
        res["extension_route_twisted_moduli"] = twisted_moduli()
        json.dump(res, open(OUT, "w"), indent=1, default=str)
        print(res["extension_route_twisted_moduli"])
        return
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
