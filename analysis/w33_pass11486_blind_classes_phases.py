"""Pass 11486: the magnitude-blind (q = 0) classes of the magic-axis law, closed by phases.

Setting (Passes 11457, 11487): U = W(a) V_M T1, M^2 z1 = z1, v = z1 + M z1 != 0, z1 in R_M = Im(M - I), Rad = R_M cap R_M^perp,
and the twisted moduli |sum_j tau_j omega^(j s + q j^2)| constant (q = 0), so magnitudes certify nothing.

LEMMA 1 (pure phases).  For q = 0 the coefficients of G = V_M T1 are pure phases on R_M:
    c_G(r) = const * omega^(Q(r)) * zeta^(s(r)^3)        (zeta = e^(2 pi i/9), Q quadratic, s = l + kappa affine on R_M)
because sum_j tau_j omega^(j s) is the diagonal of T, zeta^(s^3) for s in {0, 1, -1}.  [checked per class: the exponent
function e = arg(c)/(2 pi/9) admits e = 3Q + s^3 mod 9]
LEMMA 2 (exact F3 conditions).  By Pass 11252 (c_U(Lp) = mu omega^<b,Lp> c_U(p)), U is reversible iff some anti-symplectic L
with L R_M = R_M and t = L a - a in R_M satisfies e(L r + t) - e(r) in 3 Aff + const.  Reducing mod 3 (s^3 = s) gives
    (i)  l(L r) = l(r) on R_M,
and mod 9, with c0 = l(t) and (s + c0)^3 - s^3 = 3 c0 s^2 + 3 c0^2 s + c0^3,
    (ii) B(L r, L r') = B(r, r') - 2 c0 l(r) l(r')          (B the polar form of Q's quadratic part).
[checked: the (L, t) satisfying the exact criterion reproduce the exhaustive verdicts on every n = 2 blind class]

THEOREM.  If l(v) = 0 and l = mu B(v, .) on R_M with mu != 0, then every frame a with omega(a, v) != 0 and
omega(a, Rad cap rad(B|R_M)) = 0 is VIOLATING.
PROOF.  Put r = v in (ii): l(v) = 0 gives B(Lv, L r') = B(v, r').  By (i), B(v, L r') = mu^-1 l(L r') = mu^-1 l(r') = B(v, r').
Subtracting, B(Lv - v, L r') = 0 for all r', and L is bijective on R_M, so r0 := Lv - v is in rad(B|R_M); also r0 is in Rad
(L preserves R_M and R_M^perp).  Then -omega(a, v) = omega(L a, L v) = omega(a + t, v + r0) = omega(a, v) + omega(a, r0)
= omega(a, v), so 2 omega(a, v) = 0: contradiction.  QED.
At n = 2 the hypotheses hold on all 168 blind classes and B|R_M is nondegenerate, so F' is proved on every one of them.
"""

from __future__ import annotations

import itertools
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from sympy import GF
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11350_linear_decider as L  # noqa: E402
import w33_pass11421_four_qutrits as F4  # noqa: E402
import w33_pass11457_magic_axis_theorem as TH  # noqa: E402
import w33_pass11487_magic_axis_extension as X  # noqa: E402

OUT = ROOT / "data" / "w33_pass11486_blind_classes_phases.json"
LIFT = {0: 0, 1: 1, 2: -1}


def rank3(A):
    A = np.asarray(A) % 3
    if A.size == 0:
        return 0
    return DomainMatrix([[GF(3)(int(x)) for x in row] for row in A], A.shape, GF(3)).rank()


def null3(A):
    """row basis of {x : A x = 0} over F3"""
    sol = L.solve_affine(np.asarray(A) % 3, np.zeros(len(A), dtype=np.int64))
    return np.array(sol[1], dtype=np.int64) % 3 if sol and sol[1] else np.zeros((0, np.asarray(A).shape[1]), np.int64)


class Coeffs:
    def __init__(self, n):
        self.lw = F4.LightWeyl(n)
        lab = self.lw.labels.astype(np.int64)
        self.ph, self.k = self.lw.act(lab[:, None, :], self.lw.digits[None, :, :])

    def __call__(self, G):
        return np.einsum('pj,pj->p', np.conj(self.ph), G[self.k, np.arange(self.lw.D)[None, :]]) / self.lw.D


def phase_data(D, CO, M):
    """None if not a blind class; else dict with R, Rad, coordinates, l, B, v-coordinates, decomposition flag"""
    z, Om = D.z1, D.wl.Om
    lab = CO.lw.labels.astype(np.int64)
    pow3 = CO.lw.pow3
    if ((M @ M @ z - z) % 3).any() or not ((z + M @ z) % 3).any():
        return None
    R, Rad = TH.radical(M, Om)
    if not TH.in_span(z, R):
        return None
    c = CO(D.weil(M) @ D.T1)
    if X.level_hyperplane(np.abs(c), lab, R) is not None:
        return None
    r = len(R)
    coords = np.array(list(itertools.product(range(3), repeat=r)), dtype=np.int64)
    pts = (coords @ R) % 3
    cp = c[pts @ pow3]
    if np.abs(np.abs(cp) - np.abs(cp[0])).max() > 1e-9 or np.abs(cp[0]) < 1e-9:      # Lemma 1: pure phases on R_M
        return dict(R=R, Rad=Rad, decomposed=False, reason="moduli not constant on R_M")
    e = np.angle(cp / cp[0]) * 9 / (2 * np.pi)
    if np.abs(e - np.rint(e)).max() > 1e-6:
        return dict(R=R, Rad=Rad, decomposed=False, reason="phases not ninth roots")
    ev = np.rint(e).astype(np.int64) % 9
    mons = [()] + [(i,) for i in range(r)] + [(i, j) for i in range(r) for j in range(i, r)]
    A = np.array([[np.prod([x[i] for i in m]) if m else 1 for m in mons] for x in coords]) % 3
    Aff = np.concatenate([np.ones((len(coords), 1), np.int64), coords], 1) % 3
    for e0 in range(9):
        f = (ev + e0) % 9
        s = f % 3
        sl = L.solve_affine(Aff, s)
        if sl is None:
            continue
        g = (f - np.array([LIFT[int(x)] ** 3 for x in s])) % 9
        if (g % 3).any():
            continue
        qs = L.solve_affine(A, (g // 3) % 3)
        if qs is None:
            continue
        lvec = sl[0][1:] % 3
        Bm = np.zeros((r, r), np.int64)
        for m, cf in zip(mons, qs[0]):
            if len(m) == 2:
                i, j = m
                Bm[i, j] += cf if i != j else 2 * cf
                if i != j:
                    Bm[j, i] += cf
        Bm %= 3
        v = (z + M @ z) % 3
        vc = L.solve_affine(R.T % 3, v)[0] % 3
        return dict(R=R, Rad=Rad, coords=coords, l=lvec, B=Bm, vc=vc, v=v, decomposed=True)
    return dict(R=R, Rad=Rad, decomposed=False)


def theorem_mask(D, pdat, require_lv0=True):
    """(mask of frames the theorem proves violating, hypotheses dict); require_lv0=False for q != 0, where the distinct
    moduli force c0 = 0 and (ii) holds without the l(v) term"""
    lab = D.wl.labels.astype(np.int64)
    Om = D.wl.Om
    R, Rad, lvec, Bm, vc, v = (pdat[k] for k in ("R", "Rad", "l", "B", "vc", "v"))
    Bv = (vc @ Bm) % 3
    lv = int(lvec @ vc % 3)
    mu = [m for m in (1, 2) if ((lvec - m * Bv) % 3 == 0).all()]
    hyp = dict(l_of_v_zero=lv == 0, l_prop_Bv=bool(mu), rank_B=rank3(Bm), dim_R=len(R), dim_Rad=len(Rad))
    if not ((hyp["l_of_v_zero"] or not require_lv0) and hyp["l_prop_Bv"]):
        return None, hyp
    radB = (null3(Bm) @ R) % 3 if hyp["rank_B"] < len(R) else np.zeros((0, R.shape[1]), np.int64)
    S = X.intersect(Rad, TH.span_basis(radB.T) if len(radB) else radB, R.shape[1])
    hyp["dim_Rad_cap_radB"] = len(S)
    wv = (lab @ Om @ v) % 3 != 0
    orth = ((lab @ Om @ S.T) % 3 == 0).all(axis=1) if len(S) else np.ones(len(lab), bool)
    return wv & orth, hyp


def hyperplane_phase_data(D, CO, M):
    """q != 0 classes with z1 in R_M: c_G = const omega^Q tau~_q(s) on R_M, s linear (Pass 11487's level function).
    Returns None unless it is a hyperplane class whose Pass 11487 certificate is only partial (Rad cap K != 0)."""
    z, Om = D.z1, D.wl.Om
    lab = CO.lw.labels.astype(np.int64)
    pow3 = CO.lw.pow3
    if ((M @ M @ z - z) % 3).any() or not ((z + M @ z) % 3).any():
        return None
    R, Rad = TH.radical(M, Om)
    if not TH.in_span(z, R):
        return None
    c = CO(D.weil(M) @ D.T1)
    K = X.level_hyperplane(np.abs(c), lab, R)
    if K is None or len(X.intersect(Rad, K, R.shape[1])) == 0:
        return None
    taus = TH.taus()
    w = np.exp(2j * np.pi / 3)

    def tt(q, s):
        return sum(taus[j] * w ** (j * s + q * j * j) for j in range(3))

    r = len(R)
    coords = np.array(list(itertools.product(range(3), repeat=r)), dtype=np.int64)
    cp = c[((coords @ R) % 3) @ pow3] * np.sqrt(3 ** r)
    best = None
    for q in (1, 2):
        mods = [abs(tt(q, s)) for s in range(3)]
        sv = np.array([int(np.argmin([abs(abs(x) - m) for m in mods])) for x in cp])
        err = max(abs(abs(x) - mods[s]) for x, s in zip(cp, sv))
        if best is None or err < best[0]:
            best = (err, q, sv)
    err, q, sv = best
    assert err < 1e-8
    ratio = cp / np.array([tt(q, s) for s in sv])
    eq = np.angle(ratio / ratio[0]) * 3 / (2 * np.pi)
    assert np.abs(eq - np.rint(eq)).max() < 1e-6                         # the rest is a pure omega-phase
    Qv = np.rint(eq).astype(np.int64) % 3
    mons = [()] + [(i,) for i in range(r)] + [(i, j) for i in range(r) for j in range(i, r)]
    A = np.array([[np.prod([x[i] for i in m]) if m else 1 for m in mons] for x in coords]) % 3
    Aff = np.concatenate([np.ones((len(coords), 1), np.int64), coords], 1) % 3
    qs, sl = L.solve_affine(A, Qv), L.solve_affine(Aff, sv % 3)
    if qs is None or sl is None:
        return dict(R=R, Rad=Rad, decomposed=False)
    Bm = np.zeros((r, r), np.int64)
    for m, cf in zip(mons, qs[0]):
        if len(m) == 2:
            i, j = m
            Bm[i, j] += cf if i != j else 2 * cf
            if i != j:
                Bm[j, i] += cf
    v = (z + M @ z) % 3
    return dict(R=R, Rad=Rad, l=sl[0][1:] % 3, B=Bm % 3, vc=L.solve_affine(R.T % 3, v)[0] % 3, v=v, decomposed=True, q=q)


def run_partial_n2():
    """Pass 11487's 48 partial hyperplane classes: phases close them (c0 = 0 forced by the distinct moduli)"""
    import w33_pass11252_exact_reversibility as RV
    import w33_pass11330_orbit_census as O
    L.CAP = 3 ** 11
    D = L.Decider(2)
    CO = Coeffs(2)
    lab = D.wl.labels.astype(np.int64)
    rng = np.random.default_rng(0)
    st = Counter()
    for M in np.array(O.all_symplectic(D.wl)[0]) % 3:
        pdat = hyperplane_phase_data(D, CO, M)
        if pdat is None:
            continue
        st["partial hyperplane classes"] += 1
        st["c = omega^Q tau~_q(s) with Q quadratic, s affine"] += pdat["decomposed"]
        mask, hyp = theorem_mask(D, pdat, require_lv0=False)
        st["l prop B(v,.)"] += hyp["l_prop_Bv"]
        st["B|R nondegenerate"] += hyp["rank_B"] == hyp["dim_R"]
        g = D.good_frames(M)
        if g is None:
            V = D.weil(M)
            g = np.array([RV.decide(D.wl.W[a] @ V @ D.T1, 2, rng)[0] is True for a in range(len(lab))])
        wv = (lab @ D.wl.Om @ pdat["v"]) % 3 != 0
        if mask is not None:
            st["theorem sound"] += bool((~g)[mask].all())
            st["F' fully proved"] += bool(mask[wv].all())
    return dict(st)


def exact_criterion_n2(D, CO, M, pdat):
    """Lemma 2 checked directly: frames reversible by some (L, t) satisfying the exact phase criterion"""
    lab = CO.lw.labels.astype(np.int64)
    pow3 = CO.lw.pow3
    R = pdat["R"]
    c = CO(D.weil(M) @ D.T1)
    supp = np.flatnonzero(np.abs(c) > 1e-7)
    E = -np.ones(len(lab), np.int64)
    E[supp] = np.rint(np.angle(c[supp] / c[supp[0]]) * 9 / (2 * np.pi)).astype(np.int64) % 9
    inR = np.zeros(len(lab), bool)
    inR[supp] = True
    anti = (D.J[None] @ EXACT_MS) % 3
    keep = inR[(np.einsum('aij,rj->ari', anti, R) % 3) @ pow3].all(axis=1)
    Rs = lab[supp]
    ei = E[supp]
    Aff = np.concatenate([np.ones((len(Rs), 1), np.int64), Rs], 1) % 3
    pred = np.zeros(len(lab), bool)
    fixes_v = True
    for Lm in anti[keep]:
        img = (Rs @ Lm.T) % 3
        for t in Rs:
            f = (E[((img + t) % 3) @ pow3] - ei) % 9
            if len(set((f % 3).tolist())) != 1:
                continue
            if L.solve_affine(Aff, ((f - (f[0] % 3)) % 9) // 3) is None:
                continue
            pred |= (((lab @ Lm.T) - lab - t) % 3 == 0).all(axis=1)
            fixes_v &= bool(((Lm @ pdat["v"] - pdat["v"]) % 3 == 0).all())
    return pred, fixes_v


EXACT_MS = None


def run_n2():
    global EXACT_MS
    import w33_pass11252_exact_reversibility as RV
    import w33_pass11330_orbit_census as O
    L.CAP = 3 ** 11
    D = L.Decider(2)
    CO = Coeffs(2)
    Ms = np.array(O.all_symplectic(D.wl)[0]) % 3
    EXACT_MS = Ms
    lab = D.wl.labels.astype(np.int64)
    rng = np.random.default_rng(0)
    st = Counter()
    for M in Ms:
        pdat = phase_data(D, CO, M)
        if pdat is None:
            continue
        st["blind classes"] += 1
        st["e = 3Q + s^3"] += pdat["decomposed"]
        mask, hyp = theorem_mask(D, pdat)
        for k in ("l_of_v_zero", "l_prop_Bv"):
            st[k] += hyp[k]
        st["B|R nondegenerate"] += hyp["rank_B"] == hyp["dim_R"]
        g = D.good_frames(M)
        if g is None:
            V = D.weil(M)
            g = np.array([RV.decide(D.wl.W[a] @ V @ D.T1, 2, rng)[0] is True for a in range(len(lab))])
        wv = (lab @ D.wl.Om @ pdat["v"]) % 3 != 0
        st["F' true"] += bool((~g)[wv].all())
        if mask is not None:
            st["theorem sound"] += bool((~g)[mask].all())
            st["F' fully proved"] += bool(mask[wv].all())
        pred, fixes = exact_criterion_n2(D, CO, M, pdat)
        st["Lemma 2 criterion = exhaustive verdicts"] += bool((pred == g).all())
        st["every admissible L fixes v"] += fixes
    return dict(st)


def run_n3():
    D = L.Decider(3)
    CO = Coeffs(3)
    lab = D.wl.labels.astype(np.int64)
    st = Counter()
    mass = Counter()
    for M, size in TH.n3_reps():
        pdat = phase_data(D, CO, M)
        if pdat is None:
            continue
        st["blind orbits"] += 1
        mass["blind"] += size
        st["e = 3Q + s^3"] += pdat["decomposed"]
        mask, hyp = theorem_mask(D, pdat)
        st["l(v) = 0"] += hyp["l_of_v_zero"]
        st["l prop B(v,.)"] += hyp["l_prop_Bv"]
        wv = (lab @ D.wl.Om @ pdat["v"]) % 3 != 0
        key = "hypotheses fail" if mask is None else ("F' fully proved" if mask[wv].all() else
                                                       f"partial (dim Rad cap rad B = {hyp['dim_Rad_cap_radB']})")
        st[key] += 1
        mass[key] += size
        g = D.good_frames(M)
        if g is not None and mask is not None:
            st["decided"] += 1
            st["theorem sound on decided"] += bool((~g)[mask].all())
    tot = mass.pop("blind")
    return dict(counts=dict(st), mass_fraction_of_blind={k: v / tot for k, v in mass.items()})


def transvection_case(D, M):
    """THEOREM (transvection case of Pass 11457's partial classes).  If z1 is not in R_M and R_M = span(r1) is a line (M a
    symplectic transvection x -> x + c omega(x, r1) r1), then M^2 z1 = z1 forces omega(z1, r1) = 0, so M z1 = z1 and
    v = -z1; any reversing L fixes the three cosets j z1 + R_M (distinct |tau_j|), so L z1 = z1 + kappa r1, and on each
    coset the coefficient phase is omega^(q1 lam^2) with q1 = chi(r1, r1) the Wall form, nonzero.  Matching the j*lam cross
    term of the exact criterion forces 2 q1 kappa = 0, so kappa = 0, L v = v, and every frame with omega(a, v) != 0
    violates.  Returns the hypotheses as checked (None if not in the case)."""
    z, Om = D.z1, D.wl.Om
    if ((M @ M @ z - z) % 3).any() or not ((z + M @ z) % 3).any():
        return None
    R, Rad = TH.radical(M, Om)
    if len(R) != 1 or TH.in_span(z, R):
        return None
    r1 = R[0]
    N2 = len(z)
    sol = L.solve_affine((M - np.eye(N2, dtype=np.int64)) % 3, r1)       # (M - I) x = r1
    x = sol[0] % 3
    wall = int(x @ Om @ r1 % 3)                                          # chi(r1, r1) = omega(x, r1)
    return dict(Mz1_eq_z1=bool(((M @ z - z) % 3 == 0).all()), omega_z1_r1=int(z @ Om @ r1 % 3), wall_q1=wall,
                holds=bool(((M @ z - z) % 3 == 0).all() and wall != 0))


def run_transvection_n2():
    import w33_pass11252_exact_reversibility as RV
    import w33_pass11330_orbit_census as O
    L.CAP = 3 ** 11
    D = L.Decider(2)
    lab = D.wl.labels.astype(np.int64)
    rng = np.random.default_rng(0)
    st = Counter()
    for M in np.array(O.all_symplectic(D.wl)[0]) % 3:
        z = D.z1
        if ((M @ M @ z - z) % 3).any() or not ((z + M @ z) % 3).any():
            continue
        v = (z + M @ z) % 3
        wv = (lab @ D.wl.Om @ v) % 3 != 0
        old = TH.theorem_frames(M, z, D.wl.Om, lab)
        if old is None or old[0][wv].all():
            continue
        st["Pass 11457 partial classes"] += 1
        h = transvection_case(D, M)
        if h is None:
            st["not a transvection"] += 1
            continue
        st["transvection, Mz1 = z1, Wall q1 != 0"] += h["holds"]
        g = D.good_frames(M)
        if g is None:
            V = D.weil(M)
            g = np.array([RV.decide(D.wl.W[a] @ V @ D.T1, 2, rng)[0] is True for a in range(len(lab))])
        st["F' true (theorem sound)"] += bool((~g)[wv].all())
    return dict(st)


def classify_all(D, CO, M):
    """(key, mask or None): which theorem proves F' for the class, trying 11457, 11487, then the phase theorems"""
    lab = D.wl.labels.astype(np.int64)
    z = D.z1
    v = (z + M @ z) % 3
    wv = (lab @ D.wl.Om @ v) % 3 != 0
    old = TH.theorem_frames(M, z, D.wl.Om, lab)
    if old is not None and old[0][wv].all():
        return "fully proved by Pass 11457", old[0], wv
    ext, _ = X.extension_frames(D, M)
    if ext is not None and ext[wv].all():
        return "fully proved by Pass 11487", ext, wv
    pdat = phase_data(D, CO, M)
    req = True
    if pdat is None:
        pdat = hyperplane_phase_data(D, CO, M)
        req = False
    tv = transvection_case(D, M) if pdat is None else None
    if tv is not None and tv["holds"]:
        return "transvection: fully proved by Pass 11486", wv, wv
    if pdat is None or not pdat.get("decomposed"):
        return ("outside the phase theorems (z1 not in R_M, Rad != 0)" if pdat is None else "no decomposition"), None, wv
    mask, hyp = theorem_mask(D, pdat, require_lv0=req)
    kind = "blind (q = 0)" if req else "partial hyperplane (q != 0)"
    if mask is None:
        return f"{kind}: hypotheses fail", None, wv
    return (f"{kind}: fully proved by Pass 11486" if mask[wv].all() else f"{kind}: partly proved by Pass 11486"), mask, wv


def run_n3_all():
    """every F'-setting orbit at n = 3 (Pass 11373's orbits in the two bad cells): which theorem covers it, by orbit mass,
    with soundness against the decider wherever it reaches"""
    D = L.Decider(3)
    CO = Coeffs(3)
    st, mass = Counter(), Counter()
    sound = Counter()
    for M, size in TH.n3_reps():
        z = D.z1
        if ((M @ M @ z - z) % 3).any() or not ((z + M @ z) % 3).any():
            continue
        key, mask, wv = classify_all(D, CO, M)
        st[key] += 1
        mass[key] += size
        if mask is not None:
            g = D.good_frames(M)
            if g is not None:
                sound["decided"] += 1
                sound["sound"] += bool((~g)[mask].all())
    tot = sum(mass.values())
    proved = sum(v for k, v in mass.items() if "fully proved" in k) / tot
    return dict(classes=dict(st), mass_fraction={k: v / tot for k, v in mass.items()}, soundness=dict(sound),
                fully_proved_mass=proved)


def run_n4():
    """the Stab(z1) classes the n = 4 decider could not reach (Pass 11433): how many does the phase theorem close?
    (no exhaustive verdicts exist here; the certificate is the theorem itself, hypotheses checked per class)"""
    D = F4.SparseDecider(4)
    CO = Coeffs(4)
    lab = D.wl.labels.astype(np.int64)
    slow = []
    for M, size in TH.n4_classes():
        dims = [len(sol[1]) for sol in (D.symplectic_solutions(M, k) for k in range(3)) if sol is not None]
        if dims and max(3 ** d for d in dims) > 3 ** 11:
            slow.append((M, size))
    st, mass = Counter(), Counter()
    for M, size in slow:
        z = D.z1
        v = (z + M @ z) % 3
        wv = (lab @ D.wl.Om @ v) % 3 != 0
        old = TH.theorem_frames(M, z, D.wl.Om, lab)
        if old is not None and old[0][wv].all():
            key = "fully proved by Pass 11457"
        else:
            ext, _ = X.extension_frames(D, M)
            if ext is not None and ext[wv].all():
                key = "fully proved by Pass 11487"
            else:
                pdat = phase_data(D, CO, M)
                req = True
                if pdat is None:
                    pdat = hyperplane_phase_data(D, CO, M)
                    req = False
                tv = transvection_case(D, M) if pdat is None else None
                if tv is not None and tv["holds"]:
                    key = "transvection: fully proved by Pass 11486"
                elif pdat is None or not pdat.get("decomposed"):
                    key = "outside the phase theorems (z1 not in R_M, Rad != 0)" if pdat is None else "no decomposition"
                else:
                    mask, hyp = theorem_mask(D, pdat, require_lv0=req)
                    kind = "blind (q = 0)" if req else "partial hyperplane (q != 0)"
                    if mask is None:
                        key = f"{kind}: hypotheses fail"
                    elif mask[wv].all():
                        key = f"{kind}: fully proved by Pass 11486"
                    else:
                        key = f"{kind}: partly proved by Pass 11486"
        st[key] += 1
        mass[key] += size
    tot = sum(mass.values())
    return dict(count=len(slow), classes=dict(st), mass_fraction={k: v / tot for k, v in mass.items()})


def main():
    res = json.load(open(OUT)) if OUT.exists() else dict(pass_id=11486)
    stages = [a for a in sys.argv[1:] if a in ("n2", "n3", "partial", "n4", "transvection", "n3all")] or [
        "n2", "n3", "partial", "transvection", "n4", "n3all"]
    if "n3all" in stages:
        res["n3_all_theorems"] = run_n3_all()
        print(res["n3_all_theorems"], flush=True)
        json.dump(res, open(OUT, "w"), indent=1, default=str)
    if "transvection" in stages:
        res["n2_transvection"] = run_transvection_n2()
        print(res["n2_transvection"], flush=True)
        json.dump(res, open(OUT, "w"), indent=1, default=str)
    if "partial" in stages:
        res["n2_partial_hyperplane"] = run_partial_n2()
        print(res["n2_partial_hyperplane"], flush=True)
        json.dump(res, open(OUT, "w"), indent=1, default=str)
    if "n2" in stages:
        res["n2"] = run_n2()
        print(res["n2"], flush=True)
        json.dump(res, open(OUT, "w"), indent=1, default=str)
    if "n3" in stages:
        res["n3"] = run_n3()
        print(res["n3"], flush=True)
        json.dump(res, open(OUT, "w"), indent=1, default=str)
    if "n4" in stages:
        res["n4_unreached"] = run_n4()
        print(res["n4_unreached"], flush=True)
        json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
