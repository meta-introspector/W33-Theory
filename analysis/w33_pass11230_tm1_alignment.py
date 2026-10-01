#!/usr/bin/env python3
"""Pass 11230: can a flavon potential pick TM1, and what TM1 then predicts for delta and theta23.

Pass 11225 found that the only symmetry column of Aut(W33) / Sp(4,3) within ~1 sigma of the data is TM1,
(2/3, 1/6, 1/6), from the S4 quotient of the line stabiliser 3^3:S4.  TM1 needs a residual Z3 on the charged leptons
and a specific Z2 on the neutrinos.  Residual symmetries are left over by flavon vacua.  This pass asks, by explicit
computation, whether the most general S4-invariant potential puts flavons there.

Part A (sum rules).  For each surviving first column (TM1, the qutrit-Clifford column (2/3)sin^2(k pi/9), the golden
column) the unitarity relation |U_mu1|^2 = target fixes cos(delta) as a function of (theta13, theta23).  Over the
NuFIT 6.0 3-sigma box we tabulate cos(delta)(sin^2 theta23).  TM1: cos(delta) = 0 exactly at maximal theta23, and the
sign of cos(delta) is the octant of theta23.  The Clifford column only fits a narrow theta23 window with delta ~ pi.

Part B (vacuum alignment).  In the 3-dim basis where the Z3 is the cyclic permutation, the TM1 involution is
g = diag(1,-1,-1) P_(23) (computed below: its isolated eigenvector (0,1,1) has Fourier weights (2/3, 1/6, 1/6)).  A
flavon fixed by exactly <g> is the (0,1,-1) direction of the triplet 3 (or (0,1,1) of 3').  Random S4-invariant
potentials (the Reynolds average of random polynomials, so every invariant is generic) are minimised globally and the
residual symmetries of the vacuum are classified:
  * one triplet, renormalisable (degree <= 4): the global minimum is never the (0,1,-1)-type direction;
  * one triplet with degree-6 terms: it sometimes is;
  * two triplets (phi_e, phi_nu) with every allowed cross-coupling, with or without a shaping Z2 x Z2 that keeps the two
    flavons apart: how often does the vacuum pair give TM1 (G_e contains a nondegenerate Z3, G_nu contains the TM1
    involution relative to it)?
Outcome (certificate): TM1 is not a renormalisable vacuum.  It needs non-renormalisable (sextic) terms or a different
mechanism (F-term / driving-field alignment, as in the CSD / littlest-seesaw literature).  The geometry hosts TM1 but
does not, by an S4 potential alone, select it at renormalisable level.
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "w33_pass11230_tm1_alignment.json"

# NuFIT 6.0 normal-ordering 3-sigma box, as quoted in Pass 11225
BOX = dict(s12=(0.275, 0.345), s13=(0.02030, 0.02388), s23=(0.430, 0.596))
W = np.exp(2j * np.pi / 3)
FOURIER = np.array([[W ** (j * k) for k in range(3)] for j in range(3)]) / np.sqrt(3)   # columns: Z3 eigenvectors


# ---------------------------------------------------------------- Part A: sum rules
def cos_delta(col, s13, s23):
    """first PMNS column fixed to col = (|U_e1|^2, |U_mu1|^2, |U_tau1|^2); returns (sin^2 theta12, cos delta)"""
    e, mu, _ = col
    c13 = 1 - s13
    s12 = 1 - e / c13
    if not 0 <= s12 <= 1:
        return None
    c12, c23 = 1 - s12, 1 - s23
    # |U_mu1|^2 = s12 c23 + c12 s23 s13 + 2 sqrt(s12 c12 s23 c23 s13) cos delta   (squares of sines/cosines)
    cd = (mu - s12 * c23 - c12 * s23 * s13) / (2 * np.sqrt(s12 * c12 * s23 * c23 * s13))
    return s12, cd


def sum_rule_table(col, label, n=34):
    rows = []
    s13 = 0.02195                                        # NuFIT 6.0 central value region; the s13 spread is tabulated too
    for s23 in np.linspace(*BOX["s23"], n):
        vals = [cos_delta(col, x, s23) for x in (BOX["s13"][0], s13, BOX["s13"][1])]
        cds = [v[1] for v in vals]
        ok = [abs(c) <= 1 for c in cds]
        rows.append(dict(s23=round(float(s23), 4), s12=round(vals[1][0], 4), cos_delta=[round(float(c), 4) for c in cds],
                         allowed=bool(any(ok))))
    allowed = [r["s23"] for r in rows if r["allowed"]]
    return dict(column=label, target=[round(float(c), 6) for c in col], rows=rows,
                s23_allowed=[min(allowed), max(allowed)] if allowed else None)


def tm1_exact_check():
    """TM1: cos delta = 0 iff theta23 maximal (for every theta13) -- the TM1 'atmospheric sum rule'"""
    out = []
    for s13 in (0.0203, 0.022, 0.0239):
        _, cd = cos_delta((2 / 3, 1 / 6, 1 / 6), s13, 0.5)
        out.append(abs(cd) < 1e-12)
    return all(out)


# ---------------------------------------------------------------- the S4 triplet
def s4_triplet(prime=False):
    """3: D P with D = diag(+-1), det D = 1, P any permutation matrix.  3' = sign(P) * (D P)."""
    out = []
    for perm in itertools.permutations(range(3)):
        P = np.zeros((3, 3))
        for i, j in enumerate(perm):
            P[j, i] = 1
        sgn = round(np.linalg.det(P))
        for d in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)):
            M = np.diag(d) @ P
            out.append(sgn * M if prime else M)
    return np.array(out)


G3, G3P = s4_triplet(False), s4_triplet(True)


def tm1_involution():
    """the involutions of S4 whose isolated eigenvector has Fourier weights (2/3,1/6,1/6) relative to the cyclic Z3"""
    C = np.roll(np.eye(3), 1, axis=0)
    U = FOURIER
    D = U.conj().T @ C @ U
    assert np.allclose(D, np.diag(np.diag(D)))               # the Fourier basis diagonalises the cyclic Z3
    found = []
    for g in G3:
        if np.allclose(g @ g, np.eye(3)) and not np.allclose(g, np.eye(3)):
            ev, V = np.linalg.eigh(g)
            iso = V[:, 0] if np.sum(ev < 0) == 1 else V[:, 2]
            weights = np.sort(np.abs(U.conj().T @ iso) ** 2)[::-1]
            found.append((g, weights))
    return found


# ---------------------------------------------------------------- random invariant potentials
def reynolds(T, R):
    """average a symmetric tensor over the group (rows of R act on every index)"""
    out = np.zeros_like(T)
    k = T.ndim
    for g in R:
        X = T
        for i in range(k):
            X = np.moveaxis(np.tensordot(g.T, X, axes=([1], [i])), 0, i)   # T_{..i..} g_{i j}
        out += X
    return out / len(R)


def sym_random(rng, dim, k):
    return _sym_big(rng.normal(size=(dim,) * k))


def _sym_big(T):
    k = T.ndim
    perms = list(itertools.permutations(range(k)))
    return sum(np.transpose(T, p) for p in perms) / len(perms)


def contract(T, x):
    for _ in range(T.ndim):
        T = T @ x
    return T


class Potential:
    def __init__(self, tensors):
        self.T = tensors           # dict degree -> symmetric invariant tensor

    def __call__(self, x):
        return sum(contract(T, x) for T in self.T.values())

    def grad(self, x):
        g = 0
        for k, T in self.T.items():
            X = T
            for _ in range(k - 1):
                X = X @ x
            g = g + k * X
        return g


def group_block(Ra, Rb):
    return np.array([np.block([[a, np.zeros((3, 3))], [np.zeros((3, 3)), b]]) for a, b in zip(Ra, Rb)])


def random_potential(rng, R, degrees, shaping=None):
    """generic invariant: Reynolds average of random tensors; negative mass terms; stabilised by +|x|^top"""
    dim = R.shape[1]
    group = R
    if shaping is not None:          # extra Z2 x Z2 (phi -> -phi, chi -> -chi): average over it as well
        group = np.array([s @ g for s in shaping for g in R])
    tens = {}
    for k in degrees:
        T = reynolds(sym_random(rng, dim, k), group)
        if k == 2:
            T = T - 2.0 * np.eye(dim) * (1 + rng.random())      # tachyonic mass so the symmetry breaks
        tens[k] = T
    top = max(degrees)
    # stabiliser lambda |x|^top with lambda large enough that the potential is bounded
    I = np.eye(dim)
    stab = I
    for _ in range(top // 2 - 1):
        stab = np.multiply.outer(stab, I)
    stab = _sym_big(stab)
    n = rng.normal(size=(4000, dim))
    n /= np.linalg.norm(n, axis=1, keepdims=True)
    T = tens[top]
    vals = np.array([contract(T, v) for v in n[:800]])
    tens[top] = T + (1.5 * np.abs(vals).max() + 0.5) * stab
    return Potential(tens)


def global_min(V, dim, rng, starts=24):
    best = None
    for _ in range(starts):
        x0 = rng.normal(size=dim)
        r = minimize(V, x0, jac=V.grad, method="BFGS", options=dict(gtol=1e-10, maxiter=2000))
        if best is None or r.fun < best.fun - 1e-9:
            best = r
    return best.x


def stabiliser(R, v, tol=1e-5):
    n = max(np.linalg.norm(v), 1e-12)
    return [i for i, g in enumerate(R) if np.linalg.norm(g @ v - v) < tol * max(n, 1)]


def orbit_type(v, tol=1e-4):
    """the S4 stratum of a triplet vev: (1,0,0), (1,1,1), (1,1,0)-type or generic"""
    a = np.sort(np.abs(v))[::-1] / max(np.linalg.norm(v), 1e-12)
    if a[0] < tol:
        return "zero"
    if a[1] < tol:
        return "(1,0,0)"
    if abs(a[0] - a[2]) < tol:
        return "(1,1,1)"
    if a[2] < tol and abs(a[0] - a[1]) < tol:
        return "(1,1,0)"
    return "other"


def mixing_outcome(Re, Rn, ve, vn, tol=1e-3):
    """residual symmetries of the pair of vevs -> TM1 / other column / full pattern / no prediction.

    G_e = stab(phi_e) must contain an order-3 element (nondegenerate: it fixes the charged-lepton basis); G_nu =
    stab(phi_nu) must be a Z2 (one predicted column) or a Klein group (a full pattern).  A vanishing vev leaves all of
    S4 and predicts nothing."""
    if np.linalg.norm(ve) < tol or np.linalg.norm(vn) < tol:
        return "a flavon has zero vev"
    Ge = [Re[i] for i in stabiliser(Re, ve)]
    Gn = [Rn[i] for i in stabiliser(Rn, vn)]
    z3 = [g for g in Ge if np.allclose(np.linalg.matrix_power(g, 3), np.eye(3)) and not np.allclose(g, np.eye(3))]
    if not z3:
        return "G_e has no Z3"
    invol = [g for g in Gn if np.allclose(g @ g, np.eye(3)) and not np.allclose(g, np.eye(3))
             and not np.allclose(g, -np.eye(3))]
    if len(Gn) > 4 or any(len(set(np.round(np.linalg.eigvals(g), 6))) == 3 for g in Gn):
        return "G_nu too large (degenerate neutrinos)"
    if not invol:
        return "G_nu trivial (no prediction)"
    _, U = np.linalg.eig(z3[0])
    U = np.linalg.qr(U)[0]
    cols = []
    for g in invol:
        ev, V = np.linalg.eigh(g)
        iso = V[:, 0] if np.sum(ev < 0) == 1 else V[:, 2]
        cols.append(tuple(np.round(np.sort(np.abs(U.conj().T @ iso) ** 2)[::-1], 4)))
    if len(invol) >= 2:
        return "full pattern (Klein G_nu): " + ("TBM, theta13 = 0" if (0.6667, 0.1667, 0.1667) in cols else
                                                  "columns " + str(sorted(set(cols))))
    c = cols[0]
    name = {(0.6667, 0.1667, 0.1667): "TM1", (0.3333, 0.3333, 0.3333): "TM2", (0.5, 0.5, 0.0): "theta13 = 0 column"}
    return name.get(c, "column " + str(c))


def single_scan(prime, degrees, samples, rng):
    R = G3P if prime else G3
    counts = {}
    for _ in range(samples):
        V = random_potential(rng, R, degrees)
        x = global_min(V, 3, rng)
        t = f"{orbit_type(x)} |stab|={len(stabiliser(R, x))}"
        counts[t] = counts.get(t, 0) + 1
    return dict(rep="3'" if prime else "3", degrees=list(degrees), samples=samples, global_minimum_types=counts)


def pair_scan(prime_e, prime_nu, degrees, samples, rng, shaping):
    Re = G3P if prime_e else G3
    Rn = G3P if prime_nu else G3
    R = group_block(Re, Rn)
    sh = None
    if shaping:
        sh = [np.diag([a] * 3 + [b] * 3) for a in (1, -1) for b in (1, -1)]
    counts, types = {}, {}
    zero = 0
    while sum(counts.values()) < samples:
        V = random_potential(rng, R, degrees, shaping=sh)
        x = global_min(V, 6, rng)
        ve, vn = x[:3], x[3:]
        if min(np.linalg.norm(ve), np.linalg.norm(vn)) < 1e-3:      # a vanishing vev breaks nothing: resample
            zero += 1
            continue
        out = mixing_outcome(Re, Rn, ve, vn)
        counts[out] = counts.get(out, 0) + 1
        key = f"|stab e|={len(stabiliser(Re, ve))} x |stab nu|={len(stabiliser(Rn, vn))}"
        types[key] = types.get(key, 0) + 1
    return dict(phi_e="3'" if prime_e else "3", phi_nu="3'" if prime_nu else "3", degrees=list(degrees), shaping_Z2xZ2=shaping,
                samples=samples, discarded_zero_vev=zero, outcomes=counts, vev_types=types)


def run(samples=120, seed=11230):
    rng = np.random.default_rng(seed)
    res = dict(pass_id=11230)
    cols = {"TM1 (line stabiliser S4)": (2 / 3, 1 / 6, 1 / 6),
            "Clifford (2/3)sin^2(k pi/9), mu=k1": (0.646564, 0.077985, 0.275451),
            "Clifford, mu=k2": (0.646564, 0.275451, 0.077985),
            "golden A5, mu=0.0955": (0.654508, 0.095492, 0.25),
            "golden A5, mu=0.25": (0.654508, 0.25, 0.095492)}
    res["sum_rules"] = [sum_rule_table(c, k) for k, c in cols.items()]
    res["tm1_cos_delta_zero_at_maximal_theta23"] = tm1_exact_check()
    inv = tm1_involution()
    res["tm1_involutions"] = [dict(matrix=g.astype(int).tolist(), weights=np.round(w, 6).tolist()) for g, w in inv
                              if np.allclose(w, [2 / 3, 1 / 6, 1 / 6])]
    res["involution_weight_types"] = sorted({tuple(np.round(w, 4)) for _, w in inv})
    res["single_flavon"] = [single_scan(False, (2, 3, 4), samples, rng), single_scan(True, (2, 4), samples, rng),
                            single_scan(False, (2, 3, 4, 6), samples, rng), single_scan(True, (2, 4, 6), samples, rng)]
    res["two_flavons"] = []
    for pe, pn in ((False, False), (False, True), (True, False), (True, True)):
        res["two_flavons"].append(pair_scan(pe, pn, (2, 3, 4), samples, rng, False))
        res["two_flavons"].append(pair_scan(pe, pn, (2, 4), samples, rng, True))
        res["two_flavons"].append(pair_scan(pe, pn, (2, 4, 6), samples // 2, rng, True))
    return res


def main():
    samples = int(sys.argv[1]) if len(sys.argv) > 1 else 120
    res = run(samples)
    OUT.write_text(json.dumps(res, indent=1))
    for s in res["sum_rules"]:
        print(s["column"], "s23 allowed:", s["s23_allowed"])
    print("TM1 involutions:", res["tm1_involutions"])
    print("weights:", res["involution_weight_types"])
    for s in res["single_flavon"]:
        print("single", s["rep"], s["degrees"], s["global_minimum_types"])
    for s in res["two_flavons"]:
        print("pair", s["phi_e"], s["phi_nu"], s["degrees"], s["shaping_Z2xZ2"], s["outcomes"], s["vev_types"])


if __name__ == "__main__":
    main()
