#!/usr/bin/env python3
"""Pass 11225: can the geometry's own finite symmetry fix the mixing angles?  An exhaustive residual-symmetry scan.

The residual-symmetry method of flavour physics: if the three generations carry a 3-dimensional representation rho of a
finite group, and the charged-lepton and neutrino mass terms are left invariant by subgroups G_e and G_nu, then the
mixing matrix is fixed by group theory alone: U_PMNS = U_e^dagger U_nu, with U_e the joint eigenbasis of rho(G_e) and
U_nu that of rho(G_nu).  For Majorana neutrinos G_nu is a Klein group Z2 x Z2 of involutions (full prediction) or a
single Z2 (one predicted column); for Dirac neutrinos any abelian G_nu with nondegenerate joint spectrum.  A Z2 on the
charged-lepton side predicts one row.

This pass applies the method to EVERY 3-dimensional irreducible representation of EVERY subgroup (up to conjugacy) of
  * W(E6) = Aut(W33), order 51840 (350 subgroup classes), and
  * Sp(4,3), the linear symmetry of two qutrits (its centre acts on the Weil representation),
exported by `analysis/gap/w33_pass11225_three_dim_images.g` (gzipped in data/).  Mixing patterns depend only on the image group in U(3) up
to unitary conjugacy, so images are deduplicated by (order, multiset of (trace, det, order)).  GAP's representations
are first made unitary (conjugation by the square root of the invariant Hermitian form): a first run without this step
reported "viable" full patterns for the qutrit Clifford group that were not even doubly stochastic, an artefact of a
non-orthogonal eigenbasis; every pattern is now asserted to be unistochastic.

Data (as quoted, NuFIT 6.0, normal ordering, 3 sigma): sin^2 th12 in [0.275, 0.345], sin^2 th13 in [0.02030, 0.02388],
sin^2 th23 in [0.430, 0.596], delta unconstrained at 3 sigma.  A pattern is viable when some assignment of generations
(row and column permutations) and some angles and phase inside that box reproduce it.  CKM: |V_us| in
[0.2230, 0.2271] (PDG 2024 0.22501 +- 0.00068, 3 sigma).
"""
from __future__ import annotations

import itertools
import json
import re
import sys
from cmath import exp, pi
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
GAP_OUT = ROOT / "data" / "w33_pass11225_three_dim_images.txt.gz"
OUT = ROOT / "data" / "w33_pass11225_flavour_from_symmetry.json"
BOX = dict(s12=(0.275, 0.345), s13=(0.02030, 0.02388), s23=(0.430, 0.596))
VUS = (0.2230, 0.2271)
LEGACY = dict(s12=4 / 13, s23=7 / 13, s13=2 / 91)
TOL = 1e-8


# ---------------------------------------------------------------- parsing GAP cyclotomics
def E(n):
    return exp(2j * pi / n)


def cyc(s):
    s = s.strip().replace("^", "**")
    return complex(eval(s, {"E": E, "__builtins__": {}}))


def parse(path=GAP_OUT):
    import gzip
    text = gzip.open(path, "rt").read().replace("\\\n", "")      # GAP's line continuations
    lines = []
    for raw in text.splitlines():                                      # GAP also wraps printed lines silently
        if raw.startswith(("M ", "REP ", "GROUP ")) or not lines:
            lines.append(raw)
        else:
            lines[-1] += raw
    reps, cur = [], None
    for line in lines:
        if line.startswith("REP "):
            meta = dict(kv.split("=", 1) for kv in line.split()[2:])
            meta["family"] = line.split()[1]
            cur = dict(meta=meta, mats=[])
            reps.append(cur)
        elif line.startswith("M ") and cur is not None:
            vals = [cyc(t) for t in line[2:].split(";")]
            cur["mats"].append(np.array(vals).reshape(3, 3))
    return reps


def unitarise(mats):
    """GAP's Dixon representations need not be unitary.  Conjugate by T = H^(1/2), H = sum_g rho(g)^dagger rho(g)
    (the invariant positive form), so that every matrix is unitary; eigenbases are then orthonormal."""
    H = sum(M.conj().T @ M for M in mats)
    w, W = np.linalg.eigh(H)
    T = W @ np.diag(np.sqrt(w)) @ W.conj().T
    Ti = W @ np.diag(1 / np.sqrt(w)) @ W.conj().T
    out = [T @ M @ Ti for M in mats]
    assert all(np.allclose(M.conj().T @ M, np.eye(3), atol=1e-9) for M in out), "unitarisation failed"
    return out


def fingerprint(mats):
    keys = []
    for M in mats:
        o = element_order(M)
        keys.append((round(np.trace(M).real, 6), round(np.trace(M).imag, 6),
                     round(np.linalg.det(M).real, 6), round(np.linalg.det(M).imag, 6), o))
    return (len(mats), tuple(sorted(keys)))


def element_order(M, cap=240):
    P = M.copy()
    for k in range(1, cap + 1):
        if np.allclose(P, np.eye(3), atol=1e-9):
            return k
        P = P @ M
    return -1


# ---------------------------------------------------------------- residual bases
def nondegenerate(M):
    ev = np.linalg.eigvals(M)
    return all(abs(ev[i] - ev[j]) > 1e-6 for i in range(3) for j in range(i + 1, 3))


def eigbasis(M):
    """orthonormal eigenbasis of a unitary matrix with distinct eigenvalues (columns)"""
    ev, V = np.linalg.eig(M)
    V = V / np.linalg.norm(V, axis=0)
    assert np.allclose(V.conj().T @ V, np.eye(3), atol=1e-8), "eigenbasis not orthonormal"
    return V


def basis_key(U):
    """a basis up to phases and order = its set of rank-one projectors"""
    return frozenset(tuple(np.round(np.outer(U[:, i], U[:, i].conj()).ravel(), 6)) for i in range(3))


def unique_bases(lst):
    seen, out = set(), []
    for kind, U in lst:
        k = basis_key(U)
        if k not in seen:
            seen.add(k)
            out.append((kind, U))
    return out


def full_bases(mats):
    """eigenbases of nondegenerate elements and joint eigenbases of Klein groups of involutions (deduplicated)"""
    bases, klein = [], []
    for M in mats:
        if nondegenerate(M):
            bases.append(("cyclic", eigbasis(M)))
    inv = [M for M in mats if np.allclose(M @ M, np.eye(3), atol=1e-9) and not np.allclose(M, np.eye(3), atol=1e-9)
           and not np.allclose(M, -np.eye(3), atol=1e-9)]
    for a, b in itertools.combinations(range(len(inv)), 2):
        A, B = inv[a], inv[b]
        if np.allclose(A @ B, B @ A, atol=1e-9) and not np.allclose(A @ B, np.eye(3), atol=1e-9):
            C = A + np.pi * B
            if nondegenerate(C):
                klein.append(("klein", eigbasis(C)))
    return unique_bases(bases), unique_bases(klein), inv


def z2_vectors(inv):
    """the singled-out eigenvector of each involution with eigenvalues {1,1,-1} or {-1,-1,1} (deduplicated)"""
    out, seen = [], set()
    for M in inv:
        ev, V = np.linalg.eig(M)
        r = np.round(ev.real).astype(int)
        for s in (1, -1):
            if (r == s).sum() == 1:
                v = V[:, list(r).index(s)]
                v = v / np.linalg.norm(v)
                k = tuple(np.round(np.outer(v, v.conj()).ravel(), 6))
                if k not in seen:
                    seen.add(k)
                    out.append(v)
    return out


def dedupe_mags(mats_list, perm_invariant=True):
    pre, uniq = set(), []
    for P in mats_list:
        k = tuple(sorted(np.round(P, 6).ravel()))
        k2 = (k, tuple(sorted(tuple(sorted(np.round(r, 6))) for r in P)))
        if k2 not in pre:
            pre.add(k2)
            uniq.append(P)
    mats_list = uniq
    seen, out = set(), []
    for P in mats_list:
        key = canonical(P) if perm_invariant else tuple(np.round(P, 6).ravel())
        if key not in seen:
            seen.add(key)
            out.append(P)
    return out


def canonical(P):
    best = None
    for r in itertools.permutations(range(P.shape[0])):
        for c in itertools.permutations(range(P.shape[1])):
            key = tuple(np.round(P[np.ix_(r, c)], 6).ravel())
            if best is None or key < best:
                best = key
    return best


def canonical_vec(v):
    return tuple(sorted(np.round(v, 6)))


# ---------------------------------------------------------------- comparison with data
def pmns(s12, s13, s23, d):
    c12, c13, c23 = np.sqrt(1 - s12), np.sqrt(1 - s13), np.sqrt(1 - s23)
    a12, a13, a23 = np.sqrt(s12), np.sqrt(s13), np.sqrt(s23)
    e = np.exp(1j * d)
    return np.array([[c12 * c13, a12 * c13, a13 * np.conj(e)],
                     [-a12 * c23 - c12 * a23 * a13 * e, c12 * c23 - a12 * a23 * a13 * e, a23 * c13],
                     [a12 * a23 - c12 * c23 * a13 * e, -c12 * a23 - a12 * c23 * a13 * e, c23 * c13]])


def in_box(s12, s13, s23):
    return (BOX["s12"][0] <= s12 <= BOX["s12"][1] and BOX["s13"][0] <= s13 <= BOX["s13"][1]
            and BOX["s23"][0] <= s23 <= BOX["s23"][1])


def full_viable(P):
    """P = |U|^2 (3x3).  Viable iff some row/column permutation has standard-parametrisation angles in the box
    (any unistochastic P comes from some unitary, and delta is free at 3 sigma)"""
    hits = []
    for r in itertools.permutations(range(3)):
        for c in itertools.permutations(range(3)):
            Q = P[np.ix_(r, c)]
            s13 = Q[0, 2]
            if s13 >= 1 - 1e-12:
                continue
            s12, s23 = Q[0, 1] / (1 - s13), Q[1, 2] / (1 - s13)
            if in_box(s12, s13, s23):
                hits.append((r, c, (s12, s13, s23)))
    return hits


def column_viable(col):
    """col = predicted |U_{alpha i}|^2 for one column (up to row permutation): is there a point of the box whose
    column i matches?  Returns the best residual over i, row permutations and the box."""
    best = (np.inf, None)
    lo = np.array([BOX["s12"][0], BOX["s13"][0], BOX["s23"][0], 0.0])
    hi = np.array([BOX["s12"][1], BOX["s13"][1], BOX["s23"][1], 2 * np.pi])
    for i in range(3):
        for perm in set(itertools.permutations(np.round(col, 12))):
            target = np.array(perm)

            def f(x):
                U = pmns(*x)
                return float(np.sum((np.abs(U[:, i]) ** 2 - target) ** 2))
            for start in np.linspace(0.05, 0.95, 4):
                x0 = lo + start * (hi - lo)
                r = minimize(f, x0, bounds=list(zip(lo, hi)), method="L-BFGS-B")
                if r.fun < best[0]:
                    best = (r.fun, dict(column=i, target=[float(t) for t in target], angles=[float(t) for t in r.x]))
    return best


def row_viable(row):
    best = (np.inf, None)
    lo = np.array([BOX["s12"][0], BOX["s13"][0], BOX["s23"][0], 0.0])
    hi = np.array([BOX["s12"][1], BOX["s13"][1], BOX["s23"][1], 2 * np.pi])
    for a in range(3):
        for perm in set(itertools.permutations(np.round(row, 12))):
            target = np.array(perm)

            def f(x):
                U = pmns(*x)
                return float(np.sum((np.abs(U[a, :]) ** 2 - target) ** 2))
            for start in np.linspace(0.05, 0.95, 4):
                x0 = lo + start * (hi - lo)
                r = minimize(f, x0, bounds=list(zip(lo, hi)), method="L-BFGS-B")
                if r.fun < best[0]:
                    best = (r.fun, dict(row=a, target=[float(t) for t in target], angles=[float(t) for t in r.x]))
    return best


# ---------------------------------------------------------------- the scan
def analyse_image(mats):
    bases, klein, inv = full_bases(mats)
    allfull = bases + klein
    z2 = z2_vectors(inv)
    res = dict(order=len(mats), nondegenerate_elements=len(bases), klein_bases=len(klein), involutions=len(inv))
    # full patterns: charged leptons any full basis; Majorana nu: Klein; Dirac nu: any full basis
    maj, dirac = [], []
    for _, Ue in allfull:
        for _, Un in klein:
            maj.append(np.abs(Ue.conj().T @ Un) ** 2)
        for _, Un in allfull:
            dirac.append(np.abs(Ue.conj().T @ Un) ** 2)
    for P in maj + dirac:
        assert np.allclose(P.sum(0), 1, atol=1e-8) and np.allclose(P.sum(1), 1, atol=1e-8), "pattern not unistochastic"
    maj, dirac = dedupe_mags(maj), dedupe_mags(dirac)
    res["majorana_full_patterns"] = len(maj)
    res["dirac_full_patterns"] = len(dirac)
    res["full_viable"] = [dict(kind=k, P=np.round(P, 6).tolist(), hits=len(full_viable(P)))
                          for k, lst in (("majorana", maj), ("dirac", dirac)) for P in lst if full_viable(P)]
    # one predicted column (Majorana Z2 on the neutrino side) and one predicted row (Z2 on the charged-lepton side)
    cols = {}
    for _, Ue in allfull:
        for v in z2:
            c = np.abs(Ue.conj().T @ v) ** 2
            assert abs(c.sum() - 1) < 1e-8
            cols[canonical_vec(c)] = c
    rows = {}
    for v in z2:
        for _, Un in klein + bases:
            r = np.abs(v.conj() @ Un) ** 2
            assert abs(r.sum() - 1) < 1e-8
            rows[canonical_vec(r)] = r
    res["columns"] = []
    for key, c in cols.items():
        val, info = column_viable(c)
        res["columns"].append(dict(column=[float(x) for x in key], viable=bool(val < 1e-10), residual=val, fit=info))
    res["rows"] = []
    for key, r in rows.items():
        val, info = row_viable(r)
        res["rows"].append(dict(row=[float(x) for x in key], viable=bool(val < 1e-10), residual=val, fit=info))
    # CKM: Dirac on both sides; Cabibbo-only and full patterns
    cab = []
    for P in dirac:
        V = np.sqrt(np.clip(P, 0, 1))
        for r in itertools.permutations(range(3)):
            for c in itertools.permutations(range(3)):
                Q = V[np.ix_(r, c)]
                if VUS[0] <= Q[0, 1] <= VUS[1] and Q[0, 0] > 0.9 and Q[2, 2] > 0.9:
                    cab.append(np.round(Q, 6).tolist())
    res["cabibbo_patterns"] = cab[:5]
    res["cabibbo_any"] = bool(cab)
    # the legacy numerology: is (4/13, 2/91, 7/13) the angle set of any full pattern?
    leg = []
    for P in maj + dirac:
        for r in itertools.permutations(range(3)):
            for c in itertools.permutations(range(3)):
                Q = P[np.ix_(r, c)]
                s13 = Q[0, 2]
                if s13 < 1 - 1e-9:
                    s12, s23 = Q[0, 1] / (1 - s13), Q[1, 2] / (1 - s13)
                    if abs(s12 - LEGACY["s12"]) < 1e-6 and abs(s13 - LEGACY["s13"]) < 1e-6 \
                            and abs(s23 - LEGACY["s23"]) < 1e-6:
                        leg.append(np.round(Q, 6).tolist())
    res["legacy_pattern_found"] = bool(leg)
    return res


def column_predictions(columns):
    """for every viable column, every PMNS column index i and every assignment of its entries to (e, mu, tau): the
    predicted sin^2 th12 over the th13 box (i = 1, 2; for i = 3 the e-entry is sin^2 th13 itself) and the range of
    cos(delta) that the mu-entry forces over the th13 x th23 box.  Only assignments with a solution are kept."""
    out = []
    s13 = np.linspace(*BOX["s13"], 25)
    s23 = np.linspace(*BOX["s23"], 25)
    deltas = np.linspace(0, np.pi, 1441)
    for col in columns:
        for i in range(3):
            for perm in sorted(set(itertools.permutations(col))):
                e, mu, tau = perm
                if i == 0:
                    s12 = 1 - e / (1 - s13)
                elif i == 1:
                    s12 = e / (1 - s13)
                else:
                    if not BOX["s13"][0] <= e <= BOX["s13"][1]:
                        continue
                    s12 = np.linspace(*BOX["s12"], 25)
                ok12 = (s12 >= BOX["s12"][0]) & (s12 <= BOX["s12"][1])
                if not ok12.any():
                    continue
                cosd = []
                pairs = zip(s13, s12) if i < 2 else ((e, x) for x in s12)
                for a13, a12 in pairs:
                    if not BOX["s12"][0] <= a12 <= BOX["s12"][1]:
                        continue
                    for a23 in s23:
                        vals = np.array([abs(pmns(a12, a13, a23, d)[1, i]) ** 2 for d in deltas])
                        k = int(np.argmin(np.abs(vals - mu)))
                        if abs(vals[k] - mu) < 1e-3:
                            cosd.append(float(np.cos(deltas[k])))
                if cosd:
                    out.append(dict(column=[float(x) for x in col], index=i + 1, assignment=dict(e=e, mu=mu, tau=tau),
                                    sin2_th12=[float(np.min(s12[ok12])), float(np.max(s12[ok12]))] if i < 2 else None,
                                    cos_delta=[min(cosd), max(cosd)]))
    return out


def run():
    reps = parse()
    for rep in reps:
        rep["mats"] = unitarise(rep["mats"])
    images = {}
    for rep in reps:
        fp = fingerprint(rep["mats"])
        images.setdefault(fp, dict(mats=rep["mats"], sources=[]))["sources"].append(
            {k: rep["meta"][k] for k in ("family", "sub", "order", "struct", "irr", "image_struct")})
    out = dict(pass_id=11225, representations=len(reps), distinct_images=len(images), images=[])
    for fp, im in sorted(images.items(), key=lambda kv: kv[0][0]):
        a = analyse_image(im["mats"])
        a["image_struct"] = im["sources"][0]["image_struct"]
        a["sources"] = im["sources"][:6]
        a["n_sources"] = len(im["sources"])
        out["images"].append(a)
        print(a["image_struct"], a["order"], "full viable:", len(a["full_viable"]),
              "viable columns:", [c["column"] for c in a["columns"] if c["viable"]],
              "viable rows:", [r["row"] for r in a["rows"] if r["viable"]],
              "cabibbo:", a["cabibbo_any"], "legacy:", a["legacy_pattern_found"], flush=True)
    out["any_full_viable"] = any(im["full_viable"] for im in out["images"])
    out["viable_columns"] = sorted({tuple(c["column"]) for im in out["images"] for c in im["columns"] if c["viable"]})
    out["viable_rows"] = sorted({tuple(r["row"]) for im in out["images"] for r in im["rows"] if r["viable"]})
    out["any_cabibbo"] = any(im["cabibbo_any"] for im in out["images"])
    out["legacy_numerology_realised"] = any(im["legacy_pattern_found"] for im in out["images"])
    out["column_predictions"] = column_predictions([np.array(c) for c in out["viable_columns"]])
    return out


def main():
    res = run()
    OUT.write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps({k: v for k, v in res.items() if k != "images"}, indent=1, default=float))


if __name__ == "__main__":
    main()
