"""Pass 11236: the cubic (magic) gate as a residual flavour symmetry -- full mixing matrices including delta.

The parallel session's Pass 11225 scanned every residual-symmetry pattern available inside the geometry's FINITE
symmetry groups (W(E6), Sp(4,3) and their subgroups, which include the single-qutrit Clifford group): no complete PMNS
pattern and no Cabibbo angle; single columns survive (TM1; an order-9 Clifford column with delta ~ pi).  The cubic gate
T = diag(1, zeta9, zeta9^8) is not in those groups -- it is the non-Clifford resource (Passes 11213, 11227) -- so the
patterns it generates were never tested.

Scan (one qutrit = three generations):
  * charged-lepton residual group G_e = <C T C^dag>, C in the 216 Cliffords (T has non-degenerate spectrum 1, zeta9,
    zeta9^8, so U_e = eigenbasis of C T C^dag);
  * neutrino residual group G_nu = <g>, g any single-qutrit Clifford with non-degenerate spectrum (Dirac neutrinos:
    every abelian residual group with non-degenerate spectrum fixes U_nu), or a Klein group of Clifford involutions
    (Majorana);
  * U = U_e^dag U_nu up to row and column permutations; extract sin^2 theta_12, theta_13, theta_23, J and delta and
    compare with NuFIT 6.0 (3 sigma, normal ordering, as quoted in Pass 11225) and with |V_us|.
"""

from __future__ import annotations

import itertools
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11213_cubic_t_violation as P  # noqa: E402

OUT = ROOT / "data" / "w33_pass11236_cubic_flavour.json"
S12, S13, S23 = (0.275, 0.345), (0.02030, 0.02388), (0.430, 0.596)
VUS = (0.2230, 0.2271)


def eigbasis(M):
    w, v = np.linalg.eig(M)
    if min(abs(w[i] - w[j]) for i in range(3) for j in range(i + 1, 3)) < 1e-6:
        return None
    q, _ = np.linalg.qr(v)          # orthonormalise (eigenvectors of a unitary are orthogonal already)
    return v / np.linalg.norm(v, axis=0)


def angles(U):
    a = np.abs(U) ** 2
    s13 = a[0, 2]
    if s13 > 1 - 1e-12:
        return None
    s12 = a[0, 1] / (1 - s13)
    s23 = a[1, 2] / (1 - s13)
    J = float(np.imag(U[0, 0] * U[1, 1] * np.conj(U[0, 1]) * np.conj(U[1, 0])))
    c12, c13, c23 = np.sqrt(1 - s12), np.sqrt(1 - s13), np.sqrt(1 - s23)
    den = np.sqrt(s12) * c12 * np.sqrt(s23) * c23 * c13 ** 2 * np.sqrt(s13)
    sind = J / den if den > 1e-12 else 0.0
    return float(s12), float(s13), float(s23), J, float(np.clip(sind, -1, 1))


def finite_order(W, kmax=60):
    M = np.eye(3, dtype=complex)
    for k in range(1, kmax + 1):
        M = M @ W
        if np.allclose(M / M[0, 0] if abs(M[0, 0]) > 1e-9 else M, np.eye(3), atol=1e-8) and abs(M[0, 0]) > 1e-9:
            return k
    return None


def basis_key(V):
    """projector set of an eigenbasis, rounded (phase- and order-free)"""
    return tuple(sorted(tuple(np.round(np.outer(V[:, i], V[:, i].conj()).ravel(), 5).view(float)) for i in range(3)))


def depth1_bases(C):
    """eigenbases of finite-order one-cubic-gate words C1 T C2 (non-degenerate spectrum), deduplicated"""
    out, orders = {}, Counter()
    for c1 in C:
        for c2 in C:
            W = c1 @ P.T1 @ c2
            o = finite_order(W)
            if o is None:
                continue
            V = eigbasis(W)
            if V is None:
                continue
            k = basis_key(V)
            if k not in out:
                out[k] = V
                orders[o] += 1
    return list(out.values()), orders


def scan_pairs(E, N):
    hits, cab, mods = [], [], Counter()
    seen = set()
    for Ve in E:
        for Vn in N:
            U = Ve.conj().T @ Vn
            mods[tuple(sorted(np.round((np.abs(U) ** 2).ravel(), 4)))] += 1
            for rp in itertools.permutations(range(3)):
                for cp in itertools.permutations(range(3)):
                    W = U[np.ix_(rp, cp)]
                    a = angles(W)
                    if a is None or not all(np.isfinite(a)):
                        continue
                    s12, s13, s23, J, sd = a
                    key = (round(s12, 4), round(s13, 4), round(s23, 4), round(abs(J), 4))
                    if key in seen:
                        continue
                    seen.add(key)
                    if S12[0] <= s12 <= S12[1] and S13[0] <= s13 <= S13[1] and S23[0] <= s23 <= S23[1]:
                        hits.append(dict(s12=s12, s13=s13, s23=s23, J=J, sin_delta=sd))
                    vus = abs(W[0, 1])
                    if VUS[0] <= vus <= VUS[1] and abs(W[0, 2]) < 0.02 and abs(W[1, 2]) < 0.06:
                        cab.append(dict(Vus=float(vus), Vub=float(abs(W[0, 2])), Vcb=float(abs(W[1, 2]))))
    return hits, cab, mods


def run():
    C = P.clifford1()
    res = dict(pass_id=11236)
    # neutrino side: Clifford elements with non-degenerate spectrum
    nu = []
    for g in C:
        V = eigbasis(g)
        if V is not None:
            nu.append(V)
    # charged-lepton side: Clifford conjugates of T
    e = []
    for c in C:
        V = eigbasis(c @ P.T1 @ c.conj().T)
        if V is not None:
            e.append(V)
    res["n_neutrino_symmetries"] = len(nu)
    res["n_charged_symmetries"] = len(e)
    pmns_hits, cab_hits, patterns = [], [], Counter()
    seen = set()
    for Ve in e:
        for Vn in nu:
            U = Ve.conj().T @ Vn
            for rp in itertools.permutations(range(3)):
                for cp in itertools.permutations(range(3)):
                    W = U[np.ix_(rp, cp)]
                    a = angles(W)
                    if a is None:
                        continue
                    s12, s13, s23, J, sd = a
                    key = (round(s12, 5), round(s13, 5), round(s23, 5), round(abs(J), 5))
                    if key in seen:
                        continue
                    seen.add(key)
                    patterns[key] += 1
                    if S12[0] <= s12 <= S12[1] and S13[0] <= s13 <= S13[1] and S23[0] <= s23 <= S23[1]:
                        pmns_hits.append(dict(s12=s12, s13=s13, s23=s23, J=J, sin_delta=sd))
                    vus = abs(W[0, 1])
                    if VUS[0] <= vus <= VUS[1] and abs(W[0, 2]) < 0.02 and abs(W[1, 2]) < 0.06:
                        cab_hits.append(dict(Vus=float(vus), Vub=float(abs(W[0, 2])), Vcb=float(abs(W[1, 2]))))
    res["distinct_patterns"] = len(patterns)
    res["pmns_complete_hits"] = pmns_hits[:20]
    res["n_pmns_complete_hits"] = len(pmns_hits)
    res["cabibbo_hits"] = cab_hits[:20]
    res["n_cabibbo_hits"] = len(cab_hits)
    # which |U|^2 matrices arise at all (rounded, sorted)
    mods = Counter()
    for Ve in e:
        for Vn in nu:
            a = np.abs(Ve.conj().T @ Vn) ** 2
            mods[tuple(sorted(np.round(a.ravel(), 4)))] += 1
    res["distinct_modulus_patterns"] = [list(k) for k in mods][:30]
    res["n_distinct_modulus_patterns"] = len(mods)
    # depth-1 magic words: genuinely non-stabilizer eigenbases
    D1, orders = depth1_bases(C)
    res["depth1_finite_order_bases"] = len(D1)
    res["depth1_orders"] = {str(k): v for k, v in sorted(orders.items())}
    cliff = []
    seenb = set()
    for V in nu:
        k = basis_key(V)
        if k not in seenb:
            seenb.add(k)
            cliff.append(V)
    res["clifford_bases"] = len(cliff)
    for name, (E, N) in {"magic_e_x_clifford_nu": (D1, cliff), "clifford_e_x_magic_nu": (cliff, D1),
                         "magic_e_x_magic_nu": (D1, D1)}.items():
        hits, cab, mods = scan_pairs(E, N)
        res[name] = dict(pmns_hits=len(hits), examples=hits[:10], cabibbo_hits=len(cab), cab_examples=cab[:5],
                         modulus_patterns=len(mods))
        print(name, len(hits), len(cab), len(mods), flush=True)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=float)
    print(json.dumps({k: v for k, v in res.items() if k != "distinct_modulus_patterns"}, indent=1, default=float))
    print("modulus patterns:", res["n_distinct_modulus_patterns"])
    for m in res["distinct_modulus_patterns"][:12]:
        print(m)


if __name__ == "__main__":
    main()
