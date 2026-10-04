"""Pass 11370: the witness-degree law 10/6/4 (Pass 11357) is a twisted Frobenius-Schur count.

Pass 11355: U is substrate-reversible iff U^T is Clifford-conjugate to U, and J_2t vanishes identically iff every
Clifford-conjugation-invariant polynomial of degree (t,t) is invariant under tau: f(U) -> f(U^T).

Peter-Weyl.  Degree-(t,t) polynomials on PU(d) are spanned by the matrix coefficients U -> tr(M lambda(U)) of the
rational irreps lambda = (alpha, beta), |alpha| = |beta| = k <= t, l(alpha) + l(beta) <= d.  In a real basis
lambda(U^T) = lambda(U)^T, so tau acts on coefficients by M -> M^T and Ad_C by M -> lambda(C) M lambda(C)^dag.  The
trace of M -> P M^T Q on End(V) is tr(P Q^T), hence

    #(invariants in lambda)          = (1/|Cl|) sum_C |chi_lambda(C)|^2,
    #(tau-even) - #(tau-odd)         = (1/|Cl|) sum_C chi_lambda(conj(C) C)       (twisted Frobenius-Schur indicator of
                                                                                   Kawanaka-Matsuyama 1990, twist = K)
    n_lambda^-  = #(tau-odd invariants) = (1/2) [ first - second ].

THEOREM (exact reduction).  J_2t is not identically zero  <=>  some lambda with |alpha| = |beta| <= t has n^- > 0.
Restricting lambda to the Clifford group, lambda = sum_rho m_rho rho, gives n^- = (1/2) sum_rho m_rho (m_rho - nu_K(rho))
with nu_K in {-1, 0, 1}: odd invariants appear when a Clifford irrep repeats, or appears with nu_K != 1.

Qubits.  conj(C) = Y C Y up to phase, so the twisted indicator is the ordinary Frobenius-Schur indicator of the octahedral
group (all irreps real):  n^- = sum_rho m_rho (m_rho - 1) / 2.  Spin j restricted to O first repeats an irrep at j = 5
(11 = E + 2 T1 + T2), so the first qubit witness has degree 2 * 5 = 10.

Characters by Jacobi-Trudi (h_k from power sums); the counts are checked to be integers, and the predicted first
degrees are compared with Pass 11357 (d = 2, 3, 5) and PREDICTED, then tested by polynomial identity, for d = 7.
"""

from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11357_jarlskog_degree_by_dimension as P7  # noqa: E402

OUT = ROOT / "data" / "w33_pass11370_twisted_indicator_degree_law.json"


def partitions(k, maxlen):
    def rec(n, m, ln):
        if n == 0:
            yield ()
            return
        if ln == 0:
            return
        for p in range(min(n, m), 0, -1):
            for rest in rec(n - p, p, ln - 1):
                yield (p,) + rest
    return list(rec(k, k, maxlen))


def mixed_weights(t, d):
    """highest weights (alpha, 0..0, -beta reversed) of the PU(d) irreps in V^(x)t (x) Vbar^(x)t"""
    out = []
    for k in range(t + 1):
        for a in partitions(k, d):
            for b in partitions(k, d - len(a)):
                w = list(a) + [0] * (d - len(a) - len(b)) + [-x for x in reversed(b)]
                out.append((k, a, b, tuple(w)))
    return out


def h_from_eigs(x, kmax):
    """complete homogeneous symmetric polynomials h_0..h_kmax of the eigenvalues x (batched: x has shape (N, d))"""
    p = [None] + [np.sum(x ** i, axis=1) for i in range(1, kmax + 1)]
    h = [np.ones(x.shape[0], dtype=complex)]
    for k in range(1, kmax + 1):
        h.append(sum(p[i] * h[k - i] for i in range(1, k + 1)) / k)
    return h


def characters(eigs, w):
    """chi_w at elements with eigenvalues eigs (N, d): s_{w + m}(x) / det(x)^m, m = -min(w)"""
    m = -min(w)
    lam = [wi + m for wi in w]
    ell = sum(1 for v in lam if v > 0)
    if ell == 0:
        return np.ones(eigs.shape[0], dtype=complex)
    h = h_from_eigs(eigs, lam[0] + ell)
    N = eigs.shape[0]
    Mt = np.zeros((N, ell, ell), dtype=complex)
    for i in range(ell):
        for j in range(ell):
            idx = lam[i] - i + j
            Mt[:, i, j] = h[idx] if 0 <= idx < len(h) else 0
    return np.linalg.det(Mt) / np.prod(eigs, axis=1) ** m


def counts(Cl, t):
    eC = np.linalg.eigvals(Cl)
    eK = np.linalg.eigvals(np.einsum('cij,cjk->cik', Cl.conj(), Cl))
    rows = []
    for k, a, b, w in mixed_weights(t, Cl.shape[1]):
        ch = characters(eC, w)
        tot = float(np.mean(np.abs(ch) ** 2))
        tw = complex(np.mean(characters(eK, w)))
        nm = (tot - tw.real) / 2
        assert abs(tw.imag) < 1e-8 and abs(tot - round(tot)) < 1e-7 and abs(nm - round(nm)) < 1e-7, (w, tot, tw)
        rows.append(dict(k=k, alpha=a, beta=b, invariants=round(tot), twisted=round(tw.real), odd=round(nm)))
    return rows


def first_degree(Cl, tmax):
    for t in range(1, tmax + 1):
        rows = [r for r in counts(Cl, t) if r["k"] == t]
        odd = [r for r in rows if r["odd"]]
        if odd:
            return t, odd
    return None, []


def octahedral_spin_decomposition():
    """qubit check: spin j restricted to the octahedral rotation group (character table, classes 1, 8C3, 3C2, 6C4, 6C2')"""
    sizes = np.array([1, 8, 3, 6, 6])
    angles = [0, 2 * np.pi / 3, np.pi, np.pi / 2, np.pi]
    table = {"A1": [1, 1, 1, 1, 1], "A2": [1, 1, 1, -1, -1], "E": [2, -1, 2, 0, 0],
             "T1": [3, 0, -1, 1, -1], "T2": [3, 0, -1, -1, 1]}
    out = {}
    for j in range(0, 8):
        chi = np.array([2 * j + 1 if a == 0 else np.sin((j + 0.5) * a) / np.sin(a / 2) for a in angles])
        mult = {r: int(round(float(np.sum(sizes * chi * np.array(v)) / 24))) for r, v in table.items()}
        out[j] = dict(mult={r: m for r, m in mult.items() if m}, odd=sum(m * (m - 1) // 2 for m in mult.values()))
    return out


def identity_test(d, Cl, ts, samples=12, seed=0):
    rng = np.random.default_rng(seed)
    best = {t: 0.0 for t in ts}
    for _ in range(samples):
        Z = (rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))) / np.sqrt(2)
        Q, R = np.linalg.qr(Z)
        U = Q * (np.diag(R) / np.abs(np.diag(R)))
        for t in ts:
            best[t] = max(best[t], abs(P7.J(U, t, Cl)))
    return best


def run():
    res = dict(pass_id=11370)
    cert = json.load(open(ROOT / "data" / "w33_pass11357_jarlskog_degree_by_dimension.json"))
    res["octahedral_spin_restrictions"] = octahedral_spin_decomposition()
    table = {}
    for d, tmax in ((2, 6), (3, 4), (5, 3), (7, 3)):
        Cl = P7.clifford_group(d)
        t, odd = first_degree(Cl, tmax)
        per_t = {}
        for tt in range(1, (t or tmax) + 1):
            rows = counts(Cl, tt)
            per_t[tt] = dict(invariants=sum(r["invariants"] for r in rows), odd=sum(r["odd"] for r in rows))
        table[d] = dict(clifford_order=len(Cl), first_odd_degree_t=t, witness_degree=2 * t if t else None,
                        odd_irreps_at_first_degree=[dict(alpha=r["alpha"], beta=r["beta"], invariants=r["invariants"],
                                                         twisted=r["twisted"], odd=r["odd"]) for r in odd],
                        cumulative=per_t)
        print(d, table[d]["witness_degree"], per_t, flush=True)
    res["degree_law"] = table
    # comparison with Pass 11357's polynomial identity tests
    res["matches_pass_11357"] = all(table[d]["witness_degree"] == {2: 10, 3: 6, 5: 4}[d] for d in (2, 3, 5))
    # d = 7: prediction above, now tested by identity testing (Clifford group of order 7^2 * 336 = 16464)
    Cl7 = P7.clifford_group(7)
    t7 = table[7]["first_odd_degree_t"]
    best = identity_test(7, Cl7, list(range(1, t7 + 1)), samples=8)
    res["d7_identity_test_max_abs_J"] = {f"J{2 * t}": best[t] for t in best}
    res["d7_prediction_confirmed"] = all(best[t] < 1e-10 for t in range(1, t7)) and best[t7] > 1e-6
    _ = cert
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in res.items() if k != "octahedral_spin_restrictions"}, indent=1, default=str))


if __name__ == "__main__":
    main()
