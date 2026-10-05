"""Pass 11491: the Hilbert-Molien series of time-even and time-odd Clifford invariants of a qutrit state, exactly.

A ray invariant of degree k is a polynomial p(psi, psibar) of bidegree (k, k) fixed by the 216 projective Cliffords.  Time
reversal (complex conjugation, composed with any Clifford) splits them into even and odd ones.  With h_k the complete
homogeneous symmetric polynomial,

    D_k  = (1/216) sum_g |h_k(eig g)|^2          (all invariants)
    Tw_k = (1/216) sum_g h_k(eig(conj(g) g))     (trace of time reversal on them)
    E_k = (D_k + Tw_k)/2,  O_k = (D_k - Tw_k)/2.

Everything is computed exactly: eigenvalues are identified as roots of unity of order dividing 72 and every sum is reduced
modulo the 72nd cyclotomic polynomial over the integers.  The four closed forms below are then PROVED, not fitted: both sides
are quasi-polynomials in k of period P (computed) and degree <= 4, so agreement on 5P consecutive values is equality.
"""

from __future__ import annotations

import json
import sys
from math import lcm
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11357_jarlskog_degree_by_dimension as P7  # noqa: E402

OUT = ROOT / "data" / "w33_pass11491_time_odd_molien.json"
N = 72
KMAX = 72
t = sp.symbols("t")

# closed forms (numerators over the common primary denominator (1-t)(1-t^2)(1-t^3)(1-t^4)(1-t^6))
PRIM = (1 - t) * (1 - t**2) * (1 - t**3) * (1 - t**4) * (1 - t**6)
FORMS = {
    "odd": t**6 * (1 + t) / PRIM,
    "even": sp.expand((1 - t + t**5 - t**6 + t**7 - t**11 + t**12) * (1 + t)) / PRIM,
    "total": sp.expand((1 - t + t**5 + t**7 - t**11 + t**12) * (1 + t)) / PRIM,
    "twisted": (1 + t**5) / ((1 - t) * (1 - t**3) * (1 - t**4)),
}


def exponents(vals):
    """Exact exponents e with vals = exp(2 pi i e / N); asserts the identification."""
    e = np.rint(np.angle(vals) / (2 * np.pi) * N).astype(np.int64) % N
    assert np.abs(np.exp(2j * np.pi * e / N) - vals).max() < 1e-9 and np.abs(np.abs(vals) - 1).max() < 1e-9
    return e


def h_counts(e, k):
    """N_k(r) = #{(a,b,c) : a+b+c = k, a e0 + b e1 + c e2 = r mod N}, so h_k(zeta^e) = sum_r N_k(r) zeta^r."""
    a, b = np.meshgrid(np.arange(k + 1), np.arange(k + 1), indexing="ij")
    m = a + b <= k
    r = (a[m] * e[0] + b[m] * e[1] + (k - a[m] - b[m]) * e[2]) % N
    return np.bincount(r, minlength=N)


def reduce_exact(coeffs):
    """sum_d coeffs[d] zeta_N^d as an exact element of Z[zeta_N]; returns the integer if it is rational."""
    x = sp.symbols("x")
    poly = sp.Poly(sum(int(c) * x**d for d, c in enumerate(coeffs) if c), x)
    rem = poly.rem(sp.Poly(sp.cyclotomic_poly(N, x), x))
    cs = rem.all_coeffs()
    assert all(c == 0 for c in cs[:-1]), "not rational"
    return int(cs[-1]) if cs else 0


def series(form, K):
    num, den = sp.fraction(sp.together(form))
    nc = sp.Poly(sp.expand(num), t).all_coeffs()[::-1]
    dc = sp.Poly(sp.expand(den), t).all_coeffs()[::-1]
    assert dc[0] in (1, -1)
    out = []
    for k in range(K):
        s = (nc[k] if k < len(nc) else 0) - sum(dc[j] * out[k - j] for j in range(1, min(k, len(dc) - 1) + 1))
        out.append(int(s / dc[0]))
    return out


def run():
    G = P7.clifford_group(3)
    assert len(G) == 216
    ev_ratio, ev_twist = [], []
    for g in G:
        lam = np.linalg.eigvals(g)
        ev_ratio.append(exponents(lam / lam[0]))                  # D depends only on ratios (projective)
        ev_twist.append(exponents(np.linalg.eigvals(np.conj(g) @ g)))  # conj(g) g is phase-independent
    period_D = lcm(*[N // np.gcd.reduce(np.append(e, N)) for e in ev_ratio])
    period_T = lcm(*[N // np.gcd.reduce(np.append(e, N)) for e in ev_twist])
    D, Tw = [], []
    for k in range(KMAX + 1):
        cD = np.zeros(N, dtype=object)
        cT = np.zeros(N, dtype=object)
        for er, et in zip(ev_ratio, ev_twist):
            n = h_counts(er, k).astype(object)
            # |h|^2 = sum_{r,s} n_r n_s zeta^(r-s)
            for r in np.nonzero(n)[0]:
                cD += n[r] * np.roll(n[::-1], r + 1)                # index d = r - s mod N
            cT += h_counts(et, k).astype(object)
        sD, sT = reduce_exact(cD), reduce_exact(cT)
        assert sD % 216 == 0 and sT % 216 == 0
        D.append(sD // 216)
        Tw.append(sT // 216)
    assert all((d + w) % 2 == 0 for d, w in zip(D, Tw))
    E = [(d + w) // 2 for d, w in zip(D, Tw)]
    O = [(d - w) // 2 for d, w in zip(D, Tw)]
    data = {"total": D, "twisted": Tw, "even": E, "odd": O}
    # proof: quasi-polynomials of period P and degree <= deg agree everywhere once they agree on (deg+1) P values
    P = lcm(period_D, period_T, 12)
    need = 5 * P
    assert KMAX + 1 >= need, (KMAX, need)
    proved = {}
    for name, form in FORMS.items():
        den_roots_ok = all(sp.Pow(r, P).equals(1) for r in sp.roots(sp.fraction(sp.together(form))[1], t))
        proper = sp.degree(sp.fraction(sp.together(form))[0], t) < sp.degree(sp.fraction(sp.together(form))[1], t)
        proved[name] = bool(series(form, KMAX + 1) == data[name] and den_roots_ok and proper)
    res = dict(
        pass_id=11491,
        group_order=len(G),
        periods=dict(ratio_eigenvalues=period_D, twisted_eigenvalues=period_T, used=P, values_needed=need,
                     values_checked=KMAX + 1),
        sequences={k: v[:41] for k, v in data.items()},
        closed_forms={k: str(sp.factor(v)) for k, v in FORMS.items()},
        closed_forms_display=dict(
            odd="t^6 (1+t) / ((1-t)(1-t^2)(1-t^3)(1-t^4)(1-t^6))",
            twisted="(1+t^5) / ((1-t)(1-t^3)(1-t^4))",
        ),
        proved=proved,
        lowest_odd_degree=next(k for k, v in enumerate(O) if v),
        odd_not_multiples_of_h6=dict(
            O=O[6:16], h6_times_even=E[:10],
            excess=[O[6 + j] - E[j] for j in range(10)],
            first_new_odd_degree=6 + next(j for j in range(1, 20) if O[6 + j] > E[j]),
        ),
        design_check=dict(D1=D[1], D2=D[2], D3=D[3]),
    )
    print(json.dumps(res, indent=1), flush=True)
    return res


def explicit_generators(n_states=14, seed=3):
    """independent lower bounds: odd Reynolds averages R_m = (1/432) sum_g sgn(g) m(g psi) of stabiliser-probability
    monomials (Pass 11434's construction) span exactly O_6, O_7, O_8 = 1, 2, 3 dimensions on random rays; g7 = R of
    p_a^4 p_b^2 p_c over three different MUBs is new at degree 7 for every one of the 24 oriented triples."""
    import itertools
    import w33_pass11434_h6_explicit as H
    Cl = P7.clifford_group(3)
    st, mub = H.stabiliser_states()
    rng = np.random.default_rng(seed)
    psis = [v / np.linalg.norm(v) for v in (rng.normal(size=(n_states, 3)) + 1j * rng.normal(size=(n_states, 3)))]
    A = np.array([[np.abs(st.conj() @ (C @ p)) ** 2 for C in Cl] for p in psis])
    B = np.array([[np.abs(st.conj() @ (C @ p.conj())) ** 2 for C in Cl] for p in psis])

    def R(m):
        m = list(m)
        return (np.prod(A[:, :, m], axis=2) - np.prod(B[:, :, m], axis=2)).sum(axis=1) / (2 * len(Cl))

    def mono(states, pat):
        return tuple(sum(([s] * e for s, e in zip(states, pat)), []))

    def cands(deg):
        pats = [p for p in itertools.product(range(deg + 1), repeat=4) if sum(p) == deg and sorted(p)[-1] <= deg - 1]
        out = {mono((0, 3, 6, 9), p) for p in pats}
        out |= {mono((0, 1, 3, 6), p) for p in pats if p[1] > 0}             # two states of the same MUB
        return sorted(out)

    out = {}
    for deg in (6, 7, 8):
        V = np.stack([R(m) for m in cands(deg)], 1)
        sv = np.linalg.svd(V, compute_uv=False)
        rank = int((sv > 1e-9 * sv[0]).sum())
        out[str(deg)] = dict(candidates=V.shape[1], rank=rank, last_kept_sv=float(sv[rank - 1]),
                             first_dropped_sv=float(sv[rank]) if rank < len(sv) else None)
    h6 = R(mono((0, 3, 6), (3, 2, 1)))
    res_g7 = []
    for tri in itertools.permutations(range(4), 3):
        v = R(mono(tuple(3 * t for t in tri), (4, 2, 1)))
        c = (v @ h6) / (h6 @ h6)
        res_g7.append(float(np.linalg.norm(v - c * h6) / np.linalg.norm(v)))
    out["g7_p4p2p_residual_vs_h6_all_24_oriented_triples"] = [min(res_g7), max(res_g7)]
    return out


def main():
    res = run()
    assert all(res["proved"].values())
    res["explicit_generators"] = explicit_generators()
    print(res["explicit_generators"], flush=True)
    assert [res["explicit_generators"][k]["rank"] for k in ("6", "7", "8")] == [res["sequences"]["odd"][k] for k in (6, 7, 8)]
    json.dump(res, open(OUT, "w"), indent=1)


if __name__ == "__main__":
    main()
