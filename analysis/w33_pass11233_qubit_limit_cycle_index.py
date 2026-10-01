#!/usr/bin/env python3
"""Pass 11233: the qubit protection limit lim E_n[c] -- exact unitary part, exact unipotent data to n = 6, rigorous bracket.

Pass 11220 computed E_n[c] for random Clifford ticks on n qubits exactly for n <= 6 from GAP's classes and
extrapolated lim E_n[c] ~ 0.665166 (geometric fit; not a proof).  Master's Pass 11215 obtained the qutrit limit from
Fulman's cycle index; for qubits the unipotent factor needs Hesselink's char-2 classes, which is why 11220 summed classes.

Here the char-2 cycle index is assembled from what is exact:
  * c = m_1(1)/2 + chi_2 m_2(1) + m_1(x^2+x+1)  (Pass 11217; additive over the primary parts);
  * the x^2+x+1 part is the only self-reciprocal irreducible of degree 2 over F_2: its factor is the unitary sum
    sum_lambda u^{|lambda|} / c_U(lambda) with Q = 2 (Pass 11215's unitary_types with q = 2);
  * the unipotent (x+1) factor is P(u) = sum_m u^m w_m, w_m = 2^{2m^2}/|Sp(2m,2)| (Steinberg), which is the q-binomial
    product P(u) = prod_{i>=0} 1/(1 - u 2^{-1-2i}) -- so 1/P is entire and the unipotent contribution to the limit is
    exactly F(1)/P(1), F(u) = sum_m u^m w_m F_m, F_m = E[m_1/2 + chi_2 m_2 | u unipotent in Sp(2m,2)];
  * F_m for m <= 6 is exact from the GAP class data (Pass 11208/11220 files);
  * F_m <= log2(3 - 2^{1-2m}) for every m (rigorous: f <= dim ker(u-1) = d, Jensen, and E[2^d] over unipotents is
    1 + (2^{2m}-1) 2^{1-2m} = 3 - 2^{1-2m}, counting unipotents in the vector stabiliser 2^{2m-1}:Sp(2m-2,2)).
Validation: Steinberg's counts, E[2^d] on the class data, and E_n[c] for n = 2..6 reproduced exactly from the factors.
Monte Carlo F_7, F_8 (random transvection walks, unipotent rejection; validated on m = 5, 6 against the exact F_m) give
a statistical estimate inside the rigorous bracket.
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11208_arrow_universality as U  # noqa: E402
import w33_pass11215_protected_qutrit_limit as Q3  # noqa: E402
import w33_pass11217_qubit_arrow_law as Q  # noqa: E402
import w33_pass11220_qubit_mean_protection as M  # noqa: E402

OUT = ROOT / "data" / "w33_pass11233_qubit_limit_cycle_index.json"


def sp2_order(m):
    return Q3.sp_order(2 * m, 2)


def w(m):
    return Fraction(2 ** (2 * m * m), sp2_order(m))


def unitary_types_q2(N):
    out = []
    for size in range(0, N + 1):
        for lam in Q3.partitions(size):
            mm = Q3.mult(lam)
            c = Fraction(2) ** (Q3.conj_sq_sum(lam) - sum(x * x for x in mm.values()))
            for x in mm.values():
                c *= Q3.u_order(x, 2)
            out.append((size, 1 / c, mm.get(1, 0)))
    return out


def unipotent_data(m):
    """exact (sum over unipotent classes of 1/|C|, of f/|C|, of 2^d/|C|) from the GAP classes of Sp(2m,2)"""
    if m == 1:      # Sp(2,2) = SL(2,2): every element is its own class representative, size 1
        els = [np.array([[a, b], [c, d]]) for a in range(2) for b in range(2) for c in range(2) for d in range(2)
               if (a * d - b * c) % 2 == 1]
        cls, head = [dict(S=e, size=1) for e in els], (6, 6)
    else:
        cls, head = U.load_classes(M.class_file(m), 2, m)
    order = head[0]
    I = np.eye(2 * m, dtype=np.int64)
    tot = ftot = ptot = Fraction(0)
    for c in cls:
        N = (c["S"] + I) % 2
        if Q.dimker(Q.mpow(N, 2 * m)) != 2 * m:
            continue
        wt = Fraction(c["size"], order)
        tot += wt
        ftot += wt * Q.c_formula(c["S"] % 2, m)
        ptot += wt * 2 ** Q.dimker(N)
    return tot, ftot, ptot


def series_inverse(a, N):
    return Q3.inverse(a, N)


def exact_part(N=6, NU=40):
    res = {}
    W = [w(m) for m in range(N + 1)]
    Fm, checks = [Fraction(0)], []
    for m in range(1, N + 1):
        tot, ftot, ptot = unipotent_data(m)
        checks.append(dict(m=m, steinberg=tot == W[m], E_2_to_d=str(ptot / tot),
                           E_2_to_d_formula_ok=ptot / tot == 3 - Fraction(2, 4 ** m)))
        Fm.append(ftot / tot)
    res["unipotent_checks"] = checks
    res["F_m_exact"] = {str(m): dict(value=str(Fm[m]), float=float(Fm[m])) for m in range(1, N + 1)}
    # finite-n reconstruction of E_n[c], n <= N
    Fu = [W[m] * Fm[m] for m in range(N + 1)]
    a = Q3.mul(Fu, series_inverse(W, N), N)
    uni = unitary_types_q2(NU)
    PU = Q3.series(uni, NU, False)
    FU = Q3.series(uni, NU, True)
    b = Q3.mul(FU, series_inverse(PU, NU), NU)
    ref = {r["n"]: Fraction(r["E_c"]) for r in json.loads(M.OUT.read_text())["means"]}
    recon = {n: sum(a[: n + 1]) + sum(b[: n + 1]) for n in range(2, N + 1)}
    res["reconstructs_pass11220"] = {str(n): recon[n] == ref[n] for n in ref if n <= N}
    res["unitary_steinberg"] = all(sum(x for k, x, _ in uni if k == n) * Q3.u_order(n, 2) == 2 ** (n * (n - 1))
                                   for n in range(0, 8))
    # the unitary (x^2+x+1) contribution to the limit: partial sums of b converge geometrically
    from mpmath import mp, mpf
    mp.dps = 40
    f = lambda x: mpf(x.numerator) / x.denominator
    Bpart = [f(sum(b[: n + 1])) for n in range(NU + 1)]
    res["unitary_limit"] = str(Bpart[NU])
    res["unitary_limit_change_last5"] = float(abs(Bpart[NU] - Bpart[NU - 5]))
    # unipotent contribution: F(1)/P(1); P(1) = prod 1/(1 - 2^{-1-2i})
    P1 = mp.mpf(1)
    for i in range(200):
        P1 /= (1 - mpf(2) ** (-1 - 2 * i))
    known = sum(f(Fu[m]) for m in range(N + 1))
    tail_w = P1 - sum(f(W[m]) for m in range(N + 1))
    tail_max = sum(f(w(m)) * mp.log(3 - mpf(2) ** (1 - 2 * m), 2) for m in range(N + 1, 120))
    lo, hi = known / P1, (known + tail_max) / P1
    res["P1"] = str(P1)
    res["unipotent_known_part"] = str(lo)
    res["unipotent_tail_weight"] = str(tail_w / P1)
    res["limit_bracket"] = [str(lo + Bpart[NU]), str(hi + Bpart[NU])]
    # point estimate: F_m = F_N for m > N (F_4, F_5, F_6 differ by 4.4e-4, 4.6e-5: geometric ratio ~ 1/10)
    NG = 40
    Wg = [w(m) for m in range(NG + 1)]
    Fg = [Wg[m] * (Fm[m] if m <= N else Fm[N]) for m in range(NG + 1)]
    ag = Q3.mul(Fg, series_inverse(Wg, NG), NG)
    En = [f(sum(ag[: n + 1])) + Bpart[min(n, NU)] for n in range(NG + 1)]
    res["E_n_if_F_m_flat"] = {str(n): str(En[n]) for n in range(2, 13)}
    res["limit_if_F_m_flat"] = str(En[NG])
    res["sensitivity_per_unit_F_tail"] = str(tail_w / P1)
    res["E_n_not_monotone"] = bool(En[7] < En[6])
    res["_floats"] = dict(Fm=[float(x) for x in Fm], W=[float(x) for x in W], P1=float(P1), known=float(known),
                          unitary=float(Bpart[NU]))
    return res


# ------------------------------------------------------------------ Monte Carlo F_m
def symplectic_form(m):
    return U.form(m, 2)


def random_sp(m, batch, rng, steps=None):
    """batch of (approximately uniform) elements of Sp(2m,2): products of random transvections x -> x + w(x,v) v"""
    D = 2 * m
    J = symplectic_form(m).astype(np.float32)
    S = np.broadcast_to(np.eye(D, dtype=np.float32), (batch, D, D)).copy()
    for _ in range(steps or (12 * D + 40)):
        v = rng.integers(0, 2, size=(batch, D)).astype(np.float32)
        r = np.einsum("bi,ij,bjk->bk", v, J, S) % 2          # row (v^T J S)
        S = (S + v[:, :, None] * r[:, None, :]) % 2
    return S


def mc_F(m, batch, rounds, seed):
    rng = np.random.default_rng(seed)
    D = 2 * m
    J = symplectic_form(m)
    vals, total, pow2 = [], 0, []
    for _ in range(rounds):
        S = random_sp(m, batch, rng)
        N = (S + np.eye(D, dtype=np.float32)) % 2
        P = N.copy()
        k = 1
        while k < D:
            P = np.matmul(P, P) % 2
            k *= 2
        unip = np.flatnonzero(~P.reshape(batch, -1).any(axis=1))
        total += batch
        for i in unip:
            Si = S[i].astype(np.int64)
            assert not ((Si.T @ J @ Si - J) % 2).any()
            vals.append(Q.c_formula(Si, m))
            pow2.append(2 ** Q.dimker((Si + np.eye(D, dtype=np.int64)) % 2))
    vals = np.array(vals, float)
    return dict(m=m, samples=total, unipotent=len(vals), p_unipotent=len(vals) / total, p_unipotent_exact=float(w(m)),
                F_mean=float(vals.mean()), F_stderr=float(vals.std(ddof=1) / np.sqrt(len(vals))),
                E_2_to_d=float(np.mean(pow2)), E_2_to_d_exact=float(3 - Fraction(2, 4 ** m)))


def run(mc=True):
    res = dict(pass_id=11233)
    res.update(exact_part())
    if mc:
        rounds = {5: 6, 6: 10, 7: 24, 8: 40}
        res["monte_carlo"] = [mc_F(m, 20000, rounds[m], 11233 + m) for m in (5, 6, 7, 8)]
        from mpmath import mpf
        fl = res["_floats"]
        est = fl["known"]
        var = 0.0
        for r in res["monte_carlo"]:
            if r["m"] >= 7:
                est += fl_w(r["m"]) * r["F_mean"]
                var += (fl_w(r["m"]) * r["F_stderr"]) ** 2
        # m >= 9: use the m = 8 estimate as a flat guess (weight 2^-9-ish), error = its full bracket
        tail9 = sum(fl_w(m) for m in range(9, 120))
        F8 = res["monte_carlo"][-1]["F_mean"]
        est += tail9 * F8
        res["mc_estimate"] = dict(value=(est / fl["P1"]) + fl["unitary"],
                                  stat_err=float(np.sqrt(var)) / fl["P1"],
                                  m9_tail_weight=tail9 / fl["P1"],
                                  pass11220_extrapolation=json.loads(M.OUT.read_text())["extrapolated_limit"])
    return res


def fl_w(m):
    return float(w(m))


def main():
    res = run(mc="--no-mc" not in sys.argv)
    OUT.write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
