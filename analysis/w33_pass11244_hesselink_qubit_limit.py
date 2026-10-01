#!/usr/bin/env python3
"""Pass 11244: the char-2 unipotent factor in closed form -- and the qubit protection limit to seven digits.

Pass 11233 pinned the qubit limit up to the unknown unipotent means F_m (m >= 7) of Sp(2m,2) and noted that the
char-2 (Hesselink) classes block Fulman's cycle index.  Here the exact class data of Sp(2m,2), m <= 6 (all unipotent
classes from the GAP files of Passes 11208/11220), are compared with two closed forms:

  IDENTITY 1 (Jordan types).  The total weight sum_{u unipotent, Jordan type lambda} 1/|Sp(2m,2)| equals the odd-q
     symplectic-partition weight of Pass 11215 evaluated at q = 2 (sum over the +- signs of the even parts):
        c(lambda, signs) = q^{dimC - dimR} prod_{i odd} |Sp(m_i,q)| prod_{i even} |O^{+-}(m_i,q)|.
     Checked exactly on all 92 Jordan types with m <= 6.
  IDENTITY 2 (Hesselink index at 2).  Among unipotents of Jordan type lambda, the fraction with chi_2 = 0 (the form
     omega(v, (u-1)v) vanishing on ker (u-1)^2) is 2^{-m_2} when m_2 = m_2(lambda) is even and 0 when it is odd,
     independently of the other parts.  Checked exactly on all 54 Jordan types with a part 2, m <= 6.
     (For m_2 even this is the chance that a random nondegenerate symmetric F_2-form of rank m_2 is alternating.)

Together they give F_m = E[m_1/2 + chi_2 m_2 | unipotent of Sp(2m,2)] for every m as a finite sum over partitions:
  F_m = sum_lambda w(lambda) (m_1/2 + (1 - [m_2 even] 2^{-m_2}) m_2) / sum_lambda w(lambda),
and the qubit limit = F(1)/P(1) + (unitary part of Pass 11233).  Status: the identities are verified, not proved, so
the limit below is CONDITIONAL on them; Pass 11233's bracket [0.6626, 0.6777] remains the unconditional statement.
An independent Monte Carlo of F_7 (random transvection walks, Pass 11233's sampler) tests the m = 7 prediction.
"""
from __future__ import annotations

import itertools
import json
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11208_arrow_universality as U  # noqa: E402
import w33_pass11215_protected_qutrit_limit as Q3  # noqa: E402
import w33_pass11217_qubit_arrow_law as Q  # noqa: E402
import w33_pass11220_qubit_mean_protection as M  # noqa: E402
import w33_pass11233_qubit_limit_cycle_index as P33  # noqa: E402

OUT = ROOT / "data" / "w33_pass11244_hesselink_qubit_limit.json"
QQ = 2


def lam_weight(lam):
    """Pass 11215's symplectic-partition weight at q = 2, summed over the signs of the even parts"""
    m = Q3.mult(lam)
    evens = [i for i in m if i % 2 == 0]
    expo2 = (Q3.conj_sq_sum(lam) + sum(mi for i, mi in m.items() if i % 2)
             - sum(mi * (mi + 1) for i, mi in m.items() if i % 2) - sum(mi * (mi - 1) for i, mi in m.items() if i % 2 == 0))
    odd = 1
    for i, mi in m.items():
        if i % 2:
            odd *= Q3.sp_order(mi, QQ)
    tot = Fraction(0)
    for s in range(1 << len(evens)):
        c = Fraction(odd)
        for b, i in enumerate(evens):
            c *= Q3.o_order(m[i], 1 if (s >> b) & 1 else -1, QQ)
        c *= Fraction(QQ) ** (expo2 // 2)
        tot += 1 / c
    return tot


def symplectic_partitions(size):
    for lam in Q3.partitions(size):
        m = Q3.mult(lam)
        if not any(i % 2 == 1 and mi % 2 == 1 for i, mi in m.items()):
            yield lam


def p_chi0(m2):
    return Fraction(1, 2 ** m2) if m2 % 2 == 0 else Fraction(0)


def f_value(lam):
    m = Q3.mult(lam)
    m1, m2 = m.get(1, 0), m.get(2, 0)
    return Fraction(m1, 2) + (1 - p_chi0(m2)) * m2


# ----------------------------------------------------------------------------- exact class data, m <= 6
def jordan(N, D):
    k = [0] + [Q.dimker(Q.mpow(N, j)) for j in range(1, D + 1)]
    lam = []
    for s in range(1, D + 1):
        lam += [s] * ((k[s] - k[s - 1]) - ((k[s + 1] - k[s]) if s + 1 <= D else 0))
    return tuple(sorted(lam, reverse=True))


def chi2(N, J):
    K = Q.nullspace(Q.mpow(N, 2))
    for c in itertools.product(range(2), repeat=len(K)):
        if any(c):
            v = (np.array(c) @ K) % 2
            if (v @ J @ N @ v) % 2:
                return 1
    return 0


def class_data(maxm=6):
    by_lam, by_chi = defaultdict(Fraction), defaultdict(Fraction)
    for m in range(1, maxm + 1):
        if m == 1:
            els = [np.array([[a, b], [c, d]]) for a in range(2) for b in range(2) for c in range(2) for d in range(2)
                   if (a * d - b * c) % 2 == 1]
            cls, order = [dict(S=e, size=1) for e in els], 6
        else:
            cls, head = U.load_classes(M.class_file(m), 2, m)
            order = head[0]
        J = U.form(m, 2)
        I = np.eye(2 * m, dtype=np.int64)
        for c in cls:
            N = (c["S"] + I) % 2
            if Q.dimker(Q.mpow(N, 2 * m)) != 2 * m:
                continue
            lam = jordan(N, 2 * m)
            w = Fraction(c["size"], order)
            by_lam[lam] += w
            if 2 in lam:
                by_chi[(lam, chi2(N, J))] += w
    return by_lam, by_chi


def verify(maxm=6):
    by_lam, by_chi = class_data(maxm)
    id1 = {lam: by_lam[lam] == lam_weight(lam) for lam in by_lam}
    allsp = {lam for size in range(2, 2 * maxm + 1, 2) for lam in symplectic_partitions(size)}
    id2 = {}
    for lam in {k[0] for k in by_chi}:
        tot = by_chi[(lam, 0)] + by_chi[(lam, 1)]
        id2[lam] = by_chi[(lam, 0)] / tot == p_chi0(Q3.mult(lam)[2])
    return dict(identity1_checked=len(id1), identity1_all=all(id1.values()),
                identity1_covers_all_symplectic_partitions=set(by_lam) == allsp,
                identity2_checked=len(id2), identity2_all=all(id2.values()))


def F_series(N):
    W, Fw = [Fraction(0)] * (N + 1), [Fraction(0)] * (N + 1)
    for k in range(N + 1):
        for lam in symplectic_partitions(2 * k):
            w = lam_weight(lam)
            W[k] += w
            Fw[k] += w * f_value(lam)
    return W, Fw


def run(N=20, mc=True):
    from mpmath import mp, mpf
    mp.dps = 40
    res = dict(pass_id=11244)
    res.update(verify())
    W, Fw = F_series(N)
    res["steinberg_all_m"] = all(W[m] == P33.w(m) for m in range(N + 1))
    exact = json.loads(P33.OUT.read_text())["F_m_exact"]
    res["matches_pass11233_F_m"] = all(Fw[m] / W[m] == Fraction(exact[str(m)]["value"]) for m in range(1, 7))
    res["F_m"] = {str(m): dict(value=str(Fw[m] / W[m]), float=float(Fw[m] / W[m])) for m in range(1, N + 1)}
    f = lambda x: mpf(x.numerator) / x.denominator
    P1 = mpf(1)
    for i in range(400):
        P1 /= 1 - mpf(2) ** (-1 - 2 * i)
    known = sum(f(Fw[m]) for m in range(N + 1))
    tailw = P1 - sum(f(W[m]) for m in range(N + 1))
    unitary = mpf(json.loads(P33.OUT.read_text())["unitary_limit"])
    Finf = f(Fw[N] / W[N])
    res["unipotent_limit_part"] = str(known / P1)
    res["tail_weight"] = str(tailw / P1)
    res["limit_conditional"] = str(known / P1 + Finf * tailw / P1 + unitary)
    res["limit_conditional_bracket"] = [str(known / P1 + unitary), str((known + mp.log(3, 2) * tailw) / P1 + unitary)]
    res["pass11233_flat_estimate"] = json.loads(P33.OUT.read_text())["limit_if_F_m_flat"]
    if mc:
        res["monte_carlo_F7"] = P33.mc_F(7, 20000, 12, 11244)
        r = res["monte_carlo_F7"]
        res["mc_F7_z"] = (r["F_mean"] - float(Fw[7] / W[7])) / r["F_stderr"]
    return res


def main():
    res = run(mc="--no-mc" not in sys.argv)
    OUT.write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps({k: v for k, v in res.items() if k != "F_m"}, indent=1, default=str))
    print({m: res["F_m"][m]["float"] for m in list(res["F_m"])[:10]})


if __name__ == "__main__":
    main()
