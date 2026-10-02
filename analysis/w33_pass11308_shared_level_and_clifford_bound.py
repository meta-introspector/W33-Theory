"""Pass 11308: the shared cyclotomic level, derived -- and how close-to-Clifford bounds time-reversal violation.

Codex's balanced-G26 note (PASS20261001_BALANCED_G26_TMAGIC_VACUUM.md) records that
((1 + 2cos 2pi/9)/3)^2 = 0.712386... is both the largest qutrit-stabilizer probability of the T-magic vacuum ray
|T3> = T|+> and the exact two-qutrit reversal-fidelity level of Passes 11238/11251, 'a shared algebraic invariant, not
yet a derived physical identity'.

DERIVATION (exact).  T = diag(zeta^{x^3}) has spectrum zeta^{0, 1, -1}.  (i) <+|T|+> = tr T / 3 = (1 + zeta + zeta^-1)/3
= F_min, so |<+|T3>|^2 = |tr T/3|^2 = F_min^2, and |+> attains the maximum over the 12 stabilizer states.  (ii) Pass
11251: the optimal residual of (I(x)T)SUM(I(x)T^2) has the spectrum of T (x) T (the sumset {-1,0,1}+{-1,0,1}), so its
fidelity is |tr(T (x) T)|/9 = |tr T/3|^2.  Both numbers are the squared normalised trace of the cubic gate.  (iii) The
same number is the Clifford fidelity F_Cl(U) = max_C |tr(C^dag U)|/d of that word: the three are one invariant.

BOUND (proof).  Let C be a Clifford with |tr(C^dag U)| = d F_Cl, U = C R.  Every Clifford has an anti-unitary Clifford
reversal V C^* V^dag = lambda C^dag (no Clifford T-violation; Wonenburger, cited in the repo).  Then
tr(V U^* V^dag U) = lambda tr(A R) with A = C^dag V R^* V^dag C, |tr A| = |tr R|.  Split A, R into trace and traceless
parts; Cauchy-Schwarz on the traceless parts (Frobenius norms^2 = d(1 - F_Cl^2)) gives
        F_T(U) >= 2 F_Cl(U)^2 - 1.
Tested numerically: the sharper F_T >= F_Cl^2 on every level, and F_T >= F_Cl on two qutrits (false on one qutrit).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11213_cubic_t_violation as P1  # noqa: E402
import w33_pass11227_two_qutrit_cp_mixing as P2  # noqa: E402
import w33_pass11228_t_witness as W1  # noqa: E402
import w33_pass11251_exact_reversal as E  # noqa: E402
import w33_pass11253_cubic_field_levels as L53  # noqa: E402

OUT = ROOT / "data" / "w33_pass11308_shared_level_and_clifford_bound.json"
FMIN = (1 + 2 * np.cos(2 * np.pi / 9)) / 3


def exact_identities():
    # tr T / 3 and its square, in Z[zeta]
    g = E.add(E.add(E.ONE, E.zpow(1)), E.zpow(8))                      # tr T = 1 + z + z^8
    g2 = E.mul(g, g)
    sumset = E.ZERO
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            sumset = E.add(sumset, E.zpow(i + j))
    w = np.exp(2j * np.pi / 3)
    stab = [np.eye(3)[i] for i in range(3)] + [np.array([1, w ** a, w ** b]) / np.sqrt(3) for a in range(3) for b in range(3)]
    plus = np.ones(3) / np.sqrt(3)
    T3 = P1.T1 @ plus
    probs = [abs(np.vdot(s, T3)) ** 2 for s in stab]
    return dict(trT_over3_equals_Fmin=bool(abs(E.to_complex(g) / 3 - FMIN) < 1e-14),
                T_tensor_T_trace_equals_sumset=(sumset == g2),
                T3_max_stabilizer_probability=max(probs),
                T3_max_attained_at_plus=bool(np.argmax(probs) == int(np.argmax([abs(np.vdot(s, plus)) for s in stab]))),
                equals_Fmin_squared=bool(abs(max(probs) - FMIN ** 2) < 1e-14))


def fcl_one(U, Cf):
    return max(abs(np.trace(C.conj().T @ U)) / 3 for C in Cf)


def fcl_two(U, reps):
    best = 0.0
    for s in range(0, len(reps), 4000):
        B = np.einsum('rji,jk->rik', reps[s:s + 4000].conj(), U)
        best = max(best, float(np.abs(np.einsum('aji,rji->ra', P2.PA.conj(), B)).max() / 9))
    return best


def run():
    res = dict(pass_id=11308, identities=exact_identities())
    print(res["identities"], flush=True)
    Cf = P1.clifford1()
    rows1 = []
    for d, fac in L53.one_qutrit_stream():
        U = L53.prod(fac)
        f = W1.fidelity(U, Cf)
        if f > 1 - 1e-9:
            continue
        rows1.append((d, f, fcl_one(U, Cf)))
    reps = np.load(P2.CACHE)
    rows2 = []
    for d, fac in L53.two_qutrit_stream():
        U = L53.prod(fac)
        f = float(E.float_maximisers(U, reps).max())
        if f > 1 - 1e-9:
            continue
        rows2.append((d, f, fcl_two(U, reps)))
    out = {}
    for name, rows in (("one_qutrit", rows1), ("two_qutrit", rows2)):
        a = np.array([[f, s] for _, f, s in rows])
        out[name] = dict(
            violating_words=len(rows),
            min_FT_minus_2Fcl2_plus1=float((a[:, 0] - (2 * a[:, 1] ** 2 - 1)).min()),
            min_FT_minus_Fcl2=float((a[:, 0] - a[:, 1] ** 2).min()),
            min_FT_minus_Fcl=float((a[:, 0] - a[:, 1]).min()),
            equality_FT_eq_Fcl=int(np.sum(np.abs(a[:, 0] - a[:, 1]) < 1e-9)),
            Fcl_values=sorted({round(float(s), 9) for s in a[:, 1]}),
            counterexamples_FT_lt_Fcl=sorted({(round(float(f), 7), round(float(s), 7)) for f, s in a if f < s - 1e-9}),
        )
        print(name, out[name], flush=True)
    res["bounds"] = out
    res["proved"] = "F_T >= 2 F_Cl^2 - 1 (any n)"
    res["holds_FT_ge_2Fcl2_minus_1"] = all(v["min_FT_minus_2Fcl2_plus1"] > -1e-9 for v in out.values())
    res["holds_FT_ge_Fcl2_on_sample"] = all(v["min_FT_minus_Fcl2"] > -1e-9 for v in out.values())
    res["holds_FT_ge_Fcl_two_qutrit_sample"] = out["two_qutrit"]["min_FT_minus_Fcl"] > -1e-9
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in res.items() if k != "bounds"}, indent=1, default=str))


if __name__ == "__main__":
    main()
