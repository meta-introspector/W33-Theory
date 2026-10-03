"""Pass 11355: a substrate Jarlskog invariant -- a degree-six closed-form witness of time-reversal violation.

THEOREM 1 (exact).  A tick U is substrate-time-reversible (an anti-unitary Clifford V K with V Ubar V^dag = lambda U^-1
exists) iff U^T is Clifford-conjugate to U up to a phase.  Proof: V Ubar V^dag = lambda U^-1 <=> Ubar = lambda V^dag U^-1 V
<=> U^T = (Ubar)^-1 = lambda-bar V^dag U V; conversely read backwards.  (Pass 11252's criterion, restated.)

COROLLARY.  Every Clifford-conjugation-invariant, phase-invariant function f satisfies f(U^T) = f(U) on reversible ticks,
so f(U) != f(U^T) witnesses T-violation.  The invariant polynomials of degree (t, t) in (U, Ubar) are tr(X pi_t(U)),
pi_t(U) = U^(x)t (x) Ubar^(x)t, X in the commutant of the Clifford group; all are captured by
A_t(U) = avg_C pi_t(C U C^dag) (the orthogonal projection of pi_t(U) onto the commutant), and since
<pi_t(X), pi_t(Y)> = |tr(X^dag Y)|^(2t):
        J_2t(U) := avg_C |tr(U^dag C U C^dag)|^(2t) - avg_C |tr(U^dag C U^T C^dag)|^(2t)  =  1/2 ||A_t(U) - A_t(U^T)||^2 >= 0,
with J_2t(U) > 0 iff some degree-(t,t) invariant separates U from U^T.
THEOREM 2.  J_2 = 0 identically for any unitary 2-design (degree-(1,1) invariants: the commutant of C (x) Cbar is
spanned by the identity and the maximally entangled projector, both transpose-symmetric).  J_4 also vanishes
identically for qutrits -- NOT because of a design property (the degree-(t,t) invariants need the frame potential at 2t;
an earlier draft wrongly derived J_4 = 0 from the 2-design property), but as a verified polynomial identity (Pass 11357:
Haar identity test; J_4 is NOT identically zero for d = 5).  So for qutrits the lowest witness degree is (3,3) -- cubic,
like the Jarlskog invariant Im tr[H_u, H_d]^3 of the quark sector.
FINDING (computer-verified, one qutrit).  J_6 is COMPLETE on every word checked: J_6(U) > 0 exactly for the violators --
per word against the exact Weyl criterion for all words with 1, 2, 3 cubic gates (coset-reduced: 216, 5184, 124,416), and
by exact count for 4 gates (2,985,984 words: the number with J_6 > 0 equals the 2,077,650 violators of Pass 11312).
The identity J_6 = 1/2||A_3(U) - A_3(U^T)||^2 is checked on explicit 729 x 729 matrices.
"""

from __future__ import annotations

import itertools
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11213_cubic_t_violation as P1  # noqa: E402
import w33_pass11252_exact_reversibility as R  # noqa: E402
import w33_pass11312_depth_law as DL  # noqa: E402

OUT = ROOT / "data" / "w33_pass11355_substrate_jarlskog.json"
CF = np.array(P1.clifford1())
TOL = 1e-9


def J(U, t=3, Cs=CF):
    X = np.einsum('cij,jk,clk->cil', Cs, U, Cs.conj())
    Y = np.einsum('cij,jk,clk->cil', Cs, U.T, Cs.conj())
    a = np.abs(np.einsum('ji,cji->c', U.conj(), X)) ** (2 * t)
    b = np.abs(np.einsum('ji,cji->c', U.conj(), Y)) ** (2 * t)
    return float(a.mean() - b.mean())


def det1(U):
    return U / np.linalg.det(U) ** (1 / 3)


def A3(U):
    S = np.zeros((729, 729), complex)
    for C in CF:
        V = det1(C @ U @ C.conj().T)
        V3 = np.kron(np.kron(V, V), V)
        S += np.kron(V3, V3.conj())
    return S / 216


_S = {}


def _init():
    R.WEYL[1] = R.Weyl(1)
    _S["R"] = DL.coset_reps(list(CF))


def _block(args):
    """words R_k T ... R_2 T C_1 T (Pass 11312 reduction); per-word J6 and (if check) the exact verdict"""
    k, prefix, check = args
    Rr, T = _S["R"], P1.T1
    rng = np.random.default_rng(0)
    pos = neg = mism = 0
    jmin = np.inf
    for rest in itertools.product(range(24), repeat=k - 1 - len(prefix)):
        W = np.eye(3, dtype=complex)
        for r in list(prefix) + list(rest):
            W = W @ Rr[r] @ T
        for C in CF:
            U = W @ C @ T
            j = J(U)
            pos += j > TOL
            if j > TOL:
                jmin = min(jmin, j)
            if check:
                viol = R.decide(U, 1, rng)[0] is False
                mism += viol != (j > TOL)
    return pos, mism, jmin


def run():
    res = dict(pass_id=11355)
    rng = np.random.default_rng(11355)
    # identity check on explicit matrices (violators and reversible ticks)
    R.WEYL[1] = R.Weyl(1)
    dev = 0.0
    for _ in range(6):
        U = det1(CF[rng.integers(216)] @ P1.T1 @ CF[rng.integers(216)] @ P1.T1)
        dev = max(dev, abs(0.5 * np.linalg.norm(A3(U) - A3(U.T)) ** 2 - J(U)))
    res["identity_J6_equals_half_norm_sq_max_dev"] = dev
    res["J2_J4_max_abs_on_random_words"] = max(abs(J(U, t)) for U in
                                               [CF[rng.integers(216)] @ P1.T1 @ CF[rng.integers(216)] @ P1.T1
                                                for _ in range(200)] for t in (1, 2))
    exact = {1: 18, 2: 2106, 3: 68202, 4: 2077650}
    with Pool(11, initializer=_init) as pool:
        out = {}
        for k in (1, 2, 3, 4):
            check = k <= 3
            if k == 1:
                jobs = [(1, (), check)]
            elif k <= 3:
                jobs = [(k, (i,), check) for i in range(24)]
            else:
                jobs = [(k, (i, j), False) for i in range(24) for j in range(24)]
            r = pool.map(_block, jobs)
            pos = sum(x[0] for x in r)
            mism = sum(x[1] for x in r)
            jmin = min(x[2] for x in r)
            out[str(k)] = dict(words=216 * 24 ** (k - 1), J6_positive=int(pos), exact_violators=exact[k],
                               counts_equal=pos == exact[k], per_word_checked=check, per_word_mismatches=int(mism),
                               min_positive_J6=float(jmin))
            print(k, out[str(k)], flush=True)
    res["one_qutrit"] = out
    res["J6_complete_on_all_checked"] = all(v["counts_equal"] and v["per_word_mismatches"] == 0 for v in out.values())
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in res.items() if k != "one_qutrit"}, indent=1, default=str))


if __name__ == "__main__":
    main()
