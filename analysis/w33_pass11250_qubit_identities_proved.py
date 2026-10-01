"""Pass 11250: the two identities behind the qubit protection limit are theorems in the literature -- and a q = 4 check.

Pass 11244 reduced the qubit limit lim E_n[c] = 0.66516065 to two identities about unipotent elements of Sp(2m,2),
verified (not proved) on every Jordan type with m <= 6; Pass 11245 reduced identity 2 to a uniformity statement.

IDENTITY 1 (Jordan-type weights = the odd-q formula at q = 2, summed over signs).  For Sp_{2n} in characteristic 2,
two unipotent classes lie in the same Lusztig 'unipotent piece' iff they have the same Jordan type (Xue, 'On unipotent
and nilpotent pieces for classical groups', arXiv:0912.3820, Lemma 2.7, citing Lusztig).  The number of F_q-points of a
piece is a polynomial in q with integer coefficients independent of the characteristic (Lusztig, 'Unipotent elements
in small characteristic', Transform. Groups 10 (2005), property P5, proved for Sp in section 3; restated for F_{p^s}
in Xue, section 6).  In odd characteristic the piece of type lambda is the single geometric class lambda, whose
F_q-points are the rational classes, i.e. the odd-q formula summed over the signs of the even parts.  Evaluating the
same polynomial at q = 2 gives identity 1.

IDENTITY 2 (the Hesselink index at 2).  In Lusztig's proof of P5-P8 for Sp, p = 2 (op. cit., section 3, the
identification of E with triples (w, alpha, j)), the fibre is a product in which each j_n ranges independently over
ALL nondegenerate symmetric bilinear forms on the multiplicity space of the size-n blocks, and the classes inside a
piece are distinguished exactly by which j_n are symplectic (alternating).  Pass 11244's chi_2 = 0 (omega(v,(u-1)v)
vanishes on ker(u-1)^2) says j_2 is alternating: on ker N^2 only the size-2 blocks contribute, since for a block of
size k >= 3 the pairing omega(N^{k-2}w, N^{k-1}w) involves N^{2k-3} w = 0.  So the chi_2 = 0 fraction equals
#{nondegenerate alternating} / #{nondegenerate symmetric} forms of rank m_2, which Pass 11245 computes:
2^{-m_2} for m_2 even and 0 for m_2 odd.

CONSEQUENCE.  The qubit limit 0.66516065 of Pass 11244 (bracket [0.6651604993, 0.6651614201] from the tail bound) is
no longer conditional on unproved identities; it rests on the cited theorems.

INDEPENDENT CHECK IN A NEW FIELD.  The literature statement is characteristic-free, so it must also hold at q = 4,
where q - 1 = 3 (Lusztig's P5 is stated for split structures with q - 1 sufficiently divisible; q = 2 was already
checked for m <= 6 by Pass 11244).  Here all 979200 elements of Sp(4,4) are enumerated, the unipotents are sorted by
Jordan type and chi_2, and both identities are tested in their general-q form.
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11215_protected_qutrit_limit as Q3  # noqa: E402

OUT = ROOT / "data" / "w33_pass11250_qubit_identities_proved.json"

# F_4 = {0, 1, a, a+1} coded 0, 1, 2, 3; addition = XOR
MUL = np.array([[0, 0, 0, 0], [0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2]], dtype=np.uint8)
J = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]], dtype=np.uint8)
I4 = np.eye(4, dtype=np.uint8)


def matmul(A, B):
    """batched F_4 matrix product: A (N,4,4), B (4,4) or (N,4,4)"""
    if B.ndim == 2:
        P = MUL[A[:, :, :, None], B[None, None, :, :]]
    else:
        P = MUL[A[:, :, :, None], B[:, None, :, :]]
    return np.bitwise_xor.reduce(P, axis=2)


def keys(A):
    flat = A.reshape(len(A), 16).astype(np.uint64)
    return (flat << (2 * np.arange(16, dtype=np.uint64))).sum(axis=1)


def transvection(v, lam):
    v = np.array(v, dtype=np.uint8)
    Jv = (MUL[J, v[None, :]]).astype(np.uint8)
    Jv = np.bitwise_xor.reduce(Jv, axis=1)                      # J v
    outer = MUL[MUL[lam, v][:, None], Jv[None, :]]              # lam v (Jv)^T
    return I4 ^ outer


def enumerate_sp44():
    gens = [transvection(v, lam) for v in ([1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [1, 1, 0, 0],
                                            [1, 0, 1, 0], [0, 1, 0, 1], [1, 0, 0, 1]) for lam in (1, 2)]
    for g in gens:
        assert (matmul(matmul(g.T[None].copy(), J), g)[0] == J).all()          # g^T J g = J
    seen = {int(k) for k in keys(I4[None])}
    frontier = I4[None]
    allm = [frontier]
    while len(frontier):
        new = []
        for g in gens:
            C = matmul(frontier, g)
            k = keys(C)
            _, first = np.unique(k, return_index=True)
            C, k = C[first], k[first]
            mask = np.array([int(x) not in seen for x in k])
            C, k = C[mask], k[mask]
            seen.update(int(x) for x in k)
            new.append(C)
        frontier = np.concatenate(new) if new else np.zeros((0, 4, 4), np.uint8)
        if len(frontier):
            allm.append(frontier)
    G = np.concatenate(allm)
    return G


def rank_f4(A):
    A = A.copy()
    r = 0
    rows, cols = A.shape
    inv = {1: 1, 2: 3, 3: 2}
    for c in range(cols):
        p = next((i for i in range(r, rows) if A[i, c]), None)
        if p is None:
            continue
        A[[r, p]] = A[[p, r]]
        A[r] = MUL[inv[int(A[r, c])], A[r]]
        for i in range(rows):
            if i != r and A[i, c]:
                A[i] ^= MUL[int(A[i, c]), A[r]]
        r += 1
    return r


VECS = np.array(np.meshgrid(*[range(4)] * 4, indexing="ij")).reshape(4, -1).T.astype(np.uint8)   # all 256 vectors


def apply(M, V):
    """M (4,4) applied to all vectors V (K,4)"""
    return np.bitwise_xor.reduce(MUL[M[None, :, :], V[:, None, :]], axis=2)


def jordan_and_chi2(N):
    k = [0] + [4 - rank_f4(_pow(N, j)) for j in range(1, 5)]
    lam = []
    for s in range(1, 5):
        lam += [s] * ((k[s] - k[s - 1]) - ((k[s + 1] - k[s]) if s + 1 <= 4 else 0))
    lam = tuple(sorted(lam, reverse=True))
    chi = None
    if 2 in lam:
        N2 = _pow(N, 2)
        ker = VECS[(apply(N2, VECS) == 0).all(axis=1)]
        Nv = apply(N, ker)
        JNv = apply(J, Nv)
        q = np.bitwise_xor.reduce(MUL[ker, JNv], axis=1)             # omega(v, N v) = v^T J N v
        chi = int((q != 0).any())
    return lam, chi


def _pow(N, j):
    P = N.copy()
    for _ in range(j - 1):
        P = matmul(P[None], N)[0]
    return P


def lam_weight(lam, QQ):
    """Pass 11244's weight (odd-q symplectic-partition formula summed over signs) at a general q"""
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


def gl_order(k, q):
    out = 1
    for i in range(k):
        out *= q ** k - q ** i
    return out


def gauss_binom(d, k, q):
    num = den = 1
    for i in range(k):
        num *= q ** (d - i) - 1
        den *= q ** (i + 1) - 1
    return num // den


def nondeg_symmetric(d, q):
    """Lusztig's recursion: f(d) = q^{d(d+1)/2} - sum_{d'>=1} [d, d']_q f(d - d')"""
    f = [1]
    for n in range(1, d + 1):
        f.append(q ** (n * (n + 1) // 2) - sum(gauss_binom(n, k, q) * f[n - k] for k in range(1, n + 1)))
    return f[d]


def p_chi0(m2, q):
    if m2 % 2:
        return Fraction(0)
    return Fraction(gl_order(m2, q) // Q3.sp_order(m2, q), nondeg_symmetric(m2, q))


def run():
    G = enumerate_sp44()
    order = 4 ** 4 * (4 ** 2 - 1) * (4 ** 4 - 1)
    res = dict(pass_id=11250, sp44_order_enumerated=len(G), sp44_order_formula=order)
    assert len(G) == order
    N = G ^ I4
    N4 = matmul(matmul(matmul(N, N), N), N)
    unip = N[(N4 == 0).reshape(len(N), 16).all(axis=1)]
    res["unipotents"] = len(unip)
    res["steinberg_q8"] = 4 ** 8
    by_lam, by_chi = defaultdict(int), defaultdict(int)
    for M in unip:
        lam, chi = jordan_and_chi2(M)
        by_lam[lam] += 1
        if chi is not None:
            by_chi[(lam, chi)] += 1
    id1 = {}
    for lam, cnt in by_lam.items():
        pred = lam_weight(lam, 4) * order
        id1[",".join(map(str, lam))] = dict(count=cnt, predicted=str(pred), equal=(pred == cnt))
    id2 = {}
    for lam in {k[0] for k in by_chi}:
        tot = by_chi[(lam, 0)] + by_chi[(lam, 1)]
        frac = Fraction(by_chi[(lam, 0)], tot)
        pred = p_chi0(Q3.mult(lam)[2], 4)
        id2[",".join(map(str, lam))] = dict(chi0=by_chi[(lam, 0)], total=tot, fraction=str(frac), predicted=str(pred),
                                            equal=(frac == pred))
    res["identity1_q4"] = id1
    res["identity2_q4"] = id2
    res["identity1_all"] = all(v["equal"] for v in id1.values())
    res["identity2_all"] = all(v["equal"] for v in id2.values())
    # the q = 2 counting lemma of Pass 11245, re-derived from the same general-q formulas
    res["alt_over_sym_q2"] = {k: str(p_chi0(k, 2)) for k in range(1, 9)}
    res["alt_over_sym_q4"] = {k: str(p_chi0(k, 4)) for k in range(1, 7)}
    res["qubit_limit"] = "0.66516065 (Pass 11244), now resting on Lusztig 2005 section 3 and Xue 2.7 instead of verification"
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
