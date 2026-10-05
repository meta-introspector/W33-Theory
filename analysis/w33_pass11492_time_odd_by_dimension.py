"""Pass 11492: the lowest time-odd Clifford invariant of a qudit STATE, for d = 2, 3, 5 -- the state-level companion of
Pass 11357's Jarlskog degrees (10, 6, 4 for unitaries).

For the projective Clifford group G_d of one qudit (orders 24, 216, 3000) and k = 0..KMAX:
    D_k  = (1/|G|) sum_g |h_k(eig g)|^2,   Tw_k = (1/|G|) sum_g h_k(eig(conj(g) g)),   O_k = (D_k - Tw_k)/2
(Pass 11491).  Exact: eigenvalues are identified as N-th roots of unity, h_k(zeta^e) is an integer vector over Z/N
(complete homogeneous counts by exponent), and the cyclotomic sums are evaluated at 60 digits and asserted to be integers
divisible by |G|.  Conjugacy is used only to group identical eigenvalue multisets.
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from math import lcm
from pathlib import Path

import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11357_jarlskog_degree_by_dimension as P7  # noqa: E402

OUT = ROOT / "data" / "w33_pass11492_time_odd_by_dimension.json"
mp.mp.dps = 60


def find_N(angles_list, cands):
    for N in cands:
        ok = all(np.abs(a * N / (2 * np.pi) - np.rint(a * N / (2 * np.pi))).max() < 1e-7 for a in angles_list)
        if ok:
            return N
    raise AssertionError("no N")


def h_vectors(e, N, K):
    """H[k] = integer vector over Z/N with h_k(zeta^e) = sum_r H[k][r] zeta^r, via h^(i)_k = h^(i-1)_k + shift(h^(i)_(k-1))"""
    prev = np.zeros((K + 1, N), dtype=object)
    prev[0, 0] = 1                                          # zero variables: h_0 = 1, h_k = 0
    for ei in e:
        cur = np.zeros((K + 1, N), dtype=object)
        cur[0] = prev[0]
        for k in range(1, K + 1):
            cur[k] = prev[k] + np.roll(cur[k - 1], int(ei))
        prev = cur
    return prev


def zsum(vec, N):
    z = [mp.exp(2j * mp.pi * r / N) for r in range(N)]
    return mp.fsum(int(c) * z[r] for r, c in enumerate(vec) if c)


def clifford_proj(p, seed=1):
    """projective Clifford group of one qudit (p odd prime) as {W(a) V_S}, V_S from the Weyl twirl
    sum_q W(Sq) A W(q)^dag (an intertwiner, unique up to scale by Schur); symmetric Weyl W(a,b) = omega^(ab/2) X^a Z^b.
    Independent of Pass 11357's BFS construction (which it reproduces at p = 3, 5)."""
    import itertools
    w = np.exp(2j * np.pi / p)
    X = np.roll(np.eye(p), 1, axis=0)
    Zm = np.diag([w ** k for k in range(p)])
    h = pow(2, -1, p)
    Ws = {(a, b): w ** ((h * a * b) % p) * np.linalg.matrix_power(X, a) @ np.linalg.matrix_power(Zm, b)
          for a in range(p) for b in range(p)}
    A = np.random.default_rng(seed).normal(size=(p, p)) + 1j * np.random.default_rng(seed + 1).normal(size=(p, p))
    out = []
    for a, b, c, dd in itertools.product(range(p), repeat=4):
        if (a * dd - b * c) % p != 1:
            continue
        S = np.array([[a, b], [c, dd]])
        V = sum(Ws[tuple(int(x) for x in (S @ np.array(q)) % p)] @ A @ Ws[q].conj().T for q in Ws)
        V = V / np.sqrt(np.trace(V.conj().T @ V).real / p)
        assert np.allclose(V.conj().T @ V, np.eye(p), atol=1e-8)
        out.append(V)
    return np.array([Ws[q] @ V for q in Ws for V in out])


def series(d, K, G=None):
    if G is None:
        G = P7.clifford_group(d) if d in (2, 3) else clifford_proj(d)
    lam = [np.linalg.eigvals(g) for g in G]
    rat = [np.angle(x / x[0]) for x in lam]
    tw = [np.angle(np.linalg.eigvals(np.conj(g) @ g)) for g in G]
    cands = sorted({a * b for a in (1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 24, 30, 40, 60, 120) for b in (1, d, d * d)})
    N = find_N(rat + tw, cands)
    cr = Counter(tuple(sorted(int(x) for x in np.rint(a * N / (2 * np.pi)) % N)) for a in rat)
    ct = Counter(tuple(sorted(int(x) for x in np.rint(a * N / (2 * np.pi)) % N)) for a in tw)
    return series_from_classes(cr, ct, N, len(G), K, d)


def series_from_classes(cr, ct, N, order, K, d):
    """the exact series from eigenvalue-multiset class counts: cr (ratio exponents of g), ct (exponents of conj(g) g)"""
    Dsum = np.zeros((K + 1, N), dtype=object)
    for e, m in cr.items():
        H = h_vectors(e, N, K)
        for k in range(K + 1):
            v = np.array(H[k], dtype=np.int64)
            assert v.max() < 2 ** 26
            full = np.convolve(v, v[::-1])                      # index (r - s) + N - 1
            acc = full[N - 1:].copy()
            acc[1:] += full[:N - 1]                             # fold r - s < 0 into Z/N
            Dsum[k] += m * acc.astype(object)
    Tsum = np.zeros((K + 1, N), dtype=object)
    for e, m in ct.items():
        Tsum += m * h_vectors(e, N, K)
    D, Tw = [], []
    for k in range(K + 1):
        for vec, out in ((Dsum[k], D), (Tsum[k], Tw)):
            x = zsum(vec, N)
            xi = int(mp.nint(mp.re(x)))
            assert abs(x - xi) < mp.mpf(10) ** -25 and xi % order == 0, (d, k, x)
            out.append(xi // order)
    O = [(a - b) // 2 for a, b in zip(D, Tw)]
    E = [(a + b) // 2 for a, b in zip(D, Tw)]
    assert all((a - b) % 2 == 0 for a, b in zip(D, Tw))
    per_r = lcm(*[N // np.gcd.reduce(np.array(e + (N,))) for e in cr])
    per_t = lcm(*[N // np.gcd.reduce(np.array(e + (N,))) for e in ct])
    return dict(d=d, group_order=order, N=N, eigen_classes=len(cr), periods=[per_r, per_t], total=D, twisted=Tw,
                even=E, odd=O, lowest_odd_degree=next((k for k, v in enumerate(O) if v), None),
                design_strength=next(k for k in range(1, K + 1) if D[k] > 1) - 1)


def two_qutrit_classes():
    """eigenvalue-multiset classes of the projective TWO-qutrit Clifford group {W(a) V_M}: 81 x |Sp(4,3)| = 4,199,040
    elements, streamed one Weyl translation at a time (V_M from Pass 11350's Weil construction)"""
    import w33_pass11330_orbit_census as O
    import w33_pass11350_linear_decider as L
    D = L.Decider(2)
    Ms, _ = O.all_symplectic(D.wl)
    V = np.array([D.weil(np.array(M) % 3) for M in Ms])
    assert len(V) == 51840
    cands = [n for n in range(1, 6481) if 6480 % n == 0]
    probe = V[:2000]
    lam = np.linalg.eigvals(probe)
    tw = np.linalg.eigvals(np.conj(probe) @ probe)
    N = find_N([np.angle(lam / lam[:, :1]).ravel(), np.angle(tw).ravel()], cands)
    cr, ct = Counter(), Counter()
    for a in range(81):
        G = np.einsum('ij,mjk->mik', D.wl.W[a], V)
        lam = np.linalg.eigvals(G)
        er = np.angle(lam / lam[:, :1]) * N / (2 * np.pi)
        et = np.angle(np.linalg.eigvals(np.conj(G) @ G)) * N / (2 * np.pi)
        if np.abs(er - np.rint(er)).max() > 1e-6 or np.abs(et - np.rint(et)).max() > 1e-6:
            N2 = find_N([np.angle(lam / lam[:, :1]).ravel(), np.angle(np.linalg.eigvals(np.conj(G) @ G)).ravel()], cands)
            raise AssertionError(("N too small", N, N2))
        er = np.sort(np.rint(er).astype(np.int64) % N, axis=1)
        et = np.sort(np.rint(et).astype(np.int64) % N, axis=1)
        for arr, cnt in ((er, cr), (et, ct)):
            u, c = np.unique(arr, axis=0, return_counts=True)
            for row, m in zip(u, c):
                cnt[tuple(int(x) for x in row)] += int(m)
    assert sum(cr.values()) == sum(ct.values()) == 81 * 51840
    return cr, ct, N


def two_qutrits(K=8):
    cr, ct, N = two_qutrit_classes()
    r = series_from_classes(cr, ct, N, 81 * 51840, K, 9)
    r["register"] = "two qutrits (C^3 x C^3), full two-qutrit Clifford group"
    return r


def restriction_to_one_qutrit(n_states=2, seed=11):
    """the unique degree-6 odd invariant of the two-qutrit register on psi (x) s, s a stabilizer state.  With O_6 = 1 the
    frame-potential difference Delta2(Psi) = avg_C |<Psi|C Psi>|^12 - avg_C |<Psi|C conj(Psi)>|^12 (all 4,199,040 projective
    two-qutrit Cliffords) is proportional to |h^(2)(Psi)|^2; compare with the one-qutrit Delta6 (Pass 11419: 2|h6|^2)."""
    import w33_pass11330_orbit_census as O
    import w33_pass11350_linear_decider as L
    D = L.Decider(2)
    V = np.array([D.weil(np.array(M) % 3) for M in O.all_symplectic(D.wl)[0]])
    CF1 = P7.clifford_group(3)

    def delta2(Psi):
        Psi = Psi / np.linalg.norm(Psi)
        VP, VQ = V @ Psi, V @ Psi.conj()
        a1 = sum(float((np.abs((Wa @ VP.T).T @ Psi.conj()) ** 12).sum()) for Wa in D.wl.W)
        a2 = sum(float((np.abs((Wa @ VQ.T).T @ Psi.conj()) ** 12).sum()) for Wa in D.wl.W)
        return (a1 - a2) / (81 * len(V))

    def delta1(psi):
        psi = psi / np.linalg.norm(psi)
        return float(np.mean(np.abs(CF1 @ psi @ psi.conj()) ** 12) - np.mean(np.abs(CF1 @ psi.conj() @ psi.conj()) ** 12))

    rng = np.random.default_rng(seed)
    w = np.exp(2j * np.pi / 3)
    rows = []
    for stab in (np.array([1, 0, 0], complex), np.array([1, w, w], complex) / np.sqrt(3)):
        for _ in range(n_states):
            v = rng.normal(size=3) + 1j * rng.normal(size=3)
            d1 = delta1(v)
            rows.append(dict(ratio_psi_s=delta2(np.kron(v, stab)) / d1, ratio_s_psi=delta2(np.kron(stab, v)) / d1,
                             delta1=d1))
    prod = []
    for _ in range(2):
        a = rng.normal(size=3) + 1j * rng.normal(size=3)
        b = rng.normal(size=3) + 1j * rng.normal(size=3)
        prod.append(dict(delta2=delta2(np.kron(a, b)), delta1_a=delta1(a), delta1_b=delta1(b)))
    return dict(rows=rows, max_dev_from_1_over_108=max(abs(108 * r[k] - 1) for r in rows for k in ("ratio_psi_s", "ratio_s_psi")),
                generic_products=prod)


def derive_108(seed=5):
    """DERIVATION of Delta2(psi x s) = Delta6(psi)/108.  Compress every two-qutrit Clifford C to the slot psi x |0>:
    A_C = (1 x <0|) C (1 x |0>).  Then <Psi|C Psi> = <psi|A_C psi> and <Psi|C conj Psi> = <psi|A_C conj psi>.  Counted over all
    4,199,040 projective Cliffords: A_C is a unitary on 1/120 of them, 3^(-1/2) x unitary on 27/40, rank one on the rest; the
    induced unitaries cover the one-qutrit Clifford group uniformly, and the rank-one terms cancel in the odd difference.
    Hence Delta2 = (1/120 + (27/40) 3^(-6)) Delta6 = (1/120 + 1/1080) Delta6 = Delta6/108."""
    from fractions import Fraction
    import w33_pass11330_orbit_census as O
    import w33_pass11350_linear_decider as L
    D = L.Decider(2)
    Ms = np.array(O.all_symplectic(D.wl)[0]) % 3
    CF1 = P7.clifford_group(3)
    Iso = np.kron(np.eye(3), np.array([[1.0], [0.0], [0.0]]))
    rng = np.random.default_rng(seed)
    psi = rng.normal(size=3) + 1j * rng.normal(size=3)
    psi /= np.linalg.norm(psi)
    counts = Counter()
    hits = {1.0: np.zeros(216, np.int64), 3 ** -0.5: np.zeros(216, np.int64)}
    rank_one_odd = 0.0
    rank_one_scale = 0.0
    for M in Ms:
        V = D.weil(M)
        A = np.einsum('ip,apq,qj->aij', Iso.T, np.einsum('apq,qr->apr', D.wl.W, V), Iso)     # 81 compressions
        sv = np.linalg.svd(A, compute_uv=False)
        for Ai, s_ in zip(A, sv):
            if s_[0] > 1e-9 and abs(s_[2] - s_[0]) < 1e-9:
                key = 1.0 if abs(s_[0] - 1) < 1e-9 else 3 ** -0.5
                assert abs(s_[0] - key) < 1e-9
                U = Ai / s_[0]
                k = int(np.argmax(np.abs(np.einsum('cij,ij->c', CF1.conj(), U))))
                assert abs(abs(np.trace(CF1[k].conj().T @ U)) - 3) < 1e-8          # U is a one-qutrit Clifford
                hits[key][k] += 1
                counts["unitary" if key == 1.0 else "unitary / sqrt3"] += 1
            elif s_[0] > 1e-9:
                assert s_[1] < 1e-9
                counts["rank one"] += 1
                x, y = abs(np.vdot(psi, Ai @ psi)) ** 12, abs(np.vdot(psi, Ai @ psi.conj())) ** 12
                rank_one_odd += x - y
                rank_one_scale += x
            else:
                counts["zero"] += 1
    n = 81 * len(Ms)
    f1 = Fraction(counts["unitary"], n)
    f3 = Fraction(counts["unitary / sqrt3"], n)
    return dict(counts=dict(counts), fraction_unitary=str(f1), fraction_unitary_over_sqrt3=str(f3),
                uniform_over_Cl1={("1" if k == 1.0 else "1/sqrt3"): bool(v.min() == v.max()) for k, v in hits.items()},
                rank_one_odd_sum_relative=rank_one_odd / rank_one_scale,
                predicted=str(f1 + f3 * Fraction(1, 729)), equals_1_over_108=(f1 + f3 * Fraction(1, 729)) == Fraction(1, 108))


def closed_form(seq, d, period, fit=36):
    """rational generating function, fitted on the first `fit` terms and PROVED on all of them when the quasi-polynomial
    bound (2d - 1) * period <= len(seq) holds (coefficients of degree <= 2d - 2, period `period`)"""
    import sympy as sp
    from sympy.concrete.guess import guess_generating_function_rational
    t = sp.symbols("t")
    m = next(i for i, x in enumerate(seq) if x)                  # the guesser needs a nonzero first term
    g = guess_generating_function_rational([int(x) for x in seq[m:m + fit]], X=t)
    if g is None:
        return None
    g = t ** m * g
    num, den = sp.fraction(sp.together(g))
    ser = sp.Poly(sp.series(num / den, t, 0, len(seq)).removeO(), t).all_coeffs()[::-1]
    ser = [int(x) for x in ser] + [0] * (len(seq) - len(ser))
    roots_ok = all(sp.simplify(r ** period - 1) == 0 for r in sp.roots(sp.Poly(den, t)))
    proper = sp.degree(num, t) < sp.degree(den, t)
    return dict(form=str(sp.factor(g)), matches=ser == list(seq), proved=bool(ser == list(seq) and roots_ok and proper
                                                                              and (2 * d - 1) * period <= len(seq)))


def run(K=None):
    if "--derive-108" in sys.argv:
        res = json.load(open(OUT))
        res["derivation_108"] = derive_108()
        print(res["derivation_108"], flush=True)
        return res
    res = dict(pass_id=11492)
    for d, k in ((2, 72), (3, 72), (5, 30), (7, 16), (11, 9)):
        r = series(d, K or k)
        if d in (2, 3):
            P = lcm(*r["periods"])
            r["closed_forms"] = {nm: closed_form(r[nm], d, P, fit=30 if d == 2 else 40) for nm in ("odd", "total", "twisted")}
            print(r["closed_forms"], flush=True)
        res[f"d{d}"] = r
        print(d, r["group_order"], "N", r["N"], "lowest odd", r["lowest_odd_degree"], "odd", r["odd"][:16], flush=True)
    res["lowest_odd_degree_by_d"] = {d: res[f"d{d}"]["lowest_odd_degree"] for d in (2, 3, 5, 7, 11)}
    res["count_at_lowest_odd_degree"] = {d: res[f"d{d}"]["odd"][res[f"d{d}"]["lowest_odd_degree"]] for d in (2, 3, 5, 7, 11)}
    if "--skip-n2" not in sys.argv:
        res["two_qutrits"] = two_qutrits()
        print("two qutrits", res["two_qutrits"]["lowest_odd_degree"], res["two_qutrits"]["odd"], res["two_qutrits"]["total"],
              flush=True)
    if "--skip-n2" not in sys.argv:
        res["two_qutrit_restriction"] = restriction_to_one_qutrit()
        print("restriction", res["two_qutrit_restriction"]["max_dev_from_1_over_108"], flush=True)
        res["derivation_108"] = derive_108()
    res["control_p5_two_constructions_agree"] = series(5, 12, G=P7.clifford_group(5))["total"] == res["d5"]["total"][:13]
    res["jarlskog_degree_by_d_Pass11357"] = {2: 10, 3: 6, 5: 4}
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
