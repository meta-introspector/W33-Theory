"""Pass 11371: the flavour Jarlskog of a substrate tick, and why it cannot see the arrow of time; the arrow as the
Altland-Zirnbauer class of the tick (az_census, run with --az).

The map (answering the open 'explicit full-rank flavour map' of the Codex flavour audit).  A one-qutrit tick U is a
3 x 3 unitary.  Read it as a mixing matrix between the clock basis (where the magic gate T is diagonal, the 'mass basis')
and its image:  H_u = diag(m_u),  H_d = U diag(m_d) U^dag.  Then the flavour CP invariant is
        det[H_u, H_d] = 2 i J_CKM(U) prod_{i<j} (m_u,i - m_u,j)(m_d,i - m_d,j)    (checked numerically, sign included),
        J_CKM(U) = Im(U_11 U_22 conj(U_12) conj(U_21)).

THEOREM 1 (exact, elementary).  J_CKM(U^T) = J_CKM(U),  J_CKM(conj U) = J_CKM(U^dag) = -J_CKM(U), and
J_CKM(C U C^dag) = J_CKM(U) for every MONOMIAL C (phases times a permutation sigma: rows and columns move by the same
sigma, sgn(sigma)^2 = 1).  Hence
  * flavour CP is EVEN under the substrate time reversal U -> U^T: no function of J_CKM can witness the arrow;
  * if U has an antiunitary symmetry  conj U ~ C U C^dag  or a unitary reversal  U^dag ~ C U C^dag  with C in the
    Borel (clock-basis-preserving) Clifford subgroup (54 of 216), then J_CKM(U) = 0.

LEMMA 2 (exact).  {1, T: U -> U^T, K: U -> conj U, I: U -> U^dag} is a Klein four-group, and the set of operations X with
X(U) Clifford-conjugate to U (up to phase) is a SUBGROUP of it (the Clifford group is closed under conj and C^T =
conj(C)^-1).  So 'any two of {reversible, antiunitary-symmetric, unitarily reversible} imply the third'; a tick has one
of five symmetry types.

COMPUTED (exact Clifford search on every coset-reduced word with 1, 2, 3 cubic gates, Pass 11312):
  * the joint distribution of symmetry type and J_CKM;
  * whether the substrate arrow and flavour CP are correlated.
"""

from __future__ import annotations

import itertools
import json
import sys
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11213_cubic_t_violation as P1  # noqa: E402
import w33_pass11312_depth_law as DL  # noqa: E402
import w33_pass11355_substrate_jarlskog as SJ  # noqa: E402

OUT = ROOT / "data" / "w33_pass11371_flavour_vs_arrow.json"
CF = SJ.CF
JMAX = 1 / (6 * np.sqrt(3))
TOL = 1e-9
_S = {}


def jcp(U):
    return float(np.imag(U[0, 0] * U[1, 1] * np.conj(U[0, 1]) * np.conj(U[1, 0])))


def monomial(C):
    return bool(((np.abs(C) > 1e-9).sum(axis=0) == 1).all())


BOREL = np.array([monomial(C) for C in CF])


def symmetric_via(U, V):
    """indices of Cliffords C with V ~ C U C^dag (up to phase)"""
    X = np.einsum('cij,jk,clk->cil', CF, U, CF.conj())
    ov = np.abs(np.einsum('cij,ij->c', X.conj(), V)) / 3
    return np.flatnonzero(ov > 1 - 1e-9)


def classify(U):
    T = symmetric_via(U, U.T)
    K = symmetric_via(U, U.conj())
    I = symmetric_via(U, U.conj().T)
    typ = "".join(s for s, x in (("T", T), ("K", K), ("I", I)) if len(x)) or "none"
    return typ, bool(BOREL[K].any()), bool(BOREL[I].any()), bool(BOREL[T].any())


def flavour_identity(rng, trials=20):
    dev = 0.0
    for _ in range(trials):
        U = SJ.det1(CF[rng.integers(216)] @ P1.T1 @ CF[rng.integers(216)] @ P1.T1)
        mu, md = rng.normal(size=3), rng.normal(size=3)
        Hu, Hd = np.diag(mu), U @ np.diag(md) @ U.conj().T
        lhs = np.linalg.det(Hu @ Hd - Hd @ Hu)
        pr = np.prod([(mu[i] - mu[j]) * (md[i] - md[j]) for i, j in ((0, 1), (0, 2), (1, 2))])
        dev = max(dev, abs(lhs - 2j * jcp(U) * pr))
    return dev


def theorem1_checks(rng, trials=200):
    dev = 0.0
    for _ in range(trials):
        Z = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
        U, _ = np.linalg.qr(Z)
        C = CF[rng.choice(np.flatnonzero(BOREL))]
        dev = max(dev, abs(jcp(U.T) - jcp(U)), abs(jcp(U.conj()) + jcp(U)), abs(jcp(U.conj().T) + jcp(U)),
                  abs(jcp(C @ U @ C.conj().T) - jcp(U)))
    return dev


def _init():
    _S["R"] = DL.coset_reps(list(CF))


def _block(args):
    k, prefix = args
    Rr, T = _S["R"], P1.T1
    cnt = Counter()
    for rest in itertools.product(range(24), repeat=k - 1 - len(prefix)):
        W = np.eye(3, dtype=complex)
        for r in list(prefix) + list(rest):
            W = W @ Rr[r] @ T
        for C in CF:
            U = W @ C @ T
            typ, kb, ib, tb = classify(U)
            j = jcp(U)
            jl = "0" if abs(j) < TOL else ("max" if abs(abs(j) - JMAX) < 1e-9 else "other")
            # Theorem 1 consequence, checked on every word
            assert not ((kb or ib) and jl != "0"), "Borel-realised K or I symmetry with J_CKM != 0"
            cnt[(typ, jl, kb or ib)] += 1
    return cnt


def run():
    rng = np.random.default_rng(11371)
    res = dict(pass_id=11371)
    res["det_commutator_identity_max_dev"] = flavour_identity(rng)
    res["theorem1_max_dev"] = theorem1_checks(rng)
    res["borel_cliffords"] = int(BOREL.sum())
    res["clifford_J_CKM_values"] = dict(Counter(
        "0" if abs(jcp(C)) < TOL else ("+max" if jcp(C) > 0 else "-max") if abs(abs(jcp(C)) - JMAX) < 1e-9 else "other"
        for C in CF))
    res["nonmonomial_cliffords_all_maximal_J"] = all(abs(abs(jcp(C)) - JMAX) < 1e-9 for C in CF[~BOREL])
    out = {}
    with Pool(11, initializer=_init) as pool:
        for k in (1, 2, 3):
            jobs = [(1, ())] if k == 1 else [(k, (i,)) for i in range(24)]
            tot = Counter()
            for c in pool.map(_block, jobs):
                tot.update(c)
            n = sum(tot.values())
            table = Counter()
            for (typ, jl, b), v in tot.items():
                table[f"{typ} | J_CKM={jl}"] += v
            rev = {"0": Counter(), "max": Counter(), "other": Counter()}
            for (typ, jl, b), v in tot.items():
                rev[jl]["reversible" if "T" in typ else "violating"] += v
            out[str(k)] = dict(words=n, types=dict(sorted(table.items())),
                               arrow_by_flavour_cp={jl: dict(c) for jl, c in rev.items()},
                               borel_realised_K_or_I=sum(v for (typ, jl, b), v in tot.items() if b),
                               types_seen=sorted({typ for (typ, jl, b) in tot}))
            print(k, json.dumps(out[str(k)]), flush=True)
    res["words"] = out
    res["only_subgroup_types"] = all(set(v["types_seen"]) <= {"none", "T", "K", "I", "TKI"} for v in out.values())
    return res


AZ = {"none": "A", "T": "AI", "I": "AIII", "K": "D", "TKI": "BDI"}


def az_census(kmax=2):
    """ALTLAND-ZIRNBAUER READING.  With U = exp(-iHt):
         antiunitary Theta, Theta U Theta^-1 ~ U^-1   <=>  Theta H Theta^-1 = +H   (time reversal T: the arrow channel)
         antiunitary Theta, Theta U Theta^-1 ~ U      <=>  Theta H Theta^-1 = -H   (particle-hole C: the K channel)
         unitary S,         S U S^-1 ~ U^-1           <=>  S H S^-1 = -H           (chiral S: the I channel)
    and S = T C is Lemma 2.  From U^T ~ C U C^dag the reversal is Theta = conj(C) K with Theta^2 = conj(C) C; from
    conj(U) ~ C U C^dag the particle-hole operator is C^dag K with square (conj(C) C)^dag.  In odd dimension a scalar
    square must be +1 (det(V conj V) = |det V|^2 > 0).  So the five symmetry types are the AZ classes
    A, AI, AIII, D, BDI.  Checked here on every coset-reduced word with up to kmax cubic gates:
      * the square of each Clifford-realised T and C (scalar +1, or a non-scalar unitary = an extra symmetry);
      * chiral spectral pairing e^{i th} -> lambda e^{-i th} for every chiral (I-type) tick."""
    import w33_pass11312_depth_law as DL
    Rr = DL.coset_reps(list(CF))
    out = {}
    for k in range(1, kmax + 1):
        words = [C @ P1.T1 for C in CF] if k == 1 else [Rr[r] @ P1.T1 @ C @ P1.T1 for r in range(24) for C in CF]
        cnt = Counter()
        squares = Counter()
        pairing_fail = 0
        rev_with_plus = rev = 0
        for U in words:
            typ = classify(U)[0]
            if "T" in typ:
                rev += 1
                rev_with_plus += any(np.abs(CF[c].conj() @ CF[c] - np.eye(3)).max() < 1e-9
                                     for c in symmetric_via(U, U.T))
            for ch, V in (("T", U.T), ("C", U.conj())):
                for c in symmetric_via(U, V):
                    M = CF[c].conj() @ CF[c]
                    if np.abs(M - M[0, 0] * np.eye(3)).max() < 1e-9:
                        squares[f"{ch}^2 = {M[0, 0].real:+.0f}"] += 1
                    else:
                        squares[f"{ch}^2 non-scalar"] += 1
            for c in symmetric_via(U, U.conj().T):
                S = CF[c]
                X = S @ U @ S.conj().T
                lam = np.vdot(U.conj().T, X) / 3
                ev = np.linalg.eigvals(U)
                img = lam * np.conj(ev)
                pairing_fail += max(min(abs(e - i) for i in img) for e in ev) > 1e-9
            cnt[AZ[typ]] += 1
        out[str(k)] = dict(words=len(words), az_classes=dict(cnt), squares=dict(squares),
                           chiral_pairing_failures=int(pairing_fail), reversible_words=rev,
                           reversible_words_with_a_reversal_squaring_to_plus_one=int(rev_with_plus))
        print(k, out[str(k)], flush=True)
    return out


def spectrally_paired(U, tol=1e-7):
    """is the eigenphase multiset symmetric under theta -> c - theta for some c (necessary for a chiral symmetry)?"""
    th = np.angle(np.linalg.eigvals(U))
    for i in range(len(th)):
        for j in range(len(th)):
            img = np.angle(np.exp(1j * (th[i] + th[j] - th)))
            if all(min(abs(np.angle(np.exp(1j * (x - y)))) for y in th) < tol for x in img):
                return True
    return False


def two_qutrit_chirality(samples=3000, seed=5):
    """two qutrits, one magic gate: is every violator chiral?  Exact verdicts from the F3-linear decider (Pass 11350)
    on uniformly sampled (class, frame); chirality tested by the necessary spectral-pairing condition"""
    import w33_pass11330_orbit_census as O
    import w33_pass11350_linear_decider as L
    D = L.Decider(2)
    Ms, _ = O.all_symplectic(D.wl)
    rng = np.random.default_rng(seed)
    cnt = Counter()
    for _ in range(samples):
        i = int(rng.integers(len(Ms)))
        g = D.good_frames(Ms[i])
        if g is None:
            continue
        a = int(rng.integers(81))
        U = D.wl.W[a] @ D.weil(Ms[i]) @ D.T1
        cnt[f"{'violating' if not g[a] else 'reversible'}, {'paired' if spectrally_paired(U) else 'unpaired'}"] += 1
    return dict(cnt)


def main():
    if "--az2" in sys.argv:
        res = json.load(open(OUT))
        res["two_qutrit_one_gate_chirality"] = two_qutrit_chirality()
        print(res["two_qutrit_one_gate_chirality"])
        json.dump(res, open(OUT, "w"), indent=1, default=str)
        return
    if "--az" in sys.argv:
        res = json.load(open(OUT))
        res["altland_zirnbauer"] = az_census()
        json.dump(res, open(OUT, "w"), indent=1, default=str)
        return
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
