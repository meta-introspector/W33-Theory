# Pass 11213 — T-violation enters exactly with the cubic layer: the cubic phase times a cyclic shift

Producer: `analysis/w33_pass11213_cubic_t_violation.py`
Certificate: `data/w33_pass11213_cubic_t_violation.json`
Regression: `tests/test_w33_pass11213_11216.py`

**Question (from Pass 11206).** Every Clifford tick is inverted by a substrate time reversal (an anti-unitary Clifford
Θ = V K; Wonenburger 1966), so the Clifford skeleton has no T-violation. Every unitary is inverted by *some*
anti-unitary, because U^* and U^† are unitarily similar. The physical question is therefore whether the substrate's
own time reversals still suffice once the non-Clifford resource is switched on. That resource is the cubic phase gate
T = diag(1, ζ₉, ζ₉⁸) = ζ₉^{x³}, the qutrit form of the degree-3 (E6-cubic) layer of the universal gate set.

**Result (single qutrit; all 216 Cliffords mod phase; 3000 random words per depth).**

| cubic gates d | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|---|
| fraction with **no** substrate time reversal | 0 | 1/12 (exact) | 0.40 | 0.55 | 0.69 | 0.78 | 0.84 | 0.90 | 0.92 |

* **The Clifford skeleton is T-symmetric; the cubic layer breaks T generically.** The fraction climbs towards 1 with
  depth.
* **The minimal T-violating ticks are exactly T·X^a Z^b S^c with a ≠ 0.** These are 18 of the 216 one-cubic-gate
  ticks: the cubic phase composed with a **cyclic shift of the three levels**. The cubic phase alone is inverted by
  complex conjugation. So is any cubic phase dressed only by diagonal Cliffords.
* All 18 minimal violators have the same Jarlskog-type value J = Im tr(U Z U† X U Z† U† X†) = 3√3/2. J is **not**
  a clean witness, however: it is also non-zero on many invertible words.

**Reading.** In the Standard Model, CP (equivalently T) violation needs a complex phase *and* mixing among three
generations. On the substrate, T-violation needs the cubic phase ζ₉^{x³} *and* a cyclic permutation of the three
levels. Neither alone suffices. This is a structural analogy on one qutrit, not a derivation of the CKM matrix.

**Scope.** Single qutrit, time reversals restricted to anti-unitary Cliffords. The two-qutrit extension (anti-unitary
two-qutrit Clifford group, 51 840 × 81 elements) is not run here.
