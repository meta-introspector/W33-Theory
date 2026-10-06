# Pass 11531 — why the cubed-phase test is exact on words but not on PU(3): only T-count ≤ 1 words ever carry a non-affine relation

Producer: `analysis/w33_pass11531_cubed_test_words_vs_pu3.py`
Certificate: `data/w33_pass11531_cubed_test_words_vs_pu3.json`
Regression: `tests/test_w33_pass11531_11538.py`

**Setting.**
* A *relation* (L, f) for U means c(Lp) = μ·ω^{f(p)}·c(p) at every p, where L is anti-symplectic on F₃² and c the Weyl
  coefficients.
* U is reversible iff some relation has f affine on the support (Pass 11252).
* The cubed-phase test only sees that some relation exists, whatever f is.
* Pass 11512 found non-reversible unitaries passing it.
* Passes 11500 and 11513 found it exact on every word through depth 10.

## (A) The counterexamples form a 3-parameter family

* Take L = −swap, (a, b) ↦ (−b, −a), and the quadratic f(a, b) = δ_{b,2} + 2δ_{a,2} = a² − b² + 2a + b (mod 3), or its
  double.
* The unitaries in the 6-dimensional eigenspace form a smooth family of real dimension **4**: Jacobian rank, at 599 of
  600 sampled points. That is a **3-parameter family in PU(3)** after the global phase.
* **All 600 sampled points are non-reversible.** Their maximal Clifford overlap ranges from 2.39 to just below 3, so the
  family approaches the reversible set without meeting it.

## (B) The words: non-affine relations occur only at T-count ≤ 1

Every distinct one-qutrit Clifford+T operator of depth ≤ 7 was scanned: 5 567 184 operators at depth 7, all 24 L's.

| depth | distinct operators | affine relations | non-affine: Clifford | non-affine: T-count 1 | non-affine: T-count ≥ 2 |
|---|---|---|---|---|---|
| 1 | 216 | 270 | 0 | 324 | **0** |
| 3 | 5 616 | 3 567 | 1 350 | 1 314 | **0** |
| 5 | 166 752 | 20 748 | 2 376 | 2 556 | **0** |
| 6 | 959 040 | 48 426 | 2 376 | 2 592 | **0** |
| 7 | 5 567 184 | 120 423 | 2 376 | 2 592 | **0** |

* Non-affine relations on words occur **only** in two kinds, and both have full Weyl support:
  * **Cliffords:** 2376, a constant count. Their coefficients have flat moduli and quadratic phases, so many L satisfy a
    quadratic relation.
  * **T-count-1 operators C₁T^{±1}C₂:** the count saturates at 2592 = 216 × 12.
* Both kinds are reversible, through an affine relation as well.
* **No operator of T-count ≥ 2 carries a non-affine relation**, through depth 7.

**The mechanism.**
* The cubed test discards only one thing: whether the relation's phase is affine.
* That matters only for a U whose relations are all non-affine.
* On words such relations exist only for Cliffords and single-T operators. Those always also have an affine relation, so
  the discarded information never decides anything.
* The PU(3) counterexample family consists of non-flat, full-support unitaries. A word reaches it only by having T-count
  ≥ 2 and a non-affine relation, and that never happens.

**Scope.**
* "T-count ≥ 2 never carries a non-affine relation" is verified exhaustively through depth 7. It is not proved for all
  depths.
* The natural route to a proof is the arithmetic of the coefficients in ℤ[ζ₉, 1/√−3]. It is not attempted here.
* The saturation at depth 6 strongly suggests the set of relation-carrying words is finite up to depth: all of them have
  T-count ≤ 1.
