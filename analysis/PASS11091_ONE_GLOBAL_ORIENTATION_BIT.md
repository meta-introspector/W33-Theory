# Pass 11091 — the substrate has exactly one global orientation bit, and it is the chirality of Pass 346

Producer: `analysis/w33_pass11091_one_global_orientation_bit.py`
Certificate: `data/w33_pass11091_one_global_orientation_bit.json`
Regression: `tests/test_w33_pass11091_one_global_orientation_bit.py`
Cross-track join. It cites:
* **Pass 346** (this track): the chirality no-go;
* **Passes 11072, 11073, 11084, 11085** (the other track): the internal Borel deck involution, the outer
  similitude and the two-bit orientation codec;
* **Holotrade e92a047**: frame chirality is reversed only by antiunitaries.

## Prior art, stated first

The Forty Points paper's **Corollary 9.7** already shows the following. In the stabiliser of a Bell line,
the elements outside PSp(4,3) form one coset. It is exactly the set that reverses the history graph's
orientation (Theorem 4.4), the set that acts antiunitarily (Pass 5730), and the set that swaps the
half-spins (Theorem 9.5 / e92a047). This pass does not re-derive that. It adds three things:
1. the other track's cubic-tick chirality χ is the **same** coset;
2. a **uniqueness** theorem, which makes the coincidence forced rather than observed;
3. the other track's chamber-sheet sign σ is **not** a global orientation.

## The question

The other track found **two** Z₂ orientation labels over every W(3,3) flag:
* **σ**, the tetracode or chamber-sheet sign. It is flipped by the internal Borel deck involution
  T = diag(1,−1,1,−1) ∈ Sp(4,3), and by the outer similitude S = diag(1,−1,−1,1), which has
  SᵀJS = −J.
* **χ**, the cubic-tick chirality. It is flipped by S and by TS, not by T.

Pass 346 proved that the half-spin chirality S⁺ ↔ S⁻ of the D5 structure is exchanged by the substrate's
outer controller, which has det = −1 on the D5 module. So no invariant of the full automorphism group
selects it. Are these the same bit, or three?

## The theorem (checked by explicit construction of GSp(4,3))

1. **Orders.** |Sp(4,3)| = 51,840 and |GSp(4,3)| = 103,680. With the scalars {I, 2I} removed,
   |PGSp(4,3)| = 51,840 = |W(E6)|.
2. **The multiplier character.** The similitude multiplier μ (MᵀJM = μJ) is a homomorphism GSp(4,3) →
   {±1}. It is surjective and trivial on scalars, so it descends to PGSp(4,3) with kernel PSp(4,3), of
   order 25,920.
3. **Uniqueness.** Sp(4,3) is **perfect**: its commutators generate all 51,840 elements. Hence PSp(4,3) is
   perfect and the abelianization of PGSp(4,3) is Z₂. **μ is the only nontrivial homomorphism
   PGSp(4,3) → Z₂.**
4. **The two codec bits.** On the other track's deck group, μ = (1: +1, T: +1, S: −1, TS: −1). That is
   exactly their χ. σ is flipped by T, which has μ = +1, so **σ is not the restriction of any global
   character**. It is a label defined relative to a flag.

Four Z₂ data are each a nontrivial homomorphism on PGSp(4,3) ≅ W(E6):
* the half-spin swap of Pass 346;
* the cubic-tick reversal of Pass 11073;
* the sign character of W(E6), i.e. det on the reflection representation;
* antiunitarity of the frame action (e92a047).

By step 3 they are **one and the same character, μ**. Corollary 9.7 established the coincidence for three
of these, one coset at a time. Step 3 makes it a necessity: any Z₂-valued orientation of the substrate that
is a group homomorphism, including ones not yet found, is either trivial or μ.

## Reading

* **Exactly one global orientation bit.** The substrate carries one, μ: the orientation of the symplectic
  form. It is also the choice of primitive cube root of unity (the coefficient conjugation inside Pass
  333's T), the handedness of the half-spins, the direction of the cubic tick, and unitary versus
  antiunitary.
* **Pass 346's no-go transfers to χ.** No datum invariant under the full automorphism group selects the
  other track's tick chirality χ. Selecting it means choosing μ, i.e. breaking PGSp(4,3) to PSp(4,3).
* **σ is local, not a second global orientation.** Pass 11085's statement that "the temporal orientation
  requires two Z₂ labels" is correct *per flag*. Only one of the two labels, χ = μ, is a property of the
  substrate as a whole; σ depends on the chosen flag.

This sharpens the paper's open problem 4 ("chirality from outside"). The one bit that must come from
outside is μ, and it is simultaneously handedness, tick direction and time reversal.

Rediscovery guard. The hook flagged the compound "e6 + tetracode" against BT927 (E8-lift artifact
reconciliation) and Passes 9185–9196 (Golay/tetracode glue bifurcation). Both were read. Neither concerns
orientation characters, chirality or the similitude multiplier; the overlap is vocabulary only.
