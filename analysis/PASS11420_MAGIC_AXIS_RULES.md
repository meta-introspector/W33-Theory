# Pass 11420 — the magic-axis rules have one law behind them: when M² fixes the magic axis, every frame not orthogonal to the M-fixed vector z₁ + Mz₁ breaks the arrow

Producer: `analysis/w33_pass11420_magic_axis_rules.py`
Certificate: `data/w33_pass11420_magic_axis_rules.json`
Regression: `tests/test_w33_pass11418_11422.py`

**Background.** Pass 11373 found four magic-axis rules, exact at n = 2 and 3. Pass 11421 tested them at n = 4 on
every decided targeted class.
* Mz₁ = z₁: some frame violates.
* Mz₁ = −z₁: every frame is reversible.
* Same line, M²z₁ = z₁: some frame violates.
* Same line, M²z₁ = −z₁: every frame is reversible.

## Conjecture F′ (computer-verified at n = 2, 3, 4): one law for both "bad" cells

> **If M²z₁ = z₁, then v = z₁ + Mz₁ is fixed by M, and a frame a violates whenever ω(a, v) ≠ 0.**

* **Fixed axis.** v = 2z₁ = −z₁, and the condition is a_{x₁} ≠ 0 (2/3 of all frames).
* **Same line.** v is a new M-fixed vector, again giving 2/3 of the frames.
* **Reversed axis (Mz₁ = −z₁).** v = 0 and the law predicts nothing, consistent with that cell being all good.

| n | Mz₁ = z₁ | same line, M²z₁ = z₁ |
|---|---|---|
| 2 (all classes) | 648/648 (624 with equality, 24 with extra bad frames) | **648/648, all with equality** |
| 3 (every orbit of Pass 11373) | 97 fast + 29 slow orbits, no exception | 75 fast + 6 slow orbits, no exception |
| 4 (targeted, sparse decider) | 18/18 with equality (2 undecided) | 16/16 with equality (4 undecided) |

* "With equality" means the violating frames are **exactly** {a : ω(a, v) ≠ 0}.
* The slow n = 3 classes are decided on frame orbits of C_H(M). C_H(M) fixes z₁ and Mz₁, so it preserves ω(a, v)
  (Lemma L4), and every orbit lies entirely inside or outside the predicate.

## Proved lemmas

* **L1.** Mz₁ = ±z₁ ⟺ V_M T₁ V_M† = T₁^{±1}. The map k ↦ k³ is odd. Checked on all 1296 such two-qutrit classes.
* **L2.** If Mz₁ = z₁, then U³ is Clifford for every frame.
  * For a_{x₁} = 0, everything commutes with T₁, and T₁³ = Z₁.
  * For a_{x₁} = c ≠ 0, U cycles the Z₁-eigenspaces and the cubic phases sum to
    Σⱼ(k + jc)³ ≡ 3k³ + 6kc² (mod 9). Then ζ raised to that sum is ω^{k(1 + 2c²)} = 1.
  * Checked on 2592 frames. A violating fixed-axis tick is a cube root of a Clifford.
* **L3.** If Mz₁ = −z₁ and dim ker(M + I) = 1, there is a reversing anti-symplectic involution σ with σz₁ = −z₁.
  This gives the solution Q = σJ of equation (S) at k = 0.
  * Proof: Wonenburger's τ reverses M. Then τz₁ lies in ker(M⁻¹ + I) = ker(M + I), so τz₁ = ±z₁. Take ±τ.
  * The condition holds for 414 of the 648 reversed classes at n = 2, and for 8 049 618 of 12 597 120 (64%) at
    n = 3.
  * At n = 2 such a σ exists for all 648 (Pass 11373, computed).
* **L4.** C_H(M) preserves ω(a, z₁ + Mz₁).

## Status of the four rules

* **"Bad" rules.** They follow from F′, which is computer-verified at n = 2 and 3 exactly and at n = 4 on samples. F′
  itself is not proved.
* **"Good" rules.** They reduce to an existence lemma, proved by L3 in the generic case, plus a frame-covering lemma
  (open).

**Open.**
* A proof of F′. The fixed vector z₁ + Mz₁ suggests an obstruction carried by the Pauli W(z₁ + Mz₁), which commutes
  with V_M.
* The extra bad frames of the 24 two-qutrit fixed-axis classes.
