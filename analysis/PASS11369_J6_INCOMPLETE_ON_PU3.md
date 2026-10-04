# Pass 11369 — J₆ is complete on Clifford+T words but NOT on PU(3); the hierarchy is complete with a finite degree, and J₈ is the conjectured complete witness

Producer: `analysis/w33_pass11369_j6_completeness.py`
Certificate: `data/w33_pass11369_j6_completeness.json`
Regression: `tests/test_w33_pass11369_11373.py`

**Question (Pass 11355 next step).** Prove that J₆ is a complete arrow witness for one qutrit: J₆(U) = 0 ⟺ U is
reversible ⟺ Uᵀ is Clifford-conjugate to U.

**Answer: false on PU(3).** Pass 11355's finding is unaffected, since it was scoped to words: J₆ is exactly complete
on all one-qutrit words with up to 4 cubic gates.

## Theorem A — the witness hierarchy is complete, with a finite degree

* Put P(U) = U ⊗ Ū. The projective Clifford group, of order 216, acts **linearly** on P by conjugation, and
  P(Uᵀ) = P(U)ᵀ.
* Invariants of a finite group separate its orbits, and a separating set exists in degree ≤ |G| (Derksen–Kemper,
  *Computational Invariant Theory*; Kemper 2009).
* The linear invariant ⟨Φ|P|Φ⟩ = tr(UU†)/3 = 1 lifts every lower-degree invariant to degree t. So A_t(U) carries all
  invariants of bidegree ≤ (t,t).

> **U reversible ⟺ J_{2t}(U) = 0 for every t ≥ t\*, with 3 ≤ t\* ≤ 216.** Completeness is monotone in t.

## Finding 1 (computer-certified): t\* > 3. J₆ has a 4-dimensional family of non-reversible zeros

* **How it was found.** Gradient descent on J₆ from 2000 Haar-random starts reached 1936 zeros.
  * **313** of them are **certified non-reversible**: their distance to the reversible set is up to 0.86, and the next
    witness J₈ is > 10⁻⁶ there, up to 4.7.
  * Any nonzero J_{2t} certifies a violator, because every J_{2t} vanishes at reversible points.
* **One point refined to rounding level.**
  * Gauss–Newton on the full 729 × 729 residual F = A₃(U) − A₃(Uᵀ) converges quadratically: ‖F‖ goes
    4×10⁻⁸ → 7×10⁻¹³ → 8×10⁻¹⁵.
  * There J₈ = 0.305, J₁₀ = 9.29 and the distance to the reversible set is 0.316.
  * dF has rank exactly 4 there (singular values 1, 0.36, 0.17, 0.026, then below 10⁻⁸). So the spurious zero set is
    locally a smooth **4-dimensional** submanifold of the 8-dimensional PU(3).
  * The point has no antiunitary symmetry and no unitary reversal either (Pass 11371 channels; best Clifford overlaps
    0.60). Its structure is open.
* **Two populations (scratch probe, 49 spurious points from fresh descents; not in the certificate).**
  * 32 have a generic spectrum, with the smallest eigenphase gap > 0.039 · 2π, like the refined point.
  * **17 (35%) have a nearly degenerate spectrum**, gap < 0.003 · 2π. For Haar-random unitaries that rate is 0.01%.
  * The low-dimensional reversible strata (below) also have a degenerate spectrum, and J₆ is flat there. So either a
    second spurious component exists, or the 4-dimensional family accumulates on the degenerate-spectrum locus.
    Unresolved.
* **Why it can happen.** At degree 3 there are only **7** independent T-odd invariants (Pass 11370). Their common
  zero set contains the reversible set but is larger.
* **64 runs** stopped with 10⁻¹² ≤ J₆ < 10⁻⁶. They are unconverged, not positive minima: no local minimum with
  J₆ > 10⁻⁶ was found (8 such runs for J₈).
* **26 further zeros** lie at distance > 10⁻³ but have J₈ < 10⁻⁶. They are unresolved near-misses: slow convergence
  onto the low-dimensional strata below is the likely cause, but this is not verified.

## Finding 2 (computer-verified, not a proof): J₈ shows no spurious zeros

* Gradient descent on J₈ from 1000 Haar-random starts: 992 reached zeros.
* None is certified non-reversible. The two farthest lie at distances 0.0035 and 0.0014, with J₁₀ ≈ 10⁻¹¹, i.e. on
  the reversible set to the precision of the descent.
* **Conjecture: t\* = 4.** The degree-8 witness is complete on PU(3).

## Local structure of the reversible set

The reversible set is a union of strata R_C = {U : Uᵀ ∝ CUC†}, one per twisted class C ~ D̄CD⁻¹. There are 4 twisted
classes, of sizes 36, 72, 54 and 54.

| twisted class | C̄C | stratum dim | codim | Hessian rank of J₆ | of J₈ |
|---|---|---|---|---|---|
| 36 (symmetric-unitary type) | 1 | **5** | 3 | **3** = codim ✓ | **3** ✓ |
| 72 | order 3 | 1 | 7 | 1 | 1 |
| 54, 54 | order 4 | 1 | 7 | 0 | 0 |

The certificate field `hessian_rank` uses a relative threshold. Where it reads 8 on a 54-class, every eigenvalue is
below 2×10⁻⁸: the Hessian vanishes there, so the rank is 0.

* On the generic stratum (U symmetric up to Clifford), J₆ has full normal rank: near generic reversible points the
  zero set of J₆ is exactly the stratum.
* For C̄C ≠ 1, Uᵀ ∝ CUC† forces U to commute projectively with C̄C. Those strata are 1-dimensional, have a degenerate
  spectrum, and J₆ is flat to second order there.

## Control (non-vacuous search)

The same pipeline on d = 5 with J₄, which Pass 11357 found incomplete on words:
* **all 400 of 400** descents end at zeros of J₄ that are certified non-reversible (J₆ ≥ 0.50, distance up to 2.0).
* So the search does find spurious zeros when they exist. J₈'s clean result for qutrits is therefore not an artefact of
  a blind search.

**Reading.**
* Pass 11355's "J₆ is complete" is true on the Clifford+T words where it was tested. On all of PU(3) it is false.
* A fixed-degree polynomial witness detects the arrow exactly only from some degree t\* on. For one qutrit
  t\* ∈ [4, 216], and the computations point to t\* = 4.
* The Clifford+T words avoid J₆'s spurious family on every word checked. Why they do is open: the family is
  4-dimensional, the words are countable, and nothing forbids a deep word from landing on it.

## Deep words: still complete, but nearly blind

Random words of depth k, each decided exactly by the Weyl criterion, 200 000 per depth (`--deep`; see the
certificate's `deep_words`):

* Depths 5, 6, 8, 10, 12, 16, 20 (1 399 996 decided words): **no word with J₆ = 0 that violates**, and no reversible
  word with J₆ > 0.
* Four words at depths 16 and 20 had floating-point coefficient magnitudes too close for the exact decider. They are
  reported as numerically ambiguous, not guessed.
* The smallest J₆ of a violator shrinks with depth: 1.5×10⁻⁴ at depth 5, 1.19×10⁻⁷ at depths 6–10, 2.3×10⁻⁸ at depth
  12, 9.1×10⁻⁹ at depth 16.
* **The words come very close to the spurious family.** A depth-6 violator, also found again at depth 8 (the same
  group element), lies at distance 0.217 from the reversible set. There J₆ = 1.19×10⁻⁷, J₈ = 2.0×10⁻⁵ and
  J₁₀ = 1.2×10⁻³. J₆ is positive, far above rounding level (10⁻¹³), but three to four orders of magnitude weaker than
  the higher witnesses.
* **That approach is forced, not surprising.** Clifford+T words are dense in PU(3). So violating words come
  arbitrarily close to any non-reversible point of J₆'s spurious family, and inf J₆ over violating words is **0**.
  There is no positive lower bound on the cubic witness over words. What remains open is only whether some word lies
  **exactly** on the family.
* **Operational reading.** The witness of degree 2t uses t copies of the tick. Detecting the arrow on all of PU(3)
  needs at least 4 copies (t\* ≥ 4). Three copies suffice on every Clifford+T word checked, but can be extremely
  weak.

**Correction recorded.** Pass 11355's "next step" phrasing ("prove J₆ completeness") was an over-read. This pass
shows that statement is false beyond words.
