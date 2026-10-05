# Pass 11490 — re-audit of the fixed numerical cuts after the 2.999 artefact: the rest hold, one margin is thin

Producer: `analysis/w33_pass11490_cut_audit.py`
Certificate: `data/w33_pass11490_cut_audit.json`
Regression: `tests/test_w33_pass11486_11492.py`

**Why.** Pass 11459 found that the overlap cut 2.999 failed at depth (79% false reversibles at k = 30).
* Clifford+T words are dense, so violators approach the reversible value continuously.
* A cut that separates at shallow depth can therefore stop separating deep down.
* Every **sampled** decision with a hard cut is exposed to the same mechanism. This pass re-checks them at depths beyond
  those originally used, on a fresh sample of 24 000 uniform words per depth.

## A. The deep-word decider (used by Passes 11312 and 11369)

* **Two deciders compared.**
  * R.decide is Pass 11252's Weyl-coefficient criterion.
  * The overlap decider is Pass 11355's Theorem 1 at the tight cut 3 − 10⁻⁷ (Pass 11459).
  * They are independent constructions, so their agreement is a genuine cross-check.

| k | 12 | 16 | 20 | 24 | 30 |
|---|---|---|---|---|---|
| reversible (overlap) | 521 | 136 | 39 | 10 | 1 |
| deciders agree | 24 000 | 24 000 | 24 000 | 24 000 | 23 998 (2 refused as ambiguous) |
| max deficit, reversible | 2.0e-14 | 2.4e-14 | 2.5e-14 | 2.5e-14 | 3.4e-14 |
| min deficit, violator | 1.2e-3 | 2.4e-4 | 2.6e-4 | 2.2e-4 | 6.5e-4 |

**Verdict.** No disagreement in 119 998 decided words. The overlap gap spans about ten decades at every depth.

## B. Pass 11369's J₆ cut (J₆ < 10⁻⁹ counted as a zero)

| k | 12 | 16 | 20 | 24 | 30 |
|---|---|---|---|---|---|
| max \|J₆\|, reversible | 5.3e-14 | 2.5e-14 | 2.8e-14 | 8.9e-15 | 1.4e-14 |
| min J₆, violator | 2.0e-6 | 1.2e-6 | 6.6e-8 | **7.0e-9** | 1.8e-7 |
| violators with J₆ < 10⁻⁵ | 4 | 9 | 8 | 9 | 9 |

**Verdict.**
* **11369's count stands.** No violator has J₆ < 10⁻⁹ at any depth to 30, and no reversible word has |J₆| > 10⁻¹³.
* **The upper margin is thin.** The smallest violating J₆ is 7 × 10⁻⁹ (k = 24, here) and 9.1 × 10⁻⁹ (k = 16, in
  11369), less than one decade above the cut. The lower side has about five decades.
* The left tail of violating J₆ is fed by density, exactly as for the overlap. A larger sample at depth will eventually
  put a violator below any fixed cut.
* **Recommendation for future runs:** decide zeros of J₆ against the rounding scale of the reversible words (about 10⁻¹³),
  or decide reversibility by the overlap first and report J₆ only on violators. Do not use a fixed 10⁻⁹.
* **Scope.** J₆'s completeness (Pass 11355) is an algebraic statement. A sampled cut can only fail to find
  counterexamples; it was never evidence of more than that.

## C. Pass 11353's degenerate-spacing cut (s < 10⁻⁹)

* Tested on 4000 two-qutrit magic circuits at depths 4–32.
* 28 spacings fall below the cut, all at most 2.4 × 10⁻¹⁵ (exact symmetry degeneracies).
* The next spacing is 0.0116, and none lies between 10⁻¹² and 10⁻⁶.
* **Clean**, with about thirteen decades of separation.

## Not re-audited (exact, finite)

* **Cuts that compare finitely many exact objects are not exposed** to the density mechanism, because there are no deep
  words to approach the cut:
  * 11357 and 11372 (Clifford labelling by |tr| = 1 or 3);
  * 11360 (J on depth ≤ 2 words);
  * 11330 and 11331 (equality of finite Clifford products).
* **The exact enumerations of Passes 11422 and 11436 (k ≤ 7) were already re-checked in Pass 11459.** Their largest
  violating overlap is 2.99861.

**Lesson, kept.** Report both one-sided extremes of every cut at the deepest depth used. If the gap on either side shrinks
with depth, the cut is living on borrowed time.
