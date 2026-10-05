# Pass 11456 — the magic-axis law at four qutrits: what is settled, and the same-line cell

Producer: `analysis/w33_pass11456_n4_closure.py`
Certificate: `data/w33_pass11456_n4_closure.json`
Regression: `tests/test_w33_pass11456_11460.py`

## Fixed-axis cell = Stab(z₁) ⊂ Sp(8,3), all 549 conjugacy classes

| status | classes | share of H (by class size) |
|---|---|---|
| decided by the sparse decider (Pass 11433): F′ holds, 0 failures | 232 | 97.60% |
| beyond the decider, but **proved** by the magnitude theorem (Pass 11457) | 44 | 0.61% |
| open | 273 | 1.79% |

**Reading.**
* At four qutrits, F′ is established on every fixed-axis class except the "open" ones. Those are small,
  high-symmetry classes with z₁ ∈ Im(M − I) or a degenerate Im(M − I).
* Their status would follow from Pass 11457's extension route (q ≠ 0), or from phase-based arguments for the
  magnitude-blind classes.

## Same-line cell (Mz₁ on a line through z₁, M²z₁ = z₁): 200 targeted classes

* **Construction.** An exact n = 2 class of the cell, embedded as A₂ ⊕ B with B ∈ Sp(4,3), then conjugated by a random
  element of Stab(z₁).

| | decider: F′ with equality | decider: undecided |
|---|---|---|
| theorem proves F′ fully | 45 | **2** |
| theorem proves part of the frames | 19 | 6 |
| theorem silent (z₁ ∈ Im(M − I)) | 102 | 26 |

* **Decided classes.** All **166** decided classes satisfy F′ with exact equality (violating frames =
  {ω(a, z₁ + Mz₁) ≠ 0}). None fails.
* **The theorem adds two.** It proves F′ for 2 classes the decider cannot reach.
* **Still open:** 32 of the 200 sampled classes.
