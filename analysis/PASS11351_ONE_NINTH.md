# Pass 11351 — the 1/9 is (2/9) × (1/2): a twisted-centralizer count times the parity −I

Producer: `analysis/w33_pass11351_one_ninth.py` (uses the linear decider of Pass 11350)
Certificate: `data/w33_pass11351_one_ninth.json`
Regression: `tests/test_w33_pass11350_11354.py`

On the cell where M sends the magic axis z₁ to a point **not collinear** with it in W(3,3) (34 992 classes), exactly
3888 = 1/9 are bad (Pass 11331). The equation (S) of Pass 11350 explains why.

**Count of solutions of (S) against verdict (exhaustive over the cell).** Only k = 0 has solutions on this cell.

| symplectic solutions of (S) | 1 | 3 | 4 | 6 | 24 |
|---|---|---|---|---|---|
| classes | 23 328 | 7776 | 972 | 2592 | 324 |
| bad | **0** | **3888 (half)** | 0 | 0 | 0 |

* **Only the three-solution classes can be bad.** In all 7776 of them, the three solutions differ by a
  **transvection**: they form a coset of a transvection subgroup.
* **The half is a sign.** M is bad ⇔ **−M is good** (7776/7776). Badness is unchanged under M ↦ M⁻¹ and M ↦ JMJ.
* So **1/9 = (7776/34 992) × ½ = (2/9) × (1/2).**

**The parity −I across all cells.** −I is qutrit charge conjugation, the central involution of Sp(4,3).

| cell | (bad M, bad −M) |
|---|---|
| z₁ fixed | (1,0) × 648 |
| z₁ reversed | (0,1) × 648 |
| non-collinear | (0,0) × 27 216, (1,0) × 3888, (0,1) × 3888 |
| collinear, different lines | (0,0) × 11 664 |
| same line, M²z₁ = z₁ | (1,1) × 648 |
| same line, M²z₁ = −z₁ | (0,0) × 648 |
| same line, other | (0,0) × 1296, (1,1) × 1296 |

* −I exchanges reversible and violating classes exactly on the non-line cells. It preserves the verdict on the
  single-line cells, because (−M)² = M².
* One qutrit shows the same exchange: the unit shears are bad, and the parity-twisted shears (−I times a shear) are
  good (Pass 11266).

**Status.** The 1/8 at n = 2 is now accounted for cell by cell: 648 + 648 + 1296 + 3888 = 6480. Pass 11352 shows the
fraction is **not** 1/8 at n = 3 (0.1333 ± 0.0019).
