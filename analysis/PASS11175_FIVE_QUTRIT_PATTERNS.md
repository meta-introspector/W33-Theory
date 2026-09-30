# Pass 11175 — perfect five-qutrit gates: five possible orientation patterns, a rarity bound, and an exact F₉ census that realises only two

Producer: `analysis/w33_pass11175_five_qutrit_patterns.py`
Scans: `analysis/w33_pass11175_scan_sp10_3.py` (`data/w33_pass11175_sp10_3_sample.json`) and
`analysis/w33_pass11175_scan_five_qutrit.py` (`data/w33_pass11175_five_qutrit_scan.json`)
Regression: `tests/test_w33_pass11175_five_qutrit_patterns.py`

**Pattern theorem.** Each column and row of block determinants of a perfect S ∈ Sp(10,3) sums to 1 mod 3 (Pass 11158),
so it holds five or two orientation reversals (−1's). Up to relabelling there are **exactly five** patterns (exhaustive
over the 11⁵ row choices):

| pattern | labelled matrices |
|---|---|
| 10-cycle | 1440 |
| 4-cycle + 6-cycle | 600 |
| full row and column + permutation | 600 |
| two full rows and columns | 100 |
| all −1 | 1 |

**Criterion.** A gate is perfect iff all 25 blocks and all 100 two-party 4×4 minors are invertible. The 3- and 4-party
minors follow from det S[I,J] = ±det S[Iᶜ,Jᶜ], which is Jacobi's identity for symplectic S.

**Rarity.**
* There are **no perfect gates in 7 000 000 samples** of Sp(10,3), so their fraction is below 4.3·10⁻⁷ (95%). For
  comparison, the fraction is 128/4095 for three qutrits and 4/15 for two.
* One local orbit is far smaller still. The explicit circulant of Pass 11170 has local stabiliser 4, so its orbit
  (24¹⁰/4 ≈ 1.6·10¹³ gates) is a fraction 1.0·10⁻¹³ of Sp(10,3). Sampling cannot see it.

**Exact F₉-linear census.** Unitary 5×5 matrices over F₉ with every entry and 2×2 minor nonzero:
* Method: orbit–stabiliser over the column monomial unitaries, with exact controls (2048 row sets for 3×3, 0 for 4×4).
* **2 642 411 520 = 2²³·3²·5·7** of them, a fraction 0.00256 of U(5,F₉). The corresponding fractions for n = 2, 3, 4 are
  2/3, 32/63 and 0.
* They realise **only two of the five patterns**: full row and column + permutation (15 728 640 row sets) and the
  10-cycle (6 291 456).

**Graph states.** Hill climbing for AME(10,3) weighted graph states found none in 6 × 25 minutes. They exist
(Danielsen, arXiv:1106.2428) but are not reachable this way.

**Open.** Whether the 4-cycle + 6-cycle, double-cross and all-(−1) patterns occur for five qutrits.
