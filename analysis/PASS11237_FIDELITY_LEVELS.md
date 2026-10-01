# Pass 11237 — the time-reversal fidelity levels are algebraic: the minimal quantum, its square, and a cubic field

Producer: `analysis/w33_pass11237_fidelity_levels.py`
Certificate: `data/w33_pass11237_fidelity_levels.json`
Regression: `tests/test_w33_pass11236_11239.py`

**Question.** Passes 11228, 11235 and 11238 found a discrete spectrum of best time-reversal fidelities F_T for
Clifford+cubic dynamics. Which numbers are they?

**Method (sound identification).**
* A first attempt ran findpoly on double-precision values with coefficients up to 10⁷. It "identified" every level with
  a degree-2 polynomial, but such fits have more free digits than 15-digit data. They were spurious and are discarded.
* Here the 216 single-qutrit Cliffords are rebuilt **exactly** in mpmath at 120 digits.
* Words realising each level are regenerated from the same random stream, and F_T is recomputed over all 216 reversals.
* A minimal polynomial of F_T² is accepted only if it has small coefficients (|c| ≤ 10⁵, degree ≤ 12) and vanishes to
  10⁻¹⁰⁵.

**Result.**

| level F_T | minimal polynomial of F_T² | field |
|---|---|---|
| 0.84402962874598535680… = (1 + 2cos 2π/9)/3 | 81x³ − 81x² + 18x − 1 | Q(cos 2π/9) (PSLQ relation) |
| 0.81009548548684328962… | 59049x³ − 39366x² + 405x − 1 (59049 = 3¹⁰) | Q(cos 2π/9) (PSLQ relation) |
| 0.93926257709015054796… | not of degree ≤ 12 with \|c\| ≤ 10⁵ | — |
| 0.72578765402643956338… | not of degree ≤ 12 with \|c\| ≤ 10⁵ | — |

* **The minimal quantum is cubic, as its closed form requires**, which also checks the method.
* **A second level, 0.8100955, lies in the same cubic field Q(cos 2π/9).**
* **The level 0.712386014 equals F_min² = 0.71238601420108587…**, the very number whose minimal polynomial is the
  cubic in the first row. This level arises only for two-qutrit words, so it is checked in double precision: the best
  reversal of (I⊗T)·SUM·(I⊗T²) over all 51 840 × 81 anti-unitary Cliffords gives 0.7123860142010856, which differs from
  F_min² by 2×10⁻¹⁶. So this level's F is the minimal level's F², to machine precision. It is the fidelity of every target-sector
  violator and of four both-sector words in Pass 11238.
  Numerically the stronger violation reads as **two 2π/9 quanta composed multiplicatively**. An exact proof, and a
  mechanism (why the best reversal factorises over the two cubic gates), are open.
* The levels 0.939 and 0.7258 have no minimal polynomial of degree ≤ 12 with |c| ≤ 10⁵. Their fields are not
  identified here.
