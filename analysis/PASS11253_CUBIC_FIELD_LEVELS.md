# Pass 11253 — every Clifford+T time-reversal fidelity has F² in the cubic field Q(cos 2π/9)

Producer: `analysis/w33_pass11253_cubic_field_levels.py` (exact Z[ζ₉] arithmetic from Pass 11251)
Certificate: `data/w33_pass11253_cubic_field_levels.json`
Regression: `tests/test_w33_pass11250_11254.py`

## Theorem (any number n of qutrits)

Let U be any word in qutrit Cliffords and T = diag(ζ^{x³}), with ζ = e^{2πi/9}, and let V be any Clifford.
* Up to global phases, which cancel in V U* V†, all entries lie in Q(ζ) times a real power of √3.
* The powers of √3 pair up into rational powers of 3, so t = tr(V U* V† U) ∈ Q(ζ).
* |t|² is fixed by complex conjugation, so it lies in Q(ζ)⁺ = Q(cos 2π/9), which is **cubic**.
* The best reversal fidelity is attained, so **F_T² ∈ Q(cos 2π/9)**: rational or a cubic irrationality, with
  denominators powers of 3.

**Correction to Pass 11237.** The levels 0.939 and 0.7258 were "not identified at degree ≤ 12, |c| ≤ 10⁵". They are
cubic. Their minimal polynomials simply have coefficients up to 3¹⁰ and beyond, outside that search bound. The failure
was the bound, not the degree.

## Computation

* The random word streams of Pass 11235 (two qutrits) and Pass 11237 (one qutrit) are replayed with the same seeds.
* For every violating word, U and an optimal reversal are rebuilt exactly over Z[ζ], and F² is written as
  a + b·c + e·c² with c = 2cos(2π/9) (c³ = 3c − 1).
* **148 distinct levels**: 103 one-qutrit and 45 two-qutrit. All lie in Q(c): 147 are cubic and 1 is rational
  (F = 2/3). Float and exact values agree to 5.9×10⁻¹⁵.

Selected exact levels:

| F | F² minimal polynomial | F itself |
|---|---|---|
| 0.9392625771 | x³ − (40/27)x² + (1264/2187)x − 23104/531441 | (−4 + 2c + 4c²)/9 |
| 0.8440296287 = F_min | x³ − x² + (2/9)x − 1/81 | (1 + c)/3 |
| 0.8100954855 | x³ − (2/3)x² + (5/729)x − 1/59049 | (−2 + 3c + 2c²)/9 |
| 0.7257876540 | x³ − x² + (638/2187)x − 11881/531441 | **(5 + c)/9** |
| 0.7123860142 | x³ − (5/9)x² + (2/81)x − 1/6561 | (1 + 2c + c²)/9 = F_min² |
| 0.6666666667 | x − 4/9 | **2/3** |

**F itself.**
* For 143 of the 148 levels, F itself lies in Q(c).
* For **5** levels it does not (0.849189, 0.841168, 0.836103 on one qutrit; 0.611065, 0.594859 on two). This was
  checked by an exact search over all square roots. The search is rigorous because Z[c] is the full ring of
  integers: disc(c³ − 3c + 1) = 81 equals the field discriminant.
* So the theorem's statement, about F² and not F, is sharp.
