# Pass 11152 — the best spooky correlation W(3,3)'s own measurements allow is an exact algebraic number

Producer: `analysis/w33_pass11152_pauli_bell_charpoly.py`
Regression: `tests/test_w33_pass11152_pauli_bell_charpoly.py`

With Pauli measurements only, the best two-qutrit state reaches 2.6016791 on the CGLMP inequality (Pass 11149). The
optimal Pauli strategy is the simplest one: both parties measure in the **Z and X eigenbases**, with the natural
labelling.

**Exact characteristic polynomial.** The spectral projectors of Z and X have entries in ℚ(ω), so each coefficient of the
Bell operator's characteristic polynomial is a rational with denominator dividing 9⁹. A 60-digit Faddeev–LeVerrier
computation therefore identifies the coefficients exactly (deviation below 10⁻⁵⁸):

    char(x) = x⁹ − 12x⁷ + 45x⁵ − 6x⁴ − 57x³ + 18x² + 6x − 2 = (x³ − 6x − 2)(x³ − 3x + 1)².

**What it gives.**
* **The Pauli-only maximum** is the largest root of x³ − 6x − 2 = 0:
  **2√2·cos(⅓·arccos(1/(2√2))) = 2.6016791**, above the local bound 2. It needs the *state* to carry magic: the
  maximising state is not a stabilizer state (mana 0.392).
* **The rest of the spectrum** is 2cos(2πk/9) for k = 1, 2, 4 (the roots of x³ − 3x + 1, each twice): ninth roots of
  unity in a qutrit problem. Why ninth roots appear is **not** explained here and is left open.
