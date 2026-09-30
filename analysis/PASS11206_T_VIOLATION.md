# Pass 11206 — T-violation on the substrate: none, but half the dynamics can only be reversed anti-unitarily

Producer: `analysis/w33_pass11206_t_violation.py` (GAP: `analysis/gap/w33_pass11206_time_reversal_invertible.g`)
Certificate: `data/w33_pass11206_t_violation.json` (GAP output `data/w33_pass11206_gap_time_reversal.txt`)
Regression: `tests/test_w33_pass11206_t_violation.py`

**Question.** CP violation in the Standard Model is T-violation. Is the substrate's Clifford dynamics time-reversal
symmetric: for each tick S, is there a time reversal T (anti-unitary, i.e. anti-symplectic) with T S T⁻¹ = S⁻¹?

**Result (GAP, every class; independent Python check for two qutrits).**
* **Every tick is T-invertible**, for two and three qutrits. There is no T-violation at the Clifford level. This is an
  instance of **Wonenburger's theorem** (M. J. Wonenburger, 1966): every symplectic map is a product of two
  skew-symplectic involutions. Prior art, re-derived here on the substrate.
* **Real ticks**, those inverted by some *unitary*, make up 85/162 of two-qutrit and 73451/157464 of three-qutrit
  ticks.
* So **77/162 and 84013/157464 of the dynamics can be run backwards only anti-unitarily**. For these ticks Wigner's
  anti-unitarity of time reversal is forced by the dynamics, not by convention.
* All of these ticks have order divisible by 3: orders 3, 6, 9, 12 for two qutrits, and 3 to 36 for three qutrits.
  They are exactly where the qutrit phase e^{2πi/3} enters.

**Reading.** Any T-violation (and so CP violation) in a W(3,3) theory cannot come from its Clifford dynamics. It must
come from the non-Clifford resource, e.g. the signed E6 cubic (the magic layer). This is stated as an open question,
not answered here.
