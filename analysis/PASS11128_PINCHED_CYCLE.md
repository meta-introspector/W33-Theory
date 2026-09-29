# Pass 11128 — the pinched cycle: an exact linear winding law down to the −½ ground state, with no island and no new vacuum

Producer: `analysis/w33_pass11128_pinched_cycle.py`
Scan: `analysis/w33_pass11128_scan_small_radius.py` (Pass 11107 enumerator, winding ranges scaled as 2.5/y)
Frozen: `data/w33_pass11128_small_radius_six.json`
Regression: `tests/test_w33_pass11128_pinched_cycle.py`

Setup: the neutral Wilson-line torus is shrunk from y = 1.6 to 0.2, with the other at 3i, B = 0, and the family torus at ρ.

**Exact law (all six, to 10⁻¹⁰).** Δ = −½ + p_R²/2 with:
* p_R² = y₁/√3 in 36621 and 17224;
* p_R² = y₁/√3 + 1/(3√3·y₂) in 46043, 24165, 40521 and 5904.

This is a unit winding on the A2 torus, plus, in the golden models, a Wilson-line-shifted momentum ⅓ on the other torus.
Level matching then fixes l² = 1.

**The golden radius, derived.** On the diagonal y₁ = y₂ = y, the onset p_R² = 1 is the root of y² − √3·y + ⅓ = 0:

    y_c = (√3 + √(5/3))/2 = (3 + √5)/(2√3) = φ²/√3.

The golden ratio in Pass 11123 enters only through this quadratic, not through any symmetry.

**No tachyon-free island.**
* The lowest state is SM-neutral all the way down; no charged tachyon appears on this line.
* It deepens monotonically toward Δ = −½.
* The number of tachyonic states grows: 4 for y ≥ 0.5, 6 for 0.25–0.4, 8 at 0.2.
* This is the winding tower turning into the momentum tower of a decompactifying T-dual circle, carrying the ground-state
  tachyon of a tachyonic non-supersymmetric heterotic string. The SO(16)×SO(16) string on a small circle with Wilson lines
  is T-dual to the tachyonic 10D strings; see the map of 9D non-supersymmetric heterotic moduli space by
  Fraiman–Graña–Parra De Freitas–Sethi, arXiv:2307.13745.

Reading: shrinking the cycle does not reach a nearby tachyon-free orbifold. The perturbative picture points to the T-dual
tachyonic string, whose own instability takes over. Together with Pass 11123 (no contact quartic), the endpoint of the
neutral condensation is non-perturbative.

Scope:
* B = 0, with one torus shrunk at a time.
* The gauge group of the T-dual string is not identified. Its tachyon has l² = 1, the vector class shared by the SO(32)
  and SO(16)×E8 tachyonic strings.
