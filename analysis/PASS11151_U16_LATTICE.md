# Pass 11151 — the zero-radius vector lattice is the U(16) lattice, identified by genus and root system

Producer: `analysis/w33_pass11151_u16_lattice.py`
Scan: `analysis/w33_pass11151_scan_zero_radius_lattice.py`
Frozen: `data/w33_pass11151_zero_radius_lattice_21.json`
Regression: `tests/test_w33_pass11151_u16_lattice.py`

The lattice: L0 = {π + N·W : π ∈ E8⊕E8, π·V0 ∈ ℤ, π·W + N·W²/2 ∈ ℤ}. It is built exactly as the kernel of two linear
characters on ℤ¹⁷, and its Gram matrix, determinant, parity and Smith form are computed in exact arithmetic.

* **Validation.** With W = 0 the same construction gives Γ_v⁽⁰⁾ = D8⊕D8 + (s,s), the O(16)×O(16) vector lattice: rank 16,
  even, determinant 4, discriminant ℤ₂×ℤ₂.
* **Zero radius, 21/21 models.** Rank 16, even, determinant 4, discriminant ℤ₂×ℤ₂, the same invariants.
* **Same genus.** The nonzero discriminant classes take the values {0, 1} mod 2 (checked explicitly in 36621, 40521 and
  2233, for the validation and zero-radius lattices alike). That is the even hyperbolic form u: the other even form v
  takes only the value 1, and the odd forms take half-integers. This is consistent with Milgram (signature 16 ≡ 0 mod
  8). By Nikulin, L0 is in the genus of Γ_v⁽⁰⁾.
* **Which member.** Fraiman–Graña–Parra De Freitas–Sethi (arXiv:2307.13745) show this genus is exhausted by the six vector
  lattices Γ_v⁽ᵖ⁾ of the non-supersymmetric heterotic strings, which differ by root system. L0's root system is A₁₅
  (Pass 11140), so **L0 = Γ_v⁽¹⁾, the U(16) string's lattice.**
* **Wilson-line data.** Every model's data (W² ∈ {10/3, …, 38/3}, W·V0 ∈ {0, 2}) lands on the same lattice.

This makes Pass 11140's identification, which counted roots and tachyons, a lattice identification.

A bug caught on the way: the first version of the scan wrote the Wilson-line offset as quadratic in N. That set is not
closed under addition. The correct offset is linear in N, as in the validated tachyon enumerator.
