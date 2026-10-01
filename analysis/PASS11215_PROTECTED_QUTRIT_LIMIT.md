# Pass 11215 — the mean number of protected qutrits, exactly for every n and in the limit: 0.72553535079905834944…

Producer: `analysis/w33_pass11215_protected_qutrit_limit.py`
Certificate: `data/w33_pass11215_protected_qutrit_limit.json`
Regression: `tests/test_w33_pass11213_11216.py`

**Background.** Pass 11207 (merged from the parallel session) proved A(S) = n − c(S), with the Jordan closed form
c = Σ_{λ=±1}(m₁(λ)/2 + m₂(λ)) + m₁^{F₉}(x²+1). Pass 11208 computed E[c] exactly for n ≤ 5 by summing over conjugacy
classes, observed E[c] → 0.7255… by sampling, and noted that "a closed form would follow from cycle-index methods
(Fulman)". This pass carries that out.

**Method.**
* Fulman's cycle index of Sp(2n,q) factorises over the polynomials φ.
* c is additive over φ = z − 1, z + 1 and z² + 1. So Σₙ uⁿ E_n[c] = (1 − u)⁻¹ Σ_φ F_φ(u)/P_φ(u):
  * P_φ sums 1/|centraliser| over the Jordan types at φ (symplectic signed partitions for ±1, unitary partitions over
    F₉ for z² + 1);
  * F_φ weights each type by its contribution to c.
* The centraliser order is |C| = q^{dim C − dim R}·|R|, with dim C = ½Σλ′ᵢ² + ½Σ_{i odd} mᵢ and R the reductive part
  ∏ Sp(mᵢ) × ∏ O^±(mᵢ).
* A first version omitted the −dim R term. It failed the built-in Steinberg check and gave E₁[c] ≠ 1, and was fixed.

**Validation.**
* Steinberg's counts are reproduced exactly: q^{2n²} unipotent elements in Sp(2n,3) and q^{m(m−1)} in U(m,3), for
  n, m ≤ 6.
* The finite-n values reproduce the parallel session's exact class sums for n = 2, 3 and 4: 133/180, 106927/147420
  and 142080247/195832728.
* The n = 5 value, 764515771313563/1053726755353272 = 0.7255351, matches their 940-class value 0.725535.

**Result.**
* **E_n[c] is now exact for every n.** The fractions to n = 8 are in the certificate.
* **Limit (32 digits):** lim E_n[c] = 0.72553535079905834944908732725662.
  * It splits as 0.26397770268859568932465903977598 from each eigenvalue ±1, plus 0.19757994542186697079976924770466
    from the z² + 1 blocks.
  * The n = 12 and n = 14 values differ by 10⁻³¹.
* So a large Clifford dynamics has an unavoidable arrow of **n − 0.72553535…** trits per tick on average.
* The limit is not identified in closed form. Simple q-series in 3 and −1/3 do not match.
* **Correction (same day).** An earlier version of this note said that evaluating the φ-factors at u = 1 is invalid
  because the unitary factor's coefficients tend to a constant. That was wrong. The coefficients are
  3^{−m}/∏(1 − (−1/3)^j) and 3^{−n}/∏(1 − 3^{−2j}), so both factors converge at u = 1 (radius 3). The limit equals
  Σ_φ F_φ(1)/P_φ(1) exactly. The direct evaluation had only differed by its truncation error (~3^{−12}).
  * Exact identity: P_{Sp}(1) = ∏_{r≥1}(1 − 3^{−(2r−1)})^{−1} = ∏_{r≥1}(1 + 3^{−r}) = 1.5649340185670115…
  * Split of the ±1 contribution: E[m₁/2] = 0.0281360246692…, E[m₂] = 0.2358416780193…; z² + 1 part:
    E[m₁] = 0.1975799454218…, with P_U(1) = 1.3891204763652…
  * No closed form was found among simple q-series in 3 and −1/3. The exponent ½Σλ′ᵢ² couples part sizes, so the
    weighted sums do not factor part by part.
