# Pass 10978 — the E6 lattice fixes the admissible R-rule of Z12-I; the last 9 models are decided

Producer: `analysis/w33_pass10978_z12_e6_admissible_r_rule.py`
Certificate: `data/w33_pass10978_z12_e6_admissible_r_rule.json`
Regression: `tests/test_w33_pass10978_z12_e6_admissible_r_rule.py`

Pass 10974 reclassified 9 Z12-I verdicts as undetermined: no R-rule is established for Z12-I on the
non-factorizable E6 lattice. This pass settles them by asking which R-rules the lattice can support at
all.

## A. Which symmetries commute with the twist

The Z12-I twist on the E6 root lattice is the Coxeter element c (order 12, eigenphases
1, 4, 5, 7, 8, 11 in units of 2π/12). The producer generates W(E6), all 51,840 elements, exactly as
integer matrices:

* the centraliser of c in W(E6) is **⟨c⟩**, of order 12;
* **−1 ∉ W(E6)**, so Aut(E6) = W × {±1} and the full centraliser is ⟨c⟩ × ⟨−1⟩;
* −c⁶ has eigenvalues (−1, −1, 1, 1, 1, 1): a π-rotation of the **order-3 plane only**.

An R-symmetry of the worldsheet instantons comes from a lattice symmetry commuting with the twist
(Bizet et al., arXiv:1301.2322). So the only R-symmetry beyond the orbifold projection is
**a Z2^R on the order-3 plane**, with Σ R³ ≡ −1 mod 2. The per-plane Z12 × Z12 × Z3 rules used
tentatively in Pass 10968 are not supported by the lattice.

## B. Bracketing every verdict

Each of the 14 models with a parity-viable FI-cancelling vacuum is evaluated twice: with the weakest
rules (gauge + point group) and with the strongest admissible ones (+ Z2^R).

* **Exotic d-triplets.** The all-orders lattice rank over-allows couplings, so a rank deficiency is
  robust.
* **Higgs doublets.** Explicit monomials through order 9 give a lower bound on the rank, so "no light
  pair" is robust. A lattice rank below what is needed means extra massless doublets at all orders.

A model is dead only if it fails in a way that holds under every admissible rule.

(see the certificate for the per-model table; summary in the docs card)

## C. The one candidate, Z12I_1063_303

Under Z2^R its 5-direction vacuum does well at first:

* the exotics are full rank (4/4);
* the Higgs matrix has **structural rank 3 of 4 at all orders**, so exactly one Higgs pair is
  protected;
* W vanishes identically on the vacuum.

But an outside singlet n₃₆ appears linearly at **order 4** (F ~ ε³M_s², supersymmetry broken near the
string scale). Further linear singlets appear at orders 7, 7, 7 and 11. Absorbing **any** of them into
the vacuum lifts the protected pair (Higgs rank 4/4). The symmetry that keeps the Higgs light is exactly
the one the order-4 F-term breaks: flat or light, never both.

## D. Audit of Pass 10967

Pass 10967 stated that Z6II_23 has no u^c d^c d^c at orders 3 and 4. That statement holds under
orbifolder's rules, under the established prime-plane rules alone, and with the γ-corrected G2-plane
charge: 0 and 0 in every case. Gauge and space-group rules alone forbid them, so the claim never
depended on the incomplete R-rule.

## Scope

* The Z2^R is taken without γ-type corrections for its action on fixed points.
* Bizet et al.'s additional instanton "Rule 6" is not implemented. It can only forbid more couplings,
  which could revive a model killed by the μ-problem but not one killed by massless exotics.
