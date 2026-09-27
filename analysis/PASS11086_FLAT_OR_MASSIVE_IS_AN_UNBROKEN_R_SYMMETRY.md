# Pass 11076 — "flat-or-massive" is an unbroken R-symmetry, not a law

Producer: `analysis/w33_pass11076_flat_or_massive_is_an_unbroken_r_symmetry.py`
Certificate: `data/w33_pass11076_flat_or_massive_is_an_unbroken_r_symmetry.json`
Regression: `tests/test_w33_pass11076_flat_or_massive_is_an_unbroken_r_symmetry.py`
Uses: the exact lattice solver `analysis/w33_exact_monomial_orders.py` (Pass 11024).

Since Pass 10968 the programme has recorded a recurring pattern. Parity-preserving vacua on which the
superpotential vanishes (W|_S ≡ 0) keep an exotic vector-like colour triplet massless. The manuscript
states the causal version: the triplet is massless *because* W vanishes on the vacuum. This pass tests
the mechanism exactly.

## The test

On every parity vacuum with a decided exotic spectrum, three things are classified.
* **W.** Either W|_S has a monomial, or it is forbidden. A forbidden W is *lattice-forbidden* when no
  integer exponents of any sign exist, so an exact unbroken symmetry of the vacuum carries W's charge (an
  unbroken discrete R-symmetry). Otherwise the obstruction is *holomorphy*: only non-negativity fails.
* **The exotic d–d̄ mass matrix.** Its holomorphic rank, and its rank ignoring holomorphy (lattice).
* **R-neutral pairs.** These are d·d̄ pairs whose total charge lies in the vacuum lattice. Their mass
  term carries exactly W's charge, so it is forbidden iff W is.

## Result

| vacua | W on the vacuum | exotic triplets |
|---|---|---|
| Z12-I, weakest rules (14) | nonzero in 14 | massive 3 · massless 11 (8 by an unbroken symmetry, 3 by holomorphy) |
| Z12-I, admissible Z₂^R (14) | **zero by an unbroken symmetry in 14** | massless 12 (10 by symmetry, 2 by holomorphy) · **massive 2** |
| Z3×Z3 minimal supports (36) | zero by an unbroken symmetry in 36 | massless by an unbroken symmetry in 36 |
| Z3×Z3 extended vacua (28) | nonzero in 28 | massive in 28 |

## Reading

* **Where both W and the exotic masses vanish, they share a cause.** In 46 of those 48 cases an exact
  unbroken symmetry of the vacuum forbids W: an R-symmetry, since W is charged under it. The same symmetry
  forbids the exotic masses. In the other 2, the masses vanish by holomorphy while W is still forbidden
  by symmetry. This is the familiar field-theory fact that an unbroken R-symmetry forbids the
  superpotential and the masses of R-neutral vector-like pairs together. It is the mechanism behind the
  pattern.
* **Neither implies the other.**
  * With the admissible Z₂^R, **2 Z12-I vacua have W forbidden and the exotics massive**. One of them is
    Z12I_1063, which Pass 10978 excluded for other reasons.
  * With the weakest rules, **11 have W ≠ 0 and the exotics massless**.

  So "W|_S ≡ 0 ⇒ exotic triplet massless" is false as a statement, and so is its converse.
* **Correction to the manuscript.** "An exotic triplet massless *because* the superpotential vanishes"
  becomes: *both are forbidden by the same unbroken R-symmetry of the vacuum*. No exclusion in the
  programme relied on the causal reading; every exotic-mass verdict was computed directly.

## Scope

The classification covers the vacua listed above, under the rules stated. Z6-II's Pass 10968 vacua are
not included: the rule that produced them was withdrawn in Pass 10974. The classification says nothing
about how large a mass is, only whether a mass term is allowed.
