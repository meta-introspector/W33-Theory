# Pass 11111 — three generations cannot carry the two-qutrit Pauli group, alone or glued to colour

Producer: `analysis/w33_pass11111_no_243_family_group.py`
Certificate: `data/w33_pass11111_no_243_family_group.json`
Regression: `tests/test_w33_pass11111_no_243_family_group.py`

Pass 11105 found exactly one Wilson-line-free torus, and so one qutrit of flavour symmetry, in all 104 models. This pass
shows that the count is forced. It also closes the last route by which W(3,3)'s group of order 243 could act on the
generations.

## A. Representation theory (computed)

The k-qutrit Pauli group 3^(1+2k) is built as explicit 3^k × 3^k clock/shift matrices. Its structure:
* 3^(2k) linear characters;
* exactly **two** irreducible representations on which the centre acts nontrivially, each of dimension **3^k**. The
  defining representation has ⟨χ, χ⟩ = 1, and the degrees satisfy 3^(2k) + 2·9^k = |G|.

The degree statement is textbook (Stone–von Neumann for extraspecial p-groups); what is new is its application.
In the orbifold the centre is the point-group twist, which acts on a θ^l field as ω^l ≠ 1 on twisted matter. So
twisted generations that carry the group faithfully come in multiples of 3^k:

**three generations force k = 1: one qutrit, one free torus.**

Two free tori (k = 2) would give twisted families in 9s. Untwisted families (l = 0) see the group only through its
abelian quotient 3^4, so the symplectic form, and with it the W(3,3) geometry, is invisible on them.

## B. The colour loophole (tested)

Quarks carry generation ⊗ colour: 3 × 3 = 9, the faithful dimension of 3^(1+4). SU(3)_colour contains its own clock/shift
H27, whose centre is Z(SU(3)). The central product (family H27) ∘ (colour H27) would then be the 243 group acting on
quarks. That needs the point-group phase to equal a Standard-Model gauge-centre element on every quark:

    ω^l = ω^{a t} e^{2πi b Y} (−1)^{2 s T3}     (t = colour triality).

In all 12 survivors the light quarks are Q (3, 1/6, l = 2), ū (3̄, −2/3, l = 2) and d̄ (3̄, 1/3, l = 2). **No (a, b, s)
solves the condition**: doubling the Q condition contradicts the d̄ condition mod 1.

The control is non-vacuous: with l = 0 the same search finds the Standard Model's own Z6 centre elements. The orbifold
twist is not a gauge-centre element on quarks, so the family H27 cannot be glued to colour.

## Reading

The W(3,3) group of order 243 (the other track's E8 trinification construction, `analysis/2026-09-21_e8_trinification_two_qutrit_pauli243.md`)
is excluded as a family symmetry on both routes:
* **as a pure family group**, by representation theory;
* **as family ⊗ colour**, by the sector assignment of the quarks.

The flavour symmetry of three string generations is at most one qutrit, Δ(54), exactly as found in Pass 11105.
