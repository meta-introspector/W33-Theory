# E6 cubic tensor: a concrete ±1 sign structure from E8

This repo already verifies that the **support** of the E6-invariant cubic on the 27 is exactly the set of 45 **tritangent planes** (triangles in the meet graph) of the 27 lines on a smooth cubic surface.

The next question is: can we recover a **consistent sign / phase** assignment, not just the set of nonzero monomials?

## What’s implemented here

- `tools/e8_lattice_cocycle.py`
  - Implements a deterministic bimultiplicative cocycle `ε(α,β) ∈ {±1}` on the E8 **root lattice** (using the repo’s Bourbaki simple roots as a Z-basis).
  - Intended as a concrete, reproducible sign convention for experiments with Chevalley-type structure constants.

- `tools/solve_e6_cubic_sign_gauge.py`
  - Uses `ε(α,β)` to assign “raw” signs to ordered mixed triples `(α,β,γ)` in the `(27,3)` sector with `α+β+γ=0`.
  - Solves the GF(2) gauge problem for phases `σ_a(i)∈{±1}` on the three SU(3) colors and symmetric coefficients `d_{ijk}∈{±1}` on the 45 tritangent planes:
    - forces consistency across all 270 ordered triples
    - determines a concrete sign assignment up to gauge.
  - Writes `artifacts/e6_cubic_sign_gauge_solution.json`.

## Key computed outcome

With this cocycle convention:

- A consistent sign gauge **does** exist (the linear system is consistent).
- There is **no gauge** in which all 45 coefficients are `+1` (the “all-plus” system is inconsistent).
- The induced sign split in the current saved solution is `22` positive and `23` negative tritangent-plane
  coefficients. Multiplying the entire cubic by `-1` gives the equivalent `23` positive and `22` negative convention;
  the invariant statement is the `23/22` imbalance and the absence of an all-plus gauge.

This is a precise, testable “phase obstruction” result: **support-only** data is not enough; the cubic invariant carries genuine sign/phase structure.

## Important follow-up: Weyl equivariance happens at the *root-structure-constant* level

When you include the correct Weyl-action signs on root vectors (the `μ` factors from the Chevalley `n_α(1)` elements),
the mixed structure constants satisfy the exact equivariance identity:

`μ(α) μ(β) N(wα,wβ) = N(α,β) μ(w(-γ))` for `α+β+γ=0`.

This is verified computationally for all 6 E6 simple reflections across all 270 mixed triples:
- tool: `tools/verify_mixed_structure_constant_equivariance.py`
- artifact: `artifacts/mixed_structure_constant_equivariance.json`

The remaining “hard part” is to choose a **canonical gauge** that makes the induced 45-term cubic sign tensor
manifestly W(E6)-equivariant in a *single fixed 27-weight basis* (this involves diagonal torus factors and is subtler
than a pure signed-permutation model).

## New: canonical SU(3) color gauge (across all 27 weights)

We can now fix the SU(3) “color” phases **canonically** (and compatibly with the E6 cubic sign tensor) using the A2
singleton roots as ladder operators:

- `tools/solve_canonical_su3_gauge_and_cubic.py`
  - Solves one combined GF(2) system for:
    - sign phases on **all six** mixed 27-orbits (the `3` and the `3̄` triangles),
    - the 45-term symmetric E6 cubic tensor signs `d_{ijk}`,
    - SU(3) antisymmetry bits `su3_eps(oa,ob)` for ordered orbit-pairs in the `3` triangle,
    - uniform ladder-constant bits per singleton-root orbit transition (does **not** force them all to be `+1`).
  - Output: `artifacts/canonical_su3_gauge_and_cubic.json`.

- `tools/verify_canonical_su3_gauge_and_cubic.py`
  - Re-checks every equation instance written by the solver.
  - Output: `artifacts/verify_canonical_su3_gauge_and_cubic.json`.

This is the first truly “canonical” piece of the E8→E6×SU3 factorization in this repo: it identifies the three
mixed 27-orbits as **27 copies of the SU(3) fundamental** with a globally consistent phase convention.

## Where the 45 triples come from (references)

- 27 lines on a cubic surface and their classical subconfigurations:
  - https://en.wikipedia.org/wiki/27_lines_on_a_cubic_surface
  - https://en.wikipedia.org/wiki/Schl%C3%A4fli_graph
  - https://en.wikipedia.org/wiki/Double_six
- A readable notes-style reference that explicitly enumerates the 45 “tritangent trios” and states W(E6) transitivity:
  - https://dept.math.lsa.umich.edu/~idolga/CAG-2.pdf
- A convenient modern reference explicitly connecting tritangent planes and the Schläfli fan:
  - https://link.springer.com/article/10.1007/s00454-020-00215-x

## A classical reference on signs (E6 triads)

Vavilov discusses the 45 “triads” (tritangent planes) and an explicit sign normalization of the cubic form:
- https://iopscience.iop.org/article/10.1070/SM2004v195n09ABEH000825/pdf

One accessible modern-algebra reference point for the same invariant is the Albert algebra / Freudenthal viewpoint:
- https://arxiv.org/abs/1201.0291

## Why a cocycle shows up at all (references)

If you build a Lie algebra (or a lattice VOA) from an even lattice, you must choose a ±1-valued cocycle to fix signs consistently.
This is a standard ingredient in the “lattice construction” literature:

- https://en.wikipedia.org/wiki/Vertex_operator_algebra
- https://en.wikipedia.org/wiki/Lattice_vertex_operator_algebra

## How to reproduce

- `python3 tools/solve_e6_cubic_sign_gauge.py`
- `python3 tools/solve_canonical_su3_gauge_and_cubic.py`
- `python3 tools/verify_canonical_su3_gauge_and_cubic.py`
- `./.venv_test/bin/python -m pytest -q tests/test_e8_lattice_cocycle.py tests/test_e6_cubic_sign_gauge.py`
