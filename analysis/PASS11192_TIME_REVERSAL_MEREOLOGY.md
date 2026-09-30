# Pass 11192 — time reversal on the substrate: the E₆ reflections are the Kramers time reversals, and they never respect a qutrit

Producer: `analysis/w33_pass11192_time_reversal_mereology.py`
GAP: `analysis/gap/w33_pass11192_antisymplectic_classes.g` → `data/w33_pass11192_gap_outer_classes.txt`
Frozen: `data/w33_pass11192_time_reversal_mereology.json`
Regression: `tests/test_w33_pass11192_time_reversal_mereology.py`

**Setup.** A time reversal of n qutrits is an antiunitary Clifford operation, taken mod Paulis and phases. It acts on
phase space by an **anti-symplectic** map θ: ω(θx, θy) = −ω(x, y). Complex conjugation is τ = (x, z) ↦ (x, −z).
Projectively, time reversals form the outer coset of PGSp(2n, 3) = PSp(2n, 3).2. For two qutrits this group is W(E₆),
and its odd (det −1) half is the set of time reversals. A time reversal is *local* in a split F if it maps every plane of
F onto a plane of F: a product of single-qutrit antiunitaries, possibly composed with a relabelling of the qutrits.

**Theorems (every n).**
1. **θ² = +1 ⇒ θ is a local complex conjugation in some split.** Its eigenspaces E₊ and E₋ are transverse Lagrangians,
   since ω(x, y) = −ω(x, y) on each. Dual bases eᵢ ∈ E₊, fᵢ ∈ E₋ with ω(eᵢ, fⱼ) = δᵢⱼ give θ-invariant planes ⟨eᵢ, fᵢ⟩,
   on each of which θ = diag(1, −1) = τ. So θ²=+1 time reversals correspond to ordered pairs of transverse Lagrangians:
   1080 matrices for n = 2 and 816 480 for n = 3 (540 and 408 240 projectively).
2. **θ² = −1 (Kramers type: T² is the parity operator) ⇒ n is even.** Let h(x, y) = ω(θx, y) + i·ω(x, y), with i acting
   as θ. This is an F₉-bilinear, alternating, nondegenerate form on V viewed as an F₉-space. Moreover, a single qutrit
   has no anti-symplectic map with square −1: a 2×2 map of determinant −1 squares to a scalar only if its trace
   vanishes, and then it squares to +1 (checked: 12 of the 24 maps square to +1, none to −1). So a Kramers time
   reversal is local only in splits where it **pairs** the qutrits, i.e. permutes them by a fixed-point-free involution.
   For odd n it is local in no split at all.

**Two qutrits (exact: all 51 840 anti-symplectic matrices, 10 outer classes of W(E₆)).**

| class | θ² | size (projective) | local in | as |
|---|---|---|---|---|
| **the 36 reflections of E₆** | −1 | 36 | 15 of 45 splits | all 15 swap the qutrits |
| products of 3 orthogonal reflections | +1 | 540 | 7 splits | 6 local complex conjugations + 1 swap |

* **The reflections of E₆ are exactly the Kramers time reversals.** None of them respects any qutrit. Each is local only
  as an antiunitary SWAP, and it is so in exactly 15 of the 45 splits (tritangent planes).
* **19/45 of all time reversals are local in no split** (10 944 projective, orders 6, 10, 12). This is the same fraction
  as the intrinsically entangling ticks (19/45, Pass 11180). Equivalently, W(E₆) acting on the 45 tritangent planes has
  exactly 10 944 derangements in each coset of W(E₆)⁺. We record this coincidence without an explanation.

**Three qutrits (26 outer classes of PGSp(6, 3)).**
* **No Kramers class exists**, as Theorem 2 requires.
* The single involution class (408 240 elements) is a local complex conjugation in 234 splits and a conjugation with a
  transposition in 117 more.
* **Time reversals are more often intrinsically non-local than ticks: 8377/12285** are local in no split, against
  7922/12285 for ticks (Pass 11181). They have orders 12, 20, 24, 26, 28 and 40.

**Reading.** Time reversal is not one operation on the substrate. Those with θ² = +1 are complex conjugation, relative to
suitable subsystems. The E₆ reflections are a genuinely two-qutrit, Kramers-type time reversal that exists only for an
even number of qutrits and always exchanges them. Pass 11189 showed that τ is what swaps the enantiomeric split relations.

**Prior art.** The count 15 is classical. The reflection s_α of a double-six D swaps aᵢ ↔ bᵢ and fixes every c_ij,
so it fixes exactly the 15 tritangent planes disjoint from D (the all-c_ij syntheme slice, recorded in
`analysis/PASS7217_7232_double_six_doily_spread_code.md`). Under the splits = tritangent planes dictionary (Pass 11177),
the 15 splits in which a Kramers time reversal is local are therefore that slice. What is new here is the reading: θ² = −1,
the parity lemma, and that locality always comes with a swap.

**Scope.**
* The two-qutrit statements are exhaustive. The three-qutrit statements are exact over GAP's complete list of outer
  classes; class sizes sum to |PSp(6, 3)|.
* The class-size identification with W(E₆) is by counting: the involutions of the odd coset of W(E₆) are exactly the
  classes of sizes 36 and 540.
