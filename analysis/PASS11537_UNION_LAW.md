# Pass 11537 — the union law: one-gate time reversal is exactly linear algebra on the fixed vectors orthogonal to the magic axis

Producer: `analysis/w33_pass11537_union_law.py`
Certificate: `data/w33_pass11537_union_law.json`
Regression: `tests/test_w33_pass11531_11538.py`

## Statement

* **The setting.** U = W(a)V_M T₁ on n qutrits. A reversal is a symplectic Q solving (S): M s^{−k} Q (JMJ) = Q with
  Qz₁ = z₁ (Pass 11350), and A_Q = QJ.
* **The space K₀.** Let K₀ = {w : Mw = w, ω(w, z₁) = 0}, the fixed vectors orthogonal to the magic axis.
* **The space N_Q.** For each solution Q, N_Q = (I − A_Q⁻¹)K₀.

> **Union law.** Frame a is reversible **if and only if** some solution Q of (S) has a ⊥ N_Q, i.e.
> ω((I − A_Q⁻¹)w, a) = 0 for every w ∈ K₀. If (S) has no solution, no frame is reversible.

**"Only if" is a theorem for every n.**
* Pair the frame condition (F) of a reversing (Q, k) with ω(w, ·) for w ∈ K₀.
* M fixes w. The shear s^k fixes w, because w_{x₁} = ω(w, z₁) = 0.
* ω(w, t_k) = 0, since t_k ∈ span z₁. And r_Q = 0, since Qz₁ = z₁.
* What remains is ω(w, a) + ω(w, A_Q a) = 0, i.e. ω((I − A_Q⁻¹)w, a) = 0. ∎

**"If" is verified exhaustively.**

| | classes / orbits checked | union law exact |
|---|---|---|
| n = 2 (all classes the decider settles) | 51 838 | **51 838** |
| n = 3 (orbits of Pass 11373, enumerable solution spaces) | 2 117 | **2 117** |

* Every cell is covered at both sizes: fixed axis, reversed axis, the same-line cells, collinear-different-lines and
  non-collinear.
* At n = 3, 109 orbits are beyond the decider and 82 have solution spaces too large to enumerate.
* **The checked orbits carry 99.970% of Sp(6, 3)** (28 296 133 / 28 304 640). The unchecked mass is 0.020% beyond the
  decider and 0.010% beyond enumeration.

## What it unifies

The union law contains every earlier result as a special case:
* **The magic-axis law (Pass 11498).** v = z₁ + Mz₁ ∈ K₀, and every A_Q sends it to −v. So (I − A_Q⁻¹)v = −v, and every
  reversible frame satisfies ω(v, a) = 0.
* **The good-rule criterion (Passes 11499 and 11511).** In the good cells, every fixed vector is automatically ⊥ z₁, so
  K₀ = ker(M − I). The criterion there is exactly the union law, and its sufficiency was proved directly.
* **The one-gate W-law (Pass 11533).** v ∈ W means (I − A_Q⁻¹)v = −v for every Q, so v ∈ N_Q for all Q. Violators
  therefore include {ω(v, a) ≠ 0}.
* **The W-law's exceptions.** These are the 24 transvection classes at n = 2 and 9 orbits at n = 3. Their violating sets
  are affine pieces, intersections of the N_Q^⊥ that no single vector explains. The union law gets them exactly.

## Per-solution exactness, and where the union is needed

* **The test** (stage `perq`): for each single solution (Q, k), are the frames reversible through it exactly N_Q^⊥?
  It was run on every 7th class of Sp(4, 3), with all of their solutions.
* **The generic cells: per-solution exactness.** In the non-collinear, collinear-different-lines and all same-line cells,
  every solution has its frame set **exactly** N_Q^⊥: 13 735 (Q, k) pairs, no exception. There the extra k-dependent
  condition below is automatically implied.
* **The eigenvector cells: the union is needed.** In Mz₁ = z₁ and Mz₁ = −z₁, some solutions have frame sets that are
  proper subsets of N_Q^⊥ (684 and 558 of the pairs). The union over several Q still recovers the full set.
* **Proof targets.**
  * In the generic cells: per-solution sufficiency.
  * In the eigenvector cells: that other solutions cover the gaps. There Pass 11498 (fixed axis) and Pass 11511 (good
    cells) already prove what is needed for the class verdicts.

## What a proof of sufficiency needs

* (F) for (Q, k) is solvable iff its right-hand side pairs correctly with every w in K′ = {w : (I − M)w ∈ span z₁}.
* For symplectic M, K^⊥ = Im(M − I), so exactly one of two things happens:
  1. z₁ ∉ Im(M − I): some fixed vector has ω(w, z₁) ≠ 0. Its pairing picks up the shear and the T-twist.
  2. z₁ ∈ Im(M − I): then K = K₀, but K′ = K + span(w₀) with (I − M)w₀ = z₁. That adds the affine condition
     ω((A⁻¹ − I)w₀ + (1 − kc₀)z₁, a) = τ_k c₀ − k, where c₀ = ω(w₀, z₁).
* So sufficiency amounts to showing that the freedom in (Q, k) always absorbs this one extra condition when the K₀
  conditions hold. The data say it always does. The argument is open.

**Consequence if proved.** One-gate reversibility at every n would reduce to K₀ and the affine solution spaces of (S): the
union of the subspaces N_Q^⊥ over the solutions. That replaces the decider's frame-by-frame search. The union can still
range over many Q when the solution spaces are large.
