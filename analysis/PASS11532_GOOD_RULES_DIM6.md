# Pass 11532 — the good rules extend to blocks of dimension 6: the four indecomposable cases, each settled by witnesses

Producer: `analysis/w33_pass11532_good_rules_dim6.py`
Certificate: `data/w33_pass11532_good_rules_dim6.json`
Regression: `tests/test_w33_pass11531_11538.py`

**Where Pass 11511 left it.**
* The good rules reduce to two block statements:
  * **(U)** the reversers of the unipotent part on V₁ cover every frame;
  * **(E)** some reverser of M on z₁'s block negates z₁.
* Both were proved for blocks of dimension ≤ 4 by exhaustion. Sp(6, 3) has 3¹⁸ unipotent elements, too many to exhaust.

## Two reductions

1. **Conjugation invariance.** Both statements are invariant under conjugation in Sp, so one representative per class
   suffices.
2. **Orthogonal sums.** If u = u_a ⊥ u_b, block-diagonal reversers A_a ⊕ A_b act blockwise.
   * (U) for u follows from (U) on each block: choose A_a for a's component in V_a and A_b for its component in V_b.
   * Likewise for (E): split z = z_a + z_b.
   * By the standard orthogonal decomposition of unipotent classes in odd characteristic, every unipotent element of
     Sp(2m, 3) is an orthogonal sum of indecomposables. These are a single even Jordan block J_{2k}, or a pair J_k ⊕ J_k
     of odd blocks.
   * In dimension 6 the indecomposables of full dimension are J₆, in two classes told apart by the square class of
     ω(v, N⁵v), and J₃ ⊕ J₃. All others split into blocks of dimension ≤ 4, which Pass 11511 settled.
   * For the x² + 1 type, M = s·u with s² = −1 makes V a 3-dimensional Hermitian space over F₉. The only indecomposable of
     real dimension 6 is the regular unipotent of U(3, 9).

## Witnesses

* **How a witness works.** For a given element, the reversers are A₀·C(u). Here A₀ is one anti-symplectic reverser and
  C(u) the symplectic centraliser, both found by F₃ linear algebra plus a symplecticity filter.
  * Reversers are drawn until their criteria cover all 729 frames: that proves (U).
  * For (E), some drawn reverser must negate each admissible z.
* **Sampling.** Random elements of Sp(6, 3) were drawn until each indecomposable class had appeared 4 times: 172 samples.

| class | (U) | (E) |
|---|---|---|
| J₆, form class 1 | 4/4 proved | 4/4 proved (for −u) |
| J₆, form class 2 | 4/4 proved | 4/4 proved (for −u) |
| J₃ ⊕ J₃ | 4/4 proved | 4/4 proved (for −u) |
| regular unipotent of U(3, 9) (x² + 1 type) | — | 4/4 proved |

> **Theorem.** For every n, the good rules hold whenever every indecomposable piece of M's generalised 1-eigenspace and of
> z₁'s block has dimension ≤ 6. Equivalently, at eigenvalue 1 the even Jordan blocks have size ≤ 6 and the odd ones size
> ≤ 3, and z₁'s block pieces are of the corresponding types.

**What is left for every n.** Indecomposables of dimension ≥ 8:
* J_{2k} for 2k ≥ 8;
* J_k ⊕ J_k for odd k ≥ 5;
* regular unipotents of U(k, 9) for k ≥ 4.

These are single classes, or pairs of classes, in each dimension, so the witness method extends dimension by dimension.
An inductive argument would close the general case.

**Scope.**
* The witnesses are proofs for the classes they represent.
* That the listed indecomposables are all of them rests on the classification of unipotent classes in finite symplectic
  and unitary groups (Springer–Steinberg; Wall).
