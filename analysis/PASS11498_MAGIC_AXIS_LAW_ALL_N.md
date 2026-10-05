# Pass 11498 — the magic-axis law is a theorem for every number of qutrits

Producer: `analysis/w33_pass11498_magic_axis_law_all_n.py`
Certificate: `data/w33_pass11498_magic_axis_law_all_n.json`
Regression: `tests/test_w33_pass11498_11505.py`

## Statement

> **Theorem.** Let U = W(a) V_M T₁ be a one-magic-gate tick on any number n of qutrits, with z₁ the magic axis. If
> M²z₁ = z₁ and v = z₁ + Mz₁ ≠ 0, then every frame a with ω(v, a) ≠ 0 is violating.

This is Pass 11420's conjecture F′. Passes 11433, 11457, 11487 and 11486 proved it piece by piece:
* exhaustively at n = 2;
* class by class under F₃ hypotheses;
* 79.7% of the mass at n = 3;
* all but 0.24% of Stab(z₁) at n = 4.

**This pass proves it uniformly, for every n, in five lines.** The class-by-class theorems become special cases.

## Proof

The variables are those of Pass 11350's linear decider. U is reversible iff, for some k ∈ F₃, a symplectic Q solves

* **(S)** M s^{−k} Q (J M J) = Q and Q z₁ = z₁,
* **(F)** g(e) + s^{−k} r_Q = M^{−1}(e − a) − s^{−k} Q J a, with e_{x₁} = k.

Here J = diag(1, −1, …), s is the unit shear (x₁, z₁) ↦ (x₁, z₁ + x₁), g(e) = e + t_k, and r_Q is the frame of
T₁ V_Q T₁⁻¹ V_Q†.

1. **The magic axis is isotropic to its image.** ω(Mz₁, z₁) = ω(M²z₁, Mz₁) = ω(z₁, Mz₁) = −ω(Mz₁, z₁), so
   ω(Mz₁, z₁) = 0. Since ω(y, z₁) = y_{x₁}:
   * z₁ and Mz₁ have no x₁-component;
   * s fixes v;
   * ω(v, z₁) = 0.
2. **A := QJ reverses v.** A is anti-symplectic, and (S) reads A = M s^{−k} A M, i.e. AM = s^k M⁻¹ A.
   * Jz₁ = −z₁ and Qz₁ = z₁ give Az₁ = −z₁.
   * M⁻¹z₁ = Mz₁, so AMz₁ = s^k M⁻¹ Az₁ = −s^k Mz₁ = −Mz₁.
   * Hence **Av = −v**.
3. **The T-twists are harmless.**
   * Qz₁ = z₁ means V_Q fixes W(z₁). T₁ is diagonal on qutrit 1, hence a function of W(z₁). So V_Q commutes with T₁,
     and r_Q = 0.
   * T X^k T⁻¹ = X^k times a quadratic phase, so t_k ∈ span(z₁).
4. **Pair (F) with ω(v, ·).**
   * M and s^k fix v, so ω(v, M⁻¹y) = ω(v, s^{−k}y) = ω(v, y).
   * ω(v, t_k) = 0 by step 1, and r_Q = 0 by step 3.
   * What remains is ω(v, a) + ω(v, Aa) = 0.
   * A is anti-symplectic with A⁻¹v = −v, so ω(v, Aa) = −ω(A⁻¹v, a) = ω(v, a).
   * Hence 2ω(v, a) = 0, i.e. ω(v, a) = 0. ∎

**What the proof rests on.** It uses only Pass 11350's normal form (S)/(F). That form was validated on all 51 840
two-qutrit classes against the exhaustive census of Pass 11330.

## Verification (each step, and the conclusion)

| check | n = 1 | n = 2 | n = 3 | n = 4 |
|---|---|---|---|---|
| t_k ∈ span(z₁), k = 0, 1, 2 | ✓ | ✓ | ✓ | ✓ |
| r_Q = 0 for random symplectic Q with Qz₁ = z₁ | 30/30 | 30/30 | 30/30 | 4/4 |

**n = 2 (all 1296 F′-setting classes).**
* ω(Mz₁, z₁) = 0 on all of them.
* **QJv = −v for all 23 328 enumerated solutions of (S).** Three classes exceed the enumeration cap.
* On the 1295 classes the decider settles, every reversible frame has ω(v, a) = 0.

**n = 3 (Pass 11373's 207 orbit representatives in the F′ setting).**
* All 16 020 enumerated solutions have QJv = −v.
* On the 172 orbits the decider reaches, every reversible frame has ω(v, a) = 0.

**Step (2) on whole solution spaces, for n = 2, 3, 4 (stage `linear`).**
* The proof never uses the symplecticity of Q, only (S) and Qz₁ = z₁. Those cut out an affine space x₀ + span(bᵢ).
* QJv = −v holds on all of it iff x₀Jv = −v and every bᵢJv = 0, which is a linear check with no enumeration.

| n | classes | affine solution spaces (k = 0, 1, 2) | total dimension | QJv = −v on the whole space |
|---|---|---|---|---|
| 2 | 1296 | 2848 | 8 802 | **2848/2848** |
| 3 | 207 orbits | 471 | 4 154 | **471/471** |
| 4 | all 549 Stab(z₁) classes | 1365 | 19 988 | **1365/1365** |

* This includes the n = 4 classes that were beyond every earlier decider cap.
* The identity is a theorem; this checks the implementation of (S) against it.

**How much of the one-gate landscape the theorem governs.**
* Mz₁ is uniform over the 3^{2n} − 1 nonzero vectors.
* Sp(2n, 3) is transitive on pairs (z₁, y) with a given ω(z₁, y). For M²z₁ = z₁, y = Mz₁ must satisfy ω(z₁, y) = 0, and
  the pair can be swapped.
* So:
  * **Pr(M²z₁ = z₁, v ≠ 0) = 2/(3^{2n} − 1)** (the y = z₁ term plus the isotropic independent y; the y = −z₁ term gives
    v = 0);
  * Pr(Mz₁ = −z₁) = 1/(3^{2n} − 1);
  * Pr(M²z₁ = −z₁) = 3/(3^{2n} − 1).
* **Check at n = 2:** 51 840 · 2/80 = **1296** F′-setting classes and 51 840 · 4/80 = **2592** good-cell classes, exactly
  the counts found.
* **Scope.** The theorem is a structural law of one cell, with share 2/(3^{2n} − 1), which shrinks with n. It does
  **not** explain the bulk of the ≈ 1/8 bad-class fraction, which lives in the collinear and non-collinear cells.

**What this changes.**
* The partial coverage tables of Passes 11486 and 11487 (n = 3: 79.7%; n = 4: 0.24% open) are **superseded**. Their
  theorems remain correct but are no longer needed.
* The magnitude and phase arguments were explaining F′ through the Weyl coefficients. The decider's own variables
  make it a two-line consequence of (S): **every reversing map anti-fixes v**.
* Pass 11433's grading lemma (Π = W(v) commutes with V_M T₁; U has Π-degree ω(v, a)) suggests an operator form of the
  same fact: a reversal would send the degree d to −d, forcing d = 0. That reading is not worked out here; the proof
  above does not use it.

**Prior art.**
* The law is this repository's (Pass 11420, conjectured).
* The normal form is Passes 11331 and 11350. No external source is known.
* A search of the corpus for the identities QJv = −v and "Av = −v" found nothing.
