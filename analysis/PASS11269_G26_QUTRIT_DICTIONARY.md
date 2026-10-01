# Pass 11269 — the G26 Cartan dictionary: mirrors are the qutrit SIC and stabilizer states; the invariants are SIC moments

Producer: `analysis/w33_pass11269_g26_qutrit_dictionary.py` (uses the parallel session's
`w33_pass11255_11259_g26_common.py` coordinates and invariants)
Certificate: `data/w33_pass11269_g26_qutrit_dictionary.json`
Regression: `tests/test_w33_pass11265_11269.py`

**Cited, not re-derived.**
* G26 = C₂ × G25 with G25 = 3^{1+2}.SL(2,3) (`data/w33_extended_clifford_g26_no_go.json`, Magma handbook).
* The collineation image is the Hessian group of order 216, which is the projective one-qutrit Clifford group ASL(2,3)
  (`analysis/2026-09-23_extended_clifford_hesse_null_cone.md`).
* The Hesse configuration's nine flexes form the qutrit SIC (Hughston 2007; Bengtsson et al.).
* The anti-linear extended Clifford group is **not** G26 (the repo's no-go stands).
* The algebraic identities below are very likely classical (Maschke 1889 invariant theory of the Hessian group).
  What is new is their use as a reading of the parallel session's vacuum data.

## Read Codex's Cartan coordinates (x, y, z) as qutrit amplitudes. Then, exactly:

1. **The reflection Jacobian is a product of SIC and MUB overlaps.**
   det ∂(u₆, u₁₂, u₉²)/∂(x, y, z) = **19440 · ∏_{9 SIC} (n_a·v) · ∏_{12 stabilizer} (s·v)²**
   * Proved symbolically in Q[ω]/(ω²+ω+1)[x,y,z]; both sides reduce to the same 42-term polynomial.
   * Degree count: 9·1 + 12·2 = 33 = 5 + 11 + 17.
2. **The 21 mirrors are the qutrit's distinguished states.**
   * The 9 order-2 mirror normals form a **SIC**: every pairwise |overlap|² is 1/4.
   * The 12 order-3 mirror normals are the **4 MUBs**: overlaps 0 or 1/3. These are all qutrit stabilizer states.
3. **The reflections are qutrit Cliffords and generate exactly G26.** They generate a group of order **1296** whose
   collineation image has order **216**, and every element is a qutrit Clifford.
4. **Every basic invariant is a SIC moment** (exact symbolic identities), with n_a the nine SIC normals:
   * u₆ = (1/6) Σ_a (n_a·v)⁶
   * u₁₂ = (1/6) Σ_a (n_a·v)¹²
   * u₉ = ∏_a (n_a·v)

   Stabilizer power sums do **not** reproduce u₆ (least-squares residual 0.54).
5. **Codex's control point (1, 2, 3)** has every SIC and stabilizer overlap nonzero, so it is regular.

## What this gives the G26 vacuum programme (Passes 11255–11263 and 11269)

* **A vacuum is a qutrit state.** Its G26 invariant coordinates are the 6th and 12th moments of its nine SIC amplitudes
  and the square of their product.
* **Degeneration.** A Cartan vacuum is non-regular (reflection ramification, the Hessian zero modes of Passes
  11257/11259) **exactly when its qutrit is orthogonal to a SIC vector or to a stabilizer state**.
* **The two kinds of ramification are physically distinct.**
  * Order-2 mirrors (SIC) are where u₉ vanishes.
  * Order-3 mirrors (stabilizer states) are where the Jacobian vanishes to second order.
* The control state of Pass 11258's CP interface (built on (T⊗T)SUM from Pass 11227) is a regular point.
