# Pass 11457 — a proof of the magic-axis law whenever Im(M − I) is nondegenerate, from Weyl-coefficient magnitudes

Producer: `analysis/w33_pass11457_magic_axis_theorem.py`
Certificate: `data/w33_pass11457_magic_axis_theorem.json`
Regression: `tests/test_w33_pass11456_11460.py`

**The law F′ (Pass 11420).** Let U = W(a)V_M T₁ and v = z₁ + Mz₁. If M²z₁ = z₁, every frame with ω(a, v) ≠ 0 violates.
* M²z₁ = z₁ forces Mz₁ to be collinear with z₁, since ω(z₁, Mz₁) = ω(Mz₁, M²z₁) = −ω(z₁, Mz₁).
* So the hypothesis covers exactly the fixed-axis cell and the same-line cell with M²z₁ = z₁.

## Lemmas

* **L1 (support).** The Weyl coefficients of V_M are supported exactly on R_M = Im(M − I), with constant modulus
  |R_M|^{−1/2}. This is standard; checked on 400 two-qutrit and 60 three-qutrit classes.
* **L2 (magnitudes are a reversal invariant).**
  * C = W(b)V_S rescales coefficients by phases and relabels p ↦ Sp. Transposition relabels p ↦ −Jp.
  * So if CUC† ∝ Uᵀ, the modulus m_U = |c_U| is invariant under the **anti-symplectic linear** map L = −JS.
* **L3 (three levels).** T₁ = Σ_j τ_j W(jz₁) with τ_j = (1 + ζ^{1−3j} + ζ^{8−6j})/3. The moduli |τ_j| = 0.844, 0.449
  and 0.293 are **distinct**.
  * Hence c_U(p) = phase · Σ_j τ_j c_{V_M}(p − a − jz₁).
  * If z₁ ∉ R_M, m_U takes the three distinct values |τ_j|·|R_M|^{−1/2} exactly on the cosets R_M + jz₁ + a.

## Theorem (proved)

> Let M²z₁ = z₁, v = z₁ + Mz₁ ≠ 0, z₁ ∉ R_M, and Rad = R_M ∩ R_M^⊥. Then every frame a with ω(a, v) ≠ 0 **and**
> ω(a, Rad) = 0 violates. In particular, **if Im(M − I) is nondegenerate, F′ holds for M.**
> (Nondegeneracy forces z₁ ∉ R_M: otherwise v ∈ R_M ∩ R_M^⊥ = 0.)

**Proof.**
1. Suppose U is reversible, and take L from L2. L preserves the values of m_U, so it fixes each level coset
   R_M + jz₁ + a.
2. Hence L(R_M) = R_M, La − a ∈ R_M and Lz₁ − z₁ ∈ R_M.
3. L is anti-symplectic and preserves R_M, so it preserves R_M^⊥ ∋ v, because Mv = v and ker(M − I) = R_M^⊥.
4. Also Lv ∈ 2Lz₁ + R_M = v + R_M. So r = Lv − v lies in Rad.
5. Then −ω(a, v) = ω(La, Lv) = ω(a, v) + ω(a, r). The last term vanishes for a ⊥ Rad.
6. So 2ω(a, v) ≡ 0, i.e. ω(a, v) ≡ 0 (mod 3), a contradiction. ∎

## Soundness and coverage

| | classes in the two cells | hypotheses hold | F′ fully proved (Rad = 0) | partly proved (Rad ≠ 0) | z₁ ∈ R_M (argument silent) |
|---|---|---|---|---|---|
| n = 2 (all classes) | 1296 | 376 | **352** | 24 | 920 |
| n = 3 (all orbits) | 207 | 58 | **38** orbits, **21.2%** of the cells | 20 | 149 orbits (75.2%) |
| n = 4 (fixed-axis cell only: all 549 classes of Stab(z₁)) | 549 | 141 | **72** classes, **21.4%** of H | 69 | 408 (75.0%) |

* **Soundness at n = 2:** at every one of the 376 classes where the theorem applies, all frames it claims are violating
  are violating in the exhaustive census.
* **Pass 11433's unreachable four-qutrit classes:** of the 317 classes the decider could not reach, the theorem
  **fully settles 44** (25% of their mass). The unresolved fixed-axis mass at n = 4 falls from 2.40% to about 1.79% of H.

**Reading.**
* The magic-axis law is now a theorem on the "nondegenerate" quarter of each bad cell, at every number of qutrits.
* The proof uses only the magnitudes of the Weyl coefficients and the three distinct moduli of T's Weyl expansion.
* The remaining three quarters have z₁ ∈ Im(M − I): the three levels interfere on a single coset and the magnitude
  argument is silent. There F′ is verified (n = 2, 3 exhaustively; n = 4 on 97.6% of the fixed cell) but not proved.

**Extension route (open).**
* When z₁ ∈ R_M, the phases of c_{V_M} are quadratic. On R_M + a the modulus becomes |Σ_j τ_j ω^{js + qj²}|, where
  s = ω(u, p) for a linear functional and q is the quadratic phase of V_M along z₁.
* Computed:
  * q = 0 gives the three moduli 1, 1, 1. The magnitudes are uniform and the argument is genuinely silent.
  * q = ±1 gives **0.778, 0.508, 1.462**, three distinct values.
* So the silent case splits. For q ≠ 0 the level sets are the hyperplanes {ω(u, p) = s} of R_M + a, which an
  anti-symplectic L must preserve: the same kind of constraint that closes the main proof.
* Completing that case, and counting how many silent classes have q ≠ 0, is left open.

## How far can magnitudes go? (`--magnitude`, n = 2)

* **The test.** For U = W(a)G we have m_U(p) = m_G(p − a). An anti-symplectic L preserves m_U iff the affine map
  q ↦ Lq + (La − a) preserves m_G.
* **The computation.** For each class, all affine magnitude automorphisms of m_G were found by exhaustive search over
  the 51 840 anti-symplectic maps and 81 translations. A frame F′ calls violating is **certified by magnitudes alone**
  when no automorphism matches it.
* **Sample:** 60 random classes of the two bad cells.

| class type | all F′ frames certified | some | none |
|---|---|---|---|
| Rad = 0 (the theorem) | **16/16** | — | — |
| z₁ ∈ R_M (silent for the theorem) | **33/42** | 2 | **7** |
| Rad ≠ 0, z₁ ∉ R_M | — | 2 | — |

* **82% of the sampled classes are fully provable from magnitudes.** That includes most of the classes where the
  theorem as stated is silent, as the q ≠ 0 extension predicts.
* About **12% are magnitude-blind**, presumably the uniform q = 0 case. A proof there must use the phases of the Weyl
  coefficients, not only their moduli.
