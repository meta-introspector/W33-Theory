# Pass 11217: the arrow law A(S) = n − c(S) for qubits (characteristic 2)

Producer: `analysis/w33_pass11217_qubit_arrow_law.py`
Certificate: `data/w33_pass11217_qubit_arrow_law.json`
Regression: `tests/test_w33_pass11217_qubit_arrow_law.py`

**Background.**
* Pass 11207 proved the intrinsic-arrow law **A(S) = n − c(S)** for every odd local dimension.
  * A(S) is the forgotten record: 2n minus the largest total overlap Σ dim(P ∩ SP) over splits into n planes.
  * c(S) is the largest number of mutually orthogonal S-invariant nondegenerate planes.
* Pass 11208 checked the qubit case, Sp(2n,2), class by class for n ≤ 5, but left it as a conjecture.
* The odd-characteristic proof fails for qubits. It diagonalises b(x,y) = ω(x,Sy) on a Lagrangian L that is also
  Lagrangian for ω(·,(S+S⁻¹)·) (a *bi-Lagrangian*). In characteristic 2, b can be alternating, and an alternating form
  has no orthogonal basis.
* This pass closes the gap: **the law holds for every number of qubits.**

## 1. A criterion that is exactly right in characteristic 2

Let q(x) = ω(x,Sx). In characteristic 2 its polar form is ω(x,(S+S⁻¹)y), so **q is additive (F₂-linear) on every
bi-Lagrangian**.

**Lemma (criterion).** V has a split into n planes that each meet their image if and only if some bi-Lagrangian L
satisfies both of the following:
* (a) L ∩ SL is spanned by eigenvectors of S (over F₂, by fixed vectors);
* (b) q|L ≠ 0, or L = L ∩ SL.

*Proof of ⇐.*
* The radical of b on L is L ∩ (SL)^⊥ = L ∩ SL = R, and q vanishes on R.
* Take a complement L′ of R in L. By (b), b on L′ is nondegenerate and non-alternating, or L′ = 0.
* Over a perfect field of characteristic 2, such a form has an orthogonal basis l_i.
  * The producer builds this basis explicitly. When an alternating remainder appears, a hyperbolic pair (u, v) and an
    earlier x₀ are traded for x₀+u, x₀+v and x₀+u+v.
* The planes span(l_i, Sl_i) are nondegenerate, mutually orthogonal and each meets its image.
* Each eigenvector e_j spanning R is paired with a dual vector x_j in the ω-complement, which contains R as a Lagrangian.
  Each plane span(e_j, x_j) then contains Se_j.

*Proof of ⇒.*
* In each plane pick l_i ∈ P_i with Sl_i ∈ P_i. Then L = span(l_i) is Lagrangian.
* ω(l_i, Sl_j) = 0 for i ≠ j, and ω(l_i, Sl_i) + ω(Sl_i, l_i) = 0 in characteristic 2. So L is bi-Lagrangian.
* b is diagonal on the basis l_i, so its radical R = L ∩ SL is spanned by the eigenvector l_i's.
* q(l_i) ≠ 0 for every cyclic plane. So (b) holds.

**Reduction.**
* The criterion is additive over orthogonal S-invariant sums, because q|L is the sum of its pieces.
* The lower bound A ≥ n − c holds in every characteristic. A split contains at most c invariant planes, and every other
  plane meets its image in at most a line.
* Take a maximal orthogonal family of c invariant planes. Its complement W has no invariant plane. So W decomposes into
  orthogonally indecomposable S-modules of dimension ≥ 4, and the law follows once each such piece meets the criterion.

## 2. Every indecomposable piece meets the criterion

In the ring models below, σ is the involution x ↦ x⁻¹, and ℓ is a linear functional whose kernel contains no nonzero
ideal.

| piece | construction | radical | q ≠ 0 because |
|---|---|---|---|
| U ⊕ U*, U = F[x]/(fᵐ), f ≠ f* | graph of multiplication by a unit g in the model R ⊕ R, with S = (x, x⁻¹) and ω = ℓ(ab′) + ℓ(a′b) | 0 (1 − x⁻² is a unit) | q = ℓ(τ g a²) with τ = x + x⁻¹; units span R, so some g has ℓ(τg) ≠ 0 |
| F[x]/(fᵐ), f = f* ≠ x+1 (R/R₀ étale) | L = R₀η | 0 (R = R₀ ⊕ R₀ξ) | ℓ(R₀) = 0 because the trace is onto; q = λ(a² N(η)) with λ(m) = ℓ(mξ), the norm is onto R₀^×, and units span R₀ |
| **V(2k)**, one Jordan block, k ≥ 3 | R = F[t]/(t^{2k}), σ(t) = t/(1+t), s = tσ(t) = t²/(1+t), P = F[s], L = Pη | 0 (P ∩ Pt = 0; valuations in P are even) | q(aη) = λ(a² N(η)) with λ(m) = ℓ(mt); the span of {a² N(η)} contains every sʲ with j ≠ 1, and λ(s^{k−1}) = ℓ(t^{2k−1}) ≠ 0 |
| **W(k)** = J_k ⊕ J_k*, k ≥ 4 | L_g = P(1,0) + P(t,g) in R_k² | ann(g): 0 (k even, g a unit) or F(t^{k−1},0), a fixed vector (k odd, v(g) = 1) | see below |
| V(4), W(2), W(3) | exhaustive search (15, 15 and 135 Lagrangians) | 1, 2 (= L, fixed), 1 | found |

**W(k) in detail.**
* L_g is automatically bi-Lagrangian: it is a P-module and S + S⁻¹ = s ∈ P.
* L_g is Lagrangian exactly when y = σ(g) ∈ P^⊥, where ⊥ is taken for the pairing ⟨u,v⟩ = ℓ(uv).
* On L_g, q = ℓ(c² s t y). So q ≢ 0 exactly when y ∉ Q^⊥, with Q = span{c² s t}.
* A suitable y exists because a vector space over any field is not the union of two proper subspaces. The two needed
  failures of containment come from valuation parity:
  * k even: P^⊥ ⊄ tR, because (tR)^⊥ = ann(t) = F t^{k−1} and t^{k−1} ∉ P (odd valuation);
  * k odd: P^⊥ ⊄ t²R, because t^{k−2} ∉ P;
  * both parities: P^⊥ ⊄ Q^⊥, because st has valuation 3 < k, so st ∉ P.
* The certificate checks these three containments for k = 4…12 and builds the split each time.

**Theorem (qubit arrow law).** For every n and every Clifford tick S ∈ Sp(2n,2), A(S) = n − c(S).

The table and the reduction cover every orthogonally indecomposable symplectic isometry in characteristic 2:
* pieces for f ≠ x+1 follow Wall (1963); the unitary-type pieces need only an element of trace one, which exists because
  R/R₀ is étale;
* the unipotent pieces V(2k) and W(k) follow Hesselink (1980), and Liebeck–Seitz, *Unipotent and Nilpotent Classes in
  Simple Algebraic Groups and Lie Algebras* (AMS 2012), Lemma 6.2.

The constructions for V(2k), W(k) and the f ≠ x+1 pieces work over every field of characteristic 2. The three finite
exceptions are certified over F₂ only, so the theorem is stated for qubits.

## 3. What the certificate checks

* **Constructions, with the split built and verified plane by plane.** In every case the split has n nondegenerate,
  mutually orthogonal planes, each meeting its image.
  * V(2k) for k = 3…7, over **every** admissible form (4, 8, 16, 32, 64 functionals ℓ).
  * W(k) for k = 4…12.
  * f ≠ f*: x³+x+1 (m = 1, 2), x⁴+x+1, x⁵+x²+1 and x⁴+x³+1 (m = 2).
  * f = f* ≠ x+1: x²+x+1 (m = 2, 3, 4), x⁴+x³+x²+x+1 (m = 1, 2) and x⁶+x³+1.
  * The exceptions V(4), W(2) and W(3).
* **Verifier control.** 200 random splits of W(4) are all rejected by the same verifier.
* **Classes.**
  * For every class of Sp(2n,2) with n ≤ 4 (11, 30 and 81 classes), A = n − c.
  * Every class with c = 0 (5, 7 and 22 classes) carries a criterion Lagrangian, and its split is built and verified.
* **Qubit closed form for c**, exact on every class for n ≤ 5 (11 + 30 + 81 + 198 = 320 classes):

  c(S) = m₁(1)/2 + χ₂·m₂(1) + m₁(x²+x+1),

  where:
  * m_s(f) is the number of Jordan blocks of size s at f;
  * χ₂ = 1 exactly when some v ∈ ker(S+1)² has ω(v,(S+1)v) = 1, i.e. the size-2 blocks can be taken as transvection
    planes V(2). This is Hesselink's index at 2.

  The closed form is verified (Tier C), not proved here. The law itself does not depend on it: c is defined geometrically.

## 4. Reading

* The law is now a theorem for every prime local dimension: the finite-field symplectic groups in odd characteristic
  (Pass 11207) and in characteristic 2 (this pass).
* For qubits, the protected planes are of three kinds:
  * identity planes, two 1-blocks;
  * transvection planes, one 2-block, allowed only when χ₂ = 1;
  * order-3 planes, x²+x+1.
* The unavoidable forgotten record is one bit per unprotected qubit per tick. Pass 11207's Landauer reading therefore
  carries over with k_B T ln 2 in place of k_B T ln 3.
* **How much a random qubit tick protects.** From Pass 11208's class census (E[c] = n − E[A], exact):

  | n | E[c] |
  |---|---|
  | 2 | 29/40 |
  | 3 | 5981/8960 = 0.6675 |
  | 4 | 12945139/19496960 = 0.66396 |
  | 5 | 6791923649683/10212039720960 = 0.66509 |

  * Uniform samples scored with the closed form give 0.654 ± 0.012 at n = 5 (the exact value is 0.6651), 0.653 ± 0.014
    at n = 8, 0.665 ± 0.017 at n = 12, 0.674 ± 0.020 at n = 16, and 0.654 ± 0.028 at n = 24.
  * Qubits therefore protect about 0.665 planes per tick, fewer than the qutrits' 0.72553535… (master's Pass 11215).
  * P(c = 0), the maximal arrow, is ≈ 0.49, against ≈ 0.45 for qutrits.
  * The exact qubit limit needs Fulman's cycle index refined by Hesselink's index. That is left open.
* The obstruction that blocked the odd-characteristic argument (b alternating on V(4) and W(2)) is real: those pieces
  have no transverse non-alternating bi-Lagrangian. It is bypassed, not removed:
  * W(2) uses its fixed Lagrangian (L = L ∩ SL);
  * V(4) uses a radical of one fixed vector.

**Prior art checked.**
* `RESULTS_INDEX.md` and the analysis notes contain no qubit proof. Pass 11208 lists it as open, and master's Passes
  11213–11216 do not touch characteristic 2.
* The individual ingredients are classical:
  * Hesselink's V/W normal forms;
  * Wall's classification;
  * diagonalisation of non-alternating symmetric forms (Albert).
* The criterion and the F[s]-lattice constructions are new here as far as the corpus shows.
