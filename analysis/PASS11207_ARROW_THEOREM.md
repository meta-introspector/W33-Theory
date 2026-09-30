# Pass 11207 — the intrinsic arrow of time for every number of qutrits: A(S) = n − c(S)

Producers:
* `analysis/w33_pass11207_arrow_theorem.py` (theorem checks);
* `analysis/w33_pass11207_arrow_by_planes.py` (part 1: plane formula; exact n ≤ 4);
* `analysis/w33_pass11207_half_moving_splits.py` (part 2: n = 5, 6 search; semisimple case).

GAP: `analysis/gap/w33_pass11207_class_reps{,_n4,_n5}.g`
Frozen: `data/w33_pass11207_{arrow_theorem,arrow_by_planes,half_moving_splits}.json` and `data/w33_pass11207_gap_class_reps*.txt`
Regression: `tests/test_w33_pass11207_{arrow_theorem,arrow_by_planes,half_moving_splits}.py`

**Numbering.** This work was first committed on branch `claude/gallant-sagan-i33c7e` as Passes 11188, 11191 and 11199,
from 19:10 UTC on 2026-09-30. Master then received the parallel session's own Passes 11188–11192 (20:49 UTC) and a
reservation of 11199–11206. Master is canonical, so the branch material is renumbered 11207–11211 and cites master's
passes wherever the results overlap.

**What master already has.**
* Pass 11183 defined the intrinsic arrow A(S): the loss that no choice of subsystems avoids, i.e. the minimum over all
  splits of the number of trits exported.
* Pass 11188 (master) computed A exactly on every class for n = 2, 3: A ∈ {0, 2} and {0, 2, 3}, never 1 and never more
  than n.
* Its reserved item "a second law from counting (arrow-free fraction vs n)" is answered exactly by the theorem below:
  the arrow-free fraction is P(c = n).

**What this pass adds.**
* The plane formula A(S) = 2n − max Σ_P dim(P ∩ SP) over orthogonal families of nondegenerate planes (part 1). It gives
  exact A for n = 4 (all 278 classes).
* The general lower bound, A ≤ n − 1 ⇒ an invariant plane.
* The theorem below, for every n.

## Theorem

Let **c(S)** be the largest number of mutually orthogonal, S-invariant, nondegenerate planes: the qutrits the dynamics
can keep to themselves. For every n and every S ∈ Sp(2n, 3),

    A(S) = n − c(S).

Consequences:
* **A ≤ n.** The unavoidable arrow is at most one trit per qutrit.
* **A = n** exactly when S fixes no nondegenerate plane.
* **A is never 1.** Here c = n − 1 forces c = n, because the last plane is the orthogonal complement of the others.
* **The arrow counts the qutrits that are not dynamically protected**: each qutrit the tick can keep invariant saves one
  trit, and each one it cannot costs exactly one trit.

**Reading: half the record.** In an optimal split, each protected qutrit exports nothing. Each unprotected qutrit exports
**exactly one** of its two trits, since its plane meets its own image in exactly a line. Theorem 4.7 of the paper prices
a forgotten record at two trits. The intrinsic price is half that: one trit, which by Landauer's principle is
k_BT ln 3, per unprotected qutrit per tick. (Bridge reading; the mathematics is the theorem.)

## Proof

**Lower bound.** In any split, let a be the number of planes that are invariant (d = 2), so a ≤ c; the rest have
d ≤ 1. Then E = 2n − Σd ≥ 2n − 2a − (n − a) = n − a ≥ n − c.

**Upper bound.**

*Step 1: reduction.* Split off a maximum orthogonal family of c invariant planes; each exports 0. The complement W is
S-invariant and has no invariant nondegenerate plane. A is subadditive over orthogonal S-invariant pieces, since splits of
the pieces combine. So it is enough to show that every tick without an invariant nondegenerate plane has a
**half-moving split**, one in which every plane meets its image, so that E = n.

*Step 2: decomposition.* W splits orthogonally and S-invariantly (Wall 1963; Milnor 1969):
* by irreducible factors f of the characteristic polynomial, into V_f ⊕ V_{f*} when f ≠ f* (reciprocal), and into V_f
  when f = f*;
* a component with f = x ∓ 1 is ±u with u unipotent symplectic. In odd characteristic it splits into blocks V(2m) (one
  Jordan block of even size) and W(k) (two dual Jordan blocks of odd size k);
* a component with f = f* of degree 2d ≥ 2 is a Hermitian space over E = F_{3^{2d}} ⊃ F = F_{3^d}. There S = μU, with
  μ^{q+1} = 1, μ ∉ F, and U unipotent unitary, and it splits into single unitary Jordan blocks.

Blocks of dimension 2 (V(2), W(1), and the Hermitian block with d = k = 1) are invariant planes, so they do not occur in W.

*Step 3: the bi-Lagrangian / radical lemma.* Let T = S + S⁻¹, which is ω-self-adjoint, and let L be a Lagrangian for both
ω and ω_T = ω(·, T·).
* Then b(x, y) = ω(x, Sy) is symmetric on L, with radical L ∩ SL.
* If L ∩ SL = 0: diagonalize b to get eᵢ with b(eᵢ, eᵢ) ≠ 0. The planes ⟨eᵢ, Seᵢ⟩ are nondegenerate and mutually
  orthogonal, and each meets its image.
* If L ∩ SL = ⟨e⟩ with Se = ±e: diagonalize b modulo e. The resulting n − 1 planes are orthogonal to e. The plane
  complementary to them contains e = ±Se, so it meets its image too.

*Step 4: a construction for each block.*

| block | Lagrangian L | why it works |
|---|---|---|
| **(δ)** V_f ⊕ V_{f*}, f ≠ f*: S = A ⊕ A^{−T} on U ⊕ U* | L = {(x, Φx)}, Φ symmetric with ΦA = AᵀΦ (Taussky–Zassenhaus) | L ∩ SL = ker(A² − 1) = 0, since ±1 are not roots of f |
| **(α)** V(2m) | N = (u−1)(u+1)⁻¹ is skew and nilpotent; L = ⟨v, N²v, …, N^{2m−2}v⟩, v cyclic | isotropic (ω(v, N^{even}v) = 0); T = u + u⁻¹ is even in N, so L is T-stable; L ∩ uL = 0 by parity (uf ∈ L forces Nf ∈ L ∩ NL = 0, and ker N ∩ L = 0) |
| **(β)** W(k), k odd | reduce by the eigenvector e ∈ U: e^⊥/e ≅ J_{k−1} ⊕ J_{k−1}^{−T}, where the parity Lagrangian ⟨N^{even}v, N^{even}w⟩ works | its lift has L ∩ uL = ⟨e⟩; apply the radical lemma. Without the eigenvector plane it is impossible: all 972 optimal splits of W(3) use one |
| **(γ)** μ·J_k, unitary | the real form V_R = F-span{εⱼ Nʲv}, with εⱼ = 1 (j even) and κ (j odd), κ̄ = −κ | the skew-Hermitian form is imaginary on V_R, so V_R is ω-Lagrangian; T maps V_R into itself; V_R ∩ SV_R = 0 because μ ∉ F |

∎

## Closed form

Jordan types add over orthogonal S-invariant sums. An invariant nondegenerate plane is one of three 2-dimensional
indecomposables:
* a pair of size-1 Jordan blocks at λ = ±1;
* one size-2 block at λ = ±1;
* one size-1 F₉-block of the factor x² + 1.

So

    c(S) = Σ_{λ=±1} ( m₁(λ)/2 + m₂(λ) ) + m₁^{F₉}(x² + 1),     A(S) = n − c(S),

where m_s is the number of Jordan blocks of size s, read off from the ranks of (S − λ)^s and (S² + 1)^s. **The intrinsic
arrow of any qutrit Clifford dynamics is computable in polynomial time** (`arrow(S, n)` in the producer). The formula is
checked on all 372 classes (n ≤ 4) and on all 240 random ticks at n = 5, 6.

## Checks (all in the producer)

* **All conjugacy classes, n ≤ 4.** A = n − c(S) on 20/20, 74/74 and 278/278 classes of PSp(4,3), PSp(6,3), PSp(8,3), with A
  taken from exhaustive enumeration over all splits (part 1 for n = 4; n = 2, 3 agree with master's Pass 11188).
* **Every construction built and checked**, testing nondegenerate, mutually orthogonal planes each meeting its image:
  * V(2m) for m ≤ 8, both form classes;
  * W(k) for k ≤ 11;
  * unitary blocks over F₉ (k ≤ 5) and F₈₁ (k ≤ 3);
  * f ≠ f* pairs x²+x+2, (x²+x+2)², (x²+x+2)³, x³+2x+1, (x³+2x+1)².
* **Random ticks at n = 5, 6.** Take a maximum invariant family, then find a half-moving split of its complement. This
  gives a split with E = n − c, which meets the lower bound, so A is computed exactly without using the theorem. It
  equals n − c on every sample.

**What this closes.**
* The conjecture "A ≤ n for all n", which master's Pass 11188 proves for n = 2, 3.
* The semisimple special case, first proved in the branch's part 2.

A(S) = n − c(S) is now an exact formula for the intrinsic arrow of any qutrit Clifford dynamics, together with a
polynomial-time closed form.

**Scope and prior art.**
* The structure theory used is classical: orthogonal decomposition of isometries of symplectic spaces by characteristic
  polynomial, and unipotent classes of Sp and U in odd characteristic (Wall, *J. Austral. Math. Soc.* 3 (1963); Milnor,
  *Topology* 8 (1969); Springer–Steinberg).
* Taussky–Zassenhaus (*Pacific J. Math.* 9 (1959)): every square matrix is similar to its transpose via a symmetric matrix.
* The arrow A is Pass 11183; exact n = 2, 3 is master's Pass 11188; c(S), the plane formula and the all-n theorem are this pass. The operational-mereology idea of minimising scrambling over subsystem structures is Zanardi, Dallas, Andreadakis & Lloyd, Quantum 8, 1406 (2024), cited in Pass 11188.
* The proof is written for F₃. Every step uses only odd characteristic, so it extends verbatim to qudits of any odd prime
  power dimension. That extension is not checked numerically here.
