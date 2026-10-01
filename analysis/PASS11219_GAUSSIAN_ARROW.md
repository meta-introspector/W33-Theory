# Pass 11219: the arrow law over the reals (Gaussian dynamics)

Producer: `analysis/w33_pass11219_gaussian_arrow.py`
Certificate: `data/w33_pass11219_gaussian_arrow.json`
Regression: `tests/test_w33_pass11219_gaussian_arrow.py`

**Background.**
* A(S) = n − c(S) is now proved for every finite symplectic group: Pass 11207 for odd q, Pass 11217 for q = 2.
* Every step of the odd-characteristic proof works over any field of characteristic ≠ 2:
  * the characteristic-free lower bound;
  * the Wall/Milnor orthogonal decomposition;
  * the bi-Lagrangian radical lemma;
  * Taussky–Zassenhaus graphs for f ≠ f* (which include real λ ≠ ±1 and loxodromic quartets);
  * the Cayley transform for V(2m);
  * the eigenvector reduction for W(k), k odd;
  * the Hermitian real form for unitary-type blocks (over ℝ: e^{±iθ}, κ = i);
  * diagonalisation of a symmetric form.
* So the law holds for **S ∈ Sp(2n, ℝ)**, the symplectic matrices of Gaussian unitaries, i.e. of free bosonic dynamics.
  This pass states what the law says there and checks it numerically.

## What the arrow is for Gaussian dynamics

For a mode decomposition V = P₁ ⊕ … ⊕ P_n into symplectic planes,

    dim(P ∩ SP) = 2 − rank S[P^⊥ ← P].

Here S[P^⊥ ← P] is the block of S that couples mode P to the other modes. So

    A(S) = min over mode decompositions of Σ_k rank S[P_k^⊥ ← P_k]

is the smallest number of **quadrature channels** through which the modes must leak into one another per tick. The law
gives:

    A(S) = n − c(S),   c(S) = #(elliptic size-1 pairs e^{±iθ}) + #(hyperbolic size-1 pairs λ, 1/λ) + Σ_{±1}(m₁/2 + m₂).

The elliptic and hyperbolic pairs are the normal modes, and c counts them. For semisimple S without eigenvalues ±1:

    **A(S) = 2 × (number of loxodromic quartets r e^{±iθ}, r⁻¹ e^{±iθ}).**

* **Stable free dynamics has no intrinsic arrow.** A positive quadratic Hamiltonian has only elliptic modes, so c = n and
  A = 0. The Williamson normal modes are exactly the subsystems the dynamics keeps to itself.
* **Pure squeezing has none either.** Hyperbolic pairs are invariant planes, so A = 0 even though their
  Kolmogorov–Sinai entropy is positive.
* **The arrow lives in complex instability and Jordan defects.**
  * Each loxodromic quartet costs exactly two channels. Its best split couples each of its two modes to the other by a
    rank-one block, never more.
  * So does each non-semisimple elliptic or hyperbolic pair, and V(4) at ±1.
* **Krein's theorem, read as a statement about the arrow.**
  * Two elliptic modes of *opposite* Krein signature collide at a Hamiltonian–Hopf point. There the pair is
    non-semisimple, and A jumps 0 → 2. It stays 2 on the loxodromic side.
  * Modes of the *same* signature never destabilise, and an exact same-signature degeneracy is semisimple, so A = 0
    throughout.

## Checks (`data/w33_pass11219_gaussian_arrow.json`)

1. **Generic mixtures.**
   * Setup: 33 random S = M S₀ M⁻¹, with M a random symplectic matrix and S₀ built from e elliptic, h hyperbolic and l
     loxodromic blocks, n ≤ 6 (types (2,0,0), (0,2,0), (1,1,0), (0,0,1), (1,0,1), (0,1,1), (2,1,1), (0,0,2), (1,1,2),
     (0,0,3) and (2,2,1), three draws each).
   * Construction: invariant planes plus a Taussky–Zassenhaus split for each quartet.
   * In every case the result is a symplectic split whose total coupling rank is **exactly 2l = n − c**.
2. **Non-semisimple pieces.**
   * Pieces: the elliptic Jordan pair at a Krein collision (Meyer–Hall normal form), V(4) at +1
     (H = p₁x₂ + p₂²/2), and a hyperbolic Jordan pair.
   * Method: a transverse bi-Lagrangian is found numerically (least squares on ω(·,(S+S⁻¹)·) = 0 over graph
     Lagrangians), and b is diagonalised.
   * Each split has total coupling rank 2 = n − c.
3. **Krein families.** Modes with ω₁ = 1 and ω₂ = 1.2, coupled by ε(x₁x₂ − p₁p₂):
   * with opposite signature the pair is semisimple (A = 0) for ε < 0.1, collides at ε = |ω₁ − ω₂|/2 = 0.1, and is
     loxodromic (A = 2) beyond it, with growth rate √(ε² − 0.01) to 10⁻⁹;
   * with equal signature it is semisimple (A = 0) for every ε tested up to 0.5;
   * at the exact degeneracies, the opposite-signature point is non-semisimple (A = 2) and the equal-signature one is
     semisimple (A = 0).

The lower bound A ≥ n − c is the characteristic-free argument of Pass 11207, so these constructions give equality.

## Scope and what is not claimed

* The theorem is about **dimensions**: coupling ranks of the symplectic matrix, a property of the operator.
* A first attempt to read it as entanglement growth of the vacuum was dropped, because it is wrong as stated. For
  S = diag(A, A^{−T}) with A = r·R(θ), the vacuum stays a product state in the standard split even though S is
  loxodromic: equal squeezing commutes with the passive rotation. An entropic version would have to quantify over
  initial states or use operator entanglement. It is left open.
* Bianchi, Hackl and Yokomizo (JHEP 03 (2018) 025) relate entanglement growth in unstable quadratic systems to Lyapunov
  exponents and are the natural comparison for such a version.
* Classical background, cited rather than rederived:
  * Williamson normal modes;
  * Krein and Gelfand–Lidskii strong stability;
  * Laub–Meyer normal forms of real symplectic matrices.
* Within the corpus, the statement that the minimal number of inter-mode channels equals n − c, with each loxodromic
  quartet costing exactly two, is new. `RESULTS_INDEX.md` and the notes have no Gaussian or Krein version of the arrow;
  the "Krein" hits are association-scheme Krein parameters.

## Reading

* The arrow law is the same theorem for the finite (Clifford) substrate and for its continuum (Gaussian) shadow:
  symplectic linear algebra over any field.
* On the substrate a random tick is almost never decomposable: c is about 0.73 for qutrits and about 0.665 for qubits,
  out of n.
* In the continuum, stable free fields are fully decomposable, and the arrow needs complex instability.
