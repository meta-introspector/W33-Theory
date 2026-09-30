# Pass 11190 — every AME(10,3) stabilizer state carries one sign structure; the five-qutrit conjecture is proved

Producer: `analysis/w33_pass11190_ame10_sign_structure.py`
Certificate: `data/w33_pass11190_ame10_sign_structure.json`
Regression: `tests/test_w33_pass11190_ame10_sign_structure.py`

**What was open.** Pass 11186 found that the 17 892 perfect five-qutrit Clifford gates from 71 AME(10,3) graph states
realise only two of the five block-determinant patterns the determinant law admits: 10-cycle and cross+permutation. It
conjectured that C4+C6, double cross and all-(−1) never occur.

**1. Signs (proved).** For a 4-set K of the ten parties:
* The stabilisers trivial on K form a 2-dim space W_K.
* Restriction to each of the six other parties is an isomorphism, since anything also trivial there has weight ≤ 5.
* Pulling back each party's symplectic form gives signs c_p = ±1 with Σ_p c_p ≡ 0 (mod 3) (isotropy).
* So every 6-set of parties is **uniform (6|0) or split 3|3**.
* A gate block determinant is d_ij = −c_i c_j, computed in U = outputs ∪ {j}. Checked against all 17 892 gates.

**2. Local graphs (proved).** For a 3-set A:
* W_A is 4-dim. The seven restriction kernels are pairwise skew lines of PG(3,3), and their Plücker points satisfy
  Σ λ_p ℓ_p = 0 on the Klein quadric.
* G_pq = λ_pλ_q B(ℓ_p, ℓ_q) is a symmetric ±1 matrix whose row d is the split of the 6-set A^c ∖ {d}.
* All degrees of its (−1)-graph lie in {0, 3, 6}. An exhaustive search over 2²¹ graphs finds 1052 such graphs, which
  up to complement are empty, K4+3K1, K33+K1 and prism+K1.
* Empty and K33+K1 have rank 6 with an *elliptic* O⁻(6,3) Gram matrix, impossible inside the hyperbolic Klein quadric
  O⁺(6,3). K4+3K1 is O⁺; prism+K1 has rank 4.

**3. Duality (proved).** For disjoint 4-sets K₁, K₂ with K₁ ∪ K₂ = complement of {i, j}, isotropy between W_{K₁} and
W_{K₂} gives det φ^{K₁}_{j←i} = det φ^{K₂}_{j←i}. There are 0 violations in 111 825 checks.

**4. CP-SAT (exact).** The model has 210 six-set states (11 values each), the 120 local-graph tables and the 1575
duality equalities.
* "Some local graph is K4+3K1" is **INFEASIBLE**, also with a single worker and K4 forced at A = {0,1,2}. S₁₀ is
  transitive on 3-sets, so forcing one triple loses nothing.
* "prism+K1 everywhere" is feasible.

**Theorem.** In every AME(10,3) stabilizer state, every local graph is a prism plus an isolated vertex. So each 3-set
of parties has exactly one uniform extension, and **the 30 uniform 4-sets form a Steiner system S(3,4,10)**, the
inversive plane of order 3.

**Corollary (the conjecture).**
* A 5-set contains at most one block of an S(3,4,10), so a gate has at most one all-(−1) column. Double cross and
  all-(−1) are impossible.
* In a prism, two columns with two (−1)s never share both, so there is no 4-cycle and C4+C6 is impossible.
* Hence **cross+permutation ⇔ the input set contains a block** (180 of 252 splits), and a **10-cycle otherwise** (72).
  This is exactly the 12 780 / 5 112 census.
* CP-SAT confirms directly: C4+C6, double cross and all-(−1) are INFEASIBLE, while C10 and cross+perm are feasible.

**Uniqueness.**
* A fixed S(3,4,10) admits exactly **2** sign structures.
* The sign structure of the first graph has automorphism group of order **720**, with 81 involutions and elements of
  order 8 and 10 but none of order 6. That is **PGL(2,9)**, of index 2 in the PΓL(2,9) (order 1440) of its Steiner
  system, which therefore swaps the two structures.
* So every AME(10,3) stabilizer state carries **the same** sign structure up to relabelling, and all 71 graphs check
  isomorphic. The ten parties behave as the projective line over F₉.

**Scope.** Stabilizer (graph) AME(10,3) states only. Non-stabilizer AME(10,3) states are not covered.

**Prior art.**
* AME(n,q) ↔ [[n,0,n/2+1]]_q codes and graph-state constructions are standard (Helwig et al.; Huber–Wyderka table;
  Grassl tables). The GF(9) MDS construction is Pass 11175.
* Context in the corpus: the inversive plane of order 3 is the Miquelian plane of a regular spread of PG(3,3)
  (J. A. Thas, *Symplectic spreads in PG(3,q), inversive planes and projective planes*, 1997). The repo's
  regular-spread packets are `analysis/BT2088_BT2092_regular_spread_quadratic_structure_controller.md` and Pass 2064.
  Whether this S(3,4,10) is that spread's plane, i.e. whether the ten parties can be identified with ten spread lines
  of PG(3,3) compatibly with the local skew-line configurations of step 2, is **open**. No such map is built here.
* A limited search (3 web queries) found no link between AME sign structures and S(3,4,10), inversive planes or
  PGL(2,9). The dedicated prior-art agent for this topic hit a session limit, so this search is weaker than usual.
